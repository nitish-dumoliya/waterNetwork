import time
import sys
from tabulate import tabulate
from amplpy import AMPL
import contextlib
import os
import numpy as np
import math

class WaterNetworkSolver:
    def __init__(self, model_file, solver_name, data_file, data_number):
        self.model_file = model_file
        self.solver_name = solver_name
        self.data_file = data_file
        self.data_number = data_number
        self.ampl = AMPL()
        
        # To store solutions
        self.q_init = {}
        self.h_init = {}
        self.l_init = {}
        self.eps_init = {}
        self.q1_init = {}
        self.q2_init = {}
        self.eps = {}
        
    def read_model_and_data(self):
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)
        
        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')
        
        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
 

    def compute_adaptive_eps(self, demand):
        """
        Adaptive epsilon calculation based on the minimum demand value.
        """
        epsilon = 1e-14 + 1e-10*(demand)**2
        #if max_demand <= 0.1:
            #    epsilon =  1e-13
            #elif max_demand >= 0.1 and max_demand <= 0.2:
            #epsilon = 1e-14
        #elif max_demand >= 0.2:
        #    epsilon = 1e-6
        
        return epsilon

    def constraint_violations(self, q_values, h_values, l_values, epsilon, solver):
        total_absolute_constraint_violation = 0
        total_relative_constraint_violation = 0
         
        con1_gap = {}
        if self.data_number==5:
            q1 = self.ampl.get_variable('q1').get_values().to_dict()
            q2 = self.ampl.get_variable('q2').get_values().to_dict()
            for i in self.nodes:
                if i not in self.source:
                    con1_rhs = self.D[i]
                    incoming_flow = sum(q1[j, i] + q2[j,i] for j in self.nodes if (j, i) in self.arcs)
                    outgoing_flow = sum(q1[i, j] + q2[i,j] for j in self.nodes if (i, j) in self.arcs)
                    con1_lhs = incoming_flow - outgoing_flow
                    con1_violation = con1_lhs - con1_rhs
                    con1_gap[f"{i}"] = con1_violation
                    total_absolute_constraint_violation += abs(con1_violation)
        else:
            for i in self.nodes:
                if i not in self.source:
                    con1_rhs = self.D[i]
                    incoming_flow = sum(q_values[j, i] for j in self.nodes if (j, i) in self.arcs)
                    outgoing_flow = sum(q_values[i, j] for j in self.nodes if (i, j) in self.arcs)
                    con1_lhs = incoming_flow - outgoing_flow
                    con1_violation = con1_lhs - con1_rhs
                    con1_gap[f"{i}"] = con1_violation
                    total_absolute_constraint_violation += abs(con1_violation)
        #print("con1_gap:", con1_gap) 
        con2_original_gap = {}
        con2_approx_gap = {}
        absolute_violations = {}
        relative_violations = {}
        con2_absolute_constraint_violation = 0
        con2_relative_constraint_violation = 0
        con2_original_violation = 0
        con2_approx_violation = 0

        if self.data_number==5:
            #q1 = self.ampl.get_variable('q1').get_values().to_dict()
            #q2 = self.ampl.get_variable('q2').get_values().to_dict()
            self.exdiam = self.ampl.getParameter('exdiam').to_dict()
            for (i, j) in q1.keys():
                # Original constraint value
                lhs = 2*(h_values[i] - h_values[j])
                alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[i,j] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                original_rhs = q1[i, j] * (abs(q1[i, j])) ** 0.852 * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + q2[i, j] * (abs(q2[i, j])) ** 0.852 * alpha_rhs  
                # original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs
                # Approximated constraint value
                # approx_rhs = (q1[i, j]**3 * ((q1[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q1[i,j]**2 + 0.426*epsilon[i,j]**2)) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j]**3 * ((q2[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q2[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                approx_rhs = (q1[i, j] * ((q1[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j] * ((q2[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

                #approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs
                con2_original_gap[f"{i},{j}"] = lhs - original_rhs
                con2_original_violation += abs(lhs - original_rhs) 
                con2_approx_gap[f"{i},{j}"] = lhs - approx_rhs 
                total_absolute_constraint_violation += abs(lhs - approx_rhs)    
                con2_approx_violation += abs(lhs - approx_rhs) 

                # Compute absolute violation
                absolute_violation =  original_rhs - approx_rhs
                absolute_violations[f"{i},{j}"] = absolute_violation
                con2_absolute_constraint_violation += abs(absolute_violation)

                # Compute relative violation between original_rhs and approx_rhs
                relative_violation = (original_rhs - approx_rhs) / (original_rhs+1e-14)
                relative_violations[f"{i},{j}"] = relative_violation
                con2_relative_constraint_violation += abs(relative_violation)

        elif self.data_number==6:
            self.fixarcs = self.ampl.getSet('fixarcs')
            self.fix_r = self.ampl.getParameter('fix_r').to_dict()
            self.exdiam = self.ampl.getParameter('fixdiam').to_dict()
            for (i, j) in self.arcs:
                if (i,j) not in self.fixarcs:
                    # Original constraint value
                    lhs = h_values[i] - h_values[j]
                    alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                    #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                    original_rhs = q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs  
                    #original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs 
                    # Approximated constraint value
                    # approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                    approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

                    #approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs
                    con2_original_gap[f"{i},{j}"] = lhs - original_rhs
                    con2_approx_gap[f"{i},{j}"] = lhs - approx_rhs

                    total_absolute_constraint_violation += abs(lhs - approx_rhs)    
                    con2_original_violation += abs(lhs - original_rhs) 
                    con2_approx_violation += abs(lhs - approx_rhs) 

                    # Compute absolute violation
                    absolute_violation =  original_rhs - approx_rhs
                    absolute_violations[f"{i},{j}"] = absolute_violation
                    con2_absolute_constraint_violation += abs(absolute_violation)

                    # Compute relative violation between original_rhs and approx_rhs
                    relative_violation = (original_rhs - approx_rhs) / (original_rhs+1e-14)
                    relative_violations[f"{i},{j}"] = relative_violation
                    con2_relative_constraint_violation += abs(relative_violation)

            #print("con2_gap:", con2_gap)
            for (i, j) in self.fixarcs:
                # Original constraint value
                lhs = h_values[i] - h_values[j]
                alpha_rhs = 10.67 * self.L[i, j] / ((self.fix_r[i,j] ** 1.852) * ((self.exdiam[i,j]) ** 4.87))
                #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                original_rhs = q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * 10.67 * self.L[i,j]/(self.fix_r[i,j]**1.852 * self.exdiam[i,j]**4.87) 
                #original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs

                # Approximated constraint value
                # approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * 10.67 * self.L[i,j]/(self.fix_r[i,j]**1.852 * self.exdiam[i,j]**4.87) 
                approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * 10.67 * self.L[i,j]/(self.fix_r[i,j]**1.852 * self.exdiam[i,j]**4.87) 

                #approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs

                #con2_original_violation =  lhs - original_rhs
                con2_original_gap[f"{i},{j}"] = lhs - original_rhs
                con2_original_violation += abs(lhs - original_rhs) 

                #con2_approx_violation =  lhs - approx_rhs
                con2_approx_gap[f"{i},{j}"] = lhs - approx_rhs

                total_absolute_constraint_violation += abs(lhs - approx_rhs)    
                con2_approx_violation += abs(lhs - approx_rhs) 

                # Compute absolute violation
                absolute_violation =  original_rhs - approx_rhs
                absolute_violations[f"{i},{j}"] = absolute_violation
                con2_absolute_constraint_violation += abs(absolute_violation)

                # Compute relative violation between original_rhs and approx_rhs
                relative_violation = (original_rhs - approx_rhs) / (original_rhs+1e-14)
                relative_violations[f"{i},{j}"] = relative_violation
                con2_relative_constraint_violation += abs(relative_violation)
        else:
            for (i, j) in q_values.keys():
                # Original constraint value
                lhs = h_values[i] - h_values[j]
                alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs
                
                # Approximated constraint value
                # approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

                #approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs

                con2_original_gap[f"{i},{j}"] = lhs - original_rhs
                con2_original_violation += abs(lhs - original_rhs) 
                
                con2_approx_gap[f"{i},{j}"] = lhs - approx_rhs
                
                total_absolute_constraint_violation += abs(lhs - approx_rhs)    
                con2_approx_violation += abs(lhs - approx_rhs) 

                 # Compute absolute violation
                absolute_violation =  original_rhs - approx_rhs
                absolute_violations[f"{i},{j}"] = absolute_violation
                con2_absolute_constraint_violation += abs(absolute_violation)

                # Compute relative violation between original_rhs and approx_rhs
                relative_violation = (original_rhs - approx_rhs) / (original_rhs + 1e-14)
                relative_violations[f"{i},{j}"] = relative_violation
                con2_relative_constraint_violation += abs(relative_violation)
           
        #print("con2_gap:", con2_gap)
        
        con3_gap = {}
        if self.data_number==6:
            self.fixarcs = self.ampl.getSet('fixarcs')
            for (i,j) in self.arcs:
                if (i,j) not in self.fixarcs:
                    con3_rhs = self.L[i,j]
                    con3_lhs = sum(l_values[i,j,k] for k in self.pipes) 
                    con3_violation = con3_lhs - con3_rhs
                    con3_gap[f"{i},{j}"] = con3_violation 
                    total_absolute_constraint_violation += abs(con3_violation)
        else:
            for (i,j) in self.arcs:
                con3_rhs = self.L[i,j]
                con3_lhs = sum(l_values[i,j,k] for k in self.pipes) 
                con3_violation = con3_lhs - con3_rhs
                con3_gap[f"{i},{j}"] = con3_violation 
                total_absolute_constraint_violation += abs(con3_violation)
        #print("con3_gap:", con3_gap)

        con4_gap = {}
        for (i,j) in self.arcs:
            for k in self.pipes:
                #con4_rhs = self.L[i,j]
                #con4_lhs = l_values[i,j,k]
                con4_violation = max(0,l_values[i,j,k]-self.L[i,j])
                con4_gap[f"{i},{j},{k}"] = con4_violation 
                total_absolute_constraint_violation += abs(con4_violation)
        #print("con4_gap:", con4_gap)
        
        con5_gap = {}
        for j in self.source:
            con5_rhs = self.E[j]
            con5_lhs = h_values[j]
            con5_violation = con5_lhs - con5_rhs
            con5_gap[f"{j}"] = con5_violation 
            total_absolute_constraint_violation += abs(con5_violation)
        #print("con5_gap:", con5_gap)

        con6_gap = {}
        for j in self.nodes:
            if j not in self.source:
                #con6_rhs = self.E[j] + self.P[j]
                #con6_lhs = h_values[j]
                con6_violation = max(0, -h_values[j] + self.E[j] + self.P[j])
                con6_gap[f"{j}"] = con6_violation 
                total_absolute_constraint_violation += abs(con6_violation)
        #print("con6_gap:", con6_gap)

        print("*******************************************************************************\n")
        print("Constraints violation:\n")

        table_data = []
        for constraint, vio in con1_gap.items():
               table_data.append([constraint, f"{con1_gap[constraint]:.8f}"])
        for constraint, vio in con2_approx_gap.items():
               table_data.append([constraint, f"{con2_approx_gap[constraint]:.8f}"])
        for constraint, vio in con3_gap.items():
               table_data.append([constraint, f"{con3_gap[constraint]:.8f}"])
        for constraint, vio in con4_gap.items():
               table_data.append([constraint, f"{con4_gap[constraint]:.8f}"])
        for constraint, vio in con5_gap.items():
               table_data.append([constraint, f"{con5_gap[constraint]:.8f}"])
        for constraint, vio in con6_gap.items():
               table_data.append([constraint, f"{con6_gap[constraint]:.8f}"])

        #headers = ["Constraint ID", "Violation"]
        #print(tabulate(table_data, headers=headers, tablefmt="grid"))
        print("\nSum of constraints violation:", total_absolute_constraint_violation)

        #print("*******************************************************************************\n")
        #table_data = []
        #for constraint, vio in con2_original_gap.items():
        #       table_data.append([constraint, f"{con2_original_gap[constraint]:.8f}",  f"{con2_approx_gap[constraint]:.8f}"])

        print("*******************************************************************************\n")
        #print("Constraint 2 violations:\n")
        #headers = ["Constraint ID", "Original Con Violation", "Approx Con Violation"]
        #print(tabulate(table_data, headers=headers, tablefmt="grid"))
        
        #print("\nSum of violation of original con2:", con2_original_violation) 
        #print("Sum of violation of approx con2:", con2_approx_violation)


        table_data = []
        for constraint, vio in relative_violations.items():
            i_str, j_str = constraint.split(',')
            i, j = int(i_str), int(j_str)

            table_data.append([constraint, q_values[i,j], f"{con2_original_gap[constraint]:.8f}",  f"{con2_approx_gap[constraint]:.8f}", f"{absolute_violations[constraint]:.8f}", f"{relative_violations[constraint]:.8f}"])

        print("*******************************************************************************\n")
        print("Absolute and relative violations between original and approximation constraint 2:\n")
        headers = ["Constraint ID", "flow value", "Original Con Violation", "Approx Con Violation", "Absolute Violation", "Relative Violation"]
        print(tabulate(table_data, headers=headers, tablefmt="grid"))
        print("\nSum of violation of original headloss constraint:", con2_original_violation) 
        print("Sum of violation of approx headloss constraint:", con2_approx_violation)
        print("\nCon2 sum of absolute violation between original function and approximate function:", con2_absolute_constraint_violation)
        print("Con2 sum of relative violation between original function and approximate function:", con2_relative_constraint_violation)

        # Print total violations
        #print("\nTotal absolute constraint violation:", total_absolute_constraint_violation)
        #print("Total relative constraint violation:", total_relative_constraint_violation)

        print("*******************************************************************************\n")


    def kkt_condition_violation(self, l, q, h, eps, lam, x, y, u, v, w, tol=1e-6):
        #  Stationary Condition Violation
        violations = {}
        # dL/dl = 0
        for (i,j) in self.arcs:
            for k in self.pipes:
                lhs = self.C[k] - 10.67*q[i,j]*abs(q[i,j])**0.852 * x[i,j] / (self.R[k]**1.852 * self.d[k]**4.87) + y[i,j] - w[i,j,k]
                # lhs = self.C[k] - (10.67 * q[i,j]**3 * (q[i,j]**2 + eps[i,j]**2)**0.426 / (q[i,j]**2 + 0.426 * eps[i,j]**2)) * x[i,j] / (self.R[k]**1.852 * self.d[k]**4.87) + y[i,j] - w[i,j,k]
                if abs(lhs) > tol:
                    violations.setdefault("con1", []).append(((i,j,k), lhs))
        # dL/dq = 0
        for (i,j) in self.arcs:
            lhs = lam[j] - lam[i] - 1.852 * x[i,j] * abs(q[i,j])**0.852 * sum(10.67*l[i,j,k]/(self.R[k]**1.852 * self.d[k]**4.87) for k in self.pipes)
            if abs(lhs) > tol:
                violations.setdefault("con2", []).append(((i,j), lhs))
        # dL/dh_s = 0
        for s in self.source:
            lhs = -sum(x[s,j] for j in self.nodes if (s,j) in self.arcs) + u[s]
            if abs(lhs) > tol:
                violations.setdefault("con4", []).append((s, lhs))

        # dL/dh_j = 0
        for j in self.nodes:
            if j not in self.source:
                lhs = -sum(x[i,j] for i in self.nodes if (i,j) in self.arcs) + sum(x[j,i] for i in self.nodes if (j,i) in self.arcs) - v[j]
                if abs(lhs) > tol:
                    violations.setdefault("con3", []).append((j, lhs))
        
        # Dual Feasibility violation
        # con5
        for j in self.nodes: 
            if j not in self.source:
                if v[j] <= -tol:  # should be >= 0
                    violations.setdefault("con5", []).append((j, v[j]))
        
        # Complementary Slackness Violations
        # con6
        for j in self.nodes:
            if j not in self.source:
                lhs = (-h[j] + self.E[j] + self.P[j]) * v[j]
                if abs(lhs) > tol:
                    violations.setdefault("con6", []).append((j, lhs))

        # Primal Feasibility constraints violation
        for (i,j) in self.arcs:
            lhs = h[i] - h[j] - q[i,j] * abs(q[i,j])**0.852 * sum(10.67*l[i,j,k]/(self.R[k]**1.852 * self.d[k]**4.87) for k in self.pipes)
            if abs(lhs) > tol:
                violations.setdefault("con2", []).append(((i,j), lhs))

        return violations



    def solve_ipopt(self):
        """
        First solve with IPOPT to get a good starting point.
        """
        print(f"\n-------------------------------- Solving with IPOPT --------------------------")
        self.ampl.option['solver'] = 'ipopt' 
        # self.ampl.option["ipopt_options"] = "outlev = 0  bound_push = 0.01 bound_frac = 0.01" 
        self.ampl.option["ipopt_options"] = "outlev = 1 expect_infeasible_problem = no bound_relax_factor=0 tol = 1e-9 constr_viol_tol = 1e-9 acceptable_constr_viol_tol = 1e-9 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = no halt_on_ampl_error = yes"
        # self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 tol = 1e-10 constr_viol_tol = 1e-10 acceptable_constr_viol_tol = 1e-10 honor_original_bounds = yes acceptable_tol = 1e-10 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = no halt_on_ampl_error = yes"
        
        # self.ampl.option["ipopt_options"] = (
        #     "outlev=0 "
        #     "expect_infeasible_problem=yes "
        #     "bound_relax_factor=1e-8 "  # small relaxation to avoid eps slightly below lower bound
        #     "tol=1e-8 "
        #     "constr_viol_tol=1e-8 "
        #     "acceptable_constr_viol_tol=1e-8 "
        #     "honor_original_bounds=yes "
        #     "acceptable_tol=1e-8 "
        #     "bound_push=1e-8 "
        #     "bound_frac=1e-8 "
        #     "warm_start_init_point=no "
        #     "halt_on_ampl_error=yes"
        # )
        self.ampl.option["presolve_eps"] = "8.53e-15"

        #min_demand = self.ampl.getParameter('D_min').getValues().to_list()[0]
        #max_demand = self.ampl.getParameter('D_max').getValues().to_list()[0]
        #max_flow = self.ampl.getParameter('Q_max').getValues().to_list()[0]

        #print("min_demand:", min_demand)
        #print("max_demand:", max_demand)
        #print("max_flow:", max_flow)
        #d_max = self.ampl.getParameter('d_max').getValues().to_list()[0]
        #d_min = self.ampl.getParameter('d_min').getValues().to_list()[0]
        #max_L = max(self.L[i,j] for (i,j) in self.arcs)
        #R_min = min(self.R[k] for k in self.pipes)
        #MaxK = 10.67 / ((R_min ** 1.852) * ((d_min) ** 4.87))
        
        #epsilon = (10**(-6)/(0.07508*MaxK))**(1/0.926)
        #epsilon = (10**(-6)/(0.04001571*MaxK))**(1/1.852)
        #epsilon = self.compute_adaptive_eps(min_demand/1000)
        #epsilon = 1e-5
 
        #print("eps:", epsilon,"\n")
        
        
        #eps = ampl.getParameter('eps').to_list()
        #for (i,j) in eps.items():
            #eps[i,j].setValue(epsilon)
        #self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = {epsilon};")


        print("Ipopt solver outputs: \n")
        self.ampl.solve()
        # self.ampl.eval('write gsol "warmstart.sol";')

        total_cost = self.ampl.getObjective("total_cost").value()
        print("total_cost:", total_cost, "\n")

        #l_init = self.ampl.getVariable('l').getValues().to_dict()
        #q_init = self.ampl.getVariable('q').getValues().to_dict()
        #h_init = self.ampl.getVariable('h').getValues().to_dict()
        # eps = self.ampl.getParameter('eps').getValues().to_dict()
        #eps = self.ampl.getVariable('eps').getValues().to_dict()

        #print(l_init)
        #print(q_init)
        #print(h_init)

        #print("eps:",eps, "\n")

        print("*******************************************************************************\n")
        #print("Print the decision variables value:\n")
        #ampl.eval("display l;")
        #ampl.eval("display q;")
        # self.ampl.eval("display eps;")
        
        self.q = self.ampl.get_variable('q').get_values().to_dict()
        self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.l = self.ampl.get_variable('l').get_values().to_dict()
        self.eps = self.ampl.getParameter('eps').get_values().to_dict()
        # self.eps = self.ampl.get_variable('eps').get_values().to_dict()
        #self.ampl.eval("display eps;")
        self.ampl.eval("display q;")
        # self.ampl.eval("display h;")
        # self.ampl.eval("display q1;")
        # self.ampl.eval("display q2;")
        self.ampl.eval("display l;")
        self.ampl.eval("display eps;")
        # self.ampl.eval("display dvar;")
        # print("eps: ", self.eps[next(iter(self.eps))])
        for (i,j) in self.arcs:
            if np.abs(self.q[i,j]) <=1e-3:
                print(f"q[{i},{j}]:",self.q[i,j])
        self.constraint_violations(self.q, self.h, self.l, self.eps, self.solver_name)

        #self.constraint_violations(q_init, h_init, l_init, eps, "ipopt")
        #ampl.eval("display con1.body;")
        #ampl.eval("display con2.body;")
        #ampl.eval("display con3.body;")
        #ampl.eval("display con4.body;")
        #ampl.eval("display con5.body;")
        #ampl.eval("display con6_.body;")
        #ampl.eval("display con7.body;")
        #ampl.eval("display con8.body;")

        solve_time = self.ampl.get_value('_solve_elapsed_time')
        total_cost = self.ampl.getObjective("total_cost").value()

        print(f"Total cost using ipopt:", total_cost)
        print(f"IPOPT solve time: {solve_time:.2f} seconds")

        # Extract solutions
        self.q_init = self.ampl.get_variable('q').get_values().to_dict()
        self.h_init = self.ampl.get_variable('h').get_values().to_dict()
        self.l_init = self.ampl.get_variable('l').get_values().to_dict()
        # self.eps_init = self.ampl.get_variable('eps').get_values().to_dict()
        if self.data_number ==5:
             self.q1_init = self.ampl.get_variable('q1').get_values().to_dict()
             self.q2_init = self.ampl.get_variable('q2').get_values().to_dict()
        # Save initial points
        #for idx in q_sol.keys():
        #    self.q_init[idx] = q_sol[idx]
        #for idx in h_sol.keys():
        #    self.h_init[idx] = h_sol[idx]
        #for idx in l_sol.keys():
        #    self.l_init[idx] = l_sol[idx]
        #for idx in eps_sol.keys():
        #    self.eps_init[idx] = eps_sol[idx]

        print("*******************************************************************************\n")

    @contextlib.contextmanager
    def suppress_output(self):
        # Open devnull to suppress the output
        with open(os.devnull, 'w') as devnull:
            # Redirect stdout and stderr
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            sys.stdout = devnull
            sys.stderr = devnull
            try:
                yield
            finally:
                # Restore original stdout and stderr
                sys.stdout = old_stdout
                sys.stderr = old_stderr


    def reduced_diameter(self):
        arc_max_dia = {}
        for (i, j, d), val in self.l.items():
            if val > 1e-6:
                if (i, j) not in arc_max_dia:
                    arc_max_dia[(i, j)] = d
                else:
                    arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)
        print(arc_max_dia)
        sorted_arcs = sorted(list(self.arcs), key=lambda arc: abs(self.q[arc if arc in self.arcs else (arc[1], arc[0])]), reverse=True)
        print(sorted_arcs)
        for (i,j) in sorted_arcs:
            print("arc:",(i,j))
            ampl = AMPL()
            ampl.reset()
            ampl.read("reduced_wdnmodel.mod")
            ampl.read_data(self.data_file)
            #ampl.set['arc_max_dia'] = arc_max_dia
            new_arcs = [arc for arc in self.arcs if arc != (i, j)]
            ampl.eval(f"set new_arcs := {{{set(new_arcs)}}};")
            
            #for (x, y, k), val in self.l.items():
            #    ampl.eval(f'let l[{x},{y},{k}] := {val};')
            #for (x, y), val in self.q.items():
            #    ampl.eval(f'let q[{x},{y}] := {val};')
            #for x, val in self.h.items():
            #    ampl.eval(f'let h[{x}] := {val};')
            #self.ampl.eval(f"subject to flow_2_4: q[2,4] >= 0;")
            #self.ampl.eval(f"subject to flow_3_5: q[3,5] >= 0;")
            #self.ampl.eval(f"subject to flow_4_5: q[4,5] >= 0;")
            ampl.eval(f"subject to flow_4_6: q[4,6] >= 0;")
            #ampl.eval(f"subject to flow_6_7: q[6,7] <= 0;")
                   
            #ampl.eval(f"subject to flow_7_5: q[5,7] <= 0;")
            #ampl.eval(f"subject to flow_7_5: h[7] - h[5] <= 0;")
            ampl.eval(f"subject to con3{{(i,j) in new_arcs}}: sum{{k in pipes}} l[i,j,k] = L[i,j];")
            ampl.eval(f"subject to con3_{i}_{j}: sum{{k in pipes: k <=  {arc_max_dia[i,j]-1}}} l[{i},{j},k] = L[{i},{j}];")
            #ampl.eval(f"subject to con3_{i}_{j}: l[{i},{j},{arc_max_dia[i,j]-2}] + l[{i},{j},{arc_max_dia[i,j]-1}] = L[{i},{j}];")
            
            #ampl.eval(f"subject to con2{{(i,j) in new_arcs}}: h[i] - h[j]  = (q[i,j])^3 *((((q[i,j])^2 + eps[i,j])^0.426) /((q[i,j])^2 + 0.426*eps[i,j]))  * sum{{k in pipes}} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));")
            #ampl.eval(f"subject to con2_{i}_{j}: h[{i}] - h[{j}]  = (q[{i},{j}])^3 *((((q[{i},{j}])^2 + eps[{i},{j}])^0.426) /((q[{i},{j}])^2 + 0.426*eps[{i},{j}]))  * sum{{k in pipes: k <= {arc_max_dia[i,j]-1}}} (omega * l[{i},{j},k] / ( (R[k]^1.852) * (d[k])^4.87));")
            
            ampl.option['solver'] = "ipopt" 
            ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.001 bound_frac = 0.001 warm_start_init_point = yes halt_on_ampl_error = yes "
            with self.suppress_output():
                ampl.solve()
        
            l = ampl.getVariable('l').getValues().to_dict()
            q = ampl.getVariable('q').getValues().to_dict()
            h = ampl.getVariable('h').getValues().to_dict()
            total_cost = ampl.getObjective("total_cost").value()
            if ampl.solve_result == "solved":
                print(f"Total cost using ipopt:", total_cost)
                if total_cost < self.total_cost:
                    print(f"New optimal solution:", total_cost)
                    self.total_cost = total_cost
                    ampl.eval("display {(i,j) in arcs, k in pipes: l[i,j,k]>1e-6}: l[i,j,k];")
 
        print("best optimal solution:", self.total_cost)

    def second_solve(self,model_file, solver_name):
        """
        Second solve with user-specified solver, using IPOPT solution as a starting point.
        """
        print(f"\n-------------------------------- Solving with {self.solver_name} --------------------------")

        # self.ampl.reset()
        #self.read_model_and_data()
        self.ampl.reset()
        self.ampl.read(model_file)
        self.ampl.read_data(self.data_file)
        
        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')
        
        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
        d_min = self.ampl.getParameter('d_min').getValues().to_list()
        r_min = min(self.ampl.getParameter('R').getValues().to_list())
        #epsilon = max(ampl.getParameter('eps').getValues().to_list())
        MaxK = {}
        for (i,j) in self.arcs:
            MaxK[i,j] = 10.67*self.L[i,j]/((r_min[0]**1.852) * (d_min[0]**4.87)) 
            # print(MaxK[i,j])
        eps = {}
        # MaxK = self.ampl.getParameter('MaxK').getValues().to_dict()
        for (i,j) in self.arcs:
            val = ((10**(-1))/(0.36061*MaxK[i,j]))**(1/1.852)
            # val = ((10**(-1))/(0.07508 *MaxK[i,j]))**(1/1.852)
            order = math.floor(math.log10(abs(val)))   # find exponent
            coeff = round(val / (10 ** order), 1)
            eps[i, j] = coeff * (10 ** order)
            #eps[i,j] = 0.0048 * (10**(-2))

        # print("epsilon:", eps)
        # self.ampl.param["eps"] = eps
        # Set initial values
        #q_var = self.ampl.get_variable('q')
        #h_var = self.ampl.get_variable('h')
        #l_var = self.ampl.get_variable('l')
        #eps_var = self.ampl.get_variable('eps')

        #for idx in self.q_init:
        #   q_var[idx].set_value(self.q_init[idx])
        #for idx in self.h_init:
        #   h_var[idx].set_value(self.h_init[idx])
        #for idx in self.l_init:
        #   l_var[idx].set_value(self.l_init[idx])
        #for idx in self.eps_init:
        #   eps_var[idx].set_value(self.eps_init[idx])

        # for (i, j, k), val in self.l_init.items():
        #    self.ampl.eval(f'let l[{i},{j},{k}] := {val};')
        # for (i, j), val in self.q_init.items():
        #    self.ampl.eval(f'let q[{i},{j}] := {val};')
        #    if self.data_number ==5:
        #        self.ampl.eval(f'let q1[{i},{j}] := {self.q1_init[i,j]};')
        #        self.ampl.eval(f'let q2[{i},{j}] := {self.q2_init[i,j]};')
        # for i, val in self.h_init.items():
        #    self.ampl.eval(f'let h[{i}] := {val};')
        
        # for (i, j), val in self.eps_init.items():
        #     self.ampl.eval(f'let eps[{i},{j}] := {val};')

        #for (i, j, k), val in self.l_init.items():
        #    self.ampl.var["l"][i,j,k].fix(val)
        #for (i, j), val in self.q_init.items():
        #    self.ampl.var["q"][i,j].fix(val)
        #    if self.data_number ==5:
        #        self.ampl.var["q1"][i,j].fix(self.q1[i,j])
        #        self.ampl.var["q2"][i,j].fix(self.q2[i,j])
        #for i, val in self.h_init.items():
        #    self.ampl.var["h"][i].fix(val)

        # Change solver and solve
        self.ampl.option['solver'] = solver_name

        self.ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"

        self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = no bound_relax_factor=0 tol = 1e-9 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = no halt_on_ampl_error = yes"

        #ampl.set_option("ipopt_options", "outlev = 0 expect_infeasible_problem = yes bound_push = 0.001 bound_frac = 0.001 nlp_scaling_method = gradient-based  warm_start_init_point = yes halt_on_ampl_error = yes warm_start_bound_push=1e-9 warm_start_mult_bound_push=1e-9")   #max_iter = 1000
        self.ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 2 warm_start_init_point = no bonmin.num_resolve_at_root = 10 tol = 1e-9 expect_infeasible_problem = yes bound_relax_factor = 0 bound_push = 0.01 bound_frac = 0.01"
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 2 nodemethod = 2 concurrentmethod = 3 nonconvex = 2  warmstart = 1 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 600 warmstart = 0 method = 1  mipgapabs = 1e-6 mipgap = 1e-9 barconvtol = 1e-9 sol:chk:feastol = 1e-5 sol:chk:feastolrel = 1e-9 NumericFocus = 1 tech:optionfile = gurobiOpt.prm" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        self.ampl.option["gurobi_options"] = "outlev 1 outlev_mp = 1 presolve 1 aggregate = 1 timelimit 600 alg:numericfocus = 1 obbt = 1 pre:scale = 1 method = 2 nodemethod = 1 varbranch = 3  nonconvex = 2 mipfocus = 1 nlpheur = 1 varbranch 0  mipgapabs = 1e-5 mipgap = 1e-6 alg:feastol = 1e-5 pre:feastol = 1e-5 pre:feastolrel = 1e-6 chk:feastol = 1e-5 chk:feastolrel = 1e-6" 
        # self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 4 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 feastol = 1e-6" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 0 timelimit 3600 NumericFocus = 1" # iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 basis = 1 premiqcpform = 2 preqlin = 2"# intfeastol = 1e-5 feastol = 1e-6 chk:epsrel = 1e-6 checkinfeas chk:inttol = 1e-5 scale = 3 aggregate = 1 intfocus = 1  BarHomogeneous = 1  startnodelimit = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        
        # self.ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 version objbound wantsol = 2 iisfind = 4 threads = 8 epsr = 1e-9" # lsolver = conopt
        self.ampl.option["baron_options"]= "optfile = optfile version objbound wantsol = 2 outlev = 2 barstats" # lsolver = conopt
        #self.ampl.option["baron_options"]= "optfile = optfile" # lsolver = conopt
        # self.ampl.option["scip_options"] = (
        #                                     "outlev=1 "
        #                                     "timelimit=10 "
        #                                     "lim:gap=1e-9 "
        #                                     "chk:feastol=1e-8 "
        #                                     "chk:feastolrel=1e-8 "
        #                                     "param:read=scip.set"
        #                                      )
        # self.ampl.option["scip_options"] = "outlev 1 timelimit 20 param:read = scip.set" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        self.ampl.option["scip_options"] = "outlev 1 timelimit 100 heu:settings = 0 method = p lim:absgap=1e-5 lim:gap = 1e-6 chk:feastol = 1e-5 chk:feastolrel=1e-6 param:read = scip.set" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        # self.ampl.option["scip_options"] = "outlev  1 "
        self.ampl.option["knitro_options"] = "maxtime_real = 3600 outlev = 4 opttol_abs=1e-5 opttol = 1e-6 feastol_abs = 1.0e-5 feastol = 1.0e-5  ms_enable = 1 ms_maxsolves = 10"
        #self.ampl.option["conopt_options"]= "outlev = 4"
        self.ampl.option["presolve"] = "1"
        self.ampl.option["presolve_eps"] = "8.53e-15"
        
        # self.ampl.option["presolve_eps"] = "3.75e-14"
        #print(f"{self.solver_name} solver outputs:\n")
        
        #min_demand = self.ampl.getParameter('D_min').getValues().to_list()[0]
        #max_demand = self.ampl.getParameter('D_max').getValues().to_list()[0]
        #max_flow = self.ampl.getParameter('Q_max').getValues().to_list()[0]

        #print("min_demand:", min_demand)
        #print("max_demand:", max_demand)
        #print("max_flow:", max_flow)
        #d_max = self.ampl.getParameter('d_max').getValues().to_list()[0]
        #d_min = self.ampl.getParameter('d_min').getValues().to_list()[0]
        #max_L = max(self.L[i,j] for (i,j) in self.arcs)
        #R_min = min(self.R[k] for k in self.pipes)
        #MaxK = 10.67*max_L/((R_min**1.852) * (d_min**4.87))
 
        #epsilon = ((10**(-6))/(0.07508*MaxK))**(1/0.926)
        #epsilon = (10**(-6)/(0.04001571*MaxK))**(1/1.852)

        #eps = self.ampl.getParameter('eps').get_values().to_dict()
        #print("eps:",eps) 

        #epsilon = 1e-3
        #epsilon = self.compute_adaptive_eps(min_demand/1000)
        
        #print("eps:", epsilon,"\n")
        
        
        #eps = ampl.getParameter('eps').to_list()
        #for (i,j) in eps.items():
            #eps[i,j].setValue(epsilon)
        #self.ampl.eval(f"param MaxK{{arcs}};")
        #self.ampl.eval(f"subject to MaxK{{(i,j) in arcs}}: MaxK[i,j] = 10.67 * L[i,j]/({R_min**1.852}*d_min^4.87);")
        #self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = ((10^(-6))/(0.07508*MaxK[i,j]))^(1/0.0926);")
        #self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = 1e-6 + q[i,j]^2;")
        
        #self.ampl.eval(f"subject to flow_1_2: q[1,2] >= 0;")
        #self.ampl.eval(f"subject to flow_2_3: q[2,3] >= 0;")
        #self.ampl.eval(f"subject to flow_2_4: q[2,4] >= 0;")
        #self.ampl.eval(f"subject to flow_3_5: q[3,5] >= 0;")
        #self.ampl.eval(f"subject to flow_4_5: q[4,5] >= 0;")
        #self.ampl.eval(f"subject to flow_4_6: q[4,6] >= 0;")
        #self.ampl.eval(f"subject to flow_6_7: q[6,7] <= 0;")
        #self.ampl.eval(f"subject to flow_7_5: q[5,7] <= 0;")
        #self.ampl.eval(f"subject to flow_7_5: h[7] - h[5] <= 0;")
        #self.ampl.eval(f"subject to flow_7: abs(q[6,7]) - abs(q[5,7]) = D[7];")
        #self.ampl.eval(f"subject to flow_5: abs(q[5,7]) + abs(q[3,5]) + abs(q[4,5]) = D[5];")
        # Define MaxK as a computed parameter
        #self.ampl.eval(f"""param MaxK{{(i,j) in arcs}} := 10.67 * L[i,j] / ({R_min}^1.852 * {d_min}^4.87);""")
        
        # Now use MaxK in your constraint
        #self.ampl.eval("""subject to Epsilon_Selection{(i,j) in arcs}:eps[i,j] = (1e-6 / (0.07508 * MaxK[i,j]))^(1 / 0.0926);""")
        #self.ampl.eval("display {i in 1.._ncons} (_conname[i]);")
        #self.ampl.eval("expand ;")
        # self.ampl.eval(f"solution {warmstart.sol};")
        self.ampl.solve()

        # Check constraint violations
        self.q = self.ampl.get_variable('q').get_values().to_dict()
        self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.l = self.ampl.get_variable('l').get_values().to_dict()
        self.eps = self.ampl.getParameter('eps').get_values().to_dict()
        # self.eps = self.ampl.get_variable('eps').get_values().to_dict()
        # self.ampl.eval("display eps;")
        # self.ampl.eval("display q;")
        # self.ampl.eval("display h;")
        # self.ampl.eval("display q1;")
        # self.ampl.eval("display q2;")
        # self.ampl.eval("display eps;")
        # print("eps: ", self.eps[next(iter(self.eps))])
        # for (i,j) in self.arcs:
        #     if np.abs(self.q[i,j]) <=1e-3:
        #         print(f"q[{i},{j}]:",self.q[i,j])

        # current_duals = {}
        # v = {}
        # for con_name, val in self.ampl.get_constraints():
        #     dual_values = val.get_values()
        #     current_duals[con_name] = dual_values
        #     if con_name == "con1":
        #         lam = dual_values.to_dict()
        #         # print(lam)
        #     if con_name == "con2":
        #         x = dual_values.to_dict()
        #         print(x)
        #     if con_name == "con3":
        #         y = dual_values.to_dict()
        #         print(y)
        #     if con_name == "con5":
        #         w = dual_values.to_dict()
        #     if con_name == "con6":
        #         u = dual_values.to_dict() 
        #     if con_name == "con7":
        #         v = dual_values.to_dict()
        #         # for j in val.keys():
        #         #     v[j] = val[j]
        #         # print(v)

        # tol = 1e-6
        # self.ampl.eval("display con3;")
        # self.ampl.eval("display con3.dual;")
        
        def dualmodelsolve(file):
            ampl = AMPL()
            ampl.reset()
            # ampl.read("lagDualWdnModel.mod")
            ampl.read(file)
            ampl.read_data(self.data_file)

            for (i, j, k), value in self.l.items():
                ampl.param['l'][i, j, k] = value
            for (i, j), value in self.q.items():
                ampl.param['q'][i, j] = value
            for i, value in self.h.items():
                ampl.param['h'][i] = value
            for (i,j), value in self.eps.items():
                ampl.param['eps'][i,j] = value
            
            ampl.option['solver'] = "cplex"
            ampl.option["presolve_eps"] = "7.583e-13"
            ampl.solve()
            lam = ampl.get_variable('lam').get_values().to_dict()
            x = ampl.get_variable('x').get_values().to_dict()
            y = ampl.get_variable('y').get_values().to_dict()
            u = ampl.get_variable('u').get_values().to_dict()
            v = ampl.get_variable('v').get_values().to_dict()
            w = ampl.get_variable('w').get_values().to_dict()
            # s = ampl.get_variable('s').get_values().to_dict()

            # print("lam: ", lam)
            # print("x: ", x)
            # print("y: ", y)
            # print("u: ", u)
            # print("v: ", v)
            # print("w: ", w)
            # print("s: ", s)
            return lam, x, y, u, v, w
        # lam, x, y, u, v, w = dualmodelsolve("approximatekkt.mod")
        # dualmodelsolve("originalkkt.mod")
    
        self.constraint_violations(self.q, self.h, self.l, self.eps, self.solver_name)
        # violations = self.kkt_condition_violation(self.l, self.q, self.h, self.eps, lam, x, y, u, v, w)
        # print("violations: ", violations)
        solve_time = self.ampl.get_value('_solve_elapsed_time')
        self.total_cost = self.ampl.getObjective("total_cost").value()

        print(f"Total cost using {self.solver_name}:", self.total_cost)
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")

        # Extract solutions
        #q_sol = self.ampl.get_variable('q').get_values().to_dict()
        #h_sol = self.ampl.get_variable('h').get_values().to_dict()
        #l_sol = self.ampl.get_variable('l').get_values().to_dict()

        #self.ampl.eval("display q;")
        #self.ampl.eval("display q2;")
        #self.ampl.eval("display l;")
        #self.ampl.eval("display {(i,j) in arcs}: h[i] - h[j];")
        #self.ampl.eval("display {(i,j) in arcs, k in pipes: l[i,j,k]>1e-6}: l[i,j,k];")
        #self.ampl.eval("display {(i,j) in arcs} h[i] - h[j] - q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));")

        #self.ampl.eval("display {(i,j) in arcs} h[i] - h[j]  - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));")

        #self.ampl.eval("display {(i,j) in arcs} sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));")
 
        #self.ampl.eval("display sum {(i,j) in arcs} (q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)));")
        
        #self.ampl.eval("display sum {(i,j) in arcs} (q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)))/(q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)));")

        print("*******************************************************************************\n")


    def solve_content_model(self):
        self.ampl.reset()
        self.ampl.read("wdn_content_model.mod")
        self.ampl.read_data(self.data_file)
        
        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')
        
        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
 
        for cor, val in self.l.items():
            #ampl.eval(f"s.t. l_val{cor[0]}_{cor[1]}_{cor[2]}: l[{cor[0]},{cor[1]},{cor[2]}] = {val};")
            self.ampl.param['l'][cor[0], cor[1], cor[2]] = val


        # Set initial values
        #q_var = self.ampl.get_variable('q')
        #h_var = self.ampl.get_variable('h')
        #l_var = self.ampl.get_variable('l')

        #for idx in self.q_init:
        #    q_var[idx].set_value(self.q_init[idx])
        #for idx in self.h_init:
        #    h_var[idx].set_value(self.h_init[idx])
        #for idx in self.l_init:
        #    l_var[idx].set_value(self.l_init[idx])

        # Change solver and solve
        self.ampl.option['solver'] = self.solver_name
        self.ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"
        self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = no halt_on_ampl_error = yes "
        self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 warmstart = 1 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0 mipgap = 1e-9 NumericFocus = 1" 
        self.ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 barstats version objbound" # lsolver = conopt
        self.ampl.option["scip_options"] = "outlev  1 timelimit 3600 lim:gap = 1e-9 chk:feastol = 1e-5 chk:feastolrel=0 " #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        self.ampl.option["knitro_options"]= "maxtime_real = 3600 outlev = 4 threads=8 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10"
        #self.ampl.option["conopt_options"]= "outlev = 4"
        self.ampl.option["presolve"] = "1"
        self.ampl.option["presolve_eps"] = "8.53e-15"
        self.ampl.solve()
        #self.q = self.ampl.get_variable('q').get_values().to_dict()
        #self.h = self.ampl.get_variable('h').get_values().to_dict()
        #self.l = self.ampl.get_variable('l').get_values().to_dict()
        #self.eps = self.ampl.getParameter('eps').get_values().to_dict()
        #self.ampl.eval("display eps;")
        #self.ampl.eval("display q;")
        #self.ampl.eval("display q1;")
        #self.ampl.eval("display q2;")
        #self.ampl.eval("display eps;")

        #self.constraint_violations(self.q, self.h, self.l, self.eps, self.solver_name)

        solve_time = self.ampl.get_value('_solve_elapsed_time')
        self.total_cost = self.ampl.getObjective("total_cost").value()

        print(f"Total cost using {self.solver_name}:", self.total_cost)
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")

        self.ampl.eval("display sum{(i,j) in arcs, k in pipes} l[i,j,k]*C[k];")
        # Extract solutions
        #q_sol = self.ampl.get_variable('q').get_values().to_dict()
        #h_sol = self.ampl.get_variable('h').get_values().to_dict()
        #l_sol = self.ampl.get_variable('l').get_values().to_dict()

        self.ampl.eval("display q;")
        self.ampl.eval("display con1.dual;")
        self.ampl.eval("display E[1] - sum{k in pipes} (10.67*l[1,2,k]*(q[1,2]^1.852)/(R[k]^1.852 * d[k]^4.87)) ;")
        #self.ampl.eval("display q2;")
        #self.ampl.eval("display l;")
        #self.ampl.eval("display {(i,j) in arcs}: h[i] - h[j];")
        #self.ampl.eval("display {(i,j) in arcs, k in pipes: l[i,j,k]>1e-6}: l[i,j,k];")

        # Extract duals (node potentials)
        con1 = self.ampl.getConstraint("con1")
        duals = {}
        
        print("\nOriginal duals (up to constant):")
        for n in self.nodes:
            duals[n] = con1[n].dual()
            print(f"  {n}: {duals[n]:.6g}")
        
        source = list(self.source)[0]
        H_source = self.E[source]
        
        lambda_source = duals[source]
        print(lambda_source)
        
        shift = H_source - lambda_source
        
        print("\nShifted duals (physical heads):")
        for n in self.nodes:
            head = duals[n] + shift
            print(f"  {n}: {head:.6g}")
        
        print(f"\nAt source node '{source}': shifted dual = {duals[source] + shift} (should equal H_source = {H_source})")
        

    def run(self):
        # self.read_model_and_data()
        # First solve: IPOPT
        # self.solve_ipopt()
        if self.data_number==5:
            model_file = "newyork_model.mod"
        elif self.data_number==6:
            model_file = "blacksburg_model.mod"
        else:
            model_file = sys.argv[1]
            print(model_file)
        # Second solve: self.solver_name
        # self.second_solve(model_file, "ipopt")
        self.second_solve(self.model_file, self.solver_name)
        #self.reduced_diameter()
        #self.solve_content_model()
if __name__ == "__main__":
    model = sys.argv[1]

    data_list = [
        "d1_bessa",
        "d2_shamir",
        "d3_hanoi",
        "d4_double_hanoi",
        "d5_triple_hanoi",
        "d6_newyork",
        "d7_blacksburg",
        "d8_fossolo_iron",
        "d9_fossolo_poly_0",
        "d10_fossolo_poly_1",
        "d11_kadu",
        "d12_pescara",
        "d13_modena",
        "d14_balerma"
    ]

    # Select the data number here (0 to 18)
    data_number = int(sys.argv[3]) -1
    data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{data_list[(data_number)]}.dat"
    #data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{sys.argv[3]}"
    #data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{sys.argv[3]}"
    #data = f"/home/nitishdumoliya/waterNetwork/data/minlplib_data/{data_list1[(data_number)]}.dat"
    #print("Water Network:", f"{data_list[(data_number)]}.dat")
    print("Water Network:", f"{sys.argv[3]}")

    print("Results smooth approximation 1 of Hazen--Williams headloss constraint\n")
    # print("Results smooth approximation 2 of Hazen--Williams headloss constraint\n")

    print("*******************************************************************************")

    solver = sys.argv[2] 

    solver_instance = WaterNetworkSolver(model, solver, data, data_number)
    solver_instance.run()
