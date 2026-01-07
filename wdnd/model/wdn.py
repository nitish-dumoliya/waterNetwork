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
        self.ampl.setOption('hsllib', '/usr/local/lib/libma57.so')        
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

        n_vars = self.ampl.getValue("_nvars")
        n_cons = self.ampl.getValue("_ncons")
        # n_nl_cons = self.ampl.getValue("_snl")
        # n_nl_obj  = self.ampl.getValue("_nlobj")

        print("Variables:", n_vars)
        print("Constraints:", n_cons)
        # print("Nonlinear constraints:", n_nl_cons)
        # print("Nonlinear objective:", n_nl_obj)


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
                approx_rhs = (q1[i, j]**3 * ((q1[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q1[i,j]**2 + 0.426*epsilon[i,j]**2)) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j]**3 * ((q2[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q2[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                # approx_rhs = (q1[i, j] * ((q1[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j] * ((q2[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

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
                    approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                    # approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

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
                approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * 10.67 * self.L[i,j]/(self.fix_r[i,j]**1.852 * self.exdiam[i,j]**4.87) 
                # approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * 10.67 * self.L[i,j]/(self.fix_r[i,j]**1.852 * self.exdiam[i,j]**4.87) 

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
                approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                # approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs
                # approx_rhs = self.a[i,j]*q_values[i, j]*abs(q_values[i,j]) + self.b[i,j]*q_values[i,j]

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

    def second_solve(self,model_file, solver_name):
        print(f"\n-------------------------------- Solving with {self.solver_name} --------------------------")
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
        self.ampl.option['solver'] = solver_name
        self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = no bound_relax_factor=0 tol = 1e-9 bound_push = 0.1 bound_frac = 0.1 warm_start_init_point = no halt_on_ampl_error = yes "
        self.ampl.option["presolve"] = "1"
        self.ampl.option["presolve_eps"] = "8.53e-15" 
        self.ampl.solve()
        self.q = self.ampl.get_variable('q').get_values().to_dict()
        if self.data_number == 5:
            self.q1 = self.ampl.get_variable('q1').get_values().to_dict()
            self.q2 = self.ampl.get_variable('q2').get_values().to_dict()
        self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.l = self.ampl.get_variable('l').get_values().to_dict()
        self.eps = self.ampl.getParameter('eps').get_values().to_dict()
        solve_time = self.ampl.get_value('_solve_elapsed_time')
        self.total_cost = self.ampl.getObjective("total_cost").value()
        print(f"Total cost using {self.solver_name}:", self.total_cost)
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")
        print("*******************************************************************************\n") 

    def third_solve(self,model_file, solver_name):
        print(f"\n-------------------------------- Solving with {self.solver_name} --------------------------")
        # arc_max_dia = {}
        # if self.data_number == 6:
        #     self.fixarcs = self.ampl.getSet('fixarcs')
        #     #print("fixarcs:",self.fixarcs)
        #     for (i, j, d), val in self.l.items():
        #         if (i,j) not in self.fixarcs or (j,i) not in self.fixarcs:
        #             if val > 1e-3:
        #                 if (i, j) not in arc_max_dia:
        #                     arc_max_dia[(i, j)] = d
        #                 else:
        #                     arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)
        # else:
        #     for (i, j, d), val in self.l.items():
        #         if val > 1e-3:
        #             if (i, j) not in arc_max_dia:
        #                 arc_max_dia[(i, j)] = d
        #             else:
        #                 arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)
        ampl = AMPL()
        ampl.reset()
        ampl.read(model_file)
        ampl.read_data(self.data_file) 
        # for (x, y, k), val in self.l.items():
        #     ampl.eval(f'let l[{x},{y},{k}] := {val};')
        # for (x, y), val in self.q.items():
        #     ampl.eval(f'let q[{x},{y}] := {val};')
        #     if self.data_number ==5:
        #         ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
        #         ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
        # for x, val in self.h.items():
        #     ampl.eval(f'let h[{x}] := {val};') 
        edge = min(self.con2_violation, key=self.con2_violation.get) 
        print(edge, self.con2_violation[edge[0], edge[1]])
        # current_duals = {}
        # for con_name, val in self.ampl.get_constraints():
        #     dual_values = val.get_values()
        #     current_duals[con_name] = dual_values
        #
        # # Initialize dual values for all constraints
        # for con_name, dual_values in self.all_duals.items():
        #     if con_name in current_duals:
        #         # Initialize dual values for each constraint
        #         self.ampl.get_constraint(con_name).set_values(dual_values)               
        for k in self.pipes:
            if k>=self.arc_max_dia[edge[0], edge[1]]:
                ampl.eval(f"subject to con3__{edge[0]}_{edge[1]}_{k}: l[{edge[0]},{edge[1]},{k}] = 0;")
            elif k == self.arc_max_dia[edge[0], edge[1]]-1:
                ampl.eval(f"subject to con3__{edge[0]}_{edge[1]}_{k}: l[{edge[0]},{edge[1]},{k}] = L[{edge[0]}, {edge[1]}];")

        ampl.option['solver'] = solver_name
        ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"
        ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = no bound_relax_factor=0 tol = 1e-9 bound_push = 0.1 bound_frac = 0.1 warm_start_init_point = yes halt_on_ampl_error = yes "
        ampl.option["presolve"] = "1"
        ampl.option["presolve_eps"] = "8.53e-15" 
        ampl.solve()
        # Check constraint violations
        q = ampl.get_variable('q').get_values().to_dict()
        h = ampl.get_variable('h').get_values().to_dict()
        l = ampl.get_variable('l').get_values().to_dict()
        eps = ampl.getParameter('eps').get_values().to_dict()
        solve_time = ampl.get_value('_solve_elapsed_time')
        total_cost = ampl.getObjective("total_cost").value()
        print(f"Total cost using {self.solver_name}:", total_cost)
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")
        print("*******************************************************************************\n")

    def solve_reducred_model(self):
        self.arc_max_dia = {}
        if self.data_number == 6:
            self.fixarcs = self.ampl.getSet('fixarcs')
            #print("fixarcs:",self.fixarcs)
            for (i, j, d), val in self.l.items():
                if (i,j) not in self.fixarcs or (j,i) not in self.fixarcs:
                    if val > 1e-3:
                        if (i, j) not in self.arc_max_dia:
                            self.arc_max_dia[(i, j)] = d
                        else:
                            self.arc_max_dia[(i, j)] = max(self.arc_max_dia[(i, j)], d)
        else:
            for (i, j, d), val in self.l.items():
                if val > 1e-3:
                    if (i, j) not in self.arc_max_dia:
                        self.arc_max_dia[(i, j)] = d
                    else:
                        self.arc_max_dia[(i, j)] = max(self.arc_max_dia[(i, j)], d)
        # print("\n*********************************************************************************************")

        self.all_duals = {}
        for con_name, val in self.ampl.get_constraints():
            self.all_duals[con_name] = val.getValues()
        sorted_arcs = []
        dual_dict = self.all_duals["con2"].to_dict()
        print("dual_value:", dual_dict)
        sorted_duals = dict(sorted(dual_dict.items(), key=lambda kv: abs(kv[1]), reverse=True))
        sorted_arcs = list(sorted_duals.keys())
        # print("\nsorted_arcs:", sorted_arcs)
        # if self.data_number == 6:
        #     sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fixarcs]
        # sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fix_arc_set]
        sorted_arcs = [arc for arc in sorted_arcs if self.arc_max_dia[arc[0], arc[1]] != 1]
        print("sorted_arcs:", sorted_arcs)
        self.cost = {}
        for (i,j) in sorted_arcs:
            print("Arc:", (i,j))
            ampl = AMPL()
            ampl.reset()
            ampl.read("reduced_wdnmodel.mod")
            # ampl.read("wdn_content_model.mod")
            ampl.read_data(self.data_file) 
            for cor, val in self.l.items():
                if (cor[0],cor[1]) != (i,j): 
                    #ampl.eval(f"s.t. l_val{cor[0]}_{cor[1]}_{cor[2]}: l[{cor[0]},{cor[1]},{cor[2]}] = {val};")
                    ampl.param['l'][cor[0], cor[1], cor[2]] = val
                # else:
                # self.ampl.param['l'][cor[0], cor[1], arc_max_dia[i,j] - 1] = self.L[i,j]
            for k in self.pipes:
                if k>=self.arc_max_dia[i,j]:
                    ampl.param['l'][i, j, k] = 0
                elif k == self.arc_max_dia[i,j]-1:
                    ampl.param['l'][i, j, k] = self.L[i,j]
                else:
                    ampl.param['l'][i, j, k] = 0 
            # Change solver and solve
            ampl.option['solver'] = self.solver_name
            # self.ampl.option["hsllib"]= "/usr/local/lib/libma57.so"        
            ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.1 bound_frac = 0.1 warm_start_init_point = no halt_on_ampl_error = yes"
            ampl.option["presolve"] = "1"
            ampl.option["presolve_eps"] = "2.73e-13" 
            with self.suppress_output():
                ampl.solve() 
            q = ampl.get_variable('q').get_values().to_dict()
            h = ampl.get_variable('h').get_values().to_dict()
            # solve_time = ampl.get_value('_solve_elapsed_time')
            # self.violation = ampl.getObjective("violation").value()
            # self.con2_violation[i,j] = self.violation
            # print(f"Total cost using {self.solver_name}:", self.violation)
            # ampl.eval("display {(i,j) in arcs} h[i] - h[j] - q[i,j] * abs(q[i,j])^0.852 * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));")
            # ampl.eval("display {i in nodes} E[i] + P[i]; display {i in nodes} h[i];")
            # ampl.eval("display sum {i in nodes} abs(h[i] - E[i] - P[i]);")
            # print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")
            h_vio = 0 
            for i in self.nodes:
                if h[i]- self.E[i] - self.P[i]>=0:
                    # pass
                    h_vio = h_vio + 0
                else:
                    h_vio = h_vio + np.abs(h[i]- self.E[i] - self.P[i])
            print("h_vio: ",h_vio)
            # ampl = AMPL()
            # ampl.reset()
            # # ampl.read("reduced_wdnmodel.mod")
            # ampl.read("lp_model2.mod")
            # ampl.read_data(self.data_file) 
            # for cor, val in q.items():
            #     ampl.param['q'][cor[0], cor[1]] = val
            # # Change solver and solve
            # ampl.option['solver'] = "cplex"
            # # self.ampl.option["hsllib"]= "/usr/local/lib/libma57.so"        
            # ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.1 bound_frac = 0.1 warm_start_init_point = no halt_on_ampl_error = yes"
            # ampl.option["presolve"] = "1"
            # ampl.option["presolve_eps"] = "2.73e-13" 
            # with self.suppress_output():
            #     ampl.solve() 
            self.solve_result = ampl.solve_result
            print("solver_result:", self.solve_result)
            solve_time = ampl.get_value('_solve_elapsed_time')
            self.total_cost = ampl.getObjective("total_cost").value()
            self.cost[i,j] = self.total_cost
            print(f"Total cost:", self.total_cost)
            # ampl.eval("display {(i,j) in arcs} h[i] - h[j] - q[i,j] * abs(q[i,j])^0.852 * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));")
            print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")

        edge = min(self.cost, key=self.cost.get) 
        print("minimum cost:", self.cost[edge[0], edge[1]])

    # Function to suppress output
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

    def run(self):
        # First solve: IPOPT
        # self.solve_ipopt()
        if self.data_number==5:
            model_file = "newyork_model.mod"
        elif self.data_number==6:
            model_file = "blacksburg_model.mod"
        else:
            model_file = "wdnmodel.mod"
            print(model_file)
        # self.read_model_and_data()
        # Second solve: self.solver_name
        # self.second_solve(model_file, "ipopt")
        self.second_solve(self.model_file, self.solver_name)
        #self.reduced_diameter()
        self.solve_reducred_model()
        # self.third_solve(self.model_file, self.solver_name)

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
    print("Results smooth approximation 2 of Hazen--Williams headloss constraint using relative error\n")
    print("*******************************************************************************")
    solver = sys.argv[2] 
    solver_instance = WaterNetworkSolver(model, solver, data, data_number)
    solver_instance.run()
