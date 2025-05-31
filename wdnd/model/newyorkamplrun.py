import time
import sys
from tabulate import tabulate
from amplpy import AMPL

class WaterNetworkSolver:
    def __init__(self, model_file, solver_name, data_file):
        self.model_file = model_file
        self.solver_name = solver_name
        self.data_file = data_file
        self.ampl = AMPL()
        
        # To store solutions
        self.q_init = {}
        self.h_init = {}
        self.l_init = {}
        self.eps = {}
        
    def read_model_and_data(self):
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
        self.exdiam = self.ampl.getParameter('exdiam').to_dict()
 

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

    def constraint_violations(self, q_values, q1, q2, h_values, l_values, epsilon, solver):
        total_absolute_constraint_violation = 0
        total_relative_constraint_violation = 0
         
        con1_gap = {}
        for i in self.nodes:
            if i not in self.source:
                con1_rhs = self.D[i]
                incoming_flow = sum(q_values[j, i] for j in self.nodes if (j, i) in self.arcs)
                outgoing_flow = sum(q_values[i, j] for j in self.nodes if (i, j) in self.arcs)
                con1_lhs = incoming_flow - outgoing_flow
                con1_violation = con1_lhs - con1_rhs
                con1_gap[f"con1_{i}"] = con1_violation
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
        for (i, j) in q_values.keys():
            # Original constraint value
            lhs = 2*(h_values[i] - h_values[j])
            alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[i,j] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
            #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
            original_rhs = q1[i, j] * (abs(q1[i, j])) ** 0.852 * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87)+ q2[i, j] * (abs(q2[i, j])) ** 0.852 * alpha_rhs  
            #original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs
            
            # Approximated constraint value
            approx_rhs = (q1[i, j]**3 * ((q1[i, j]**2 + epsilon[i,j]) ** 0.426)/(q1[i,j]**2 + 0.426*epsilon[i,j])) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j]**3 * ((q2[i, j]**2 + epsilon[i,j]) ** 0.426)/(q2[i,j]**2 + 0.426*epsilon[i,j])) * alpha_rhs

            #approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs

            con2_original_violation =  lhs - original_rhs
            con2_original_gap[f"con2_{i}_{j}"] = con2_original_violation
            con2_original_violation += abs(con2_original_violation) 
            
            con2_approx_violation =  lhs - approx_rhs
            con2_approx_gap[f"con2_{i}_{j}"] = con2_approx_violation
            
            total_absolute_constraint_violation += abs(con2_approx_violation)    
            con2_approx_violation += abs(con2_approx_violation) 

             # Compute absolute violation
            absolute_violation =  original_rhs - approx_rhs
            absolute_violations[f"con2_{i},{j}"] = absolute_violation
            con2_absolute_constraint_violation += abs(absolute_violation)

            # Compute relative violation between original_rhs and approx_rhs
            relative_violation = (original_rhs - approx_rhs) / (original_rhs+1e-14)
            relative_violations[f"con2_{i},{j}"] = relative_violation
            con2_relative_constraint_violation += abs(relative_violation)
           
        #print("con2_gap:", con2_gap)
        
        con3_gap = {}
        for (i,j) in self.arcs:
            con3_rhs = self.L[i,j]
            con3_lhs = sum(l_values[i,j,k] for k in self.pipes) 
            con3_violation = con3_lhs - con3_rhs
            con3_gap[f"con3_{i}"] = con3_violation 
            total_absolute_constraint_violation += abs(con3_violation)
        #print("con3_gap:", con3_gap)

        con4_gap = {}
        for (i,j) in self.arcs:
            for k in self.pipes:
                #con4_rhs = self.L[i,j]
                #con4_lhs = l_values[i,j,k]
                con4_violation = max(0,l_values[i,j,k]-self.L[i,j])
                con4_gap[f"con4_{i}_{j}_{k}"] = con4_violation 
                total_absolute_constraint_violation += abs(con4_violation)
        #print("con4_gap:", con4_gap)
        
        con5_gap = {}
        for j in self.source:
            con5_rhs = self.E[j]
            con5_lhs = h_values[j]
            con5_violation = con5_lhs - con5_rhs
            con5_gap[f"con5_{j}"] = con5_violation 
            total_absolute_constraint_violation += abs(con5_violation)
        #print("con5_gap:", con5_gap)

        con6_gap = {}
        for j in self.nodes:
            if j not in self.source:
                #con6_rhs = self.E[j] + self.P[j]
                #con6_lhs = h_values[j]
                con6_violation = max(0, -h_values[j] + self.E[j] + self.P[j])
                con6_gap[f"con6_{j}"] = con6_violation 
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

        print("*******************************************************************************\n")
        table_data = []
        for constraint, vio in con2_original_gap.items():
               table_data.append([constraint, f"{con2_original_gap[constraint]:.8f}",  f"{con2_approx_gap[constraint]:.8f}"])

        print("*******************************************************************************\n")
        print("Constraint 2 violations:\n")
        #headers = ["Constraint ID", "Original Con Violation", "Approx Con Violation"]
        #print(tabulate(table_data, headers=headers, tablefmt="grid"))
        
        print("\nSum of violation of original con2:", con2_original_violation)
        print("Sum of violation of approx con2:", con2_approx_violation)


        table_data = []
        for constraint, vio in relative_violations.items():
               table_data.append([constraint, f"{absolute_violations[constraint]:.8f}", f"{relative_violations[constraint]:.8f}"])

        print("*******************************************************************************\n")
        print("Absolute and relative violations between original and approximation constraint 2:\n")
        #headers = ["Constraint ID", "Absolute Violation", "Relative Violation"]
        #print(tabulate(table_data, headers=headers, tablefmt="grid"))
        print("\nCon2 sum of absolute violation:", con2_absolute_constraint_violation)
        print("Con2 sum of relative violation:", con2_relative_constraint_violation)

        # Print total violations
        #print("\nTotal absolute constraint violation:", total_absolute_constraint_violation)
        #print("Total relative constraint violation:", total_relative_constraint_violation)

        print("*******************************************************************************\n")



    def solve_ipopt(self):
        """
        First solve with IPOPT to get a good starting point.
        """
        print(f"\n-------------------------------- Solving with IPOPT --------------------------")
        self.ampl.option['solver'] = 'ipopt' 
        self.ampl.option["ipopt_options"] = "outlev = 0  expect_infeasible_problem = yes tol = 1e-9 bound_relax_factor=0  bound_push = 0.01 bound_frac = 0.01 nlp_scaling_method = none" 
        self.ampl.option["presolve_eps"] = "8.53e-15"

        min_demand = self.ampl.getParameter('D_min').getValues().to_list()[0]
        max_demand = self.ampl.getParameter('D_max').getValues().to_list()[0]
        max_flow = self.ampl.getParameter('Q_max').getValues().to_list()[0]

        print("min_demand:", min_demand)
        print("max_demand:", max_demand)
        print("max_flow:", max_flow)
        d_max = self.ampl.getParameter('d_max').getValues().to_list()[0]
        d_min = self.ampl.getParameter('d_min').getValues().to_list()[0]
        #max_L = max(self.L[i,j] for (i,j) in self.arcs)
        R_min = min(self.R[k] for k in self.pipes)
        MaxK = 10.67 / ((R_min ** 1.852) * ((d_min) ** 4.87))
        
        epsilon = (10**(-6)/(0.07508*MaxK))**(1/0.926)
        #epsilon = (10**(-6)/(0.04001571*MaxK))**(1/1.852)
        #epsilon = self.compute_adaptive_eps(min_demand/1000)
        #epsilon = 1e-5
 
        print("eps:", epsilon,"\n")
        
        
        #eps = ampl.getParameter('eps').to_list()
        #for (i,j) in eps.items():
            #eps[i,j].setValue(epsilon)
        self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = {epsilon};")


        print("Ipopt solver outputs: \n")
        self.ampl.solve()

        total_cost = self.ampl.getObjective("total_cost").value()
        print("total_cost:", total_cost, "\n")

        l_init = self.ampl.getVariable('l').getValues().to_dict()
        q_init = self.ampl.getVariable('q').getValues().to_dict()
        h_init = self.ampl.getVariable('h').getValues().to_dict()
        eps = self.ampl.getVariable('eps').getValues().to_dict()
        #eps = self.ampl.getParameter('eps').getValues().to_dict()

        #print(l_init)
        #print(q_init)
        #print(h_init)

        #print("eps:",eps, "\n")

        print("*******************************************************************************\n")
        #print("Print the decision variables value:\n")
        #ampl.eval("display l;")
        #ampl.eval("display q;")
        #ampl.eval("display h;")

        self.constraint_violations(q_init, h_init, l_init, eps, "ipopt")
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
        q_sol = self.ampl.get_variable('q').get_values().to_dict()
        h_sol = self.ampl.get_variable('h').get_values().to_dict()
        l_sol = self.ampl.get_variable('l').get_values().to_dict()

        # Save initial points
        for idx in q_sol.keys():
            self.q_init[idx] = q_sol[idx]
        for idx in h_sol.keys():
            self.h_init[idx] = h_sol[idx]
        for idx in l_sol.keys():
            self.l_init[idx] = l_sol[idx]
        print("*******************************************************************************\n")


    def second_solve(self):
        """
        Second solve with user-specified solver, using IPOPT solution as a starting point.
        """
        print(f"\n-------------------------------- Solving with {self.solver_name} --------------------------")

        self.ampl.reset()
        self.read_model_and_data()

        # Set initial values
        q_var = self.ampl.get_variable('q')
        h_var = self.ampl.get_variable('h')
        l_var = self.ampl.get_variable('l')

        for idx in self.q_init:
            q_var[idx].set_value(self.q_init[idx])
        for idx in self.h_init:
            h_var[idx].set_value(self.h_init[idx])
        for idx in self.l_init:
            l_var[idx].set_value(self.l_init[idx])

        # Change solver and solve
        self.ampl.option['solver'] = self.solver_name

        self.ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"
        
        self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = no halt_on_ampl_error = yes "
        
        #ampl.set_option("ipopt_options", "outlev = 0 expect_infeasible_problem = yes bound_push = 0.001 bound_frac = 0.001 nlp_scaling_method = gradient-based  warm_start_init_point = yes halt_on_ampl_error = yes warm_start_bound_push=1e-9 warm_start_mult_bound_push=1e-9")   #max_iter = 1000
        #self.ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 2 warm_start_init_point = no bonmin.num_resolve_at_root = 10 "
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 2 nodemethod = 2 concurrentmethod = 3 nonconvex = 2  warmstart = 1 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 method = 2 warmstart = 1 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0 mipgap = 1e-9" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 " 
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 4 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 feastol = 1e-6" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600" # iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 basis = 1 premiqcpform = 2 preqlin = 2"# intfeastol = 1e-5 feastol = 1e-6 chk:epsrel = 1e-6 checkinfeas chk:inttol = 1e-5 scale = 3 aggregate = 1 intfocus = 1  BarHomogeneous = 1  startnodelimit = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        
        #self.ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 version objbound wantsol = 2 iisfind = 4 threads = 8 epsr = 1e-9" # lsolver = conopt
        self.ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 barstats version objbound" # lsolver = conopt
        #self.ampl.option["baron_options"]= "optfile = optfile" # lsolver = conopt
        self.ampl.option["scip_options"] = "outlev  1 timelimit 300 wantsol lpmethod = b" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        self.ampl.option["knitro_options"]= "maxtime_real = 300 outlev = 4 threads=8 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10"
        #self.ampl.option["conopt_options"]= "outlev = 4"
        self.ampl.option["presolve"] = "1"
        self.ampl.option["presolve_eps"] = "8.53e-15"
        
        #print(f"{self.solver_name} solver outputs:\n")
        
        min_demand = self.ampl.getParameter('D_min').getValues().to_list()[0]
        max_demand = self.ampl.getParameter('D_max').getValues().to_list()[0]
        max_flow = self.ampl.getParameter('Q_max').getValues().to_list()[0]

        print("min_demand:", min_demand)
        print("max_demand:", max_demand)
        print("max_flow:", max_flow)
        d_max = self.ampl.getParameter('d_max').getValues().to_list()[0]
        d_min = self.ampl.getParameter('d_min').getValues().to_list()[0]
        #max_L = max(self.L[i,j] for (i,j) in self.arcs)
        R_min = min(self.R[i,j] for (i,j) in self.arcs)
        MaxK = 10.67/((R_min**1.852) * (d_min**4.87))
 
        epsilon = ((10**(-6))/(0.07508*MaxK))**(1/0.926)
        #epsilon = (10**(-6)/(0.04001571*MaxK))**(1/1.852) 

        #epsilon = 1e-3
        #epsilon = self.compute_adaptive_eps(min_demand/1000)
        
        #print("eps:", epsilon,"\n")
        
        
        #eps = ampl.getParameter('eps').to_list()
        #for (i,j) in eps.items():
            #eps[i,j].setValue(epsilon)
        self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = {epsilon};")
        #self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = 1e-6 + q[i,j]^2;")
        
        #self.ampl.eval("display {i in 1.._ncons} (_conname[i]);")
        #self.ampl.eval("expand ;")
        self.ampl.solve()

        # Check constraint violations
        q = self.ampl.get_variable('q').get_values().to_dict()
        q1 = self.ampl.get_variable('q1').get_values().to_dict()
        q2 = self.ampl.get_variable('q2').get_values().to_dict()
        h = self.ampl.get_variable('h').get_values().to_dict()
        l = self.ampl.get_variable('l').get_values().to_dict()
        eps = self.ampl.get_variable('eps').get_values().to_dict()
        #self.ampl.eval("display eps;")
        #self.ampl.eval("display q;")
        #self.ampl.eval("display q1;")
        #self.ampl.eval("display q2;")
        #self.ampl.eval("display z;")

        self.constraint_violations(q, q1, q2, h, l, eps, self.solver_name)

        solve_time = self.ampl.get_value('_solve_elapsed_time')
        total_cost = self.ampl.getObjective("total_cost").value()

        print(f"Total cost using {self.solver_name}:", total_cost)
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")

        # Extract solutions
        q_sol = self.ampl.get_variable('q').get_values().to_dict()
        h_sol = self.ampl.get_variable('h').get_values().to_dict()
        l_sol = self.ampl.get_variable('l').get_values().to_dict()

        #self.ampl.eval("display l;")
        self.ampl.eval("display q1;")
        self.ampl.eval("display q2;")
        #self.ampl.eval("display h;")

        #self.ampl.eval("display {(i,j) in arcs} h[i] - h[j] - q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));")

        #self.ampl.eval("display {(i,j) in arcs} h[i] - h[j]  - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));")

        #self.ampl.eval("display {(i,j) in arcs} sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));")
 
        #self.ampl.eval("display sum {(i,j) in arcs} (q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)));")
        
        #self.ampl.eval("display sum {(i,j) in arcs} (q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)))/(q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)));")

        
        print("*******************************************************************************\n")

    def run(self):
        self.read_model_and_data()

        # First solve: IPOPT
        #self.solve_ipopt()

        # Second solve: self.solver_name
        self.second_solve()

if __name__ == "__main__":
    model = sys.argv[1]

    data_list = [
        "d1_shamir",
        "d2_hanoi",
        "d3_fossolo_iron",
        "d4_fossolo_poly_0",
        "d5_fossolo_poly_1",
        "d6_modena",
        "d7_pescara"
    ]

    # Select the data number here (0 to 18)
    #data_number = int(sys.argv[3]) -1
    #data = f"/home/nitishdumoliya/waterNetwork/data/{data_list[(data_number)]}.dat"
    data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{sys.argv[3]}"
    #data = f"/home/nitishdumoliya/waterNetwork/data/minlplib_data/{data_list1[(data_number)]}.dat"
    #print("Water Network:", f"{data_list[(data_number)]}.dat")
    print("Water Network:", f"{sys.argv[3]}")

    #print("Results of first order approximation 1 of head loss constraint\n")
    #print("Results of first order approximation 2 of head loss constraint\n")
    #print("Results of smooth approximation of head loss constraint\n")

    print("*******************************************************************************")

    solver = sys.argv[2] 

    solver_instance = WaterNetworkSolver(model, solver, data)
    solver_instance.run()
