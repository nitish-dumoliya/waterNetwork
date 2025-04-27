import sys
from tabulate import tabulate
from amplpy import AMPL
ampl = AMPL()

data_list = [
    "d1_Sample_input_cycle_twoloop",
    "d2_Sample_input_cycle_hanoi",
    "d3_Sample_input_double_hanoi",
    "d4_Sample_input_triple_hanoi",
    "d5_Taichung_input",
    "d6_HG_SP_1_4",
    "d7_HG_SP_2_3",
    "d8_HG_SP_3_4",
    "d9_HG_SP_4_2",
    "d10_HG_SP_5_5",
    "d11_HG_SP_6_3",
    "d12",
    "d13",
    "d14_NewYork",
    "d15_foss_poly_0",
    "d16_foss_iron",
    "d17_foss_poly_1",
    "d18_pescara",
    "d19_modena"
]

# Select the data number here (0 to 18)
data_number = int(sys.argv[3]) -1
input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data_list[(data_number)]}.dat"
print("Water Network:", f"{data_list[(data_number)]}.dat")

#print("Results of first order approximation 1 of head loss constraint\n")
#print("Results of first order approximation 2 of head loss constraint\n")
#print("Results of smooth approximation of head loss constraint\n")

print("*******************************************************************************")


def compute_all_constraint_violations(ampl, solver):
    constraints = ampl.getConstraints()
    total_absolute_violation = 0
    total_relative_violation = 0

    absolute_violations = {}
    relative_violations = {}
    table_data = []

    for con_name in constraints:
        constraint = constraints[con_name]
        body_vals = constraint.getValues().to_dict()
        dual_vals = constraint.getDuals().to_dict()

        for index, body_value in body_vals.items():
            # Constraint name formatting
            idx_str = ",".join(str(i) for i in index) if isinstance(index, tuple) else str(index)
            full_name = f"{con_name}[{idx_str}]"

            # Constraint upper and lower bounds
            con_obj = constraint[index]
            lb = con_obj.lb()
            ub = con_obj.ub()

            # Compute violation
            violation = 0
            if lb is not None and body_value < lb:
                violation = lb - body_value
            elif ub is not None and body_value > ub:
                violation = body_value - ub

            # Relative violation (avoid zero division)
            rel_violation = abs(violation) / (abs(body_value) + 1e-10)

            if violation > 1e-8:  # You can change this tolerance
                absolute_violations[full_name] = violation
                relative_violations[full_name] = rel_violation
                table_data.append([full_name, f"{violation:.8f}", f"{rel_violation:.8f}"])

                total_absolute_violation += abs(violation)
                total_relative_violation += rel_violation

    if table_data:
        print("\nConstraint Violations:\n")
        print(tabulate(table_data, headers=["Constraint", "Absolute Violation", "Relative Violation"], tablefmt="grid"))
    else:
        print("\nNo constraint violations detected (within tolerance).\n")

    print(f"\nTotal Absolute Violation using {solver}: {total_absolute_violation:.8f}")
    print(f"Total Relative Violation using {solver}: {total_relative_violation:.8f}")


def constraint_violations(q_values, h_values, l_values, R_values, d_values, pipes, epsilon, solver):

    total_relative_constraint_violation = 0
    total_absolute_constraint_violation = 0
    relative_violations = {}
    absolute_violations = {}

    #total_relative_approx_violation = 0
    #total_absolute_approx_violation = 0
    #relative_approx_violations = {}
    #absolute_approx_violations = {}

    for (i, j) in q_values.keys():
        # Original constraint value
        original_lhs = h_values[i] - h_values[j]
        #original_rhs = q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * (0.001 ** 1.852) * \
        #               sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
        #                   for k in pipes)

        original_rhs = q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * (0.00000277971326776) 
        #original_rhs = q_values[i, j] * (abs(q_values[i, j])) ** 0.852  

        original_value = original_rhs

        # Approximated constraint value
        approx_lhs = h_values[i] - h_values[j]

        #approx_rhs1 = (0.001**1.852)*(q_values[i,j]*(abs(q_values[i,j])+148*epsilon[i,j]) /(abs(q_values[i,j])+1000*epsilon[i,j])**0.148)* \
        #             sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
        #                 for k in pipes)

        approx_rhs1 = (0.001**1.852)*(q_values[i,j]*(abs(q_values[i,j])+148*epsilon[i,j]) /(abs(q_values[i,j])+1000*epsilon[i,j])**0.148)

        #approx_rhs2 = q_values[i, j] * ((abs(q_values[i, j]) + 1000 * epsilon[i,j]) ** 0.852) * \
        #             (abs(q_values[i, j]) / (abs(q_values[i, j]) + 852 * epsilon[i,j])) * (0.001 ** 1.852) * \
        #             sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
        #                 for k in pipes)

        approx_rhs2 = q_values[i, j] * ((abs(q_values[i, j]) + 1000 * epsilon[i,j]) ** 0.852) * \
                     (abs(q_values[i, j]) / (abs(q_values[i, j]) + 852 * epsilon[i,j])) * (0.001 ** 1.852)
 
        #approx_rhs3 = ((0.001 ** 1.852)*(q_values[i, j] * (abs(q_values[i, j]) + 1000*epsilon[i,j]) ** 0.852) - \
        #            (0.002368316*epsilon[i,j] * q_values[i,j]/(abs(q_values[i,j]) + 1000*epsilon[i,j])**0.148) +\
        #            (0.175255362*(epsilon[i,j])**2) * q_values[i,j]/((abs(q_values[i,j])+1000*epsilon[i,j])**1.148)) * \
        #            sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
        #                 for k in pipes)

        approx_rhs3 = ((0.001 ** 1.852)*(q_values[i, j] * (abs(q_values[i, j]) + 1000*epsilon[i,j]) ** 0.852) - \
                    (0.002368316*epsilon[i,j] * q_values[i,j]/(abs(q_values[i,j]) + 1000*epsilon[i,j])**0.148) +\
                    (0.175255362*(epsilon[i,j])**2) * q_values[i,j]/((abs(q_values[i,j])+1000*epsilon[i,j])**1.148)) 
 
        #approx_rhs4 = (q_values[i, j] * ((abs(q_values[i, j]) + 1000 * epsilon[i,j]) ** 0.852) * \
        #             (abs(q_values[i, j]) / (abs(q_values[i, j]) + 852 * epsilon[i,j])) * (0.001 ** 1.852) +\
        #            (0.175255362*(epsilon[i,j])**2) * q_values[i,j]/((abs(q_values[i,j])+1000*epsilon[i,j])**1.148)) * \
        #            sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
        #                 for k in pipes)
        approx_rhs4 = (q_values[i, j] * ((abs(q_values[i, j]) + 1000 * epsilon[i,j]) ** 0.852) * \
                     (abs(q_values[i, j]) / (abs(q_values[i, j]) + 852 * epsilon[i,j])) * (0.001 ** 1.852) +\
                    (0.175255362*(epsilon[i,j])**2) * q_values[i,j]/((abs(q_values[i,j])+1000*epsilon[i,j])**1.148)) 

        #approx_rhs5 = (0.001**1.852)*q_values[i, j]**3 * ((q_values[i, j]**2 + 1000**2 * epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 426000*epsilon[i,j]**2)
        approx_rhs5 = (0.00000277971326776)*q_values[i, j]**3 * ((q_values[i, j]**2 + 0.0001) ** 0.426)/(q_values[i,j]**2 + 0.000426)
        #approx_rhs5 = q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)
 
        approx_value = approx_rhs5
        
        # Compute relative violation
        relative_violation = abs(original_value - approx_value) / (abs(original_value)+1e-10)
        relative_violations[f"con2_{i},{j}"] = relative_violation
        total_relative_constraint_violation += relative_violation

        # Compute absolute violation
        absolute_violation =  original_value - approx_value
        absolute_violations[f"con2_{i},{j}"] = absolute_violation
        
        total_absolute_constraint_violation += abs(absolute_violation)

        # Compute violations for the approximation function itself
        #approx_constraint_lhs = h_values[i] - h_values[j]
        #approx_constraint_rhs = approx_rhs4
        #approx_constraint_violation = approx_constraint_rhs - approx_constraint_lhs
        
        #relative_approx_violation = abs(approx_constraint_violation) / (abs(approx_constraint_rhs) + 1e-10)
        #relative_approx_violations[f"approx_con2_{i},{j}"] = relative_approx_violation
        #total_relative_approx_violation += abs(relative_approx_violation)

        #absolute_approx_violations[f"approx_con2_{i},{j}"] = approx_constraint_violation
        #total_absolute_approx_violation += abs(approx_constraint_violation)
        

    # Prepare data for tabulation
    table_data = []
    for constraint, vio in relative_violations.items():
        table_data.append([constraint, f"{absolute_violations[constraint]:.8f}", f"{relative_violations[constraint]:.8f}"])


    #print("*******************************************************************************\n")
    print("Violations between approximation constraint and original constraint:\n")
    # Print table
    headers = ["Constraint ID", "Absolute Violation", "Relative Violation"]
    print(tabulate(table_data, headers=headers, tablefmt="grid"))

    # Print total violations
    print(f"\nTotal absolute constraint violation using {solver}:", total_absolute_constraint_violation)
    print(f"Total relative constraint violation using {solver}:", total_relative_constraint_violation)


    #approx_table_data = []
    #for constraint, vio in absolute_approx_violations.items():
    #    approx_table_data.append([constraint, f"{absolute_approx_violations[constraint]:.8f}"])
    
    #print("*******************************************************************************\n")

    #print("Constraint violations for approximation function:\n")
    #print(tabulate(approx_table_data, headers=headers, tablefmt="grid"))


    print("*******************************************************************************\n")


def compute_adaptive_eps(min_demand):

    min_demand = min_demand/1000
    if min_demand < 1e-4:
        #print("min_demand1", min_demand,"\n")
        return 1e-3

    elif min_demand < 1:
        #print("min_demand2", min_demand,"\n")
        return 1e-6
    else:
        #print("min_demand3", min_demand,"\n")
        return 1e-3


ipopt_run = 0

if ipopt_run == 1:
            ampl.read(sys.argv[1])
            ampl.read_data(input_data_file)
            ampl.option["solver"] = "ipopt"
            ampl.set_option("ipopt_options", "outlev = 0  expect_infeasible_problem = yes tol = 1e-9 bound_relax_factor=0  bound_push = 0.01 bound_frac = 0.01 nlp_scaling_method = none")   #max_iter = 1000

            ampl.option["presolve_eps"] = "8.53e-15"

            
            min_demand = ampl.getParameter('D_min').getValues().to_list()[0]
            max_demand = ampl.getParameter('D_max').getValues().to_list()[0]
            max_flow = ampl.getParameter('Q_max').getValues().to_list()[0]

            print("min_demand:", min_demand/1000)
            print("max_demand:", max_demand/1000)
            print("max_flow:", max_flow/1000)
            
            epsilon = compute_adaptive_eps(min_demand)
            
            print("eps:", epsilon)
            
            
            #eps = ampl.getParameter('eps').to_list()
            #for (i,j) in eps.items():
                #eps[i,j].setValue(epsilon)
            ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = {epsilon};")


            print("Ipopt solver outputs: \n")
            ampl.solve()

            total_cost = ampl.getObjective("total_cost").value()
            print("total_cost:", total_cost, "\n")
            
            #ampl.eval(" let {(i,j) in arcs, k in pipes} l[i,j,k] := (if abs(l[i,j,k]) < 1e-6 then 0 else precision(l[i,j,k], 16));")

            l_init = ampl.getVariable('l').getValues().to_dict()
            q_init = ampl.getVariable('q').getValues().to_dict()
            h_init = ampl.getVariable('h').getValues().to_dict()
            eps = ampl.getVariable('eps').getValues().to_dict()
            #eps = ampl.getParameter('eps').getValues().to_dict()
            
            #print(l_init)
            #print(q_init)
            #print(h_init)

            #print("eps:",eps, "\n")
            nodes = ampl.getSet('nodes')
            source = ampl.getSet('Source')
            arcs = ampl.getSet('arcs')
            pipes = ampl.getSet('pipes')
            
            L = ampl.getParameter('L').to_dict()
            D = ampl.getParameter('D').to_dict()
            C = ampl.getParameter('C').to_dict()
            P = ampl.getParameter('P').to_dict()
            R = ampl.getParameter('R').to_dict()
            E = ampl.getParameter('E').to_dict()
            d = ampl.getParameter('d').to_dict()
            
            print("*******************************************************************************\n")
            #print("Print the decision variables value:\n")
            #ampl.eval("display l;")
            #ampl.eval("display q;")
            #ampl.eval("display h;")
            
            constraint_violations(q_init, h_init, l_init, R, d, pipes, eps, "ipopt")

            #ampl.eval("display con1.body;")
            #ampl.eval("display con2.body;")
            #ampl.eval("display con3.body;")
            #ampl.eval("display con4.body;")
            #ampl.eval("display con5.body;")
            #ampl.eval("display con6_.body;")
            #ampl.eval("display con7.body;")
            #ampl.eval("display con8.body;")
            
            solve_time = ampl.get_value('_solve_elapsed_time')
            total_cost = ampl.getObjective("total_cost").value()
            
            print(f"total_cost using ipopt:", total_cost)
            print(f"solve_time using ipopt:", solve_time)


            ampl.close()

#***************************************************************************************
ampl = AMPL()
ampl.read(sys.argv[1])
ampl.read_data(input_data_file)

#eps = ampl.getParameter('eps').getValues().to_list()[0]

print("*******************************************************************************\n")

if ipopt_run == 1:
            for (i, j, k), val in l_init.items():
                ampl.eval(f'let l[{i},{j},{k}] := {val};')
            for (i, j), val in q_init.items():
                ampl.eval(f'let q[{i},{j}] := {val};')
            for i, val in h_init.items():
                ampl.eval(f'let h[{i}] := {val};')

            #ampl.getVariable("q").setValues(q_init)
            #ampl.getVariable("h").setValues(h_init)
            #ampl.getVariable("l").setValues(l_init)

ampl.option["solver"]= sys.argv[2]
ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"

ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes tol = 1e-8 bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = no halt_on_ampl_error = yes "

#ampl.set_option("ipopt_options", "outlev = 0 expect_infeasible_problem = yes bound_push = 0.001 bound_frac = 0.001 nlp_scaling_method = gradient-based  warm_start_init_point = yes halt_on_ampl_error = yes warm_start_bound_push=1e-9 warm_start_mult_bound_push=1e-9")   #max_iter = 1000
ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 2 warm_start_init_point = no bonmin.num_resolve_at_root = 10 "
ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 feastol = 1e-6" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
#ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 4 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 feastol = 1e-6" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
#ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600" # iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 basis = 1 premiqcpform = 2 preqlin = 2"# intfeastol = 1e-5 feastol = 1e-6 chk:epsrel = 1e-6 checkinfeas chk:inttol = 1e-5 scale = 3 aggregate = 1 intfocus = 1  BarHomogeneous = 1  startnodelimit = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10

ampl.option["baron_options"]= "maxtime = 3600  outlev = 1 iisfind = 4 lpsolver = cplex lsolver = conopt threads = 8" # lsolver = conopt
#ampl.option["baron_options"]= "maxtime = 3600  outlev = 1 lsolver = ipopt threads = 8" # lsolver = conopt
ampl.option["scip_options"] = "outlev  1 timelimit 300 wantsol lpmethod = b" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
ampl.option["knitro_options"]= "maxtime_real = 300 outlev = 4 threads=8 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10"
#ampl.option["conopt_options"]= "outlev = 4"
ampl.option["presolve"] = "1"
ampl.option["presolve_eps"] = "8.53e-15"

print(f"{sys.argv[2]} solver outputs:\n")

min_demand = ampl.getParameter('D_min').getValues().to_list()[0]
max_demand = ampl.getParameter('D_max').getValues().to_list()[0]
Q_max = ampl.getParameter('Q_max').getValues().to_list()[0]

print("min_demand:", min_demand/1000)
print("max_demand:", max_demand/1000)
print("max_flow:", Q_max/1000)

epsilon = compute_adaptive_eps(min_demand)

print("eps:", epsilon, "\n")

ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = {epsilon};")

ampl.solve()

# ampl.eval("expand ;")

# ampl.eval("display {i in 1.._ncons: _con[i].slack < -1e-3} (_conname[i], _con[i].slack);")

l = ampl.getVariable('l').getValues().to_dict()
q = ampl.getVariable('q').getValues().to_dict()
h = ampl.getVariable('h').getValues().to_dict()
eps = ampl.getVariable('eps').getValues().to_dict()

#print("eps:",eps, "\n")

nodes = ampl.getSet('nodes')
source = ampl.getSet('Source')
arcs = ampl.getSet('arcs')
pipes = ampl.getSet('pipes')

L = ampl.getParameter('L').to_dict()
D = ampl.getParameter('D').to_dict()
C = ampl.getParameter('C').to_dict()
P = ampl.getParameter('P').to_dict()
R = ampl.getParameter('R').to_dict()
E = ampl.getParameter('E').to_dict()
d = ampl.getParameter('d').to_dict()

print("*******************************************************************************\n")
#print("Print the decision variables value:\n")
#ampl.eval("display l;")
#ampl.eval("display q;")
#ampl.eval("display h;")

#ampl.eval("display con2.body;")

#constraint_violations(q, h, l, R, d, pipes, eps, sys.argv[2])

ampl.eval("display _conname, _con.slack, _con.sstatus;")
solve_time = ampl.get_value('_solve_elapsed_time')
total_cost = ampl.getObjective("total_cost").value()

print(f"total_cost using {sys.argv[2]}:", total_cost)
print(f"solve_time using {sys.argv[2]}:", solve_time,"\n")
print("min_demand:", min_demand/1000, "\n")
print("Q_max:", Q_max/1000, "\n")
print("eps:", epsilon, "\n")
compute_all_constraint_violations(ampl, solver="IPOPT")
ampl.close()

print("*******************************************************************************\n")

