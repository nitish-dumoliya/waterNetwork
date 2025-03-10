import sys
from tabulate import tabulate
from amplpy import AMPL
ampl = AMPL()

print("Water Network:", sys.argv[3],"\n")

#print("Results of first order approximation 1 of head loss constraint\n")
#print("Results of first order approximation 2 of head loss constraint\n")
print("Results of second order approximation of head loss constraint\n")

print("*******************************************************************************")

def constraint_violations(q_values, h_values, l_values, R_values, d_values, pipes, epsilon, solver):
    
    total_relative_constraint_violation = 0
    total_absolute_constraint_violation = 0
    relative_violations = {}
    absolute_violations = {}

    for (i, j) in q_values.keys():
        # Original constraint value
        original_lhs = h_values[i] - h_values[j]
        original_rhs = q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * (0.001 ** 1.852) * \
                       sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
                           for k in pipes)
        original_value = original_rhs

        # Approximated constraint value
        approx_lhs = h_values[i] - h_values[j]

        approx_rhs1 = (0.001**1.852)*(q_values[i,j]*(abs(q_values[i,j])+148*epsilon[i,j]) /(abs(q_values[i,j])+1000*epsilon[i,j])**0.148)* \
                     sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
                         for k in pipes)


        approx_rhs2 = q_values[i, j] * ((abs(q_values[i, j]) + 1000 * epsilon[i,j]) ** 0.852) * \
                     (abs(q_values[i, j]) / (abs(q_values[i, j]) + 852 * epsilon[i,j])) * (0.001 ** 1.852) * \
                     sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
                         for k in pipes)

        approx_rhs3 = ((0.001 ** 1.852)*(q_values[i, j] * (abs(q_values[i, j]) + 1000*epsilon[i,j]) ** 0.852) - \
                    (0.002368316*epsilon[i,j] * q_values[i,j]/(abs(q_values[i,j]) + 1000*epsilon[i,j])**0.148) +\
                    (0.175255362*(epsilon[i,j])**2) * q_values[i,j]/((abs(q_values[i,j])+1000*epsilon[i,j])**1.148)) * \
                    sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
                         for k in pipes)
        approx_rhs4 = (q_values[i, j] * ((abs(q_values[i, j]) + 1000 * epsilon[i,j]) ** 0.852) * \
                     (abs(q_values[i, j]) / (abs(q_values[i, j]) + 852 * epsilon[i,j])) * (0.001 ** 1.852) +\
                    (0.175255362*(epsilon[i,j])**2) * q_values[i,j]/((abs(q_values[i,j])+1000*epsilon[i,j])**1.148)) * \
                    sum(10.67 * l_values[i, j, k] / ((R_values[k] ** 1.852) * ((d_values[k] / 1000) ** 4.87))
                         for k in pipes)




        approx_value = approx_rhs4
        
        # Compute relative violation
        relative_violation = abs(original_value - approx_value) / (abs(original_value)+1e-10)
        relative_violations[f"con2_{i},{j}"] = relative_violation
        total_relative_constraint_violation += abs(relative_violation)

        # Compute absolute violation
        absolute_violation =  original_value - approx_value
        absolute_violations[f"con2_{i},{j}"] = absolute_violation
        
        total_absolute_constraint_violation += abs(absolute_violation)

    # Prepare data for tabulation
    table_data = []
    for constraint, vio in relative_violations.items():
        table_data.append([constraint, f"{absolute_violations[constraint]:.8f}", f"{relative_violations[constraint]:.8f}"])
    
    #print("*******************************************************************************\n")
    print("Constraint violations:\n")
    # Print table
    headers = ["Constraint ID", "Absolute Violation", "Relative Violation"]
    print(tabulate(table_data, headers=headers, tablefmt="grid"))

    # Print total violations
    print(f"\nTotal absolute constraint violation using {solver}:", total_absolute_constraint_violation)
    print(f"Total relative constraint violation using {solver}:", total_relative_constraint_violation)

    print("*******************************************************************************\n")




ipopt_run = 1

if ipopt_run == 1:
            ampl.read(sys.argv[1])
            ampl.read_data(sys.argv[3])
            ampl.option["solver"] = "ipopt"

            ampl.set_option("ipopt_options", "outlev = 0  expect_infeasible_problem = yes bound_push = 0.01 bound_frac = 0.01 nlp_scaling_method = gradient-based ")   #max_iter = 1000

            ampl.option["presolve_eps"] = "8.53e-15"
            
            print("Ipopt solver outputs: \n")
            ampl.solve()

            total_cost = ampl.getObjective("total_cost").value()
            print("total_cost:", total_cost, "\n")

            l_init = ampl.getVariable('l').getValues().to_dict()
            q_init = ampl.getVariable('q').getValues().to_dict()
            h_init = ampl.getVariable('h').getValues().to_dict()
            eps = ampl.getVariable('eps').getValues().to_dict()

            print("eps:",eps, "\n")
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
            
            solve_time = ampl.get_value('_solve_elapsed_time')
            total_cost = ampl.getObjective("total_cost").value()
            
            print(f"total_cost using ipopt:", total_cost)
            print(f"solve_time using ipopt:", solve_time)


            ampl.close()

#***************************************************************************************
ampl = AMPL()
ampl.read(sys.argv[1])
ampl.read_data(sys.argv[3])

#eps = ampl.getParameter('eps').getValues().to_list()[0]

print("*******************************************************************************\n")

if ipopt_run == 1:
            ampl.getVariable("q").setValues(q_init)
            ampl.getVariable("h").setValues(h_init)
            ampl.getVariable("l").setValues(l_init)

ampl.option["solver"]= sys.argv[2]
ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"

ampl.option["ipopt_options"] = "outlev = 3 expect_infeasible_problem = yes bound_push = 0.001 bound_frac = 0.01 warm_start_init_point = yes halt_on_ampl_error = yes "

#ampl.set_option("ipopt_options", "outlev = 0 expect_infeasible_problem = yes bound_push = 0.001 bound_frac = 0.001 nlp_scaling_method = gradient-based  warm_start_init_point = yes halt_on_ampl_error = yes warm_start_bound_push=1e-9 warm_start_mult_bound_push=1e-9")   #max_iter = 1000
ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 2 warm_start_init_point = no bonmin.num_resolve_at_root = 10 "
ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 300 iisfind = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 basis = 1 " #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1

ampl.option["baron_options"]= "maxtime = 300  outlev = 1 iisfind = 4 lpsolver = cplex lsolver = conopt threads = 8" # lsolver = conopt
ampl.option["scip_options"] = "outlev  1 timelimit 300 wantsol lpmethod = b" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
ampl.option["knitro_options"]= "maxtime_real = 300 outlev = 4 threads=8 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10"
ampl.option["presolve"] = "1"
ampl.option["presolve_eps"] = "8.53e-15"

print(f"{sys.argv[2]} solver outputs:\n")

ampl.solve()

# ampl.eval("expand ;")

# ampl.eval("display {i in 1.._ncons: _con[i].slack < -1e-3} (_conname[i], _con[i].slack);")

l = ampl.getVariable('l').getValues().to_dict()
q = ampl.getVariable('q').getValues().to_dict()
h = ampl.getVariable('h').getValues().to_dict()
eps = ampl.getVariable('eps').getValues().to_dict()

print("eps:",eps, "\n")


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

constraint_violations(q, h, l, R, d, pipes, eps, sys.argv[2])

solve_time = ampl.get_value('_solve_elapsed_time')
total_cost = ampl.getObjective("total_cost").value()

print(f"total_cost using {sys.argv[2]}:", total_cost)
print(f"solve_time using {sys.argv[2]}:", solve_time)

ampl.close()

print("*******************************************************************************\n")

