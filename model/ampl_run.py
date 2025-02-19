# # Install dependencies
import sys
# import pandas as pd
# import csv
# import argparse

# pip install -q amplpy pandas numpy
from amplpy import AMPL
ampl = AMPL()

# data_list=[  #"d1_Sample_input_cycle_twoloop.dat",
            #  "d2_Sample_input_cycle_hanoi.dat",
            # "d14_NewYork.dat",
        #    "d15_foss_poly_0.dat",
        #    "d16_foss_iron.dat",
        #    "d17_foss_poly_1.dat",
        #    ]

# OBJ = []
# csv_data = []
# for data in data_list:
            
ipopt_run = sys.argv[4]

if ipopt_run == 1:
            ampl.reset()
            # ampl.read("non_convex_multiple1.mod")
            ampl.read(sys.argv[1])
            ampl.read_data(sys.argv[3])
            # ampl.read_data(f"/home/nitishdumoliya/minotaur/examples/water-network/Data/{data}")
            #ampl.eval("minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];")
            # ampl.option["solver"]= "/home/nitishdumoliya/minotaur/build-d/bin/mqg"
            # ampl.set_option("mqg_options","--presolve 1,--log_level 6, --nlp_engine IPOPT, --eval_within_bnds 1")

            # Set solver to IPOPT and solve
            ampl.option["solver"] = "ipopt"

            ampl.set_option("ipopt_options", "outlev = 0")   #max_iter = 1000

            ampl.solve()

            total_cost = ampl.getObjective("total_cost").value()
            print("total_cost:", total_cost, "\n")
            print("*******************************************************************************")


            l_init = ampl.getVariable('l').getValues().to_dict()
            q_init = ampl.getVariable('q').getValues().to_dict()
            h_init = ampl.getVariable('h').getValues().to_dict()

ampl.reset()
# ampl.read("non_convex_multiple1.mod")
ampl.read(sys.argv[1])
ampl.read_data(sys.argv[3])

if ipopt_run == 1:
            ampl.getVariable("q").setValues(q_init)
            ampl.getVariable("h").setValues(h_init)
            ampl.getVariable("l").setValues(l_init)


# ampl.option["solver"]= "/home/nitishdumoliya/minotaur/build/bin/mmultistart"

ampl.option["solver"]= sys.argv[2]
ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"
#ampl.option["gurobi_options"] = "outlev 1"
# ampl.option["presolve_eps"] = "1.09e-12"

# ampl.option["solver"] = "knitro"
# ampl.option["solver"] = "/home/nitishdumoliya/dist/bin/ipopt"
#ampl.option["ipopt_options"] = "print_level 3"
# ampl.set_option("ipopt_options", "outlev = 1 expect_infeasible_problem = yes bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = yes max_iter = 3000 halt_on_ampl_error = yes")   #max_iter = 1000
#ampl.option["ipopt_options"] = "outlev = 6 expect_infeasible_problem = yes bound_push = 0.01 bound_frac = 0.01 nlp_scaling = gradient-based warm_start_init_point = yes max_iter = 3000 halt_on_ampl_error = yes"

ampl.set_option("ipopt_options", "outlev = 0 expect_infeasible_problem = yes bound_push = 0.001 bound_frac = 0.001 nlp_scaling_method = gradient-based  warm_start_init_point = yes halt_on_ampl_error = yes warm_start_bound_push=1e-9 warm_start_mult_bound_push=1e-9")   #max_iter = 1000
ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 2 warm_start_init_point = no bonmin.num_resolve_at_root = 10 "
# ampl.eval("option gurobi_auxfiles rc;")
ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iisfind = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 " #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1
#ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 NumericFocus = 1 iisfind = 1 iismethod 0 networkalg = 1 networkcuts = 1 seed 5 varbranch = 3 cuts = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1

ampl.option["baron_options"]= "maxtime = 3600  outlev = 1 iisfind = 2 lpsolver = cplex lsolver = conopt" # lsolver = conopt
ampl.option["scip_options"] = "outlev  1 timelimit 3600 wantsol lpmethod = b" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
ampl.option["knitro_options"]= "maxtime_real = 3600 outlev = 4 threads=8 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10"
#ampl.option["knitro_options"]= "maxtime_real = 3600 outlev = 1  feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 0 ms_maxsolves = 0"
ampl.option["presolve"] = "1"
ampl.option["presolve_eps"] = "8.53e-15"



# Provide initial values
#for (i, j), value in q_init.items():
#    ampl.getVariable("q").setValues({(i, j): value})

#for i, value in h_init.items():
#    ampl.getVariable("h").setValues({i: value})

#for (i, j, k), value in l_init.items():
#    ampl.getVariable("l").setValues({(i, j, k): value})



ampl.solve()


# ampl.eval("expand ;")
# Define the file path where you want to save the results
# output_file = f"/home/nitishdumoliya/minotaur/examples/water-network/Data/multi_out/{data}.txt"

# ampl.setOption('tee', output_file)

# Write the results to the file
# with open(output_file, 'w') as file:
#     sys.stdout = file
#     # ampl.solve()
#     totalcost = ampl.get_objective("total_cost")
#     # totalcost = ampl.get_objective("Objective")
#     print("Objective is:", totalcost.value())    
#     OBJ.append(totalcost.value())
#     ampl.eval("display total_cost;")
#     # ampl.eval('display l;')
#     ampl.eval("display {(i,j) in arcs, k in pipes : l[i,j,k]>0.1}: l[i,j,k];")
#     ampl.eval('display q;')
#     ampl.eval('display h;')
#     ampl.eval("display {(i,j) in arcs} 1000*q[i,j]/((3.14/4)*((sum{k in pipes} d[k]*l[i,j,k])/L[i,j])^2);")
#     ampl.eval("display {(i,j) in arcs} 1000*q[i,j]/(sum{k in pipes} (3.14/4)*(d[k]*X[i,j,k])^2);")

# csv_data.append([data_list, OBJ])
# fields = ["Instances","Objective"]

# rows = csv_data
# filename = "/home/nitishdumoliya/minotaur/examples/water-network/Data/multi_out/Water_instances_mmultistart_results.csv"

# with open(filename, 'w') as csvfile:

#     csvwriter = csv.writer(csvfile)

#     csvwriter.writerow(fields)

#     csvwriter.writerows(rows)

# df = pd.read_csv(filename)

# excel_file = pd.ExcelWriter('/home/nitishdumoliya/minotaur/examples/water-network/Data/multi_out/Water_instances_mmultistart_results.xlsx')
# df.to_excel(excel_file, index = False)

# excel_file.save()
# excel_file = pd.ExcelWriter('/home/nitishdumoliya/minotaur/examples/water-network/Data/multi_out/Water_instances_mmultistart_results.xlsx')

# ampl.close()

# df1 = ampl.get_variable("l")
# print("Pipe length : ", df1.get_values())
# ampl.eval("display X;")
# ampl.eval("display {(i,j) in arcs, k in pipes : X[i,j,k]=1}: X[i,j,k];")
# ampl.eval("display {(i,j) in arcs, k in pipes }: l[i,j,k];")

# ampl.eval("display l;")
# ampl.eval("display q;")
# ampl.eval("display {(i,j) in arcs}: q1[i,j]+q2[i,j];")
# ampl.eval("display h;")
# ampl.eval("display total_cost;")
# ampl.eval("display x;")
# df2 = ampl.get_variable("q")
# print("Flow Values : ",df2.get_values())
# df3 = ampl.get_variable("h")
# print("head loss : ", df3.get_values())
# ampl.eval("display {i in 1.._ncons: _con[i].slack < -1e-3} (_conname[i], _con[i].slack);")
# ampl.eval("display {(i,j) in arcs} 1000*q[i,j]/(sum{k in pipes} (3.14/4)*(d[k]*X[i,j,k])^2);")
# ampl.eval("display {(i,j) in arcs} 1000*q[i,j]/(sum{k in pipes: l[i,j,k] >1} (3.14/4)*(d[k]*X[i,j,k])^2);")
 
# ampl.eval("display {(i,j) in arcs} 1000*q[i,j]/((3.14/4)*((sum{k in pipes} d[k]*l[i,j,k])/L[i,j])^2);")

# ampl.eval("display {(i,j) in arcs }: vmax[i,j]*(sum{k in pipes } (3.14/4)*(d[k])^2);")

def constraint_relative_gap(data_file, q_values, h_values, l_values):
    ampl = AMPL()

    ampl.read("original-nlp.mod")
    ampl.readData(sys.argv[3])

    approx_solution = {
        'q': q_values,  
        'h': h_values,  
        'l': l_values,  
    }

    for (i, j), value in approx_solution['q'].items():
        ampl.getVariable("q").setValues({(i, j): value})

    for i, value in approx_solution['h'].items():
        ampl.getVariable("h").setValues({i: value})

    for (i, j, k), value in approx_solution['l'].items():
        ampl.getVariable("l").setValues({(i, j, k): value})

    relative_gaps = {}

    epsilon = 1e-6  # Small value to prevent division by zero

    for (i, j) in ampl.getSet("arcs"):
        # Original constraint value from AMPL
        original_lhs = ampl.getValue(f"h[{i}]") - ampl.getValue(f"h[{j}]")
        original_rhs = ampl.getValue(f"q[{i},{j}]") * (abs(ampl.getValue(f"q[{i},{j}]"))) ** 0.852 * (0.001 ** 1.852) * \
                       sum(10.67 * ampl.getValue(f"l[{i},{j},{k}]") / ((ampl.getParameter("R")[k] ** 1.852) *
                                                                       ((ampl.getParameter("d")[k] / 1000) ** 4.87))
                           for k in ampl.getSet("pipes"))
        original_value = original_rhs

        # Approximated constraint value using provided values
        approx_lhs = approx_solution['h'][i] - approx_solution['h'][j]
        
        approx_rhs = approx_solution['q'][i, j] * ((abs(approx_solution['q'][i, j])+1000*epsilon) ** 0.852) *(abs(approx_solution['q'][i,j])/(abs(approx_solution['q'][i,j]) + 852*epsilon)) * (0.001 ** 1.852) *  sum(10.67 * approx_solution['l'][i, j, k] / ((ampl.getParameter("R")[k] ** 1.852) * ((ampl.getParameter("d")[k] / 1000) ** 4.87)) for k in ampl.getSet("pipes"))
        approx_value = approx_rhs

        # Compute relative gap
        relative_gap = (original_value - approx_value) / (original_value + epsilon)
        relative_gaps[f"con2_{i},{j}"] = relative_gap

    # Display relative gaps
    for constraint, gap in relative_gaps.items():
        print(f"{constraint}: {gap:.6f}")

l = ampl.getVariable('l').getValues().to_dict()
q = ampl.getVariable('q').getValues().to_dict()
h = ampl.getVariable('h').getValues().to_dict()
        
#constraint_relative_gap(sys.argv[3], q, h, l)

solve_time = ampl.get_value('_solve_elapsed_time')
total_cost = ampl.getObjective("total_cost").value()
print("total_cost:", total_cost)
print("solve_time:", solve_time)

