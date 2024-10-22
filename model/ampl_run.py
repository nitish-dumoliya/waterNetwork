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
ampl.reset()
# ampl.read("non_convex_multiple1.mod")
ampl.read(sys.argv[1])
ampl.read_data(sys.argv[2])
# ampl.read_data(f"/home/nitishdumoliya/minotaur/examples/water-network/Data/{data}")

# ampl.option["solver"]= "/home/nitishdumoliya/minotaur/build-d/bin/mqg"
# ampl.set_option("mqg_options","--presolve 1,--log_level 6, --nlp_engine IPOPT, --eval_within_bnds 1")
# ampl.option["solver"]= "baron"
ampl.option["solver"]= "/home/nitishdumoliya/minotaur/build/bin/mmultistart"
ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1"
# ampl.option["presolve_eps"] = "1.09e-12"

# ampl.option["solver"] = "knitro"
# ampl.option["solver"] = "/home/nitishdumoliya/dist/bin/ipopt"
ampl.option["ipopt_options"] = "print_level 3"
ampl.option["baron_options"]= "maxtime = 3600  outlev = 1 lsolver=conopt barstats  objbound  prloc =0 prtime = 100"
ampl.option["knitro_options"]= "maxtime_real = 3600 outlev = 4 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 12 "
ampl.option["presolve"] = "1"
# ampl.option["presolve_eps"] = "1.09e-12"
ampl.solve()

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

ampl.eval("display l;")
ampl.eval("display q;")
ampl.eval("display h;")
# ampl.eval("display total_cost;")
total_cost = ampl.getObjective("total_cost").value()
print("total_cost:", total_cost)
# ampl.eval("display x;")
# df2 = ampl.get_variable("q")
# print("Flow Values : ",df2.get_values())
# df3 = ampl.get_variable("h")
# print("head loss : ", df3.get_values())

# ampl.eval("display {(i,j) in arcs} 1000*q[i,j]/(sum{k in pipes} (3.14/4)*(d[k]*X[i,j,k])^2);")
# ampl.eval("display {(i,j) in arcs} 1000*q[i,j]/(sum{k in pipes: l[i,j,k] >1} (3.14/4)*(d[k]*X[i,j,k])^2);")
 
# ampl.eval("display {(i,j) in arcs} 1000*q[i,j]/((3.14/4)*((sum{k in pipes} d[k]*l[i,j,k])/L[i,j])^2);")

# ampl.eval("display {(i,j) in arcs }: vmax[i,j]*(sum{k in pipes } (3.14/4)*(d[k])^2);")
