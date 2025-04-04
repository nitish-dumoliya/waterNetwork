import networkx as nx
import time
import sys

from amplpy import AMPL

data_list=[
        # "d1_Sample_input_cycle_twoloop",
        # "d2_Sample_input_cycle_hanoi",
        # "d3_Sample_input_double_hanoi",
        # "d4_Sample_input_triple_hanoi",
        # "d5_Taichung_input",
        # "d6_HG_SP_1_4",
        # "d7_HG_SP_2_3",
        # "d8_HG_SP_3_4",
        # "d9_HG_SP_4_2",
        # "d10_HG_SP_5_5",
        # "d11_HG_SP_6_3",
        # "d12",
        # "d13",
        # "d14_NewYork",
        # "d15_foss_poly_0",
        # "d16_foss_iron",
        # "d17_foss_poly_1",
        # "d18_pescara",
        # "d19_modena"
        ]

TREE_COST = []
LOOP_COST = []
time_list = []

start_time = time.time()

print("Water Network File : ",sys.argv[1])
print("  ")

print("######################## Results Minimum Cost Spanning Tree##################")    
print(" ")

tree_ampl = AMPL()
tree_ampl.reset()
tree_ampl.read("new_spi_tree.mod")
# tree_ampl.read("spi_tree_lp.mod")
input_data_file = f"/home/nitishdumoliya/minotaur/examples/water-network/Data/{sys.argv[1]}"
# input_data_file1 = f"/home/nitish/minotaur/examples/water-network/Data/Sample_input_cycle_twoloop.dat"
# input_data_file1 = f"/home/nitish/minotaur/examples/water-network/Data/Sample_input_cycle_hanoi.dat"

tree_ampl.read_data(input_data_file)
nodes_list = []

for i in tree_ampl.getSet('nodes'):
    nodes_list.append(i)
# print("Nodes :",nodes_list)
edge_set = tree_ampl.getSet('arcs')

edges_list = tree_ampl.getParameter('L').getValues()
# print("Edges :",edges_list)

uwg = nx.Graph()
uwg.add_nodes_from(nodes_list)
# print(uwg.nodes())

uwg.add_weighted_edges_from(edges_list)
# print(uwg.edges())

##################################################################################################

# Calculate a minimum spanning tree of an undirected weighted graph with the Prim algorithm

# mst = nx.minimum_spanning_tree(uwg, algorithm='prim')
# # print(sorted(mst.edges(data=True)))
# print("Minimum spanning tree is",mst.edges() )
# print(" ")
# for (i,j) in edge_set:
#     if (i,j) in mst.edges():
#         tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 1;")
#     elif (j,i) in mst.edges():
#         tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 1;")
#     else:
#         print((i,j))
#         tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 0;")

##################################################################################################

cycle_basis = nx.cycle_basis(uwg)
print(cycle_basis)
print(" ")

tree_ampl.eval(f"s.t. tree_const: sum {{(i,j) in arcs}} x[i,j] = {len(nodes_list)-1};")

count = 0 
for cycle in cycle_basis:
    # print(cycle)
    tree_ampl.eval(f'set cycle{count};')
    tree_ampl.set[f'cycle{count}'] = cycle
    # print(tree_ampl.getSet(f'cycle{count}'))
    # print("s.t. cycle_basis"f{count}": sum{(i,j) in arcs : i,j in "f{cycle}"} x[i,j] <= "f{len(cycle)-1}";")
    # print(f"s.t. cycle_basis{count}: sum{{(i,j) in arcs : i in cycle{count} and j in cycle{count}}} x[i,j] <= {len(cycle)-1};" )  
    tree_ampl.eval(f"s.t. cycle_basis{count}: sum{{(i,j) in arcs : i in cycle{count} and j in cycle{count}}} x[i,j] <= {len(cycle)-1};" ) 
    count = count +1

####################################################################################################
print("======================Solver Results====================")
tree_ampl.option["solver"] = "knitro"
# tree_ampl.option["solver"] = "/home/nitish/minotaur/build/bin/mmultistart"
# tree_ampl.set_option("mmultistart_options","--presolve 1,--log_level 6,--eval_within_bnds 1")
# tree_ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 0 "
# tree_ampl.option["ipopt_options"] = " outlev = 0"
# tree_ampl.option["knitro_options"] = "outlev = 1 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 20 ms_maxtime_real = 50"
tree_ampl.option["knitro_options"] = "outlev =5 ms_enable 1  ms_maxsolves 10 mip_multistart 1 ms_maxtime_real = 10"
tree_ampl.option["presolve_eps"]="  6.82e-14"
# tree_ampl.set_option("baron_options","maxtime = 200  outlev = 1 lsolver=conopt firstloc 1 barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")

tree_ampl.option["presolve"]="1"
tree_ampl.solve()
# tree_ampl.eval("show;")
# tree_ampl.eval("expand;")

# tree_ampl.eval("display l;")
# tree_ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>1} l[i,j,k];")
# tree_ampl.eval("display q;")
# tree_ampl.eval("display h;")
# tree_ampl.eval("display {(i,j) in arcs} h[i]-h[j];")

tree_ampl.eval("display x;")
tree_ampl.eval("display total_cost;")
totalcost = tree_ampl.get_objective("total_cost")
print("Objective for minimum spanning tree is:", totalcost.value())
TREE_COST.append(totalcost.value())
print("==========================================================")
# break

print(" ")
print("############################# Results for Loop Network #######################")
print(" ")
loop_ampl= AMPL()
loop_ampl.reset()
loop_ampl.read("m1_basic.mod")
loop_ampl.read_data(input_data_file)

###############################################################################################

l_tree = tree_ampl.getVariable("l").getValues().toDict()

x_tree = tree_ampl.getVariable("x").getValues().toDict()

set_pipe = loop_ampl.getSet('pipes')
print(set_pipe)
max_k = []
for k in set_pipe:
    max_k.append(k)
max_K = max(max_k)
print(max_K)

for (i,j), value in x_tree.items():
    # print((i,j), value)
    if value == 1:
        PipeNo=[]
        for k in tree_ampl.getSet('pipes'):
            if l_tree[(i,j,k)]>1:
                PipeNo.append(k)
                # print(k)
            else:
                continue
        # print(PipeNo)
        loop_ampl.eval(f'set Pipeset_{i}_{j};')
        loop_ampl.set[f'Pipeset_{i}_{j}'] = PipeNo
        # print(loop_ampl.getSet(f'Pipeset_{i}_{j}'))  
        loop_ampl.eval(f"s.t. fixed_len{i}_{j}: sum {{s in Pipeset_{i}_{j}}} l[{i},{j},s] = L[{i},{j}];")
    
        # print(PipeNo)
        # K = max(PipeNo)
        # loop_ampl.eval(f"s.t. fixed_len{i}_{j}: l[{i},{j},{k}] = L[{i},{j}];")
    # else:
    #     print((i,j))
    #     loop_ampl.eval(f"s.t. fixed_len{i}_{j}: l[{i},{j},1] = L[{i},{j}];")

            
##################################################################################################

# loop_ampl.option["solver"] = "/home/nitishdumoliya/minotaur/build/bin/mmultistart"
# loop_ampl.set_option("mmultistart_options","--presolve 0,--log_level 6,--eval_within_bnds 1")
print("======================Solver Results====================")

loop_ampl.option["solver"] = "/home/nitishdumoliya/dist/bin/ipopt"
# loop_ampl.option["solver"] = "knitro"
loop_ampl.option["ipopt_options"] = "outlev 1"
# loop_ampl.set_option("knitro_options","outlev = 1 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 100 ms_maxtime_real = 10 mip_multistart=1 maxtime_real =60")

loop_ampl.option["presolve_eps"]=" 6e-06 "
loop_ampl.option["presolve"]="1"
loop_ampl.solve()
# ampl.eval("display l;")
# loop_ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>0} l[i,j,k];")
# loop_ampl.eval("display q;")
# loop_ampl.eval("display h;")
# loop_ampl.eval("display {(i,j) in arcs} h[i]-h[j];")

loop_ampl.eval("display total_cost;")
totalcost = loop_ampl.get_objective("total_cost")
print("Objective for loop network is:", totalcost.value())
# loop_ampl.display("_varname", "_var")
LOOP_COST.append(totalcost.value())

loop_ampl.close()
    
end_time = time.time()
elapsed_time = end_time - start_time
time_list.append(elapsed_time)
print("==========================================================")

print("Tree Objective Value list :", TREE_COST)
print(" ")
print("Loop Objective Value list :", LOOP_COST)
print(" ")
print("Solve Time list : ", time_list)
