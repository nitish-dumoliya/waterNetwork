import networkx as nx

from amplpy import AMPL

data_list=[
        "data1",
        "d1_Sample_input_cycle_twoloop",
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

networkFile = data_list[0]
print(networkFile)
dataFilePath = f"../../data/{networkFile}.dat"

print(dataFilePath)
print("Water Network File : ",networkFile)
print("  ")

print("#************************ Results Minimum Cost Spanning Tree********************#")    
print(" ")

def arc12(dataFilePath):
    ampl = AMPL()
    #ampl.reset()
    ampl.read("arc12.mod")
    ampl.read_data(dataFilePath)
    ampl.option["solver"] = "cplexamp"
    return ampl

ampl12 = arc12(dataFilePath)
ampl12.solve()
ampl12.eval("display total_cost;")

ampl12.eval("display E[1]-(0.001^1.852)*(Qmax^1.852)*(sum{k in pipes}omega*l[1,2,k]/(R[k]^1.852 * (d[k]/1000)^4.87));")
ampl12.eval("display {k in pipes} : l[1,2,k];")


def spaningTree(dataFilePath):
    ampl = AMPL()
    ampl.reset()
    ampl.read("new_spi_tree.mod")
    ampl.read_data(dataFilePath)
    ampl.option["solver"] = "cplexamp"
    return ampl

# spanningTreeAmpl = spaningTree(dataFilePath)

# spanningTreeAmpl.eval("s.t. cycle: x[2,3]+x[2,4]+x[3,4] = 2;")

# spanningTreeAmpl.eval("s.t. con12: x[1,2] = 1;")
# spanningTreeAmpl.eval("s.t. con23: x[2,3] = 1;")
# spanningTreeAmpl.eval("s.t. con24: x[2,4] = 1;")
 # spanningTreeAmpl.eval("s.t. con34: x[3,4] = 0;")

# spanningTreeAmpl.solve()
# spanningTreeAmpl.eval("display total_cost;")
# spanningTreeAmpl.eval("display l;")
# spanningTreeAmpl.eval("display h;")
# spanningTreeAmpl.eval("display h[2]+h[3]+h[4];")



"""
tree_ampl = AMPL()
tree_ampl.reset()
tree_ampl.read("new_tree_model.mod")
input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data}.dat"
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

tree_ampl.eval("s.t. fix_x_12: x[1,2] = 1;")
tree_ampl.eval("s.t. fix_x_23: x[2,3] = 1;")
tree_ampl.eval("s.t. fix_x_24: x[2,4] = 1;")
tree_ampl.eval("s.t. fix_x_35: x[3,5] = 0;")
tree_ampl.eval("s.t. fix_x_45: x[4,5] = 1;")
tree_ampl.eval("s.t. fix_x_46: x[4,6] = 1;")
tree_ampl.eval("s.t. fix_x_67: x[6,7] = 1;")
tree_ampl.eval("s.t. fix_x_75: x[7,5] = 0;")

####################################################################################################
print("======================Solver Results====================")
tree_ampl.option["solver"] = "cplexamp"
# tree_ampl.option["solver"] = "/home/nitishdumoliya/minotaur/build/bin/mmultistart"
# tree_ampl.set_option("mmultistart_options","--presolve 1,--log_level 6,--eval_within_bnds 1")
# tree_ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 0 "
# tree_ampl.option["ipopt_options"] = " outlev = 0"
tree_ampl.option["knitro_options"] = "outlev = 0 threads=1 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 100 ms_maxtime_real = 30 mip_multistart=1 maxtime_real =60"
# tree_ampl.option["knitro_options"] = "outlev = 3 ms_enable 1  ms_maxsolves 2 "
# tree_ampl.option["presolve_eps"]="  6.82e-14 "
# tree_ampl.set_option("baron_options","maxtime = -1  outlev = 1 lsolver=conopt  barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")

# tree_ampl.option["presolve"]="1"
tree_ampl.solve()
# tree_ampl.eval("show;")
# tree_ampl.eval("expand;")

tree_ampl.eval("display l;")
# tree_ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>1} l[i,j,k];")
tree_ampl.eval("display q;")
tree_ampl.eval("display h;")
# tree_ampl.eval("display x;")
# tree_ampl.eval("display total_cost;")
totalcost = tree_ampl.get_objective("total_cost")
print("Objective for minimum spanning tree is:", totalcost.value())
print("==========================================================")

print(" ")
print("#************************** Results for Loop Network ****************************#")
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
min_k = min(max_k)

for (i,j), value in x_tree.items():
    # print((i,j), value)
    if value ==1:
        PipeNo=[]
        for k in tree_ampl.getSet('pipes'):
            if l_tree[(i,j,k)]>1:
                PipeNo.append(k)
            else:
                continue
        # print(PipeNo)
        loop_ampl.eval(f'set Pipeset_{i}_{j};')
        loop_ampl.set[f'Pipeset_{i}_{j}'] = PipeNo
        # print(loop_ampl.getSet(f'Pipeset_{i}_{j}'))  
        loop_ampl.eval(f"s.t. fixed_len{i}_{j}: sum {{s in Pipeset_{i}_{j}}} l[{i},{j},s] = L[{i},{j}];")
    
        # print(PipeNo)
        # K = max(PipeNo)
        # loop_ampl.eval(f"s.t. fixed_len{i}_{j}: l[{i},{j},{K}] = L[{i},{j}];")
    # else:
    #     print((i,j))
    #     loop_ampl.eval(f"s.t. fixed_len{i}_{j}: l[{i},{j},{max_k[2]}] = L[{i},{j}];")
                
##################################################################################################

# loop_ampl.option["solver"] = "/home/nitishdumoliya/minotaur/build/bin/mmultistart"
# loop_ampl.set_option("mmultistart_options","--presolve 0,--log_level 6,--eval_within_bnds 1")
print("======================Solver Results====================")

# loop_ampl.option["solver"] = "/home/nitishdumoliya/Nitish/ampl.linux-intel64/ipopt"
loop_ampl.option["solver"] = "knitro"
# loop_ampl.option["solver"] = "baron"
# loop_ampl.option["ipopt_options"] = "outlev 3"
loop_ampl.set_option("knitro_options","outlev = 0 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 200 ms_maxtime_real = 60 mip_multistart=1 maxtime_real =60")
# loop_ampl.set_option("baron_options","maxtime = -1  outlev = 1 lsolver=conopt  barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")

loop_ampl.option["presolve_eps"]=" 6e-06 "
loop_ampl.option["presolve"]="1"
loop_ampl.solve()
# loop_ampl.eval("display l;")
loop_ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>0} l[i,j,k];")
# loop_ampl.eval("display q;")
# loop_ampl.eval("display h;")
loop_ampl.eval("display total_cost;")
totalcost = loop_ampl.get_objective("total_cost")
print("Objective for loop network is:", totalcost.value())
# loop_ampl.display("_varname", "_var")
LOOP_COST.append(totalcost.value())
    

print("==========================================================")

print("Tree Objective Value list :", TREE_COST)
print(" ")
print("Loop Objective Value list :", LOOP_COST)

"""
