import networkx as nx
import time
import sys

from amplpy import AMPL,Environment

data_list=[
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

start_time = time.time()

# print("Water Network File : ",sys.argv[1])
print("  ")

print("######################## Results Minimum Cost Spanning Tree##################")    
print(" ")

ampl = AMPL()
ampl.reset()
ampl.read("potential_based_flow.mod")
# ampl.read("new_spi_tree.mod")
# ampl.read("spi_tree_lp.mod")
# input_data_file = f"/home/nitishdumoliya/minotaur/examples/water-network/Data/{sys.argv[1]}"
input_data_file = "/home/nitishdumoliya/waterNetwork/data/d1_Sample_input_cycle_twoloop.dat"
# input_data_file1 = f"/home/nitish/minotaur/examples/water-network/Data/Sample_input_cycle_hanoi.dat"

ampl.read_data(input_data_file)
nodes_list = []

for i in ampl.getSet('nodes'):
    nodes_list.append(i)
# print("Nodes :",nodes_list)
edge_set = ampl.getSet('arcs')
# print("edges:",edge_set)
edges_list = ampl.getParameter('L').getValues()
# print("Edges :",edges_list)

uwg = nx.DiGraph()
uwg.add_nodes_from(nodes_list)
print(uwg.nodes())

uwg.add_weighted_edges_from(edges_list)
print(uwg.edges())

# Define sets for incoming  and outgoing arcs to nodes in AMPL
for node in uwg.nodes():
    print(node)
    incoming_edges = list(uwg.in_edges(node, data=False))
    # print(f"Edges incoming to node {node}:", incoming_edges)
    edge_set_str1 = ", ".join([f"({u},{v})" for u, v in incoming_edges])
    # print("edge_set_str1:",edge_set_str1)
    ampl.eval(f"set delta_neg{node} := {{{edge_set_str1}}};")

    outgoing_edges = list(uwg.out_edges(node, data=False))
    # print(f"Edges outgoing to node {node}:", outgoing_edges)
    edge_set_str = ", ".join([f"({u},{v})" for u, v in outgoing_edges])
    ampl.eval(f"set delta_po{node} := {{{edge_set_str}}};")

    delta_set_neg = ampl.getSet(f'delta_neg{node}')
    print(delta_set_neg)
    delta_set_po = ampl.getSet(f'delta_po{node}')
    print(delta_set_po)
    degree = uwg.degree(node)

    if delta_set_po != {}:
        print(delta_set_po)
        for edge in delta_set_po:
            if node != 1:    
            #     # print("edge : ",edge)
                (i,j) = edge 
                ampl.eval(f"s.t. z_posi_{node}_{i}_{j}:  z1[{i},{j}] <= (sum{{(i,j) in delta_po{node} }} z2[i,j]) - z2[{i},{j}] + (sum{{(i,j) in delta_neg{node}}} z1[i,j]);")
                # if degree == 2:
                #     ampl.eval(f"s.t. z_posi1_{node}_{i}_{j}: z2[{i},{j}] <= (sum{{(i,j) in delta_po{node}}} z1[i,j]) - z1[{i},{j}]  + (sum{{(i,j) in delta_neg{node}}} z2[i,j]);")
    if delta_set_neg != {}:
        for edge in delta_set_neg:
            if node != 1:    
                # print("edge : ",edge)
                (i,j) = edge   
                ampl.eval(f"s.t. z_neg_{node}_{i}_{j}: z2[{i},{j}] <= (sum{{(i,j) in delta_neg{node}}} z1[i,j]) - z1[{i},{j}] + (sum{{(i,j) in delta_po{node}}} z2[i,j]);")
                # if degree == 2:
                #     ampl.eval(f"s.t. z_neg1_{node}_{i}_{j}: z1[{i},{j}] <= (sum{{(i,j) in delta_neg{node}}} z2[i,j]) - z2[{i},{j}] + (sum{{(i,j) in delta_po{node}}} z1[i,j]);")

##################################################################################################

ampl.eval("s.t. cycle1: z1[2,3] + z1[3,5] + z2[4,5] + z2[2,4] <= 3;")
ampl.eval("s.t. cycle2: z1[4,5] + z2[7,5] + z2[6,7] + z2[4,6] <= 3;")

# ampl.eval("s.t. v1: z1[1,2] = 1;")
# ampl.eval("s.t. v2: z1[2,3] = 1;")
# ampl.eval("s.t. v3: z1[3,5] = 1;")
# ampl.eval("s.t. v4: z1[4,5] = 1;")
# ampl.eval("s.t. v5: z1[4,6] = 1;")
# ampl.eval("s.t. v6: z1[6,7] = 1;")
# ampl.eval("s.t. v7: z1[7,5] = 1;")

########################## exhibit the model that has been built ###################################

# ampl.eval("show;")
# ampl.eval("expand;")

####################################################################################################
print("======================Solver Results====================")
ampl.option["solver"] = "knitro"
# ampl.option["solver"] = "/home/nitishdumoliya/Nitish/minotaur/build/bin/mmultistart"
ampl.set_option("mmultistart_options","--presolve 1,--log_level 6,--eval_within_bnds 1")
# ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 0 "
# ampl.option["ipopt_options"] = " outlev = 0"
# ampl.option["knitro_options"] = "outlev = 1 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 20 ms_maxtime_real = 50"
ampl.option["knitro_options"] = "outlev = 4 ms_enable 1  ms_maxsolves 5 mip_multistart 1 "
# ampl.option["presolve_eps"]="  6.82e-14"
# ampl.set_option("baron_options","maxtime = 200  outlev = 1 lsolver=knitro firstloc 1 barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")
ampl.set_option("baron_options","maxtime = -1  outlev = 1 ")

# ampl.eval("expand cycle_basis0;")
# ampl.eval("expand cycle_basis1;")
# ampl.option["presolve"]="1"
ampl.solve()
# ampl.eval("show;")
# ampl.eval("expand;")

# ampl.eval("display l;")
# ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>1} l[i,j,k];")
ampl.eval("display q;")
ampl.eval("display h;")
# ampl.eval("display {(i,j) in arcs} h[i]-h[j];")
ampl.eval("display z1;")
ampl.eval("display z2;")
ampl.eval("display total_cost;")
totalcost = ampl.get_objective("total_cost")
print("Objective:", totalcost.value())
print("==========================================================")
# break

end_time = time.time()
elapsed_time = end_time - start_time
print("==========================================================")

print("elapsed_time : ", elapsed_time)