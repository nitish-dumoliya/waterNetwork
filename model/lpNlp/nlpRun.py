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

print("  ")
print("#**********************************************************************************#")    
print("Water Network File : ",data_list[0])
print(" ")

ampl = AMPL()
ampl.reset()
ampl.read("nlp_model.mod")
input_data_file = f"/home/nitishdumoliya//waterNetwork/data/{data_list[0]}.dat"

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

##################################################################################################
ampl.eval("s.t. length{(i,j) in arcs}: l[i,j,14]=L[i,j] ;")

########################## exhibit the model that has been built ###################################

# ampl.eval("show;")
# ampl.eval("expand;")

####################################################################################################
print("======================Solver Results====================")
ampl.option["solver"] = "knitro"
# ampl.option["solver"] = "/home/nitish/minotaur/build/bin/mmultistart"
# ampl.set_option("mmultistart_options","--presolve 1,--log_level 6,--eval_within_bnds 1")
# ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 0 "
# ampl.option["ipopt_options"] = " outlev = 0"
# ampl.option["knitro_options"] = "outlev = 1 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 20 ms_maxtime_real = 50"
ampl.option["knitro_options"] = "outlev =0 ms_enable 1  ms_maxsolves 10 mip_multistart 1 "
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
# ampl.eval("display z1;")
# ampl.eval("display z2;")
ampl.eval("display con1.dual;")
ampl.eval("display con2.dual;")
ampl.eval("display con3.dual;")
ampl.eval("display con4.dual;")
ampl.eval("display total_cost;")

totalcost = ampl.get_objective("total_cost")
print("Objective:", totalcost.value())
print("==========================================================")

lp_ampl = AMPL()
lp_ampl.reset()
lp_ampl.read("updatedLp.mod")
input_data_file = f"/home/nitishdumoliya//waterNetwork/data/{data_list[0]}.dat"
lp_ampl.read_data(input_data_file)

q_lp = ampl.getVariable("q").getValues().toDict()

for (i, j), value in q_lp.items():
    lp_ampl.param['q_lp'][i, j] = value

lp_ampl.option["presolve_eps"] = "6.82e-14"
lp_ampl.option["solver"] = "knitro"
lp_ampl.solve()
lp_ampl.eval("display h_lp;")
lp_ampl.eval("display con1.dual;")
lp_ampl.eval("display con2.dual;")
lp_ampl.eval("display con3.dual;")
lp_ampl.eval("display con4.dual;")
# lp_ampl.eval("display con5.dual;")
###################################################################

end_time = time.time()
elapsed_time = end_time - start_time
print("elapsed_time : ", elapsed_time)
print("#*******************************************************************************#")
