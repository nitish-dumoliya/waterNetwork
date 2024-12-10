import networkx as nx
import time
import sys

from amplpy import AMPL,Environment

data_list=[
        "data1",
        "data2",
        "twoloop",
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

datanumber = 3

print("  ")
print("****************************** Solver input *************************************")    
print("Water Network File : ",data_list[datanumber])
print(" ")

ampl = AMPL()
ampl.reset()
ampl.read("m1Basic.mod")
# ampl.read("new_spi_tree.mod")
# ampl.read("spi_tree_lp.mod")
input_data_file = f"../data/{data_list[datanumber]}.dat"
ampl.read_data(input_data_file)

########################## exhibit the model that has been built ###################################

# ampl.eval("show;")
# ampl.eval("expand;")

####################################################################################################
print("======================Solver Results====================")
ampl.option["solver"] = "knitro"
# ampl.option["solver"] = "/home/nitishdumoliya/Nitish/minotaur/build/bin/mmultistart"
# ampl.set_option("mmultistart_options","--presolve 1,--log_level 6,--eval_within_bnds 1")
# ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 0 "
# ampl.option["ipopt_options"] = " outlev = 0"
# ampl.option["knitro_options"] = "outlev = 1 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 20 ms_maxtime_real = 50"
ampl.option["knitro_options"] = "outlev = 4 ms_enable 1  ms_maxsolves 10 mip_multistart 1 "
# ampl.option["presolve_eps"]="  6.82e-14"
# ampl.set_option("baron_options","maxtime = 200  outlev = 1 lsolver=knitro firstloc 1 barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")
ampl.set_option("baron_options","maxtime = -1  outlev = 1 ")

# ampl.eval("expand cycle_basis0;")
# ampl.eval("expand cycle_basis1;")
# ampl.option["presolve"]="1"
ampl.solve()
# ampl.eval("show;")
# ampl.eval("expand;")

ampl.eval("display l;")
ampl.eval("display {(i,j) in arcs} sum{k in pipes} C[k]* l[i,j,k];")
ampl.eval("display q;")
ampl.eval("display h;")
ampl.eval("display {(i,j) in arcs} h[i]-h[j];")
# ampl.eval("display total_cost;")
ampl.eval("display con1.dual;")
ampl.eval("display con2.dual;")
ampl.eval("display con5.dual;")
ampl.eval("display con6.dual;")
totalcost = ampl.get_objective("total_cost")
print("Objective:", totalcost.value())
# break

end_time = time.time()
elapsed_time = end_time - start_time
print("elapsed_time : ", elapsed_time)
print("==========================================================")

