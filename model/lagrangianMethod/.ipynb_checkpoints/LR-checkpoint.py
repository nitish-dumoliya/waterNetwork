import networkx as nx
import time
import sys

from amplpy import AMPL, Environment

data_list = [
    "data1",
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
print("****************************** Solver input ***********************************")
print("Water Network File : ", data_list[1])
print(" ")

ampl = AMPL()
ampl.reset()
ampl.read("lagrangianRelaxationProblem.mod")
# ampl.read("new_spi_tree.mod")
# ampl.read("spi_tree_lp.mod")
input_data_file = f"../../data/{data_list[1]}.dat"
ampl.read_data(input_data_file)

######################## exhibit the model that has been built #########################


my_set = ampl.getSet('nodes')

#for [j] in my_set.getValues():
#    ampl.eval(f"s.t. u_{j}: u[{j}] = 0;")

# ampl.eval("show;")
# ampl.eval("expand;")

########################################################################################
print("=============================Solver Results====================================")
ampl.option["solver"] = "ipopt"
ampl.option["ipopt_options"] = " outlev = 4"
# ampl.option["presolve_eps"]="  6.82e-14"
# ampl.option["presolve"]="1"

ampl.solve()

ampl.eval("display l;")
# ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>1} l[i,j,k];")
ampl.eval("display q;")
ampl.eval("display h;")
ampl.eval("display x;")
# ampl.eval("display g;")
# ampl.eval("display total_cost;")
totalcost = ampl.get_objective("total_cost") 
print("Objective:", totalcost.value())

print("========================Solve the Content Model================================")

nlp_ampl = AMPL()
nlp_ampl.reset()
nlp_ampl.read("/home/nitishdumoliya/waterNetwork/model/potentialBasedFlow/content_model.mod")
# nlp_ampl.read("/home/nitishdumoliya/waterNetwork/model/lpNlp/NLP.mod")
nlp_ampl.read_data(f"/home/nitishdumoliya/waterNetwork/data/{data_list[1]}.dat")
nlp_ampl.option["solver"] = "ipopt"

nlp_ampl.option["ipopt_options"] = "outlev 1"
l_nlp = ampl.getVariable("l").getValues().toDict()

my_set = ampl.getSet('arcs')

for (i, j, k) in l_nlp.keys():
   if l_nlp[i, j, k] >= 1:
       nlp_ampl.eval(f"s.t. fix_length_{i}_{j}_{
                     k}: l[{i},{j},{k}]={l_nlp[i, j, k]};")

nlp_ampl.solve()
nlp_ampl.eval("display l;")
nlp_ampl.eval("display sum{(i,j) in arcs} (sum{ k in pipes} l[i,j,k]*C[k]);")
nlp_ampl.eval("display q;")

print("===========================Solve the LP Problem================================")

lp_ampl = AMPL()
lp_ampl.reset()
lp_ampl.read("../lpNlp/lp_model.mod")
input_data_file = f"/home/nitishdumoliya//waterNetwork/data/{data_list[1]}.dat"
lp_ampl.read_data(input_data_file)

q_lp = nlp_ampl.getVariable("q").getValues().toDict()

for (i, j), value in q_lp.items():
    lp_ampl.param['q_lp'][i, j] = value

#lp_ampl.eval("s.t. length{(i,j) in arcs}: l_lp[i,j,14]=L[i,j] ;")
lp_ampl.option["presolve_eps"] = "6.82e-14"
lp_ampl.option["solver"] = "cplex"
lp_ampl.solve()
lp_ampl.eval("display h_lp;")


print("Lower Bound: ", ampl.getObjective("total_cost").value())
print("Upper Bound: ", lp_ampl.getObjective("total_cost").value())

# Solve the Surrogate Lagrangian Relaxation Problem

s_ampl = AMPL()
s_ampl.reset()
s_ampl.read("/home/nitishdumoliya/waterNetwork/model/lagrangianMethod/lagrangianRelaxationProblem.mod")
s_ampl.read_data(f"/home/nitishdumoliya/waterNetwork/data/{data_list[1]}.dat")

u = ampl.getVariable("u").getValues().toDict()
q = nlp_ampl.getVariable("q").getValues().toDict()
l = lp_ampl.getVariable("l_lp").getValues().toDict()
h = lp_ampl.getVariable("h_lp").getValues().toDict()
d = ampl.getParameter("d").getValues().toDict()
R = ampl.getParameter("R").getValues().toDict()
E = ampl.getParameter("E").getValues().toDict()
P = ampl.getParameter("P").getValues().toDict()

set_pipes = ampl.getSet("pipes")

omega = 10.67
set_nodes = ampl.getSet("nodes")
s_ampl.eval("var g{nodes};")
for j in h.keys():
    s_ampl.eval(f"s.t. fix_g1_{j} : g[{j}] = {E[j]+P[j]-h[j]};") 

s_ampl.eval("param theta{nodes}:=1;")

for j in u.keys():
    s_ampl.eval(f"s.t. u_{j}: u[{j}] = {u[j]} - theta[{j}]*g[{j}];")

s_ampl.option["solver"]= "ipopt"
s_ampl.option["ipopt_options"] = "outlev = 1 "

s_ampl.solve()
s_ampl.eval("display total_cost;")
s_ampl.eval("display l;")
s_ampl.eval("display q;")
s_ampl.eval("display h;")
s_ampl.eval("display u;")

end_time = time.time()
elapsed_time = end_time - start_time
print("elapsed_time : ", elapsed_time)
print("===============================================================================")
