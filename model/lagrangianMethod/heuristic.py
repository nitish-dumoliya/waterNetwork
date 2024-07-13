import networkx as nx
import time
import sys
import numpy as np
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
n = 3
print("  ") 
print("****************************** Solver input ***********************************")
print("Water Network File : ", data_list[n])
print(" ")


def lagrangianRelaxationModel(n):
    ampl = AMPL()
    ampl.reset()
    ampl.read("lagrangianRelaxationModel.mod")
    input_data_file = f"../../data/{data_list[n]}.dat"
    ampl.read_data(input_data_file)
    ampl.option["solver"] = "ipopt"
    ampl.option["ipopt_options"] = "outlev 0"
    return ampl


def contentModel(n):
    content_ampl = AMPL()
    content_ampl.reset()
    content_ampl.read(
        "/home/nitishdumoliya/waterNetwork/model/potentialBasedFlow/content_model.mod")
    input_data_file = f"../../data/{data_list[n]}.dat"
    content_ampl.read_data(input_data_file)
    content_ampl.option["solver"] = "ipopt"
    content_ampl.option["ipopt_options"] = "outlev 0"
    content_ampl.option["presolve_eps"] = "1.41e-07"
    return content_ampl


def lpModel(n):
    lp_ampl = AMPL()
    lp_ampl.reset()
    lp_ampl.read("../lpNlp/lp_model.mod")
    input_data_file = f"../../data/{data_list[n]}.dat"
    lp_ampl.read_data(input_data_file)
    lp_ampl.option["solver"] = "cplex"
    lp_ampl.option["presolve_eps"] = "1.08e-08"
    return lp_ampl

iter = 1
print("Iteration: ",iter)
print("******************** Solve the Lagrangian Relaxation **************************")

ampl = lagrangianRelaxationModel(n)

my_set = ampl.getSet("nodes")
for [j] in my_set.getValues():
    ampl.eval(f"s.t. u_{j}: u[{j}] = 0;")

# ampl.eval("show;")
# ampl.eval("expand;")

ampl.solve()
ampl.eval("display l;")
ampl.eval("display q;")
ampl.eval("display h;")
ampl.eval("display u;")
totalcost = ampl.get_objective("total_cost")
print("Objective:", totalcost.value())

print("========================Solve the Content Model================================")

content_ampl = contentModel(n)

l_content = ampl.getVariable("l").getValues().toDict()

for (i, j, k) in l_content.keys():
    if l_content[i, j, k] >= 1:
        content_ampl.eval(f"s.t. fix_length_{i}_{j}_{
                          k}: l[{i},{j},{k}]={l_content[i, j, k]};")

content_ampl.solve()
#content_ampl.eval("display l;")
content_ampl.eval(
    "display sum{(i,j) in arcs} (sum{ k in pipes} l[i,j,k]*C[k]);")
#content_ampl.eval("display q;")

print("===========================Solve the LP Problem================================")

lp_ampl = lpModel(n)

q_lp = content_ampl.getVariable("q").getValues().toDict()

for (i, j), value in q_lp.items():
    lp_ampl.param['q_lp'][i, j] = value

lp_ampl.solve()
#lp_ampl.eval("display h_lp;")
lb = ampl.getObjective("total_cost").value()
ub = lp_ampl.getObjective("total_cost").value()

print("Lower Bound: ", lb)
print("Upper Bound: ", ub)

lowerBound = lb
upperBound = ub

iter = iter +1

while upperBound-lowerBound >= 0.001:
    print(" ")
    print("Iteration: ",iter)
    print("==============Solve the Surrogate Lagrangian Relaxation Problem================")
    u = ampl.getVariable("u").getValues().toDict()
    h = ampl.getVariable("h").getValues().toDict()
    d = ampl.getParameter("d").getValues().toDict()
    R = ampl.getParameter("R").getValues().toDict()
    E = ampl.getParameter("E").getValues().toDict()
    P = ampl.getParameter("P").getValues().toDict()

    ampl = lagrangianRelaxationModel(n)
    set_pipes = ampl.getSet("pipes")
    set_nodes = ampl.getSet("nodes")

    ampl.eval("var g{nodes};")
    g1=0
    for j in u.keys():
        if j !=1:
            g1=g1+(E[i]+P[j]-h[j])**2
    #steplength = 0.5*(upperBound-lb)/(g1)**0.5
    steplength = 0.5*(upperBound-lb)/(g1)
    for j in u.keys():
        u[j] = max(0,u[j] + steplength*(E[j]+P[j]-h[j]))
        #ampl.eval(f"s.t. u_{j}: u[{j}] = {u[j] + steplength*(E[j]+P[j]-h[j])};")
        ampl.eval(f"s.t. u_{j}: u[{j}] = {u[j]};")
    ampl.solve()
    #ampl.eval("display u;")
    #ampl.eval("display h;")
    #ampl.eval("display sum{(i,j) in arcs}(sum{k in pipes} l[i,j,k]*C[k]);")
    lb = ampl.getObjective("total_cost").value()
    print(lb)
    print(" ")
    print("==========================Solve the Content model==============================")
    content_ampl = contentModel(n)
    l_content = ampl.getVariable("l").getValues().toDict()

    for (i, j, k) in l_content.keys():
        if l_content[i, j, k] >= 1:
            content_ampl.eval(f"s.t. fix_length_{i}_{j}_{
                              k}: l[{i},{j},{k}]={l_content[i, j, k]};")
    content_ampl.solve()

    print(" ")
    print("=============================Solve the LP problem==============================")
    lp_ampl = lpModel(n)
    q_lp = content_ampl.getVariable("q").getValues().toDict()
    for (i, j), value in q_lp.items():
        lp_ampl.param['q_lp'][i, j] = value

    lp_ampl.solve()
    #lp_ampl.eval("display h_lp;")

    ub = lp_ampl.getObjective("total_cost").value()
    upperBound = min(ub, upperBound)
    lowerBound = max(lb, lowerBound)

    print("Best Lower Bound: ", lowerBound)
    print("Best Upper Bound: ", upperBound)
    print("diff: ", upperBound-lowerBound)
    iter = iter + 1
end_time = time.time()
elapsed_time = end_time - start_time
print("elapsed_time : ", elapsed_time)

