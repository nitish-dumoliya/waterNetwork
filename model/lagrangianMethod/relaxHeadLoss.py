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
n = 1
print("  ") 
print("****************************** Solver input ***********************************")
print("Water Network File : ", data_list[n])
print(" ")


def m1BasicModel(n):
    basicAmpl = AMPL()
    basicAmpl.reset()
    basicAmpl.read("../m1Basic.mod")
    input_data_file = f"../../data/{data_list[n]}.dat"
    basicAmpl.read_data(input_data_file)
    basicAmpl.option["solver"] = "ipopt"
    basicAmpl.option["ipopt_options"] = "outlev 0"
    basicAmpl.option["knitro_options"] = "outlev = 0 ms_enable 1  ms_maxsolves 1 mip_multistart 1 "
    basicAmpl.option["presolve_eps"] = "8.53e-15"
    return basicAmpl

def lagrangianRelaxationModel(n):
    ampl = AMPL()
    ampl.reset()
    ampl.read("lagrangianRelaxationModel.mod")
    input_data_file = f"../../data/{data_list[n]}.dat"
    ampl.read_data(input_data_file)
    ampl.option["solver"] = "knitro"
    ampl.option["ipopt_options"] = "outlev 0"
    ampl.option["knitro_options"] = "outlev = 0 ms_enable 1  ms_maxsolves 1 mip_multistart 1 "
    ampl.option["presolve_eps"] = "8.53e-15"
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

def nlpModel(n):
    nlp_ampl = AMPL()
    nlp_ampl.reset()
    nlp_ampl.read("../lpNlp/nlp_model.mod")
    input_data_file = f"../../data/{data_list[n]}.dat"
    nlp_ampl.read_data(input_data_file)
    nlp_ampl.option["solver"] = "knitro"
    nlp_ampl.option["presolve_eps"] = "1.08e-08"
    #nlp_ampl.option["knitro_options"] = "outlev 0 "
    nlp_ampl.option["knitro_options"] = "outlev = 0 ms_enable 1  ms_maxsolves 5 mip_multistart 1 "
    return nlp_ampl


iter = 1
print("Iteration: ",iter)
print("******************** Solve the Lagrangian Relaxation **************************")

ampl = lagrangianRelaxationModel(n)

#arc_set = ampl.getSet("arcs")
#for (i,j) in arc_set.getValues():
#    ampl.eval(f"s.t. x_{i}_{j}: x[{i},{j}] = 0;")

ampl.solve()
ampl.eval("display l;")
ampl.eval("display q;")
ampl.eval("display h;")
ampl.eval("display x;")
totalcost = ampl.get_objective("total_cost")
print("Objective:", totalcost.value())

lb = ampl.getObjective("total_cost").value()

print("===========================Solve the LP Problem================================")

lp_ampl = lpModel(n)

q_lp = ampl.getVariable("q").getValues().toDict()

for (i, j), value in q_lp.items():
    lp_ampl.param['q_lp'][i, j] = value

lp_ampl.solve()
lp_ampl.eval("display h_lp;")
ub = lp_ampl.getObjective("total_cost").value()


print("Lower Bound: ", lb)
print("Upper Bound: ", ub)

lowerBound = lb
upperBound = ub

iter = iter +1
bub = ub
while upperBound-lowerBound >= 0.01:
    print(" ")
    print("Iteration: ",iter)
    print("==============Solve the Surrogate Lagrangian Relaxation Problem================")
    u = ampl.getVariable("u").getValues().toDict()
    x = ampl.getVariable("x").getValues().toDict()
    l = ampl.getVariable("l").getValues().toDict()
    q = ampl.getVariable("q").getValues().toDict()
    h = ampl.getVariable("h").getValues().toDict()
    d = ampl.getParameter("d").getValues().toDict()
    R = ampl.getParameter("R").getValues().toDict()
    E = ampl.getParameter("E").getValues().toDict()
    P = ampl.getParameter("P").getValues().toDict()

    ampl = lagrangianRelaxationModel(n)
    set_pipes = ampl.getSet("pipes")
    set_nodes = ampl.getSet("nodes")

   
    arc_sum=0
    for (i,j) in x.keys():
        g = 0
        for [k] in set_pipes.getValues():
            g = g + 10.67*l[i,j,k]/((R[k]**1.852)*d[k]**4.87)
        arc_sum = arc_sum + (h[i]-h[j]-q[i,j]*(abs(q[i,j])**0.852)*g)    
        
    #steplength = 0.5*(upperBound-lb)/(arc_sum)
    steplength = 2

    for (i,j) in x.keys():
        g = 0
        for [k] in set_pipes.getValues():
             g = g + 10.67*l[i,j,k]/((R[k]**1.852)*d[k]**4.87)
        x[i,j] = x[i,j] + steplength*(h[i]-h[j]-q[i,j]*abs(q[i,j])**(0.852) * g)
        # x[i,j] = max(0,x[i,j])
        ampl.eval(f"s.t. x_{i}_{j}: x[{i},{j}] = {x[i,j]};")
    
    ampl.solve()
    #ampl.eval("display u;")
    ampl.eval("display h;")
    ampl.eval("display q;")
    ampl.eval("display x;")
    ampl.eval("display total_cost;")
    lb = ampl.getObjective("total_cost").value()
    print(lb)

    print(" ")
   
    print(" ")
    print("=============================Solve the LP problem==============================")
    lp_ampl = lpModel(n)
    q_lp = ampl.getVariable("q").getValues().toDict()
    for (i, j), value in q_lp.items():
        lp_ampl.param['q_lp'][i, j] = value
    
   
    lp_ampl.solve()
    lp_ampl.eval("display h_lp;")
    
    ub = lp_ampl.getObjective("total_cost").value()

    upperBound = min(ub, upperBound)
    lowerBound = max(lb, lowerBound)
    
    if upperBound != bub:
        optimalAmpl = lp_ampl
    bub = upperBound

    print("Best Lower Bound: ", lowerBound)
    print("Best Upper Bound: ", upperBound)
    print("Gap: ", upperBound-lowerBound)

    iter = iter + 1
print(" ")
print("********************************Optimal Solution*******************************")
#optimalAmpl.eval("display l_lp;")
#optimalAmpl.eval("display q_lp;")
#optimalAmpl.eval("display h_lp;")
optimalAmpl.eval("display total_cost;")

end_time = time.time()
elapsed_time = end_time - start_time
print("elapsed_time : ", elapsed_time)
