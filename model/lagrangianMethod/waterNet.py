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

print("  ") 
print("*******************************************************************************")
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
    ampl.option["solver"] = "ipopt"
    ampl.option["ipopt_options"] = "outlev 0"
    ampl.option["knitro_options"] = "outlev = 0 ms_enable 1  ms_maxsolves 5 mip_multistart 1 "
    ampl.option["presolve_eps"] = "8.53e-15"
    return ampl

def contentModel(n):
    content_ampl = AMPL()
    content_ampl.reset()
    content_ampl.read(
        "../potentialBasedFlow/content_model.mod")
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

n = 1
print("Water Network File : ", data_list[n])
print(" ")

iter = 1
print("Iteration: ",iter)
print("************************* Solve the Lagrangian Relaxation *********************")

ampl = lagrangianRelaxationModel(n)

my_set = ampl.getSet("nodes")
for [j] in my_set.getValues():
    if j !=1:
        ampl.eval(f"s.t. u_{j}: u[{j}] = 0;")

ro_value = 1
ro = ampl.getParameter('ro')
ro.set(ro_value)

#ro = ampl.eval("s.t. mu_param:mu=1;")
 
# ampl.eval("show;")
# ampl.eval("expand;")

ampl.solve()
ampl.eval("display l;")
ampl.eval("display q;")
ampl.eval("display h;")
ampl.eval("display u;")
totalcost = ampl.get_objective("total_cost")
lb = ampl.getObjective("total_cost").value()
print("lower bound:", lb)

print("===========================Solve the LP Problem================================")

lp_ampl = lpModel(n)

q_lp = ampl.getVariable("q").getValues().toDict()

for (i, j), value in q_lp.items():
    lp_ampl.param['q_lp'][i, j] = value

lp_ampl.solve()
lp_ampl.eval("display h_lp;")
ub = lp_ampl.getObjective("total_cost").value()
print("upper bound:", ub)

print("==============================Update the bounds================================")
print("Lower Bound: ", lb)
print("Upper Bound: ", ub)

lowerBound = lb
upperBound = ub

optimalAmpl = lp_ampl

iter = iter +1
pub = ub
plb=lb

#while upperBound-lb >= 0.0001:
while iter<=1000:
    print(" ")
    print("Iteration: ",iter)
    print("================Solve the Lagrangian Relaxation Problem====================")
    u = ampl.getVariable("u").getValues().toDict()
    h = ampl.getVariable("h").getValues().toDict()
    d = ampl.getParameter("d").getValues().toDict()
    R = ampl.getParameter("R").getValues().toDict()
    E = ampl.getParameter("E").getValues().toDict()
    P = ampl.getParameter("P").getValues().toDict()
    #ro = ampl.getParameter("ro").getValues().toDict()
    
    ampl = lagrangianRelaxationModel(n)
    set_pipes = ampl.getSet("pipes")
    set_nodes = ampl.getSet("nodes")
   
    ro = ampl.getParameter('ro')
    ro.set(ro_value)

 
    g1=0
    for j in u.keys():
        if j !=1:
            g1=g1+(E[j]+P[j]-h[j])**2
    #for j in u.keys():
    #     print(E[j]+P[j]-h[j])   
    #steplength = 2*(upperBound-lb)/(g1)
    #steplength = 1/(iter)
    #print(steplength)
    #steplength = 1000/(iter)**0.5

    #for j in u.keys():
    #    if j != 1:
    #       u[j] = u[j] + steplength*(max(0,E[j]+P[j]-h[j]))
    #       u[j] = max(0,u[j])
    #       ampl.eval(f"s.t. u_{j}: u[{j}] = {u[j]};")
    for j in u.keys():
        if j != 1:
            #print(j, E[j]+P[j]-h[j])
            if (E[j]+P[j]-h[j]>0):
                #steplength = 0.4*(upperBound-lb)/(g1)
                #ro_value = iter
                u[j] = u[j] + ro_value*(E[j]+P[j]-h[j])
            else: 
                #steplength = 0.5*(upperBound-lb)/(g1)
                #ro_value = iter
                u[j] = u[j] + ro_value*(E[j]+P[j]-h[j])
                #u[j] = 0
            u[j] = max(0,u[j])
            #ampl.eval(f"s.t. u_{j}: u[{j}] = {u[j] + steplength*(E[j]+P[j]-h[j])};")
            ampl.eval(f"s.t. u_{j}: u[{j}] = {u[j]};")
    ampl.solve()
    #ampl.eval("display l;")
    ampl.eval("display u;")
    ampl.eval("display h;")
    ampl.eval("display q;")
    #ampl.eval("display {(i,j) in arcs}: h[i]-h[j];")
    #ampl.eval("display sum{(i,j) in arcs}(sum{k in pipes} l[i,j,k]*C[k]);")
    #ampl.eval("display total_cost;")
    lb = ampl.getObjective("total_cost").value()
    print("lower bound:",lb)
    
    print(" ")
    print("===========================Solve the LP problem============================")
    lp_ampl = lpModel(n)
    q_lp = ampl.getVariable("q").getValues().toDict()
    for (i, j), value in q_lp.items():
        lp_ampl.param['q_lp'][i, j] = value
    
    lp_ampl.solve()
    lp_ampl.eval("display h_lp;")
    
    ub = lp_ampl.getObjective("total_cost").value()
    print("upper bound:",ub)
    print("============================Update the bounds==============================")

    upperBound = min(ub, upperBound)
    lowerBound = max(lb, lowerBound)
    
    if upperBound < pub:
        optimalAmpl = lp_ampl
        optimalAmpl.eval("display con4.dual;")
    pub = upperBound

    #print("Lower Bound: ", lb)
    print("Best Lower Bound: ", lowerBound)
    print("Best Upper Bound: ", upperBound)
    print("Gap: ", upperBound-lowerBound)
    #if lb<=upperBound:
    #    if abs((upperBound-lb)/lb) < 0.01:
    #        lp_ampl.eval("display con4.dual;")
    #        break
    if (upperBound-lowerBound)/lowerBound < 0.001:
        break

    plb = lowerBound

    iter = iter + 1
print(" ")
print("********************************Optimal Solution*******************************")
optimalAmpl.eval("display l_lp;")
optimalAmpl.eval("display q_lp;")
optimalAmpl.eval("display h_lp;")
optimalAmpl.eval("display con1.dual;")
optimalAmpl.eval("display total_cost;")

end_time = time.time()
elapsed_time = end_time - start_time
print("elapsed_time : ", elapsed_time)
print("*******************************************************************************")
