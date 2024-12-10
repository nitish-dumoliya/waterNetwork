from amplpy import AMPL
import time

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

def nlpModel(data):
    nlp_ampl = AMPL()
    nlp_ampl.reset()    
    nlp_ampl.read("../lpNlp/NLP.mod")
    nlp_ampl.read_data(f"../../data/{data}.dat")
    nlp_ampl.option["solver"] = "knitro"
    nlp_ampl.option["knitro_options"] = "outlev 0"
    nlp_ampl.option["presolve_eps"] = "3.41e-14"
    return nlp_ampl

def lpModel(data):
    lp_ampl = AMPL()
    lp_ampl.reset()
    lp_ampl.read("../lpNlp/lp_model.mod")
    lp_ampl.read_data(f"../../data/{data}.dat")
    lp_ampl.option["solver"] = "cplexamp"
    #lp_ampl.option["cplex_options"] = "display = 2 mipbasis=1 basis_cond=1 logfile = cplex.out endbasis = foo sensitivity"
    lp_ampl.option["presolve_eps"]= "1.28e-14"
    return lp_ampl


def dual_lpModel(data):
    lp_ampl = AMPL()
    lp_ampl.reset()
    lp_ampl.read("dual_lp.mod")
    lp_ampl.read_data(f"../../data/{data}.dat")
    lp_ampl.option["solver"] = "cplexamp"
    lp_ampl.option["presolve_eps"]= "1.28e-14"
    return lp_ampl


def flowModel(data):
    flow_ampl = AMPL()
    flow_ampl.reset()
    flow_ampl.read("flowCon.mod")
    flow_ampl.read_data(f"../../data/{data}.dat")
    flow_ampl.option["solver"] = "cplexamp"
    flow_ampl.option["presolve_eps"]= "1.44e-9"
    return flow_ampl


start = time.time()
data = data_list[0]

nlp_ampl = nlpModel(data)

max_dia=int(nlp_ampl.getSet('pipes').getValues().to_list()[-1])
print(max_dia)
nlp_ampl.eval("s.t. fix_length{(i,j) in arcs}:"f" l[i,j,{max_dia}]=L[i,j];")
# nlp_ampl.eval("s.t. fix_length1{(i,j) in arcs, k in 1...5}: l[i,j,k]=0;")

nlp_ampl.solve()
nlp_ampl.eval("display l;")
# nlp_ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>0.0001} l[i,j,k];")
nlp_ampl.eval("display q;")
nlp_ampl.eval("display h;")
# nlp_ampl.eval("display _var.dual;")
# nlp_ampl.eval("display total_cost;")
totalcost = nlp_ampl.get_objective("total_cost")
print("Objective is:", totalcost.value())

# dual_lp_ampl = dual_lpModel(data)
#
# q_lp = nlp_ampl.getVariable("q").getValues().toDict()
#
# for (i, j), value in q_lp.items():
#     dual_lp_ampl.param['q'][i, j] = value
#
# # lp_ampl.option["ipopt_options"]="outlev 1"
# dual_lp_ampl.solve()
# # lp_ampl.eval("display total_cost;"
# totalcost = dual_lp_ampl.get_objective("total_cost")
# print("Objective is:", totalcost.value())
#
# dual_lp_ampl.eval("display lambda;")
# dual_lp_ampl.eval("display beta;")
# dual_lp_ampl.eval("display gamma;")
# dual_lp_ampl.eval("display con1.dual;")

lp_ampl = lpModel(data)

q_lp = nlp_ampl.getVariable("q").getValues().toDict()
for (i, j), value in q_lp.items():
    lp_ampl.param['q_lp'][i, j] = value

#lp_ampl.eval("s.t. l_238: l_lp[2,3,8]=1000;")
#lp_ampl.eval("s.t. l_759: l_lp[7,5,9]=1000;")

lp_ampl.solve()
#lp_ampl.eval("display l_lp;")
# lp_ampl.eval("display h_lp;")
# lp_ampl.eval("display q_lp;")
lp_ampl.eval("display con1.dual;")
#lp_ampl.eval("display con2.dual;")
# lp_ampl.eval("display con3.dual;")
# lp_ampl.eval("display con4.dual;")
# lp_ampl.eval("show;")
# lp_ampl.eval("expand;")
# lp_ampl.eval("display _var;")
# lp_ampl.eval("display _con;")
# lp_ampl.eval("display _obj;")
# lp_ampl.eval("display {j in 1.._nvars} (_varname[j],_var[j],_var[j].ub,_var[j].rc);")
# lp_ampl.eval("display {j in 1.._ncons} (_conname[j],_con[j],_con[j].lb, _con[j].ub);")

print("Iteration 1")

flow_ampl = flowModel(data)
q_lp = lp_ampl.getParameter('q_lp').getValues().toDict()

#**************************************************************************#
# print(q_lp)

for (i, j), value in q_lp.items():
    flow_ampl.param['q_lp'][i, j] = value

delta = 66

#flow_ampl.eval(f"s.t. fix_q_del45: q_del[4,5]={-delta};")
flow_ampl.eval(f"s.t. fix_q_del75: q_del[7,5]={delta};")
#flow_ampl.eval(f"s.t. fix_q_del46: q_del[4,6]={delta};")
#flow_ampl.eval(f"s.t. fix_q_del67: q_del[6,7]={delta};")

flow_ampl.solve()
flow_ampl.eval("display q_del;")
flow_ampl.eval("display {(i,j) in arcs} q_lp[i,j]+q_del[i,j];")

lp_ampl = lpModel(data)

q_del = flow_ampl.getVariable("q_del").getValues().toDict()
q_lp = flow_ampl.getParameter("q_lp").getValues().toDict()
for (i, j) in q_lp.keys():
    lp_ampl.param['q_lp'][i, j] = q_lp[i,j]+q_del[i,j]

lp_ampl.solve()
lp_ampl.eval("display total_cost;")
lp_ampl.eval("display q_lp;")
lp_ampl.eval("display h_lp;")
lp_ampl.eval("display con1.dual;")

#**************************************************************************#
print("Iteration 2")

flow_ampl = flowModel(data)
q_lp = lp_ampl.getParameter('q_lp').getValues().toDict()

#print(q_lp)

for (i, j), value in q_lp.items():
    flow_ampl.param['q_lp'][i, j] = value

delta = -42.055

#flow_ampl.eval(f"s.t. fix_q_del23: q_del[2,3]={-delta};")
flow_ampl.eval(f"s.t. fix_q_del24: q_del[2,4]={delta};")
#flow_ampl.eval(f"s.t. fix_q_del45: q_del[4,5]={delta};")
#flow_ampl.eval(f"s.t. fix_q_del35: q_del[3,5]={-delta};")

flow_ampl.solve()
flow_ampl.eval("display q_del;")
flow_ampl.eval("display {(i,j) in arcs} q_lp[i,j]+q_del[i,j];")


lp_ampl = lpModel(data)

q_del = flow_ampl.getVariable("q_del").getValues().toDict()
q_lp = flow_ampl.getParameter("q_lp").getValues().toDict()
for (i, j) in q_lp.keys():
    lp_ampl.param['q_lp'][i, j] = q_lp[i,j]+q_del[i,j]

lp_ampl.solve()
lp_ampl.eval("display total_cost;")
lp_ampl.eval("display q_lp;")
lp_ampl.eval("display h_lp;")
lp_ampl.eval("display con1.dual;")

#**************************************************************************#
print("Iteration 3")

flow_ampl = flowModel(data)
q_lp = lp_ampl.getParameter('q_lp').getValues().toDict()

#print(q_lp)

for (i, j), value in q_lp.items():
    flow_ampl.param['q_lp'][i, j] = value

# delta = 0.11001
delta = 0.109 

#flow_ampl.eval(f"s.t. fix_q_del45: q_del[4,5]={-delta};")
#flow_ampl.eval(f"s.t. fix_q_del75: q_del[7,5]={delta};")
flow_ampl.eval(f"s.t. fix_q_del46: q_del[4,6]={delta};")
#flow_ampl.eval(f"s.t. fix_q_del67: q_del[6,7]={delta};")

flow_ampl.solve()
flow_ampl.eval("display q_del;")
flow_ampl.eval("display {(i,j) in arcs} q_lp[i,j]+q_del[i,j];")


lp_ampl = lpModel(data)

q_del = flow_ampl.getVariable("q_del").getValues().toDict()
q_lp = flow_ampl.getParameter("q_lp").getValues().toDict()
for (i, j) in q_lp.keys():
    lp_ampl.param['q_lp'][i, j] = q_lp[i,j]+q_del[i,j]

lp_ampl.solve()
lp_ampl.eval("display total_cost;")
lp_ampl.eval("display q_lp;")
lp_ampl.eval("display h_lp;")
lp_ampl.eval("display con1.dual;")
#**************************************************************************#
print("Iteration 4")

flow_ampl = flowModel(data)
q_lp = lp_ampl.getParameter('q_lp').getValues().toDict()

#print(q_lp)

for (i, j), value in q_lp.items():
    flow_ampl.param['q_lp'][i, j] = value

# delta = -0.10
delta = -0.109

#flow_ampl.eval(f"s.t. fix_q_del23: q_del[2,3]={-delta};")
flow_ampl.eval(f"s.t. fix_q_del24: q_del[2,4]={delta};")
#flow_ampl.eval(f"s.t. fix_q_del45: q_del[4,5]={delta};")
#flow_ampl.eval(f"s.t. fix_q_del35: q_del[3,5]={-delta};")

flow_ampl.solve()
flow_ampl.eval("display q_del;")
flow_ampl.eval("display {(i,j) in arcs} q_lp[i,j]+q_del[i,j];")


lp_ampl = lpModel(data)

q_del = flow_ampl.getVariable("q_del").getValues().toDict()
q_lp = flow_ampl.getParameter("q_lp").getValues().toDict()
for (i, j) in q_lp.keys():
    lp_ampl.param['q_lp'][i, j] = q_lp[i,j]+q_del[i,j]

lp_ampl.solve()
lp_ampl.eval("display total_cost;")
lp_ampl.eval("display q_lp;")
lp_ampl.eval("display h_lp;")
lp_ampl.eval("display con1.dual;")

#**************************************************************************#
end = time.time()

print("Total solve time :",end-start)

