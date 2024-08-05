from amplpy import AMPL

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
    lp_ampl.option["cplexamp_options"] = "mipbasis 1"
    return lp_ampl


def dual_lpModel(data):
    lp_ampl = AMPL()
    lp_ampl.reset()
    lp_ampl.read("dual_lp.mod")
    lp_ampl.read_data(f"../../data/{data}.dat")
    lp_ampl.option["solver"] = "cplexamp"
    return lp_ampl


def flowModel(data):
    flow_ampl = AMPL()
    flow_ampl.reset()
    flow_ampl.read("flowCon.mod")
    flow_ampl.read_data(f"../../data/{data}.dat")
    flow_ampl.option["solver"] = "cplexamp"
    return flow_ampl

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
# nlp_ampl.eval("display total_cost;")
totalcost = nlp_ampl.get_objective("total_cost")
print("Objective is:", totalcost.value())

dual_lp_ampl = dual_lpModel(data)

q_lp = nlp_ampl.getVariable("q").getValues().toDict()

for (i, j), value in q_lp.items():
    dual_lp_ampl.param['q'][i, j] = value

# lp_ampl.option["ipopt_options"]="outlev 1"
dual_lp_ampl.solve()
# lp_ampl.eval("display total_cost;"
totalcost = dual_lp_ampl.get_objective("total_cost")
print("Objective is:", totalcost.value())

dual_lp_ampl.eval("display lambda;")
dual_lp_ampl.eval("display beta;")
dual_lp_ampl.eval("display gamma;")

lp_ampl = lpModel(data)

q_lp = nlp_ampl.getVariable("q").getValues().toDict()
for (i, j), value in q_lp.items():
    lp_ampl.param['q_lp'][i, j] = value

lp_ampl.solve()
#lp_ampl.eval("display l_lp;")
lp_ampl.eval("display h_lp;")
lp_ampl.eval("display q_lp;")
lp_ampl.eval("display con1.dual;")
#lp_ampl.eval("display con2.dual;")
lp_ampl.eval("display con3.dual;")
lp_ampl.eval("display con4.dual;")

flow_ampl = flowModel(data)
q_lp = lp_ampl.getParameter('q_lp').getValues().toDict()

#print(q_lp)

for (i, j), value in q_lp.items():
    flow_ampl.param['q_lp'][i, j] = value

flow_ampl.eval(f"s.t. fix_q_del67: q_del[7,5]={50};")
#flow_ampl.eval("s.t. fix_q_del45: q_del[4,5]=-66;")
flow_ampl.solve()
flow_ampl.eval("display q_del;")
flow_ampl.eval("display {(i,j) in arcs} q_lp[i,j]+q_del[i,j];")

#flow_ampl.param['q_lp'][1,2] = q_lp[1,2]
#flow_ampl.param['q_lp'][2,3] = q_lp[2,3]
#flow_ampl.param['q_lp'][2,4] = q_lp[2,4]
#flow_ampl.param['q_lp'][3,5] = q_lp[3,5]
#flow_ampl.param['q_lp'][4,5] = q_lp[4,5]
#flow_ampl.param['q_lp'][4,6] = q_lp[4,6]
#flow_ampl.param['q_lp'][6,7] = q_lp[6,7]
#flow_ampl.param['q_lp'][7,5] = q_lp[7,5]

lp_ampl = lpModel(data)

q_lp = flow_ampl.getParameter("q_lp").getValues().toDict()
q_del = flow_ampl.getVariable("q_del").getValues().toDict()
for (i, j) in q_lp.keys():
    lp_ampl.param['q_lp'][i, j] = q_lp[i,j] + q_del[i,j]

lp_ampl.solve()
#lp_ampl.eval("display l_lp;")
lp_ampl.eval("display h_lp;")
lp_ampl.eval("display con1.dual;")

