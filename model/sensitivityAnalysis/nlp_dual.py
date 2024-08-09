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

def nlpDualModel(data):
    nlp_dual_ampl = AMPL()
    nlp_dual_ampl.reset()
    nlp_dual_ampl.read("nlp_dual.mod")
    nlp_dual_ampl.read_data(f"../../data/{data}.dat")
    nlp_dual_ampl.option["solver"] = "knitro"
    nlp_dual_ampl.option["knitro_options"] = "outlev = 1 ms_enable = 1 ms_maxsolves = 10"
    nlp_dual_ampl.option["baron_options"] = "outlev = 1 "
    return nlp_dual_ampl

data = data_list[0]

nlp_ampl = nlpDualModel(data)

nlp_ampl.solve()
#nlp_ampl.eval("display l;")
#nlp_ampl.eval("display sum{(i,j) in arcs}(sum{k in nodes} con1.dual[i,j,k]*C[k]);")
nlp_ampl.eval("display q;")
nlp_ampl.eval("display con1.dual;")
nlp_ampl.eval("display con2.dual;")
nlp_ampl.eval("display con3.dual;")
nlp_ampl.eval("display total_cost;")
totalcost = nlp_ampl.get_objective("total_cost")
print("Objective is:", totalcost.value())

