import time

data_list=[
           "d1_Sample_input_cycle_twoloop.dat",
        #    "d2_Sample_input_cycle_hanoi.dat",
        #    "d3_Sample_input_double_hanoi.dat",
        #    "d4_Sample_input_triple_hanoi.dat",
        #    "d5_Taichung_input.dat",
        #    "d6_HG_SP_1_4.dat",
        #    "d7_HG_SP_2_3.dat",
        #    "d8_HG_SP_3_4.dat",
        #    "d9_HG_SP_4_2.dat",
        #    "d10_HG_SP_5_5.dat",
        #    "d11_HG_SP_6_3.dat",
        #    "d12.dat",
        #    "d13.dat",
        #    "d14_NewYork.dat",
        #    "d15_foss_poly_0.dat",
        #    "d16_foss_iron.dat",
        #    "d17_foss_poly_1.dat",
           ]


OBJ=[]
TIME=[]

for data in data_list:
    print("##########################################################################")
    print(f"Results for {data} ")
    start_time = time.time()
    # cost =0
    COST=[]
    
    from amplpy import AMPL
    # add_to_path(r"/home/nitish/Nitish/ampl.linux-intel64")
    nlp_ampl = AMPL()
    nlp_ampl.reset()

    nlp_ampl.read("/home/nitishdumoliya/gittwo/model/lpNlp/NLP.mod")
    nlp_ampl.read_data(f"/home/nitishdumoliya/gittwo/data/{data}")
    nlp_ampl.option["solver"] = "knitro"
    # nlp_ampl.option["ipopt_options"] = "outlev 1"
    
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
    # COST.append(totalcost.value())
    
    break_outer = False

    while break_outer==False:
        lp_ampl = AMPL()
        lp_ampl.reset()
        lp_ampl.read("/home/nitishdumoliya/gittwo/model/lpNlp/lp_model.mod")
        lp_ampl.read_data(f"/home/nitishdumoliya/gittwo/data/{data}")
        q_lp = nlp_ampl.getVariable("q").getValues().toDict()

        for (i, j), value in q_lp.items():
            lp_ampl.param['q_lp'][i, j] = value

        lp_ampl.option["solver"] = "cplex"
        # lp_ampl.option["ipopt_options"]="outlev 1"
        lp_ampl.solve()
        # lp_ampl.eval("display total_cost;")
        totalcost = lp_ampl.get_objective("total_cost")
        print("Objective is:", totalcost.value())

        
        for t_cost in COST:
            if totalcost.value()==t_cost:
                print(t_cost)
                break_outer=True
                break
            else:
                continue
        # if break_outer:
        #     cost=totalcost.value()
        #     COST.append(cost)
        #     break
        
        cost=totalcost.value()
        COST.append(cost)
        lp_ampl.eval("display {(i,j) in arcs, k in pipes:l_lp[i,j,k]>0} l_lp[i,j,k];")
        lp_ampl.eval("display q_lp;")
        lp_ampl.eval("display h_lp;")

        nlp_ampl = AMPL()
        nlp_ampl.reset()
        nlp_ampl.read("/home/nitishdumoliya/gittwo/model/lpNlp/NLP.mod")
        nlp_ampl.read_data(f"/home/nitishdumoliya/gittwo/data/{data}")
        nlp_ampl.option["solver"] = "knitro"
        # nlp_ampl.option["ipopt_options"]="outlev 1"

        l_nlp = lp_ampl.getVariable("l_lp").getValues().toDict()

        my_set = lp_ampl.getSet('arcs')

        for (i,j) in my_set.getValues():
            i = int(i)
            j = int(j)
            # print((i,j))
            K = []
            for k in nlp_ampl.getSet('pipes').getValues().to_list():
                k = int(k)
                if (i,j,k) in list(l_nlp.keys()):
                    # print(l_nlp.get((i,j,k)))
                    if l_nlp.get((i,j,k))>0:
                        K.append(k)
            k = max(K)
            # print(k)
            nlp_ampl.eval(f"s.t. fix_length_{i}_{j}_{k}: l[{i},{j},{k}]=L[{i},{j}];")

        nlp_ampl.solve()
        nlp_ampl.eval("display l;")
        # nlp_ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>0.0001} l[i,j,k];")
        nlp_ampl.eval("display q;")
        nlp_ampl.eval("display h;")
   
    end_time = time.time()
    elapsed_time = end_time - start_time
    # print(cost)
    print(COST)
    print(min(COST))
    TIME.append(elapsed_time)
    OBJ.append(min(COST))