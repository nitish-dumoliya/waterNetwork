import networkx as nx

from amplpy import AMPL

data_list=[
        # "data1",
        # "twoloop",
        # "d1_Sample_input_cycle_twoloop",
        "d2_Sample_input_cycle_hanoi",
        # "d3_Sample_input_double_hanoi",
        # "d4_Sample_input_triple_hanoi",
        # "d5_Taichung_input",
        # "d6_HG_SP_1_4",
        # "d7_HG_SP_2_3",
        # "d8_HG_SP_3_4",
        # "d9_HG_SP_4_2",
        # "d10_HG_SP_5_5",
        # "d11_HG_SP_6_3",
        # "d12",
        # "d13",
        # "d14_NewYork",
        # "d15_foss_poly_0",
        # "d16_foss_iron",
        # "d17_foss_poly_1",
        # "d18_pescara",
        # "d19_modena"
        ]

TREE_COST = []
LOOP_COST = []
for data in data_list:
    print("Water Network File : ",data)
    print("  ")

    print("#************************ Results Minimum Cost Spanning Tree********************#")    
    print(" ")

    tree_ampl = AMPL()
    tree_ampl.reset()
    # tree_ampl.read("spanningTreeMultiple.mod")
    # tree_ampl.read("spanningTreeSingle.mod")
    tree_ampl.read("manualPathCont.mod")
    input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data}.dat"
    tree_ampl.read_data(input_data_file)
    nodes_list = []

    for i in tree_ampl.getSet('nodes'):
        nodes_list.append(i)
    # print("Nodes :",nodes_list)
    edge_set = tree_ampl.getSet('arcs')

    edges_list = tree_ampl.getParameter('L').getValues()
    G = edges_list.toDict()
    # print(G)
    for v in sorted(G.items(), key = lambda item:item[1]):
        print(v)
    # {k: v for k, v in sorted(G.items(), key=lambda item: item[1])}
    uwg = nx.Graph()
    uwg.add_nodes_from(nodes_list)
    # print(uwg.nodes())

    uwg.add_weighted_edges_from(edges_list)
    # print(uwg.edges())

    ##################################################################################################
    
    # Calculate a minimum spanning tree of an undirected weighted graph with the kruskal algorithm
    mst = nx.minimum_spanning_tree(uwg, algorithm='kruskal')
    print(mst)
    # G = mst.edges(data=True)
    # print(G)
    # {k: v for k, v in sorted(G.items(), key=lambda item: item[1])}
    # print(sorted(mst.edges(data=True)))
    print("Minimum spanning tree is",mst.edges() )
    print(" ")
    for (i,j) in edge_set:
        if (i,j) in mst.edges():
            print(f"tree_ampl.eval(s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 1;)")
            tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 1;")
        elif (j,i) in mst.edges():
            tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 1;")
        else:
            print((i,j))
            tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 0;")
            tree_ampl.eval(f"s.t. h_{i}: h[{i}] >= E[{i}] + P[{i}];")
            tree_ampl.eval(f"s.t. h_{j}: h[{j}] >= E[{j}] + P[{j}];")
    #voilate = [16,17,20,23,28,29]
    # for i in range(32):
    #     if (i+1 not in voilate):
    #         print(f"s.t. h_{i+1}: h[{i+1}] >= E[{i+1}] + P[{i+1}];")

    tree_ampl.eval("s.t. q1_2: q[1,2] = sum{j in nodes diff Source}  D[j]*v[j];")

    # tree_ampl.eval("s.t. v1: v[1]= 1;")
    # tree_ampl.eval("s.t. v2: v[2]= 1;")
    # tree_ampl.eval("s.t. v3: v[3]= 1;")
    # tree_ampl.eval("s.t. v4: v[4]= 1;")
    # tree_ampl.eval("s.t. v5: v[5]= 1;")
    # tree_ampl.eval("s.t. v6: v[6]= 1;")
    # tree_ampl.eval("s.t. v7: v[7]= 1;")
    # 
    # tree_ampl.eval("s.t. fixed_binary_1_2: x[1,2]= 1;")
    # tree_ampl.eval("s.t. fixed_binary_2_3: x[2,3]= 1;")
    # tree_ampl.eval("s.t. fixed_binary_2_4: x[2,4]= 1;")
    # tree_ampl.eval("s.t. fixed_binary_3_5: x[3,5]= 1;")
    # tree_ampl.eval("s.t. fixed_binary_4_5: x[4,5]= 0;")
    # tree_ampl.eval("s.t. fixed_binary_4_6: x[4,6]= 1;")
    # tree_ampl.eval("s.t. fixed_binary_6_7: x[6,7]= 0;")
    # tree_ampl.eval("s.t. fixed_binary_7_5: x[7,5]= 1;")

    # tree_ampl.eval("show;")
    # tree_ampl.eval("expand;")

    #tree_ampl.eval("s.t. v1: v[1]=   1;")
    #tree_ampl.eval("s.t. v2: v[2]=   1;")
    #tree_ampl.eval("s.t. v3: v[3]=   1;")
    #tree_ampl.eval("s.t. v4: v[4]=   1;")
    #tree_ampl.eval("s.t. v5: v[5]=   1;")
    #tree_ampl.eval("s.t. v6: v[6]=   1;")
    #tree_ampl.eval("s.t. v7: v[7]=   1;")
    #tree_ampl.eval("s.t. v8: v[8]=   1;")
    #tree_ampl.eval("s.t. v9: v[9]=   1;")
    #tree_ampl.eval("s.t. v10: v[10]= 1;")
    #tree_ampl.eval("s.t. v11: v[11]= 1;")
    #tree_ampl.eval("s.t. v12: v[12]= 1;")
    #tree_ampl.eval("s.t. v13: v[13]= 1;")
    #tree_ampl.eval("s.t. v14: v[14]= 1;")
    #tree_ampl.eval("s.t. v15: v[15]= 1;")
    #tree_ampl.eval("s.t. v16: v[16]= 1;")
    #tree_ampl.eval("s.t. v17: v[17]= 1;")
    #tree_ampl.eval("s.t. v18: v[18]= 1;")
    #tree_ampl.eval("s.t. v19: v[19]= 1;")
    #tree_ampl.eval("s.t. v20: v[20]= 1;")
    #tree_ampl.eval("s.t. v21: v[21]= 1;")
    #tree_ampl.eval("s.t. v22: v[22]= 1;")
    #tree_ampl.eval("s.t. v23: v[23]= 1;")
    #tree_ampl.eval("s.t. v24: v[24]= 1;")
    #tree_ampl.eval("s.t. v25: v[25]= 1;")
    #tree_ampl.eval("s.t. v26: v[26]= 1;")
    #tree_ampl.eval("s.t. v27: v[27]= 1;")
    #tree_ampl.eval("s.t. v28: v[28]= 1;")
    #tree_ampl.eval("s.t. v29: v[29]= 0;")
    #tree_ampl.eval("s.t. v30: v[30]= 1;")
    #tree_ampl.eval("s.t. v31: v[31]= 1;")
    #tree_ampl.eval("s.t. v32: v[32]= 1;")
 
    #tree_ampl.eval("s.t. fixed_binary_1_2: x[1,2]=     1;")
    #tree_ampl.eval("s.t. fixed_binary_2_3: x[2,3]=     1;")
    #tree_ampl.eval("s.t. fixed_binary_3_4: x[3,4]=     1;")
    #tree_ampl.eval("s.t. fixed_binary_4_5: x[4,5]=     1;")
    #tree_ampl.eval("s.t. fixed_binary_5_6: x[5,6]=     1;")
    #tree_ampl.eval("s.t. fixed_binary_6_7: x[6,7]=     1;")
    #tree_ampl.eval("s.t. fixed_binary_7_8: x[7,8]=     0;")
    #tree_ampl.eval("s.t. fixed_binary_8_9: x[8,9]=     1;")
    #tree_ampl.eval("s.t. fixed_binary_9_10: x[9,10]=   1;")
    #tree_ampl.eval("s.t. fixed_binary_10_11: x[10,11]= 1;")
    #tree_ampl.eval("s.t. fixed_binary_11_12: x[11,12]= 1;")
    #tree_ampl.eval("s.t. fixed_binary_12_13: x[12,13]= 1;")
    #tree_ampl.eval("s.t. fixed_binary_10_14: x[10,14]= 1;")
    #tree_ampl.eval("s.t. fixed_binary_14_15: x[14,15]= 1;")
    #tree_ampl.eval("s.t. fixed_binary_15_16: x[15,16]= 1;")
    #tree_ampl.eval("s.t. fixed_binary_17_16: x[17,16]= 1;")

    #tree_ampl.eval("s.t. fixed_binary_18_17: x[18,17]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_19_18: x[19,18]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_3_19:  x[3,19]=  1 ;")
    #tree_ampl.eval("s.t. fixed_binary_3_20:  x[3,20]=  1 ;")
    #tree_ampl.eval("s.t. fixed_binary_20_21: x[20,21]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_21_22: x[21,22]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_20_23: x[20,23]= 1 ;")

    #tree_ampl.eval("s.t. fixed_binary_23_24: x[23,24]= 0 ;")
    #tree_ampl.eval("s.t. fixed_binary_24_25: x[24,25]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_26_25: x[26,25]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_27_26: x[27,26]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_16_27: x[16,27]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_23_28: x[23,28]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_28_29: x[28,29]= 0 ;")

    #tree_ampl.eval("s.t. fixed_binary_29_30: x[29,30]= 0 ;")
    #tree_ampl.eval("s.t. fixed_binary_30_31: x[30,31]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_32_31: x[32,31]= 1 ;")
    #tree_ampl.eval("s.t. fixed_binary_25_32: x[25,32]= 1 ;")
    
    # tree_ampl.eval("show;")
    # tree_ampl.eval("expand;")
    
    
    
    # tree_ampl.eval("s.t. fix_x_12: x[1,2] = 1;")
    # tree_ampl.eval("s.t. fix_x_23: x[2,3] = 1;")
    # tree_ampl.eval("s.t. fix_x_24: x[2,4] = 1;")
    # tree_ampl.eval("s.t. fix_x_35: x[3,5] = 1;")
    # tree_ampl.eval("s.t. fix_x_45: x[4,5] = 0;")
    # tree_ampl.eval("s.t. fix_x_46: x[4,6] = 1;")
    # tree_ampl.eval("s.t. fix_x_67: x[6,7] = 1;")
    # tree_ampl.eval("s.t. fix_x_75: x[7,5] = 0;")

    ##################################################################################################
    
    # cycle_basis = nx.cycle_basis(uwg)
    # print(cycle_basis)
    # print(" ")

    # tree_ampl.eval(f"s.t. tree_const: sum {{(i,j) in arcs}} x[i,j] = {len(nodes_list)-1};")

    # count = 0 
    # for cycle in cycle_basis:
    #     # print(cycle)
    #     tree_ampl.eval(f'set cycle{count};')
    #     tree_ampl.set[f'cycle{count}'] = cycle
    #     # print(tree_ampl.getSet(f'cycle{count}'))
    #     # print("s.t. cycle_basis"f{count}": sum{(i,j) in arcs : i,j in "f{cycle}"} x[i,j] <= "f{len(cycle)-1}";")
    #     # print(f"s.t. cycle_basis{count}: sum{{(i,j) in arcs : i in cycle{count} and j in cycle{count}}} x[i,j] <= {len(cycle)-1};" )  
    #     tree_ampl.eval(f"s.t. cycle_basis{count}: sum{{(i,j) in arcs : i in cycle{count} and j in cycle{count}}} x[i,j] <= {len(cycle)-1};" ) 
    #     count = count +1 

    ####################################################################################################
    print("======================Solver Results====================")
    tree_ampl.option["solver"] = "cplexamp"
    tree_ampl.option["cplex_options"] = "iisfind=1 mipdisplay = 5"
    # tree_ampl.option["solver"] = "/home/nitishdumoliya/minotaur/build/bin/mmultistart"
    # tree_ampl.set_option("mmultistart_options","--presolve 1,--log_level 6,--eval_within_bnds 1")
    # tree_ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 0 "
    # tree_ampl.option["ipopt_options"] = " outlev = 0"
    tree_ampl.option["knitro_options"] = "outlev = 0 threads=1 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10 ms_maxtime_real = 30 mip_multistart=1 maxtime_real =60"
    # tree_ampl.option["knitro_options"] = "outlev = 3 ms_enable 1  ms_maxsolves 2 "
    tree_ampl.option["presolve_eps"]="1.64e-12"
    # tree_ampl.set_option("baron_options","maxtime = -1  outlev = 1 lsolver=conopt  barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")

    # tree_ampl.option["presolve"]="1"
    tree_ampl.solve()
    # tree_ampl.eval("show;")
    # tree_ampl.eval("expand;")

    # tree_ampl.eval("display X;")
    # tree_ampl.eval("display x;")
    tree_ampl.eval("display l;")
    #tree_ampl.eval("display {(i,j) in arcs} : sum{k in pipes}l[i,j,k]*C[k] ;")
    # tree_ampl.eval("display q;")
    tree_ampl.eval("display h;")
    # tree_ampl.eval("display {j in 1.._ncons} (_conname[j]);")
    # tree_ampl.eval("display x;")
    # tree_ampl.eval("display sum{i in nodes} h[i];")
    # tree_ampl.eval("display {(i,j) in arcs} h[i]-h[j];")
    # tree_ampl.eval("display sum{(i,j) in arcs} abs(h[i]-h[j]);")
    # tree_ampl.eval("display total_cost;")
    totalcost = tree_ampl.get_objective("total_cost")
    print("Objective for minimum spanning tree is:", totalcost.value())
    TREE_COST.append(totalcost.value())
    print("==========================================================")
    break

    print(" ")
    print("#************************** Results for Loop Network ****************************#")
    print(" ")
    
    loop_ampl= AMPL()
    loop_ampl.reset()
    loop_ampl.read("m1_basic.mod")
    loop_ampl.read_data(input_data_file)

    ###############################################################################################

    l_tree = tree_ampl.getVariable("l").getValues().toDict()

    x_tree = tree_ampl.getVariable("x").getValues().toDict()

    set_pipe = loop_ampl.getSet('pipes')
    print(set_pipe)
    max_k = []
    for k in set_pipe:
        max_k.append(k)
    max_K = max(max_k)
    print(max_K)
    min_k = min(max_k)

    for (i,j), value in x_tree.items():
        # print((i,j), value)
        if value ==1:
            PipeNo=[]
            for k in tree_ampl.getSet('pipes'):
                if l_tree[(i,j,k)]>1:
                    PipeNo.append(k)
                else:
                    continue
            # print(PipeNo)
            loop_ampl.eval(f'set Pipeset_{i}_{j};')
            loop_ampl.set[f'Pipeset_{i}_{j}'] = PipeNo
            # print(loop_ampl.getSet(f'Pipeset_{i}_{j}'))  
            loop_ampl.eval(f"s.t. fixed_len{i}_{j}: sum {{s in Pipeset_{i}_{j}}} l[{i},{j},s] = L[{i},{j}];")
        
            # print(PipeNo)
            # K = max(PipeNo)
            # loop_ampl.eval(f"s.t. fixed_len{i}_{j}: l[{i},{j},{K}] = L[{i},{j}];")
        # else:
        #     print((i,j))
        #     loop_ampl.eval(f"s.t. fixed_len{i}_{j}: l[{i},{j},{max_k[2]}] = L[{i},{j}];")
                    
    ##################################################################################################

    # loop_ampl.option["solver"] = "/home/nitishdumoliya/minotaur/build/bin/mmultistart"
    # loop_ampl.set_option("mmultistart_options","--presolve 0,--log_level 6,--eval_within_bnds 1")
    print("======================Solver Results====================")

    # loop_ampl.option["solver"] = "/home/nitishdumoliya/Nitish/ampl.linux-intel64/ipopt"
    loop_ampl.option["solver"] = "knitro"
    # loop_ampl.option["solver"] = "baron"
    # loop_ampl.option["ipopt_options"] = "outlev 3"
    loop_ampl.set_option("knitro_options","outlev = 0 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 200 ms_maxtime_real = 60 mip_multistart=1 maxtime_real =60")
    # loop_ampl.set_option("baron_options","maxtime = -1  outlev = 1 lsolver=conopt  barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")

    loop_ampl.option["presolve_eps"]=" 6e-06 "
    loop_ampl.option["presolve"]="1"
    loop_ampl.solve()
    # loop_ampl.eval("display l;")
    loop_ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>0} l[i,j,k];")
    # loop_ampl.eval("display q;")
    # loop_ampl.eval("display h;")
    loop_ampl.eval("display total_cost;")
    totalcost = loop_ampl.get_objective("total_cost")
    print("Objective for loop network is:", totalcost.value())
    # loop_ampl.display("_varname", "_var")
    LOOP_COST.append(totalcost.value())
    

print("==========================================================")

print("Tree Objective Value list :", TREE_COST)
print(" ")
print("Loop Objective Value list :", LOOP_COST)
