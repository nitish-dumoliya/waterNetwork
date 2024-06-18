import networkx as nx
import time
from scipy.optimize import bisect

def loop1(x):
    return (0.001**1.852)*(10.68*1000/((278.28**1.852) *(1.016)**4.87))*((q_value[17,16]+x)*abs(q_value[17,16]+x)**0.852 + (q_value[18,17]+x)*abs(q_value[18,17]+x)**0.852 + (q_value[19,18]+x)*abs(q_value[19,18]+x)**0.852 + (q_value[3,19]+x)*abs(q_value[3,19]+x)**0.852 - (q_value[15,16]-x)*abs(q_value[15,16]-x)**0.852 - (q_value[14,15]-x)*abs(q_value[14,15]-x)**0.852 - (q_value[10,14]-x)*abs(q_value[10,14]-x)**0.852 - (q_value[9,10]-x)*abs(q_value[9,10]-x)**0.852 - (q_value[8,9]-x)*abs(q_value[8,9]-x)**0.852 - (q_value[7,8]-x)*abs(q_value[7,8]-x)**0.852 - (q_value[6,7]-x)*abs(q_value[6,7]-x)**0.852 - (q_value[5,6]-x)*abs(q_value[5,6]-x)**0.852 - (q_value[4,5]-x)*abs(q_value[4,5]-x)**0.852 - (q_value[3,4]-x)*abs(q_value[3,4]-x)**0.852 )

def loop2(x):
    return (0.001**1.852)*(10.68*1000/((278.28**1.852) *(1.016)**4.87))*((q_value[3,20]+x)*abs(q_value[3,20]+x)**0.852 + (q_value[20,23]+x)*abs(q_value[20,23]+x)**0.852 + (q_value[23,24]+x)*abs(q_value[23,24]+x)**0.852 + (q_value[24,25]+x)*abs(q_value[24,25]+x)**0.852 - (q_value[26,25]-x)*abs(q_value[26,25]-x)**0.852 - (q_value[27,26]-x)*abs(q_value[27,26]-x)**0.852 - (q_value[16,27]-x)*abs(q_value[16,27]-x)**0.852 - (q_value[17,16]-x)*abs(q_value[17,16]-x)**0.852 - (q_value[18,17]-x)*abs(q_value[18,17]-x)**0.852 - (q_value[19,18]-x)*abs(q_value[19,18]-x)**0.852 - (q_value[3,19]-x)*abs(q_value[3,19]-x)**0.852 )

def loop3(x):
	return (0.001**1.852)*(10.68*1000/((278.28**1.852) *(1.016)**4.87))*( (q_value[23,28]+x)*abs(q_value[23,28]+x)**0.852 + (q_value[28,29]+x)*abs(q_value[28,29]+x)**0.852 + (q_value[29,30]+x)*abs(q_value[29,30]+x)**0.852 + (q_value[30,31]+ x)*abs(q_value[30,31]+x)**0.852 - (q_value[32,31]-x)*abs(q_value[32,31]-x)**0.852 - (q_value[25,32]-x)*abs(q_value[25,32]-x)**0.852 - (q_value[24,25]-x)*abs(q_value[24,25]-x)**0.852 - (q_value[23,24]-x)*abs(q_value[23,24]-x)**0.852)


from amplpy import AMPL    
ampl = AMPL()
ampl.reset()
ampl.read("fea_flow.mod")
# tree_ampl.read("spi_tree_lp.mod")
input_data_file = f"/home/nitishdumoliya/minotaur/examples/water-network/Data/d2_Sample_input_cycle_hanoi.dat"
# input_data_file = f"/home/nitish/minotaur/examples/water-network/model/Fea_Flow/one_loop.dat"
ampl.read_data(input_data_file)

# ampl.eval("s.t. fix_length{(i,j) in arcs}:l[i,j,14]=L[i,j];")

ampl.option["solver"] = "cplexamp"
# ampl.option["presolve_eps"]= "6.82e-14"
ampl.solve()
# ampl.eval("display q;")
# ampl.eval("display h;")
# ampl.eval("display l;")

iter = 1
print("Iteration : ", iter)
print("Find the initial feasible flow values in every links :")

q_value = ampl.getVariable("q").getValues().toDict()
for (i, j), value in q_value.items():
    print((i,j),round(value,2))

print("Total Head Loss in Loop 1 :",round(loop1(0),8))
print("Total Head Loss in Loop 2 :",round(loop2(0),8))
print("Total Head Loss in Loop 3 :",round(loop3(0),8))

print("We need to change the flow in loops so that total head loss in loops should be zero.")
print(" ")

root_1 = bisect(loop1, -5538.87, 5538.87)

print(root_1)

break_loop12 = False



while break_loop12 == False:
    print("Iteration : ",iter+1)
    root1 = bisect(loop1, -5538.87, 5538.87)
    # print("Approximate root (x where loop1(x) = 0):", root1)

    q_value[17,16]  = q_value[17,16] + root1
    q_value[18,17] = q_value[18,17] + root1
    q_value[19,18] = q_value[19,18] + root1
    q_value[3,19] =  q_value[3,19] + root1
    q_value[15,16] = q_value[15,16] - root1
    q_value[14,15] = q_value[14,15] - root1
    q_value[10,14]  = q_value[10,14] - root1
    q_value[9,10]  = q_value[9,10]  - root1
    q_value[8,9] =   q_value[8,9] - root1
    q_value[7,8] =   q_value[7,8] - root1
    q_value[6,7] =   q_value[6,7] - root1
    q_value[5,6] =   q_value[5,6] - root1
    q_value[4,5] =   q_value[4,5] - root1
    q_value[3,4] =   q_value[3,4] - root1
    
    # print("Flow values in links : ")
    # for (i, j), value in q_value.items():
    #     print((i,j),round(value,2))
    
    # print(" ")
    # print("Total Head Loss in Loop 1 :",round(loop1(0),8))
    # print("Total Head Loss in Loop 2 :",round(loop2(0),8))  
    # print("Total Head Loss in Loop 3 :",round(loop3(0),8))  
    # print(" ")
    
    root2 = bisect(loop2, -5538.87, 5538.87)
    # print("Approximate root (x where loop2(x) = 0):", root2)
    
    q_value[3,20] =  q_value[3,20] + root2
    q_value[20,23] = q_value[20,23] + root2
    q_value[23,24] = q_value[23,24] + root2
    q_value[24,25] = q_value[24,25] + root2
    q_value[26,25] = q_value[26,25] - root2
    q_value[27,26] = q_value[27,26] - root2
    q_value[16,27] = q_value[16,27] - root2
    q_value[17,16]  = q_value[17,16] - root2
    q_value[18,17] = q_value[18,17] - root2
    q_value[19,18] = q_value[19,18] - root2
    q_value[3,19] =  q_value[3,19] - root2
    
    # print(" ")
    # print("Flow values in links : ")
    # for (i, j), value in q_value.items():
    #     print((i,j),round(value,2))
    print(root2)
    if round(root1,8)==0:
        if round(root2,8)==0:
            break_loop23 = False
            root3 = bisect(loop3, -5538.87, 5538.87)
            print("Approximate root (x where loop3(x) = 0):", root3)
            if round(root3,8)==0:
                break
            while break_loop23==False:
                root3 = bisect(loop3, -5538.87, 5538.87)
                # print("Approximate root (x where loop3(x) = 0):", root3)
                # if round(root3,8)==0:
                #     print("Approximate root (x where loop3(x) = 0):", root3)
                #     break_loop23=True
                #     break_loop12=True
                
                q_value[23,28] =   q_value[23,28] + root3
                q_value[28,29] =   q_value[28,29] + root3
                q_value[29,30] =   q_value[29,30] + root3
                q_value[30,31] =   q_value[30,31] + root3
                q_value[32,31] =   q_value[32,31] - root3
                q_value[25,32] =   q_value[25,32] - root3
                q_value[24,25] =   q_value[24,25] - root3
                q_value[23,24] =   q_value[23,24] - root3
                
                # print("Flow values in links : ")
                # for (i, j), value in q_value.items():
                #     print((i,j),round(value,2))
                
                # print(" ")
                # print("Total Head Loss in Loop 1 :",round(loop1(0),8))
                # print("Total Head Loss in Loop 2 :",round(loop2(0),8))  
                # print("Total Head Loss in Loop 3 :",round(loop3(0),8))  
                # print(" ")

                root2 = bisect(loop2, -5538.87, 5538.87)
                # print("Approximate root (x where loop2(x) = 0):", root2)
                
                q_value[3,20] =  q_value[3,20] + root2
                q_value[20,23] = q_value[20,23] + root2
                q_value[23,24] = q_value[23,24] + root2
                q_value[24,25] = q_value[24,25] + root2
                q_value[26,25] = q_value[26,25] - root2
                q_value[27,26] = q_value[27,26] - root2
                q_value[16,27] = q_value[16,27] - root2
                q_value[17,16]  = q_value[17,16] - root2
                q_value[18,17] = q_value[18,17] - root2
                q_value[19,18] = q_value[19,18] - root2
                q_value[3,19] =  q_value[3,19] - root2
                
                # print(" ")
                # print("Flow values in links : ")
                # for (i, j), value in q_value.items():
                #     print((i,j),round(value,2))
                    
                # print(" ")
                # print("Total Head Loss in Loop 1 :",round(loop1(0),8))
                # print("Total Head Loss in Loop 2 :",round(loop2(0),8))  
                # print("Total Head Loss in Loop 3 :",round(loop3(0),8))  
                # print(" ")
                
                if round(root3,8) == 0:
                    print("Approximate root (x where loop3(x) = 0):", root3)
                    if round(root2,8)==0:
                        print("Approximate root (x where loop2(x) = 0):", root2)
                        break_loop23 == True
                        # root1 = bisect(loop1, -5538.87, 5538.87)
            print("Approximate root (x where loop1(x) = 0):", root1)
            if round(root1,8)==0:
                # print("Approximate root (x where loop1(x) = 0):", root1)
                print("Flow value :")
                break_loop12 = True
                
                for (i, j), value in q_value.items():
                    print((i,j),round(value,2))
                    
                print(" ")
                print("Total Head Loss in Loop 1 :",round(loop1(0),8))
                print("Total Head Loss in Loop 2 :",round(loop2(0),8))  
                print("Total Head Loss in Loop 3 :",round(loop3(0),8))  
                print(" ")
                # break
            # else:
            #     break_loop12 = False
    print("Approximate root (x where loop1(x) = 0):", root1)        
    root2 = bisect(loop2, -5538.87, 5538.87)
    print("Approximate root (x where loop2(x) = 0):", root2)
    root3 = bisect(loop3, -5538.87, 5538.87)
    print("Approximate root (x where loop3(x) = 0):", root3)
    iter = iter + 1
    
ampl = AMPL()
ampl.reset()
# ampl.read("fea_flow.mod")
ampl.read("/home/nitishdumoliya/minotaur/examples/water-network/model/NLP.mod")
input_data_file = f"/home/nitishdumoliya/minotaur/examples/water-network/Data/d2_Sample_input_cycle_hanoi.dat"
# input_data_file = f"/home/nitish/minotaur/examples/water-network/model/Fea_Flow/one_loop.dat"
ampl.read_data(input_data_file)

ampl.eval("s.t. fix_length{(i,j) in arcs}:l[i,j,6]=L[i,j];")

for (i, j), value in q_value.items():
    # print(f"s.t. fix_q_{i}_{j}: q[{i},{j}]= {value};")
    ampl.eval(f"s.t. fix_q_{i}_{j}: q[{i},{j}]= {value};")

ampl.option["solver"] = "cplexamp"
# ampl.option["presolve_eps"]= "1.71e-14"
ampl.option["presolve_eps"]= "8.19e-10"
# ampl.option["presolve"]= "1"

ampl.solve()
# ampl.eval("display q;")
ampl.eval("display h;")
# ampl.eval("display l;"