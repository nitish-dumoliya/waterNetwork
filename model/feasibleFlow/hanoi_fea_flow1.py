import networkx as nx
import time
from scipy.optimize import bisect

def loop1(x):
    return (0.001**1.852)*(10.68*1000/((278.28**1.852) *(1.016)**4.87))*((q_value[17,16]+x)*abs(q_value[17,16]+x)**0.852 + (q_value[18,17]+x)*abs(q_value[18,17]+x)**0.852 + (q_value[19,18]+x)*abs(q_value[19,18]+x)**0.852 + (q_value[3,19]+x)*abs(q_value[3,19]+x)**0.852 - (q_value[15,16]-x)*abs(q_value[15,16]-x)**0.852 - (q_value[14,15]-x)*abs(q_value[14,15]-x)**0.852 - (q_value[10,14]-x)*abs(q_value[10,14]-x)**0.852 - (q_value[9,10]-x)*abs(q_value[9,10]-x)**0.852 - (q_value[8,9]-x)*abs(q_value[8,9]-x)**0.852 - (q_value[7,8]-x)*abs(q_value[7,8]-x)**0.852 - (q_value[6,7]-x)*abs(q_value[6,7]-x)**0.852 - (q_value[5,6]-x)*abs(q_value[5,6]-x)**0.852 - (q_value[4,5]-x)*abs(q_value[4,5]-x)**0.852 - (q_value[3,4]-x)*abs(q_value[3,4]-x)**0.852 )

def loop2(x):
    return (0.001**1.852)*(10.68*1000/((278.28**1.852) *(1.016)**4.87))*((q_value[3,20]+x)*abs(q_value[3,20]+x)**0.852 + (q_value[20,23]+x)*abs(q_value[20,23]+x)**0.852 + (q_value[23,24]+x)*abs(q_value[23,24]+x)**0.852 + (q_value[24,25]+x)*abs(q_value[24,25]+x)**0.852 - (q_value[26,25]-x)*abs(q_value[26,25]-x)**0.852 - (q_value[27,26]-x)*abs(q_value[27,26]-x)**0.852 - (q_value[16,27]-x)*abs(q_value[16,27]-x)**0.852 - (q_value[17,16]-x)*abs(q_value[17,16]-x)**0.852 - (q_value[18,17]-x)*abs(q_value[18,17]-x)**0.852 - (q_value[19,18]-x)*abs(q_value[19,18]-x)**0.852 - (q_value[3,19]-x)*abs(q_value[3,19]-x)**0.852 )

def loop3(x):
	return (0.001**1.852)*(10.68*1000/((278.28**1.852) *(1.016)**4.87))*( (q_value[23,28]+x)*abs(q_value[23,28]+x)**0.852 + (q_value[28,29]+x)*abs(q_value[28,29]+x)**0.852 + (q_value[29,30]+x)*abs(q_value[29,30]+x)**0.852 + (q_value[30,31]+ x)*abs(q_value[30,31]+x)**0.852 - (q_value[32,31]-x)*abs(q_value[32,31]-x)**0.852 - (q_value[25,32]-x)*abs(q_value[25,32]-x)**0.852 - (q_value[24,25]-x)*abs(q_value[24,25]-x)**0.852 - (q_value[23,24]-x)*abs(q_value[23,24]-x)**0.852)

def loop4(x):
    return (0.001**1.852)*(10.68*1000/((278.28**1.852) *(1.016)**4.87))*((q_value[3,20]+x)*abs(q_value[3,20]+x)**0.852 + (q_value[20,23]+x)*abs(q_value[20,23]+x)**0.852 + (q_value[23,28]+x)*abs(q_value[23,28]+x)**0.852 + (q_value[28,29]+x)*abs(q_value[28,29]+x)**0.852 + (q_value[29,30]+x)*abs(q_value[29,30]+x)**0.852 + (q_value[30,31]+ x)*abs(q_value[30,31]+x)**0.852 - (q_value[32,31]-x)*abs(q_value[32,31]-x)**0.852 - (q_value[25,32]-x)*abs(q_value[25,32]-x)**0.852  - (q_value[26,25]-x)*abs(q_value[26,25]-x)**0.852 - (q_value[27,26]-x)*abs(q_value[27,26]-x)**0.852 - (q_value[16,27]-x)*abs(q_value[16,27]-x)**0.852 - (q_value[17,16]-x)*abs(q_value[17,16]-x)**0.852 - (q_value[18,17]-x)*abs(q_value[18,17]-x)**0.852 - (q_value[19,18]-x)*abs(q_value[19,18]-x)**0.852 - (q_value[3,19]-x)*abs(q_value[3,19]-x)**0.852 )

from amplpy import AMPL    
ampl = AMPL()
ampl.reset()
ampl.read("fea_flow.mod")
# tree_ampl.read("spi_tree_lp.mod")
input_data_file = f"/home/nitishdumoliya/Nitish/minotaur/examples/water-network/Data/d2_Sample_input_cycle_hanoi.dat"
# input_data_file = f"/home/nitishdumoliya/Nitish/minotaur/examples/water-network/model/Fea_Flow/one_loop.dat"
ampl.read_data(input_data_file)

# ampl.eval("s.t. fix_length{(i,j) in arcs}:l[i,j,14]=L[i,j];")
ampl.option["solver"] = "cplex"
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
    print((i,j),round(value,4))

print("Total Head Loss in Loop 1 :",round(loop1(0),4))
print("Total Head Loss in Loop 2 :",round(loop2(0),4))
print("Total Head Loss in Loop 3 :",round(loop3(0),4))

print("We need to change the flow in loops so that total head loss in loops should be zero.")
print(" ")

root_1 = bisect(loop1, -5538.87, 5538.87)
root_2 = bisect(loop2, -5538.87, 5538.87)
root_3 = bisect(loop3, -5538.87, 5538.87)

# print(root_1)
# print(root_2)
# print(root_3)

break_loop12 = False

while break_loop12 == False:
    print("Iteration : ",iter+1)
    root1 = bisect(loop1, -5538.87, 5538.87)
    print("Approximate root (x where loop1(x) = 0):", root1)

    q_value[17,16]  = q_value[17,16] + root1
    q_value[18,17]  = q_value[18,17] + root1
    q_value[19,18]  = q_value[19,18] + root1
    q_value[3,19]   = q_value[3,19]  + root1
    q_value[15,16]  = q_value[15,16] - root1
    q_value[14,15]  = q_value[14,15] - root1
    q_value[10,14]  = q_value[10,14] - root1
    q_value[9,10]   = q_value[9,10]  - root1
    q_value[8,9]    = q_value[8,9]   - root1
    q_value[7,8]    = q_value[7,8]   - root1
    q_value[6,7]    = q_value[6,7]   - root1
    q_value[5,6]    = q_value[5,6]   - root1
    q_value[4,5]    = q_value[4,5]   - root1
    q_value[3,4]    = q_value[3,4]   - root1
    
    # print("Flow values in links : ")
    # for (i, j), value in q_value.items():
    #     print((i,j),round(value,2))

    # print(" ")
    # print("Total Head Loss in Loop 1 :",round(loop1(0),8))
    # print("Total Head Loss in Loop 2 :",round(loop2(0),8))  
    # print("Total Head Loss in Loop 3 :",round(loop3(0),8))  
    # print(" ")

    root2 = bisect(loop2, -5538.87, 5538.87)
    print("Approximate root (x where loop2(x) = 0):", root2)
    if round(root2,11)==0:
        root3 = bisect(loop3,-5538.87, 5538.87)
        print("Approximate root (x where loop3(x) = 0):", root3)
        if round(root3,11)==0:
            break_loop12 = True
            print("root1",root1)
            print("root2",root2)
            print("root3",root3)
            
            print("Total Head Loss in Loop 1 :",round(loop1(0),4))
            print("Total Head Loss in Loop 2 :",round(loop2(0),4))  
            print("Total Head Loss in Loop 3 :",round(loop3(0),4))  
            # print("Total Head Loss in Loop 4 :",round(loop4(0),8))  
            
            # root4 = bisect(loop4,-5538.87, 5538.87)
            # print("root4", root4)

        else:
            break_loop23 = False
            while break_loop23==False:
                root3 = bisect(loop3,-5538.87, 5538.87)
                print("Approximate root (x where loop3(x) = 0):", root3)
                q_value[23,28] =   q_value[23,28] + root3
                q_value[28,29] =   q_value[28,29] + root3
                q_value[29,30] =   q_value[29,30] + root3
                q_value[30,31] =   q_value[30,31] + root3
                q_value[32,31] =   q_value[32,31] - root3
                q_value[25,32] =   q_value[25,32] - root3
                q_value[24,25] =   q_value[24,25] - root3
                q_value[23,24] =   q_value[23,24] - root3
                
                root2 = bisect(loop2, -5538.87, 5538.87)
                print("Approximate root (x where loop2(x) = 0):", root2)
                if round(root2,11)==0:
                    break_loop23=True
                else:
                    q_value[3,20]  = q_value[3,20]  + root2
                    q_value[20,23] = q_value[20,23] + root2
                    q_value[23,24] = q_value[23,24] + root2
                    q_value[24,25] = q_value[24,25] + root2
                    q_value[26,25] = q_value[26,25] - root2
                    q_value[27,26] = q_value[27,26] - root2
                    q_value[16,27] = q_value[16,27] - root2
                    q_value[17,16] = q_value[17,16] - root2
                    q_value[18,17] = q_value[18,17] - root2
                    q_value[19,18] = q_value[19,18] - root2
                    q_value[3,19]  = q_value[3,19]  - root2
    else:           
        q_value[3,20]  = q_value[3,20]  + root2
        q_value[20,23] = q_value[20,23] + root2
        q_value[23,24] = q_value[23,24] + root2
        q_value[24,25] = q_value[24,25] + root2
        q_value[26,25] = q_value[26,25] - root2
        q_value[27,26] = q_value[27,26] - root2
        q_value[16,27] = q_value[16,27] - root2
        q_value[17,16] = q_value[17,16] - root2
        q_value[18,17] = q_value[18,17] - root2
        q_value[19,18] = q_value[19,18] - root2
        q_value[3,19]  = q_value[3,19]  - root2

    # print(" ")
    # print("Flow values in links : ")
    # for (i, j), value in q_value.items():
    #     print((i,j),round(value,2))
    iter = iter + 1

print("Flow values in links : ")
for (i, j), value in q_value.items():
    print((i,j),round(value,4))

# q_value[23,24] = 590.1542578702351-0.190163

ampl1 = AMPL()
ampl1.reset()
# ampl.read("fea_flow.mod")
ampl1.read("/home/nitishdumoliya/Nitish/minotaur/examples/water-network/model/NLP.mod")
input_data_file = f"/home/nitishdumoliya/Nitish/minotaur/examples/water-network/Data/d2_Sample_input_cycle_hanoi.dat"
# input_data_file = f"/home/nitishdumoliya/Nitish/minotaur/examples/water-network/model/Fea_Flow/one_loop.dat"
ampl1.read_data(input_data_file)

ampl1.eval("s.t. fix_length{(i,j) in arcs}:l[i,j,6]=L[i,j];")

for (i, j), value in q_value.items():
    print(f"s.t. fix_q_{i}_{j}: q[{i},{j}]= {value};")
    ampl1.eval(f"s.t. fix_q_{i}_{j}: q[{i},{j}]= {round(value,4)};")

ampl1.option["solver"] = "knitro"
# ampl1.option["presolve_eps"]= "1.36e-13"
ampl1.option["presolve_eps"]= "8.19e-10"
# ampl1.option["presolve"]= "1"

ampl1.solve()
# ampl1.eval("display q;")
# ampl1.eval("display h;")
# ampl1.eval("display l;"
# ampl1.eval("expand con2;")
# ampl1.eval("show;")

# q_value = ampl1.getVariable("q").getValues().toDict()
# for (i, j), value in q_value.items():
#     print((i,j),round(value,2))

# # Get the total number of constraints

# all_constraints = ampl1.getConstraints()
# all_constraints = list(all_constraints)

# Print all constraints
# for constraint_name in all_constraints:
#     constraint_expr = ampl1.getConstraint(constraint_name[0])
#     print(constraint_name, ':', constraint_expr)  

# Print variable values
# print("Variable Values:")
# for var in ampl1.getVariables():
#     var_value = ampl1.getVariable(var[0]).getValues()
#     print(f"{var[0]}: {var_value}")

var_value = ampl1.getVariable('q').getValues()
print(f"q: {var_value}")

var_value = ampl1.getVariable('h').getValues()
print(f"h: {var_value}")

# var_value = ampl1.getVariable('l').getValues()
# print(f"l: {var_value}")

objective_value = ampl1.getObjective('total_cost').value()
print("\nObjective Value:", objective_value)

ampl1.eval("display {(i,j) in arcs} (q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * sum{k in pipes } omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87);")
