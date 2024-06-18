import networkx as nx
import time
from scipy.optimize import bisect

def loop1(x):
	return (0.001**1.852)*(10.68*1000/((130**1.852) *(0.6096)**4.87))*( (q_value[2,4]+x)*abs(q_value[2,4]+x)**0.852 + (q_value[4,5]+x)*abs(q_value[4,5]+x)**0.852 -(q_value[3,5]-x)*abs(q_value[3,5]-x)**0.852 - (q_value[2,3]- x)*abs(q_value[2,3]-x)**0.852 )

def loop2(x):
    return (0.001**1.852)*(10.68*1000/((130**1.852) *(0.6096)**4.87))*( (q_value[4,6]+x)*abs(q_value[4,6]+x)**0.852 + (q_value[6,7]+x)*abs(q_value[6,7]+x)**0.852 + (q_value[7,5]+ x)*abs(q_value[7,5]+x)**0.852 - (q_value[4,5]-x)*abs(q_value[4,5]+x)**0.852)

def loop3(x):
     return (0.001**1.852)*(10.68*1000/((130**1.852) *(0.6096)**4.87))*( (q_value[2,4]+x)*abs(q_value[2,4]+x)**0.852 + (q_value[4,6]+x)*abs(q_value[4,6]+x)**0.852 + (q_value[6,7]+x)*abs(q_value[6,7]+x)**0.852 + (q_value[7,5]+ x)*abs(q_value[7,5]+x)**0.852 -(q_value[3,5]-x)*abs(q_value[3,5]-x)**0.852 - (q_value[2,3]- x)*abs(q_value[2,3]-x)**0.852 )

# root = bisect(func, -311, 311)
# print(root)

from amplpy import AMPL, add_to_path 
add_to_path("/home/nitishdumoliya/Nitish/ampl.linux-intel64/")    
ampl = AMPL()
ampl.reset()
ampl.read("fea_flow.mod")
# tree_ampl.read("spi_tree_lp.mod")
input_data_file = f"/home/nitishdumoliya/Nitish/minotaur/examples/water-network/Data/d1_Sample_input_cycle_twoloop.dat"
# input_data_file = f"/home/nitish/minotaur/examples/water-network/model/Fea_Flow/one_loop.dat"
ampl.read_data(input_data_file)

# ampl.eval("s.t. fix_length{(i,j) in arcs}:l[i,j,14]=L[i,j];")
ampl.option["presolve"]= "0"
ampl.option["solver"] = "cplex"
ampl.option["presolve_eps"]= "6.82e-14"
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

print("Total Head Loss in Loop 1 :",round(loop1(0),4))
print("Total Head Loss in Loop 2 :",round(loop2(0),4))

print("We need to change the flow in loops so that total head loss in loops should be zero.")
print(" ")

root = bisect(loop1, -311, 311)

root_3 = bisect(loop3, -311, 311)

print(root_3)



q_value[2,3]= q_value[2,3] - root_3
q_value[2,4]= q_value[2,4] + root_3
# q_value[4,5]= q_value[4,5] + root_1
q_value[3,5]= q_value[3,5] - root_3

q_value[4,6]= q_value[4,6]+root_3
q_value[6,7]= q_value[6,7]+root_3
q_value[7,5]= q_value[7,5]+root_3


for (i, j), value in q_value.items():
    print((i,j),round(value,4))


root_1 = round(root,1)

# while round(root_f,2) !=0:
break_loop = False 

while break_loop == False:
    print("Iteration : ",iter+1)
    root_1 = bisect(loop1, -311, 311)
    print("Approximate root (x where f(x) = 0):", root_1)

    q_value[2,3]= q_value[2,3] - root_1
    q_value[2,4]= q_value[2,4] + root_1
    q_value[4,5]= q_value[4,5] + root_1
    q_value[3,5]= q_value[3,5] - root_1

    print("Flow values in links : ")
    for (i, j), value in q_value.items():
        print((i,j),round(value,4))
    print(" ")
    print("Total Head Loss in Loop 1 :",round(loop1(0),4))
    print("Total Head Loss in Loop 2 :",round(loop2(0),4))  
    print(" ")
    
    root_2 = bisect(loop2, -311, 311)
    print("Approximate root (x where g(x) = 0):", root_2)
    
    if round(root_2,4)==0:
        break_loop=True
    else:
        q_value[4,6]= q_value[4,6]+root_2
        q_value[6,7]= q_value[6,7]+root_2
        q_value[7,5]= q_value[7,5]+root_2
        q_value[4,5]= q_value[4,5]-root_2
       
        print(" ")
        print("Flow values in links : ")
        for (i, j), value in q_value.items():
            print((i,j),value)

    print(" ")
    print("Total Head Loss in Loop 1 :",round(loop1(0),8))
    print("Total Head Loss in Loop 2 :",round(loop2(0),8))  
    print(" ")
    

    iter = iter + 1 
    print(" ")


ampl = AMPL()
ampl.reset()
# ampl.read("fea_flow.mod")
ampl.read("/home/nitishdumoliya/Nitish/minotaur/examples/water-network/model/NLP.mod")
input_data_file = f"/home/nitishdumoliya/Nitish/minotaur/examples/water-network/Data/d1_Sample_input_cycle_twoloop.dat"
# input_data_file = f"/home/nitish/minotaur/examples/water-network/model/Fea_Flow/one_loop.dat"
ampl.read_data(input_data_file)

ampl.eval("s.t. fix_length{(i,j) in arcs}:l[i,j,14]=L[i,j];")
for (i, j), value in q_value.items():
    print(f"s.t. fix_q_{i}_{j}: q[{i},{j}]= {value};")
    ampl.eval(f"s.t. fix_q_{i}_{j}: q[{i},{j}]= {value};")

ampl.option["solver"] = "cplex"
# ampl.option["knitro_options"] = "outlev 0"
ampl.option["presolve_eps"]= "3.41e-14"
ampl.option["presolve"]= "1"

ampl.solve()
ampl.eval("display q;")
ampl.eval("display h;")
# ampl.eval("display l;")  

