import networkx as nx
import time
import sys
import numpy as np
from amplpy import AMPL, Environment
import contextlib
import os
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


def lagrangianRelaxationModel(n):
    ampl = AMPL()
    ampl.reset()
    ampl.read("lagrangianRelaxationModel.mod")
    input_data_file = f"../../data/{data_list[n]}.dat"
    ampl.read_data(input_data_file)
    ampl.option["solver"] = "ipopt"
    ampl.option["ipopt_options"] = "outlev 0"
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
    lp_ampl.option["solver"] = "cplexamp"
    lp_ampl.option["presolve_eps"] = "8.53e-15"
    return lp_ampl

def m1BasicModel(n):
    basic_ampl = AMPL()
    basic_ampl.reset()
    basic_ampl.read("../m1Basic.mod")
    input_data_file = f"../../data/{data_list[n]}.dat"
    basic_ampl.read_data(input_data_file)
    basic_ampl.option["solver"] = "ipopt"
    basic_ampl.option["ipopt_options"] = "outlev 0"
    basic_ampl.option["presolve_eps"] = "8.53e-15"
    return basic_ampl

# Function to suppress output
@contextlib.contextmanager
def suppress_output():
    # Open devnull to suppress the output
    with open(os.devnull, 'w') as devnull:
        # Redirect stdout and stderr
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            # Restore original stdout and stderr
            sys.stdout = old_stdout
            sys.stderr = old_stderr   
 

iter = 1
# print("Iteration: ",iter)
# print("******************** Solve the Lagrangian Relaxation **************************")

ampl = m1BasicModel(n)

with suppress_output():
    ampl.solve()
# ampl.eval("display l;")
# ampl.eval("display q;")
# ampl.eval("display h;")
# ampl.eval("display u;")
totalcost = ampl.get_objective("total_cost")
# print("Objective:", totalcost.value())

# print("========================Solve the Content Model================================")
l_solution = ampl.getVariable("l").getValues().toDict()
h_solution = ampl.getVariable("h").getValues().toDict()
q_solution = ampl.getVariable("q").getValues().toDict()

basic_ampl = m1BasicModel(n)

for (i, j, k), val in l_solution.items():
    basic_ampl.eval(f'let l[{i},{j},{k}] := {val};')
for (i, j), val in q_solution.items():
    basic_ampl.eval(f'let q[{i},{j}] := {val};')
for i, val in h_solution.items():
    basic_ampl.eval(f'let h[{i}] := {val};')

with suppress_output():
    basic_ampl.solve()

lb = ampl.getObjective("total_cost").value()
ub = basic_ampl.getObjective("total_cost").value()

# print("Upper Bound: ", ub)

lowerBound = lb
upperBound = ub

print(iter,"LowerBound: ", lowerBound, "UpperBound:",upperBound, "ReGap:", round(100*(lowerBound-upperBound)/upperBound, 4))

end_time = time.time()
elapsed_time = end_time - start_time
print("elapsed_time : ", elapsed_time)
