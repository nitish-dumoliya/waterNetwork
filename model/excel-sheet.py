import array
import os
#import commands
import time
import math
import sys
import string
import fileinput
import re
import subprocess
import traceback
from typing import List, Tuple, Union
import shutil

def find_float(arr0, st0, fl0):
    for line in arr0:
        find = re.search(st0, line)
        if (find):
            t = re.search('-*[0-9]*\.*[0-9]+', line)
            if t is not None:
                fl0 = float(t.group())
                return 1, fl0
    return -1, fl0


def search_best(arr0, st0, st1, ab, minmax, fl0):
    if minmax == "min":
        fl0 = 1e+20
    else:
        fl0 = -1e+20
    for line in arr0:
        find = re.search(st0, line)
        if (find):
            find = re.search(st1, line)
            if (find):
                split_line = line.split(' ')
                st1_split = st1.split(' ')
                words = len(st1_split)
                i = 0
                for word in split_line:
                    if word == st1_split[i]:
                        i += 1
                    if i == words:
                        if ab == "a":
                            ind = split_line.index(word) + 2
                            fl1 = split_line[ind]
                            if "inf" in fl1:
                                fl1 = 1e+20 if minmax == "min" else -1e+20
                            else:
                                fl1 = float(fl1)
                            if minmax == "min" and fl1 < fl0:
                                fl0 = fl1
                            elif minmax == "max" and fl1 > fl0:
                                fl0 = fl1
                        else:
                            ind = split_line.index(word) - words
                            fl1 = float(split_line[ind])
                            if minmax == "min" and fl1 < fl0:
                                fl0 = fl1
                            elif minmax == "max"  and fl1 > fl0:
                                fl0 = fl1
                        break
    return 1, fl0


def InstanceOutput(instance,output):
    file_data = output.read().split('\n')
    time = 0
    find, time = find_float(file_data, "MultiStart: wall clock time used", time)
    if find < 0:
        time = 3600
        
    status = ""
    for line in file_data:
        find = re.search("MultiStart: status of branch-and-bound:", line)
        if (find):
                status = line.partition("MultiStart: status of branch-and-bound: ")[2]
                if status == "Optimal solution found":
                    status = "OPT"
                elif status == "Reached time limit":
                    status = "TL"
                elif status == "Detected infeasibility":
                    status = "INFEAS"
                else:
                    print(instance, "THIS SHOULD NOT HAPPEN")
                    sys.exit(0)
    ##UB
    ub = 0
    if status == "":
        find, ub = search_best(file_data, "BranchAndBound:", "ub", "a", "min", ub)
        if ub == 1e+20:
            ub = "INFTY"
    elif status == "INFEAS":
            ub = "INFEAS"
    else:
        find, ub = find_float(file_data, "MultiStart: best solution value", ub)

    ##LB
    lb = 0
    if status == "":
        find, lb = search_best(file_data, "BranchAndBound:", "lb", "a", "max", lb)
        if lb == -1e+20:
            lb = "-INFTY"
    elif status == "INFEAS":
        lb = "INFEAS"
    elif status == "OPT":
        lb = ub
    else:
        find, lb = find_float(file_data, "MultiStart: best bound estimate from remaining nodes", lb)

    ##Gap
    gap = 0
    if status == "":
        find, gap = search_best(file_data, "BranchAndBound:", "gap%", "a", "min", gap)
        if gap == 1e+20:
            gap = "INFTY"
    elif status == "INFEAS":
            gap = "INFEAS"
    elif status == "OPT":
            gap = 0.00
    else:
        find, gap = find_float(file_data, "MultiStart: gap percentage =", gap)
    
    if status == "":
        status = "TL"
    
    # print("model instance:",instance)
    # print("Objective:",ub)
    # print("Time:",time) 
    return ub, time,ub


def BaronInstanceOutput(instance, output):
    
    # Extract total relative constraint violation
    output_content = output.read() 

    violation_match_con = re.search(r'Sum of constraints violation:\s*([\d.eE+-]+)', output_content)
    violation_match_abs = re.search(r'Con2 sum of absolute violation:\s*([\d.eE+-]+)', output_content)
    violation_match_rel = re.search(r'Con2 sum of relative violation:\s*([\d.eE+-]+)', output_content)
    #print(violation_match)
    
    violation_con = None
    if violation_match_con:
        try:
            violation_con = float(violation_match_con.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_con.group(1)}' to float.")
    violation_abs = None
    if violation_match_abs:
        try:
            violation_abs = float(violation_match_abs.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_abs.group(1)}' to float.")
    violation_rel = None
    if violation_match_rel:
        try:
            violation_rel = float(violation_match_rel.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_rel.group(1)}' to float.")
 
 
    #file_data = out.read().split('\n')
    file_data = output_content.split('\n')
    time = 0
    find, time = find_float(file_data, "baron solve time: ", time)
    Objective = 0
    find, Objective = find_float(file_data, "Objective ", Objective)
    # Objective = 0
    
    for line in file_data:
        find1 = re.search("Problem is infeasible", line)
        if (find1):
            print("Problem is infeasible")
            Objective = "Infeasible "
        find2 = re.search("No feasible solution was found", line)
        if (find2):
            print("No feasible solution was found")
            Objective = 1e+51
    return Objective, time, violation_con, violation_abs, violation_rel

def KnitroInstanceOutput(instance, output):
    file_data = output.read().split('\n')
    time = 0
    find, time = find_float(file_data, "Total program time", time)
    Objective = 0
    find, Objective = find_float(file_data, "objective  ", Objective)
    # Objective = 0
    FEA_ERR =0 
    for line in file_data:
        find1 = re.search("objective", line)
        if (find1):
            Objective = line.partition("objective")[2]
            # print(Objective)
            Objective = Objective.partition(";")[0]
            # FEA_ERR= line.partition("feasibility error")[2]
            # print(FEA_ERR)
            # a,b = FEA_ERR.split("e")
            # FEA_ERR = float(a)*(10**float(b))
            # print(FEA_ERR)
            # FEA_ERR = float(FEA_ERR)
            # if FEA_ERR >0:
            #     Objective = "infeasible"

    # print("model instance:",instance)
    # print("Objective:",Objective)
    # print("Lower bound:",lb1)
    # print("Uper bound:",ub1)
    # print("Time:",time) 
    # print(" ")
    return Objective, time


def IpoptInstanceOutput(instance, output):
    output_content = output.read()  # Read the entire output content

    # Extract total relative constraint violation
    violation_match_con = re.search(r'Sum of constraints violation:\s*([\d.eE+-]+)', output_content)
    violation_match_abs = re.search(r'Con2 sum of absolute violation:\s*([\d.eE+-]+)', output_content)
    violation_match_rel = re.search(r'Con2 sum of relative violation:\s*([\d.eE+-]+)', output_content)
    #
    violation_con = float(violation_match_con.group(1)) if violation_match_con else None
    violation_abs = float(violation_match_abs.group(1)) if violation_match_abs else None
    violation_rel = float(violation_match_rel.group(1)) if violation_match_rel else None

    # Extract objective value
    objective_match = re.search(r'Total cost using ipopt:\s*([\d.eE+-]+)', output_content)
    objective = float(objective_match.group(1)) if objective_match else None

    # Extract solve time
    time_match = re.search(r'IPOPT solve time:\s*([\d.eE+-]+)', output_content)
    solve_time = float(time_match.group(1)) if time_match else None

    return objective, solve_time, violation_con, violation_abs, violation_rel


def BonminInstanceOutput(instance, output):
    file_data = output.read().split('\n')
    time = 0
    # find, time = find_float(file_data, "Total CPU secs in IPOPT", time)
    Objective = 0
    find, Objective = find_float(file_data, "total_cost using baron", Objective)
    find, time = find_float(file_data, "solve_time using baron", time)
    
    # T = []
    # for line in file_data:
    #     find1 = re.search("NLP0014I", line)
    #     if (find1):
    #         find2 = re.search("OPT",line)
    #         if find2:
    #             time = line.partition("OPT")[2]
    #             time = str(time).split()[2]
    #             time = float(time)
    #             # print(time)
    #             T.append(time)

    #         find2 = re.search("INFEAS",line)
    #         if find2:
    #             time = line.partition("INFEAS")[2]
    #             time = str(time).split()[2]
    #             time = float(time)
    #             T.append(time)
    #             # print(time)

    #             Objective = "INFEAS"
    # time = sum(T)
    return Objective, time

def GurobiInstanceOutput(instance, output):
    # Match best objective using a robust regex pattern
    #obj_match = re.search(r'Best objective\s+([-\d\.eE]+)', output)
    obj_match = re.search(r'Best objective\s+([-]?\d+\.\d+e[+-]?\d+)', output)
    best_bound = re.search(r"best bound ([\d.e+-]+)", output)
    # Alternative objective format from "optimal solution" line
    obj_alt_match = re.search(r'optimal solution; objective\s+([-\d\.eE]+)', output, re.IGNORECASE)

    # Extract objective safely
    objective = None
    extracted_obj = obj_match.group(1) if obj_match else (obj_alt_match.group(1) if obj_alt_match else None)
    best_bound = best_bound.group(1) if best_bound else None

    no_feasible_solution = "Solution count 0" in output
    
    if no_feasible_solution:
        objective = 'no feasible solution'
    if extracted_obj:
        try:
            # Ensure no truncation and remove any extra spaces or unwanted characters
            extracted_obj = extracted_obj.strip().replace(",", "")
            objective = float(extracted_obj)  # Convert to float
            best_bound = best_bound.strip().replace(",", "")
            best_bound = float(best_bound)  # Convert to float
        except ValueError:
            print(f"Warning: Could not convert extracted objective '{extracted_obj}' to float.")

    # Match solver time (ensure correct extraction)
    time_match = re.search(r'gurobi solve time:\s*([\d\.]+)', output)

    time_limit_reached = "Time limit reached" in output

    solver_time = None
    if time_limit_reached:
        solver_time = 3600
    elif time_match:
        try:
            solver_time = float(time_match.group(1).strip())
        except ValueError:
            print(f"Warning: Could not convert solver time '{time_match.group(1)}' to float.")
   
    # Extract total relative constraint violation
    violation_match_con = re.search(r'Sum of constraints violation:\s*([\d.eE+-]+)', output)
    violation_match_abs = re.search(r'Con2 sum of absolute violation:\s*([\d.eE+-]+)', output)
    violation_match_rel = re.search(r'Con2 sum of relative violation:\s*([\d.eE+-]+)', output)
    #print(violation_match)
    
    violation_con = None
    if violation_match_con:
        try:
            violation_con = float(violation_match_con.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_con.group(1)}' to float.")
    violation_abs = None
    if violation_match_abs:
        try:
            violation_abs = float(violation_match_abs.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_abs.group(1)}' to float.")
    violation_rel = None
    if violation_match_rel:
        try:
            violation_rel = float(violation_match_rel.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_rel.group(1)}' to float.")
  
    return objective, best_bound, solver_time, violation_con, violation_abs, violation_rel

def SCIPInstanceOutput(instance, output): 
    # Match best objective using a robust regex pattern
    #obj_match = re.search(r'Best objective\s+([-\d\.eE]+)', output)
    obj_match = re.search(r"Primal Bound\s*:\s*([+-]?\d+\.\d+e[+-]?\d+)", output)
    dual_bound = re.search(r"Dual Bound\s*:\s*([+-]?\d+\.\d+e[+-]?\d+)", output)
    lp_error = re.search(r"SCIP Error \(-6\): error in LP solver", output)


    # Extract objective safely
    objective = None

    #no_feasible_solution = "Solution count 0" in output
    
    if lp_error:
        objective = 'LP error'
    if obj_match:
        try:
            # Ensure no truncation and remove any extra spaces or unwanted characters
            objective = float(obj_match.group(1))
            dual_bound = float(dual_bound.group(1))
        except ValueError:
            print(f"Warning: Could not convert extracted objective '{extracted_obj}' to float.")
    
    # Match solver time (ensure correct extraction)
    #time_match = re.search(r"Solving Time \(sec\)\s*:\s*([\d\.]+)", output)
    #time_match = re.search(r'Solving Time (sec) :\s*([\d\.]+)', output)

    time_match = re.search(r'scip solve time:\s*([\d\.]+)', output)
    
    time_limit_reached = "time limit reached" in output

    solver_time = None

    if time_match:
        try:
            solver_time = float(time_match.group(1).strip())
        except ValueError:
            print(f"Warning: Could not convert solver time '{time_match.group(1)}' to float.")
    elif time_limit_reached:
        solver_time = 3600
 
    # Extract total relative constraint violation
    violation_match_con = re.search(r'Sum of constraints violation:\s*([\d.eE+-]+)', output)
    violation_match_abs = re.search(r'Con2 sum of absolute violation:\s*([\d.eE+-]+)', output)
    violation_match_rel = re.search(r'Con2 sum of relative violation:\s*([\d.eE+-]+)', output)
    #print(violation_match)
    
    violation_con = None
    if violation_match_con:
        try:
            violation_con = float(violation_match_con.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_con.group(1)}' to float.")
    violation_abs = None
    if violation_match_abs:
        try:
            violation_abs = float(violation_match_abs.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_abs.group(1)}' to float.")
    violation_rel = None
    if violation_match_rel:
        try:
            violation_rel = float(violation_match_rel.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match_rel.group(1)}' to float.")   
    
    return objective, dual_bound, solver_time, violation_con, violation_abs, violation_rel




def HeuristicInstanceOutput(instance, output):
    
    output_content = output.read() 

    # Extract total relative constraint violation
    violation_match = re.search(r'Total relative constraint violation:\s*([\d.eE+-]+)', output_content)
    
    violation = None
    if violation_match:
        try:
            violation = float(violation_match.group(1).strip())  # Convert extracted value to float
        except ValueError:
            print(f"Warning: Could not convert extracted violation '{violation_match.group(1)}' to float.")

    #file_data = output.read().split('\n')
    file_data = output_content.split('\n')
    time = 0
    # find, time = find_float(file_data, "Total CPU secs in IPOPT", time)
    Objective = 0
    find, Objective = find_float(file_data, "Final best objective:", Objective)
    find, time = find_float(file_data, "Heuristic elapsed time:", time)

  
    return Objective, time, violation

data_list=["d1_Sample_input_cycle_twoloop",
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
           "d15_foss_poly_0",
           "d16_foss_iron",
           "d17_foss_poly_1"
           ]

# !python3 ampl_run.py m1_basic.mod d1_Sample_input_cycle_twoloop.dat > output/d1_Sample_input_cycle_twoloop.baro_out

import pandas as pd
import csv

Ins = []
Mmultistart_Objective = []
Mmultistart_Time_taken = []
Baron_Objective = []
Baron_Time_taken = []
Baron_Con_vio = []
Baron_Abs_vio = []
Baron_Rel_vio = []
Knitro_Objective = []
Knitro_Time_taken = []
Ipopt_Objective = []
Ipopt_Time_taken = []
Bonmin_Objective = []
Bonmin_Time_taken = []
Gurobi_Objective = []
Gurobi_Best_Bound = []
Gurobi_Time_taken = []
Gurobi_Con_vio = []
Gurobi_Abs_vio = []
Gurobi_Rel_vio = []
Scip_Objective = []
Scip_Dual_Bound = []
Scip_Time_taken = []
Scip_Con_vio = []
Scip_Abs_vio = []
Scip_Rel_vio = []
Heuristic_Objective = []
Heuristic_Time_taken = []
Heuristic_Rel_vio = []
InitialIpopt_Objective = []
InitialIpopt_Time_taken = []
InitialIpopt_Con_vio = []
InitialIpopt_Abs_vio = []
InitialIpopt_Rel_vio = []

approx_folder = "approx_no_initial"

print("****************************Results of Mmultistart Solver *********************************")

for ins in data_list:
    with open(f"../output/m_out/{ins}.m_out") as output:
        print("Model Name:",ins)
        obj, time,ub = InstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Ins.append(ins)
        Mmultistart_Objective.append(obj)
        Mmultistart_Time_taken.append(time)

print("******************************Results of Initial Ipopt Solve Result ************************************")
for ins in data_list:
    with open(f"../output/baron_out/{approx_folder}/{ins}.baron_out") as output:
        print("Model Name:",ins)
        out = output
        obj, time, violation_con,violation_abs, violation_rel = IpoptInstanceOutput(ins,out)
        print("Objective :",obj)
        print("Time :",time)
        print("Constraint violation :",violation_rel)
        print("Absolute violation :",violation_abs)
        print("Relative violation :",violation_rel)

        print(" ")

        InitialIpopt_Objective.append(obj)
        InitialIpopt_Time_taken.append(time)
        InitialIpopt_Con_vio.append(violation_con)
        InitialIpopt_Abs_vio.append(violation_abs)
        InitialIpopt_Rel_vio.append(violation_rel)

print("******************************Results of Baron Solver ************************************")

for ins in data_list:
    with open(f"../output/baron_out/{approx_folder}/{ins}.baron_out") as output:
        print("Model Name:",ins)
        out = output
        obj, time, violation_con,violation_abs, violation_rel= BaronInstanceOutput(ins,out)
        print("Objective :",obj)
        print("Time :",time)
        print("Constraint violation :",violation_rel)
        print("Absolute violation :",violation_abs)
        print("Relative violation :",violation_rel)

        print(" ")
        Baron_Objective.append(obj)
        Baron_Time_taken.append(time)
        Baron_Con_vio.append(violation_con)
        Baron_Abs_vio.append(violation_abs)
        Baron_Rel_vio.append(violation_rel)

print("**********************Results of Knitro Solver *********************************")
 
for ins in data_list:
    with open(f"../output/knitro_out/{ins}.knitro_out") as output:
        print("Model Name:",ins)
        obj, time = KnitroInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Knitro_Objective.append(obj)
        Knitro_Time_taken.append(time)

#print("**********************Results of IPOPT Solver *********************************")
#
#for ins in data_list:
#    with open(f"../output/ipopt_out/{ins}.ipopt_out") as output:
#        print("Model Name:",ins)
#        #output = output.read()
#        obj, time, violation = IpoptInstanceOutput(ins,output)
#        print("Objective :",obj)
#        print("Time :",time)
#        print("Relative violation", violation)
#        print(" ")
#        Ipopt_Objective.append(obj)
#        Ipopt_Time_taken.append(time)
#
print("**********************Results of BONMIN Solver *********************************")
        
for ins in data_list:
    with open(f"../output/bonmin_out/{ins}.bonmin_out") as output:
        print("Model Name:",ins)
        obj, time = BonminInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Bonmin_Objective.append(obj)
        Bonmin_Time_taken.append(time)
# Solver_name = [" ","mmultistart", " ", "Baron", " ", "Knitro "," "]

print("**********************Results of GUROBI Solver *********************************")

for ins in data_list:
    with open(f"../output/gurobi_out/{approx_folder}/{ins}.gurobi_out") as output:
        print("Model Name:",ins)
        output = output.read()
        obj, best_bound, time, violation_con,violation_abs, violation_rel = GurobiInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Best bound :",best_bound)
        print("Time :",time)
        print("Constraint violation :",violation_rel)
        print("Absolute violation :",violation_abs)
        print("Relative violation :",violation_rel)
        print(" ")
        Gurobi_Objective.append(obj)
        Gurobi_Best_Bound.append(best_bound)
        Gurobi_Time_taken.append(time)
        Gurobi_Con_vio.append(violation_con)
        Gurobi_Abs_vio.append(violation_abs)
        Gurobi_Rel_vio.append(violation_rel)

print("**********************Results of SCIP Solver *********************************")
        
for ins in data_list:
    with open(f"../output/scip_out/{approx_folder}/{ins}.scip_out") as output:
        print("Model Name:",ins)
        output = output.read()
        obj, dual_bound, time, violation_con,violation_abs, violation_rel = SCIPInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Dual Bound :",dual_bound)
        print("Time :",time)
        print("Constraint violation :",violation_rel)
        print("Absolute violation :",violation_abs)
        print("Relative violation :",violation_rel)

        print(" ")
        Scip_Objective.append(obj)
        Scip_Dual_Bound.append(dual_bound)
        Scip_Time_taken.append(time)
        Scip_Con_vio.append(violation_con)
        Scip_Abs_vio.append(violation_abs)
        Scip_Rel_vio.append(violation_rel)


print("**********************Results of Heuristic *********************************")

#for ins in data_list:
#    with open(f"../output/heuristic_out/{approx_folder}/{ins}.heuristic_out") as output:
#        print("Model Name:",ins)
#        obj, time, violation = HeuristicInstanceOutput(ins,output)
#        print("Objective :",obj)
#        print("Time :",time)
#        print("Relative violation :",violation)
#        print(" ")
#        Heuristic_Objective.append(obj)
#        Heuristic_Time_taken.append(time)
#        Heuristic_Rel_vio.append(violation)

fields = ["Instances","Mmultistart Objective","Mmultistart time taken","Ini Ipopt Objective","Ini Ipopt time taken","Ini Ipopt Con Vio","Ini Ipopt Abs Vio","Ini Ipopt Rel Vio","Baron Objective","Baron time taken","Baron Con Vio","Baron Abs Vio","Baron Rel Vio", "Gurobi Objective","Gurobi Best Bound","Gurobi time taken", "Gurobi Con Vio","Gurobi Abs Vio","Gurobi Rel Vio","Scip Objective","Scip Dual Bound", "Scip time taken", "Scip Con Vio", "Scip Abs Vio", "Scip Rel Vio","Knitro Objective","Knitro time taken", "Bonmin Objective","Bonmin time taken", "Heuristic Objective","Heuristic time taken", "Heuristic Rel Vio"  ]

filename = "mmultistart_results.csv"

with open(filename, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    # csvwriter.writerow(Solver_name)
    csvwriter.writerow(fields)

import pandas as pd
csv_input = pd.read_csv(filename)
csv_input['Instances'] = Ins
#csv_input['Mmultistart Objective'] = Mmultistart_Objective
#csv_input['Mmultistart time taken'] = Mmultistart_Time_taken
csv_input['Ini Ipopt Objective'] = InitialIpopt_Objective
csv_input['Ini Ipopt time taken'] = InitialIpopt_Time_taken
csv_input['Ini Ipopt Con Vio'] = InitialIpopt_Con_vio
csv_input['Ini Ipopt Abs Vio'] = InitialIpopt_Abs_vio
csv_input['Ini Ipopt Rel Vio'] = InitialIpopt_Rel_vio
csv_input['Baron Objective'] = Baron_Objective
csv_input['Baron time taken'] = Baron_Time_taken
csv_input['Baron Con Vio'] = Baron_Con_vio
csv_input['Baron Abs Vio'] = Baron_Abs_vio
csv_input['Baron Rel Vio'] = Baron_Rel_vio
csv_input['Gurobi Objective'] = Gurobi_Objective
csv_input['Gurobi Best Bound'] = Gurobi_Best_Bound
csv_input['Gurobi time taken'] = Gurobi_Time_taken
csv_input['Gurobi Con Vio'] = Gurobi_Con_vio
csv_input['Gurobi Abs Vio'] = Gurobi_Abs_vio
csv_input['Gurobi Rel Vio'] = Gurobi_Rel_vio
csv_input['Scip Objective'] = Scip_Objective
csv_input['Scip Dual Bound'] = Scip_Dual_Bound
csv_input['Scip time taken'] = Scip_Time_taken
csv_input['Scip Con Vio'] = Scip_Con_vio
csv_input['Scip Abs Vio'] = Scip_Abs_vio
csv_input['Scip Rel Vio'] = Scip_Rel_vio
#csv_input['Knitro Objective'] = Knitro_Objective
#csv_input['Knitro time taken'] = Knitro_Time_taken
#csv_input['Ipopt Objective'] = Ipopt_Objective
#csv_input['Ipopt time taken'] = Ipopt_Time_taken
#csv_input['Bonmin Objective'] = Bonmin_Objective
#csv_input['Bonmin time taken'] = Bonmin_Time_taken
#csv_input['Heuristic Objective'] = Heuristic_Objective
#csv_input['Heuristic time taken'] = Heuristic_Time_taken
#csv_input['Heuristic Rel Vio'] = Heuristic_Rel_vio

csv_input.to_csv('output.csv', index=False)

