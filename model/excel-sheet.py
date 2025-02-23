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
    file_data = output.read().split('\n')
    time = 0
    find, time = find_float(file_data, "Total CPU time used: ", time)
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

    # print("model instance:",instance)
    # print("Objective:",Objective)
    # print("Lower bound:",lb1)
    # print("Uper bound:",ub1)
    # print("Time:",time) 
    # print(" ")
    return Objective, time

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
    file_data = output.read().split('\n')
    time = 0
    find, time = find_float(file_data, "solve_time:", time)
    Objective = 0
    find, Objective = find_float(file_data, "total_cost:", Objective)
    return Objective, time

def BonminInstanceOutput(instance, output):
    file_data = output.read().split('\n')
    time = 0
    # find, time = find_float(file_data, "Total CPU secs in IPOPT", time)
    Objective = 0
    find, Objective = find_float(file_data, "total_cost", Objective)
    find, time = find_float(file_data, "solve_time", time)
    
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

def GurobiInstanceOutput1(instance, output):
    obj_pattern = re.search(r'Optimal objective\s+([-\d\.eE]+)', output)
    time_pattern = re.search(r'Solved in\s+([\d\.]+)\s+seconds', output)

    objective = float(obj_pattern.group(1)) if obj_pattern else None
    solver_time = float(time_pattern.group(1)) if time_pattern else None

    return objective, solver_time

def GurobiInstanceOutput(instance, output):
    # Match best objective using a robust regex pattern
    #obj_match = re.search(r'Best objective\s+([-\d\.eE]+)', output)
    obj_match = re.search(r'Best objective\s+([-]?\d+\.\d+e[+-]?\d+)', output)

    # Alternative objective format from "optimal solution" line
    obj_alt_match = re.search(r'optimal solution; objective\s+([-\d\.eE]+)', output, re.IGNORECASE)

    # Extract objective safely
    objective = None
    extracted_obj = obj_match.group(1) if obj_match else (obj_alt_match.group(1) if obj_alt_match else None)

    no_feasible_solution = "Solution count 0" in output
    
    if no_feasible_solution:
        objective = 'no feasible solution'
    if extracted_obj:
        try:
            # Ensure no truncation and remove any extra spaces or unwanted characters
            extracted_obj = extracted_obj.strip().replace(",", "")
            objective = float(extracted_obj)  # Convert to float
        except ValueError:
            print(f"Warning: Could not convert extracted objective '{extracted_obj}' to float.")

    # Match solver time (ensure correct extraction)
    time_match = re.search(r'solve_time:\s*([\d\.]+)', output)

    time_limit_reached = "Time limit reached" in output

    solver_time = None
    if time_limit_reached:
        solver_time = 3600
    elif time_match:
        try:
            solver_time = float(time_match.group(1).strip())
        except ValueError:
            print(f"Warning: Could not convert solver time '{time_match.group(1)}' to float.")

    return objective, solver_time

def SCIPInstanceOutput(instance, output):
    # Match best objective using a robust regex pattern
    #obj_match = re.search(r'Best objective\s+([-\d\.eE]+)', output)
    obj_match = re.search(r"Primal Bound\s*:\s*([+-]?\d+\.\d+e[+-]?\d+)", output)
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
        except ValueError:
            print(f"Warning: Could not convert extracted objective '{extracted_obj}' to float.")

    # Match solver time (ensure correct extraction)
    #time_match = re.search(r"Solving Time \(sec\)\s*:\s*([\d\.]+)", output)
    #time_match = re.search(r'Solving Time (sec) :\s*([\d\.]+)', output)

    time_match = re.search(r'solve_time:\s*([\d\.]+)', output)
    
    time_limit_reached = "time limit reached" in output

    solver_time = None

    if time_match:
        try:
            solver_time = float(time_match.group(1).strip())
        except ValueError:
            print(f"Warning: Could not convert solver time '{time_match.group(1)}' to float.")
    elif time_limit_reached:
        solver_time = 3600
    
    return objective, solver_time




def HeuristicInstanceOutput(instance, output):
    file_data = output.read().split('\n')
    time = 0
    # find, time = find_float(file_data, "Total CPU secs in IPOPT", time)
    Objective = 0
    find, Objective = find_float(file_data, "Final best objective:", Objective)
    find, time = find_float(file_data, "Heuristic elapsed time:", time)
    return Objective, time

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
Knitro_Objective = []
Knitro_Time_taken = []
Ipopt_Objective = []
Ipopt_Time_taken = []
Bonmin_Objective = []
Bonmin_Time_taken = []
Gurobi_Objective = []
Gurobi_Time_taken = []
Scip_Objective = []
Scip_Time_taken = []
Heuristic_Objective = []
Heuristic_Time_taken = []

approx_folder = "foa1"

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

print("******************************Results of Baron Solver ************************************")

for ins in data_list:
    with open(f"../output/baron_out/{approx_folder}/{ins}.baron_out") as output:
        print("Model Name:",ins)
        obj, time= BaronInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Baron_Objective.append(obj)
        Baron_Time_taken.append(time)

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

print("**********************Results of IPOPT Solver *********************************")

for ins in data_list:
    with open(f"../output/ipopt_out/{ins}.ipopt_out") as output:
        print("Model Name:",ins)
        obj, time = IpoptInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Ipopt_Objective.append(obj)
        Ipopt_Time_taken.append(time)

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
        obj, time = GurobiInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Gurobi_Objective.append(obj)
        Gurobi_Time_taken.append(time)

print("**********************Results of SCIP Solver *********************************")
        
for ins in data_list:
    with open(f"../output/scip_out/{approx_folder}/{ins}.scip_out") as output:
        print("Model Name:",ins)
        output = output.read()
        obj, time = SCIPInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Scip_Objective.append(obj)
        Scip_Time_taken.append(time)


print("**********************Results of Heuristic *********************************")

for ins in data_list:
    with open(f"../output/heuristic_out/{ins}.heuristic_out") as output:
        print("Model Name:",ins)
        obj, time = HeuristicInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Heuristic_Objective.append(obj)
        Heuristic_Time_taken.append(time)

fields = ["Instances","Mmultistart Objective","Mmultistart time taken","Baron Objective","Baron time taken", "Gurobi Objective","Gurobi time taken", "Scip Objective","Scip time taken""Knitro Objective","Knitro time taken","Ipopt Objective","Ipopt time taken", "Bonmin Objective","Bonmin time taken", "Heuristic Objective","Heuristic time taken" ]

filename = "mmultistart_results.csv"

with open(filename, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    # csvwriter.writerow(Solver_name)
    csvwriter.writerow(fields)

import pandas as pd
csv_input = pd.read_csv(filename)
csv_input['Instances'] = Ins
csv_input['Mmultistart Objective'] = Mmultistart_Objective
csv_input['Mmultistart time taken'] = Mmultistart_Time_taken
csv_input['Baron Objective'] = Baron_Objective
csv_input['Baron time taken'] = Baron_Time_taken
csv_input['Gurobi Objective'] = Gurobi_Objective
csv_input['Gurobi time taken'] = Gurobi_Time_taken
csv_input['Scip Objective'] = Scip_Objective
csv_input['Scip time taken'] = Scip_Time_taken
csv_input['Knitro Objective'] = Knitro_Objective
csv_input['Knitro time taken'] = Knitro_Time_taken
csv_input['Ipopt Objective'] = Ipopt_Objective
csv_input['Ipopt time taken'] = Ipopt_Time_taken
csv_input['Bonmin Objective'] = Bonmin_Objective
csv_input['Bonmin time taken'] = Bonmin_Time_taken
csv_input['Heuristic Objective'] = Heuristic_Objective
csv_input['Heuristic time taken'] = Heuristic_Time_taken

csv_input.to_csv('output.csv', index=False)

from googleapiclient.http import MediaFileUpload
from googleapiclient.discovery import build
from google.oauth2 import service_account

# Load credentials from the JSON key file
SERVICE_ACCOUNT_FILE = "credentials.json"  # Ensure this file is in the same directory

SCOPES = ["https://drive.google.com/file/d/1nDIsrSl_JhvciyTIt06J9A71UY1pRlIU/view?usp=sharing"]
creds = service_account.Credentials.from_service_account_file(SERVICE_ACCOUNT_FILE, scopes=SCOPES)

# Create the Drive API service
drive_service = build("drive", "v3", credentials=creds)

# Upload the file
file_metadata = {"name": "output.csv"}  # Change name if needed
media = MediaFileUpload("output.csv", mimetype="text/csv")

file = drive_service.files().create(body=file_metadata, media_body=media, fields="id").execute()
print(f"File uploaded successfully! File ID: {file.get('id')}")

