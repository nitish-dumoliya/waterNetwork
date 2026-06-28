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


def MmultistartInstanceOutput(instance,output):
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


def find_best_solution_time(log_text, best_known=None, tol=1e-4):
    """
    Parse a BARON log and find the time when the best known solution
    (upper bound) was first achieved.

    Parameters
    ----------
    log_text   : str   - full text of the BARON log
    best_known : float - known optimal objective value (optional).
                         If None, uses the final reported upper bound.
    tol        : float - relative tolerance for matching the best known value
    """

    # --- 1. Extract the final reported objective (upper bound) ---
    final_ub_match = re.search(
        r"Objective upper bound\s*=\s*([\d.eE+\-]+)", log_text
    )
    final_obj_match = re.search(
        r"(?:Objective|Total cost using \w+):\s*([\d.eE+\-]+)", log_text
    )

    if best_known is None:
        if final_ub_match:
            best_known = float(final_ub_match.group(1))
        elif final_obj_match:
            best_known = float(final_obj_match.group(1))
        else:
            print("Could not determine best known solution from log.")
            return None

    print(f"Target best known solution : {best_known}")
    print(f"Matching tolerance         : {tol} (relative)")
    print("-" * 60)

    # --- 2. Parse iteration table rows ---
    # Format: [*] iteration   time(s)   mem   lower_bound   upper_bound   progress%
    row_pattern = re.compile(
        r"(\*)?\s+(\d+)\s+([\d.]+)\s+\S+\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+%)"
    )

    # Also capture "Preprocessing found feasible solution with value X"
    preprocess_pattern = re.compile(
        r"Preprocessing found feasible solution with value\s+([\d.eE+\-]+)"
    )

    first_time = None
    first_iter = None
    first_ub   = None

    # Check preprocessing line
    pre_match = preprocess_pattern.search(log_text)
    if pre_match:
        pre_val = float(pre_match.group(1))
        rel_diff = abs(pre_val - best_known) / max(abs(best_known), 1e-10)
        if rel_diff <= tol:
            print(f"Best solution found during PREPROCESSING")
            print(f"  Value : {pre_val}")
            print(f"  Time  : before iteration table (< first logged time)")
            first_ub = pre_val

    # Walk through iteration rows
    for m in row_pattern.finditer(log_text):
        starred   = m.group(1) == "*"
        iteration = int(m.group(2))
        time_s    = float(m.group(3))
        lb        = float(m.group(4))
        ub        = float(m.group(5))

        rel_diff = abs(ub - best_known) / max(abs(best_known), 1e-10)

        if rel_diff <= tol:
            if first_time is None:
                first_time = time_s
                first_iter = iteration
                first_ub   = ub
                marker = " (*)" if starred else ""
                print(f"Best solution first matched in iteration table:")
                print(f"  Iteration : {first_iter}{marker}")
                print(f"  Time      : {first_time:.2f} s")
                print(f"  UB value  : {first_ub}")
                print(f"  Rel. diff : {rel_diff:.2e}")
            break   # stop at first match

    if first_time is None and first_ub is None:
        print("Best known solution not found within tolerance in the log.")
        first_time = 600
    elif first_time is None and first_ub is not None:
        # matched preprocessing but no row — try to grab earliest table time
        first_row = row_pattern.search(log_text)
        if first_row:
            print(f"  (First iteration table time: {float(first_row.group(3)):.2f} s)")

    return first_time

# ===========================================================================
# GUROBI LOG PARSER
# ===========================================================================

def find_best_solution_time_gurobi(log_text, best_known=None, tol=1e-4):
    """
    Parse a Gurobi log and find the time when the best known solution
    (incumbent) was first achieved.

    Gurobi marks new incumbents with H or * in the node table:
      H <nodes> <unexpl> ... <incumbent> <bestbd> <gap> ... <time>
      * <nodes> <unexpl> ... <incumbent> <bestbd> <gap> ... <time>

    Parameters
    ----------
    log_text   : str   - full text of the Gurobi log
    best_known : float - known optimal objective value (optional).
                         If None, uses the final reported best objective.
    tol        : float - relative tolerance for matching the best known value
    """

    # --- 1. Determine best_known from log if not provided ---
    if best_known is None:
        # "Best objective 1.238826143764e+05, best bound ..."
        final_match = re.search(
            r"Best objective\s+([\d.eE+\-]+)", log_text
        )
        # fallback: "Total cost using gurobi: X"
        fallback = re.search(
            r"Total cost using \w+:\s*([\d.eE+\-]+)", log_text
        )
        if final_match:
            best_known = float(final_match.group(1))
        elif fallback:
            best_known = float(fallback.group(1))
        else:
            print("Could not determine best known solution from Gurobi log.")
            return None

    print(f"Target best known solution : {best_known}")
    print(f"Matching tolerance         : {tol} (relative)")
    print("-" * 60)

    # --- 2. Parse incumbent rows from the node table ---
    # Lines starting with H or * indicate a new incumbent:
    # H  <nodes>  <unexpl>  <obj_depth_intinf_block>  <incumbent>  <bestbd>  <gap>  <it/node>  <time>
    # The incumbent value and time are the key fields.
    #
    # Gurobi node table row format (H/* lines):
    # [H|*] nodes  unexpl | obj depth intinf | incumbent  bestbd  gap | it/node  time
    # Time field ends with 's'
    #
    # We capture: marker, incumbent, time
    # Simpler pattern: incumbent is first large float, time is last field ending in s
    incumbent_pattern = re.compile(
        r"^[H\*][\s\d]+?"          # H or * then nodes/unexpl (non-greedy)
        r"([\d]{4,}\.[\d]+)\s+"   # incumbent (>=4 digits before decimal)
        r"[\d.,eE+\-]+\s+"         # bestbd
        r"[\d.]+%\s+"               # gap
        r"[\d.\-]+\s+"             # it/node
        r"([\d.]+)s",                # time
        re.MULTILINE
    )

    first_time = None
    first_incumbent = None

    for m in incumbent_pattern.finditer(log_text):
        incumbent = float(m.group(1))
        time_s    = float(m.group(2))

        rel_diff = abs(incumbent - best_known) / max(abs(best_known), 1e-10)

        if rel_diff <= tol:
            if first_time is None:
                first_time       = time_s
                first_incumbent  = incumbent
                print(f"Best solution first matched in node table:")
                print(f"  Incumbent : {first_incumbent}")
                print(f"  Time      : {first_time:.2f} s")
                print(f"  Rel. diff : {rel_diff:.2e}")
            break

    if first_time is None:
        print("Best known solution not found within tolerance in Gurobi log.")
        first_time = 600

    return first_time

# ===========================================================================
# SCIP LOG PARSER
# ===========================================================================
 
def find_best_solution_time_scip(log_text, best_known=None, tol=1e-4):
    """
    Parse a SCIP log and find the time when the best known solution
    (primal bound) was first achieved.
 
    SCIP node table format (pipe-separated):
      [marker] time | node | left | LP iter | LP it/n | mem/heur | mdpt |
      vars | cons | rows | cuts | sepa | confs | strbr | dualbound |
      primalbound | gap | compl.
 
    New incumbents are marked with L (heuristic), * (LP), o (feasible), etc.
    Rows where primalbound = '--' are skipped (no feasible solution yet).
 
    Parameters
    ----------
    log_text   : str   - full text of the SCIP log
    best_known : float - known optimal objective value (optional).
                         If None, uses the final Primal Bound from log.
    tol        : float - relative tolerance for matching the best known value
    """
 
    # --- 1. Determine best_known from log if not provided ---
    if best_known is None:
        final_match = re.search(
            r"Primal Bound\s*:\s*[+\-]?([\d.eE+\-]+)", log_text
        )
        fallback = re.search(
            r"Total cost using \w+:\s*([\d.eE+\-]+)", log_text
        )
        if final_match:
            best_known = float(final_match.group(1))
        elif fallback:
            best_known = float(fallback.group(1))
        else:
            print("Could not determine best known solution from SCIP log.")
            return None
 
    print(f"Target best known solution : {best_known}")
    print(f"Matching tolerance         : {tol} (relative)")
    print("-" * 60)
 
    # --- 2. Walk the node table line by line, split on '|' ---
    first_time = None
    first_pb   = None
 
    for line in log_text.splitlines():
        # Must contain at least several pipe characters (node table rows)
        if line.count('|') < 10:
            continue
 
        parts = [p.strip() for p in line.split('|')]
        if len(parts) < 4:
            continue
 
        # primalbound is 3rd field from the end: | primalbound | gap | compl.
        pb_str = parts[-3].strip()
        if pb_str == '--' or not pb_str:
            continue      # no feasible solution yet on this row
 
        try:
            pb = float(pb_str)
        except ValueError:
            continue
 
        # Time is embedded in the first field, e.g. "L 0.1s" or "  0.5s"
        time_match = re.search(r'([\d.]+)s', parts[0])
        if not time_match:
            continue
        t = float(time_match.group(1))
 
        rel_diff = abs(pb - best_known) / max(abs(best_known), 1e-10)
 
        if rel_diff <= tol:
            first_time = t
            first_pb   = pb
            # Detect marker (L, *, o, space = regular)
            marker_match = re.match(r'^([A-Za-z\*o])', line.strip())
            marker = marker_match.group(1) if marker_match else ' '
            print(f"Best solution first matched in node table:")
            print(f"  Marker    : '{marker}' "
                  f"({'heuristic' if marker=='L' else 'LP/branch' if marker=='*' else 'feasible' if marker=='o' else 'regular'})")
            print(f"  Time      : {first_time:.2f} s")
            print(f"  Primal bd : {first_pb}")
            print(f"  Rel. diff : {rel_diff:.2e}")
            break
 
    if first_time is None:
        print("Best known solution not found within tolerance in SCIP log.")
        first_time=600
 
    return first_time

# def GurobiInstanceOutput(output_content, best_known=None, tol=1e-4):
#     """
#     Wrapper matching the style of BaronInstanceOutput.
#     Pass output_content (string), not a file object.
#     Returns (Objective, time, first_solution_time)
#     """
#     # Extract final objective
#     obj_match = re.search(
#         r"Total cost using \w+:\s*([\d.eE+\-]+)", output_content
#     )
#     Objective = float(obj_match.group(1)) if obj_match else None
#
#     # Extract solve time
#     time_match = re.search(r"gurobi solve time:\s*([\d.eE+\-]+)", output_content)
#     solve_time = float(time_match.group(1)) if time_match else None
#
#     if best_known is None:
#         best_known = Objective
#
#     first_time = find_best_solution_time_gurobi(output_content, best_known=best_known, tol=tol)
#
#     return Objective, solve_time, first_time



def BaronInstanceOutput(instance, output):
    
    # Extract total relative constraint violation
    output_content = output.read() 

    # violation_match_con = re.search(r'Total absolute constraints violation: \s*([\d.eE+-]+)', output_content)
    # violation_match_abs = re.search(r'Sum of absolute diff (orig vs approx): \s*([\d.eE+-]+)', output_content)
    # violation_match_rel = re.search(r'Sum of relative diff (orig vs approx): \s*([\d.eE+-]+)', output_content)
    violation_match_con = re.search(r'Total absolute constraint violation:\s*([\d.eE+-]+)',output_content)
    
    violation_match_orig = re.search(r'Sum of original head-loss violation\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_approx = re.search(r'Sum of approx\s+head-loss violation\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_abs = re.search(r'Sum of absolute diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_rel = re.search(r'Sum of relative diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output_content)
    # Regular expression to find the lower bound
    lower_bound = re.search(r"lower bound\s*=\s*([0-9.]+)", output_content)
    
    if lower_bound:
        lower_bound = float(lower_bound.group(1))
        print("Lower bound:", lower_bound)
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

    first_time = find_best_solution_time(output_content, best_known=Objective, tol=1e-4)

    return Objective, lower_bound, first_time, violation_con, violation_abs, violation_rel

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
        solver_time = 600
    elif time_match:
        try:
            solver_time = float(time_match.group(1).strip())
        except ValueError:
            print(f"Warning: Could not convert solver time '{time_match.group(1)}' to float.")
   
    # Extract total relative constraint violation
    # violation_match_con = re.search(r'Sum of constraints violation:\s*([\d.eE+-]+)', output)
    # violation_match_abs = re.search(r'Con2 sum of absolute violation between original function and approximate function:\s*([\d.eE+-]+)', output)
    # violation_match_rel = re.search(r'Con2 sum of relative violation between original function and approximate function:\s*([\d.eE+-]+)', output)
    #print(violation_match)
    violation_match_con = re.search(r'Total absolute constraint violation:\s*([\d.eE+-]+)',output)
    
    violation_match_orig = re.search(r'Sum of original head-loss violation\s*:\s*([\d.eE+-]+)',output)
    
    violation_match_approx = re.search(r'Sum of approx\s+head-loss violation\s*:\s*([\d.eE+-]+)',output)
    
    violation_match_abs = re.search(r'Sum of absolute diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output)
    
    violation_match_rel = re.search(r'Sum of relative diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output)   

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
  
    # obj, solve_time, first_time = GurobiInstanceOutput(gurobi_log, best_known=best_known, tol=1e-4)

    first_time = find_best_solution_time_gurobi(output, best_known=objective, tol=1e-3)

    return objective, best_bound, first_time, violation_con, violation_abs, violation_rel

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
        except ValueError:
            print(f"Warning: Could not convert extracted objective '{extracted_obj}' to float.")
    if dual_bound:
        try:
            # Ensure no truncation and remove any extra spaces or unwanted characters
            dual_bound = float(dual_bound.group(1))
        except ValueError:
            print(f"Warning: Could not convert extracted dual bound '{extracted_obj}' to float.")


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
    # violation_match_con = re.search(r'Sum of constraints violation:\s*([\d.eE+-]+)', output)
    # violation_match_abs = re.search(r'Con2 sum of absolute violation between original function and approximate function:\s*([\d.eE+-]+)', output)
    # violation_match_rel = re.search(r'Con2 sum of relative violation between original function and approximate function:\s*([\d.eE+-]+)', output)
    #print(violation_match)
    violation_match_con = re.search(r'Total absolute constraint violation:\s*([\d.eE+-]+)',output)
    
    violation_match_orig = re.search(r'Sum of original head-loss violation\s*:\s*([\d.eE+-]+)',output)
    
    violation_match_approx = re.search(r'Sum of approx\s+head-loss violation\s*:\s*([\d.eE+-]+)',output)
    
    violation_match_abs = re.search(r'Sum of absolute diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output)
    
    violation_match_rel = re.search(r'Sum of relative diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output)

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

    first_time = find_best_solution_time_scip(output, best_known=objective, tol=1e-1)

    
    return objective, dual_bound, first_time, violation_con, violation_abs, violation_rel

def KnitroInstanceOutput(instance, output):
    file_data = output.read().split('\n')
    time = 0
    Objective = 0
    find, Objective = find_float(file_data, "Total cost using knitro", Objective)
    find, time = find_float(file_data, "knitro solve time:", time)
    

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
    # violation_match_con = re.search(r'Sum of constraints violation:\s*([\d.eE+-]+)', output_content)
    # violation_match_abs = re.search(r'Con2 sum of absolute violation between original function and approximate function:\s*([\d.eE+-]+)', output_content)
    # violation_match_rel = re.search(r'Con2 sum of relative violation between original function and approximate function:\s*([\d.eE+-]+)', output_content)
    #
    violation_match_con = re.search(r'Total absolute constraint violation:\s*([\d.eE+-]+)',output_content)
    
    violation_match_orig = re.search(r'Sum of original head-loss violation\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_approx = re.search(r'Sum of approx\s+head-loss violation\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_abs = re.search(r'Sum of absolute diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_rel = re.search(r'Sum of relative diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output_content)

    violation_con = float(violation_match_con.group(1)) if violation_match_con else None
    violation_abs = float(violation_match_abs.group(1)) if violation_match_abs else None
    violation_rel = float(violation_match_rel.group(1)) if violation_match_rel else None

    # Extract objective value
    objective_match = re.search(r'Total cost using ipopt:\s*([\d.eE+-]+)', output_content)
    objective = float(objective_match.group(1)) if objective_match else None

    # Extract solve time
    time_match = re.search(r'ipopt solve time:\s*([\d.eE+-]+)', output_content)
    solve_time = float(time_match.group(1)) if time_match else None

    return objective, solve_time, violation_con, violation_abs, violation_rel


def BonminInstanceOutput(instance, output):
    file_data = output.read().split('\n')
    time = 0
    # find, time = find_float(file_data, "Total CPU secs in IPOPT", time)
    Objective = 0
    find, Objective = find_float(file_data, "Total cost using /home/nitishdumoliya/build/bin/bonmin", Objective)
    find, time = find_float(file_data, "/home/nitishdumoliya/build/bin/bonmin solve time:", time)
    
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




def HeuristicInstanceOutput(instance, output):    
    output_content = output.read() 
    # Extract total relative constraint violation
    # violation_match_con = re.search(r'Sum of constraints violation:\s*([\d.eE+-]+)', output_content)
    # violation_match_abs = re.search(r'Con2 sum of absolute violation between original function and approximate function:\s*([\d.eE+-]+)', output_content)
    # violation_match_rel = re.search(r'Con2 sum of relative violation between original function and approximate function:\s*([\d.eE+-]+)', output_content)
    #
    violation_match_con = re.search(r'Total absolute constraint violation:\s*([\d.eE+-]+)',output_content)
    
    violation_match_orig = re.search(r'Sum of original head-loss violation\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_approx = re.search(r'Sum of approx\s+head-loss violation\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_abs = re.search(r'Sum of absolute diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output_content)
    
    violation_match_rel = re.search(r'Sum of relative diff\s*\(orig vs approx\)\s*:\s*([\d.eE+-]+)',output_content)

    violation_con = float(violation_match_con.group(1)) if violation_match_con else None
    violation_abs = float(violation_match_abs.group(1)) if violation_match_abs else None
    violation_rel = float(violation_match_rel.group(1)) if violation_match_rel else None

 
    # Extract total relative constraint violation
    # violation_match = re.search(r'Total relative constraint violation:\s*([\d.eE+-]+)', output_content)
    
    #violation = None
    #if violation_match:
    #    try:
    #        violation = float(violation_match.group(1).strip())  # Convert extracted value to float
    #    except ValueError:
    #        print(f"Warning: Could not convert extracted violation '{violation_match.group(1)}' to float.")

    #file_data = output.read().split('\n')
    file_data = output_content.split('\n')
    time = 0
    # find, time = find_float(file_data, "Total CPU secs in IPOPT", time)
    Objective = 0
    find, Objective = find_float(file_data, "Final best objective:", Objective)
    find, time = find_float(file_data, "Total elapsed time:", time)
    number_of_nlp = 0
    find, number_of_nlp = find_float(file_data, "Number of nlp problem solved:", number_of_nlp)
    number_of_iteration = 0
    find, number_of_iteration = find_float(file_data, "Total number of iteration:", number_of_iteration)
  
    return Objective, time, violation_con,violation_abs, violation_rel, number_of_nlp, number_of_iteration

data_list = [
    "d1_bessa",
    "d2_shamir",
    "d3_hanoi",
    "d4_double_hanoi",
    "d5_triple_hanoi",
    "d6_newyork",
    "d7_blacksburg",
    "d8_fossolo_iron",
    "d9_fossolo_poly_0",
    "d10_fossolo_poly_1",
    "d11_kadu",
    "d12_pescara",
    "d13_modena",
    "d14_balerma"
]
DATA_LIST =  [
        "d1_bessa",
        "d2_shamir",
        "d3_trn",
        "d4_newyork",
        "d5_goyang",
        "d6_kadu",
        "d7_blacksburg",
        "d8_hanoi",
        "d9_bakryun",
        "d10_fossolo_iron",
        "d11_fossolo_poly_0",
        "d12_fossolo_poly_1",
        "d13_farhadgerd",
        "d14_cgn",
        "d15_double_hanoi",
        "d16_pescara",
        "d17_triple_hanoi",
        "d18_zj_network",
        "d19_small",
        "d20_modena",
        "d21_ramnagar",
        "d22_marchi_rural",
        "d23_balerma",
        "d24_kl",
        # "d25_large"
    ]




# !python3 ampl_run.py m1_basic.mod d1_Sample_input_cycle_twoloop.dat > output/d1_Sample_input_cycle_twoloop.baro_out

import pandas as pd
import csv

Ins = []
Baron_Objective = []
Baron_Lower_bound = []
Baron_Time_taken = []
Baron_Con_vio = []
Baron_Abs_error = []
Baron_Rel_error = []
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
Gurobi_Abs_error = []
Gurobi_Rel_error = []
Scip_Objective = []
Scip_Dual_Bound = []
Scip_Time_taken = []
Scip_Con_vio = []
Scip_Abs_error = []
Scip_Rel_error = []
InitialIpopt_Objective = []
InitialIpopt_Time_taken = []
InitialIpopt_Con_vio = []
InitialIpopt_Abs_error = []
InitialIpopt_Rel_error = []

Knitro_Objective = []
Knitro_Time_taken = []
Bonmin_Objective = []
Bonmin_Time_taken = []
Mmultistart_Objective = []
Mmultistart_Time_taken = []
Heuristic_Objective = []
Heuristic_Time_taken = []
Heuristic_Con_vio = []
Heuristic_Abs_error = []
Heuristic_Rel_error = []
Heuristic_Number_of_NLP = []
Heuristic_Number_of_Iteration = []

dir0 = "opt_eng_output2/model1"
#dir0 = "europt"
dir1 = "original_10min"
dir = "smooth-approx"
dir2 = "approx2"
# error = "absolute"
error = "relative"

print("******************************Results of Initial Ipopt Solve Result ************************************")
for ins in DATA_LIST:
    with open(f"../{dir0}/ipopt_out/{dir2}/{error}/{ins}.ipopt_out") as output:
        print("Model Name:",ins)
        out = output
        obj, time, violation_con,violation_abs, violation_rel = IpoptInstanceOutput(ins,out)
        print("Objective :",obj)
        print("Time :",time)
        print("Constraint violation :",violation_rel)
        print("Absolute violation :",violation_abs)
        print("Relative violation :",violation_rel)

        print(" ")
        Ins.append(ins)
        InitialIpopt_Objective.append(obj)
        InitialIpopt_Time_taken.append(time)
        InitialIpopt_Con_vio.append(violation_con)
        InitialIpopt_Abs_error.append(violation_abs)
        InitialIpopt_Rel_error.append(violation_rel)

print("******************************Results of Baron Solver ************************************")

for ins in DATA_LIST:
    with open(f"../{dir0}/baron_out/{dir1}/{ins}.baron_out") as output:
        print("Model Name:",ins)
        out = output
        obj, lower_bound, time, violation_con,violation_abs, violation_rel= BaronInstanceOutput(ins,out)
        print("Objective :",obj)
        print("Lower bound :",lower_bound)
        print("Time :",time)
        print("Constraint violation :",violation_rel)
        print("Absolute violation :",violation_abs)
        print("Relative violation :",violation_rel)

        print(" ")
        Baron_Objective.append(obj)
        Baron_Lower_bound.append(lower_bound)
        Baron_Time_taken.append(time)
        Baron_Con_vio.append(violation_con)
        Baron_Abs_error.append(violation_abs)
        Baron_Rel_error.append(violation_rel)

print("**********************Results of GUROBI Solver *********************************")

for ins in DATA_LIST:
    with open(f"../{dir0}/gurobi_out/{dir1}/{ins}.gurobi_out") as output:
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
        Gurobi_Abs_error.append(violation_abs)
        Gurobi_Rel_error.append(violation_rel)

print("**********************Results of SCIP Solver *********************************")

for ins in DATA_LIST:
    with open(f"../{dir0}/scip_out/{dir1}/{ins}.scip_out") as output:
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
        Scip_Abs_error.append(violation_abs)
        Scip_Rel_error.append(violation_rel)

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


print("**********************Results of Knitro Solver *********************************")
 
for ins in DATA_LIST:
    with open(f"../{dir0}/knitro_out/{dir1}/{ins}.knitro_out") as output:
        print("Model Name:",ins)
        obj, time = KnitroInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Knitro_Objective.append(obj)
        Knitro_Time_taken.append(time)

print("**********************Results of BONMIN Solver *********************************")

for ins in DATA_LIST:
    with open(f"../{dir0}/bonmin_out/{dir2}/{error}/{ins}.bonmin_out") as output:
        print("Model Name:",ins)
        obj, time = BonminInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print(" ")
        Bonmin_Objective.append(obj)
        Bonmin_Time_taken.append(time)

print("****************************Results of Mmultistart Solver *********************************")

# for ins in data_list:
#     with open(f"../europt/mmultistart_out/original/{ins}.mmultistart_out") as output:
#         print("Model Name:",ins)
#         obj, time, ub = MmultistartInstanceOutput(ins,output)
        # print("Objective :",obj)
        # print("Time :",time)
        # print(" ")
        # Mmultistart_Objective.append(obj)
        # Mmultistart_Time_taken.append(time)

print("**********************Results of Heuristic *********************************")

for ins in DATA_LIST:
    with open(f"../opt_eng_output/model2/heuristic_out/{dir2}/{error}/{ins}.heuristic_out") as output:
        print("Model Name:",ins)
        obj, time, violation_con,violation_abs, violation_rel, nlp, iter = HeuristicInstanceOutput(ins,output)
        print("Objective :",obj)
        print("Time :",time)
        print("Constraint violation :",violation_con)
        print("Absolute violation :",violation_abs)
        print("Relative violation :",violation_rel)

        print(" ")
        Heuristic_Con_vio.append(violation_con)
        Heuristic_Abs_error.append(violation_abs)
        Heuristic_Rel_error.append(violation_rel)

        #print("Relative violation :",violation)
        print(" ")
        Heuristic_Objective.append(obj)
        Heuristic_Time_taken.append(time)
        Heuristic_Number_of_NLP.append(nlp)
        Heuristic_Number_of_Iteration.append(iter)

fields = ["Instances","Ini Ipopt Objective","Ini Ipopt time taken","Ini Ipopt Con Vio","Ini Ipopt Abs Error","Ini Ipopt Rel Error","Baron Lower Bound", "Baron Objective","Baron Con Vio","Baron Abs Error","Baron Rel Error","Baron time taken", "Gurobi Best Bound","Gurobi Objective", "Gurobi Con Vio","Gurobi Abs Error","Gurobi Rel Error","Gurobi time taken","Scip Dual Bound","Scip Objective",  "Scip Con Vio", "Scip Abs Error", "Scip Rel Error","Scip time taken","Knitro Objective","Knitro time taken", "Bonmin Objective","Bonmin time taken", "Mmultistart Objective","Mmultistart time taken","Heuristic Objective","Heuristic time taken", "Heuristic Con Vio", "Heuristic Abs Error", "Heuristic Rel Error", "Heuristic NLP", "Heuristic Iteration"  ]

filename = "mmultistart_results.csv"

with open(filename, 'w') as csvfile:
   csvwriter = csv.writer(csvfile)
   # csvwriter.writerow(Solver_name)
   csvwriter.writerow(fields)

import pandas as pd
csv_input = pd.read_csv(filename)
csv_input['Instances'] = Ins
# csv_input['Mmultistart Objective'] = Mmultistart_Objective
# csv_input['Mmultistart time taken'] = Mmultistart_Time_taken
csv_input['Ini Ipopt Objective'] = InitialIpopt_Objective
csv_input['Ini Ipopt time taken'] = InitialIpopt_Time_taken
csv_input['Ini Ipopt Con Vio'] = InitialIpopt_Con_vio
csv_input['Ini Ipopt Abs Error'] = InitialIpopt_Abs_error
csv_input['Ini Ipopt Rel Error'] = InitialIpopt_Rel_error
csv_input['Baron Lower Bound'] = Baron_Lower_bound
csv_input['Baron Objective'] = Baron_Objective
csv_input['Baron time taken'] = Baron_Time_taken
csv_input['Baron Con Vio'] = Baron_Con_vio
csv_input['Baron Abs Error'] = Baron_Abs_error
csv_input['Baron Rel Error'] = Baron_Rel_error
csv_input['Gurobi Objective'] = Gurobi_Objective
csv_input['Gurobi Best Bound'] = Gurobi_Best_Bound
csv_input['Gurobi time taken'] = Gurobi_Time_taken
csv_input['Gurobi Con Vio'] = Gurobi_Con_vio
csv_input['Gurobi Abs Error'] = Gurobi_Abs_error
csv_input['Gurobi Rel Error'] = Gurobi_Rel_error
csv_input['Scip Objective'] = Scip_Objective
csv_input['Scip Dual Bound'] = Scip_Dual_Bound
csv_input['Scip time taken'] = Scip_Time_taken
csv_input['Scip Con Vio'] = Scip_Con_vio
csv_input['Scip Abs Error'] = Scip_Abs_error
csv_input['Scip Rel Error'] = Scip_Rel_error
csv_input['Knitro Objective'] = Knitro_Objective
csv_input['Knitro time taken'] = Knitro_Time_taken
#csv_input['Ipopt Objective'] = Ipopt_Objective
#csv_input['Ipopt time taken'] = Ipopt_Time_taken
csv_input['Bonmin Objective'] = Bonmin_Objective
csv_input['Bonmin time taken'] = Bonmin_Time_taken
csv_input['Heuristic Objective'] = Heuristic_Objective
csv_input['Heuristic time taken'] = Heuristic_Time_taken
csv_input['Heuristic Con Vio'] = Heuristic_Con_vio
csv_input['Heuristic Abs Error'] = Heuristic_Abs_error
csv_input['Heuristic Rel Error'] = Heuristic_Rel_error
csv_input['Heuristic NLP'] = Heuristic_Number_of_NLP
csv_input['Heuristic Iteration'] = Heuristic_Number_of_Iteration

csv_input.to_csv('output.csv', index=False)

