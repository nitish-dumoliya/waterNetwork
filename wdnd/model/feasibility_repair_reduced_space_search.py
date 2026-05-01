import time
import sys
import random
import math
import copy
import numpy as np
from amplpy import AMPL
import contextlib
import os
import networkx as nx
import matplotlib.pyplot as plt
from scipy.optimize import linprog

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

class WaterNetworkSolver:
    def __init__(self, model_file, solver_name, data_file, data_number):
        self.model_file = model_file
        self.solver_name = solver_name
        self.data_file = data_file
        self.data_number = data_number

        self.ampl = AMPL()
        # self.ampl.setOption('hsllib', '/usr/local/lib/libma57.so')

        self.q = {}
        self.h = {}
        self.l = {}
        self.y = {}

        self.best_q = {}
        self.best_h = {}
        self.best_l = {}
        self.best_y = {}
        self.best_cost = float("inf")
        self.best_cost = None

        self.network_graph = None
        self.solve_result = None
        self.solver_time = 0
        self.best_acyclic_flow = None
        self.number_of_nlp = 0

        # IPOPT options
        self.max_iter = 3000
        self.mu_init = 1e-2
        self.bound_push = 0.1
        self.bound_frac = 0.1

    # ============================================================
    # READ ORIGINAL MODEL
    # ============================================================
    def read_model_and_data(self):
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)

        self.nodes = list(self.ampl.getSet('nodes'))
        self.source = list(self.ampl.getSet('Source'))
        self.arcs = list(self.ampl.getSet('arcs'))
        self.pipes = list(self.ampl.getSet('pipes'))
        self.segs = list(self.ampl.getSet('segs'))

        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
        self.eps = self.ampl.getParameter('eps').to_dict()
        self.omega = float(self.ampl.getParameter('omega').value())
        self.Q_max = float(self.ampl.getParameter('Q_max').value())
        self.alpha = self.ampl.get_parameter("alpha").get_values().to_dict()
        self.slope = self.ampl.getParameter('slope').to_dict()
        self.intercept = self.ampl.getParameter('intercept').to_dict()

        # y bounds
        self.R_min = min(self.R[k] for k in self.pipes)
        self.R_max = max(self.R[k] for k in self.pipes)
        self.d_min = min(self.d[k] for k in self.pipes)
        self.d_max = max(self.d[k] for k in self.pipes)

        self.y_lb = {}
        self.y_ub = {}
        for (i, j) in self.arcs:
            self.y_lb[(i, j)] = self.omega * self.L[(i, j)] / (self.R_max**1.852 * self.d_max**4.87)
            self.y_ub[(i, j)] = self.omega * self.L[(i, j)] / (self.R_min**1.852 * self.d_min**4.87)


    def fix_leaf_arc_flow(self, ampl):
        graph = nx.Graph()
        arc_set = self.ampl.getSet('arcs').to_list()  
        graph.add_edges_from(arc_set)
        D = self.ampl.getParameter('D').getValues().to_dict()  
        source = self.ampl.getSet('Source').to_list()[0]
        fixed_arcs = set()
        Qmax = self.ampl.getParameter('Q_max').getValues().to_list()[0] 
        D[source] = -Qmax
        # print("\nPresolve the model for fixing the flow value in the leaf arcs")
        # print("Source:",self.ampl.getSet('Source').to_list())

        while True:
            leaf_nodes = [node for node in graph.nodes if graph.degree[node] == 1]
            # print("leaf_nodes:", leaf_nodes)
            if not leaf_nodes:  
                break

            for leaf in leaf_nodes:
                neighbor = next(graph.neighbors(leaf))
                if (neighbor, leaf) in arc_set:
                    edge = (neighbor, leaf)
                    if edge not in fixed_arcs:  
                        if leaf == source:
                            flow_value = D[leaf]
                            D[neighbor] = (D[neighbor]+flow_value)
                            source = neighbor
                        else:
                            flow_value = D[leaf]
                            D[neighbor] = D[neighbor] + flow_value
                        # ampl.eval(f"s.t. fix_q_{edge[0]}_{edge[1]}: q[{edge[0]},{edge[1]}] = {flow_value};")
                        # print(f"Fixing flow for arc {edge}: {flow_value}")
                        fixed_arcs.add(edge)  

                    graph.remove_node(leaf)
                else:
                    edge = (leaf, neighbor)
                    if edge not in fixed_arcs:  
                        if leaf == source:
                            flow_value = -D[leaf]
                            D[neighbor] = D[neighbor]-flow_value
                            source = neighbor
                            # ampl.eval(f"s.t. fix_q_{edge[0]}_{edge[1]}: q[{edge[0]},{edge[1]}] = {flow_value};")
                        elif neighbor == source:
                            flow_value = -D[leaf]
                            D[neighbor] = D[neighbor] - D[leaf] 
                            # ampl.eval(f"s.t. fix_q_{edge[0]}_{edge[1]}: q[{edge[0]},{edge[1]}] = {flow_value};")
                        else:
                            flow_value = -D[leaf]
                            D[neighbor] += -flow_value
                            # ampl.eval(f"s.t. fix_q_{edge[0]}_{edge[1]}: q[{edge[0]},{edge[1]}] = {flow_value};")
                        # print(f"Fixing flow for arc {edge}: {flow_value}")
                        fixed_arcs.add(edge)  
                    graph.remove_node(leaf)
        # print("All leaf arc flows have been fixed.")
        return fixed_arcs

    def is_cycle(self, graph, start_node, end_node, visited_copy, parent):
        visited_copy[start_node] = True
        # print(f"Is node {end_node} in cycle?")
        for neighbor in graph.neighbors(start_node):
            # print("visited",neighbor,visited_copy[neighbor])
            if not visited_copy[neighbor]:
                # print("neighbor of node", start_node, "is", neighbor)
                isCycle = self.is_cycle(graph, neighbor, end_node, visited_copy, start_node)
                if isCycle:
                    return True
            else:
                # print("parent:", parent)
                if parent != neighbor:
                    if end_node == neighbor:
                        # print(f"Node {end_node} is in cycle")
                        return True
        return False

    def presolve(self, graph, node, visited, parent, set_arc):
        visited_copy = visited.copy()
        # print(visited_copy)
        isCycle = self.is_cycle(graph, node, node, visited_copy, parent)
        # print(f"Is node {node} in cycle?",isCycle)
        visited[node] = True
        if isCycle:
            for neighbor in graph.neighbors(node):
                if parent!=neighbor:
                    set_arc.append((node,neighbor))
                    # print("Fix the arc", (node, neighbor))
            return set_arc
        else:
            for neighbor in graph.neighbors(node):
                if parent != neighbor:
                    set_arc.append((node,neighbor))
                    # print(set_arc)
                    # print("Fix the arc", (node, neighbor))
                    # print("neighbor:", neighbor)
                    self.presolve(graph, neighbor, visited, node, set_arc)
        return set_arc

    def fix_arc_set(self):
        graph = nx.Graph()
        arc_set = self.ampl.getSet('arcs').to_list()
        graph.add_edges_from(arc_set)
        visited = {node: False for node in graph.nodes()}
        source = self.ampl.getSet('Source').to_list()[0]
        set_arc = []
        # print("\nPresolve the model for fixing the arc direction")
        set_ = self.presolve(graph, source, visited, -1, set_arc)
        # print("fixed arc direction:",set_, "\n") 
        return set_



    # ============================================================
    # ORIGINAL NLP SOLVE
    # ============================================================
    def solve_original_model_without_init(self):
        # print("\n-------------------------------- Solving Original Model --------------------------")

        if self.data_number == 6:
            self.ampl.eval("subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
        else:
            self.ampl.eval("subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
        self.ampl.option['solver'] = 'ipopt'
        self.ampl.option["ipopt_options"] = (
            "outlev=0 "
            "expect_infeasible_problem=no "
            "bound_relax_factor=0 "
            "tol=1e-9 "
            "bound_push=0.01 "
            "bound_frac=0.01 "
            "warm_start_init_point=no "
            "halt_on_ampl_error=yes"
        )
        self.ampl.option["baron_options"]= "maxtime = 10  outlev = 2 barstats version objbound" # lsolver = conopt

        # with self.suppress_output():
        self.ampl.solve()

        self.q = self.ampl.get_variable('q').get_values().to_dict()
        self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.l = self.ampl.get_variable('l').get_values().to_dict()

        self.best_q = copy.deepcopy(self.q)
        self.best_h = copy.deepcopy(self.h)
        self.best_l = copy.deepcopy(self.l)

        self.best_y = self.compute_y_from_l(self.best_l)
        self.best_cost = self.ampl.getObjective("total_cost").value()

        solve_time = self.ampl.get_value('_solve_elapsed_time')

        print(f"Initial local solution cost: {self.best_cost:.8f}")
        print(f"Original solve time: {solve_time:.4f} sec")

    # ============================================================
    # COMPUTE y FROM l
    # ============================================================
    def compute_y_from_l(self, l_dict):
        y = {}
        for (i, j) in self.arcs:
            y[(i, j)] = sum(
                self.omega * l_dict[(i, j, k)] / (self.R[k]**1.852 * self.d[k]**4.87)
                for k in self.pipes
            )
        return y

    # ============================================================
    # REDUCED NLP: solve q,h for fixed y
    # ============================================================
    def solve_reduced_nlp(self, y_input, q_init=None, h_init=None, model_name="reduced_nlp.mod"):
        ampl = AMPL()
        # ampl.setOption('hsllib', '/usr/local/lib/libma57.so')
        ampl.read(model_name)
        ampl.read_data(self.data_file)

        for (i, j), val in y_input.items():
            ampl.param["y"][i, j] = val

        # Warm start
        if q_init is not None:
            for (i, j), val in q_init.items():
                ampl.var["q"][i, j] = val
        if h_init is not None:
            for i, val in h_init.items():
                ampl.var["h"][i] = val

        ampl.option['solver'] = 'ipopt'
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "expect_infeasible_problem=yes "
            "bound_relax_factor=0 "
            "tol=1e-8 "
            "constr_viol_tol=1e-8 "
            "acceptable_constr_viol_tol=1e-8 "
            "bound_push=0.01 "
            "bound_frac=0.01 "
            "warm_start_init_point=yes "
            "halt_on_ampl_error=yes "
            "max_iter=2000"
        )

        with self.suppress_output():
            ampl.solve()

        solve_result = ampl.get_value("solve_result")
        q_sol = ampl.get_variable('q').get_values().to_dict()
        h_sol = ampl.get_variable('h').get_values().to_dict()
        solve_time = ampl.get_value('_solve_elapsed_time')

        return ampl, solve_result, q_sol, h_sol, solve_time

    def solve_recover_model(self, y_input):
        """
        Recover l[i,j,k] analytically from Model 1 solution y*[i,j]
        using active segment detection on the PL cost envelope.
        No LP solve needed.

        y_input: dict {(i,j): y_value}  where y_value = y[i,j] (NOT multiplied by L)
        Returns same interface as original: (None, solve_result, l_sol, cost, solve_time)
        """
        import time
        tol = 1e-6
        t0 = time.time()

        # Build alpha[k] and C[k] from self data (same formula as AMPL model)
        omega = 10.67
        pipes = self.pipes          # ordered list of pipe labels
        alpha = {k: omega / (self.R[k]**1.852 * self.d[k]**4.87) for k in pipes}
        C     = {k: self.C[k] for k in pipes}
        L     = self.L              # dict {(i,j): length}

        # Build ordered pipe list by increasing alpha (= decreasing diameter)
        pipe_order = sorted(pipes, key=lambda k: alpha[k])

        l_sol  = {}   # {(i, j, k): length_value}
        cost   = 0.0
        feasible = True

        for (i, j), y_star in y_input.items():
            # y_input stores y[i,j] already (not y*L) -- adjust if needed
            # If your y_input was stored as y[i,j]*L[i,j], divide here:
            # y_star = y_star / L[i, j]

            L_ij = L[i, j]
            alpha_vals = [alpha[k] for k in pipe_order]

            # --- initialise all l to 0 for this arc ---
            for k in pipes:
                l_sol[(i, j, k)] = 0.0

            # ── Case 2 / Case 3: y* sits exactly on a pipe breakpoint ──
            matched_breakpoint = False
            for k in pipe_order:
                if abs(y_star - alpha[k]) <= tol:
                    l_sol[(i, j, k)] = L_ij
                    cost += C[k] * L_ij
                    matched_breakpoint = True
                    break

            if matched_breakpoint:
                continue

            # ── Case 1: y* strictly inside a segment ──
            placed = False
            for idx in range(len(pipe_order) - 1):
                k_lo = pipe_order[idx]       # pipe with smaller alpha (larger diameter)
                k_hi = pipe_order[idx + 1]   # pipe with larger  alpha (smaller diameter)
                a_lo = alpha[k_lo]
                a_hi = alpha[k_hi]

                if a_lo - tol <= y_star <= a_hi + tol:
                    denom   = a_hi - a_lo
                    lam_lo  = (a_hi - y_star) / denom   # weight on k_lo
                    lam_hi  = (y_star - a_lo) / denom   # weight on k_hi

                    # clamp tiny numerical noise to [0,1]
                    lam_lo = max(0.0, min(1.0, lam_lo))
                    lam_hi = max(0.0, min(1.0, lam_hi))

                    l_sol[(i, j, k_lo)] = lam_lo * L_ij
                    l_sol[(i, j, k_hi)] = lam_hi * L_ij
                    cost += (C[k_lo] * lam_lo + C[k_hi] * lam_hi) * L_ij
                    placed = True
                    break

            if not placed:
                # y* is outside [alpha_min, alpha_max] — clamp to nearest endpoint
                if y_star < alpha[pipe_order[0]]:
                    k_snap = pipe_order[0]
                else:
                    k_snap = pipe_order[-1]
                l_sol[(i, j, k_snap)] = L_ij
                cost += C[k_snap] * L_ij
                feasible = False   # flag: y* was out of bounds

        solve_time   = time.time() - t0
        solve_result = "solved" if feasible else "infeasible"
        return None, solve_result, l_sol, cost, solve_time

    # ============================================================
    # RECOVER l FROM y
    # ============================================================
    def solve_recover_model1(self, y_input, model_name="recover_wdnmodel.mod"):
        ampl = AMPL()
        # ampl.setOption('hsllib', '/usr/local/lib/libma57.so')
        ampl.read(model_name)
        ampl.read_data(self.data_file)

        for (i, j), val in y_input.items():
            ampl.param["y"][i, j] = val

        ampl.option['solver'] = 'cplex'

        # ampl.option["ipopt_options"] = (
        #     "outlev=0 "
        #     "expect_infeasible_problem=no "
        #     "bound_relax_factor=0 "
        #     "tol=1e-9 "
        #     "bound_push=0.01 "
        #     "bound_frac=0.01 "
        #     "warm_start_init_point=no "
        #     "halt_on_ampl_error=yes"
        # )

        #with self.suppress_output():
        ampl.solve()

        solve_result = ampl.get_value("solve_result")
        l_sol = ampl.get_variable('l').get_values().to_dict()
        cost = ampl.getObjective("total_cost").value()
        solve_time = ampl.get_value('_solve_elapsed_time')
        # ampl.eval("display l;")
        return ampl, solve_result, l_sol, cost, solve_time

    # ============================================================
    # ORIGINAL NLP RE-SOLVE WITH WARM START
    # ============================================================
    def solve_original_with_init(self, l_init, q_init, h_init):
        ampl = AMPL()
        # ampl.setOption('hsllib', '/usr/local/lib/libma57.so')
        ampl.read("wdnmodel.mod")
        ampl.read_data(self.data_file)

        for (i, j), val in q_init.items():
            ampl.var["q"][i, j] = val
        for i, val in h_init.items():
            ampl.var["h"][i] = val
        for (i, j, k), val in l_init.items():
            ampl.var["l"][i, j, k] = val

        if self.data_number == 6:
            ampl.eval("subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
        else:
            ampl.eval("subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];")

        ampl.option['solver'] = 'ipopt'
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "expect_infeasible_problem=no "
            "bound_relax_factor=0 "
            f"bound_push={self.bound_push} "
            f"bound_frac={self.bound_frac} "
            "warm_start_init_point=yes "
            "halt_on_ampl_error=yes "
            "max_iter=3000"
        )

        ampl.solve()

        solve_result = ampl.get_value("solve_result")
        q_sol = ampl.get_variable('q').get_values().to_dict()
        h_sol = ampl.get_variable('h').get_values().to_dict()
        l_sol = ampl.get_variable('l').get_values().to_dict()
        cost = ampl.getObjective("total_cost").value()
        solve_time = ampl.get_value('_solve_elapsed_time')

        return ampl, solve_result, q_sol, h_sol, l_sol, cost, solve_time

    # ============================================================
    # PERTURB y TO TRY LOWER COST
    # ============================================================
    def perturb_y(self, y_current, num_arcs=3, perturb_scale=0.05):
        """
        Increase y on selected arcs to encourage lower-cost pipe combinations.
        """
        y_new = copy.deepcopy(y_current)
        # chosen_arcs = random.sample(self.arcs, min(num_arcs, len(self.arcs)))
        chosen_arcs = self.sorted_arcs

        for arc in chosen_arcs:
            # factor = 1.0 + random.uniform(0.01, perturb_scale)  # increase y
            # y_new[arc] = max(self.y_lb[arc], min(self.y_ub[arc], y_current[arc] * factor))
            y_new[arc] = y_current[arc] + perturb_scale*(self.y_ub[arc] - y_current[arc])
            print("arc:",arc, "y_new:", y_new[arc])
            if y_new[arc]>=self.y_ub[arc] or y_new[arc]<=self.y_lb[arc]:
                print(y_new[arc] - self.y_ub[arc])
                print(self.y_lb[arc] - y_new[arc])

                # y_new[arc] = y_current[arc] + 0.1*perturb_scale*
                

        return y_new, chosen_arcs

    # ============================================================
    # SMOOTH PHI(q)
    # ============================================================
    def phi(self, q, eps):
        # return q**3 * (q**2 + eps**2)**0.426 / (q**2 + 0.426 * eps**2)
        return q*(np.abs(q)**0.852)
    # ============================================================
    # HEADLOSS RESIDUALS
    # ============================================================
    def compute_headloss_residuals(self, q_sol, h_sol, y_input):
        residuals = {}
        for (i, j) in self.arcs:
            phi_q = self.phi(q_sol[(i, j)], self.eps[(i, j)])
            residuals[(i, j)] = h_sol[i] - h_sol[j] - phi_q * y_input[(i, j)]
        return residuals

    def compute_constraint_violations(self, q_sol, h_sol, y_input):
        """
        Compute all reduced-model constraint residuals / violations.
    
        Returns:
            violations = {
                "headloss": {(i,j): residual, ...},
                "flow": {j: residual, ...},
                "source_head": {i: residual, ...},
                "pressure": {i: violation, ...},
                # "y_lb": {(i,j): violation, ...},
                # "y_ub": {(i,j): violation, ...},
                # "y_total": {(i,j): total_bound_violation, ...}
            }
        """
        violations = {
            "headloss": {},
            "flow": {},
            "source_head": {},
            "pressure": {},
            # "y_lb": {},
            # "y_ub": {},
            # "y_total": {}
        }
    
        # ------------------------------------------------------------
        # 1. Headloss constraint residuals
        #     h[i] - h[j] - phi(q[i,j]) * y[i,j] = 0
        # ------------------------------------------------------------
        for (i, j) in self.arcs:
            phi_q = self.phi(q_sol[(i, j)], self.eps[(i, j)])
            violations["headloss"][(i, j)] = h_sol[i] - h_sol[j] - phi_q * y_input[(i, j)]
    
        # ------------------------------------------------------------
        # 2. Flow conservation residuals
        #     sum inflow - sum outflow - D[j] = 0
        # ------------------------------------------------------------
        for j in self.nodes:
            if j not in self.source:
                inflow = sum(q_sol[(i, j)] for i in self.nodes if (i, j) in self.arcs)
                outflow = sum(q_sol[(j, i)] for i in self.nodes if (j, i) in self.arcs)
                violations["flow"][j] = inflow - outflow - self.D[j]
                print(f"violations[flow][{j}]",violations["flow"][j])

        # ------------------------------------------------------------
        # 3. Source head residuals
        #     h[i] - E[i] = 0
        # ------------------------------------------------------------
        for i in self.source:
            violations["source_head"][i] = h_sol[i] - self.E[i]

        # ------------------------------------------------------------
        # 4. Minimum pressure violations
        #     -h[i] + (E[i] + P[i]) <= 0
        #     violation = max(0, lhs)
        # ------------------------------------------------------------
        for i in self.nodes:
            if i not in self.source:
                lhs = -h_sol[i] + (self.E[i] + self.P[i])
                violations["pressure"][i] = max(0.0, lhs)
    
        # ------------------------------------------------------------
        # 5. y lower/upper bound violations
        #     y_lb <= y <= y_ub
        # ------------------------------------------------------------
        # for (i, j) in self.arcs:
        #     y_val = y_input[(i, j)]
        #
        #     lb_violation = max(0.0, self.y_lb[(i, j)] - y_val)
        #     ub_violation = max(0.0, y_val - self.y_ub[(i, j)])
        #
        #     violations["y_lb"][(i, j)] = lb_violation
        #     violations["y_ub"][(i, j)] = ub_violation
        #     violations["y_total"][(i, j)] = lb_violation + ub_violation
    
        return violations

    def summarize_violations(self, violations):
        """
        Compute max and sum of all residuals / violations.
        """
        summary = {}
    
        for key, vals in violations.items():
            abs_vals = [abs(v) for v in vals.values()]
            if len(abs_vals) == 0:
                summary[key] = {"max": 0.0, "sum": 0.0}
            else:
                summary[key] = {
                    "max": max(abs_vals),
                    "sum": sum(abs_vals)
                }
    
        # global totals
        all_vals = []
        for vals in violations.values():
            all_vals.extend([abs(v) for v in vals.values()])
    
        summary["global"] = {
            "max": max(all_vals) if all_vals else 0.0,
            "sum": sum(all_vals) if all_vals else 0.0
        }
    
        return summary

    def is_reduced_feasible(self, violations,
                            tol_headloss=1e-6,
                            tol_flow=1e-6,
                            tol_source=1e-6,
                            tol_pressure=1e-6,
                            tol_y=1e-10):
        """
        Check whether all reduced-model violations are acceptable.
        """
        summary = self.summarize_violations(violations)
    
        feasible = (
            summary["headloss"]["max"] <= tol_headloss and
            summary["flow"]["max"] <= tol_flow and
            summary["source_head"]["max"] <= tol_source and
            summary["pressure"]["max"] <= tol_pressure 
            # summary["y_total"]["max"] <= tol_y
        )
    
        return feasible, summary

    def project_y_bounds(self, y_input):
        """
        Project y onto [y_lb, y_ub].
        """
        y_proj = {}
        projected_arcs = []

        for arc in self.arcs:
            old_val = y_input[arc]
            new_val = min(self.y_ub[arc], max(self.y_lb[arc], old_val))
            y_proj[arc] = new_val

            if abs(new_val - old_val) > 1e-12:
                projected_arcs.append(arc)

        return y_proj, projected_arcs

    def repair_y_full(self, y_input, q_sol, h_sol, violations,
                      max_headloss_arcs=5, pressure_push=0.02):
        """
        Repair y using all available reduced-model violations.

        Repairs:
        1. y-bound violations by projection
        2. headloss residuals using exact arc correction
        3. pressure-violated nodes by slightly increasing incoming arc y
           (heuristic)

        Returns:
            y_repaired, repair_info
        """
        y_repaired = copy.deepcopy(y_input)
        repair_info = {
            "projected_arcs": [],
            "headloss_repaired_arcs": [],
            "pressure_repaired_arcs": []
        }

        # ============================================================
        # 1) Project onto y bounds
        # ============================================================
        y_repaired, projected_arcs = self.project_y_bounds(y_repaired)
        repair_info["projected_arcs"] = projected_arcs

        # ============================================================
        # 2) Repair headloss on most violated arcs
        #     y_hat = y + r / phi(q)
        # ============================================================
        headloss_res = violations["headloss"]

        violated_headloss = sorted(
            self.arcs,
            key=lambda a: abs(headloss_res[a]),
            reverse=True
        )

        for arc in violated_headloss[:max_headloss_arcs]:
            i, j = arc
            phi_q = self.phi(q_sol[(i, j)], self.eps[(i, j)])

            if abs(phi_q) < 1e-10:
                continue

            r_ij = headloss_res[arc]
            y_hat = y_repaired[arc] + r_ij / phi_q

            # project to bounds
            y_hat = min(self.y_ub[arc], max(self.y_lb[arc], y_hat))

            if abs(y_hat - y_repaired[arc]) > 1e-12:
                y_repaired[arc] = y_hat
                repair_info["headloss_repaired_arcs"].append(arc)

        # ============================================================
        # 3) Pressure repair heuristic
        #     If node i violates pressure, slightly increase y on arcs
        #     entering i (heuristic: lower resistance upstream)
        # ============================================================
        for node, viol in violations["pressure"].items():
            if viol > 1e-6:
                incoming_arcs = [(u, v) for (u, v) in self.arcs if v == node]

                for arc in incoming_arcs:
                    old_y = y_repaired[arc]
                    new_y = old_y * (1.0 + pressure_push)
                    new_y = min(self.y_ub[arc], max(self.y_lb[arc], new_y))

                    if abs(new_y - old_y) > 1e-12:
                        y_repaired[arc] = new_y
                        repair_info["pressure_repaired_arcs"].append(arc)

        return y_repaired, repair_info

    # ============================================================
    # REPAIR y USING HEADLOSS RESIDUAL
    # ============================================================
    def repair_y(self, y_input, q_sol, h_sol, residuals, max_repair_arcs=5):
        y_repaired = copy.deepcopy(y_input)

        # sort arcs by absolute violation
        violated = sorted(self.arcs, key=lambda a: abs(residuals[a]), reverse=True)

        repaired_arcs = []
        for arc in violated[:max_repair_arcs]:
            i, j = arc
            phi_q = self.phi(q_sol[(i, j)], self.eps[(i, j)])

            # if abs(phi_q) < 1e-10:
            #     continue

            # exact repair:
            # y_hat = y + r/phi(q)
            y_hat = y_repaired[arc] + residuals[arc] / phi_q

            if self.y_lb[arc] <= y_hat <= self.y_ub[arc]:
                y_repaired[arc] = y_hat
                repaired_arcs.append(arc)

        return y_repaired, repaired_arcs

    # ============================================================
    # MAIN IMPROVEMENT HEURISTIC
    # ============================================================
    def improve_solution_by_y_search(self, max_outer_iter=1, max_repair_iter=1,
                                     num_arcs_perturb=5, perturb_scale=0.05):

        print("\n=================== Starting y-space improvement heuristic ===================")
        incumbent_cost = self.best_cost
        incumbent_y = copy.deepcopy(self.best_y)
        incumbent_q = copy.deepcopy(self.best_q)
        incumbent_h = copy.deepcopy(self.best_h)
        incumbent_l = copy.deepcopy(self.best_l)

        print(f"Initial incumbent cost = {incumbent_cost:.8f}")

        for outer in range(max_outer_iter):
            print(f"\n---------------- OUTER ITERATION {outer+1} ----------------")

            # Step 1: perturb y
            y_trial, changed_arcs = self.perturb_y(
                incumbent_y,
                num_arcs=num_arcs_perturb,
                perturb_scale=perturb_scale
            )

            print(f"Perturbed arcs: {changed_arcs}")

            # Step 2: solve reduced NLP for q,h
            reduced_ampl, solve_result, q_trial, h_trial, t_red = self.solve_reduced_nlp(
                y_trial, q_init=incumbent_q, h_init=incumbent_h
            )
            print(f"Reduced NLP solve_result = {solve_result}, time = {t_red:.2f} sec")

            # # If infeasible, attempt repair
            # ============================================================
            # FULL FEASIBILITY REPAIR LOOP
            # ============================================================
            repaired = False

            for repair_iter in range(max_repair_iter):
                print(f"\nRepair iteration {repair_iter+1}")

                # Compute all constraint violations
                violations = self.compute_constraint_violations(q_trial, h_trial, y_trial)
                feasible, summary = self.is_reduced_feasible(violations)

                print("Violation summary:")
                print(f"  headloss max   = {summary['headloss']['max']:.4e}")
                print(f"  flow max       = {summary['flow']['max']:.4e}")
                print(f"  source head max= {summary['source_head']['max']:.4e}")
                print(f"  pressure max   = {summary['pressure']['max']:.4e}")
                # print(f"  y-bound max    = {summary['y_total']['max']:.4e}")
                print(f"  global max     = {summary['global']['max']:.4e}")
                print(f"  reduced NLP status = {solve_result}")

                # Accept only if ALL constraints are acceptable
                if feasible and ("solved" in str(solve_result).lower() or "limit" in str(solve_result).lower()):
                    repaired = True
                    print("Reduced candidate is globally feasible enough.")
                    break
            
                # Repair all relevant violations
                y_trial, repair_info = self.repair_y_full(y_trial, q_trial, h_trial, violations)
            
                print(f"Projected arcs        : {repair_info['projected_arcs']}")
                print(f"Headloss repaired arcs: {repair_info['headloss_repaired_arcs']}")
                print(f"Pressure repaired arcs: {repair_info['pressure_repaired_arcs']}")
            
                # Re-solve reduced NLP after repair
                reduced_ampl, solve_result, q_trial, h_trial, t_red = self.solve_reduced_nlp(
                    y_trial,
                    q_init=q_trial,
                    h_init=h_trial
                )
            
                print(f"Post-repair reduced NLP solve_result = {solve_result}, time = {t_red:.2f} sec")
            
            # Final check after repair loop
            if not repaired:
                violations = self.compute_constraint_violations(q_trial, h_trial, y_trial)
                feasible, summary = self.is_reduced_feasible(violations)
            
                if feasible and ("solved" in str(solve_result).lower() or "limit" in str(solve_result).lower()):
                    repaired = True
            
            if not repaired:
                print("Could not repair reduced NLP sufficiently. Skipping candidate.")
                continue

            # Step 3: recover l from y
            rec_ampl, rec_result, l_trial, recovered_cost, t_rec = self.solve_recover_model(y_trial)
            print(f"Recover model result = {rec_result}, recovered cost = {recovered_cost:.8f}, time = {t_rec:.2f} sec")

            # Step 4: accept only if reduced recovered cost improves
            if recovered_cost < incumbent_cost - 1e-6:
                print(f"Candidate improved recovered cost: {recovered_cost:.8f} < {incumbent_cost:.8f}")

                # Step 5: solve original NLP using recovered point as warm start
                orig_ampl, orig_result, q_new, h_new, l_new, final_cost, t_orig = self.solve_original_with_init(
                    l_trial, q_trial, h_trial
                )

                print(f"Original warm-start NLP result = {orig_result}, final cost = {final_cost:.8f}, time = {t_orig:.2f} sec")

                if final_cost < incumbent_cost - 1e-6:
                    print(f"*** ACCEPTED IMPROVED LOCAL SOLUTION: {final_cost:.8f} ***")
                    incumbent_cost = final_cost
                    incumbent_q = copy.deepcopy(q_new)
                    incumbent_h = copy.deepcopy(h_new)
                    incumbent_l = copy.deepcopy(l_new)
                    incumbent_y = self.compute_y_from_l(incumbent_l)

                    self.best_cost = incumbent_cost
                    self.best_q = copy.deepcopy(incumbent_q)
                    self.best_h = copy.deepcopy(incumbent_h)
                    self.best_l = copy.deepcopy(incumbent_l)
                    self.best_y = copy.deepcopy(incumbent_y)
                else:
                    print("Warm-start original NLP did not improve incumbent.")
            else:
                print("Recovered cost not better than incumbent. Reject.")

        print("\n=================== Heuristic finished ===================")
        print(f"Best cost found = {self.best_cost:.8f}")

    def plot_value_function_for_arc_old(self, arc, n_points=300):
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.optimize import linprog

        (u, v) = arc
        Lij = self.L[(u, v)]

        pipes = list(self.pipes)
    
        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipes
        }
    
        alpha_vals = np.array([alpha[k] for k in pipes])
        cost_vals  = np.array([self.C[k] for k in pipes])
    
        y_min = Lij * np.min(alpha_vals)
        y_max = Lij * np.max(alpha_vals)
    
        y_grid = np.linspace(y_min, y_max, n_points)
        z_vals = []
    
        for y in y_grid:
            c = cost_vals
            A_eq = np.vstack([alpha_vals, np.ones(len(pipes))])
            b_eq = np.array([y, Lij])
            bounds = [(0, None) for _ in pipes]
    
            res = linprog(c=c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    
            if res.success:
                z_vals.append(res.fun)
            else:
                z_vals.append(np.nan)
    
        plt.figure(figsize=(8,5))
        plt.plot(y_grid, z_vals, lw=2)
        plt.xlabel(r"$y_{ij}$")
        plt.ylabel(r"$Z_{ij}(y_{ij})$")
        plt.title(f"Objective Function Graph for Arc ({u},{v})")
        plt.savefig(f"value_function_arc_{u}_{v}.png", dpi=300, bbox_inches="tight")
        plt.grid(True)
        plt.show()
    
        return y_grid, z_vals

    def plot_value_function_for_arc(self, arc, n_points=300):
        import numpy as np
        import matplotlib.pyplot as plt
        from amplpy import AMPL

        (i, j) = arc
        Lij = self.L[(i, j)]

        pipes = list(self.pipes)
        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipes
        }

        alpha_vals = np.array([alpha[k] for k in pipes])
        y_min = Lij * np.min(alpha_vals)
        y_max = Lij * np.max(alpha_vals)

        y_grid = np.linspace(y_min, y_max, n_points)
        z_vals = []

        # Initialize AMPL
        ampl = AMPL()
        ampl.read("value_function.mod")
        # ampl.read_data(self.data_file)

        ampl.set["pipes"] = pipes

        for k in pipes:
            ampl.param["C"][k] = self.C[k]
            ampl.param["alpha"][k] = alpha[k]

        ampl.param["L"] = Lij

        ampl.setOption("solver", "cplex")   # or gurobi/highs
        ampl.setOption("solver_msg", 0)

        # Solve for each y
        for y in y_grid:
            ampl.param["y"] = float(y)

            ampl.solve()

            if ampl.get_value("solve_result") == "solved":
                z_vals.append(ampl.get_objective("total_cost").value())
            else:
                z_vals.append(np.nan)

        # Plot
        plt.figure(figsize=(8,5))
        plt.plot(y_grid, z_vals, lw=2)
        plt.xlabel(r"$y_{ij}$")
        plt.ylabel(r"$Z_{ij}(y_{ij})$")
        plt.title(f"Value Function for Arc ({i},{j})")
        plt.grid(True)
        plt.savefig(f"value_function_arc_{i}_{j}.png", dpi=300, bbox_inches="tight")
        plt.show()

        return y_grid, z_vals

    
    def plot_value_function_for_arcs(self, arcs, n_points=40, seed=42):
        pipes = list(self.pipes)
        # Compute alpha values
        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipes
        }
        # 🔥 STEP 1: GLOBAL y RANGE (valid for ALL arcs)
        alpha_vals = np.array(list(alpha.values()))
        L_vals = np.array(list(self.L.values()))
        global_y_min = np.min(L_vals) * np.min(alpha_vals)
        global_y_max = np.max(L_vals) * np.max(alpha_vals)
        # 🔥 SAME y points for everything
        np.random.seed(seed)
        y_global = np.sort(
            np.random.uniform(global_y_min, global_y_max, n_points)
        )
        print("\nCommon y points used for ALL arcs:\n", y_global)
        # Initialize AMPL
        ampl = AMPL()
        ampl.read("value_function.mod")
        ampl.set["pipes"] = pipes
        for k in pipes:
            ampl.param["C"][k] = self.C[k]
            ampl.param["alpha"][k] = alpha[k]
        ampl.setOption("solver", "cplex")
        ampl.setOption("solver_msg", 0)
        plt.figure(figsize=(9, 6))
        # Loop over arcs
        for arc in arcs:
            (i, j) = arc
            Lij = self.L[(i, j)]
            alpha_vals = np.array([alpha[k] for k in pipes])
            y_min = Lij * np.min(alpha_vals)
            y_max = Lij * np.max(alpha_vals)
            z_vals = []
            y_used = []
            ampl.param["L"] = Lij
            # 🔥 STEP 2: USE SAME y, FILTER FEASIBLE
            for y in y_global:
                if y_min <= y <= y_max:
                    ampl.param["y"] = float(y)
                    ampl.solve()
                    if ampl.get_value("solve_result") == "solved":
                        cost = ampl.get_objective("total_cost").value()
                        z_vals.append(cost)
                        y_used.append(y)
                    else:
                        z_vals.append(np.nan)
            print(f"\nArc ({i},{j}) values:")
            for y, z in zip(y_used, z_vals):
                print(f"(y, z) = ({y:.3f}, {z:.3f})")
            # Plot
            line, = plt.plot(y_used, z_vals, lw=2, label=f"Arc ({i},{j})")
            plt.scatter(y_used, z_vals, color=line.get_color(), s=40)
            # Annotate
            for y, z in zip(y_used, z_vals):
                if not np.isnan(z):
                    plt.annotate(
                        f"({y:.1f}, {z:.1f})",
                        (y, z),
                        textcoords="offset points",
                        xytext=(5, 5),
                        fontsize=8,
                        color=line.get_color()
                    )
        plt.xlabel(r"$y_{ij}$")
        plt.ylabel(r"$Z_{ij}(y_{ij})$")
        plt.title("Value Functions (Same Global y Points)")
        plt.legend()
        plt.grid(True)
        plt.savefig("value_function_same_y.png", dpi=300, bbox_inches="tight")
        plt.show()

   
    def plot_value_function_for_arcs1(self, arcs, n_points=10):

        pipes = list(self.pipes)

        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipes
        }

        # Initialize AMPL once
        ampl = AMPL()
        ampl.read("value_function.mod")

        ampl.set["pipes"] = pipes

        for k in pipes:
            ampl.param["C"][k] = self.C[k]
            ampl.param["alpha"][k] = alpha[k]

        ampl.setOption("solver", "cplex")
        ampl.setOption("solver_msg", 0)

        plt.figure(figsize=(8,5))

        # Loop over arcs
        for arc in arcs:
            (i, j) = arc
            Lij = self.L[(i, j)]

            alpha_vals = np.array([alpha[k] for k in pipes])
            y_min = Lij * np.min(alpha_vals)
            y_max = Lij * np.max(alpha_vals)

            y_grid = np.linspace(y_min, y_max, n_points)
            z_vals = []

            ampl.param["L"] = Lij

            # Solve for each y
            for y in y_grid:
                ampl.param["y"] = float(y)

                ampl.solve()

                if ampl.get_value("solve_result") == "solved":
                    cost = ampl.get_objective("total_cost").value()
                    z_vals.append(cost)
                    # print(cost)
                else:
                    print("infeasible")
                    z_vals.append(np.nan)
            print(z_vals)
            # Plot each arc
            plt.plot(y_grid, z_vals, lw=2, label=f"Arc ({i},{j})")

        # Final plot settings
        plt.xlabel(r"$y_{ij}$")
        plt.ylabel(r"$Z_{ij}(y_{ij})$")
        plt.title("Value Functions for Multiple Arcs")
        plt.legend()
        plt.grid(True)
        plt.savefig("value_function_multiple_arcs.png", dpi=300, bbox_inches="tight")
        plt.show()

    
    def plot_m_values(self):
        # Ensure pipes are ordered
        pipe_list = list(self.pipes)
    
        # Compute alpha and C
        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipe_list
        }
    
        C = {
            k: self.C[k]
            for k in pipe_list
        }
    
        # Compute m[k] and b[k]
        m = {}
        b = {}
        x_vals = []
        m_vals = []
        b_vals = []
    
        for idx in range(len(pipe_list) - 1):
            k1 = pipe_list[idx]
            k2 = pipe_list[idx + 1]
    
            denom = alpha[k1] - alpha[k2]
    
            if abs(denom) < 1e-14:
                print(f"Skipping k={k1} because alpha[{k1}] == alpha[{k2}]")
                continue
    
            m[k1] = (C[k1] - C[k2]) / denom
            b[k1] = (alpha[k1]*C[k2] - alpha[k2]*C[k1]) / denom
    
            x_vals.append(k1)
            m_vals.append(m[k1])
            b_vals.append(b[k1])
    
        # Print values
        print("\nComputed values:")
        for k in m:
            print(f"k={k}: m={m[k]:.8f}, b={b[k]:.8f}")
    
        # 🔥 Create subplots
        fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    
        # Plot m[k]
        axes[0].plot(x_vals, m_vals, marker='o', linestyle='-')
        axes[0].set_ylabel(r"$m_k$")
        axes[0].set_title(r"Slope: $m_k = \frac{C_k - C_{k+1}}{\alpha_k - \alpha_{k+1}}$")
        axes[0].grid(True)
    
        # Plot b[k]
        axes[1].plot(x_vals, b_vals, marker='s', linestyle='-')
        axes[1].set_xlabel("Pipe index k")
        axes[1].set_ylabel(r"$b_k$")
        axes[1].set_title(r"Intercept: $b_k = \frac{\alpha_k C_{k+1} - \alpha_{k+1} C_k}{\alpha_k - \alpha_{k+1}}$")
        axes[1].grid(True)
    
        # X ticks
        plt.xticks(x_vals)
    
        # Save and show
        plt.tight_layout()
        plt.savefig("m_b_values_plot.png", dpi=300, bbox_inches='tight')
        plt.show()
    
        return m, b

    def plot_m_values1(self):
        # Ensure pipes are ordered
        pipe_list = list(self.pipes)
    
        # Compute alpha and C
        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipe_list
        }
    
        C = {
            k: self.C[k]
            for k in pipe_list
        }
    
        # Compute m[k] for consecutive pairs
        m = {}
        b = {}
        x_vals = []
        y_vals = []
        z_vals = []
    
        for idx in range(len(pipe_list) - 1):
            k1 = pipe_list[idx]
            k2 = pipe_list[idx + 1]
    
            denom = alpha[k1] - alpha[k2]
    
            if abs(denom) < 1e-14:
                print(f"Skipping k={k1} because alpha[{k1}] == alpha[{k2}]")
                continue
    
            m[k1] = (C[k1] - C[k2]) / denom
            b[k1] = (alpha[k1]*C[k2] - alpha[k2]* C[k1]) / denom
            x_vals.append(k1)
            y_vals.append(m[k1])
            z_vals.append(b[k1])
    
        # Print values
        print("\nComputed m[k] values:")
        for k in m:
            print(f"m[{k}] = {m[k]:.8f}")
            print(f"b[{k}] = {b[k]:.8f}")
    
        # Plot
        plt.figure(figsize=(8, 5))
        plt.plot(x_vals, y_vals, marker='o', linestyle='-')
        plt.xlabel("Pipe index k")
        plt.ylabel("z[k]")
        plt.title(r"Plot of $m_k = \frac{C_k - C_{k+1}}{\alpha_k - \alpha_{k+1}}$")
        # plt.title(r"Plot of $b_k = \frac{\alpha_k C_{k+1} -\alpha_{k+1}  C_{k}}{\alpha_k - \alpha_{k+1}}$")
        plt.grid(True)
        plt.xticks(x_vals)

        plt.savefig("m_values_plot.png", dpi=300, bbox_inches='tight')

        plt.show()

        return m


    def plot_piecewise_lines(self, arc, n_points=200):

        (i, j) = arc
        Lij = self.L[(i, j)]

        pipe_list = list(self.pipes)

        # Compute alpha and C
        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipe_list
        }

        C = {
            k: self.C[k]
            for k in pipe_list
        }

        # Sort pipes by alpha (IMPORTANT for correct ordering)
        pipe_list = sorted(pipe_list, key=lambda k: alpha[k])

        # Compute m[k], b[k]
        m = {}
        b = {}

        for idx in range(len(pipe_list) - 1):
            print(idx)
            k1 = pipe_list[idx]
            k2 = pipe_list[idx + 1]

            denom = alpha[k1] - alpha[k2]

            if abs(denom) < 1e-14:
                continue

            m[k1] = (C[k1] - C[k2]) / denom
            b[k1] = (alpha[k1]*C[k2] - alpha[k2]*C[k1]) / denom

        # 🔥 y-range (first quadrant)
        alpha_vals = np.array([alpha[k] for k in pipe_list])
        y_min = max(0, Lij * np.min(alpha_vals))
        y_max = Lij * np.max(alpha_vals)

        y_grid = np.linspace(y_min, y_max, n_points)

        plt.figure(figsize=(15, 6))

        print("m",m)
        # 🔥 Plot all lines
        for k in m:
            print(k)
            z_vals = np.array([m[k]*y + Lij*b[k] for y in y_grid])
            # z_vals = np.maximum(z_vals, 0)
            # plt.plot(y_grid, z_vals, linestyle='-', label=f"line k={k-1}")
            # plt.plot(y_grid,z_vals,linestyle='-',label=rf"$z(y) = {m[k]:.3g}\,y + {b[k]:.3g}\,L$")
            plt.plot(y_grid,z_vals,linestyle='-',label=rf"$z(y) = m_{{{k-1}}} y + b_{{{k-1}}} L$")


        # 🔥 Plot envelope (min of lines)
        z_env = []
        for y in y_grid:
            z_env.append(max(m[k]*y + Lij*b[k] for k in m))
        # z_env = np.maximum(z_env)

        # plt.plot(y_grid, z_env, color='black', linewidth=3, label="Envelope (max)")

        # 🔥 Vertical lines at y = alpha_k * Lij
        for k in pipe_list:
            yk = alpha[k] * Lij

            if yk >= 0:
                plt.axvline(x=yk, linestyle=':', linewidth=1.5, color="gray")
                plt.text(yk, -plt.ylim()[1]*0.01, rf"$\alpha_{{{k}}} L$", rotation=0, horizontalalignment='center', fontsize=8)


        # 🔴 Plot only breakpoint points
        for k in pipe_list:
            xk = alpha[k] * Lij
            zk = C[k] * Lij
        
            plt.scatter(xk, zk, s=50)
            plt.text(
                xk,
                zk,
                rf"$(\alpha_{{{k}}}L,\;C_{{{k}}}L)$",
                fontsize=8,
                ha='left',
                va='bottom'
            )



        # Axis settings
        plt.xlim(left=0)
        plt.xticks([])
        # print(self.C)
        # print(list(self.C))
        plt.ylim(0, 1.05*self.C[next(reversed(self.C))]*Lij)

        # plt.xlim(0, 0.2 * max(alpha[k]*Lij for k in pipe_list))
        
        plt.xlabel(r"$y_{ij}$")
        plt.ylabel(r"$z_{ij}(y_{ij})$")
        plt.title(f"Piecewise Lines + Breakpoints for Arc ({i},{j})")

        plt.legend()
        # plt.grid(True)
        plt.savefig("piecewise_with_breakpoints.png", dpi=300, bbox_inches="tight")
        plt.show()


    def plot_piecewise_lines1(self, arc, n_points=100):
    
        (i, j) = arc
        Lij = self.L[(i, j)]
    
        pipe_list = list(self.pipes)
    
        # Compute alpha and C
        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipe_list
        }
    
        C = {
            k: self.C[k]
            for k in pipe_list
        }
    
        # Compute m[k], b[k]
        m = {}
        b = {}
    
        for idx in range(len(pipe_list) - 1):
            k1 = pipe_list[idx]
            k2 = pipe_list[idx + 1]
    
            denom = alpha[k1] - alpha[k2]
    
            if abs(denom) < 1e-14:
                continue
    
            m[k1] = (C[k1] - C[k2]) / denom
            b[k1] = (alpha[k1]*C[k2] - alpha[k2]*C[k1]) / denom
    
        # 🔥 First quadrant y-range
        alpha_vals = np.array([alpha[k] for k in pipe_list])
        y_min = max(0, Lij * np.min(alpha_vals))
        y_max = Lij * np.max(alpha_vals)
    
        y_grid = np.linspace(y_min, y_max, n_points)
    
        plt.figure(figsize=(8, 5))
    
        # 🔥 Plot each line (clip z ≥ 0)
        for k in m:
            z_vals = np.array([m[k]*y + Lij*b[k] for y in y_grid])
    
            # ensure first quadrant
            z_vals = np.maximum(z_vals, 0)
    
            plt.plot(y_grid, z_vals, linestyle='--', label=f"k={k}")
    
        # Axis strictly in first quadrant
        plt.xlim(left=0)
        plt.ylim(bottom=0)
    
        # Labels
        plt.xlabel(r"$y_{ij} \geq 0$")
        plt.ylabel(r"$z_{ij} \geq 0$")
        plt.title(f"Piecewise Linear Lines (First Quadrant) for Arc ({i},{j})")
    
        plt.legend()
        plt.grid(True)
    
        plt.savefig("piecewise_lines_first_quadrant.png", dpi=300, bbox_inches="tight")
        plt.show()


    def build_piecewise_reduced_cost(self, arc=None):
        """
        Build the exact piecewise reduced recovery cost z_ij(y_ij)
        using adjacent pipe pairs (2-BFS structure of LP recovery model).
    
        Returns
        -------
        segments : list of dict
            Each segment contains:
            {
                "k": pipe k,
                "kp1": pipe k+1,
                "m": slope,
                "b": intercept,
                "y_low": lower bound of segment,
                "y_high": upper bound of segment
            }
        """
    
        if arc is None:
            arc = list(self.arcs)[0]   # use first arc if none provided
    
        i, j = arc
        Lij = self.L[(i, j)]
    
        pipes = sorted(list(self.pipes))
    
        # alpha_k = 10.67 / (R_k^1.852 * d_k^4.87)
        alpha = {
            k: 10.67 / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in pipes
        }
    
        # sort pipes by alpha descending (same as increasing diameter)
        pipes_sorted = sorted(pipes, key=lambda k: alpha[k], reverse=True)
    
        segments = []
    
        for idx in range(len(pipes_sorted) - 1):
            k = pipes_sorted[idx]
            kp1 = pipes_sorted[idx + 1]
    
            alpha_k = alpha[k]
            alpha_kp1 = alpha[kp1]
    
            # slope
            m_k = (self.C[k] - self.C[kp1]) / (alpha_k - alpha_kp1)
    
            # intercept
            b_k = Lij * (self.C[kp1] - m_k * alpha_kp1)
    
            # feasible y interval for this segment
            y_low = alpha_kp1 * Lij
            y_high = alpha_k * Lij
    
            segments.append({
                "k": k,
                "kp1": kp1,
                "m": m_k,
                "b": b_k,
                "y_low": y_low,
                "y_high": y_high
            })
    
        return segments

    def print_piecewise_reduced_cost(self, arc=None):
        """
        Print the piecewise linear reduced recovery cost z_ij(y_ij)
        for a selected arc.
        """
    
        if arc is None:
            arc = list(self.arcs)[0]
    
        segments = self.build_piecewise_reduced_cost(arc)
    
        print(f"\nPiecewise reduced recovery cost for arc {arc}:\n")
    
        for s in segments:
            print(
                f"Segment ({s['k']},{s['kp1']}): "
                f"y in [{s['y_low']:.6f}, {s['y_high']:.6f}]  -->  "
                f"z(y) = {s['m']:.6f} y + {s['b']:.6f}"
            )
    def plot_piecewise_reduced_cost(self, arc=None, save_path="piecewise_reduced_cost.png"):
        """
        Plot the exact piecewise reduced recovery cost z_ij(y_ij)
        for a selected arc and save the figure.
        """
    
        if arc is None:
            arc = list(self.arcs)[0]
    
        segments = self.build_piecewise_reduced_cost(arc)
    
        import numpy as np
        import matplotlib.pyplot as plt
    
        plt.figure(figsize=(10, 6))
    
        # plot each segment
        for s in segments:
            y_vals = np.linspace(s["y_low"], s["y_high"], 100)
            z_vals = s["m"] * y_vals + s["b"]
    
            plt.plot(
                y_vals,
                z_vals,
                linewidth=2,
                label=f"({s['k']},{s['kp1']})"
            )
    
        # draw breakpoints
        breakpoints = []
        for s in segments:
            breakpoints.append(s["y_low"])
            breakpoints.append(s["y_high"])
    
        for bp in sorted(set(breakpoints)):
            plt.axvline(bp, linestyle='--', alpha=0.3)
    
        plt.xlabel(r"$y_{ij}$")
        plt.ylabel(r"$z_{ij}(y_{ij})$")
        plt.title(f"Exact Piecewise Reduced Recovery Cost for arc {arc}")
        plt.grid(True, alpha=0.3)
        plt.legend(title="Adjacent pipe pair", fontsize=8)
        plt.tight_layout()
    
        plt.savefig(save_path, dpi=300)
        plt.show()
    
        print(f"\nPlot saved as: {save_path}")

    def plot_piecewise_reduced_cost_envelope(self, arc=None, save_path="piecewise_reduced_cost_envelope.png"):
        """
        Plot the full exact continuous envelope of z_ij(y_ij)
        for a selected arc and save the figure.
        """
    
        if arc is None:
            arc = list(self.arcs)[0]
    
        segments = self.build_piecewise_reduced_cost(arc)
    
        import numpy as np
        import matplotlib.pyplot as plt
    
        y_all = []
        z_all = []
    
        for s in segments:
            y_vals = np.linspace(s["y_low"], s["y_high"], 200)
            z_vals = s["m"] * y_vals + s["b"]
    
            y_all.extend(y_vals)
            z_all.extend(z_vals)
    
        # sort by y
        pairs = sorted(zip(y_all, z_all), key=lambda t: t[0])
        y_all = [p[0] for p in pairs]
        z_all = [p[1] for p in pairs]
    
        plt.figure(figsize=(10, 6))
        plt.plot(y_all, z_all, linewidth=2)
    
        # draw breakpoints
        breakpoints = []
        for s in segments:
            breakpoints.append(s["y_low"])
            breakpoints.append(s["y_high"])
    
        for bp in sorted(set(breakpoints)):
            plt.axvline(bp, linestyle='--', alpha=0.3)
    
        plt.xlabel(r"$y_{ij}$")
        plt.ylabel(r"$z_{ij}(y_{ij})$")
        plt.title(f"Exact Piecewise Recovery Cost Envelope for arc {arc}")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
    
        plt.savefig(save_path, dpi=300)
        plt.show()
    
        print(f"\nEnvelope plot saved as: {save_path}")

    # ============================================================
    # 🔷 ORIGINAL → REDUCED
    # ============================================================
    def map_original_to_reduced(self, l, q, h):
        # ---- precompute alpha ----
        self.alpha = {
            k: self.omega / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in self.pipes
        }

        self.alpha_min = min(self.alpha.values())
        self.alpha_max = max(self.alpha.values())
        self.NP = len(self.pipes)
        self.segs = list(range(self.NP - 1))

        self.slope = {}
        self.intercept = {}

        for s in self.segs:

            k1 = self.pipes[s]
            k2 = self.pipes[s + 1]

            alpha1 = self.alpha[k1]
            alpha2 = self.alpha[k2]

            C1 = self.C[k1]
            C2 = self.C[k2]

            # ---- slope ----
            self.slope[s] = (C1 - C2) / (alpha1 - alpha2)

            # ---- intercept ----
            self.intercept[s] = (
                alpha1 * C2 - alpha2 * C1
            ) / (alpha1 - alpha2)
        
        y = {}
        z = {}

        for (i, j) in self.arcs:

            # ---- compute y ----
            y[(i, j)] = sum(
                self.alpha[k] * l[(i, j, k)]/self.L[i,j]
                for k in self.pipes
            )

            # ---- compute z (epigraph) ----
            z[(i, j)] = max(
                (self.slope[s] * y[(i, j)]
                + self.intercept[s])*self.L[i,j]
                for s in self.segs
            )

        return y, z

    # ============================================================
    # 🔷 REDUCED FEASIBILITY CHECK
    # ============================================================
    def check_reduced_feasibility(self, q, h, y, z, tol=1e-6):

        violations = {
            "flow": 0,
            "headloss": 0,
            "pressure": 0,
            "y_bounds": 0,
            "epigraph": 0
        }

        # ---- flow balance ----
        for j in self.nodes:
            if j in self.source:
                continue

            inflow = sum(q[(i, j)] for (i, j2) in self.arcs if j2 == j)
            outflow = sum(q[(j, i)] for (j2, i) in self.arcs if j2 == j)

            res = inflow - outflow - self.D[j]
            violations["flow"] = max(violations["flow"], abs(res))

        # ---- headloss ----
        for (i, j) in self.arcs:

            qij = q[(i, j)]

            phi = (
                qij**3 *
                (qij**2 + self.eps[(i, j)]**2)**0.426 /
                (qij**2 + 0.426 * self.eps[(i, j)]**2)
            )

            res = h[i] - h[j] - phi * y[(i, j)]*self.L[i,j]

            violations["headloss"] = max(
                violations["headloss"], abs(res)
            )

        # ---- pressure ----
        for i in self.nodes:
            if i in self.source:
                continue

            violations["pressure"] = max(
                violations["pressure"],
                max(0, self.E[i] + self.P[i] - h[i])
            )

        # ---- y bounds ----
        for (i, j) in self.arcs:

            lb = self.alpha_min 
            ub = self.alpha_max

            violations["y_bounds"] = max(
                violations["y_bounds"],
                max(0, lb - y[(i, j)], y[(i, j)] - ub)
            )

        # ---- epigraph ----
        for (i, j) in self.arcs:
            for s in self.segs:

                rhs = (self.slope[s] * y[(i, j)] + self.intercept[s]) * self.L[(i, j)]

                violations["epigraph"] = max(
                    violations["epigraph"],
                    max(0, rhs - z[(i, j)])
                )

        return violations

    # ============================================================
    # 🔷 REDUCED → ORIGINAL (RECONSTRUCTION)
    # ============================================================
    def reconstruct_l_for_arc(self, i, j, yij):

        pipes = list(self.pipes)
        n = len(pipes)

        c = np.zeros(n)

        A_eq = [
            [1]*n,
            [self.alpha[k] for k in pipes]
        ]

        b_eq = [self.L[(i, j)], yij]

        bounds = [(0, self.L[(i, j)]) for _ in pipes]

        res = linprog(
            c,
            A_eq=A_eq,
            b_eq=b_eq,
            bounds=bounds,
            method='highs'
        )

        if not res.success:
            return None

        return {
            pipes[k]: res.x[k]
            for k in range(n)
        }

    def map_reduced_to_original(self, y):

        l = {}

        for (i, j) in self.arcs:

            sol = self.reconstruct_l_for_arc(i, j, y[(i, j)]*self.L[i,j])

            if sol is None:
                raise ValueError(f"Infeasible reconstruction at arc {(i,j)}")

            for k in self.pipes:
                l[(i, j, k)] = sol[k]

        return l

    # ============================================================
    # 🔷 ORIGINAL FEASIBILITY CHECK
    # ============================================================
    def check_original_feasibility(self, q, h, l, tol=1e-6):

        violations = {
            "flow": 0,
            "headloss": 0,
            "length": 0,
            "pressure": 0
        }

        # ---- flow ----
        for j in self.nodes:
            if j in self.source:
                continue

            inflow = sum(q[(i, j)] for (i, j2) in self.arcs if j2 == j)
            outflow = sum(q[(j, i)] for (j2, i) in self.arcs if j2 == j)

            res = inflow - outflow - self.D[j]
            violations["flow"] = max(violations["flow"], abs(res))

        # ---- headloss ----
        for (i, j) in self.arcs:

            coeff = sum(
                self.alpha[k] * l[(i, j, k)]
                for k in self.pipes
            )

            qij = q[(i, j)]

            phi = (
                qij**3 *
                (qij**2 + self.eps[(i, j)]**2)**0.426 /
                (qij**2 + 0.426 * self.eps[(i, j)]**2)
            )

            res = h[i] - h[j] - phi * coeff

            violations["headloss"] = max(
                violations["headloss"], abs(res)
            )

        # ---- length ----
        for (i, j) in self.arcs:

            total = sum(
                l[(i, j, k)] for k in self.pipes
            )

            violations["length"] = max(
                violations["length"],
                abs(total - self.L[(i, j)])
            )

        # ---- pressure ----
        for i in self.nodes:
            if i in self.source:
                continue

            violations["pressure"] = max(
                violations["pressure"],
                max(0, self.E[i] + self.P[i] - h[i])
            )

        return violations

    def cycle_basis(self):
        root = self.ampl.getSet('Source').to_list()[0]
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = self.ampl.getSet('arcs').to_list() 
        uwg = nx.Graph()
        uwg.add_nodes_from(nodes_list)
        uwg.add_edges_from(edges_list)
        # print("Edges in the undirected graph:", edges_list)
        cycle_basis = nx.cycle_basis(uwg, root)
        print("cycle basis for given water network: ",cycle_basis)
        return cycle_basis



    def solve_exact_reduced_model(self):
        # print("\n-------------------------------- Solving Exact Reduced Model --------------------------")
        ampl = AMPL()
        # ampl.setOption('hsllib', '/usr/local/lib/libma57.so')
        ampl.read("exact_reduced_wdn.mod")
        # ampl.read("cycle_basis_wdn.mod")
        ampl.read_data(self.data_file)

        # ampl.eval(f"""minimize total_cost:sum{{(i,j) in arcs}} z[i,j];""")
        
        # cycle_ids = [f'c{k+1}' for k in range(len(self.cycles))]
        # ampl.set['cycles'] = cycle_ids
        # q_p = {(i,j): 1.0 for (i,j) in self.arcs}
        # ampl.param['q_p'] = q_p
        #
        # source_head = self.E[self.source[0]]   # example: 30 m head above elevation
        # ampl.param['source_head'] = source_head

        # for (i, j), val in self.q.items():
        #     ampl.var["q"][i, j] = val
        # for i, val in self.h.items():
        #     ampl.var["h"][i] = val
        # for (i,j) in self.arcs:
        #     ampl.var["y"][i, j] = 10.67*self.L[i,j]/(self.R_max**1.852 * self.d_max**4.87)

        # if self.data_number == 6:
        #     ampl.eval("subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
        # else:
        #     ampl.eval("subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];")

        ampl.option['solver'] = 'ipopt'
        # ampl.option['solver'] = 'baron'
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "expect_infeasible_problem=no "
            "bound_relax_factor=0 "
            "tol=1e-9 "
            f"bound_push={self.bound_push} "
            f"bound_frac={self.bound_frac} "
            "warm_start_init_point=no "
            "halt_on_ampl_error=yes "
            "max_iter=3000"
        )

        ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 barstats version objbound" # lsolver = conopt
        ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600" 
        ampl.option["scip_options"] = "outlev  1 timelimit 3600 lim:gap = 1e-9 chk:feastol = 1e-5 chk:feastolrel=0 " #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        # ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 warmstart = 0 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0 mipgap = 1e-9 NumericFocus = 1" 
        # ampl.eval("option presolve_eps  6.98e-11;")
        # with self.suppress_output():
        ampl.solve()

        solve_result = ampl.get_value("solve_result")

        # lambda_vals = ampl.getVariable('lambda').getValues().toDict()
        # q_p = ampl.getParameter('q_p').getValues().toDict()
        #
        # # Get B
        # B = ampl.getParameter('B').getValues().toDict()
        #
        # # Get arcs and cycles
        # cycles = ampl.getSet('cycles').to_list()
        #
        # # Recover q
        # q_sol = {}
        #
        # for (i, j) in self.arcs:
        #     val = q_p.get((i, j), 0.0)
        #
        #     for c in cycles:
        #         val += B.get((i, j, c), 0.0) * lambda_vals.get(c, 0.0)
        #
        #     q_sol[(i, j)] = val
        #
        # print("Recovered flows q:", q_sol)

        q_sol = ampl.get_variable('q').get_values().to_dict()
        h_sol = ampl.get_variable('h').get_values().to_dict()
        y_sol = ampl.get_variable('y').get_values().to_dict()
        z_sol = ampl.get_variable('z').get_values().to_dict()
        cost = ampl.getObjective("total_cost").value()
        solve_time = ampl.get_value('_solve_elapsed_time')
        # ampl.eval("display z;")
        # ampl.eval("display x;")
        
        self.alpha = ampl.get_parameter('alpha').get_values().to_dict()
        # obj = 0 
        # for (i,j) in self.arcs:
        #     obj = obj + z_sol[i,j]**0.5
        #
        # print("Objective", obj)
        # for (i,j) in self.arcs:
        #     for k in self.pipes:
        #         print(f"l[{i},{j},{k}]:", (y_sol[i,j]-alpha[k+1]*self.L[i,j])/alpha[k] - alpha[k+1])

        return ampl, solve_result, z_sol, q_sol, h_sol, y_sol, cost, solve_time


    def solve_with_boundary_pushing_cut(self, z_sol, y_sol, h_sol, max_iterations=30, delta=0.5, seed=0):
        """
        Boundary-pushing heuristic to explore multiple local solutions.

        Parameters
        ----------
        max_iterations : int
            Number of IPOPT resolves with cuts
        delta : float
            Strength of boundary push (larger = more aggressive change)
        seed : int
            Random seed for reproducibility

        Returns
        -------
        best_solution : tuple
            (q_sol, h_sol, y_sol, z_sol, best_cost)
        """
        np.random.seed(seed)

        # -------------------------
        # Initialize AMPL
        # -------------------------
        seg_index = {}
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)  # descending
        for (i, j) in self.arcs:
            y_val = y_sol[(i, j)]

            # --- detect segment ---
            for k in range(len(alpha_vals) - 1):
                if alpha_vals[k+1] <= y_val <= alpha_vals[k]:
                    seg_index[i, j] = k+1
                    break
        # -------------------------
        # Main loop
        # -------------------------
        for it in range(max_iterations):
            print(f"\n--- Iteration {it+1} ---")

            ampl = AMPL()
            ampl.read("exact_reduced_wdn.mod")
            ampl.read_data(self.data_file)

            ampl.option["solver"] = "ipopt"
            # ampl.option["ipopt_options"] = "outlev=0 tol=1e-9 max_iter=3000"
            ampl.option["ipopt_options"] = (
                "outlev=0 "
                "expect_infeasible_problem=no "
                "bound_relax_factor=0 "
                "tol=1e-9 "
                "bound_push=0.1 "
                "bound_frac=0.1 "
                "warm_start_init_point=yes "
                "halt_on_ampl_error=yes "
                "max_iter=3000"
            )

            # -------------------------
            # Storage
            # -------------------------
            # best_cost = float("inf")
            best_solution = self.best_solution
            # y_sol = None
            ampl.eval("""
                param seg_index{(i,j) in arcs};
            """)
            ampl.eval("""
                param y_sol{(i,j) in arcs};
            """)
            ampl.eval("""
                param x_sol{(i,j) in arcs} default 0;
            """)
            ampl.eval("""set 2_set := {1,2};""")
            ampl.eval("""var x{arcs}>=0, <=1;""")
            # -------------------------
            # Warm-start perturbation
            # -------------------------
            if it > 0 and y_sol is not None:
                y_var = ampl.get_variable("y")
                for (i, j) in self.arcs:
                    val = y_sol[(i, j)]
                    perturb = np.random.uniform(0.97, 1.03)
                    y_var[i, j].set_value(val * perturb)

            for (i, j), val in y_sol.items():
                ampl.var["y"][i, j] = val
            for i, val in h_sol.items():
                ampl.var["h"][i] = val
            for (i, j), val in z_sol.items():
                ampl.var["z"][i, j] = val

            # -------------------------
            # Solve NLP
            # -------------------------
            for (i, j) in self.arcs:
                ampl.param["seg_index"][i, j] = seg_index[i,j]
                ampl.param["y_sol"][i, j] = y_sol[(i,j)]

            # ampl.eval(f"""minimize total_cost:sum{{(i,j) in arcs}} (z[i,j] - (y[i,j] - y_sol[i,j])^2);""")
            # ampl.eval(f"""minimize total_cost:sum{{(i,j) in arcs}} z[i,j];""")

            # -------------------------
            # Build boundary cut
            # -------------------------
            expr_terms = []
            current_sum = 0.0
            s_star = 0.0

            for (i, j) in self.arcs:
                y_val = y_sol[(i, j)]

                # --- detect segment ---
                seg_low = alpha_vals[-1]
                seg_high = alpha_vals[0]

                for k in range(len(alpha_vals) - 1):
                    if alpha_vals[k+1] <= y_val <= alpha_vals[k]:
                        seg_low = alpha_vals[k+1]
                        seg_high = alpha_vals[k]
                        # ampl.param["seg_index"][i, j] = k+1
                        s_star += seg_low - y_val
                        break

                # --- normalized coordinate in segment ---
                denom = seg_high - seg_low
                if abs(denom) < 1e-12:
                    denom = 1.0  # avoid division issues

                normalized_val = (y_val - seg_low) / denom
                current_sum += normalized_val

                # --- AMPL expression ---
                term = f"(y[{i},{j}] - {seg_low}) / {denom}"
                expr_terms.append(term)

            # -------------------------
            # Add cut
            # -------------------------
            rhs = current_sum + delta
            cut_name = f"boundary_cut_{it}"
            lhs_expr = " + ".join(expr_terms)

            cut_stmt = f"subject to {cut_name}: {lhs_expr} >= {rhs};"

            try:
                # ampl.eval(cut_stmt)
                # ampl.eval(f"""
                #     subject to pwl_dual_cut{it}:
                #         sum{{(i,j) in arcs}} z[i,j] <= sum{{(i,j) in arcs}} L[i,j] *
                #             (slope[seg_index[i,j]] * y[i,j] + intercept[seg_index[i,j]]) - {0.001};
                # """)
                # ampl.eval(f"""
                #     subject to pwl_cut{it}:
                #         sum{{(i,j) in arcs}} (alpha[seg_index[i,j]] - y[i,j]) <= {s_star};
                # """)
                    # compute distance to nearest breakpoint for each arc

                d_min = float('inf')
                for (i,j) in self.arcs:
                    s      = seg_index[(i,j)]        # active segment index
                    a_up   = self.alpha[s]                # alpha_s   (upper, larger value)
                    a_lo   = self.alpha[s+1]              # alpha_{s+1} (lower, smaller value)
                    y_val  = y_sol[(i,j)]
                    # if y_sol[i,j]>a_lo:
                    dist_to_upper = (a_up - y_val)/(a_up-a_lo)     # >= 0
                    dist_to_lower = (y_val - a_lo)/(a_up-a_lo)     # >= 0
                    d_ij   = min(dist_to_upper, dist_to_lower)
                    if d_ij >= 0.001:
                        d_min  = min(d_min, d_ij)
                        # print("d_min:",d_min)

                rhs = d_min ** 2
                # print("rhs:", rhs)
                ampl.eval(f"""
                    subject to euclidean_cut_{it}:
                        sum{{(i,j) in arcs}} ((y[i,j] - y_sol[i,j])/(alpha[seg_index[i,j]] - alpha[seg_index[i,j]+1]))^2 >= {rhs};
                """)

                # ampl.eval(f"""
                #     subject to pwl_cut{it}:
                #         sum{{(i,j) in arcs}} (y[i,j] - y_sol[i,j])^2 >= min{{(i,j) in arcs}} (min{{alpha[seg_index[i,j]+1] - y_sol[i,j], alpha[seg_index[i,j]] - y_sol[i,j]}})^2;
                # """)
                # continue

                # ampl.eval(f"""
                #     subject to pwl_cut{it}:
                #         sum{{(i,j) in arcs}} (y[i,j]/alpha[seg_index[i,j]+1] + y[i,j]/alpha[seg_index[i,j]]) >= 1 + sum{{(i,j) in arcs}} (y_sol[i,j]/alpha[seg_index[i,j]+1] + y_sol[i,j]/alpha[seg_index[i,j]]);
                # """)
                # ampl.eval(f"""
                #     subject to pwl_cut{it}:
                #         sum{{(i,j) in arcs}} (x[i,j]*(y[i,j] - alpha[seg_index[i,j]])) >= 0;
                # """)

                # if it >=1:
                #     zero_arcs = [(i, j) for (i, j) in self.arcs if abs(x_sol[(i, j)]) <= 1e-3]
                #     nz_arcs   = [(i, j) for (i, j) in self.arcs if abs(x_sol[(i, j)]) > 1e-3]
                #
                #     ampl.eval("set ZERO_ARCS within {arcs};")
                #     ampl.eval("set NZ_ARCS within {arcs};")
                #
                #     ampl.set["ZERO_ARCS"] = zero_arcs
                #     ampl.set["NZ_ARCS"] = nz_arcs
                #
                #     for (i, j) in self.arcs:
                #         ampl.param["x_sol"][i, j] = x_sol[(i,j)]
                #
                #     cut_name = f"cut_{it}"
                #
                #     # ampl.eval(f"""s.t. new_cut:
                #     #        sum {{(i,j) in arcs}} (x[i,j]) >= 0.001;
                #     # """)
                #     ampl.eval(f"""
                #     subject to {cut_name}:
                #         sum{{(i,j) in ZERO_ARCS}} (x[i,j] - x_sol[i,j])
                #       + sum{{(i,j) in NZ_ARCS}} (x[i,j] - x_sol[i,j])
                #       = 1;
                #     """)
                # else:
                #     ampl.eval(f"""s.t. new_cut:
                #            sum {{(i,j) in arcs}} (x[i,j]) = 1;
                #     """)
                # ampl.eval(f"""s.t. new_cut{{(i,j) in arcs}}:
                #             (x[i,j,1] + x[i,j,2]) = 1;
                # """)
                # print(f"Added cut: sum(normalized y) >= {rhs:.4f}")

            except Exception as e:
                print(f"Failed to add cut: {e}")
                break

            with self.suppress_output():
                ampl.solve()

            if ampl.get_value("solve_result") != "solved":
                print("Solver failed. Stopping.")
                break

            # -------------------------
            # Extract solution
            # -------------------------
            # cost = ampl.get_objective("total_cost").value()
            q_sol = ampl.get_variable("q").get_values().to_dict()
            h_sol = ampl.get_variable("h").get_values().to_dict()
            y_sol = ampl.get_variable("y").get_values().to_dict()
            z_sol = ampl.get_variable("z").get_values().to_dict()
            x_sol = ampl.get_variable('x').get_values().to_dict()

            # print("x:", x_sol)
            # ampl.eval("display y;")
            # for (i, j) in self.arcs:
            #     y_val = y_sol[(i, j)]
            #
            #     # --- detect segment ---
            #     for k in range(len(alpha_vals) - 1):
            #         if alpha_vals[k+1] <= y_val <= alpha_vals[k]:
            #             seg_index[i, j] = k+1
            #             break

            cost = sum(z_sol[i,j] for (i,j) in self.arcs)

            print(f"Local cost: {cost:.6f} Best cost: {self.best_cost:.6f}")
            for (i, j) in self.arcs:
                y_val = y_sol[(i, j)]
                # --- detect segment ---
                for k in range(len(alpha_vals) - 1):
                    if alpha_vals[k+1] <= y_val <= alpha_vals[k]:
                        seg_index[i, j] = k+1
                        break
            # -------------------------
            # Update best solution
            # -------------------------
            if cost < self.best_cost:
                print("New best solution found.")
                self.best_cost = cost
                self.best_solution = (q_sol, h_sol, y_sol, z_sol, cost)
        return self.best_solution

    def generate_random_acyclic_from_solution(self, q):
        self.network_graph = nx.DiGraph()
        self.network_graph.add_nodes_from(self.nodes) 
        # q = self.ampl.getVariable('q').getValues().toDict()
        for (i,j) in self.arcs:
            if q[i,j] >= 0:
                self.network_graph.add_edge(i,j)
            else:
                self.network_graph.add_edge(j,i)
        return self.network_graph

    def format_indian_number(self,num):
        num_str = str(num)
        if len(num_str) <= 3:
            return num_str
        else:
            # Split the number into the last three digits and the rest
            last_three = num_str[-3:]
            remaining = num_str[:-3]
            # Add commas every two digits in the remaining part
            remaining = ','.join([remaining[max(i - 2, 0):i] for i in range(len(remaining), 0, -2)][::-1])
            return remaining + ',' + last_three



    def segment_cut_based_heuristic(self):
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]

        print("----------------------------------------------------------------------------------------")
        print("Iteration :",self.iteration)
        print("----------------------------------------------------------------------------------------")
        improved = False

        # --------------------------------------------------
        # Collect duals of all constraints
        # --------------------------------------------------
        self.all_duals = {}
        for con_name, con in self.ampl.get_constraints():
            self.all_duals[con_name] = con.getValues()

        # --------------------------------------------------
        # Build candidate arc list from indegree ≥ 2 nodes
        # --------------------------------------------------
        sorted_all_arcs = []
        for node in self.indegree_2_or_more:
            for (u, v) in self.network_graph.in_edges(node):
                arc = (u, v) if (u, v) in self.arcs else (v, u)
                sorted_all_arcs.append(arc)

        # Remove fixed and already-visited arcs
        # sorted_all_arcs = [
        #     arc for arc in sorted_all_arcs
        #     if arc not in self.fix_arc_set
        # ]

        sorted_arcs = [
            arc for arc in sorted_all_arcs
            if arc not in self.visited_arc_reverse
        ]

        # sorted_arcs = [arc for arc in sorted_arcs if self.seg_index[arc[0],arc[1]]>=2]

        # sorted_arcs = [arc for arc in sorted_arcs if arc_max_dia[arc[0], arc[1]] == 1]
        # --------------------------------------------------
        # Extract duals for con2 only (restricted to sorted_arcs)
        # --------------------------------------------------
        dual_dict = {}
        for con_name, dual_values in self.all_duals.items():
            if con_name == "con2":
                tmp = dual_values.to_dict()
                dual_dict = {
                    arc: val for arc, val in tmp.items()
                    if arc in sorted_arcs
                }
                break   # only con2 is needed

        # --------------------------------------------------
        # Sensitivity score computation
        # sen_score(i,j) = - dual(i,j) * (h[i] - h[j])
        # --------------------------------------------------
        self.sen_score = {}
        # for (i, j), dual_val in dual_dict.items():
        #     self.sen_score[(i, j)] = -dual_val * np.abs(self.h[i] - self.h[j])
        #     # self.sen_score[(i, j)] = -dual_val * np.abs(self.q[i,j])
        #     if self.data_number==5:
        #         self.sen_score[(i, j)] = -dual_val *(2 + 1.852*np.abs(self.q[i,j])**0.852 * sum(10.67*self.l[i,j,k]/(float(self.R[i,j])**1.852 * self.d[k]**4.87) for k in self.pipes) + np.abs(self.q[i,j])**1.852 * sum(10.67/(float(self.R[i,j])**1.852 * self.d[k]**4.87) for k in self.pipes) )
        #     else:
        #         # self.sen_score[(i, j)] = -dual_val * 1.852*np.abs(self.q[i,j])**0.852 * sum(10.67*self.l[i,j,k]/(float(self.R[k])**1.852 * self.d[k]**4.87) for k in self.pipes)
        #         self.sen_score[(i, j)] = -dual_val * (2 + 1.852*np.abs(self.q[i,j])**0.852 * sum(10.67*self.l[i,j,k]/(float(self.R[k])**1.852 * self.d[k]**4.87) for k in self.pipes) + sum(10.67/(float(self.R[k])**1.852 * self.d[k]**4.87) for k in self.pipes)*np.abs(self.q[i,j])**1.852)
        #
        # print("sen_score:", self.sen_score)

        # --------------------------------------------------
        # Rank arcs using sensitivity score (absolute value)
        # --------------------------------------------------
        # sorted_arcs = [
        #     arc for arc, _ in sorted(
        #         self.sen_score.items(),
        #         key=lambda kv: abs(kv[1]),
        #         reverse=True
        #     )
        # ]

        sorted_arcs = [
            arc for arc, _ in sorted(
                dual_dict.items(),
                key=lambda kv: abs(kv[1]),
                reverse=True
            )
        ]

        # if self.reversed_arcs:
        # sorted_arcs = [
        #    arc for arc in sorted_arcs
        #    if arc not in self.reversed_arcs
        # ]

        print("sorted_arcs:", sorted_arcs)
        # sorted_arcs = sorted_arcs[:min(10, len(sorted_arcs))]
        # print("reversed_arcs:", self.reversed_arcs)

        print("----------------------------------------------------------------------------------------")
        print(f"{'NLP':<5}{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")
        print("----------------------------------------------------------------------------------------")

        # for edge in sorted_arcs[:min(20, len(sorted_arcs))]:
        # for edge in sorted_arcs:

        while sorted_arcs:
            (i,j) = sorted_arcs[0]
            sorted_arcs.remove((i,j))
            edge = (i,j)
            self.visited_arc_reverse.append(edge)

            ampl = AMPL()
            ampl.reset()
            if self.data_number==5:
                ampl.read("newyork_model.mod")
            elif self.data_number==6:
                ampl.read("blacksburg_model.mod")
            else:
                ampl.read("exact_reduced_wdn.mod")
            ampl.read_data(self.data_file)

            # ampl.eval("""param seg_index{(i,j) in arcs};""")
            # ampl.eval("""param y_sol{(i,j) in arcs};""")

            # for (u, v) in self.arcs:
            #     ampl.param["seg_index"][u, v] = seg_index[u,v]
            #     ampl.param["y_sol"][u, v] = self.y[(u,v)]

            for (u,v), val in self.y.items():
                ampl.eval(f'let y[{u},{v}] := {val};')
            for (u, v), val in self.q.items():
                ampl.eval(f'let q[{u},{v}] := {val};')
                if self.data_number ==5:
                    ampl.eval(f'let q1[{u},{v}] := {self.q1[u,v]};')
                    ampl.eval(f'let q2[{u},{v}] := {self.q2[u,v]};')
            for u, val in self.h.items():
                ampl.eval(f'let h[{u}] := {val};') 
            for (u, v), val in self.z.items():
                ampl.eval(f'let z[{u},{v}] := {val};')

            current_duals = {}
            for con_name, val in ampl.get_constraints():
                dual_values = val.get_values()
                current_duals[con_name] = dual_values

            # Initialize dual values for all constraints
            for con_name, dual_values in self.all_duals.items():
                if con_name in current_duals:
                    # Initialize dual values for each constraint
                    ampl.get_constraint(con_name).set_values(dual_values)               

            # add next segment bound on y variable
            ampl.eval(f"""s.t. seg_cut: alpha[{self.seg_index[i,j]}] <= y[{i},{j}];""")

            ampl.option["solver"] = "ipopt"
            ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes")   #max_iter = 1000
            
            # ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes mu_strategy = adaptive recalc_y = no"
            # ampl.set_option("snopt_options", "outlev = 2")
            ampl.option["presolve_eps"] = "6.82e-14"
            ampl.option['presolve'] = 1
            with self.suppress_output():
                ampl.solve()
            #self.ampl.eval("display q;")
            self.solve_result = ampl.solve_result
            self.total_cost = ampl.get_objective("total_cost").value()
            solve_time = ampl.get_value('_solve_elapsed_time')
            self.solver_time += solve_time
            self.number_of_nlp += 1
            if self.solve_result == "solved":
                y = ampl.getVariable('y').getValues().to_dict()
                z = ampl.getVariable('z').getValues().to_dict()
                q = ampl.getVariable('q').getValues().to_dict()
                h = ampl.getVariable('h').getValues().to_dict() 

                if self.total_cost < self.current_cost:
                    print(
                        f"{self.number_of_nlp:<5}"
                        f"{str((i, j)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(self.ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                        f"{self.solve_result:<14}{'Yes':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")

                    self.current_cost = self.total_cost
                    improved = True
                    # self.is_improved_in_arc_reversal = True
                    self.ampl = ampl
                    self.network_graph = self.generate_random_acyclic_from_solution(q)

                    self.y = y 
                    self.z = z 
                    self.q = q
                    self.h = h 

                    # print(self.q)
                    # ampl.eval("display q;")

                    if self.data_number==5:
                        self.q1 = ampl.getVariable('q1').getValues().to_dict()
                        self.q2 = ampl.getVariable('q2').getValues().to_dict()
                    # self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set)
                    for con_name, con in self.ampl.get_constraints():
                        self.all_duals[con_name] = con.getValues()
                    
                    self.seg_index = {}
                    
                    # Sort alpha values in descending order
                    alpha_vals = sorted(set(self.alpha.values()), reverse=True)
                    
                    tol = 1e-8
                    
                    for (u, v) in self.arcs:
                        y_val = self.y[(u, v)]
                        active_s = None
                    
                        for k in range(len(alpha_vals) - 1):
                            a_high = alpha_vals[k]
                            a_low  = alpha_vals[k+1]
                    
                            if (a_low - tol) <= y_val <= (a_high + tol):
                                active_s = k + 1   # 1-based segment index (matches AMPL segs)
                                break
                    
                        # fallback (numerical edge cases)
                        if active_s is None:
                            if y_val > alpha_vals[0]:
                                active_s = 1
                            elif y_val < alpha_vals[-1]:
                                active_s = len(alpha_vals) - 1
                    
                        self.seg_index[(u, v)] = active_s
                    
                        # print(f"arc {(u,v)} -> y = {y_val:.6f}, segment = {active_s}")

                    # self.seg_index = {}
                    # alpha_vals = sorted(set(self.alpha.values()), reverse=True)  # descending
                    # print("alpha:",self.alpha)

                    # dual_vals = self.all_duals["exact_cost"].to_dict()
                    # print("dual_vals:", dual_vals)
                    # for (u, v) in self.arcs:
                    #     # max_dual = 0
                    #     active_s = None
                    #
                    #     for s in self.segs:
                    #         dual = dual_vals[u,v,s]
                    #         print((u,v),dual)
                    #         if dual > 1e-8:
                    #             # max_dual = dual
                    #             active_s = s
                    #             break
                    #
                    #     self.seg_index[(u, v)] = active_s
                        # print(f"arc {(u,v)} -> active segment: {active_s}, dual={dual}") 


                    # for (u,v) in sorted_arcs:
                    #     if q[u,v]*self.q[u,v]<0:
                    #         # print("reversed arc:", (u,v), "\n")
                    #         if (u,v) not in self.reversed_arcs:
                    #             self.reversed_arcs.append((u,v))
                    #             sorted_arcs.remove((u,v))

                else: 
                    print(
                        f"{self.number_of_nlp:<5}"
                        f"{str((i, j)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                        f"{self.solve_result:<14}{'No':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    # self.reversed_arcs = []
                    # for (x,y) in sorted_arcs:
                    #     if q[x,y]*self.q[x,y]<0:
                    #         # print("reversed arc:", (x,y), "\n")
                    #         if (x,y) not in self.reversed_arcs:
                    #         # if (u,v) in sorted_arcs:
                    #             self.reversed_arcs.append((x,y))
                    #             sorted_arcs.remove((x,y))
            else:
                print(
                    f"{self.number_of_nlp:<5}"
                    f"{str((i, j)):<10}"
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{self.format_indian_number(round(self.total_cost)):<14}"
                    f"{(str(round(ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                    f"{self.solve_result:<14}{'No':<10}"
                    f"{round(time.time() - self.start_time, 2)}s")
            if improved:
                self.iteration = self.iteration + 1
                self.segment_cut_based_heuristic()
                break

    @contextlib.contextmanager
    def suppress_output(self):
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
    # ============================================================
    # RUN
    # ============================================================
    def run(self):
        self.start_time = time.time()
        self.read_model_and_data()

        fix_arc_set = self.fix_leaf_arc_flow(self.ampl)
        # print("fix_arc_set:",fix_arc_set)
        self.super_source_out_arc = self.fix_arc_set()
        # print("super_source_out_arc:", self.super_source_out_arc, "\n")
        self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set)        
        # print("fix_arc_set:",self.fix_arc_set)

        self.sorted_arcs = [
            arc for arc in self.arcs
            if arc not in self.fix_arc_set
        ]
        # self.cycles = self.cycle_basis() 
        # print("sorted_arcs:", self.sorted_arcs)

        # self.plot_value_function_for_arc((1, 2))
        # self.plot_value_function_for_arcs([(19,3), (18,19)])
        # self.plot_m_values()
        # self.plot_piecewise_lines((1,2))

        # print("\n-------------------------------- Solving Original Model --------------------------")

        # self.solve_original_model_without_init()
        
        # print("-------------REDUCED MODEL FEASIBILITY CHECK------------------")
        # y, z = self.map_original_to_reduced(self.l, self.q, self.h)
        # viol_red = self.check_reduced_feasibility(self.q, self.h, y, z)
        # print(viol_red)

        # self.print_piecewise_reduced_cost(arc=(1, 2))
        # Plot segment-wise piecewise function
        # self.plot_piecewise_reduced_cost(
        #     arc=(1, 2),
        #     save_path="piecewise_reduced_cost_arc_1_2.png"
        # )

        # Plot continuous envelope
        # self.plot_piecewise_reduced_cost_envelope(
        #     arc=(1, 2),
        #     save_path="piecewise_reduced_cost_envelope_arc_1_2.png"
        # )

        print("\n-------------------------------- Solving Exact Reduced Model --------------------------")

        ampl, solve_result, z_sol, q_sol, h_sol, y_sol, cost, solve_time = self.solve_exact_reduced_model()

        print(f"Total Cost: {cost:.8f}")
        print(f"Exact model solve time: {solve_time:.4f} sec")
        # ampl.eval("display exact_cost.dual;")
        self.ampl = ampl
        self.q = q_sol
        self.h = h_sol
        self.y = y_sol
        self.z = z_sol
        self.best_cost = cost
        self.current_cost = cost
        self.best_solution = (self.q, self.h, self.y, self.z, self.best_cost)

        self.all_duals = {}
        for con_name, con in self.ampl.get_constraints():
            self.all_duals[con_name] = con.getValues()

        self.seg_index = {}
        
        # Sort alpha values in descending order
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        
        tol = 1e-8
        
        for (u, v) in self.arcs:
            y_val = self.y[(u, v)]
            active_s = None
        
            for k in range(len(alpha_vals) - 1):
                a_high = alpha_vals[k]
                a_low  = alpha_vals[k+1]
        
                if (a_low - tol) <= y_val <= (a_high + tol):
                    active_s = k + 1   # 1-based segment index (matches AMPL segs)
                    break
        
            # fallback (numerical edge cases)
            if active_s is None:
                if y_val > alpha_vals[0]:
                    active_s = 1
                elif y_val < alpha_vals[-1]:
                    active_s = len(alpha_vals) - 1
        
            self.seg_index[(u, v)] = active_s
        
            # print(f"arc {(u,v)} -> y = {y_val:.6f}, segment = {active_s}")


        # self.seg_index = {}
        # alpha_vals = sorted(set(self.alpha.values()), reverse=True)  # descending
        # print("alpha:",self.alpha)

        # dual_vals = self.all_duals["exact_cost"].to_dict()
        # print("dual_vals:", dual_vals)
        # for (u, v) in self.arcs:
        #     # max_dual = 0
        #     active_s = None
        #
        #     for s in self.segs:
        #         dual = dual_vals[u,v,s]
        #
        #         if dual > 1e-6:
        #             # max_dual = dual
        #             active_s = s
        #             break
        #
        #     self.seg_index[(u, v)] = active_s
            # print(f"arc {(u,v)} -> active segment: {active_s}, dual={dual}") 

        # print("\n-------------------------------- Segment Cut Based Heuristic --------------------------")

        # best_solution = self.solve_with_boundary_pushing_cut(z_sol, y_sol, h_sol) 
        # q_sol, h_sol, y_sol, z_sol, cost = best_solution 
        # print("\n Best Solution:", cost)

        print("\n-------------------------------- Branch and Bound based Heuristic --------------------------")
        self.network_graph = self.generate_random_acyclic_from_solution(q_sol)
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        self.best_acyclic_flow = self.network_graph.copy() 

        self.iteration = 1
        self.visited_arc_reverse = []
        self.reversed_arcs = []
        self.segment_cut_based_heuristic()

        for (i,j) in self.arcs:
            print((i,j), self.z[i,j] - self.L[i,j]*(self.slope[self.seg_index[i,j]]*self.y[i,j] + self.intercept[self.seg_index[i,j]]))

        print("\n-------------------------------- Solving Recover Model --------------------------")

        rec_ampl, rec_result, l_trial, recovered_cost, t_rec = self.solve_recover_model1(self.y)
        print(rec_result)
        print(recovered_cost)

        # print(l_trial)
        # for (i,j) in self.arcs:
        #     for k in self.pipes:
        #         if l_trial[i,j,k]>0.0001:
        #             print(f"l[{i},{j},{k}]:", l_trial[i,j,k])
        # rec_ampl.eval("display l;")
        
        # rec_ampl, rec_result, l_trial, recovered_cost, t_rec = self.solve_recover_model(self.y)
        # print(rec_result)
        # print(recovered_cost)

        # print(l_trial)
        # for (i,j) in self.arcs:
        #     for k in self.pipes:
        #         if l_trial[i,j,k]>0.001:
        #             print(f"l[{i},{j},{k}]:", l_trial[i,j,k])

        # print("-------------ORIGINAL MODEL FEASIBILITY CHECK------------------")
        # Reduced → Original
        # l_rec = self.map_reduced_to_original(y_sol)

        # viol_org = self.check_original_feasibility(q_sol, h_sol, l_trial)
        # print(viol_org)
        #
        print("\n-------------------- Solving Original Model with Initialize ---------------------")

        orig_ampl, orig_result, q_new, h_new, l_new, final_cost, t_orig = self.solve_original_with_init(l_trial, self.q, self.h)
        print(f"Original warm-start NLP result = {orig_result}, final cost = {final_cost:.8f}, time = {t_orig:.2f} sec")

        # self.improve_solution_by_y_search(max_outer_iter=1,max_repair_iter=1,num_arcs_perturb=3,perturb_scale=0.03)




if __name__ == "__main__":
    model = sys.argv[1]

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

    data_number = int(sys.argv[3]) - 1
    data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{data_list[data_number]}.dat"

    print("Water Network:", f"{sys.argv[3]}")
    print("***********************************************************************************************")

    solver = sys.argv[2]

    solver_instance = WaterNetworkSolver(model, solver, data, data_number)
    solver_instance.run()

