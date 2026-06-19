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
        self.bound_push = 0.01
        self.bound_frac = 0.01

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
        self.alpha_min = min(self.alpha[k] for k in self.segs)
        self.alpha_max = max(self.alpha[k] for k in self.segs)

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
        ampl.option['presolve_eps'] = 2.04e-10 
        # ampl.option["ipopt_options"] = (
        #     "outlev=0 "
        #     "expect_infeasible_problem=no "
        #     "bound_relax_factor=0 "
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
        y_min = max(0, np.min(alpha_vals))
        y_max = np.max(alpha_vals)

        y_grid = np.linspace(y_min, y_max, n_points)

        #plt.figure(figsize=(15, 6))
        plt.figure(figsize=(18, 7))
        print("m",m)
        # 🔥 Plot all lines
        for k in m:
            print(k)
            z_vals = np.array([Lij*(m[k]*y + b[k]) for y in y_grid])
            # z_vals = np.maximum(z_vals, 0)
            # plt.plot(y_grid, z_vals, linestyle='-', label=f"line k={k-1}")
            # plt.plot(y_grid,z_vals,linestyle='-',label=rf"$z(y) = {m[k]:.3g}\,y + {b[k]:.3g}\,L$")
            # plt.plot(y_grid,z_vals,linestyle='-',label=rf"$z_{{ij}}(y_{{ij}}) = L(m_{{{k-1}}} y_{{ij}} + b_{{{k-1}}})$")
            plt.plot(
                y_grid,
                z_vals,
                linestyle='-',
                linewidth=2.5,
                label=rf"$z_{{ij}}(y_{{ij}})=L_{{ij}}(m_{{{k-1}}}y_{{ij}}+b_{{{k-1}}})$"
            )

        # 🔥 Plot envelope (min of lines)
        z_env = []
        for y in y_grid:
            z_env.append(max(Lij*(m[k]*y + b[k]) for k in m))
        # z_env = np.maximum(z_env)

        # plt.plot(y_grid, z_env, color='black',linestyle='--', linewidth=1, label="Convex Upper Envelope (Max-affine)")
        
        plt.plot(
            y_grid,
            z_env,
            color='black',
            linestyle='--',
            linewidth=3,
            label="Convex Upper Envelope"
        )
        # 🔥 Vertical lines at y = alpha_k * Lij
        for k in pipe_list:
            yk = alpha[k]

            if yk >= 0:
                # plt.axvline(x=yk, linestyle=':', linewidth=1.5, color="gray")
                plt.axvline(
                    x=yk,
                    linestyle=':',
                    linewidth=2.5,
                    color="gray"
                )
                # plt.text(yk, -plt.ylim()[1]*0.01, rf"$\alpha_{{{k}}}$", rotation=0, horizontalalignment='center', fontsize=8)
                plt.text(
                    yk,
                    -plt.ylim()[1]*0.03,
                    rf"$\alpha_{{{k}}}$",
                    rotation=0,
                    horizontalalignment='center',
                    fontsize=18
                )


        # 🔴 Plot only breakpoint points
        for k in pipe_list:
            xk = alpha[k]
            zk = C[k] * Lij
        
            # plt.scatter(xk, zk, s=50)
            plt.scatter(xk, zk, s=140)
            plt.text(
                xk,
                zk,
                rf"$(\alpha_{{{k}}},\;C_{{{k}}}L_{{ij}})$",
                # fontsize=8,
                fontsize=18,
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
        
        # plt.xlabel(r"$y_{ij}$")
        # plt.ylabel(r"$z_{ij}(y_{ij})$")
        plt.xlabel(r"$y_{ij}$", fontsize=26)
        plt.ylabel(r"$z_{ij}(y_{ij})$", fontsize=26)
        # plt.title(rf"Convex piecewise linear cost function $z_{{ij}}(y_{{ij}})$ for Arc $(i,j)$")

        # plt.legend()
        plt.legend(
            fontsize=18,
            loc='best',
            frameon=True
        )
        # plt.grid(True)
        plt.savefig("convex_pwl_cost.png", dpi=300, bbox_inches="tight")
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
            b_k = (self.C[kp1] - m_k * alpha_kp1)
    
            # feasible y interval for this segment
            y_low = alpha_kp1
            y_high = alpha_k
    
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
    
        Lij = self.L[(i, j)]
        # plot each segment
        for s in segments:
            y_vals = np.linspace(s["y_low"], s["y_high"], 100)
            z_vals = Lij * (s["m"] * y_vals + s["b"])
    
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
        # ampl.read("dc_reformulation.mod")
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
            "bound_relax_factor=0 "
            f"bound_push={self.bound_push} "
            f"bound_frac={self.bound_frac} "
            "warm_start_init_point=no "
            "halt_on_ampl_error=yes "
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

    def lagrangian_based_heuristic2(
        self,
        q_star,
        h_star,
        y_star,
        z_star,
        seg_index,
        slack_penalty=1.0,
        slack_relax=500,
        max_outer=30,
        top_k_fraction=0.25,
        cost_improve_tol=1e-4,
        segment_change_tol=2
    ):
        """
        Active epigraph slack-relaxation heuristic.
    
        Main ideas
        ----------
        1. Detect active epigraph segments.
        2. Relax only expensive active segments.
        3. Add slack variables directly in constraints.
        4. Penalize slacks mildly in objective.
        5. Force IPOPT into different active-set basin.
        6. Project relaxed solution back to exact NLP.
        """
    
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        NP = len(alpha_vals)
    
        # ============================================================
        # utilities
        # ============================================================
    
        def find_seg(yv):
    
            for k in range(NP - 1):
    
                if alpha_vals[k + 1] - 1e-10 <= yv <= alpha_vals[k] + 1e-10:
                    return k + 1
    
            return 1 if yv > alpha_vals[0] else NP - 1
    
        def true_cost(y_sol, seg_sol):
    
            return sum(
                self.L[i, j]
                * (
                    self.slope[seg_sol[i, j]] * y_sol[i, j]
                    + self.intercept[seg_sol[i, j]]
                )
                for (i, j) in self.arcs
            )
    
        def get_expensive_arcs(z_sol):
    
            arc_cost = sorted(
                [
                    ((i, j), z_sol[i, j])
                    for (i, j) in self.arcs
                ],
                key=lambda x: -x[1]
            )
    
            k = max(1, int(top_k_fraction * len(self.arcs)))
    
            return set(a for a, _ in arc_cost[:k])
    
        # ============================================================
        # initial incumbent
        # ============================================================
    
        C_star = true_cost(y_star, seg_index)
    
        incumbent = {
            "q": q_star.copy(),
            "h": h_star.copy(),
            "y": y_star.copy(),
            "z": z_star.copy(),
            "seg": seg_index.copy(),
            "cost": C_star
        }
    
        best_global = incumbent.copy()
    
        print("\n" + "=" * 75)
        print(" ACTIVE-CONSTRAINT SLACK RELAXATION HEURISTIC")
        print(f" INITIAL COST = {C_star:.6f}")
        print("=" * 75)
    
        stagnation = 0
    
        # ============================================================
        # outer iterations
        # ============================================================
    
        for outer_it in range(1, max_outer + 1):
    
            print(f"\nOUTER ITERATION {outer_it}")
            print("-" * 75)
    
            y_inc = incumbent["y"]
            seg_inc = incumbent["seg"]
    
            # ========================================================
            # STEP 1:
            # identify expensive arcs
            # ========================================================
    
            expensive_arcs = get_expensive_arcs(incumbent["z"])
    
            print(f"Expensive arcs selected: {len(expensive_arcs)}")
    
            # ========================================================
            # STEP 2:
            # build relaxed NLP
            # ========================================================
    
            ampl = AMPL()
    
            ampl.read("exact_reduced_wdn.mod")
            ampl.read_data(self.data_file)
    
            ampl.option["solver"] = "ipopt"
    
            ampl.option["ipopt_options"] = (
                "outlev=0 "
                "max_iter=3000 "
                "warm_start_init_point=yes "
                "bound_relax_factor=0 "
            )
    
            # ========================================================
            # remove original objective
            # ========================================================
    
            ampl.eval("drop total_cost;")
    
            # ========================================================
            # add slack variables
            # ========================================================
    
            ampl.eval("""
                var relax_slack{arcs,segs} >= 0;
            """)
    
            # ========================================================
            # relaxed objective
            # ========================================================
    
            ampl.eval(f"""
    
                minimize relaxed_obj:
    
                    sum{{(i,j) in arcs}}
                        z[i,j]
    
                    +
    
                    {slack_penalty}
    
                    *
    
                    sum{{(i,j) in arcs, s in segs}}
                        relax_slack[i,j,s];
    
            """)
    
            # ========================================================
            # remove original epigraph constraints
            # ========================================================
    
            ampl.eval("drop exact_cost;")
    
            # ========================================================
            # rebuild epigraph constraints
            # ========================================================
    
            for (i, j) in self.arcs:
    
                s_cur = seg_inc[(i, j)]
    
                for s in self.segs:
    
                    # ------------------------------------------------
                    # relax only current active expensive arcs
                    # ------------------------------------------------
    
                    # if (s == s_cur):
                    if ((i, j) in expensive_arcs) and (s == s_cur):
    
                        ampl.eval(f"""
    
                            s.t. relaxed_epi_{i}_{j}_{s}:
    
                                z[{i},{j}]
                                +
                                relax_slack[{i},{j},{s}]
    
                                >=
    
                                L[{i},{j}]
                                *
                                (
                                    slope[{s}] * y[{i},{j}]
                                    +
                                    intercept[{s}]
                                );
    
                        """)
    
                    # ------------------------------------------------
                    # keep all other cuts exact
                    # ------------------------------------------------
    
                    else:
    
                        ampl.eval(f"""
    
                            s.t. epi_{i}_{j}_{s}:
    
                                z[{i},{j}]
    
                                >=
    
                                L[{i},{j}]
                                *
                                (
                                    slope[{s}] * y[{i},{j}]
                                    +
                                    intercept[{s}]
                                );
    
                        """)
    
            # ========================================================
            # warm start
            # ========================================================
    
            for (i, j) in self.arcs:
    
                ampl.var["q"][i, j] = incumbent["q"][(i, j)]
    
                ampl.var["y"][i, j] = incumbent["y"][(i, j)]
    
                ampl.var["z"][i, j] = incumbent["z"][(i, j)]
    
            for i in self.nodes:
    
                ampl.var["h"][i] = incumbent["h"][i]
    
            # ========================================================
            # solve relaxed NLP
            # ========================================================
    
            with self.suppress_output():
                ampl.solve()
    
            status = ampl.get_value("solve_result")
    
            if status != "solved":
    
                print("Relaxed NLP failed.")
    
                slack_penalty *= 0.7
    
                continue
    
            # ========================================================
            # retrieve relaxed solution
            # ========================================================
    
            q_rel = ampl.get_variable("q").get_values().to_dict()
    
            y_rel = ampl.get_variable("y").get_values().to_dict()
    
            h_rel = ampl.get_variable("h").get_values().to_dict()
    
            z_rel = ampl.get_variable("z").get_values().to_dict()
    
            seg_rel = {
                (i, j): find_seg(y_rel[(i, j)])
                for (i, j) in self.arcs
            }
    
            changed_arcs = [
    
                (i, j)
    
                for (i, j) in self.arcs
    
                if seg_rel[(i, j)] != seg_inc[(i, j)]
    
            ]
    
            seg_changes = len(changed_arcs)
    
            print(f"Segment changes = {seg_changes}")
    
            # ========================================================
            # insufficient diversification
            # ========================================================
    
            if seg_changes < segment_change_tol:
    
                slack_penalty *= 0.8
    
                stagnation += 1
    
                print("Insufficient basin escape.")
    
                continue
    
            # ========================================================
            # STEP 3:
            # exact NLP projection
            # ========================================================
    
            ampl_ex = AMPL()
    
            ampl_ex.read("exact_reduced_wdn.mod")
            ampl_ex.read_data(self.data_file)
    
            ampl_ex.option["solver"] = "ipopt"
    
            ampl_ex.option["ipopt_options"] = (
                "outlev=0 "
                "warm_start_init_point=yes "
                "bound_relax_factor=0 "
            )
    
            # ========================================================
            # warm start exact model
            # ========================================================
    
            for (i, j) in self.arcs:
    
                ampl_ex.var["q"][i, j] = q_rel[(i, j)]
    
                ampl_ex.var["y"][i, j] = y_rel[(i, j)]
    
                ampl_ex.var["z"][i, j] = z_rel[(i, j)]
    
            for i in self.nodes:
    
                ampl_ex.var["h"][i] = h_rel[i]
    
            # ========================================================
            # exact solve
            # ========================================================
    
            with self.suppress_output():
                ampl_ex.solve()
    
            status_ex = ampl_ex.get_value("solve_result")
    
            if status_ex != "solved":
    
                print("Exact projection failed.")
    
                slack_penalty *= 1.2
    
                continue
    
            # ========================================================
            # retrieve exact solution
            # ========================================================
    
            q_ex = ampl_ex.get_variable("q").get_values().to_dict()
    
            y_ex = ampl_ex.get_variable("y").get_values().to_dict()
    
            h_ex = ampl_ex.get_variable("h").get_values().to_dict()
    
            z_ex = ampl_ex.get_variable("z").get_values().to_dict()
    
            seg_ex = {
                (i, j): find_seg(y_ex[(i, j)])
                for (i, j) in self.arcs
            }
    
            cost_ex = true_cost(y_ex, seg_ex)
    
            print(
                f"Projected exact cost = {cost_ex:.6f} "
                f"Best cost = {best_global['cost']:.6f}"
            )
    
            # ========================================================
            # accept/reject
            # ========================================================
    
            if cost_ex < incumbent["cost"] - cost_improve_tol:
    
                print("*** IMPROVED LOCAL SOLUTION FOUND ***")
    
                incumbent = {
                    "q": q_ex.copy(),
                    "h": h_ex.copy(),
                    "y": y_ex.copy(),
                    "z": z_ex.copy(),
                    "seg": seg_ex.copy(),
                    "cost": cost_ex
                }
    
                if cost_ex < best_global["cost"]:
    
                    best_global = incumbent.copy()
    
                slack_penalty *= 1.1
    
                stagnation = 0
    
            else:
    
                print("No improvement.")
    
                slack_penalty *= 0.9
    
                stagnation += 1
    
            # ========================================================
            # strong diversification
            # ========================================================
    
            if stagnation >= 5:
    
                print("Strong diversification activated.")
    
                slack_penalty *= 0.5
    
                stagnation = 0
    
        # ============================================================
        # final
        # ============================================================
    
        print("\n" + "=" * 75)
        print(f" FINAL BEST COST = {best_global['cost']:.6f}")
        print("=" * 75)
    
        return (
            best_global["q"],
            best_global["h"],
            best_global["y"],
            best_global["z"]
        )


    def bilevel_segment_optimization(
        self,
        q_start,
        h_start,
        y_start,
        z_start,
        seg_start,
        max_outer_iterations=20,
        num_structures=15,
        modify_fraction=0.25,
        diversification_fraction=0.40,
        cost_improve_tol=1e-4,
        stagnation_limit=5
    ):
    
        """
        ==============================================================
        BILEVEL SEGMENT-BASED OPTIMIZATION
        ==============================================================
    
        UPPER LEVEL:
            Generate candidate segment structures
    
        LOWER LEVEL:
            Solve exact hydraulic NLP for each structure
    
        ==============================================================
    
        """
    
        import copy
        import random
    
        from amplpy import AMPL
    
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
    
        NP = len(alpha_vals)
    
        # ============================================================
        # utilities
        # ============================================================
    
        def find_seg(yv):
    
            tol = 1e-8
    
            active_segments = []
    
            for k in range(NP - 1):
    
                a_high = alpha_vals[k]
                a_low  = alpha_vals[k + 1]
    
                # ----------------------------------------------------
                # interior
                # ----------------------------------------------------
    
                if (a_low + tol) < yv < (a_high - tol):
    
                    active_segments = [k + 1]
    
                    break
    
                # ----------------------------------------------------
                # upper breakpoint
                # ----------------------------------------------------
    
                elif abs(yv - a_high) <= tol:
    
                    if k == 0:
                        active_segments = [1]
                    else:
                        active_segments = [k, k + 1]
    
                    break
    
                # ----------------------------------------------------
                # lower breakpoint
                # ----------------------------------------------------
    
                elif abs(yv - a_low) <= tol:
    
                    if k == NP - 2:
                        active_segments = [k + 1]
                    else:
                        active_segments = [k + 1, k + 2]
    
                    break
    
            # --------------------------------------------------------
            # safeguard
            # --------------------------------------------------------
    
            if not active_segments:
    
                if yv > alpha_vals[0]:
    
                    active_segments = [1]
    
                elif yv < alpha_vals[-1]:
    
                    active_segments = [NP - 1]
    
                else:
    
                    distances = [
                        abs(yv - alpha_vals[k])
                        for k in range(NP)
                    ]
    
                    idx = distances.index(min(distances))
    
                    if idx == 0:
    
                        active_segments = [1]
    
                    elif idx == NP - 1:
    
                        active_segments = [NP - 1]
    
                    else:
    
                        active_segments = [idx, idx + 1]
    
            return active_segments
    
        # ============================================================
        # exact cost
        # ============================================================
    
        def true_cost(y_sol, seg_sol):
    
            total = 0.0
    
            for (i, j) in self.arcs:
    
                s = seg_sol[(i, j)][0]
    
                total += (
                    self.L[(i, j)]
                    * (
                        self.slope[s] * y_sol[(i, j)]
                        + self.intercept[s]
                    )
                )
    
            return total
    
        # ============================================================
        # generate candidate segment structures
        # ============================================================
    
        def generate_segment_structures(current_seg, z_sol):
    
            candidate_structures = []
    
            arcs = list(self.arcs)
    
            # --------------------------------------------------------
            # expensive arcs
            # --------------------------------------------------------
    
            expensive_arcs = sorted(
                arcs,
                key=lambda a: -z_sol[a]
            )
    
            for structure_id in range(num_structures):
    
                new_structure = copy.deepcopy(current_seg)
    
                # ----------------------------------------------------
                # choose arcs to modify
                # ----------------------------------------------------
    
                num_modify = max(
                    1,
                    int(modify_fraction * len(arcs))
                )
    
                # ----------------------------------------------------
                # mix:
                # expensive + random
                # ----------------------------------------------------
    
                expensive_count = max(1, num_modify // 2)
    
                expensive_subset = expensive_arcs[:expensive_count]
    
                remaining = [
                    a for a in arcs
                    if a not in expensive_subset
                ]
    
                random_subset = random.sample(
                    remaining,
                    min(len(remaining), num_modify - expensive_count)
                )
    
                arcs_to_modify = (
                    expensive_subset + random_subset
                )
    
                # ----------------------------------------------------
                # mutate arcs
                # ----------------------------------------------------
    
                for arc in arcs_to_modify:
    
                    current_active = current_seg[arc]
    
                    s_cur = current_active[0]
    
                    # ------------------------------------------------
                    # diversification
                    # ------------------------------------------------
    
                    if random.random() <= diversification_fraction:
    
                        far_candidates = [
                            s for s in range(1, NP)
                            if abs(s - s_cur)
                            >= max(2, NP // 3)
                        ]
    
                        if not far_candidates:
    
                            far_candidates = [
                                s for s in range(1, NP)
                                if s != s_cur
                            ]
    
                        s_new = random.choice(far_candidates)
    
                    # ------------------------------------------------
                    # local mutation
                    # ------------------------------------------------
    
                    else:
    
                        step = random.choice([-1, 1])
    
                        s_new = s_cur + step
    
                        s_new = max(
                            1,
                            min(NP - 1, s_new)
                        )
    
                    # ------------------------------------------------
                    # breakpoint structure
                    # ------------------------------------------------
    
                    if random.random() <= 0.25:
    
                        if s_new < NP - 1:
    
                            new_structure[arc] = [
                                s_new,
                                s_new + 1
                            ]
    
                        else:
    
                            new_structure[arc] = [s_new]
    
                    else:
    
                        new_structure[arc] = [s_new]
    
                candidate_structures.append(
                    copy.deepcopy(new_structure)
                )
    
            return candidate_structures
    
        # ============================================================
        # LOWER LEVEL NLP SOLVER
        # ============================================================
    
        def solve_lower_level(
            seg_structure,
            q_warm,
            h_warm,
            y_warm,
            z_warm
        ):
    
            ampl = AMPL()
    
            ampl.read("exact_reduced_wdn.mod")
    
            ampl.read_data(self.data_file)
    
            ampl.option["solver"] = "ipopt"
    
            ampl.option["ipopt_options"] = (
                "outlev=0 "
                "max_iter=5000 "
                "warm_start_init_point=yes "
                "bound_relax_factor=0 "
            )
    
            # ========================================================
            # enforce segment structure
            # ========================================================
    
            for (i, j) in self.arcs:
    
                active_set = seg_structure[(i, j)]
    
                # ----------------------------------------------------
                # single segment
                # ----------------------------------------------------
    
                if len(active_set) == 1:
    
                    s = active_set[0]
    
                    y_upper = alpha_vals[s - 1]
    
                    y_lower = alpha_vals[s]
    
                    ampl.eval(f"""
    
                        s.t. seg_lb_{i}_{j}:
    
                            y[{i},{j}] >= {y_lower};
    
                        s.t. seg_ub_{i}_{j}:
    
                            y[{i},{j}] <= {y_upper};
    
                    """)
    
                # ----------------------------------------------------
                # breakpoint structure
                # ----------------------------------------------------
    
                else:
    
                    s1 = active_set[0]
    
                    s2 = active_set[1]
    
                    y_upper = alpha_vals[min(s1, s2) - 1]
    
                    y_lower = alpha_vals[max(s1, s2)]
    
                    ampl.eval(f"""
    
                        s.t. seg_lb_{i}_{j}:
    
                            y[{i},{j}] >= {y_lower};
    
                        s.t. seg_ub_{i}_{j}:
    
                            y[{i},{j}] <= {y_upper};
    
                    """)
    
            # ========================================================
            # warm start
            # ========================================================
    
            for (i, j) in self.arcs:
    
                ampl.var["q"][i, j] = q_warm[(i, j)]
    
                ampl.var["y"][i, j] = y_warm[(i, j)]
    
                ampl.var["z"][i, j] = z_warm[(i, j)]
    
            for i in self.nodes:
    
                ampl.var["h"][i] = h_warm[i]
    
            # ========================================================
            # solve NLP
            # ========================================================
    
            with self.suppress_output():
    
                ampl.solve()
    
            status = ampl.get_value("solve_result")
    
            if status != "solved":
    
                return None
    
            # ========================================================
            # retrieve solution
            # ========================================================
    
            q_sol = ampl.get_variable("q").get_values().to_dict()
    
            y_sol = ampl.get_variable("y").get_values().to_dict()
    
            h_sol = ampl.get_variable("h").get_values().to_dict()
    
            z_sol = ampl.get_variable("z").get_values().to_dict()
    
            seg_sol = {
                (i, j): find_seg(y_sol[(i, j)])
                for (i, j) in self.arcs
            }
    
            cost_sol = true_cost(y_sol, seg_sol)
    
            return {
                "q": q_sol,
                "y": y_sol,
                "h": h_sol,
                "z": z_sol,
                "seg": seg_sol,
                "cost": cost_sol
            }
    
        # ============================================================
        # INITIAL INCUMBENT
        # ============================================================
    
        initial_cost = true_cost(y_start, seg_start)
    
        incumbent = {
            "q": q_start.copy(),
            "h": h_start.copy(),
            "y": y_start.copy(),
            "z": z_start.copy(),
            "seg": copy.deepcopy(seg_start),
            "cost": initial_cost
        }
    
        best_global = copy.deepcopy(incumbent)
    
        print("\n" + "=" * 75)
    
        print(" BILEVEL SEGMENT-BASED OPTIMIZATION ")
    
        print("=" * 75)
    
        print(f"INITIAL COST = {initial_cost:.6f}")
    
        print("=" * 75)
    
        stagnation = 0
    
        # ============================================================
        # OUTER LEVEL LOOP
        # ============================================================
    
        for outer_it in range(1, max_outer_iterations + 1):
    
            print(f"\nOUTER ITERATION {outer_it}")
    
            print("-" * 75)
    
            # ========================================================
            # UPPER LEVEL:
            # generate segment structures
            # ========================================================
    
            candidate_structures = generate_segment_structures(
                incumbent["seg"],
                incumbent["z"]
            )
    
            iteration_best = None
    
            # ========================================================
            # LOWER LEVEL:
            # solve NLP for each structure
            # ========================================================
    
            for idx, structure in enumerate(candidate_structures):
    
                print(f"Solving candidate structure {idx+1}")
    
                solution = solve_lower_level(
                    structure,
                    incumbent["q"],
                    incumbent["h"],
                    incumbent["y"],
                    incumbent["z"]
                )
    
                if solution is None:
    
                    print("NLP failed.")
    
                    continue
    
                print(
                    f"Candidate cost = "
                    f"{solution['cost']:.6f}"
                )
    
                if (
                    iteration_best is None
                    or
                    solution["cost"]
                    < iteration_best["cost"]
                ):
    
                    iteration_best = copy.deepcopy(solution)
    
            # ========================================================
            # accept/reject
            # ========================================================
    
            if iteration_best is None:
    
                print("No feasible candidate.")
    
                stagnation += 1
    
                continue
    
            # ========================================================
            # improvement
            # ========================================================
    
            if (
                iteration_best["cost"]
                <
                incumbent["cost"] - cost_improve_tol
            ):
    
                print("\n*** IMPROVED SOLUTION FOUND ***")
    
                incumbent = copy.deepcopy(iteration_best)
    
                if (
                    incumbent["cost"]
                    <
                    best_global["cost"]
                ):
    
                    best_global = copy.deepcopy(incumbent)
    
                stagnation = 0
    
                # ----------------------------------------------------
                # intensification
                # ----------------------------------------------------
    
                diversification_fraction *= 0.9
    
                modify_fraction *= 0.9
    
            else:
    
                print("\nNo improvement.")
    
                stagnation += 1
    
                # ----------------------------------------------------
                # diversification
                # ----------------------------------------------------
    
                diversification_fraction = min(
                    0.95,
                    diversification_fraction + 0.1
                )
    
                modify_fraction = min(
                    0.75,
                    modify_fraction + 0.05
                )
    
            print(
                f"Current Best Cost = "
                f"{best_global['cost']:.6f}"
            )
    
            # ========================================================
            # strong diversification
            # ========================================================
    
            if stagnation >= stagnation_limit:
    
                print("\nSTRONG DIVERSIFICATION ACTIVATED")
    
                diversification_fraction = 0.95
    
                modify_fraction = min(
                    0.90,
                    modify_fraction + 0.2
                )
    
                stagnation = 0
    
        # ============================================================
        # FINAL
        # ============================================================
    
        print("\n" + "=" * 75)
    
        print(
            f"FINAL BEST COST = "
            f"{best_global['cost']:.6f}"
        )
    
        print("=" * 75)
    
        return (
            best_global["q"],
            best_global["h"],
            best_global["y"],
            best_global["z"]
        )

    """
    Sequential Continuation Segment-Cut Heuristic
    ==============================================
    
    Algorithm overview:
      1. Rank arcs by dual variable magnitude
      2. For each candidate arc, attempt sequential segment reduction: n → n-1 → n-2 → ...
      3. Accept any feasible result; track cost improvement
      4. If improved  → continue reducing the same arc
      5. If feasible but no improvement → restart outer loop from updated structure
      6. If infeasible → skip to next arc
    """
    # ──────────────────────────────────────────────────────────────────────────────
    # Main heuristic
    # ──────────────────────────────────────────────────────────────────────────────
    
    def sequential_segment_cut_heuristic(
        self,
        max_outer_iterations=10,
        cost_improve_tol=1e-4,
        max_segment_reduction_per_arc=10
    ):
        """
        Sequential Continuation Segment-Cut Heuristic
        ─────────────────────────────────────────────
        For each candidate arc (ranked by |dual|), attempts to reduce its active
        segment index by solving a constrained NLP.  Feasible improvements are
        greedily accepted; the outer loop restarts whenever a feasible (even
        non-improving) update changes the problem structure.
        """
        from contextlib import contextmanager
        from amplpy import AMPL
        
        
        # ──────────────────────────────────────────────────────────────────────────────
        # Formatting helpers
        # ──────────────────────────────────────────────────────────────────────────────
        
        def _hline(widths, char="─", cross="┼", left="├", right="┤"):
            return left + cross.join(char * w for w in widths) + right
        
        
        def _header_line(widths, char="─", cross="┬", left="┌", right="┐"):
            return left + cross.join(char * w for w in widths) + right
        
        
        def _footer_line(widths, char="─", cross="┴", left="└", right="┘"):
            return left + cross.join(char * w for w in widths) + right
        
        
        def _row(cells, widths):
            parts = []
            for cell, w in zip(cells, widths):
                text = str(cell)
                parts.append(" " + text.ljust(w - 2) + " ")
            return "│" + "│".join(parts) + "│"
        
        
        def _center(text, width):
            pad = width - len(text)
            lpad = pad // 2
            rpad = pad - lpad
            return " " * lpad + text + " " * rpad
        
        
        def print_banner(title, width=92):
            print("\n" + "═" * width)
            print(_center(title, width))
            print("═" * width)
        
        
        def print_section(title, width=92):
            print("\n" + "─" * width)
            print(f"  {title}")
            print("─" * width)
        
        
        def print_kv_table(pairs, title=None):
            """Print a two-column key-value table."""
            if title:
                print(f"\n  {title}")
            w_key = max(len(k) for k, _ in pairs) + 2
            w_val = max(len(str(v)) for _, v in pairs) + 2
            widths = [w_key, w_val]
            print("  " + _header_line(widths))
            for k, v in pairs:
                print("  " + _row([k, v], widths))
            print("  " + _footer_line(widths))
        
        
        def print_solve_table(rows, title=None):
            """
            Print the per-solve results table.
        
            Columns:
              Iter | Arc | Seg | Target | Result | Old cost | New cost | Δ cost | Time(s) | Cumul(s)
            """
            if title:
                print(f"\n  {title}")
        
            headers = [
                "Outer", "Arc", "Seg",
                "→Tgt", "Status",
                "Old cost", "New cost", "Δ cost",
                "NLP #", "Time(s)", "Cumul(s)"
            ]
            widths = [7, 12, 5, 5, 12, 14, 14, 14, 7, 9, 9]
        
            print("  " + _header_line(widths))
            print("  │" + "│".join(_center(h, w) for h, w in zip(headers, widths)) + "│")
            print("  " + _hline(widths))
        
            for r in rows:
                cells = [
                    str(r.get("outer", "")),
                    str(r.get("arc", "")),
                    str(r.get("seg", "")),
                    str(r.get("target", "")),
                    str(r.get("status", "")),
                    f"{r['old_cost']:.3f}" if r.get("old_cost") is not None else "—",
                    f"{r['new_cost']:.3f}" if r.get("new_cost") is not None else "—",
                    f"{r['delta']:+.3f}" if r.get("delta") is not None else "—",
                    str(r.get("nlp_num", "")),
                    f"{r['solve_time']:.2f}" if r.get("solve_time") is not None else "—",
                    f"{r['cumul_time']:.2f}" if r.get("cumul_time") is not None else "—",
                ]
                print("  " + _row(cells, widths))
        
            print("  " + _footer_line(widths))
        
        
        def print_segment_table(seg_index, changed_arcs=None, title=None):
            """
            Print segment assignments for all arcs.
            Highlights arcs whose segment changed (if changed_arcs provided).
            """
            if title:
                print(f"\n  {title}")
        
            changed_arcs = changed_arcs or set()
            widths = [14, 12, 8]
            print("  " + _header_line(widths))
            print("  │" + "│".join(_center(h, w)
                  for h, w in zip(["Arc", "Segments", "Changed"], widths)) + "│")
            print("  " + _hline(widths))
        
            for arc, segs in sorted(seg_index.items()):
                flag = "◀ YES" if arc in changed_arcs else ""
                cells = [str(arc), str(segs), flag]
                print("  " + _row(cells, widths))
        
            print("  " + _footer_line(widths))
        
        
        def print_candidate_arcs_table(sorted_arcs, dual_dict, seg_index, title=None):
            """Print ranked candidate arcs with their dual values and current segments."""
            if title:
                print(f"\n  {title}")
        
            widths = [5, 12, 12, 10]
            print("  " + _header_line(widths))
            print("  │" + "│".join(_center(h, w)
                  for h, w in zip(["Rank", "Arc", "Dual |val|", "Seg"], widths)) + "│")
            print("  " + _hline(widths))
        
            for rank, arc in enumerate(sorted_arcs, 1):
                dual_val = dual_dict.get(arc, 0.0)
                segs     = seg_index.get(arc, [])
                cells    = [str(rank), str(arc), f"{abs(dual_val):.3f}", str(segs)]
                print("  " + _row(cells, widths))
        
            if not sorted_arcs:
                widths_total = sum(widths) + len(widths) - 1
                print("  │" + _center("(no candidate arcs)", widths_total) + "│")
        
            print("  " + _footer_line(widths))
    
 
        # ── alpha breakpoints ────────────────────────────────────────────────────
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        NP  = len(alpha_vals)
        tol = 1e-8
    
        # ── global timing ────────────────────────────────────────────────────────
        heuristic_start = time.time()
    
        # ── log accumulator for the solve table ─────────────────────────────────
        solve_log = []
    
        # ── NLP counter (may already be initialised on self) ────────────────────
        if not hasattr(self, "number_of_nlp"):
            self.number_of_nlp = 0
    
        # =========================================================================
        # Helper: segment index
        # =========================================================================
        def compute_segment_index(y_sol):
            """Return {arc: [segment_indices]} for every arc."""
            seg_index = {}
    
            for (u, v) in self.arcs:
                y_val           = y_sol[(u, v)]
                active_segments = []
    
                for k in range(len(alpha_vals) - 1):
                    a_high = alpha_vals[k]
                    a_low  = alpha_vals[k + 1]
    
                    if (a_low + tol) < y_val < (a_high - tol):
                        # strictly interior to segment k+1
                        active_segments = [k + 1]
                        break
    
                    elif abs(y_val - a_high) <= tol:
                        # on upper breakpoint
                        active_segments = [1] if k == 0 else [k, k + 1]
                        break
    
                    elif abs(y_val - a_low) <= tol:
                        # on lower breakpoint
                        active_segments = (
                            [k + 1] if k == len(alpha_vals) - 2
                            else [k + 1, k + 2]
                        )
                        break
    
                # safeguard: y outside the alpha range or no match found
                if not active_segments:
                    if y_val > alpha_vals[0]:
                        active_segments = [1]
                    elif y_val < alpha_vals[-1]:
                        active_segments = [NP - 1]
                    else:
                        distances = [abs(y_val - alpha_vals[k]) for k in range(NP)]
                        idx = distances.index(min(distances))
                        if idx == 0:
                            active_segments = [1]
                        elif idx == NP - 1:
                            active_segments = [NP - 1]
                        else:
                            active_segments = [idx, idx + 1]
    
                seg_index[(u, v)] = active_segments
    
            return seg_index
    
        # =========================================================================
        # Helper: true piecewise-linear cost
        # =========================================================================
        def compute_true_cost(y_sol, seg_sol):
            return sum(
                self.L[(i, j)] * (
                    self.slope[seg_sol[(i, j)][0]] * y_sol[(i, j)]
                    + self.intercept[seg_sol[(i, j)][0]]
                )
                for (i, j) in self.arcs
            )

        # =========================================================================
        # Helper: warm-start and solve a single NLP
        # =========================================================================
        def build_and_solve_nlp(i, j, target_seg):
            """
            Construct the AMPL model with current warm-start values and a
            segment-cut constraint on arc (i,j), then solve with IPOPT.

            Returns (ampl_instance, solve_result, solve_time)
            """
            ampl = AMPL()
            ampl.reset()

            # load model
            model_map = {5: "newyork_model.mod", 6: "blacksburg_model.mod"}
            ampl.read(model_map.get(self.data_number, "exact_reduced_wdn.mod"))
            ampl.read_data(self.data_file)

            # warm-start primal variables
            for (u, v), val in self.y.items():
                ampl.eval(f"let y[{u},{v}] := {val};")
    
            for (u, v), val in self.q.items():
                ampl.eval(f"let q[{u},{v}] := {val};")
                if self.data_number == 5:
                    ampl.eval(f"let q1[{u},{v}] := {self.q1[u, v]};")
                    ampl.eval(f"let q2[{u},{v}] := {self.q2[u, v]};")
    
            for u, val in self.h.items():
                ampl.eval(f"let h[{u}] := {val};")
    
            for (u, v), val in self.z.items():
                ampl.eval(f"let z[{u},{v}] := {val};")

            # warm-start duals
            # current_duals = {
            #     name: con.get_values()
            #     for name, con in ampl.get_constraints()
            # }
            # for con_name, dual_values in self.all_duals.items():
            #     if con_name in current_duals:
            #         ampl.get_constraint(con_name).set_values(dual_values)

            # segment-cut constraint
            ampl.eval(f"""
                s.t. seg_cut_lower:
                      # {self.y[i,j]} <= y[{i},{j}];
                    alpha[{target_seg}] <= y[{i},{j}];
            """)

            # solver options
            ampl.option["solver"] = "ipopt"
            ampl.set_option(
                "ipopt_options",
                "outlev=0 bound_relax_factor=0 warm_start_init_point=yes"
            )
            ampl.option["presolve_eps"] = "6.82e-14"
            ampl.option["presolve"]     = 1

            with self.suppress_output():
                ampl.solve()
            self.all_duals = {}
            for con_name, con in ampl.get_constraints():
                self.all_duals[con_name] = con.getValues()

            solve_time = ampl.get_value("_solve_elapsed_time")
            return ampl, ampl.solve_result, solve_time

        # =========================================================================
        # Initialisation
        # =========================================================================
        self.seg_index   = compute_segment_index(self.y)
        self.current_cost = compute_true_cost(self.y, self.seg_index)

        best_solution = {
            "q":    copy.deepcopy(self.q),
            "h":    copy.deepcopy(self.h),
            "y":    copy.deepcopy(self.y),
            "z":    copy.deepcopy(self.z),
            "seg":  copy.deepcopy(self.seg_index),
            "cost": self.current_cost,
        }

        print_banner("SEQUENTIAL CONTINUATION SEGMENT-CUT HEURISTIC")
        print_kv_table([
            ("Initial cost",      f"{self.current_cost:.3f}"),
            ("Max outer iters",   str(max_outer_iterations)),
            ("Improve tolerance", str(cost_improve_tol)),
            ("Max seg reductions", str(max_segment_reduction_per_arc)),
            # ("Alpha breakpoints", (alpha_vals)),
            ("Total arcs",        str(len(self.arcs))),
        ], title="Configuration")

        print_segment_table(self.seg_index, title="Initial segment assignments")

        outer_it             = 0
        restart_global_search = True

        # =========================================================================
        # OUTER LOOP
        # =========================================================================
        while restart_global_search and outer_it < max_outer_iterations:
            outer_it             += 1
            restart_global_search = False

            print_section(f"OUTER ITERATION {outer_it}")

            # ── candidate arcs (max segment ≥ 2) ─────────────────────────────────
            candidate_arcs = [
                (i, j) for (i, j) in self.sorted_arcs
                if max(self.seg_index[i, j]) >= 2 
                if (i,j) not in self.visited_arc_reverse
            ]

            # ── extract duals for candidate arcs ─────────────────────────────────
            dual_dict = {}
            
            for con_name, dual_values in self.all_duals.items():
                if con_name == "con2":
                    tmp       = dual_values.to_dict()
                    dual_dict = {arc: val for arc, val in tmp.items()
                                 if arc in candidate_arcs}
                    break

            sorted_arcs = sorted(
                dual_dict, key=lambda a: abs(dual_dict[a]), reverse=False
            )

            print_candidate_arcs_table(
                sorted_arcs, dual_dict, self.seg_index,
                title=f"Candidate arcs (outer iter {outer_it})"
            )

            # ── iterate over arcs ─────────────────────────────────────────────────
            for (i, j) in sorted_arcs:

                print_section(f"Arc ({i},{j})  —  initial seg = {self.seg_index[(i,j)]}")

                self.visited_arc_reverse.append((i, j))
                current_seg       = self.seg_index[(i, j)][0]
                reduction_counter = 0

                if current_seg <= 1:
                    print("  ⊘  Already at minimum segment — skipping.")
                    continue

                # ── CONTINUATION LOOP ─────────────────────────────────────────────
                # while (current_seg > 1 ):
                # while (current_seg > 1 and reduction_counter < max_segment_reduction_per_arc):

                target_seg        = current_seg - 1
                reduction_counter += 1
                iter_start         = time.time()

                # ── solve NLP ─────────────────────────────────────────────────
                ampl_inst, solve_result, solve_time = build_and_solve_nlp(i, j, target_seg)
                self.number_of_nlp += 1
                cumul_time          = time.time() - heuristic_start

                # ── parse solution ────────────────────────────────────────────
                if solve_result == "solved":
                    y_new    = ampl_inst.getVariable("y").getValues().to_dict()
                    q_new    = ampl_inst.getVariable("q").getValues().to_dict()
                    h_new    = ampl_inst.getVariable("h").getValues().to_dict()
                    z_new    = ampl_inst.getVariable("z").getValues().to_dict()
                    seg_new  = compute_segment_index(y_new)
                    new_cost = compute_true_cost(y_new, seg_new)
                    achieved = seg_new[(i, j)][0]
                    delta    = new_cost - self.current_cost
                    status   = "IMPROVED" if delta < -cost_improve_tol else "feasible"
                else:
                    new_cost = None
                    achieved = None
                    delta    = None
                    status   = "infeasible"

                # ── log row ───────────────────────────────────────────────────
                solve_log.append({
                    "outer":       outer_it,
                    "arc":         f"({i},{j})",
                    "seg":         current_seg,
                    "target":      target_seg,
                    "status":      status,
                    "old_cost":    self.current_cost,
                    "new_cost":    new_cost,
                    "delta":       delta,
                    "nlp_num":     self.number_of_nlp,
                    "solve_time":  solve_time,
                    "cumul_time":  cumul_time,
                })

                # ── display running solve table ───────────────────────────────
                print_solve_table([solve_log[-1]], title="  Latest solve")

                # ── infeasible: stop this arc ─────────────────────────────────
                if solve_result != "solved":
                    print("  ✗  NLP infeasible — moving to next arc.")
                    break

                # ── always accept feasible structure ─────────────────────────
                prev_seg_index   = copy.deepcopy(self.seg_index)
                self.ampl        = ampl_inst
                self.y           = copy.deepcopy(y_new)
                self.q           = copy.deepcopy(q_new)
                self.h           = copy.deepcopy(h_new)
                self.z           = copy.deepcopy(z_new)
                self.seg_index   = copy.deepcopy(seg_new)
                self.network_graph = self.generate_random_acyclic_from_solution(q_new)

                # show what changed in segments
                changed_arcs = {
                    arc for arc in self.arcs
                    if self.seg_index[arc] != prev_seg_index[arc]
                }
                # print_segment_table(
                #     self.seg_index,
                #     changed_arcs=changed_arcs,
                #     title=f"  Updated segments after seg {current_seg}→{target_seg}"
                # )

                # ── improved ─────────────────────────────────────────────────
                if delta < -cost_improve_tol:
                    print(f"  ★  IMPROVEMENT  {self.current_cost:.6f} → {new_cost:.6f}"
                          f"  (Δ = {delta:+.6f})")
                    self.current_cost = new_cost
                    best_solution = {
                        "q":    copy.deepcopy(q_new),
                        "h":    copy.deepcopy(h_new),
                        "y":    copy.deepcopy(y_new),
                        "z":    copy.deepcopy(z_new),
                        "seg":  copy.deepcopy(seg_new),
                        "cost": new_cost,
                    }
                    current_seg = achieved    # continue on same arc
                    restart_global_search = True
                    # continue
                    # break

                # ── feasible but not improving ────────────────────────────────
                # print("  ↺  Feasible but no improvement — restarting outer loop.")
                # break

                # ── restart outer loop if flagged ─────────────────────────────────
                if restart_global_search:
                    break

        # =========================================================================
        # Summary
        # =========================================================================
        total_time = time.time() - heuristic_start

        print_banner("FINAL RESULTS")
        print_kv_table([
            ("Best cost",         f"{best_solution['cost']:.3f}"),
            ("Initial cost",      f"{solve_log[0]['old_cost']:.3f}" if solve_log else "—"),
            ("Total improvement", f"{best_solution['cost'] - (solve_log[0]['old_cost'] if solve_log else best_solution['cost']):+.6f}"),
            ("Outer iterations",  str(outer_it)),
            ("NLP solves",        str(self.number_of_nlp)),
            ("Total wall time",   f"{total_time:.2f}s"),
        ], title="Summary statistics")

        if solve_log:
            print_solve_table(solve_log, title="Complete solve history")

        # print_segment_table(best_solution["seg"], title="Best solution — segment assignments")

        return (
            best_solution["q"],
            best_solution["h"],
            best_solution["y"],
            best_solution["z"],
        )


    def lagrangian_based_heuristic(
        self,
        q_star,
        h_star,
        y_star,
        z_star,
        seg_index,
        theta=0.05,
        max_outer=20,
        top_k_fraction=0.25,
        distance_weight=5.0,
        cost_improve_tol=1e-4,
        segment_change_tol=2
    ):
    
        """
        Active-set dual perturbation heuristic
        with multi-active breakpoint segments.
        """
    
        from amplpy import AMPL
    
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
    
        NP = len(alpha_vals)
    
        # ============================================================
        # utilities
        # ============================================================
    
        def find_seg(yv):
    
            tol = 1e-8
    
            active_segments = []
    
            for k in range(NP - 1):
    
                a_high = alpha_vals[k]
                a_low = alpha_vals[k + 1]
    
                # ----------------------------------------------------
                # interior
                # ----------------------------------------------------
    
                if (a_low + tol) < yv < (a_high - tol):
    
                    active_segments = [k + 1]
                    break
    
                # ----------------------------------------------------
                # upper breakpoint
                # ----------------------------------------------------
    
                elif abs(yv - a_high) <= tol:
    
                    if k == 0:
                        active_segments = [1]
                    else:
                        active_segments = [k, k + 1]
    
                    break
    
                # ----------------------------------------------------
                # lower breakpoint
                # ----------------------------------------------------
    
                elif abs(yv - a_low) <= tol:
    
                    if k == NP - 2:
                        active_segments = [k + 1]
                    else:
                        active_segments = [k + 1, k + 2]
    
                    break
    
            # --------------------------------------------------------
            # safeguard
            # --------------------------------------------------------
    
            if not active_segments:
    
                if yv > alpha_vals[0]:
    
                    active_segments = [1]
    
                elif yv < alpha_vals[-1]:
    
                    active_segments = [NP - 1]
    
                else:
    
                    distances = [
                        abs(yv - alpha_vals[k])
                        for k in range(NP)
                    ]
    
                    idx = distances.index(min(distances))
    
                    if idx == 0:
                        active_segments = [1]
    
                    elif idx == NP - 1:
                        active_segments = [NP - 1]
    
                    else:
                        active_segments = [idx, idx + 1]
    
            return active_segments
    
        # ============================================================
        # exact cost
        # ============================================================
    
        def true_cost(y_sol, seg_sol):
    
            total = 0.0
    
            for (i, j) in self.arcs:
    
                s = seg_sol[(i, j)][0]
    
                total += (
                    self.L[(i, j)]
                    * (
                        self.slope[s] * y_sol[(i, j)]
                        + self.intercept[s]
                    )
                )
    
            return total
    
        # ============================================================
        # expensive arcs
        # ============================================================
    
        def get_expensive_arcs(z_sol):
    
            arc_cost = sorted(
                [
                    ((i, j), z_sol[(i, j)])
                    for (i, j) in self.arcs
                ],
                key=lambda x: -x[1]
            )
    
            k = max(1, int(top_k_fraction * len(self.arcs)))
    
            return set(a for a, _ in arc_cost[:k])
    
        # ============================================================
        # initial incumbent
        # ============================================================
    
        C_star = true_cost(y_star, seg_index)
    
        incumbent = {
            "q": q_star.copy(),
            "h": h_star.copy(),
            "y": y_star.copy(),
            "z": z_star.copy(),
            "seg": seg_index.copy(),
            "cost": C_star
        }
    
        best_global = incumbent.copy()
    
        # ============================================================
        # initialize vartheta
        # ============================================================
    
        vartheta = {}
    
        for (i, j) in self.arcs:
    
            active_set = seg_index[(i, j)]
    
            mass = 1.0 / len(active_set)
    
            for s in range(1, NP):
    
                if s in active_set:
                    vartheta[(i, j, s)] = mass
                else:
                    vartheta[(i, j, s)] = 0.0
    
        print("\n" + "=" * 75)
        print(" ACTIVE-SET LAGRANGIAN ESCAPE HEURISTIC")
        print(f" INITIAL COST = {C_star:.6f}")
        print("=" * 75)
    
        stagnation = 0
    
        # ============================================================
        # OUTER ITERATIONS
        # ============================================================
    
        for outer_it in range(1, max_outer + 1):
    
            print(f"\nOUTER ITERATION {outer_it}")
            print("-" * 75)
    
            y_inc = incumbent["y"]
            seg_inc = incumbent["seg"]
    
            # ========================================================
            # STEP 1:
            # expensive arcs
            # ========================================================
    
            expensive_arcs = get_expensive_arcs(incumbent["z"])
    
            print("Theta value:", theta)
            print(f"Expensive arcs selected: {len(expensive_arcs)}")
    
            # ========================================================
            # STEP 2:
            # update dual multipliers
            # ========================================================
    
            vartheta_new = {}
    
            for (i, j) in self.arcs:
    
                active_set = seg_inc[(i, j)]
    
                # ----------------------------------------------------
                # initialize
                # ----------------------------------------------------
    
                for s in range(1, NP):
    
                    vartheta_new[(i, j, s)] = 0.0
    
                # ----------------------------------------------------
                # current mass
                # ----------------------------------------------------
    
                mass_cur = 1.0 / len(active_set)
    
                # ----------------------------------------------------
                # neighboring cheaper segments
                # ----------------------------------------------------
    
                target_set = set()
    
                for s_cur in active_set:
    
                    if s_cur > 1:
                        target_set.add(s_cur - 1)
                    else:
                        target_set.add(1)
    
                target_set = list(target_set)
    
                mass_target = theta / len(target_set)
    
                # ----------------------------------------------------
                # retain current mass
                # ----------------------------------------------------
    
                for s_cur in active_set:
    
                    vartheta_new[(i, j, s_cur)] += (
                        (1.0 - theta) * mass_cur
                    )
    
                # ----------------------------------------------------
                # transfer perturbation mass
                # ----------------------------------------------------
    
                for s_tar in target_set:
    
                    vartheta_new[(i, j, s_tar)] += mass_target
    
            vartheta = vartheta_new
    
            # ========================================================
            # STEP 3:
            # construct coefficients
            # ========================================================
    
            coef_z = {}
            coef_y = {}
            const_term = {}
    
            for (i, j) in self.arcs:
    
                coef_z[(i, j)] = 1.0
                coef_y[(i, j)] = 0.0
                const_term[(i, j)] = 0.0
    
                for s in range(1, NP):
    
                    v = vartheta[(i, j, s)]
    
                    coef_z[(i, j)] -= v
    
                    coef_y[(i, j)] += (
                        v
                        * self.L[(i, j)]
                        * self.slope[s]
                    )
    
                    const_term[(i, j)] += (
                        v
                        * self.L[(i, j)]
                        * self.intercept[s]
                    )
    
            # ========================================================
            # STEP 4:
            # solve relaxed NLP
            # ========================================================
    
            ampl = AMPL()
    
            ampl.read("exact_reduced_wdn.mod")
            ampl.read_data(self.data_file)
    
            ampl.option["solver"] = "ipopt"
    
            ampl.option["ipopt_options"] = (
                "outlev=0 "
                "warm_start_init_point=yes "
                "bound_relax_factor=0 "
            )
    
            ampl.eval("""
    
                param coef_z{arcs};
                param coef_y{arcs};
                param const_term{arcs};
    
                param y_ref{arcs};
                param distance_weight default 1;
    
                drop total_cost;
    
                minimize relaxed_obj:
    
                    sum{(i,j) in arcs}
                    (
                        coef_z[i,j] * z[i,j]
                        + coef_y[i,j] * y[i,j]
                        + const_term[i,j]
                    );
            """)
    
            ampl.param["distance_weight"] = distance_weight
    
            for (i, j) in self.arcs:
    
                ampl.param["coef_z"][i, j] = coef_z[(i, j)]
                ampl.param["coef_y"][i, j] = coef_y[(i, j)]
                ampl.param["const_term"][i, j] = const_term[(i, j)]
    
                ampl.param["y_ref"][i, j] = y_inc[(i, j)]
    
            # ========================================================
            # warm start
            # ========================================================
    
            for (i, j) in self.arcs:
    
                ampl.var["q"][i, j] = incumbent["q"][(i, j)]
                ampl.var["y"][i, j] = incumbent["y"][(i, j)]
                ampl.var["z"][i, j] = incumbent["z"][(i, j)]
    
            for i in self.nodes:
    
                ampl.var["h"][i] = incumbent["h"][i]
    
            # ========================================================
            # relax epigraph constraints
            # ========================================================
    
            ampl.eval("drop exact_cost;")
    
            for (i, j) in self.arcs:
    
                active_set = seg_inc[(i, j)]
    
                relaxed_set = set(active_set)
    
                # ----------------------------------------------------
                # also relax neighboring cheaper segments
                # ----------------------------------------------------
    
                for s_cur in active_set:
    
                    if s_cur > 1:
                        relaxed_set.add(s_cur - 1)
    
                # ----------------------------------------------------
                # keep remaining constraints
                # ----------------------------------------------------
    
                for s in self.segs:
    
                    if s not in relaxed_set:
    
                        ampl.eval(f"""
                            s.t. epigraph_{i}_{j}_{s}:
    
                                z[{i},{j}]
                                >=
                                L[{i},{j}] *
                                (
                                    slope[{s}] * y[{i},{j}]
                                    + intercept[{s}]
                                );
                        """)
    
            # ========================================================
            # solve relaxed NLP
            # ========================================================
    
            with self.suppress_output():
                ampl.solve()
    
            status = ampl.get_value("solve_result")
    
            if status != "solved":
    
                print("Relaxed NLP failed.")
    
                theta = min(0.95, 1.5 * theta)
    
                continue
    
            # ========================================================
            # retrieve relaxed solution
            # ========================================================
    
            q_rel = ampl.get_variable("q").get_values().to_dict()
    
            y_rel = ampl.get_variable("y").get_values().to_dict()
    
            h_rel = ampl.get_variable("h").get_values().to_dict()
    
            z_rel = ampl.get_variable("z").get_values().to_dict()
    
            seg_rel = {
                (i, j): find_seg(y_rel[(i, j)])
                for (i, j) in self.arcs
            }
    
            C_relaxed_model = true_cost(y_rel, seg_rel)
    
            # ========================================================
            # segment changes
            # ========================================================
    
            changed_arcs = [
                (i, j)
                for (i, j) in self.arcs
                if set(seg_rel[(i, j)]) != set(seg_inc[(i, j)])
            ]
    
            seg_changes = len(changed_arcs)
    
            print(f"Segment changes = {seg_changes}")
    
            print("Relaxed model objective value:", C_relaxed_model)
    
            # ========================================================
            # insufficient escape
            # ========================================================
    
            if seg_changes < segment_change_tol:
    
                theta = min(0.95, theta + 0.1)
    
                distance_weight *= 1.2
    
                stagnation += 1
    
                print("Insufficient basin escape.")
    
                continue
    
            # ========================================================
            # STEP 5:
            # exact NLP projection
            # ========================================================
    
            ampl_ex = AMPL()
    
            ampl_ex.read("exact_reduced_wdn.mod")
            ampl_ex.read_data(self.data_file)
    
            ampl_ex.option["solver"] = "ipopt"
    
            ampl_ex.option["ipopt_options"] = (
                "outlev=0 "
                "warm_start_init_point=yes "
                "bound_relax_factor=0 "
            )
    
            # --------------------------------------------------------
            # warm start
            # --------------------------------------------------------
    
            for (i, j) in self.arcs:
    
                ampl_ex.var["q"][i, j] = q_rel[(i, j)]
                ampl_ex.var["y"][i, j] = y_rel[(i, j)]
                ampl_ex.var["z"][i, j] = z_rel[(i, j)]
    
            for i in self.nodes:
    
                ampl_ex.var["h"][i] = h_rel[i]
    
            # ========================================================
            # exact solve
            # ========================================================
    
            with self.suppress_output():
                ampl_ex.solve()
    
            status_ex = ampl_ex.get_value("solve_result")
    
            if status_ex != "solved":
    
                print("Exact projection failed.")
    
                theta *= 0.7
    
                continue
    
            # ========================================================
            # retrieve exact solution
            # ========================================================
    
            q_ex = ampl_ex.get_variable("q").get_values().to_dict()
    
            y_ex = ampl_ex.get_variable("y").get_values().to_dict()
    
            h_ex = ampl_ex.get_variable("h").get_values().to_dict()
    
            z_ex = ampl_ex.get_variable("z").get_values().to_dict()
    
            seg_ex = {
                (i, j): find_seg(y_ex[(i, j)])
                for (i, j) in self.arcs
            }
    
            cost_ex = true_cost(y_ex, seg_ex)
    
            print(
                f"Projected exact cost = {cost_ex:.6f} "
                f"Best cost = {best_global['cost']:.6f}"
            )
    
            # ========================================================
            # STEP 6:
            # accept/reject
            # ========================================================
    
            if cost_ex < incumbent["cost"] - cost_improve_tol:
    
                print("*** IMPROVED LOCAL SOLUTION FOUND ***")
    
                incumbent = {
                    "q": q_ex.copy(),
                    "h": h_ex.copy(),
                    "y": y_ex.copy(),
                    "z": z_ex.copy(),
                    "seg": seg_ex.copy(),
                    "cost": cost_ex
                }
    
                if cost_ex < best_global["cost"]:
    
                    best_global = incumbent.copy()
    
                theta = max(0.15, 0.5 * theta)
    
                distance_weight *= 0.8
    
                stagnation = 0
    
            else:
    
                print("No improvement.")
    
                theta = min(0.95, theta + 0.05)
    
                distance_weight *= 1.1
    
                stagnation += 1
    
            # ========================================================
            # diversification
            # ========================================================
    
            if stagnation >= 5:
    
                print("Strong diversification activated.")
    
                theta = min(0.95, theta + 0.25)
    
                distance_weight *= 1.5
    
                stagnation = 0
    
        # ============================================================
        # FINAL
        # ============================================================
    
        print("\n" + "=" * 75)
    
        print(f" FINAL BEST COST = {best_global['cost']:.6f}")
    
        print("=" * 75)
    
        return (
            best_global["q"],
            best_global["h"],
            best_global["y"],
            best_global["z"]
        )


    def lagrangian_based_heuristic3(
        self,
        q_star,
        h_star,
        y_star,
        z_star,
        seg_index,
        theta=0.0001,
        max_outer=30,
        top_k_fraction=0.25,
        distance_weight=5.0,
        cost_improve_tol=1e-4,
        segment_change_tol=2
    ):
        """
        Active-set dual perturbation heuristic.

        Main ideas:
        ----------
        1. Relax currently active epigraph planes.
        2. Penalize staying near previous basin.
        3. Push only expensive arcs.
        4. Encourage neighboring cheaper segments.
        5. Adaptive theta update.
        6. Exact NLP projection.
        """

        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        NP = len(alpha_vals)

        # ============================================================
        # utilities
        # ============================================================

        def find_seg(yv):
            for k in range(NP - 1):
                if alpha_vals[k + 1] - 1e-10 <= yv <= alpha_vals[k] + 1e-10:
                    return k + 1
            return 1 if yv > alpha_vals[0] else NP - 1

        def true_cost(y_sol, seg_sol):
            return sum(
                self.L[i, j]
                * (
                    self.slope[seg_sol[i, j]] * y_sol[i, j]
                    + self.intercept[seg_sol[i, j]]
                )
                for (i, j) in self.arcs
            )

        def get_expensive_arcs(z_sol):
            arc_cost = sorted(
                [
                    ((i, j), z_sol[i, j])
                    for (i, j) in self.arcs
                ],
                key=lambda x: -x[1]
            )

            k = max(1, int(top_k_fraction * len(self.arcs)))

            return set(a for a, _ in arc_cost[:k])

        # ============================================================
        # initial incumbent
        # ============================================================

        C_star = true_cost(y_star, seg_index)

        incumbent = {
            "q": q_star.copy(),
            "h": h_star.copy(),
            "y": y_star.copy(),
            "z": z_star.copy(),
            "seg": seg_index.copy(),
            "cost": C_star
        }

        best_global = incumbent.copy()

        # ============================================================
        # initialize vartheta multipliers
        # ============================================================

        vartheta = {}

        for (i, j) in self.arcs:
            for s in range(1, NP):
                vartheta[(i, j, s)] = (
                    1.0 if s == seg_index[(i, j)] else 0.0
                )

        print("\n" + "=" * 75)
        print(f" ACTIVE-SET LAGRANGIAN ESCAPE HEURISTIC")
        print(f" INITIAL COST = {C_star:.6f}")
        print("=" * 75)

        stagnation = 0

        # ============================================================
        # outer iterations
        # ============================================================

        for outer_it in range(1, max_outer + 1):

            print(f"\nOUTER ITERATION {outer_it}")
            print("-" * 75)

            y_inc = incumbent["y"]
            seg_inc = incumbent["seg"]

            # ========================================================
            # STEP 1:
            # identify expensive arcs
            # ========================================================

            expensive_arcs = get_expensive_arcs(incumbent["z"])
            print("Theta value:", theta)
            print(f"Expensive arcs selected: {len(expensive_arcs)}")

            # ========================================================
            # STEP 2:
            # update dual multipliers
            # ========================================================

            vartheta_new = {}

            for (i, j) in self.arcs:

                s_cur = seg_inc[(i, j)]

                # initialize all multipliers to zero
                for s in range(1, NP):
                    vartheta_new[(i, j, s)] = 0.0

                if s_cur==1:
                    vartheta_new[(i,j,s_cur)] = 1
                    continue
                # ----------------------------------------------------
                # only perturb expensive arcs
                # ----------------------------------------------------

                # if (i, j) not in expensive_arcs:
                #
                #     vartheta_new[(i, j, s_cur)] = 1.0
                #     continue

                # ----------------------------------------------------
                # choose neighboring cheaper segment
                # ----------------------------------------------------

                s_target = min(NP - 1, s_cur - 1)

                # ----------------------------------------------------
                # transfer mass from current active segment
                # to neighboring cheaper segment
                # ----------------------------------------------------

                vartheta_new[(i, j, s_cur)] = 1.0 - theta
                vartheta_new[(i, j, s_target)] = theta

            vartheta = vartheta_new
            # print("vartheta:",vartheta)
            # ========================================================
            # STEP 3:
            # construct modified objective coefficients
            # ========================================================

            coef_z = {}
            coef_y = {}
            const_term = {}

            for (i, j) in self.arcs:

                coef_z[(i, j)] = 1.0

                coef_y[(i, j)] = 0.0

                const_term[(i, j)] = 0.0

                for s in range(1, NP):

                    v = vartheta[(i, j, s)]

                    coef_z[(i, j)] -= v

                    coef_y[(i, j)] += (
                        v
                        * self.L[(i, j)]
                        * self.slope[s]
                    )

                    const_term[(i, j)] += (
                        v
                        * self.L[(i, j)]
                        * self.intercept[s]
                    )

            # ========================================================
            # STEP 4:
            # solve relaxed NLP
            # ========================================================

            ampl = AMPL()

            ampl.read("exact_reduced_wdn.mod")
            ampl.read_data(self.data_file)

            ampl.option["solver"] = "ipopt"

            ampl.option["ipopt_options"] = (
                "outlev=0 "
                "max_iter=3000 "
                "warm_start_init_point=yes "
                "bound_relax_factor=0 "
            )

            ampl.eval("""
                param coef_z{arcs};
                param coef_y{arcs};
                param const_term{arcs};

                param y_ref{arcs};
                param distance_weight default 1;

                drop total_cost;

                minimize relaxed_obj:

                    sum{(i,j) in arcs}
                    (
                        coef_z[i,j] * z[i,j]
                        + coef_y[i,j] * y[i,j]
                        + const_term[i,j]
                    )

                    # - distance_weight *
                    # sum{(i,j) in arcs}
                    # (y[i,j] - y_ref[i,j])^2
                    ;
            """)

            ampl.param["distance_weight"] = distance_weight

            for (i, j) in self.arcs:

                ampl.param["coef_z"][i, j] = coef_z[(i, j)]
                ampl.param["coef_y"][i, j] = coef_y[(i, j)]
                ampl.param["const_term"][i, j] = const_term[(i, j)]

                ampl.param["y_ref"][i, j] = y_inc[(i, j)]

            # ========================================================
            # warm start
            # ========================================================

            for (i, j) in self.arcs:

                ampl.var["q"][i, j] = incumbent["q"][(i, j)]
                ampl.var["y"][i, j] = incumbent["y"][(i, j)]
                ampl.var["z"][i, j] = incumbent["z"][(i, j)]

            for i in self.nodes:
                ampl.var["h"][i] = incumbent["h"][i]

            # ========================================================
            # relax currently active epigraph constraints
            # ========================================================

            ampl.eval("drop exact_cost;")

            # for (i, j) in self.arcs:
            #
            #     s_cur = seg_inc[(i, j)]
            #
            #     for s in self.segs:
            #
            #         if s != s_cur:
            #
            #             ampl.eval(f"""
            #                 s.t. epigraph_{i}_{j}_{s}:
            #                     z[{i},{j}]
            #                     >=
            #                     L[{i},{j}] *
            #                     (
            #                         slope[{s}] * y[{i},{j}]
            #                         + intercept[{s}]
            #                     );
            #             """)

            for (i, j) in self.arcs:
            
                s_cur = seg_inc[(i, j)]
            
                if s_cur == 1:
                    relaxed_set = {s_cur}
                else:
                    s_target = s_cur - 1
                    relaxed_set = {s_cur, s_target}
            
                for s in self.segs:
            
                    # KEEP only non-relaxed constraints
                    if s not in relaxed_set:
            
                        ampl.eval(f"""
                            s.t. epigraph_{i}_{j}_{s}:
                                z[{i},{j}]
                                >=
                                L[{i},{j}] *
                                (
                                    slope[{s}] * y[{i},{j}]
                                    + intercept[{s}]
                                );
                        """)
            # ========================================================
            # solve relaxed NLP
            # ========================================================

            with self.suppress_output():
                ampl.solve()

            status = ampl.get_value("solve_result")

            if status != "solved":

                print("Relaxed NLP failed.")

                theta = min(0.95, 1.5 * theta)

                continue

            # ========================================================
            # retrieve solution
            # ========================================================

            q_rel = ampl.get_variable("q").get_values().to_dict()
            y_rel = ampl.get_variable("y").get_values().to_dict()
            h_rel = ampl.get_variable("h").get_values().to_dict()
            z_rel = ampl.get_variable("z").get_values().to_dict()

            seg_rel = {
                (i, j): find_seg(y_rel[(i, j)])
                for (i, j) in self.arcs
            }

            C_relaxed_model = true_cost(y_rel, seg_rel)
            # ampl.eval("display y;")
            # ampl.eval("display exact_cost.dual;")
            # for (i, j) in self.arcs:
            #     s_cur = seg_inc[(i, j)]
            #     for s in self.segs:
            #         if s != s_cur:
            #             ampl.eval(f"""display epigraph_{i}_{j}_{s}.dual;""")
            # ========================================================
            # measure segment changes
            # ========================================================

            changed_arcs = [
                (i, j)
                for (i, j) in self.arcs
                if seg_rel[(i, j)] != seg_inc[(i, j)]
            ]

            seg_changes = len(changed_arcs)

            print(f"Segment changes = {seg_changes}")
            print(f"Relaxed model objective value:",C_relaxed_model)
            # ========================================================
            # if insufficient changes => stronger push
            # ========================================================

            if seg_changes < segment_change_tol:

                theta = min(0.95, theta + 0.1)

                distance_weight *= 1.2

                stagnation += 1

                print("Insufficient basin escape.")

                continue

            # ========================================================
            # STEP 5:
            # exact NLP projection
            # ========================================================

            ampl_ex = AMPL()

            ampl_ex.read("exact_reduced_wdn.mod")
            ampl_ex.read_data(self.data_file)

            ampl_ex.option["solver"] = "ipopt"

            ampl_ex.option["ipopt_options"] = (
                "outlev=0 "
                "warm_start_init_point=yes "
                "bound_relax_factor=0 "
            )

            # --------------------------------------------------------
            # warm start from relaxed solution
            # --------------------------------------------------------

            for (i, j) in self.arcs:

                ampl_ex.var["q"][i, j] = q_rel[(i, j)]
                ampl_ex.var["y"][i, j] = y_rel[(i, j)]
                ampl_ex.var["z"][i, j] = z_rel[(i, j)]

            for i in self.nodes:
                ampl_ex.var["h"][i] = h_rel[i]

            # ========================================================
            # exact solve
            # ========================================================

            with self.suppress_output():
                ampl_ex.solve()

            status_ex = ampl_ex.get_value("solve_result")

            if status_ex != "solved":

                print("Exact projection failed.")

                theta *= 0.7

                continue

            # ========================================================
            # retrieve exact solution
            # ========================================================

            q_ex = ampl_ex.get_variable("q").get_values().to_dict()

            y_ex = ampl_ex.get_variable("y").get_values().to_dict()

            h_ex = ampl_ex.get_variable("h").get_values().to_dict()

            z_ex = ampl_ex.get_variable("z").get_values().to_dict()

            seg_ex = {
                (i, j): find_seg(y_ex[(i, j)])
                for (i, j) in self.arcs
            }

            cost_ex = true_cost(y_ex, seg_ex)

            print(f"Projected exact cost = {cost_ex:.6f} Best cost = {best_global["cost"]}")

            # ========================================================
            # STEP 6:
            # accept/reject
            # ========================================================

            if cost_ex < incumbent["cost"] - cost_improve_tol:

                print(f"*** IMPROVED LOCAL SOLUTION FOUND ***")

                incumbent = {
                    "q": q_ex.copy(),
                    "h": h_ex.copy(),
                    "y": y_ex.copy(),
                    "z": z_ex.copy(),
                    "seg": seg_ex.copy(),
                    "cost": cost_ex
                }

                if cost_ex < best_global["cost"]:

                    best_global = incumbent.copy()

                # ----------------------------------------------------
                # reset perturbation after improvement
                # ----------------------------------------------------

                theta = max(0.15, 0.5 * theta)

                distance_weight *= 0.8

                stagnation = 0

            else:

                print("No improvement.")

                theta = min(0.95, theta + 0.05)

                distance_weight *= 1.1

                stagnation += 1

            # ========================================================
            # diversification
            # ========================================================

            if stagnation >= 5:

                print("Strong diversification activated.")

                theta = min(0.95, theta + 0.25)

                distance_weight *= 1.5

                stagnation = 0

        # ============================================================
        # final
        # ============================================================

        print("\n" + "=" * 75)
        print(f" FINAL BEST COST = {best_global['cost']:.6f}")
        print("=" * 75)

        return (
            best_global["q"],
            best_global["h"],
            best_global["y"],
            best_global["z"]
        )

    def lagrangian_based_heuristic1(self,q_star, h_star, y_star, z_star, seg_index, strategy='shift', theta=0.4, max_outer=30):
        """
        Update dual of active epigraph constraints of previous local solution
        and re-solve relaxed model to find different local solution.
        """
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        NP         = len(alpha_vals)

        def find_seg(yv):
            for k in range(NP-1):
                if alpha_vals[k+1]-1e-10 <= yv <= alpha_vals[k]+1e-10:
                    return k+1
            return 1 if yv > alpha_vals[0] else NP-1

        def true_cost(y_sol, seg_index):
            return sum(
                self.L[i,j]*(self.slope[seg_index[i,j]]*y_sol[i,j] + self.intercept[seg_index[i,j]])
                for (i,j) in self.arcs
            )
        # print("seg_index",seg_index) 
        C_star    = true_cost(y_star, seg_index)
        # C_star = sum(self.L[(i,j)] * (self.slope[seg_index[i,j]] * y_sol[(i,j)]+ self.intercept[seg_index[(i,j)]]) for (i,j) in self.arcs)

        incumbent = {
            "q": q_star.copy(), "h": h_star.copy(),
            "y": y_star.copy(), "z": z_star.copy(),
            "seg": seg_index.copy(), "cost": C_star
        }

        # ── initialise vartheta at current KKT solution ──
        # vartheta[i,j,s] = 1 if s == s*_ij, else 0
        vartheta = {}
        for (i,j) in self.arcs:
            for s in range(1, NP):
                vartheta[(i,j,s)] = 1.0 if s == seg_index[(i,j)] else 0.0
        print(f"\n{'='*65}")
        print(f" DUAL UPDATE SEARCH  |  C* = {C_star:.4f}")
        print(f"{'='*65}")
        for it in range(1, max_outer+1):
            y_inc  = incumbent["y"]
            seg_inc= incumbent["seg"]
            # ══════════════════════════════════════════════════
            # STEP 1: compute update direction for each arc
            # ══════════════════════════════════════════════════
            vartheta_new = {}
            for (i,j) in self.arcs:
                s_cur = seg_inc[(i,j)]
                yv    = y_inc[(i,j)]
                if strategy == 'shift':
                    # transfer theta from s* to s*+1 (cheaper segment)
                    s_next = min(NP-1, s_cur + 1)
                    for s in range(1, NP):
                        vartheta_new[(i,j,s)] = 0.0
                    vartheta_new[(i,j,s_cur)]  = 1.0 - theta
                    vartheta_new[(i,j,s_next)] = theta

            # update vartheta
            vartheta = vartheta_new
            # ══════════════════════════════════════════════════
            # STEP 2: compute modified objective coefficients
            # ══════════════════════════════════════════════════
            # Modified obj:
            # min sum z_ij + sum_{ij,s in I_ij} vartheta*_ijs * (L*(m_s*y+b_s) - z_ij)
            #
            # = min sum z_ij * (1 - sum_{s in I_ij} vartheta*_ijs)
            #     + sum_{ij,s in I_ij} vartheta*_ijs * L_ij * m_s * y_ij
            #     + constant terms (b_s contributions)
            #
            # Coefficients:
            #   coef_z[i,j]   = 1 - sum_s vartheta*_ijs
            #   coef_y[i,j]   = sum_s vartheta*_ijs * L_ij * m_s
            #   constant      = sum_s vartheta*_ijs * L_ij * b_s

            coef_z   = {}
            coef_y   = {}
            constant = {}

            for (i,j) in self.arcs:
                s_cur = seg_inc[(i,j)]
                I_ij  = {s_cur}   # active segment set

                sum_vartheta = sum(vartheta[(i,j,s)] for s in I_ij)
                coef_z[(i,j)] = 1.0 - sum_vartheta

                coef_y[(i,j)] = sum(
                    vartheta[(i,j,s)] * self.L[(i,j)] * self.slope[s]
                    for s in I_ij
                )

                constant[(i,j)] = sum(
                    vartheta[(i,j,s)] * self.L[(i,j)] * self.intercept[s]
                    for s in I_ij
                )

            # ══════════════════════════════════════════════════
            # STEP 3: build and solve relaxed NLP
            # ══════════════════════════════════════════════════
            ampl = AMPL()
            ampl.read("exact_reduced_wdn.mod")
            ampl.read_data(self.data_file)
            ampl.option["solver"] = "ipopt"
            ampl.option["ipopt_options"] = (
                "outlev=0 max_iter=2000 "
                "warm_start_init_point=yes "
            )

            # ampl.eval("""param seg_index{arcs};""")

            # declare modified objective parameters
            ampl.eval("""
                param coef_z{arcs};
                param coef_y{arcs};
                param const_term{arcs};
                param seg_index{arcs};
                param theta default 1;
                drop total_cost;
                minimize relaxed_obj:
                    sum{(i,j) in arcs} (
                        coef_z[i,j] * z[i,j]
                      + coef_y[i,j] * y[i,j]
                      + const_term[i,j]
                    );
                # minimize relaxed_obj:
                #     sum{(i,j) in arcs} (
                #         z[i,j]
                #       + theta*(slope[seg_index[i,j]] * y[i,j]
                #       + intercept[seg_index[i,j]]-z[i,j])
                #     );
            """)

            for (i,j) in self.arcs:
                ampl.param["seg_index"][i,j] = seg_inc[(i,j)]
            
            ampl.param["theta"] = theta

            for (i,j) in self.arcs:
                ampl.param["coef_z"][i,j]    = coef_z[(i,j)]
                ampl.param["coef_y"][i,j]    = coef_y[(i,j)]
                ampl.param["const_term"][i,j]= constant[(i,j)]

            # warm start from incumbent
            for (i,j) in self.arcs:
                ampl.var["q"][i,j] = incumbent["q"][(i,j)]
                ampl.var["y"][i,j] = incumbent["y"][(i,j)]
                ampl.var["z"][i,j] = incumbent["z"][(i,j)]
            for i in self.nodes:
                ampl.var["h"][i] = incumbent["h"][i]

            ampl.eval("drop exact_cost;")
            for (i,j) in self.arcs:
                for s in self.segs:
                    if s != seg_inc[i,j]:
                        ampl.eval(f"s.t. epigraph_{i}_{j}_{s}: z[{i},{j}] >= L[{i},{j}]*(slope[{s}] * y[{i},{j}] + intercept[{s}]);")

            with self.suppress_output():
                ampl.solve()

            status = ampl.get_value("solve_result")

            if status != "solved":
                print(f"[{it:02d}] solver failed — adjusting theta")
                theta = min(1.0, theta * 1.5)
                continue

            q_new  = ampl.get_variable("q").get_values().to_dict()
            y_new  = ampl.get_variable("y").get_values().to_dict()
            h_new  = ampl.get_variable("h").get_values().to_dict()
            z_new  = ampl.get_variable("z").get_values().to_dict()

            seg_new = {(i,j): find_seg(y_new[(i,j)]) for (i,j) in self.arcs}
            seg_changes = sum(
                1 for (i,j) in self.arcs
                if seg_new[(i,j)] != seg_inc[(i,j)]
            )

            # ══════════════════════════════════════════════════
            # STEP 4: check stationarity perturbation
            # ══════════════════════════════════════════════════
            # Verify the perturbation shifted the equilibrium
            print(f"\n[{it:02d}] strategy={strategy}  theta={theta:.3f}")
            print(f"      Segment changes: {seg_changes} arcs")

            # for (i,j) in self.arcs:
            #     if seg_new[(i,j)] != seg_inc[(i,j)]:
            #         print(f"        arc ({i},{j}): "
            #               f"seg {seg_inc[(i,j)]} -> {seg_new[(i,j)]}  "
            #               f"y: {y_inc[(i,j)]:.6f} -> {y_new[(i,j)]:.6f}")

            # ══════════════════════════════════════════════════
            # STEP 5: exact NLP re-solve if segment changed
            # ══════════════════════════════════════════════════
            if seg_changes > 0:
                ampl_ex = AMPL()
                ampl_ex.read("exact_reduced_wdn.mod")
                ampl_ex.read_data(self.data_file)
                ampl_ex.option["solver"] = "ipopt"
                ampl_ex.option["ipopt_options"] = (
                    "outlev=0 max_iter=3000 "
                    "warm_start_init_point=yes "
                )

                for (i,j) in self.arcs:
                    ampl_ex.var["q"][i,j] = q_new[(i,j)]
                    ampl_ex.var["y"][i,j] = y_new[(i,j)]
                    ampl_ex.var["z"][i,j] = z_new[(i,j)]
                for i in self.nodes:
                    ampl_ex.var["h"][i] = h_new[i]
                with self.suppress_output():
                    ampl_ex.solve()

                if ampl_ex.get_value("solve_result") == "solved":
                    q_ex  = ampl_ex.get_variable("q").get_values().to_dict()
                    y_ex  = ampl_ex.get_variable("y").get_values().to_dict()
                    h_ex  = ampl_ex.get_variable("h").get_values().to_dict()
                    z_ex  = ampl_ex.get_variable("z").get_values().to_dict()
                    seg_ex= {(i,j): find_seg(y_ex[(i,j)]) for (i,j) in self.arcs}
                    cost_ex = true_cost(y_ex, seg_ex)

                    print(f"      Exact NLP cost: {cost_ex:.4f}  "
                          f"(incumbent: {incumbent['cost']:.4f})")

                    if cost_ex < incumbent["cost"] - 1e-4:
                        print(f"      *** IMPROVED: {cost_ex:.4f} ***")
                        incumbent = {
                            "q": q_ex.copy(), "h": h_ex.copy(),
                            "y": y_ex.copy(), "z": z_ex.copy(),
                            "seg": seg_ex.copy(), "cost": cost_ex
                        }
                        # reset vartheta at new solution
                        for (i,j) in self.arcs:
                            for s in range(1, NP):
                                vartheta[(i,j,s)] = (
                                    1.0 if s == seg_ex[(i,j)] else 0.0
                                )
                        theta = 0.5   # reset theta after improvement
                    else:
                        # no improvement — increase theta to push harder
                        theta = min(1.0, theta + 0.1)
            else:
                # no segment change — need stronger perturbation
                theta = min(1.0, theta + 0.15)
                print(f"      No segment change — theta increased to {theta:.3f}")

        print(f"\n{'='*65}")
        print(f" FINAL BEST COST: {incumbent['cost']:.4f}")
        print(f"{'='*65}")

        return (incumbent["q"], incumbent["h"],
                incumbent["y"], incumbent["z"])



    def update_dual_and_resolve(self,q_star, h_star, y_star, z_star, seg_index, strategy='shift', theta=0.4, max_outer=30):
        """
        Update vartheta* (dual of active epigraph constraints)
        and re-solve relaxed model to find different local solution.

        Parameters
        ----------
        strategy : 'shift'     — transfer weight theta to adjacent segment
                   'subgrad'   — subgradient update on all segments
                   'target'    — directly assign weight to cheapest segment
        theta    : weight transfer amount (for 'shift' strategy)
        """
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        NP         = len(alpha_vals)

        def find_seg(yv):
            for k in range(NP-1):
                if alpha_vals[k+1]-1e-10 <= yv <= alpha_vals[k]+1e-10:
                    return k+1
            return 1 if yv > alpha_vals[0] else NP-1

        def true_cost(y_sol, seg_index):
            return sum(
                self.L[i,j]*(self.slope[seg_index[i,j]]*y_sol[i,j] + self.intercept[seg_index[i,j]])
                for (i,j) in self.arcs
            )
        # print("seg_index",seg_index) 
        C_star    = true_cost(y_star, seg_index)
        # C_star = sum(self.L[(i,j)] * (self.slope[seg_index[i,j]] * y_sol[(i,j)]+ self.intercept[seg_index[(i,j)]]) for (i,j) in self.arcs)

        incumbent = {
            "q": q_star.copy(), "h": h_star.copy(),
            "y": y_star.copy(), "z": z_star.copy(),
            "seg": seg_index.copy(), "cost": C_star
        }

        # ── initialise vartheta at current KKT solution ──
        # vartheta[i,j,s] = 1 if s == s*_ij, else 0
        vartheta = {}
        for (i,j) in self.arcs:
            for s in range(1, NP):
                vartheta[(i,j,s)] = 1.0 if s == seg_index[(i,j)] else 0.0
        print(f"\n{'='*65}")
        print(f" DUAL UPDATE SEARCH  |  C* = {C_star:.4f}")
        print(f"{'='*65}")
        for it in range(1, max_outer+1):
            y_inc  = incumbent["y"]
            seg_inc= incumbent["seg"]
            # ══════════════════════════════════════════════════
            # STEP 1: compute update direction for each arc
            # ══════════════════════════════════════════════════
            vartheta_new = {}
            for (i,j) in self.arcs:
                s_cur = seg_inc[(i,j)]
                yv    = y_inc[(i,j)]
                if strategy == 'shift':
                    # transfer theta from s* to s*+1 (cheaper segment)
                    s_next = min(NP-1, s_cur + 1)
                    for s in range(1, NP):
                        vartheta_new[(i,j,s)] = 0.0
                    vartheta_new[(i,j,s_cur)]  = 1.0 - theta
                    vartheta_new[(i,j,s_next)] = theta
                elif strategy == 'subgrad':
                    # subgradient: g_ijs = L*(m_s*y + b_s) - z
                    z_cur = self.slope[s_cur]*yv + self.intercept[s_cur]
                    alpha_step = 0.5 / np.sqrt(it)
                    raw = {}
                    for s in range(1, NP):
                        g_ijs = (self.slope[s]*yv + self.intercept[s]) - z_cur
                        raw[(i,j,s)] = max(0.0, vartheta[(i,j,s)] + alpha_step*g_ijs)
                    # normalise to sum = 1
                    total = sum(raw.values()) + 1e-12
                    for s in range(1, NP):
                        vartheta_new[(i,j,s)] = raw[(i,j,s)] / total
                elif strategy == 'target':
                    # find segment with lowest cost at current y*
                    # (excluding current active segment)
                    best_s    = s_cur
                    best_cost = float('inf')
                    for s in range(1, NP):
                        if s == s_cur:
                            continue
                        c_s = self.slope[s]*yv + self.intercept[s]
                        if c_s < best_cost:
                            best_cost = c_s
                            best_s    = s
                    for s in range(1, NP):
                        vartheta_new[(i,j,s)] = 0.0
                    vartheta_new[(i,j,best_s)] = 1.0
            # update vartheta
            vartheta = vartheta_new
            # ══════════════════════════════════════════════════
            # STEP 2: compute modified objective coefficients
            # ══════════════════════════════════════════════════
            # Modified obj:
            # min sum z_ij + sum_{ij,s in I_ij} vartheta*_ijs * (L*(m_s*y+b_s) - z_ij)
            #
            # = min sum z_ij * (1 - sum_{s in I_ij} vartheta*_ijs)
            #     + sum_{ij,s in I_ij} vartheta*_ijs * L_ij * m_s * y_ij
            #     + constant terms (b_s contributions)
            #
            # Coefficients:
            #   coef_z[i,j]   = 1 - sum_s vartheta*_ijs
            #   coef_y[i,j]   = sum_s vartheta*_ijs * L_ij * m_s
            #   constant      = sum_s vartheta*_ijs * L_ij * b_s

            coef_z   = {}
            coef_y   = {}
            constant = {}

            for (i,j) in self.arcs:
                s_cur = seg_inc[(i,j)]
                I_ij  = {s_cur}   # active segment set

                sum_vartheta = sum(vartheta[(i,j,s)] for s in I_ij)
                coef_z[(i,j)] = 1.0 - sum_vartheta

                coef_y[(i,j)] = sum(
                    vartheta[(i,j,s)] * self.L[(i,j)] * self.slope[s]
                    for s in I_ij
                )

                constant[(i,j)] = sum(
                    vartheta[(i,j,s)] * self.L[(i,j)] * self.intercept[s]
                    for s in I_ij
                )

            # ══════════════════════════════════════════════════
            # STEP 3: build and solve relaxed NLP
            # ══════════════════════════════════════════════════
            ampl = AMPL()
            ampl.read("exact_reduced_wdn.mod")
            ampl.read_data(self.data_file)
            ampl.option["solver"] = "ipopt"
            ampl.option["ipopt_options"] = (
                "outlev=0 max_iter=2000 "
                "warm_start_init_point=yes "
                "bound_relax_factor=0 "
            )

            # ampl.eval("""param seg_index{arcs};""")

            # declare modified objective parameters
            ampl.eval("""
                param coef_z{arcs};
                param coef_y{arcs};
                param const_term{arcs};
                param seg_index{arcs};
                param theta default 1;
                drop total_cost;
                minimize relaxed_obj:
                    sum{(i,j) in arcs} (
                        coef_z[i,j] * z[i,j]
                      + coef_y[i,j] * y[i,j]
                      + const_term[i,j]
                    );
                # minimize relaxed_obj:
                #     sum{(i,j) in arcs} (
                #         z[i,j]
                #       + theta*(slope[seg_index[i,j]] * y[i,j]
                #       + intercept[seg_index[i,j]]-z[i,j])
                #     );
            """)
            for (i,j) in self.arcs:
                ampl.param["seg_index"][i,j] = seg_inc[(i,j)]
            
            ampl.param["theta"] = theta

            for (i,j) in self.arcs:
                ampl.param["coef_z"][i,j]    = coef_z[(i,j)]
                ampl.param["coef_y"][i,j]    = coef_y[(i,j)]
                ampl.param["const_term"][i,j]= constant[(i,j)]

            # warm start from incumbent
            for (i,j) in self.arcs:
                ampl.var["q"][i,j] = incumbent["q"][(i,j)]
                ampl.var["y"][i,j] = incumbent["y"][(i,j)]
                ampl.var["z"][i,j] = incumbent["z"][(i,j)]
            for i in self.nodes:
                ampl.var["h"][i] = incumbent["h"][i]

            ampl.eval("drop exact_cost;")
            # for (i,j) in self.arcs:
            #     for s in self.segs:
            #         if s != seg_inc[i,j]:
            #             ampl.eval(f"s.t. epigraph_{i}_{j}_{s}: z[{i},{j}] >= L[{i},{j}]*(slope[{s}] * y[{i},{j}] + intercept[{s}]);")

            with self.suppress_output():
                ampl.solve()

            status = ampl.get_value("solve_result")

            if status != "solved":
                print(f"[{it:02d}] solver failed — adjusting theta")
                theta = min(1.0, theta * 1.5)
                continue

            q_new  = ampl.get_variable("q").get_values().to_dict()
            y_new  = ampl.get_variable("y").get_values().to_dict()
            h_new  = ampl.get_variable("h").get_values().to_dict()
            z_new  = ampl.get_variable("z").get_values().to_dict()

            seg_new = {(i,j): find_seg(y_new[(i,j)]) for (i,j) in self.arcs}
            seg_changes = sum(
                1 for (i,j) in self.arcs
                if seg_new[(i,j)] != seg_inc[(i,j)]
            )

            # ══════════════════════════════════════════════════
            # STEP 4: check stationarity perturbation
            # ══════════════════════════════════════════════════
            # Verify the perturbation shifted the equilibrium
            print(f"\n[{it:02d}] strategy={strategy}  theta={theta:.3f}")
            print(f"      Segment changes: {seg_changes} arcs")

            # for (i,j) in self.arcs:
            #     if seg_new[(i,j)] != seg_inc[(i,j)]:
            #         print(f"        arc ({i},{j}): "
            #               f"seg {seg_inc[(i,j)]} -> {seg_new[(i,j)]}  "
            #               f"y: {y_inc[(i,j)]:.6f} -> {y_new[(i,j)]:.6f}")

            # ══════════════════════════════════════════════════
            # STEP 5: exact NLP re-solve if segment changed
            # ══════════════════════════════════════════════════
            if seg_changes > 0:
                ampl_ex = AMPL()
                ampl_ex.read("exact_reduced_wdn.mod")
                ampl_ex.read_data(self.data_file)
                ampl_ex.option["solver"] = "ipopt"
                ampl_ex.option["ipopt_options"] = (
                    "outlev=0 max_iter=3000 "
                    "warm_start_init_point=yes "
                    "bound_relax_factor=0 "
                )

                for (i,j) in self.arcs:
                    ampl_ex.var["q"][i,j] = q_new[(i,j)]
                    ampl_ex.var["y"][i,j] = y_new[(i,j)]
                    ampl_ex.var["z"][i,j] = z_new[(i,j)]
                for i in self.nodes:
                    ampl_ex.var["h"][i] = h_new[i]
                with self.suppress_output():
                    ampl_ex.solve()

                if ampl_ex.get_value("solve_result") == "solved":
                    q_ex  = ampl_ex.get_variable("q").get_values().to_dict()
                    y_ex  = ampl_ex.get_variable("y").get_values().to_dict()
                    h_ex  = ampl_ex.get_variable("h").get_values().to_dict()
                    z_ex  = ampl_ex.get_variable("z").get_values().to_dict()
                    seg_ex= {(i,j): find_seg(y_ex[(i,j)]) for (i,j) in self.arcs}
                    cost_ex = true_cost(y_ex, seg_ex)

                    print(f"      Exact NLP cost: {cost_ex:.4f}  "
                          f"(incumbent: {incumbent['cost']:.4f})")

                    if cost_ex < incumbent["cost"] - 1e-4:
                        print(f"      *** IMPROVED: {cost_ex:.4f} ***")
                        incumbent = {
                            "q": q_ex.copy(), "h": h_ex.copy(),
                            "y": y_ex.copy(), "z": z_ex.copy(),
                            "seg": seg_ex.copy(), "cost": cost_ex
                        }
                        # reset vartheta at new solution
                        for (i,j) in self.arcs:
                            for s in range(1, NP):
                                vartheta[(i,j,s)] = (
                                    1.0 if s == seg_ex[(i,j)] else 0.0
                                )
                        theta = 0.5   # reset theta after improvement
                    else:
                        # no improvement — increase theta to push harder
                        theta = min(1.0, theta + 0.1)
            else:
                # no segment change — need stronger perturbation
                theta = min(1.0, theta + 0.15)
                print(f"      No segment change — theta increased to {theta:.3f}")

        print(f"\n{'='*65}")
        print(f" FINAL BEST COST: {incumbent['cost']:.4f}")
        print(f"{'='*65}")

        return (incumbent["q"], incumbent["h"],
                incumbent["y"], incumbent["z"])

    def escape_and_project(
        self,
        q_loc,
        h_loc,
        y_loc,
        z_loc,
        seg_index_loc,
        max_outer_iter=1,
        infeas_tol=1e-5,
        proj_tol=1e-6,
        cost_margin=1e-3
    ):
        """
        ------------------------------------------------------------
        MULTI-BASIN FEASIBILITY ESCAPE ALGORITHM
        ------------------------------------------------------------
        PHASE A:
            Infeasibility detection model
            -> pushes solution outside current basin
            -> enforces reduced cost
        PHASE B:
            Projection model
            -> projects infeasible point to feasible manifold
            -> preserves reduced cost
            -> avoids convergence to same local point
        PHASE C:
            Exact NLP solve
            -> obtains improved local solution
        ------------------------------------------------------------
        """
        best_cost = sum(
            z_loc[i, j]
            for (i, j) in self.arcs
        )
        incumbent = {
            "q": q_loc.copy(),
            "h": h_loc.copy(),
            "y": y_loc.copy(),
            "z": z_loc.copy(),
            "cost": best_cost
        }
        alpha_vals = sorted(
            set(self.alpha.values()),
            reverse=True
        )
        print("\n===================================================")
        print(" MULTI-BASIN ESCAPE ALGORITHM ")
        print("===================================================")
        for OUTER in range(1, max_outer_iter + 1):
            print("\n---------------------------------------------------")
            print(f" OUTER ITERATION : {OUTER}")
            print("---------------------------------------------------")
            target_cost = incumbent["cost"] - 0.00001*incumbent["cost"]
            # =========================================================
            # =========================================================
            # PHASE A:
            # INFEASIBILITY DETECTION MODEL
            # =========================================================
            # =========================================================
            ampl = AMPL()
            ampl.read("repair_feasibility_wdn.mod")
            ampl.read_data(self.data_file)
            ampl.option["solver"] = "ipopt"
            ampl.option["ipopt_options"] = (
                "outlev=0 "
                "bound_relax_factor=0 "
                "warm_start_init_point=yes "
                f"bound_push={self.bound_push} "
                f"bound_frac={self.bound_frac} "
            )
            # ---------------------------------------------------------
            # PARAMETERS
            # ---------------------------------------------------------
            ampl.eval("param alpha_arc{arcs};")
            ampl.eval("param y_ref{arcs};")
            for (i, j) in self.arcs:
                ampl.param["alpha_arc"][i, j] = self.alpha[seg_index_loc[i,j]]
                ampl.param["y_ref"][i, j] = incumbent["y"][i, j]
            # ---------------------------------------------------------
            # WARM START
            # ---------------------------------------------------------
            for (i, j), val in incumbent["q"].items():
                ampl.var["q"][i, j] = val
            for (i, j), val in incumbent["y"].items():
                ampl.var["y"][i, j] = val
            for i, val in incumbent["h"].items():
                ampl.var["h"][i] = val
            for (i, j), val in incumbent["z"].items():
                ampl.var["z"][i, j] = val
            # ---------------------------------------------------------
            # OBJECTIVE
            # ---------------------------------------------------------
            ampl.eval("""
            minimize infeasibility_detection:
                  1e-1 * sum{(i,j) in arcs}
                        delta_headloss[i,j]^2
                + 1e-1 * sum{(i,j) in arcs}
                        delta_y[i,j]^2
                - 1e-1 * sum{(i,j) in arcs}
                        (y[i,j] - y_ref[i,j])^2;
            """)
            # ---------------------------------------------------------
            # REDUCED COST CONSTRAINT
            # ---------------------------------------------------------
            ampl.eval(f"""
            subject to reduced_cost:
                sum{{(i,j) in arcs}} z[i,j]
                <= {target_cost};
            """)
            ampl.eval("""
                subject to bound_y{(i,j) in arcs}:
                    y[i,j] >= alpha_arc[i,j] - delta_y[i,j];
            """)
            # ---------------------------------------------------------
            # ESCAPE CUT
            # ---------------------------------------------------------
            factors = []
            for (i, j) in self.arcs:
                s = seg_index_loc[(i, j)]
                a_up = alpha_vals[s - 1]
                a_lo = alpha_vals[s]
                factors.append(
                    f"(({a_up} - y[{i},{j}])"
                    f"*(y[{i},{j}] - {a_lo}))"
                )
            product_expr = " * ".join(factors)
            # ampl.eval(f"""
            #
            # subject to escape_cut:
            #     {product_expr} <= 0;
            #
            # """)
            # ---------------------------------------------------------
            # SOLVE INFEASIBILITY MODEL
            # ---------------------------------------------------------
            print("\nSolving infeasibility detection model ...")
            with self.suppress_output():
                ampl.solve()
            solve_status = ampl.get_value("solve_result")
            if solve_status != "solved":
                print("Infeasibility model failed.")
                continue
            q_inf = ampl.get_variable(
                "q"
            ).get_values().to_dict()
            y_inf = ampl.get_variable(
                "y"
            ).get_values().to_dict()
            h_inf = ampl.get_variable(
                "h"
            ).get_values().to_dict()
            z_inf = ampl.get_variable(
                "z"
            ).get_values().to_dict()
            delta_head = ampl.get_variable(
                "delta_headloss"
            ).get_values().to_dict()
            delta_y = ampl.get_variable(
                "delta_y"
            ).get_values().to_dict()
            head_inf = max(
                abs(delta_head[i, j])
                for (i, j) in self.arcs
            )
            y_inf_norm = max(
                abs(delta_y[i, j])
                for (i, j) in self.arcs
            )
            det_cost = sum(
                z_inf[i, j]
                for (i, j) in self.arcs
            )
            print(
                f"Detection infeasibility = "
                f"{head_inf + y_inf_norm:.3e}"
            )
            print(
                f"Detection cost = {det_cost:.4f}"
            )
            # =========================================================
            # =========================================================
            # PHASE B:
            # PROJECTION TO FEASIBLE REGION
            # =========================================================
            # =========================================================
            ampl2 = AMPL()
            ampl2.read("exact_reduced_wdn.mod")
            ampl2.read_data(self.data_file)
            ampl2.option["solver"] = "ipopt"
            ampl2.option["ipopt_options"] = (
                "outlev=0 "
                "max_iter=3000 "
                "bound_relax_factor=0 "
                "warm_start_init_point=yes "
                f"bound_push={self.bound_push} "
                f"bound_frac={self.bound_frac} "

            )
            # ---------------------------------------------------------
            # WARM START FROM INFEASIBLE POINT
            # ---------------------------------------------------------
            for (i, j), val in q_inf.items():
                ampl2.var["q"][i, j] = val
            for (i, j), val in y_inf.items():
                ampl2.var["y"][i, j] = val
            for i, val in h_inf.items():
                ampl2.var["h"][i] = val
            for (i, j), val in z_inf.items():
                ampl2.var["z"][i, j] = val
            # ---------------------------------------------------------
            # STAY AWAY FROM PREVIOUS LOCAL SOLUTION
            # ---------------------------------------------------------
            ampl2.eval("""
            param y_prev{arcs};

            param y_inf{arcs};
            param z_inf{arcs};
            param q_inf{arcs};
            param h_inf{nodes};
            """)
            for (i, j) in self.arcs:
                ampl2.param["y_prev"][i, j] = incumbent["y"][i, j]
                ampl2.param["y_inf"][i, j] = y_inf[i, j]
                ampl2.param["q_inf"][i, j] = q_inf[i, j]
                ampl2.param["z_inf"][i, j] = z_inf[i, j]
            for i in self.nodes:
                ampl2.param["h_inf"][i] = h_inf[i]

            ampl2.eval("drop total_cost;")
            ampl2.eval("""
            minimize projection_obj:
                  sum{(i,j) in arcs} (z[i,j] + (z[i,j] - z_inf[i,j])^2 + (y[i,j] - y_inf[i,j])^2 + (q[i,j] - q_inf[i,j])^2) + sum{i in nodes} (h[i] - h_inf[i])^2;
            """)
            # ampl2.eval(f"""
            # subject to reduced_cost_proj:
            #     sum{{(i,j) in arcs}} z[i,j]
            #     <= {target_cost};
            # """)
            # ---------------------------------------------------------
            # SOLVE PROJECTION MODEL
            # ---------------------------------------------------------
            print("\nProjecting to feasible region ...")
            with self.suppress_output():
                ampl2.solve()
            solve_status = ampl2.get_value(
                "solve_result"
            )
            if solve_status != "solved":
                print("Projection failed.")
                continue
            else:
                print("Optimal Solution Found")
            q_proj = ampl2.get_variable(
                "q"
            ).get_values().to_dict()
            y_proj = ampl2.get_variable(
                "y"
            ).get_values().to_dict()
            h_proj = ampl2.get_variable(
                "h"
            ).get_values().to_dict()
            z_proj = ampl2.get_variable(
                "z"
            ).get_values().to_dict()
            proj_cost = sum(
                z_proj[i, j]
                for (i, j) in self.arcs
            )
            print(
                f"Projected feasible cost = "
                f"{proj_cost:.4f}"
            )
            # =========================================================
            # =========================================================
            # PHASE C:
            # EXACT LOCAL NLP
            # =========================================================
            # =========================================================
            ampl3 = AMPL()
            ampl3.read("exact_reduced_wdn.mod")
            ampl3.read_data(self.data_file)
            ampl3.option["solver"] = "ipopt"
            ampl3.option["ipopt_options"] = (
                "outlev=0 "
                "bound_relax_factor=0 "
                "warm_start_init_point=yes "
                f"bound_push={self.bound_push} "
                f"bound_frac={self.bound_frac} "

            )
            # ---------------------------------------------------------
            # WARM START
            # ---------------------------------------------------------
            for (i, j), val in q_proj.items():
                ampl3.var["q"][i, j] = val
            for (i, j), val in y_proj.items():
                ampl3.var["y"][i, j] = val
            for i, val in h_proj.items():
                ampl3.var["h"][i] = val
            for (i, j), val in z_proj.items():
                ampl3.var["z"][i, j] = val

            print("\nSolving exact NLP ...")

            with self.suppress_output():
                ampl3.solve()

            solve_status = ampl3.get_value(
                "solve_result"
            )
            if solve_status != "solved":
                print("Exact NLP failed.")
                continue
            q_new = ampl3.get_variable(
                "q"
            ).get_values().to_dict()
            y_new = ampl3.get_variable(
                "y"
            ).get_values().to_dict()
            h_new = ampl3.get_variable(
                "h"
            ).get_values().to_dict()
            z_new = ampl3.get_variable(
                "z"
            ).get_values().to_dict()
            new_cost = sum(
                z_new[i, j]
                for (i, j) in self.arcs
            )
            print(
                f"New local solution cost = "
                f"{new_cost:.4f}"
            )
            # =========================================================
            # ACCEPT IMPROVED SOLUTION
            # =========================================================
            if new_cost < incumbent["cost"] - 1e-6:
                print("\nImproved local solution found.")
                incumbent = {
                    "q": q_new.copy(),
                    "h": h_new.copy(),
                    "y": y_new.copy(),
                    "z": z_new.copy(),
                    "cost": new_cost
                }
                # ---------------------------------------------
                # UPDATE SEGMENTS
                # ---------------------------------------------
                for (i, j) in self.arcs:
                    yv = y_new[i, j]
                    for k in range(
                        len(alpha_vals) - 1
                    ):
                        if (
                            alpha_vals[k + 1]
                            <= yv
                            <= alpha_vals[k]
                        ):
                            seg_index_loc[(i, j)] = k + 1
                            break
                q_sol, h_sol, y_sol, z_sol = self.escape_and_project(q_new, h_new, y_new, z_new, seg_index_loc)
            else:
                print(
                    "\nNo improving feasible basin found."
                )
                break
        print("\n===================================================")
        print(" FINAL BEST COST :", incumbent["cost"])
        print("===================================================")
        
        return (
            incumbent["q"],
            incumbent["h"],
            incumbent["y"],
            incumbent["z"]
        )

    def solve_restricted_model(self, q_sol, z_sol, y_sol, h_sol, seg_index):
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        NP         = len(alpha_vals)
    
        # ── build ONE persistent AMPL instance ──
        ampl = AMPL()
        ampl.read("exact_reduced_wdn.mod")
        ampl.read_data(self.data_file)
        ampl.option["solver"] = "ipopt"
    
        # ── declare parametric bounds and cost target ──
        ampl.eval("""
            param y_lb{arcs} default alpha_min;
            param y_ub{arcs} default alpha_max;
            param cost_target;
    
            subject to local_lb{(i,j) in arcs}:
                y[i,j] >= y_lb[i,j];
            subject to local_ub{(i,j) in arcs}:
                y[i,j] <= y_ub[i,j];
            subject to cost_cut:
                sum{(i,j) in arcs} z[i,j] <= cost_target;
        """)

        # ── tracking ──
        best_cost    = self.best_cost
        best_q       = q_sol.copy()
        best_y       = y_sol.copy()
        best_h       = h_sol.copy()
        best_z       = z_sol.copy()
        best_seg     = seg_index.copy()
    
        seg_cur      = seg_index.copy()
        y_cur        = y_sol.copy()
        q_cur        = q_sol.copy()
        h_cur        = h_sol.copy()
        z_cur        = z_sol.copy()
    
        no_improve   = 0
        radius       = 1          # segment window half-width
        MAX_RADIUS   = 2
        MAX_ITER     = 60
        MAX_NO_IMP   = 8
        FEAS_TOL     = 1e-6
    
        print(f"\n{'iter':>4}  {'status':<12}  {'head_inf':>10}  "
              f"{'cost':>12}  {'best':>12}  {'radius':>6}")
    
        for it in range(1, MAX_ITER + 1):
    
            # ── set local branching window around current segments ──
            for (i,j) in self.arcs:
                s    = seg_cur[(i,j)]
                s_lo = max(1,    s - radius)
                s_hi = min(NP-1, s + radius)
                ampl.param["y_lb"][i,j] = alpha_vals[s_hi]     # lower y
                ampl.param["y_ub"][i,j] = alpha_vals[s_lo - 1] # upper y
    
            # ── set cost target slightly below current best ──
            ampl.param["cost_target"] = best_cost - 1e-3
    
            # ── warm start from current point ──
            ampl.option["ipopt_options"] = (
                "outlev=0 tol=1e-8 max_iter=800 "
                "warm_start_init_point=yes "
            )
            for (i,j) in self.arcs:
                ampl.var["q"][i,j] = q_cur[(i,j)]
                ampl.var["y"][i,j] = max(
                    ampl.param["y_lb"][i,j],
                    min(ampl.param["y_ub"][i,j], y_cur[(i,j)])
                )
            for i in self.nodes:
                ampl.var["h"][i] = h_cur[i]
            for (i,j) in self.arcs:
                ampl.var["z"][i,j] = z_cur[(i,j)]
    
            with self.suppress_output():
                ampl.solve()
    
            status = ampl.get_value("solve_result")
    
            q_it = ampl.get_variable("q").get_values().to_dict()
            y_it = ampl.get_variable("y").get_values().to_dict()
            h_it = ampl.get_variable("h").get_values().to_dict()
            z_it = ampl.get_variable("z").get_values().to_dict()

            # ── compute exact headloss residual ──
            def f_hw(q, eps):
                return q**3*(q**2+eps**2)**0.426/(q**2+0.426*eps**2)
    
            head_inf = max(
                abs(h_it[i] - h_it[j]
                    - f_hw(q_it[(i,j)], self.eps[(i,j)])
                    * y_it[(i,j)] * self.L[(i,j)])
                for (i,j) in self.arcs
            )
    
            # ── compute cost from PL envelope ──
            seg_it = {}
            for (i,j) in self.arcs:
               yv = y_sol[(i,j)]
               for k in range(len(alpha_vals)-1):
                   if (
                       alpha_vals[k+1]
                       <= yv
                       <= alpha_vals[k]
                   ):
                       seg_it[(i,j)] = k+1
                       break

            cost_it = sum(
                self.L[(i,j)] * (
                    self.slope[seg_it[i,j]] * y_it[(i,j)]
                    + self.intercept[seg_it[(i,j)]]
                )
                for (i,j) in self.arcs
            )
    
            feasible = (status == "solved") and (head_inf <= FEAS_TOL)
    
            print(f"{it:>4}  {status:<12}  {head_inf:>10.3e}  "
                  f"{cost_it:>12.4f}  {best_cost:>12.4f}  {radius:>6}")
    
            # ── update best if improved feasible solution found ──
            if feasible and cost_it < best_cost - 1e-4:
                best_cost = cost_it
                best_q    = q_it.copy()
                best_y    = y_it.copy()
                best_h    = h_it.copy()
                best_z    = z_it.copy()
                best_seg  = seg_it.copy()
                no_improve = 0
                radius     = 1   # reset radius after improvement
                print(f"  *** Improved: {cost_it:.4f} ***")
    
                # update current point to best
                q_cur   = q_it.copy()
                y_cur   = y_it.copy()
                h_cur   = h_it.copy()
                z_cur   = z_it.copy()
                seg_cur = seg_it.copy()
    
            else:
                no_improve += 1
    
                # ── adaptive response to no improvement ──
                if status == "solved" and head_inf <= FEAS_TOL:
                    # feasible but not cheaper — try expanding radius
                    radius = min(MAX_RADIUS, radius + 1)
                    # update current point even if not better cost
                    q_cur   = q_it.copy()
                    y_cur   = y_it.copy()
                    h_cur   = h_it.copy()
                    z_cur   = z_it.copy()
                    seg_cur = seg_it.copy()
    
                elif status != "solved":
                    # infeasible — use Newton correction and stay at current
                    for (i,j) in self.arcs:
                        qv    = q_it[(i,j)]
                        if abs(qv) < 1e-8:
                            continue
                        phi   = f_hw(qv, self.eps[(i,j)])
                        resid = (h_it[i] - h_it[j]
                                 - phi * y_it[(i,j)] * self.L[(i,j)])
                        step  = resid / (phi * self.L[(i,j)] + 1e-12)
                        # move y in correction direction, stay in window
                        y_new = y_cur[(i,j)] + 0.3 * step
                        lb    = ampl.param["y_lb"][i,j]
                        ub    = ampl.param["y_ub"][i,j]
                        y_cur[(i,j)] = max(lb, min(ub, y_new))
    
                if no_improve >= MAX_NO_IMP:
                    if radius < MAX_RADIUS:
                        # try wider window before giving up
                        radius += 1
                        no_improve = 0
                        print(f"  Expanding radius to {radius}")
                    else:
                        print(f"  No improvement after {MAX_NO_IMP} iters. Stopping.")
                        break

        # ════════════════════════════════════════════════════════════
        # FINAL EXACT SOLVE of original NLP from best point found
        # ════════════════════════════════════════════════════════════
        print("\n--- Final exact solve of original NLP ---")
    
        ampl_final = AMPL()
        ampl_final.read("exact_reduced_wdn.mod")   # original model, no extra constraints
        ampl_final.read_data(self.data_file)
        ampl_final.option["solver"] = "ipopt"
        ampl_final.option["ipopt_options"] = (
            "outlev=0 tol=1e-9 max_iter=3000 "
            "warm_start_init_point=yes "
        )
    
        # warm start from best point
        for (i,j) in self.arcs:
            ampl_final.var["q"][i,j] = best_q[(i,j)]
            ampl_final.var["y"][i,j] = best_y[(i,j)]
            ampl_final.var["z"][i,j] = best_z[(i,j)]
        for i in self.nodes:
            ampl_final.var["h"][i] = best_h[i]
    
        with self.suppress_output():
            ampl_final.solve()
    
        final_status = ampl_final.get_value("solve_result")
    
        q_f = ampl_final.get_variable("q").get_values().to_dict()
        y_f = ampl_final.get_variable("y").get_values().to_dict()
        h_f = ampl_final.get_variable("h").get_values().to_dict()
        z_f = ampl_final.get_variable("z").get_values().to_dict()
    
        # final segment index
        for (i,j) in self.arcs:
           yv = y_sol[(i,j)]
           for k in range(len(alpha_vals)-1):
               if (
                   alpha_vals[k+1]
                   <= yv
                   <= alpha_vals[k]
               ):
                   seg_index[(i,j)] = k+1
                   break

        # for (i,j) in self.arcs:
        #     yv = y_f[(i,j)]
        #     for k in range(NP - 1):
        #         if alpha_vals[k+1] <= yv <= alpha_vals[k]:
        #             seg_index[(i,j)] = k + 1
        #             break
    
        final_cost = sum(z_f.values())
        final_hl   = max(
            abs(h_f[i] - h_f[j]
                - f_hw(q_f[(i,j)], self.eps[(i,j)])
                * y_f[(i,j)] * self.L[(i,j)])
            for (i,j) in self.arcs
        )
    
        print(f"  Final status:      {final_status}")
        print(f"  Final cost:        {final_cost:.6f}")
        print(f"  Final headloss inf:{final_hl:.3e}")
        print(f"  Improvement over C*: {self.best_cost - final_cost:.4f}")
    
        return q_f, h_f, y_f, z_f

    def solve_restricted_model7(
        self,
        q_sol,
        z_sol,
        y_sol,
        h_sol,
        seg_index
    ):
    
        import random
        import numpy as np
    
        # ============================================================
        # PARAMETERS
        # ============================================================
        MAX_ITER = 80
    
        FEAS_TOL = 1e-6
    
        MAX_CRITICAL_ARCS = 3
    
        LOCAL_BRANCH_RADIUS = 1
    
        MAX_NO_IMPROVE = 8
    
        TRUST_REGION = 0.03 * (
            self.alpha_max - self.alpha_min
        )
    
        # ============================================================
        # INITIALIZATION
        # ============================================================
        alpha_vals = sorted(
            set(self.alpha.values()),
            reverse=True
        )
    
        best_local_cost = float("inf")
    
        best_q = q_sol.copy()
        best_y = y_sol.copy()
        best_h = h_sol.copy()
        best_z = z_sol.copy()
    
        no_improve_counter = 0
    
        # ============================================================
        # MAIN LOOP
        # ============================================================
        print(
            "\niter    head_inf      "
            "critical_arcs      "
            "cost"
        )
    
        for it in range(1, MAX_ITER + 1):
    
            # ========================================================
            # BUILD NLP MODEL
            # ========================================================
            ampl = AMPL()
    
            ampl.read("exact_reduced_wdn.mod")
    
            ampl.read_data(self.data_file)
    
            ampl.option["solver"] = "ipopt"
    
            ampl.option["ipopt_options"] = (
                "outlev=0 "
                "tol=1e-7 "
                "acceptable_tol=1e-5 "
                "max_iter=500 "
                "mu_strategy=adaptive "
                "warm_start_init_point=yes "
                "bound_push=1e-3 "
                "bound_frac=1e-3 "
                "nlp_scaling_method=gradient-based "
                "print_level=0 "
            )
    
            # ========================================================
            # WARM START
            # ========================================================
            for (i,j), val in q_sol.items():
                ampl.var["q"][i,j] = val
    
            for (i,j), val in y_sol.items():
                ampl.var["y"][i,j] = val
    
            for i, val in h_sol.items():
                ampl.var["h"][i] = val
    
            for (i,j), val in z_sol.items():
                ampl.var["z"][i,j] = val
    
            # ========================================================
            # LOCAL BRANCHING RESTRICTION
            # ========================================================
            ampl.eval("""
                param y_lb{arcs};
                param y_ub{arcs};
            """)
    
            for (i,j) in self.arcs:
    
                s = seg_index[(i,j)]
    
                s_lo = max(
                    1,
                    s - LOCAL_BRANCH_RADIUS
                )
    
                s_up = min(
                    len(alpha_vals)-1,
                    s + LOCAL_BRANCH_RADIUS
                )
    
                y_lb = alpha_vals[s_up]
                y_ub = alpha_vals[s_lo-1]
    
                ampl.param["y_lb"][i,j] = y_lb
                ampl.param["y_ub"][i,j] = y_ub
    
            ampl.eval("""
    
                subject to local_branch_lb{(i,j) in arcs}:
                    y[i,j] >= y_lb[i,j];
    
                subject to local_branch_ub{(i,j) in arcs}:
                    y[i,j] <= y_ub[i,j];
            """)
    
            # ========================================================
            # COST CUT
            # ========================================================
            ampl.eval(f"""
    
                subject to cost_cut:
    
                    sum{{(i,j) in arcs}} z[i,j]
    
                    <= {self.best_cost - 1e-5};
    
            """)
    
            # ========================================================
            # SOLVE NLP
            # ========================================================
            with self.suppress_output():
                ampl.solve()
    
            solve_status = ampl.get_value(
                "solve_result"
            )
    
            if solve_status != "solved":
    
                no_improve_counter += 1
    
                if no_improve_counter >= MAX_NO_IMPROVE:
                    break
    
                continue
    
            # ========================================================
            # EXTRACT SOLUTION
            # ========================================================
            q_sol = ampl.get_variable(
                "q"
            ).get_values().to_dict()
    
            y_sol = ampl.get_variable(
                "y"
            ).get_values().to_dict()
    
            h_sol = ampl.get_variable(
                "h"
            ).get_values().to_dict()
    
            z_sol = ampl.get_variable(
                "z"
            ).get_values().to_dict()
    
            # ========================================================
            # COMPUTE HEADLOSS VIOLATION
            # ========================================================
            delta_head = {}
    
            for (i,j) in self.arcs:
    
                qv = q_sol[(i,j)]
    
                phi = (
                    qv**3
                    * (qv**2 + self.eps[i,j]**2)**0.426
                ) / (
                    qv**2
                    + 0.426*self.eps[i,j]**2
                )
    
                delta = (
                    h_sol[i]
                    - h_sol[j]
                    - phi
                    * y_sol[(i,j)]
                    * self.L[i,j]
                )
    
                delta_head[(i,j)] = delta
    
            head_inf = max(
                abs(delta_head[a])
                for a in self.arcs
            )
    
            # ========================================================
            # UPDATE ACTIVE SEGMENTS
            # ========================================================
            for (i,j) in self.arcs:
    
                yv = y_sol[(i,j)]
    
                for k in range(len(alpha_vals)-1):
    
                    if (
                        alpha_vals[k+1]
                        <= yv
                        <= alpha_vals[k]
                    ):
    
                        seg_index[(i,j)] = k+1
                        break
    
            # ========================================================
            # COST
            # ========================================================
            cost = sum(
                self.L[i,j] * (
                    self.slope[seg_index[(i,j)]]
                    * y_sol[(i,j)]
                    +
                    self.intercept[
                        seg_index[(i,j)]
                    ]
                )
                for (i,j) in self.arcs
            )
    
            print(
                f"{it:03d}      "
                f"{head_inf:.3e}      "
                f"{MAX_CRITICAL_ARCS}      "
                f"{cost:.2f}"
            )
    
            # ========================================================
            # STORE BEST SOLUTION
            # ========================================================
            if (
                head_inf <= FEAS_TOL
                and
                cost < best_local_cost
            ):
    
                best_local_cost = cost
    
                best_q = q_sol.copy()
                best_y = y_sol.copy()
                best_h = h_sol.copy()
                best_z = z_sol.copy()
    
                no_improve_counter = 0
    
                print(
                    f"\nImproved feasible solution:"
                    f" {cost:.2f}"
                )
    
            else:
    
                no_improve_counter += 1
    
            # ========================================================
            # TERMINATION
            # ========================================================
            if no_improve_counter >= MAX_NO_IMPROVE:
    
                print(
                    "\nNo improvement."
                )
    
                break
    
            # ========================================================
            # SENSITIVITY RANKING
            # ========================================================
            ranked_arcs = []
    
            for (i,j) in self.arcs:
    
                qv = q_sol[(i,j)]
    
                phi = abs(
                    (
                        qv**3
                        * (qv**2 + self.eps[i,j]**2)**0.426
                    ) / (
                        qv**2
                        + 0.426*self.eps[i,j]**2
                    )
                )
    
                sensitivity = max(
                    1e-8,
                    phi * self.L[i,j]
                )
    
                score = abs(
                    delta_head[(i,j)]
                ) / sensitivity
    
                ranked_arcs.append(
                    (
                        score,
                        (i,j)
                    )
                )
    
            ranked_arcs.sort(reverse=True)
    
            critical_arcs = [
                arc
                for _, arc in ranked_arcs[
                    :MAX_CRITICAL_ARCS
                ]
            ]
    
            # ========================================================
            # ANALYTICAL FEASIBILITY CORRECTION
            # ========================================================
            for (i,j) in critical_arcs:
    
                qv = q_sol[(i,j)]
    
                if abs(qv) <= 1e-8:
                    continue
    
                phi = (
                    qv**3
                    * (qv**2 + self.eps[i,j]**2)**0.426
                ) / (
                    qv**2
                    + 0.426*self.eps[i,j]**2
                )
    
                step = (
                    delta_head[(i,j)]
                    / (
                        phi
                        * self.L[i,j]
                    )
                )
    
                y_new = (
                    y_sol[(i,j)]
                    + 0.25 * step
                )
    
                # ====================================================
                # TRUST REGION
                # ====================================================
                y_new = max(
                    y_sol[(i,j)] - TRUST_REGION,
                    min(
                        y_sol[(i,j)] + TRUST_REGION,
                        y_new
                    )
                )
    
                # ====================================================
                # RANDOM PERTURBATION
                # ====================================================
                if no_improve_counter >= 4:
    
                    perturb = np.random.uniform(
                        -0.02,
                        0.02
                    ) * (
                        self.alpha_max
                        - self.alpha_min
                    )
    
                    y_new += perturb
    
                # ====================================================
                # CLIP GLOBAL BOUNDS
                # ====================================================
                y_new = max(
                    self.alpha_min,
                    min(
                        self.alpha_max,
                        y_new
                    )
                )
    
                y_sol[(i,j)] = y_new
    
            # ========================================================
            # RANDOM SEGMENT JUMP
            # ========================================================
            if no_improve_counter >= 5:
    
                jump_arc = random.choice(
                    critical_arcs
                )
    
                s = seg_index[jump_arc]
    
                direction = random.choice(
                    [-1, 1]
                )
    
                s_new = max(
                    1,
                    min(
                        len(alpha_vals)-1,
                        s + direction
                    )
                )
    
                seg_index[jump_arc] = s_new
    
                y_sol[jump_arc] = 0.5 * (
                    alpha_vals[s_new-1]
                    +
                    alpha_vals[s_new]
                )
    
                print(
                    f"Jump arc {jump_arc} "
                    f"segment {s}->{s_new}"
                )
    
        # ============================================================
        # FINAL EXACT SOLVE
        # ============================================================
        print(
            "\n------------- FINAL EXACT SOLVE -------------"
        )
    
        ampl = AMPL()
    
        ampl.read("exact_reduced_wdn.mod")
    
        ampl.read_data(self.data_file)
    
        ampl.option["solver"] = "ipopt"
    
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "tol=1e-8 "
            "acceptable_tol=1e-6 "
            "max_iter=3000 "
            "mu_strategy=adaptive "
            "warm_start_init_point=yes "
            "bound_push=1e-4 "
            "bound_frac=1e-4 "
        )
    
        # ============================================================
        # WARM START BEST SOLUTION
        # ============================================================
        for (i,j), val in best_q.items():
            ampl.var["q"][i,j] = val
    
        for (i,j), val in best_y.items():
            ampl.var["y"][i,j] = val
    
        for i, val in best_h.items():
            ampl.var["h"][i] = val
    
        for (i,j), val in best_z.items():
            ampl.var["z"][i,j] = val
    
        with self.suppress_output():
            ampl.solve()
    
        q_sol = ampl.get_variable(
            "q"
        ).get_values().to_dict()
    
        y_sol = ampl.get_variable(
            "y"
        ).get_values().to_dict()
    
        h_sol = ampl.get_variable(
            "h"
        ).get_values().to_dict()
    
        z_sol = ampl.get_variable(
            "z"
        ).get_values().to_dict()
    
        final_cost = sum(
            z_sol[i,j]
            for (i,j) in self.arcs
        )
    
        print(
            "\nFinal feasible cost:",
            final_cost
        )
    
        return (
            q_sol,
            h_sol,
            y_sol,
            z_sol
        )

    def solve_restricted_model6(
        self,
        q_sol,
        z_sol,
        y_sol,
        h_sol,
        seg_index
    ):
    
        # ============================================================
        # INITIALIZATION
        # ============================================================
        ampl = AMPL()
    
        ampl.read("repair_feasibility_wdn.mod")
        ampl.read_data(self.data_file)
    
        ampl.option["solver"] = "ipopt"
    
        ampl.option["presolve_eps"] = 1e-7
    
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "max_iter=2000 "
            "warm_start_init_point=yes "
            "print_level=0 "
        )
    
        # ============================================================
        # PARAMETERS
        # ============================================================
        ampl.eval("""
            param seg_index{arcs};
            param alpha_arc{arcs};
            param y_ref{arcs};
            param w_head{arcs} >= 1 default 1;
        """)
    
        for (i,j) in self.arcs:
    
            ampl.param["seg_index"][i,j] = seg_index[(i,j)]
    
            ampl.param["alpha_arc"][i,j] = \
                self.alpha[seg_index[(i,j)]]
    
            ampl.param["y_ref"][i,j] = y_sol[(i,j)]
    
            ampl.param["w_head"][i,j] = 1.0
    
        # ============================================================
        # WARM START
        # ============================================================
        for (i,j), val in q_sol.items():
            ampl.var["q"][i,j] = val
    
        for (i,j), val in y_sol.items():
            ampl.var["y"][i,j] = val
    
        for i, val in h_sol.items():
            ampl.var["h"][i] = val
    
        for (i,j), val in z_sol.items():
            ampl.var["z"][i,j] = val
    
        # ============================================================
        # OBJECTIVE
        # ============================================================
        ampl.eval("""
    
            minimize restoration_obj:
    
                if mode = 1 then
    
                    sum{(i,j) in arcs}
                        w_head[i,j]
                        * delta_headloss[i,j]^2
    
                    + 1e-4 * sum{(i,j) in arcs}
                        delta_y[i,j]^2
    
                    - 1e-5 * sum{(i,j) in arcs}
                        (y[i,j] - y_ref[i,j])^2
    
                else
    
                    sum{(i,j) in arcs}
                        z[i,j]
    
                    + 1e-2 * sum{(i,j) in arcs}
                        delta_y[i,j]^2;
        """)
    
        # ============================================================
        # COST LIMIT
        # ============================================================
        ampl.eval(f"""
            subject to bound_cost:
    
                sum{{(i,j) in arcs}} z[i,j]
                <= {self.best_cost - 1e-5};
        """)
    
        # ============================================================
        # SEGMENT ESCAPE CUT
        # ============================================================
        alpha_vals = sorted(
            set(self.alpha.values()),
            reverse=True
        )
    
        factors = []
    
        eps_break = 1e-5
    
        for (i,j) in self.arcs:
    
            s = seg_index[(i,j)]
    
            a_up = alpha_vals[s-1]
            a_lo = alpha_vals[s]
    
            yv = y_sol[(i,j)]
    
            # ----------------------------------------
            # If on breakpoint
            # ----------------------------------------
            if (
                abs(yv - a_up) <= 1e-10
                or
                abs(yv - a_lo) <= 1e-10
            ):
    
                factors.append(
                    f"((y[{i},{j}] - ({a_lo} - {eps_break}))"
                    f"*(({a_lo} - {eps_break}) - y[{i},{j}]))"
                )
    
            # ----------------------------------------
            # Interior of segment
            # ----------------------------------------
            else:
    
                factors.append(
                    f"(({a_up} - y[{i},{j}])"
                    f"*(y[{i},{j}] - {a_lo}))"
                )
    
        product_expr = " * ".join(factors)
    
        ampl.eval(f"""
            subject to segment_escape:
                {product_expr} <= 0;
        """)
    
        # ============================================================
        # ITERATIVE RESTORATION
        # ============================================================
        max_iter = 200
    
        tol_head = 1e-6
        tol_y    = 1e-6
    
        lambda_k = 0.20
    
        prev_merit = 1e20
    
        stagnation_counter = 0
    
        print(
            "\niter  phase      "
            "head_inf      "
            "delta_y_inf      "
            "cost"
        )
    
        for it in range(1, max_iter + 1):
    
            # ========================================================
            # PHASE 1 : FEASIBILITY RESTORATION
            # ========================================================
            ampl.param["mode"] = 1
    
            for (i,j) in self.arcs:
    
                ampl.var["delta_headloss"][i,j] = 0.0
    
            with self.suppress_output():
                ampl.solve()
    
            solve_status = ampl.get_value("solve_result")
    
            if solve_status != "solved":
    
                print(f"\nRestoration NLP failed at iter {it}")
                break
    
            # ========================================================
            # EXTRACT SOLUTION
            # ========================================================
            q_sol = ampl.get_variable(
                "q"
            ).get_values().to_dict()
    
            y_sol = ampl.get_variable(
                "y"
            ).get_values().to_dict()
    
            h_sol = ampl.get_variable(
                "h"
            ).get_values().to_dict()
    
            z_sol = ampl.get_variable(
                "z"
            ).get_values().to_dict()
    
            delta_head = ampl.get_variable(
                "delta_headloss"
            ).get_values().to_dict()
    
            # ========================================================
            # INFEASIBILITY
            # ========================================================
            head_inf = max(
                abs(delta_head[i,j])
                for (i,j) in self.arcs
            )
    
            # ========================================================
            # UPDATE ACTIVE SEGMENTS
            # ========================================================
            for (i,j) in self.arcs:
    
                yv = y_sol[(i,j)]
    
                for k in range(len(alpha_vals)-1):
    
                    if (
                        alpha_vals[k+1]
                        <= yv
                        <= alpha_vals[k]
                    ):
    
                        seg_index[(i,j)] = k+1
                        break
    
            # ========================================================
            # COMPUTE COST
            # ========================================================
            cost = sum(
                self.L[i,j] * (
                    self.slope[seg_index[(i,j)]]
                    * y_sol[(i,j)]
                    +
                    self.intercept[seg_index[(i,j)]]
                )
                for (i,j) in self.arcs
            )
    
            print(
                f"{it:03d}   restore   "
                f"{head_inf:.3e}   "
                f"{0.0:.3e}   "
                f"{cost:.2f}"
            )
    
            # ========================================================
            # FEASIBLE
            # ========================================================
            if head_inf <= tol_head:
    
                print("\nRestoration feasible.")
                break
    
            # ========================================================
            # COMPUTE IMPROVED FEASIBILITY DIRECTION
            # ========================================================
            y_corr = {}
    
            for (i,j) in self.arcs:
    
                qv = q_sol[(i,j)]
    
                if abs(qv) <= 1e-8:
    
                    y_corr[(i,j)] = y_sol[(i,j)]
                    continue
    
                phi = (
                    qv**3
                    * (qv**2 + self.eps[i,j]**2)**0.426
                ) / (
                    qv**2
                    + 0.426*self.eps[i,j]**2
                )
    
                # ----------------------------------------
                # Newton-like correction
                # ----------------------------------------
                step = (
                    delta_head[(i,j)]
                    / (phi * self.L[i,j])
                )
    
                y_new = (
                    y_sol[(i,j)]
                    + lambda_k * step
                )
    
                # ----------------------------------------
                # Trust region
                # ----------------------------------------
                max_move = 0.03 * (
                    self.alpha_max
                    - self.alpha_min
                )
    
                y_new = max(
                    y_sol[(i,j)] - max_move,
                    min(
                        y_sol[(i,j)] + max_move,
                        y_new
                    )
                )
    
                # ----------------------------------------
                # Random perturbation
                # ----------------------------------------
                if stagnation_counter >= 4:
    
                    perturb = np.random.uniform(
                        -0.01,
                        0.01
                    ) * (
                        self.alpha_max
                        - self.alpha_min
                    )
    
                    y_new += perturb
    
                # ----------------------------------------
                # Global bounds
                # ----------------------------------------
                y_new = max(
                    self.alpha_min,
                    min(self.alpha_max, y_new)
                )
    
                y_corr[(i,j)] = y_new
    
            # ========================================================
            # UPDATE ACTIVE SEGMENTS
            # ========================================================
            for (i,j) in self.arcs:
    
                yv = y_corr[(i,j)]
    
                for k in range(len(alpha_vals)-1):
    
                    if (
                        alpha_vals[k+1]
                        <= yv
                        <= alpha_vals[k]
                    ):
    
                        seg_index[(i,j)] = k+1
                        break
    
            # ========================================================
            # UPDATE PARAMETERS
            # ========================================================
            for (i,j) in self.arcs:
    
                ampl.param["alpha_arc"][i,j] = \
                    y_corr[(i,j)]
    
                ampl.param["seg_index"][i,j] = \
                    seg_index[(i,j)]
    
            # ========================================================
            # PHASE 2 : EXACT PROJECTION
            # ========================================================
            ampl.param["mode"] = 2
    
            for (i,j) in self.arcs:
    
                ampl.var["delta_y"][i,j] = 0.0
    
            with self.suppress_output():
                ampl.solve()
    
            solve_status = ampl.get_value(
                "solve_result"
            )
    
            if solve_status != "solved":
    
                lambda_k = max(
                    0.05,
                    0.5 * lambda_k
                )
    
                stagnation_counter += 1
    
                continue
    
            # ========================================================
            # EXTRACT PROJECTED POINT
            # ========================================================
            q_sol = ampl.get_variable(
                "q"
            ).get_values().to_dict()
    
            y_sol = ampl.get_variable(
                "y"
            ).get_values().to_dict()
    
            h_sol = ampl.get_variable(
                "h"
            ).get_values().to_dict()
    
            z_sol = ampl.get_variable(
                "z"
            ).get_values().to_dict()
    
            delta_y = ampl.get_variable(
                "delta_y"
            ).get_values().to_dict()
    
            delta_y_inf = max(
                abs(delta_y[i,j])
                for (i,j) in self.arcs
            )
    
            # ========================================================
            # UPDATE SEGMENTS
            # ========================================================
            for (i,j) in self.arcs:
    
                yv = y_sol[(i,j)]
    
                for k in range(len(alpha_vals)-1):
    
                    if (
                        alpha_vals[k+1]
                        <= yv
                        <= alpha_vals[k]
                    ):
    
                        seg_index[(i,j)] = k+1
                        break
    
            # ========================================================
            # UPDATE PARAMETERS
            # ========================================================
            for (i,j) in self.arcs:
    
                ampl.param["alpha_arc"][i,j] = \
                    y_sol[(i,j)]
    
                ampl.param["seg_index"][i,j] = \
                    seg_index[(i,j)]
    
            # ========================================================
            # COMPUTE COST
            # ========================================================
            cost = sum(
                self.L[i,j] * (
                    self.slope[seg_index[(i,j)]]
                    * y_sol[(i,j)]
                    +
                    self.intercept[seg_index[(i,j)]]
                )
                for (i,j) in self.arcs
            )
    
            print(
                f"{it:03d}   project   "
                f"{0.0:.3e}   "
                f"{delta_y_inf:.3e}   "
                f"{cost:.2f}"
            )
    
            # ========================================================
            # MERIT FUNCTION
            # ========================================================
            merit = head_inf + delta_y_inf
    
            if merit < prev_merit:
    
                lambda_k = min(
                    0.50,
                    1.10 * lambda_k
                )
    
                stagnation_counter = 0
    
            else:
    
                lambda_k = max(
                    0.05,
                    0.50 * lambda_k
                )
    
                stagnation_counter += 1
    
            prev_merit = merit
    
            # ========================================================
            # TERMINATION
            # ========================================================
            if (
                head_inf <= tol_head
                and
                delta_y_inf <= tol_y
            ):
    
                print("\nFeasible point found.")
                break
    
        # ============================================================
        # FINAL EXACT NLP
        # ============================================================
        print(
            "\n---------------- FINAL EXACT NLP ----------------"
        )
    
        ampl = AMPL()
    
        ampl.read("exact_reduced_wdn.mod")
    
        ampl.read_data(self.data_file)
    
        ampl.option["solver"] = "ipopt"
    
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "tol=1e-8 "
            "max_iter=3000 "
            "mu_strategy=adaptive "
            "bound_relax_factor=0 "
            "warm_start_init_point=yes "
            "halt_on_ampl_error=yes "
        )
    
        # ============================================================
        # WARM START
        # ============================================================
        for (i,j), val in q_sol.items():
            ampl.var["q"][i,j] = val
    
        for (i,j), val in y_sol.items():
            ampl.var["y"][i,j] = val
    
        for i, val in h_sol.items():
            ampl.var["h"][i] = val
    
        for (i,j), val in z_sol.items():
            ampl.var["z"][i,j] = val
    
        with self.suppress_output():
            ampl.solve()
    
        # ============================================================
        # EXTRACT FINAL SOLUTION
        # ============================================================
        q_sol = ampl.get_variable(
            "q"
        ).get_values().to_dict()
    
        y_sol = ampl.get_variable(
            "y"
        ).get_values().to_dict()
    
        h_sol = ampl.get_variable(
            "h"
        ).get_values().to_dict()
    
        z_sol = ampl.get_variable(
            "z"
        ).get_values().to_dict()
    
        final_cost = sum(
            z_sol[i,j]
            for (i,j) in self.arcs
        )
    
        print(
            "\nFinal feasible cost:",
            final_cost
        )
    
        print(
            "Best cost:",
            self.best_cost
        )
    
        print(
            "\n------------- Exit Feasibility Restoration -------------\n"
        )
    
        return (
            q_sol,
            h_sol,
            y_sol,
            z_sol
        )

    def solve_restricted_model5(
        self,
        q_sol,
        z_sol,
        y_sol,
        h_sol,
        seg_index
    ):
        ampl = AMPL()
        ampl.read("repair_feasibility_wdn.mod")
        ampl.read_data(self.data_file)
        # =========================================================
        # IPOPT OPTIONS
        # =========================================================
        ampl.option["solver"] = "ipopt"
        ampl.option["presolve_eps"] = 9.43e-7
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "max_iter=1000 "
            "warm_start_init_point=yes "
        )
        # ampl.option["presolve_eps"] = 9.99e-05
        # =========================================================
        # PARAMETERS
        # =========================================================
        ampl.eval("param seg_index{arcs};")
        ampl.eval("param alpha_arc{arcs};")
        ampl.eval("param y_ref{arcs};")
        ampl.eval("param w_head{arcs} >= 1 default 1;")
        # ampl.eval("param mode integer default 1;")
        for (i,j) in self.arcs:
            ampl.param["seg_index"][i,j] = seg_index[i,j]
            ampl.param["alpha_arc"][i,j] = self.alpha[seg_index[i,j]]
            ampl.param["y_ref"][i,j] = y_sol[i,j]
            ampl.param["w_head"][i,j] = 1.0
        # =========================================================
        # WARM START
        # =========================================================
        for (i,j), val in q_sol.items():
            ampl.var["q"][i,j] = val
        for (i,j), val in y_sol.items():
            ampl.var["y"][i,j] = val
        for i, val in h_sol.items():
            ampl.var["h"][i] = val
        for (i,j), val in z_sol.items():
            ampl.var["z"][i,j] = val
        # =========================================================
        # COST LIMIT
        # =========================================================
        ampl.eval("""
            minimize restoration_obj:
                if mode = 1 then
                    sum{(i,j) in arcs}
                        w_head[i,j] * delta_headloss[i,j]^2
                    - 0 * sum{(i,j) in arcs}
                        (y[i,j] - y_ref[i,j])^2
                else
                    sum{(i,j) in arcs} (z[i,j] + delta_y[i,j]^2) 
                    - 1 * sum{(i,j) in arcs}
                        (y[i,j] - y_ref[i,j])^2;

            subject to bound_y{(i,j) in arcs}:
                y[i,j] >=
                    if mode = 1 then 
                        alpha_arc[i,j]
                    else 
                        alpha_arc[i,j] - delta_y[i,j];
        """)
        ampl.eval(f"""
            subject to bound_cost:
                sum{{(i,j) in arcs}} z[i,j]
                <= {self.best_cost} - 1e-5;
        """)
        # =========================================================
        # SEGMENT ESCAPE CUT
        # =========================================================
        # alpha_vals = sorted(
        #     set(self.alpha.values()),
        #     reverse=True
        # )
        # factors = []
        # for (i,j) in self.arcs:
        #     s = seg_index[(i,j)]
        #     a_up = alpha_vals[s-1]
        #     a_lo = alpha_vals[s]
        #     factors.append(
        #         f"(({a_up}-y[{i},{j}])"
        #         f"*(y[{i},{j}]-{a_lo}))"
        #     )
        # product_expr = " * ".join(factors)
        # ampl.eval(f"""
        #     subject to segment_escape:
        #         {product_expr} <= 0;
        # """)
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        # factors = []
        # for (i,j) in self.arcs:
        #     s = seg_index[(i,j)]
        #     a_up = alpha_vals[s-1]
        #     a_lo = alpha_vals[s]
        #     yv = y_sol[(i,j)]
        #     eps_break = 1e-4
        #     # if on breakpoint
        #     if abs(yv - a_up) <= 1e-8 or abs(yv - a_lo) <= 1e-8:
        #         factors.append(
        #             f"((y[{i},{j}] - ({a_lo} - {eps_break}))"
        #             f"*(({a_lo} - {eps_break}) - y[{i},{j}]))"
        #         )
        #     else:
        #         factors.append(
        #             f"(({a_up} - y[{i},{j}])"
        #             f"*(y[{i},{j}] - {a_lo}))"
        #         )
        # product_expr = " * ".join(factors)
        # ampl.eval(f"""
        #     subject to segment_escape:
        #         {product_expr} <= 0;
        # """)
        # =========================================================
        # FEASIBILITY RESTORATION LOOP
        # =========================================================
        max_iter = 400
        tol = 1e-6
        lambda_k = 0.25
        prev_inf = 1e20
        print(
            "\niter  phase     "
            "head_inf    "
            "delta_y_inf    "
            "cost"
        )
        for it in range(1, max_iter+1):
            # =====================================================
            # PHASE 1: RESTORATION
            # =====================================================
            ampl.param["mode"] = 1
            
            for (i,j) in self.arcs:
                ampl.var["delta_headloss"][i,j] = 0

            with self.suppress_output():
                ampl.solve()
            solve_status = ampl.get_value("solve_result")
            if solve_status != "solved":
                print(f"\nRestoration failed at iter {it}")
                break
            q_sol = ampl.get_variable(
                "q"
            ).get_values().to_dict()
            y_sol = ampl.get_variable(
                "y"
            ).get_values().to_dict()
            h_sol = ampl.get_variable(
                "h"
            ).get_values().to_dict()
            z_sol = ampl.get_variable(
                "z"
            ).get_values().to_dict()
            delta_head = ampl.get_variable(
                "delta_headloss"
            ).get_values().to_dict()
            # ampl.eval("display delta_headloss;")
            # =====================================================
            # INFEASIBILITY
            # =====================================================
            head_inf = max(
                abs(delta_head[i,j])
                for (i,j) in self.arcs
            )
            # =====================================================
            # UPDATE SEGMENTS
            # =====================================================
            for (i,j) in self.arcs:
               yv = y_sol[i,j]
               for k in range(len(alpha_vals)-1):
                   if (
                       alpha_vals[k+1]
                       <= yv
                       <= alpha_vals[k]
                   ):
                       seg_index[(i,j)] = k+1
                       break
            # =====================================================
            # COST
            # =====================================================
            cost = sum(
                self.L[i,j] * (
                    self.slope[seg_index[i,j]]
                    * y_sol[i,j]
                    + self.intercept[seg_index[i,j]]
                )
                for (i,j) in self.arcs
            )
            print(
                f"{it:03d}   restore   "
                f"{head_inf:.3e}   "
                f"{0.0:.3e}   "
                f"{cost:.2f}"
            )
            # =====================================================
            # FEASIBLE
            # =====================================================
            if head_inf <= tol:
                print("\nFeasible point restored.")
                break
            # =====================================================
            # FEASIBILITY DIRECTION
            # =====================================================
            y_corr = {}
            for (i,j) in self.arcs:
                q = q_sol[i,j]
                if abs(q) <= 1e-8:
                    y_corr[i,j] = y_sol[i,j]
                    continue
                phi = (
                    q**3
                    * (q**2 + self.eps[i,j]**2)**0.426
                ) / (
                    q**2 + 0.426*self.eps[i,j]**2
                )
                step = (
                    delta_head[i,j]
                    / (phi * self.L[i,j])
                )
                y_new = y_sol[i,j] + lambda_k * step
                # TRUST REGION
                max_step = 0.03 * (
                    self.alpha_max - self.alpha_min
                )
                y_new = max(
                    y_sol[i,j] - max_step,
                    min(y_sol[i,j] + max_step, y_new)
                )
                # GLOBAL BOUNDS
                y_new = max(
                    self.alpha_min,
                    min(self.alpha_max, y_new)
                )
                y_corr[i,j] = y_new
            # =====================================================
            # UPDATE SEGMENTS
            # =====================================================
            for (i,j) in self.arcs:
                yv = y_corr[i,j]
                for k in range(len(alpha_vals)-1):
                    if (
                        alpha_vals[k+1]
                        <= yv
                        <= alpha_vals[k]
                    ):
                        seg_index[(i,j)] = k+1
                        break
            # =====================================================
            # UPDATE PARAMETERS
            # =====================================================
            for (i,j) in self.arcs:
                ampl.param["alpha_arc"][i,j] = y_corr[i,j]
                ampl.param["seg_index"][i,j] = seg_index[(i,j)]
            # =====================================================
            # PHASE 2: EXACT PROJECTION
            # =====================================================
            ampl.param["mode"] = 2
            
            for (i,j) in self.arcs:
                ampl.var["delta_y"][i,j] = 0 
            with self.suppress_output():
                ampl.solve()
            solve_status = ampl.get_value("solve_result")
            if solve_status != "solved":
                lambda_k *= 0.5
                # continue
            # =====================================================
            # ADAPTIVE STEP
            # =====================================================
            if head_inf < prev_inf:
                lambda_k = min(0.5, 1.1*lambda_k)
            else:
                lambda_k = max(0.05, 0.5*lambda_k)
            prev_inf = head_inf

            q_sol = ampl.get_variable(
                "q"
            ).get_values().to_dict()
            y_sol = ampl.get_variable(
                "y"
            ).get_values().to_dict()
            h_sol = ampl.get_variable(
                "h"
            ).get_values().to_dict()
            z_sol = ampl.get_variable(
                "z"
            ).get_values().to_dict()

            delta_y = ampl.get_variable(
                "delta_y"
            ).get_values().to_dict()
            # ampl.eval("display delta_y;")
            # =====================================================
            # INFEASIBILITY
            # =====================================================
            delta_y_inf = max(
                abs(delta_y[i,j])
                for (i,j) in self.arcs
            ) 
            # =====================================================
            # FEASIBLE
            # =====================================================
            if delta_y_inf <= tol:
                print("\nFeasible point restored.")
                break
            for (i,j) in self.arcs:
                yv = y_sol[i,j]
                for k in range(len(alpha_vals)-1):
                    if (
                        alpha_vals[k+1]
                        <= yv
                        <= alpha_vals[k]
                    ):
                        seg_index[(i,j)] = k+1
                        break
            # =====================================================
            # UPDATE PARAMETERS
            # =====================================================
            for (i,j) in self.arcs:
                ampl.param["alpha_arc"][i,j] = y_sol[i,j]
                ampl.param["seg_index"][i,j] = seg_index[(i,j)]
            print(
                f"{it:03d}   project   "
                f"{0.0:.3e}   "
                f"{delta_y_inf:.3e}   "
                f"{cost:.2f}"
            )
        # =========================================================
        # FINAL EXACT NLP
        # =========================================================
        # ampl.param["mode"] = 2
        # ampl.solve()
        print("------------------------------FINAL EXACT NLP---------------------------------")
        ampl = AMPL()
        ampl.read("exact_reduced_wdn.mod")
        ampl.read_data(self.data_file)

        ampl.option["solver"] = "ipopt"
        # ampl.option["ipopt_options"] = "outlev=0 tol=1e-9 max_iter=3000"
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "expect_infeasible_problem=no "
            "bound_relax_factor=0 "
            "bound_push=0.1 "
            "bound_frac=0.1 "
            "warm_start_init_point=yes "
            "halt_on_ampl_error=yes "
            "max_iter=3000"
        )

        # -------------------------
        # Storage
        # -------------------------
        # Warm-start perturbation
        # -------------------------
        for (i, j), val in q_sol.items():
            ampl.var["q"][i, j] = val
        for (i, j), val in y_sol.items():
            ampl.var["y"][i, j] = val
        for i, val in h_sol.items():
            ampl.var["h"][i] = val
        for (i, j), val in z_sol.items():
            ampl.var["z"][i, j] = val

        ampl.solve()
        # ampl.eval("display total_cost;")
 
        q_sol = ampl.get_variable(
            "q"
        ).get_values().to_dict()
        y_sol = ampl.get_variable(
            "y"
        ).get_values().to_dict()
        h_sol = ampl.get_variable(
            "h"
        ).get_values().to_dict()
        z_sol = ampl.get_variable(
            "z"
        ).get_values().to_dict()

        final_cost = sum(
            z_sol[i,j]
            for (i,j) in self.arcs
        )
        print(
            "\nFinal feasible cost:",
            final_cost
        )
        print("Best cost:", self.best_cost)
        print("\n-----------------Exit the Feasibility Restoration Phase------------\n")
        return (
            q_sol,
            h_sol,
            y_sol,
            z_sol
        )

    def solve_restricted_model5(self, q_sol, z_sol, y_sol, h_sol, seg_index):
    
        ampl = AMPL()
        ampl.read("repair_feasibility_wdn.mod")
        ampl.read_data(self.data_file)
    
        ampl.option["solver"] = "ipopt"
        ampl.option["presolve_eps"] = 1e-8
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "max_iter=2000 "
            # "tol=1e-7 "
            "warm_start_init_point=yes "
            # "mu_strategy=adaptive "
        )
    
        # ============================================================
        # PARAMETERS
        # ============================================================
        ampl.eval("param seg_index{arcs};")
        ampl.eval("param alpha_arc{arcs};")
        ampl.eval("param y_ref{arcs};")
        ampl.eval("param w_head{arcs} >= 1 default 1;")
    
        for (i,j) in self.arcs:
            ampl.param["seg_index"][i,j] = seg_index[i,j]
            ampl.param["alpha_arc"][i,j] = self.alpha[seg_index[i,j]]
            ampl.param["y_ref"][i,j] = y_sol[i,j]
            ampl.param["w_head"][i,j] = 1.0
    
        # ============================================================
        # WARM START
        # ============================================================
        for (i,j), val in q_sol.items():
            ampl.var["q"][i,j] = val
    
        for (i,j), val in y_sol.items():
            ampl.var["y"][i,j] = val
    
        for i,val in h_sol.items():
            ampl.var["h"][i] = val
    
        for (i,j), val in z_sol.items():
            ampl.var["z"][i,j] = val
    
        # ============================================================
        # FEASIBILITY RESTORATION MODEL
        # ============================================================
        ampl.eval("""
    
            # minimize restoration_obj:
            #
            #     sum{(i,j) in arcs}
            #         w_head[i,j] * delta_headloss[i,j]^2
            #
            #     + 1e-2 * sum{(i,j) in arcs}
            #         delta_y[i,j]^2
            #
            #     - 1e-3 * sum{(i,j) in arcs}
            #         (y[i,j] - y_ref[i,j])^2;
            param mode integer default 1;

            minimize restoration_obj:

            if mode = 1 then

                sum{(i,j) in arcs}
                    delta_headloss[i,j]^2

                + 1e-4 * sum{(i,j) in arcs}
                    delta_y[i,j]^2

                - 1e-5 * sum{(i,j) in arcs}
                    (y[i,j] - alpha_arc[i,j])^2

            else

                sum{(i,j) in arcs}
                    z[i,j];
        """)
    
        # ampl.eval("""
        #
        #     subject to bound_y{(i,j) in arcs}:
        #         y[i,j] >= alpha_arc[i,j] - delta_y[i,j];
        #
        # """)
    
        ampl.eval(f"""
    
            subject to bound_cost:
                sum{{(i,j) in arcs}} z[i,j]
                <= {self.best_cost};
    
        """)
    
        # ============================================================
        # SEGMENT ESCAPE CUT
        # ============================================================
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
    
        factors = []
    
        for (i,j) in self.arcs:
    
            s = seg_index[(i,j)]
    
            a_up = alpha_vals[s-1]
            a_lo = alpha_vals[s]
    
            yv = y_sol[(i,j)]
    
            eps_break = 1e-4
    
            # if on breakpoint
            if abs(yv - a_up) <= 1e-8 or abs(yv - a_lo) <= 1e-8:
    
                factors.append(
                    f"((y[{i},{j}] - ({a_lo} - {eps_break}))"
                    f"*(({a_lo} - {eps_break}) - y[{i},{j}]))"
                )
    
            else:
    
                factors.append(
                    f"(({a_up} - y[{i},{j}])"
                    f"*(y[{i},{j}] - {a_lo}))"
                )
    
        product_expr = " * ".join(factors)
    
        ampl.eval(f"""
            subject to segment_escape:
                {product_expr} <= 0;
        """)
    
        # ============================================================
        # ITERATIVE FEASIBILITY JUMP
        # ============================================================
        max_iter = 100
        tol = 1e-6
        
        lambda_k = 0.3
        
        print("\niter    phase      head_inf       cost")
        
        for it in range(1, max_iter+1):
        
            # =====================================================
            # PHASE 1: RESTORATION NLP
            # =====================================================
        
            ampl.param["mode"] = 1
            ampl.eval(f"""
                subject to con2_{it}{{(i,j) in arcs}}:
                    h[i] - h[j] - ( q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426 * eps[i,j]^2)) * y[i,j] * L[i,j] = delta_headloss[i,j];

                subject to bound_y_{it}{{(i,j) in arcs}}:
                    y[i,j] >= alpha_arc[i,j];
            """)
        
            with self.suppress_output():
                ampl.solve()
        
            q_sol = ampl.get_variable("q").get_values().to_dict()
            y_sol = ampl.get_variable("y").get_values().to_dict()
            h_sol = ampl.get_variable("h").get_values().to_dict()
            z_sol = ampl.get_variable("z").get_values().to_dict()
        
            delta_head = ampl.get_variable(
                "delta_headloss"
            ).get_values().to_dict()
        
            head_inf = max(
                abs(delta_head[i,j])
                for (i,j) in self.arcs
            )
         
            alpha_vals = sorted(
                set(self.alpha.values()),
                reverse=True
            )
        
            for (i,j) in self.arcs:
        
                yv = y_sol[i,j]
        
                for k in range(len(alpha_vals)-1):
        
                    if alpha_vals[k+1] <= yv <= alpha_vals[k]:
        
                        seg_index[(i,j)] = k+1
                        break

            cost = sum(
                self.L[i,j] * (
                    self.slope[seg_index[i,j]] * y_sol[i,j]
                    + self.intercept[seg_index[i,j]]
                )
                for (i,j) in self.arcs
            )
        
            print(
                f"{it:03d}    restore    "
                f"{head_inf:.3e}    "
                f"{cost:.2f}"
            )
        
            # =====================================================
            # CONVERGED
            # =====================================================
        
            if head_inf <= tol:
        
                print("\nFeasible solution restored.")
                break
        
            # =====================================================
            # COMPUTE FEASIBILITY DIRECTION
            # =====================================================
        
            y_corr = {}
        
            for (i,j) in self.arcs:
        
                q = q_sol[i,j]
        
                if abs(q) <= 1e-8:
        
                    y_corr[i,j] = y_sol[i,j]
                    continue
        
                phi = (
                    q**3
                    * (q**2 + self.eps[i,j]**2)**0.426
                ) / (
                    q**2 + 0.426*self.eps[i,j]**2
                )
        
                step = (
                    delta_head[i,j]
                    / (phi * self.L[i,j])
                )
        
                y_new = y_sol[i,j] + lambda_k * step
        
                # trust region
                max_move = 0.05 * (
                    self.alpha_max - self.alpha_min
                )
        
                y_new = max(
                    y_sol[i,j] - max_move,
                    min(y_sol[i,j] + max_move, y_new)
                )
        
                y_new = max(
                    self.alpha_min,
                    min(self.alpha_max, y_new)
                )
        
                y_corr[i,j] = y_new
        
            # =====================================================
            # UPDATE SEGMENTS
            # =====================================================
        
            alpha_vals = sorted(
                set(self.alpha.values()),
                reverse=True
            )
        
            for (i,j) in self.arcs:
        
                yv = y_corr[i,j]
        
                for k in range(len(alpha_vals)-1):
        
                    if alpha_vals[k+1] <= yv <= alpha_vals[k]:
        
                        seg_index[(i,j)] = k+1
                        break
        
            # =====================================================
            # UPDATE PARAMETERS
            # =====================================================
        
            for (i,j) in self.arcs:
                ampl.param["alpha_arc"][i,j] = y_corr[i,j]
        
            # =====================================================
            # PHASE 2: EXACT PROJECTION NLP
            # =====================================================
        
            ampl.param["mode"] = 2
            ampl.eval(f"drop bound_y_{it};")
            ampl.eval(f"drop con2_{it};")
            ampl.eval(f"""
                subject to con2_new_{it}{{(i,j) in arcs}}:
                    h[i] - h[j] - ( q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426 * eps[i,j]^2)) * y[i,j] * L[i,j] = 0;

                subject to bound_y_new_{it}{{(i,j) in arcs}}:
                    y[i,j] >= alpha_arc[i,j] - delta_y[i,j];
            """) 

            with self.suppress_output():
                ampl.solve()
             
            ampl.eval(f"drop bound_y_new_{it};")
            ampl.eval(f"drop con2_new_{it};")
            # ampl.eval("""
            #     subject to con2{(i,j) in arcs}:
            #         h[i] - h[j] - ( q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426 * eps[i,j]^2)) * y[i,j] * L[i,j] = delta_headloss[i,j];
            #
            #     subject to bound_y_new{(i,j) in arcs}:
            #         y[i,j] >= alpha_arc[i,j];
            # """)
            q_sol = ampl.get_variable("q").get_values().to_dict()
            y_sol = ampl.get_variable("y").get_values().to_dict()
            h_sol = ampl.get_variable("h").get_values().to_dict()
            z_sol = ampl.get_variable("z").get_values().to_dict()
        
            cost = sum(
                self.L[i,j] * (
                    self.slope[seg_index[i,j]] * y_sol[i,j]
                    + self.intercept[seg_index[i,j]]
                )
                for (i,j) in self.arcs
            )
        
            print(
                f"{it:03d}    project    "
                f"{0.0:.3e}    "
                f"{cost:.2f}"
            )

        # max_iter = 500
        # tol = 1e-5
        #
        # stagnation_counter = 0
        # prev_inf = 1e20
        #
        # lambda_k = 0.2
        #
        # print("\niter   head_inf     dy_inf      total_inf      cost")
        #
        # for it in range(1, max_iter+1):
        #
        #     # ========================================================
        #     # SOLVE RESTORATION NLP
        #     # ========================================================
        #     with self.suppress_output():
        #         ampl.solve()
        #
        #     # ========================================================
        #     # EXTRACT SOLUTION
        #     # ========================================================
        #     q_sol = ampl.get_variable("q").get_values().to_dict()
        #     y_sol = ampl.get_variable("y").get_values().to_dict()
        #     h_sol = ampl.get_variable("h").get_values().to_dict()
        #     z_sol = ampl.get_variable("z").get_values().to_dict()
        #
        #     delta_head = ampl.get_variable(
        #         "delta_headloss"
        #     ).get_values().to_dict()
        #
        #     delta_y = ampl.get_variable(
        #         "delta_y"
        #     ).get_values().to_dict()
        #
        #     # ========================================================
        #     # INFEASIBILITY
        #     # ========================================================
        #     head_inf = max(
        #         abs(delta_head[i,j]) for (i,j) in self.arcs
        #     )
        #
        #     dy_inf = max(
        #         abs(delta_y[i,j]) for (i,j) in self.arcs
        #     )
        #
        #     total_inf = head_inf + dy_inf
        #
        #     # ========================================================
        #     # COMPUTE COST
        #     # ========================================================
        #     cost = sum(
        #         self.L[i,j] * (
        #             self.slope[seg_index[i,j]] * y_sol[i,j]
        #             + self.intercept[seg_index[i,j]]
        #         )
        #         for (i,j) in self.arcs
        #     )
        #
        #     print(
        #         f"{it:02d}   "
        #         f"{head_inf:.2e}   "
        #         f"{dy_inf:.2e}   "
        #         f"{total_inf:.2e}   "
        #         f"{cost:.2f}"
        #     )
        #
        #     # ========================================================
        #     # CONVERGENCE
        #     # ========================================================
        #     if total_inf <= tol:
        #
        #         print("\nFeasible point restored.")
        #         break
        #
        #     # ========================================================
        #     # STAGNATION DETECTION
        #     # ========================================================
        #     rel_improve = abs(prev_inf - total_inf) / max(1.0, prev_inf)
        #
        #     if rel_improve < 1e-3:
        #         stagnation_counter += 1
        #     else:
        #         stagnation_counter = 0
        #
        #     prev_inf = total_inf
        #
        #     # ========================================================
        #     # ADAPTIVE WEIGHTS
        #     # ========================================================
        #     for (i,j) in self.arcs:
        #
        #         w = ampl.param["w_head"][i,j]
        #
        #         if abs(delta_head[i,j]) > 1e-3:
        #             ampl.param["w_head"][i,j] = min(1e6, 1.5*w)
        #
        #     # ========================================================
        #     # COMPUTE FEASIBILITY DIRECTION
        #     # ========================================================
        #     y_corr = {}
        #
        #     for (i,j) in self.arcs:
        #
        #         q = q_sol[i,j]
        #
        #         if abs(q) <= 1e-8:
        #
        #             y_corr[i,j] = y_sol[i,j]
        #             continue
        #
        #         phi = (
        #             q**3
        #             * (q**2 + self.eps[i,j]**2)**0.426
        #         ) / (
        #             q**2 + 0.426*self.eps[i,j]**2
        #         )
        #
        #         correction = (
        #             delta_head[i,j]
        #             / (phi * self.L[i,j])
        #         )
        #
        #         # damped correction
        #         y_new = (
        #             (1-lambda_k)*y_sol[i,j]
        #             + lambda_k*(y_sol[i,j] + correction)
        #         )
        #
        #         # random perturbation if stagnating
        #         if stagnation_counter >= 3:
        #
        #             perturb = np.random.uniform(
        #                 -0.02,
        #                 0.02
        #             ) * (self.alpha_max - self.alpha_min)
        #
        #             y_new += perturb
        #
        #         y_new = max(
        #             self.alpha_min,
        #             min(self.alpha_max, y_new)
        #         )
        #
        #         y_corr[i,j] = y_new
        #
        #     # ========================================================
        #     # UPDATE SEGMENTS
        #     # ========================================================
        #     for (i,j) in self.arcs:
        #
        #         yv = y_corr[i,j]
        #
        #         for k in range(len(alpha_vals)-1):
        #
        #             if alpha_vals[k+1] <= yv <= alpha_vals[k]:
        #
        #                 seg_index[(i,j)] = k+1
        #                 break
        #
        #     # ========================================================
        #     # UPDATE PARAMETERS
        #     # ========================================================
        #     for (i,j) in self.arcs:
        #
        #         ampl.param["alpha_arc"][i,j] = y_corr[i,j]
        #         ampl.param["seg_index"][i,j] = seg_index[(i,j)]
        #         # ampl.param["y_ref"][i,j] = y_sol[i,j]
        #
        #     # ========================================================
        #     # ADAPTIVE STEP SIZE
        #     # ========================================================
        #     if stagnation_counter >= 3:
        #         lambda_k = max(0.05, 0.5*lambda_k)
        #     else:
        #         lambda_k = min(0.5, 1.1*lambda_k)
    
        # ============================================================
        # FINAL EXACT NLP
        # ============================================================
        print("\n--- Exact Feasible Solve ---")
    
        ampl.eval("drop con2;")
        ampl.eval("drop restoration_obj;")
    
        ampl.eval("""
    
            minimize total_cost:
                sum{(i,j) in arcs} z[i,j];
    
            subject to headloss_exact{(i,j) in arcs}:
    
                h[i] - h[j]
    
                - (
    
                    q[i,j]^3
                    * (q[i,j]^2 + eps[i,j]^2)^0.426
    
                    / (q[i,j]^2 + 0.426*eps[i,j]^2)
    
                ) * y[i,j] * L[i,j]
    
                = 0;
    
        """)
    
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "warm_start_init_point=yes "
        )
    
        ampl.solve()
    
        q_sol = ampl.get_variable("q").get_values().to_dict()
        y_sol = ampl.get_variable("y").get_values().to_dict()
        h_sol = ampl.get_variable("h").get_values().to_dict()
        z_sol = ampl.get_variable("z").get_values().to_dict()
    
        final_cost = ampl.get_objective(
            "total_cost"
        ).value()
    
        print("\nFinal feasible cost:", final_cost)
    
        return q_sol, h_sol, y_sol, z_sol

    def solve_restricted_model3(self, q_sol, z_sol, y_sol, h_sol, seg_index):

        ampl = AMPL()
        ampl.read("repair_feasibility_wdn.mod")
        ampl.read_data(self.data_file)

        ampl.option["solver"] = "ipopt"
        ampl.option["ipopt_options"] = (
            "outlev=0 max_iter=3000 warm_start_init_point=no"
        )

        # -------------------------
        # Parameters
        # -------------------------
        ampl.eval("param seg_index{(i,j) in arcs};")
        ampl.eval("param optimal_seg_index{(i,j) in arcs};")
        ampl.eval("param alpha_arc{arcs};")
        ampl.eval("param y_sol{arcs};")
        
        for (i,j) in self.arcs:
            ampl.param["seg_index"][i, j] = seg_index[(i, j)]
            ampl.param["optimal_seg_index"][i, j] = seg_index[(i, j)]
            ampl.param["y_sol"][i, j] = y_sol[i,j]
        
        ampl.eval("""set fixed_arcs within {i in nodes, j in nodes: i != j};""")
        ampl.get_set("fixed_arcs").set_values(self.fixed_arcs)

        # -------------------------
        # Warm start
        # -------------------------
        # for (i, j), val in q_sol.items():
        #     ampl.var["q"][i, j] = val
        # for (i, j), val in y_sol.items():
        #     ampl.var["y"][i, j] = val
        # for i, val in h_sol.items():
        #     ampl.var["h"][i] = val
        # for (i, j), val in z_sol.items():
        #     ampl.var["z"][i, j] = val

        # -------------------------
        # Relaxed constraints
        # -------------------------
        ampl.eval(f"""
            minimize infeas_obj:
                    sum{{(i,j) in arcs}} (delta_headloss[i,j]^2 + delta_y[i,j]^2 - 0.0*(y[i,j] - y_sol[i,j])^2 + 0.0*z[i,j]);

            subject to bound_y{{(i,j) in arcs diff fixed_arcs}}:
                y[i,j] >= alpha_arc[i,j] - delta_y[i,j];
            
            subject to fix_y{{(i,j) in fixed_arcs}}:
                y[i,j] >= y_sol[i,j];

            subject to bound_delta_y{{(i,j) in arcs diff fixed_arcs}}:
                delta_y[i,j] <= max(0, alpha[seg_index[i,j]] - alpha_min);
        """)

        # ampl.eval(f"""
        #     subject to bound_cost:
        #         sum{{(i,j) in arcs}} z[i,j] <= {self.best_cost};
        # """)

        d_min = float('inf')
        for (i,j) in self.arcs:
            s      = seg_index[(i,j)]        # active segment index
            a_up   = self.alpha[s]                # alpha_s   (upper, larger value)
            a_lo   = self.alpha[s+1]              # alpha_{s+1} (lower, smaller value)
            y_val  = y_sol[(i,j)]
            # if y_sol[i,j]>a_lo:
            dist_to_upper = (a_up - y_val)     # >= 0
            dist_to_lower = (y_val - a_lo)     # >= 0
            d_ij   = min(dist_to_upper, dist_to_lower)
            if d_ij >= 0.001:
                d_min  = min(d_min, d_ij)
                # print("d_min:",d_min)

        rhs = d_min ** 2
        # ampl.eval(f"""
        #     subject to pwl_cut:
        #         sum{{(i,j) in arcs}} (y[i,j] - y_sol[i,j])^2 >= {rhs};
        # """)

        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        factors = []
        delta = 1e-4
        for (i,j) in self.arcs:
            s    = seg_index[(i,j)]
            a_up = alpha_vals[s-1]
            a_lo = alpha_vals[s]

            # factor > 0 iff y[i,j] is strictly inside segment s
            # factor <= 0 iff y[i,j] has left segment s
            if y_sol[i,j]>=self.alpha[s+1]+1e-4:
                factors.append(
                    f"(({a_up} - y[{i},{j}]) * "
                    f"(y[{i},{j}] - {a_lo}))"
                )
            else:
                factors.append(
                    f"(({a_up - delta} - y[{i},{j}]) * "
                    f"(y[{i},{j}] - {a_up + delta}))"
                )

        product = " * ".join(factors)
        ampl.eval(f"""
            subject to segment_escape_{1}:
                {product} <= 0;
        """)
        print(f"  Segment escape constraint {1} added")

        # -------------------------
        # Iterative restoration
        # -------------------------
        max_iter = 30
        tol = 1e-4

        lambda_k = 0.5        # damping
        lambda_min = 1e-3
        lambda_max = 1.0

        prev_inf = 1e20
        log_fmt = "[{:02d}] {:.1e} {:.1e} {:.1e} {:.1e} {}"
        print("iter head_inf delta_y_inf total_inf lambda cost")

        for it in range(1, max_iter+1):
            if it == 1:
                for (i, j) in self.arcs:
                    if (i,j) not in self.fixed_arcs:
                        ampl.param["alpha_arc"][i, j] = self.alpha[seg_index[i,j]]
                    else:
                        ampl.param["alpha_arc"][i, j] = y_sol[i,j]

            else:
                for (i, j) in self.arcs:
                    ampl.param["alpha_arc"][i, j] = y_corr[i,j]
                    ampl.param["seg_index"][i, j] = seg_index[(i, j)]
                    # ampl.param["y_sol"][i, j] = y_sol[i,j]
                    # ampl.eval("drop bound_delta_y;")
                    # ampl.eval("""
                    #     subject to bound_y_new{(i,j) in arcs}:
                    #         y[i,j] >= alpha_arc[i,j] - delta_y[i,j];
                    # """)

            with self.suppress_output():
                ampl.solve()

            q_sol = ampl.get_variable("q").get_values().to_dict()
            y_sol = ampl.get_variable("y").get_values().to_dict()
            h_sol = ampl.get_variable("h").get_values().to_dict()
            z_sol = ampl.get_variable("z").get_values().to_dict()

            delta_head = ampl.get_variable("delta_headloss").get_values().to_dict()
            delta_y_sol = ampl.get_variable("delta_y").get_values().to_dict()
            # delta_y_minus = ampl.get_variable("delta_y_minus").get_values().to_dict()

            # for (i,j) in self.arcs: 
            #     print(f"delta_y[{i},{j}]:", delta_y_sol[i,j])
            #     print(f"delta_headloss[{i},{j}]:", delta_head[i,j])

            # ampl.eval("display y;")
            # ampl.eval("display seg_index;")
            # ampl.eval("display delta_y;")
            # ampl.eval("display delta_headloss;")
            #     print(f"delta_y_minus[{i},{j}]:", delta_y_minus[i,j])
            # -------------------------
            # Infeasibility metrics
            # -------------------------
            head_inf = max(abs(delta_head[i,j]) for (i,j) in delta_head)
            dy_inf   = max(abs(delta_y_sol[i,j]) for (i,j) in delta_y_sol)

            total_inf = head_inf + dy_inf

            alpha_vals = sorted(set(self.alpha.values()), reverse=True)

            for (i,j) in self.arcs:
                yv = y_sol[i,j]
                for k in range(len(alpha_vals)-1):
                    if alpha_vals[k+1] <= yv <= alpha_vals[k]:
                        seg_index[(i,j)] = k+1
                        break
            cost = sum(
               self.L[i,j] * (
                   self.slope[seg_index[i,j]] * y_sol[i,j]
                   + self.intercept[seg_index[i,j]]
               )
               for (i,j) in self.arcs
            )

            print(log_fmt.format(it, head_inf, dy_inf, total_inf, lambda_k, self.format_indian_number(round(cost))))
            # ampl.eval("display delta_headloss;")
            # -------------------------
            # Convergence check
            # -------------------------
            if total_inf < tol:
                print("Converged")
                break

            # -------------------------
            # Compute correction
            # -------------------------
            y_corr = {}

            for (i,j) in self.arcs:

                q = q_sol[i,j]

                # avoid division issues
                if abs(q) < 1e-8:
                    y_corr[i,j] = y_sol[i,j]
                    continue

                phi = (q**3 * (q**2 + self.eps[i,j]**2)**0.426) / \
                      (q**2 + 0.426 * self.eps[i,j]**2)

                dy = delta_head[i,j] / (phi * self.L[i,j])

                # damped update
                # y_new = y_sol[i,j] + lambda_k * dy
                y_new = y_sol[i,j] + dy

                # trust region (important)
                # max_step = 0.2 * (self.alpha_max - self.alpha_min)
                # y_new = max(y_sol[i,j] - max_step,
                #             min(y_sol[i,j] + max_step, y_new))

                # clip bounds
                y_corr[i,j] = min(self.alpha_max,
                                 max(self.alpha_min, y_new))

            # -------------------------
            # Accept / reject step
            # -------------------------
            if total_inf < prev_inf:
                lambda_k = min(lambda_max, 1.2 * lambda_k)
            else:
                lambda_k = max(lambda_min, 0.5 * lambda_k)

            prev_inf = total_inf

            # -------------------------
            # Update AMPL parameters
            # -------------------------
            # for (i,j) in self.arcs:
            #     ampl.param["alpha_arc"][i,j] = y_corr[i,j]

            # -------------------------
            # Update segment index
            # -------------------------
            alpha_vals = sorted(set(self.alpha.values()), reverse=True)

            for (i,j) in self.arcs:
                yv = y_corr[i,j]
                for k in range(len(alpha_vals)-1):
                    if alpha_vals[k+1] <= yv <= alpha_vals[k]:
                        seg_index[(i,j)] = k+1
                        break

        # -------------------------
        # Final exact solve
        # -------------------------
        print("\n--- Final exact solve ---")

        # ampl.eval("drop con2;")
        # ampl.eval("drop infeas_obj;")

        # ampl.eval("""
        #     # minimize total_cost: sum{(i,j) in arcs} z[i,j];
        #
        #     subject to headloss_exact{(i,j) in arcs}:
        #         h[i] - h[j]
        #       - ( q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426
        #         / (q[i,j]^2 + 0.426 * eps[i,j]^2))
        #         * y[i,j] * L[i,j] = 0;
        # """)

        # ampl.solve()
        #
        # q_sol = ampl.get_variable("q").get_values().to_dict()
        # h_sol = ampl.get_variable("h").get_values().to_dict()
        # y_sol = ampl.get_variable("y").get_values().to_dict()
        # z_sol = ampl.get_variable("z").get_values().to_dict()
        #
        # alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        #
        # for (i,j) in self.arcs:
        #     yv = y_sol[i,j]
        #     for k in range(len(alpha_vals)-1):
        #         if alpha_vals[k+1] <= yv <= alpha_vals[k]:
        #             seg_index[(i,j)] = k+1
        #             break
        #
        # cost = sum(
        #        self.L[i,j] * (
        #            self.slope[seg_index[i,j]] * y_sol[i,j]
        #            + self.intercept[seg_index[i,j]]
        #        )
        #        for (i,j) in self.arcs
        #     ) 
        #
        # print("Feasible solution obtained.")
        print("Final cost:", cost) 

        ampl = AMPL()
        ampl.read("exact_reduced_wdn.mod")
        ampl.read_data(self.data_file)

        ampl.option["solver"] = "ipopt"
        # ampl.option["ipopt_options"] = "outlev=0 tol=1e-9 max_iter=3000"
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "expect_infeasible_problem=no "
            "bound_relax_factor=0 "
            "bound_push=0.1 "
            "bound_frac=0.1 "
            "warm_start_init_point=yes "
            "halt_on_ampl_error=yes "
            "max_iter=3000"
        )

        # -------------------------
        # Storage
        # -------------------------
        # Warm-start perturbation
        # -------------------------
        for (i, j), val in q_sol.items():
            ampl.var["q"][i, j] = val
        for (i, j), val in y_sol.items():
            ampl.var["y"][i, j] = val
        for i, val in h_sol.items():
            ampl.var["h"][i] = val
        for (i, j), val in z_sol.items():
            ampl.var["z"][i, j] = val

        ampl.solve()
        print("Best cost:", self.best_cost)
        ampl.eval("display total_cost;")
 
        print("\n-----------------Exit the Feasibility Restoration Phase------------\n")
        return q_sol, h_sol, y_sol, z_sol


    def solve_restricted_model2(self, q_sol, z_sol, y_sol, h_sol, seg_index):
    
        ampl = AMPL()
        ampl.read("repair_feasibility_wdn.mod")
        ampl.read_data(self.data_file)
    
        ampl.option["solver"] = "ipopt"
        ampl.option["presolve_eps"] = 1.2e-08
        ampl.option["ipopt_options"] = (
            "outlev=0 max_iter=3000 "
            "warm_start_init_point=yes"
        )

        # -------------------------
        # Parameters
        # -------------------------
        ampl.eval("param seg_index{(i,j) in arcs};")

        for (i, j) in self.arcs:
            ampl.param["seg_index"][i, j] = seg_index[(i, j)]

        # -------------------------
        # Warm start
        # -------------------------
        for (i, j), val in q_sol.items():
            ampl.var["q"][i, j] = val
        for (i, j), val in y_sol.items():
            ampl.var["y"][i, j] = val
        for i, val in h_sol.items():
            ampl.var["h"][i] = val
        for (i, j), val in z_sol.items():
            ampl.var["z"][i, j] = val

        # -------------------------
        # PHASE 1: relaxed structure
        # -------------------------
        ampl.eval("param alpha_arc{arcs};")
        ampl.eval("param y_sol{arcs};")

        for (i, j) in self.arcs:
            ampl.param["alpha_arc"][i, j] = self.alpha[seg_index[(i, j)]]
            ampl.param["y_sol"][i, j] = y_sol[i,j]

        ampl.eval("""
            subject to bound_y{(i,j) in arcs}:
                y[i,j] >= alpha_arc[i,j] - delta_y[i,j];

            subject to bound_delta_y{(i,j) in arcs}:
                delta_y[i,j] <= max(0, alpha[seg_index[i,j]] - y_sol[i,j]);
        """)

        # objective includes delta_y ONLY in phase 1
        # ampl.eval("""
        #     minimize phase1_obj:
        #         sum{(i,j) in arcs} (
        #             delta_headloss[i,j]^2 + delta_y[i,j]^2
        #         );
        # """)
        ampl.eval(f"""
           subject to bound_cost:
               sum{{(i,j) in arcs}}z[i,j]<={self.best_cost};
        """)
        ampl.solve()

        y_sol = ampl.get_variable("y").get_values().to_dict()
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)

        for (i,j) in self.arcs:
            yv = y_sol[i,j]
            for k in range(len(alpha_vals)-1):
                if alpha_vals[k+1] <= yv <= alpha_vals[k]:
                    seg_index[(i,j)] = k+1
                    break

        cost = sum(
           self.L[i,j] * (
               self.slope[seg_index[i,j]] * y_sol[i,j]
               + self.intercept[seg_index[i,j]]
           )
           for (i,j) in self.arcs
        )
        print("cost:", cost)
        # delta_y_sol = ampl.get_variable("delta_y").get_values().to_dict()
        # delta_headloss_sol = ampl.get_variable("delta_headloss").get_values().to_dict()
        # for (i,j) in delta_headloss_sol:
        #     print(f"delta_y[{i},{j}]:", delta_y_sol[i,j])
        #     print(f"delta_headloss[{i},{j}]:", delta_headloss_sol[i,j],"\n")

        # -------------------------
        # PHASE 2: FIX delta_y
        # -------------------------
        # ampl.eval("param delta_y_fix{(i,j) in arcs};")

        # for (i, j) in self.arcs:
        #     ampl.param["delta_y_fix"][i, j] = delta_y_sol[(i, j)]

        # ampl.eval("""
        #     subject to fix_delta_y{(i,j) in arcs}:
        #         delta_y[i,j] = delta_y_fix[i,j];
        # """)

        # remove phase1 objective
        # ampl.eval("drop phase1_obj;")

        # -------------------------
        # PHASE 3: continuation
        # -------------------------
        rho_list = [1e2, 1e4, 1e6, 1e8]
        # ampl.eval("""
        #     minimize phase2_obj:
        #         rho * (sum{(i,j) in arcs} (
        #             delta_headloss[i,j]^2 + delta_y[i,j]^2
        #         ));
        # """)
        max_iter = 10 
        iter = 1
        total_inf = 1e+5 
        while total_inf>=1e-5:
            print(f"\n--------------------------------iter:", iter,"---------------------------------")
            # ampl.param["rho"] = rho
            # for (i, j) in self.arcs:
            #     ampl.param["alpha_arc"][i, j] = self.alpha[seg_index[(i, j)]] - delta_y_sol[i,j]
            #
            if iter>=2:
                for (i,j) in self.arcs:
                    ampl.param["alpha_arc"][i,j] = y_corr[i,j]
                    ampl.param["seg_index"][i, j] = seg_index[(i, j)]
                    ampl.param["y_sol"][i, j] = y_sol[i,j]

            with self.suppress_output():
                ampl.solve()

            q_sol = ampl.get_variable("q").get_values().to_dict()
            y_sol = ampl.get_variable("y").get_values().to_dict()
            h_sol = ampl.get_variable("h").get_values().to_dict()
            z_sol = ampl.get_variable("z").get_values().to_dict()
            delta_y_sol = ampl.get_variable("delta_y").get_values().to_dict()
            delta_head = ampl.get_variable("delta_headloss").get_values().to_dict()

            head_inf = sum(abs(delta_head[i,j]) for (i,j) in delta_head)
            delta_y_inf = sum(abs(delta_y_sol[i,j]) for (i,j) in delta_y_sol)

            # for (i,j) in delta_head:
            #     print(f"delta_headloss[{i},{j}]:", delta_head[i,j],"\n")

            max_delta_head_inf = max(delta_head[i,j] for (i,j) in delta_head)
            max_delta_y_inf = max(delta_y_sol[i,j] for (i,j) in delta_y_sol)
            total_inf = np.abs(max_delta_head_inf) + np.abs(max_delta_y_inf)

            print(f"[iter={iter}] head_inf={head_inf:.3e} delta_y_inf={delta_y_inf:.3e} total_inf={total_inf:.3e}")

            # if total_inf < 1e-8:
            #     break

            y_corr = {}
            for (i,j) in self.arcs:
                phi = (q_sol[i,j]**3 * (q_sol[i,j]**2 + self.eps[i,j]**2)**0.426) / \
                      (q_sol[i,j]**2 + 0.426 * self.eps[i,j]**2)
            
                y_corr[i,j] = y_sol[i,j] + delta_head[i,j] / (phi * self.L[i,j])

                y_corr[i,j] = min(self.alpha_max, max(self.alpha_min, y_corr[i,j]))

            alpha_vals = sorted(set(self.alpha.values()), reverse=True)

            tol = 1e-8
            
            for (u, v) in self.arcs:
                y_val = y_corr[(u, v)]
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
            
                seg_index[(u, v)] = active_s

            iter += 1
        # -------------------------
        # PHASE 4: remove slack
        # -------------------------
        print("\n---------------------------Solve without slack-----------------------")
        ampl.eval("drop con2;")
        # ampl.eval("drop bound_y;")
        ampl.eval("drop infeas_obj;")
        ampl.eval("minimize total_cost:sum{(i,j) in arcs} z[i,j];")

        ampl.eval("""
            subject to headloss_exact{(i,j) in arcs}:
                h[i] - h[j]
              - ( q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426
                / (q[i,j]^2 + 0.426 * eps[i,j]^2))
                * y[i,j] * L[i,j] = 0;
        """)

        #with self.suppress_output():
        ampl.solve()

        # -------------------------
        # Extract solution
        # -------------------------
        q_sol = ampl.get_variable("q").get_values().to_dict()
        h_sol = ampl.get_variable("h").get_values().to_dict()
        y_sol = ampl.get_variable("y").get_values().to_dict()
        z_sol = ampl.get_variable("z").get_values().to_dict()
        # delta_y_sol = ampl.get_variable("delta_y").get_values().to_dict()
        # delta_head = ampl.get_variable("delta_headloss").get_values().to_dict()

        alpha_vals = sorted(set(self.alpha.values()), reverse=True)

        tol = 1e-8
        
        for (u, v) in self.arcs:
            y_val = y_sol[(u, v)]
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

        cost = sum(
            self.L[i,j] * (
                self.slope[self.seg_index[i,j]] * y_sol[i,j]
                + self.intercept[self.seg_index[i,j]]
            )
            for (i,j) in self.arcs
        )

        print("Final feasible cost:", cost)

        ampl.eval("display y;")

        ampl.eval("display delta_headloss;")
        ampl.eval("display delta_y;")
        print("\n-----------------Exit the Feasibility Restoration Phase------------")

        return q_sol, h_sol, y_sol, z_sol

    def solve_restricted_model1(self, q_sol, z_sol, y_sol, h_sol, seg_index):
        max_iter = 10 
        iter = 1
        while iter<=max_iter:
            ampl = AMPL()
            ampl.read("repair_feasibility_wdn.mod")
            ampl.read_data(self.data_file)

            ampl.option["solver"] = "ipopt"
            # ampl.option["ipopt_options"] = "outlev=0 tol=1e-9 max_iter=3000"
            ampl.option["ipopt_options"] = (
                "outlev=0 "
                "expect_infeasible_problem=no "
                "bound_relax_factor=0 "
                "bound_push=0.1 "
                "bound_frac=0.1 "
                "warm_start_init_point=yes "
                "halt_on_ampl_error=yes "
                "max_iter=3000"
            )

            # -------------------------
            ampl.eval("""
                param seg_index{(i,j) in arcs};
            """)
            ampl.eval("""
                param y_sol{(i,j) in arcs};
            """)
            # ampl.eval("""
            #     param x_sol{(i,j) in arcs} default 0;
            # """)
            # ampl.eval("""set 2_set := {1,2};""")
            # ampl.eval("""var x{arcs}>=0, <=1;""")
            # -------------------------
            # Warm-start perturbation
            # -------------------------
            for (i, j), val in q_sol.items():
                ampl.var["q"][i, j] = val
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

            try:
                ampl.eval(f"""
                    subject to bound_y_{{(i,j) in arcs}}:
                        y[i,j]>=alpha[seg_index[i,j]]-delta_y[i,j];
                    subject to bound_delta_y{{(i,j) in arcs}}:
                        delta_y[i,j]<=alpha[seg_index[i,j]] - y_sol[i,j];
                """)
                # ampl.eval(f"""
                #    subject to bound_cost:
                #        sum{{(i,j) in arcs}}y[i,j]>={sum(y_sol[i,j] for (i,j) in self.arcs)};
                # """)
                ampl.eval(f"""
                   subject to bound_cost:
                       sum{{(i,j) in arcs}}z[i,j]<={self.best_cost};
                """)
                print("Constraints added.")

            except Exception as e:
                print(f"Failed to solve restricted model: {e}")

            #with self.suppress_output():
            ampl.solve()

            if ampl.get_value("solve_result") != "solved":
                print("Solver failed. Stopping.")

            y_star_sum = sum(y_sol[i,j] for (i,j) in self.arcs)
            # -------------------------
            # Extract solution
            # -------------------------
            # obj = ampl.get_objective("total_infeasibility").value()
            q_sol = ampl.get_variable("q").get_values().to_dict()
            h_sol = ampl.get_variable("h").get_values().to_dict()
            y_sol = ampl.get_variable("y").get_values().to_dict()
            z_sol = ampl.get_variable("z").get_values().to_dict()
            # delta_flow_balance = ampl.get_variable("delta_flow_balance").get_values().to_dict()
            delta_headloss = ampl.get_variable("delta_headloss").get_values().to_dict()
            delta_y = ampl.get_variable("delta_y").get_values().to_dict()

            # ampl.eval("display z;")
            # ampl.eval("display q;")
            # ampl.eval("display h;")
            # ampl.eval("display y;")
            # ampl.eval("display delta_y;")
            # ampl.eval("display delta_flow_balance;")
            ampl.eval("display delta_headloss;")

            alpha_vals = sorted(set(self.alpha.values()), reverse=True)

            tol = 1e-8
            
            for (u, v) in self.arcs:
                y_val = y_sol[(u, v)]
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

            print("seg_index:", self.seg_index)
            # for (i,j) in self.arcs:
            #     print((i,j), z_sol[i,j] - self.L[i,j]*(self.slope[self.seg_index[i,j]]*y_sol[i,j] + self.intercept[self.seg_index[i,j]]))

            cost = sum(self.L[i,j]*(self.slope[self.seg_index[i,j]] * y_sol[i,j] + self.intercept[self.seg_index[i,j]]) for (i,j) in self.arcs)

            # total_flow_balance_infeasibility = sum(np.abs(delta_flow_balance[i]) for i in self.nodes if i not in self.source)
            total_headloss_infeasibility = sum(np.abs(delta_headloss[i,j]) for (i,j) in self.arcs)
            # print("total_flow_balance_infeasibility:", total_flow_balance_infeasibility) 
            print("total_headloss_infeasibility:", total_headloss_infeasibility) 
            # print("total_infeasibility:", obj)
            print("infeasible cost:", cost)
            
            # seg_index = ampl.getParameter('seg_index').to_dict()
            # print(self.alpha,"\n")
            # print("\n")
            # for (i,j) in self.arcs:
            #     # print(seg_index[i,j])
            #     print(f"alpha[{i},{j}] - alpha_min:",self.alpha[seg_index[i,j]] - delta_y[i,j])
            # print("\n")
            # print(self.alpha)
            # print(index)

            # ampl.eval("display {(i,j) in arcs}: y[i,j], alpha[seg_index[i,j]] - delta_y[i,j];")
            y_hat = {}
            for (i,j) in self.arcs:
                y_hat[i,j] = y_sol[i,j] + (delta_headloss[i,j]/(self.L[i,j]*q_sol[i,j] * np.abs(q_sol[i,j])**0.852))
                # print(f"y_hat[{i},{j}]-alpha_min:",y_hat[i,j]-self.alpha_min)
                if y_hat[i,j]>=self.alpha_min and y_hat[i,j]<=self.alpha_max:
                    print(y_hat[i,j])
                else:
                    diff1 = self.alpha_min - y_hat[i,j]
                    diff2 = y_hat[i,j] - self.alpha_max
                    if diff1>=1e-8:
                        print("diff1:",diff1)
                    else:
                        print("diff2:",diff2)

            # y_hat_sum = sum(y_hat[i,j] for (i,j) in self.arcs)

            # print("diff:",y_hat_sum-y_star_sum)

        return q_sol, h_sol, y_hat, z_sol

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
        # self.fixed_arcs = self.sorted_arcs
        # self.cycles = self.cycle_basis() 
        # print("sorted_arcs:", self.sorted_arcs)

        # self.plot_value_function_for_arc((1, 2))
        # self.plot_value_function_for_arcs([(19,3), (18,19)])
        # self.plot_m_values()
        self.plot_piecewise_lines((1,2))

        # print("\n-------------------------------- Solving Original Model --------------------------")

        # self.solve_original_model_without_init()
        
        # print("-------------REDUCED MODEL FEASIBILITY CHECK------------------")
        # y, z = self.map_original_to_reduced(self.l, self.q, self.h)
        # viol_red = self.check_reduced_feasibility(self.q, self.h, y, z)
        # print(viol_red)

        # self.print_piecewise_reduced_cost(arc=(1, 2))
        # # Plot segment-wise piecewise function
        # self.plot_piecewise_reduced_cost(
        #     arc=(1, 2),
        #     save_path="piecewise_reduced_cost_arc_1_2.png"
        # )
        #
        # # Plot continuous envelope
        # self.plot_piecewise_reduced_cost_envelope(
        #     arc=(1, 2),
        #     save_path="piecewise_reduced_cost_envelope_arc_1_2.png"
        # )

        ampl, solve_result, z_sol, q_sol, h_sol, y_sol, cost, solve_time = self.solve_exact_reduced_model()

        print(f"Total Cost: {cost:.8f}")
        print(f"Exact model solve time: {solve_time:.4f} sec")
        # ampl.eval("display exact_cost.dual;")
        # ampl.eval("display flow_balance.dual;")
        # ampl.eval("display h;")
        self.ampl = ampl
        self.q = q_sol
        self.h = h_sol
        self.y = y_sol
        self.z = z_sol
        self.best_cost = cost
        self.current_cost = cost
        self.best_solution = (self.q, self.h, self.y, self.z, self.best_cost)

        # self.ampl.eval("display y;")
        self.all_duals = {}
        for con_name, con in self.ampl.get_constraints():
            self.all_duals[con_name] = con.getValues()

        self.seg_index = {}
        # descending alpha breakpoints
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
        tol = 1e-8
        for (u, v) in self.arcs:
            y_val = self.y[(u, v)]
            active_segments = []
            for k in range(len(alpha_vals) - 1):
                a_high = alpha_vals[k]
                a_low  = alpha_vals[k + 1]
                # --------------------------------------------------------
                # interior of segment
                # --------------------------------------------------------
                if (a_low + tol) < y_val < (a_high - tol):
                    active_segments = [k + 1]
                    break
                # --------------------------------------------------------
                # breakpoint at upper boundary
                # shared by segment k and k-1
                # --------------------------------------------------------
                elif abs(y_val - a_high) <= tol:
                    if k == 0:
                        active_segments = [1]
                    else:
                        active_segments = [k, k + 1]
                    break
                # --------------------------------------------------------
                # breakpoint at lower boundary
                # shared by segment k and k+1
                # --------------------------------------------------------
                elif abs(y_val - a_low) <= tol:
                    if k == len(alpha_vals) - 2:
                        active_segments = [k + 1]
                    else:
                        active_segments = [k + 1, k + 2]
                    break
            # ------------------------------------------------------------
            # fallback
            # ------------------------------------------------------------
            if not active_segments:
                if y_val > alpha_vals[0]:
                    active_segments = [1]
                elif y_val < alpha_vals[-1]:
                    active_segments = [len(alpha_vals) - 1]
                else:
                    # numerical safeguard
                    distances = [
                        abs(y_val - alpha_vals[k])
                        for k in range(len(alpha_vals))
                    ]
                    idx = distances.index(min(distances))
                    if idx == 0:
                        active_segments = [1]
                    elif idx == len(alpha_vals) - 1:
                        active_segments = [len(alpha_vals) - 1]
                    else:
                        active_segments = [idx, idx + 1]
            self.seg_index[(u, v)] = active_segments
            # print(f"arc {(u,v)} -> y = {y_val:.6f}, segment = {active_s}")
        print(self.seg_index)
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

        print("\n-------------------------------- Segment Cut using Lagrangian Relaxation Based Heuristic --------------------------")
        best_global = {
            "q": q_sol.copy(), "h": h_sol.copy(),
            "y": y_sol.copy(), "z": z_sol.copy(),
            "seg": self.seg_index.copy(), "cost": self.best_cost
        }

        # best_global = incumbent.copy()
        # best_global["q"],best_global["h"],best_global["y"],best_global["z"] = self.update_dual_and_resolve(q_sol, h_sol, y_sol, z_sol, self.seg_index, strategy = "shift")
        # best_global["q"],best_global["h"],best_global["y"],best_global["z"] = self.lagrangian_based_heuristic(q_sol, h_sol, y_sol, z_sol, self.seg_index)
                
        # print("\n-------------------------------- Segment Cut Based Heuristic --------------------------")

        # best_solution = self.solve_with_boundary_pushing_cut(z_sol, y_sol, h_sol) 

        # q_sol, h_sol, y_sol, z_sol, cost = best_solution 
        # print("\n Best Solution:", cost)

        print("\n-------------------------------- Solve Restriced Model --------------------------")

        # q_sol, h_sol, y_sol, z_sol = self.solve_restricted_model(self.q, self.z, self.y, self.h, self.seg_index)
        # q_sol, h_sol, y_sol, z_sol = self.escape_and_project(self.q, self.h, self.y, self.z, self.seg_index)

        # print(q_sol)
        # viol_red = self.check_reduced_feasibility(q_sol, h_sol, y_sol, z_sol)
        # print(viol_red)
        print("\n-------------------------------- Branch and Bound based Heuristic --------------------------")
        self.network_graph = self.generate_random_acyclic_from_solution(q_sol)
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        self.best_acyclic_flow = self.network_graph.copy() 
        #
        self.iteration = 1
        self.visited_arc_reverse = []
        self.reversed_arcs = []
        self.segment_cut_based_heuristic()

        # best_global["q"],best_global["h"],best_global["y"],best_global["z"] = self.sequential_segment_cut_heuristic()
        # for (i,j) in self.arcs:
        #     print((i,j), self.z[i,j] - self.L[i,j]*(self.slope[self.seg_index[i,j]]*self.y[i,j] + self.intercept[self.seg_index[i,j]]))

        # print("\n-------------------------------- Solving Recover Model --------------------------")
        #
        # rec_ampl, rec_result, l_trial, recovered_cost, t_rec = self.solve_recover_model1(best_global["y"])
        # print(rec_result)
        # print(recovered_cost)

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
        # print("\n-------------------- Solving Original Model with Initialize ---------------------")
        #
        # orig_ampl, orig_result, q_new, h_new, l_new, final_cost, t_orig = self.solve_original_with_init(l_trial, self.q, self.h)
        # print(f"Original warm-start NLP result = {orig_result}, final cost = {final_cost:.8f}, time = {t_orig:.2f} sec")

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

