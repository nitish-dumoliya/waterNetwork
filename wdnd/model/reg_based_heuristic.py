from matplotlib.cbook import violin_stats
import networkx as nx
from amplpy import AMPL
import matplotlib.pyplot as plt
from networkx.algorithms import cycle_basis
import numpy as np
import time
import copy
import sys
import os
import contextlib
from sklearn.decomposition import PCA
import random
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from tabulate import tabulate
import optuna
import warnings
warnings.filterwarnings("ignore")
from pyswarm import pso
import multiprocessing
from network_layout import node_position
import math
class WaterNetworkOptimizer:
    def __init__(self, model_file, data_file, data_number, data_list):
        self.ampl = AMPL()
        self.model_file = model_file
        self.data_list = data_list
        self.data_file = data_file
        self.data_number = data_number
        self.total_cost = None
        self.network_graph = None
        self.solve_result = None
        self.solver_time = 0
        self.best_acyclic_flow = None
        self.number_of_nlp = 0
        self.tol = 0

    def load_model(self):
        """Load the model and data."""
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file) 
        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')
        if self.data_number == 6:
            self.fixarcs = self.ampl.getSet('fixarcs')
        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
        #self.eps = self.ampl.getParameter('eps').to_dict()
        self.delta = 0.1
        self.p = 1.852
        self.omega = 10.67

    def generate_random_acyclic_from_solution(self, q):
        # print("Generate the acyclic network using ipopt solution") 
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

    def fix_leaf_arc_flow(self):
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
                        # self.ampl.eval(f"s.t. fix_q_{edge[0]}_{edge[1]}: q[{edge[0]},{edge[1]}] = {flow_value};")
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
                            # self.ampl.eval(f"s.t. fix_q_{edge[0]}_{edge[1]}: q[{edge[0]},{edge[1]}] = {flow_value};")
                        elif neighbor == source:
                            flow_value = -D[leaf]
                            D[neighbor] = D[neighbor] - D[leaf] 
                            # self.ampl.eval(f"s.t. fix_q_{edge[0]}_{edge[1]}: q[{edge[0]},{edge[1]}] = {flow_value};")
                        else:
                            flow_value = -D[leaf]
                            D[neighbor] += -flow_value
                            # self.ampl.eval(f"s.t. fix_q_{edge[0]}_{edge[1]}: q[{edge[0]},{edge[1]}] = {flow_value};")
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

    def diameter_reduction(self):
        improved = False
        arc_max_dia = {}
        if self.data_number == 6:
            self.fixarcs = self.ampl.getSet('fixarcs')
            #print("fixarcs:",self.fixarcs)
            for (i, j, d), val in self.l.items():
                if (i,j) not in self.fixarcs or (j,i) not in self.fixarcs:
                    if val > 1e-3:
                        if (i, j) not in arc_max_dia:
                            arc_max_dia[(i, j)] = d
                        else:
                            arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)
        else:
            for (i, j, d), val in self.l.items():
                if val > 1e-3:
                    if (i, j) not in arc_max_dia:
                        arc_max_dia[(i, j)] = d
                    else:
                        arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)
        # print("\n*********************************************************************************************")
        print("Iteration :",self.dia_red_iteration)
        # self.sen_score = {}
        # for (i,j) in self.arcs:
        #     self.sen_score[i,j] = -1.852 * (self.h[i] - self.h[j])/(np.abs(self.q[i,j])**2.852)
        # print((i,j), self.sen_score[i,j])
        # print("sen_score:", self.sen_score)
        # edge = min(self.sen_score, key=self.sen_score.get) 
        # print("minimum sen_score:",edge, self.sen_score[edge[0], edge[1]])

        self.all_duals = {}
        for con_name, val in self.ampl.get_constraints():
            self.all_duals[con_name] = val.getValues()
        sorted_arcs = []
        dual_dict = self.all_duals["con2"].to_dict()
        # print("dual_value:", dual_dict)
        sorted_duals = dict(sorted(dual_dict.items(), key=lambda kv: abs(kv[1]), reverse=True))
        sorted_arcs = list(sorted_duals.keys())
        sorted_arcs = [arc for arc in sorted_arcs if arc not in self.visited_arc]
        # print("sorted_arcs:", self.visited_arc)
        # print("\nsorted_arcs:", sorted_arcs)
        if self.data_number == 6:
            sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fixarcs]
        sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fix_arc_set]
        sorted_arcs = [arc for arc in sorted_arcs if arc_max_dia[arc[0], arc[1]] != 1]
        
        sorted_arcs = [sorted_arcs[0]]
        print("sorted_arcs:", sorted_arcs)
        
        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")
        print("----------------------------------------------------------------------------------------")

        for (i,j) in sorted_arcs:
            self.visited_arc.append((i,j))
            ampl = AMPL()
            ampl.reset()
            if self.data_number==5:
                ampl.read("newyork_model.mod")
            elif self.data_number==6:
                ampl.read("blacksburg_model.mod")
            else:
                ampl.read("wdnmodel.mod")
            ampl.read_data(self.data_file)

            for (x, y, k), val in self.l.items():
                ampl.eval(f'let l[{x},{y},{k}] := {val};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number ==5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            for x, val in self.h.items():
                ampl.eval(f'let h[{x}] := {val};') 

            # current_duals = {}
            # for con_name, val in ampl.get_constraints():
            #     dual_values = val.get_values()
            #     current_duals[con_name] = dual_values
            #
            # # Initialize dual values for all constraints
            # for con_name, dual_values in self.all_duals.items():
            #     if con_name in current_duals:
            #         # Initialize dual values for each constraint
            #         ampl.get_constraint(con_name).set_values(dual_values) 

            # for (i,j) in self.visited_arc:
            #     ampl.eval(f"s.t. arc_dia{i}_{j}: sum{{k in pipes: k <=  {arc_max_dia[i,j]}}} l[{i},{j},k] = L[{i},{j}];")
            # if arc_max_dia[i,j]>=2:
            #     ampl.eval(f"s.t. arc_dia{i}_{j}: sum{{k in pipes: k <=  {arc_max_dia[i,j]}}} l[{i},{j},k] = L[{i},{j}];")
            # ampl.eval(f"s.t. arc_dia{i}_{j}: l[{i},{j},{arc_max_dia[i,j]-1}] + l[{i}, {j}, {arc_max_dia[i,j]}] = L[{i},{j}];")

            # self.update_initial_points_with_perturbation(ampl, self.l, self.q, self.h)
            # ampl.eval(f"""subject to con3_l_:sum {{(i,j) in arcs}} sum {{k in pipes}} C[k] * l[i,j,k] <= {self.current_cost};""")
            # ampl.eval(f"subject to con3_{i}_{j}: sum{{k in pipes: k <=  {arc_max_dia[i,j]-1}}} l[{i},{j},k] = L[{i},{j}];")
            # ampl.eval(f"subject to con3_l{i}_{j}_{arc_max_dia[i,j]}: l[{i},{j},{arc_max_dia[i,j]}] <= {self.l[i,j,arc_max_dia[i,j]]} - 1e-0;"
            # ampl.eval(f"subject to con3_l{i}_{j}: sum {{k in pipes}} (omega * l[{i},{j},k] / (R[k]^1.852 * d[k]^4.87)) <= {sum(10.67 * self.l[i,j,k]/(self.R[k]**1.852 * self.d[k]**4.87) for k in self.pipes)}*(1-0.);")
            # ampl.eval(f"subject to con3_l{i}_{j}: sum {{k in pipes}} C[k]*l[{i},{j},k] <= {self.C[arc_max_dia[i,j]] * self.L[i,j]};")
            
            for k in self.pipes:
                if k>=arc_max_dia[i,j]:
                    ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = 0;")
            
            # ampl.eval(f"s.t. flow_value{i}_{j}: q[{i}, {j}]>=abs({self.q[i,j]});")
            # if self.h[i] - self.h[j]>= 0:
            #     ampl.eval(f"s.t. headloss1{i}_{j}: h[{i}]- h[{j}]>={self.h[i] - self.h[j]} + 1e-0;")
            # else:
            #     ampl.eval(f"s.t. headloss2{i}_{j}: h[{i}]- h[{j}]<={self.h[i] - self.h[j]} - 1e-0;")
            # for key, val in self.l.items():
            #     u,v,k = key
            #     if (key[0], key[1]) != (i,j):
            #         ampl.eval(f"subject to con3__{u}_{v}_{k}: l[{u},{v},{k}] = {self.l[u,v,k]};")

            # pmin= {}
            # for u in self.nodes:
            #     if u in self.source:
            #         pass
            #     else:
            #         pmin[u] = self.E[u] + self.P[u]
            # key, value = max(pmin.items(), key=lambda x: x[1]) 
            # print(key, value)

            # ampl.eval(f"s.t. head_bound_{key}: h[{key}] >= {value};")

            # for k in self.pipes:
            #     if k==arc_max_dia[i,j]-1:
            #         ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = {self.L[i,j]};")
            #     else:
            #         ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = 0;")

            ampl.option['solver'] = "ipopt" 

            ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"
            #ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no  bound_relax_factor=0 warm_start_init_point = yes halt_on_ampl_error = yes"
            #with self.suppress_output():
            #ampl.option["presolve_eps"]= "7.19e-13"

            with self.suppress_output():
                ampl.solve()
            solve_time = ampl.get_value('_solve_elapsed_time')
            self.solver_time += solve_time
            self.number_of_nlp += 1

            l1 = ampl.getVariable('l').getValues().to_dict()
            q1 = ampl.getVariable('q').getValues().to_dict()
            h1 = ampl.getVariable('h').getValues().to_dict()
            # ampl.eval("display q;")
            # ampl.eval("display h;")
            total_cost = ampl.getObjective("total_cost").value()
            if ampl.solve_result == "solved":
                # self.plot_graph(self.super_source_out_arc, total_cost, 0, q, h, self.D, (0,0), l, self.C)
                if total_cost < self.current_cost:
                    print(f"{str((i,j)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(total_cost)):<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{ampl.solve_result:<14}{'Yes':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    self.current_cost = total_cost
                    self.ampl = ampl
                    improved = True
                    self.network_graph = self.generate_random_acyclic_from_solution(q1)
                    self.best_acyclic_flow = self.network_graph.copy()
                    self.indegree_2_or_more = [node for node, indeg in self.best_acyclic_flow.in_degree() if indeg >= 2]
                    # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
                    self.l = l1 
                    self.q = q1
                    self.h = h1 
                    # ampl.eval("display l;")
                    # ampl.eval("display l.rc;")
                    if self.data_number==5:
                        self.q1 = ampl.getVariable('q1').getValues().to_dict()
                        self.q2 = ampl.getVariable('q2').getValues().to_dict()
                    headvio = 0
                    for u in self.nodes:
                        if self.h[u]>=self.E[u] + self.P[u]:
                            pass
                        else:
                            headvio = headvio + self.E[u] + self.P[u] - self.h[u]
                    print("headvio:", headvio)
                    head_vio = 0 
                    #self.sorted_nodes = sorted(self.indegree_2_or_more, key=lambda node: self.D[node], reverse=True)
                    #print("\nvisited_nodes:", self.visited_nodes)
                    #if self.visited_nodes:
                    #    self.sorted_nodes = [item for item in self.sorted_nodes if item not in self.visited_nodes]
                    #print("sorted_nodes", self.sorted_nodes) 
                    self.fix_arc_set = list(set(self.super_source_out_arc) | set(self.fix_arc_set))
                    print("----------------------------------------------------------------------------------------")
                else:
                    print(f"{str((i,j)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(total_cost)):<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{ampl.solve_result:<14}{'No':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
            else:
                print(f"{str((i,j)):<10}" 
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{self.format_indian_number(round(total_cost)):<14}" 
                    f"{(str(round(solve_time, 2)) + 's'):<12}"
                    f"{ampl.solve_result:<14}{'No':<10}"
                    f"{round(time.time() - self.start_time, 2)}s ")
            # headloss = 0
            # for (u,v) in self.arcs:
            #     headloss = headloss + np.abs(h1[u] - h1[v])
            # print("headloss:", headloss)


            if improved:
                self.dia_red_iteration = self.dia_red_iteration + 1
                self.diameter_reduction()
                break


    def run(self):
        """Main function to run the Heuristic Approach."""
        self.start_time = time.time()
        self.bound_push , self.bound_frac = (0.1, 0.1)
        self.mu_init = 0.1
        # print("NLP solve using:  smooth approximation 1, Epsilon selection using absolute error\n")
        # print("NLP solve using: smooth approximation 1, epsilon selection using relative error\n")
        # print("NLP solve using: smooth approximation 2, epsilon selection using absolute error\n")
        print("NLP Model: Smooth Approximatate WDN Design Model 2, Epsilon Selection using Relative Error")
        print("NLP Solver: Ipopt")
        print("************************Initial Solution of Approximate WDN Design Model**********************")
        self.load_model()
        fix_arc_set = self.fix_leaf_arc_flow()
        print("fix_arc_set:",fix_arc_set)
        self.super_source_out_arc = self.fix_arc_set()
        print("super_source_out_arc:", self.super_source_out_arc, "\n")

        self.solve()
        print("Objective: ",self.total_cost)
        print("Solve_result: ",self.solve_result)
        print("Solve_time:", self.ampl.get_value('_solve_elapsed_time'),"\n")
        if self.solve_result != 'solved':
            print("IPOPT did not solve the initial problem optimally. Exiting Heuristic.")
            sys.exit()
        self.current_cost = self.total_cost
        self.l = self.ampl.getVariable('l').getValues().to_dict()
        self.q = self.ampl.getVariable('q').getValues().to_dict()
        self.h = self.ampl.getVariable('h').getValues().to_dict()
        # for (i,j) in self.arcs:
        #     print(f"h[{i}] - h[{j}]:", self.h[i] - self.h[j])
        self.ampl.eval("display l;")
        self.ampl.eval("display q;")
        self.ampl.eval("display h;")
        self.ampl.eval("display {i in nodes} E[i] + P[i];")
        self.ampl.eval("display {i in nodes} h[i] - E[i] - P[i];")
        # self.ampl.eval("display l.rc;")
        # self.ampl.eval("display l;")
        # self.ampl.eval("display l.rc;")
        # headloss = 0
        # for (i,j) in self.arcs:
        #     headloss = headloss + np.abs(self.h[i] - self.h[j])
        # print("headloss:", headloss)
        if self.data_number==5:
            self.q1 = self.ampl.getVariable('q1').getValues().to_dict()
            self.q2 = self.ampl.getVariable('q2').getValues().to_dict()
        # self.plot_graph(self.super_source_out_arc, self.current_cost, 0, self.q, self.h, self.D, (0,0), self.l, self.C)
        print("*****************************Improve the Initial Solution*************************************\n")
        self.super_source_out_arc = self.fix_arc_set()
        self.network_graph = self.generate_random_acyclic_from_solution(self.q)
        # print("Fix the flow direction in optimization model and solve the updated model")
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set)        
        self.best_acyclic_flow = self.network_graph.copy() 
        self.visited_nodes = []
        self.sorted_nodes = []
        # self.visited_arc = []
        # self.plot_graph(fix_arc_set, self.total_cost, 0, self.q, self.h, self.D, (0,0), self.l, self.C)
        print("\n----------------------------Diameter Reduction Approach------------------------------------")
        self.dia_red_iteration = self.headloss_increase_iteration + 1
        self.visited_arc = []
        self.diameter_reduction()

       print("\n************************************Final Best Results*****************************************")
        print("Water Network:", self.data_list[self.data_number])
        # self.eps = self.ampl.get_variable('eps').get_values().to_dict()
        self.eps = self.ampl.getParameter('eps').getValues().to_dict()
        # cost = 0
        # for (i,j) in self.arcs:
        #     for k in self.pipes:
        #         cost = cost + self.C[k]*self.l[i,j,k]
        # print("cost:", cost)
        print(f"Final best objective: {self.current_cost}")
        #self.ampl.eval("display q;")
        print("Number of nlp problem solved:", self.number_of_nlp)
        print("Total number of iteration:", self.iteration + self.headloss_increase_iteration + self.dia_red_iteration)
        # self.constraint_violations(self.q, self.h, self.l, self.eps, "ipopt")
        elapsed_time = time.time() - self.start_time
        solver_time = self.solver_time
        print(f"Solver_time: {solver_time:.2f} seconds")
        # print(f"Heuristic elapsed time: {elapsed_time:.2f} seconds = {elapsed_time/60:.2f} minutes.\n")
        print(f"Heuristic elapsed time: {elapsed_time:.2f} seconds")
        print("***********************************************************************************************\n")

if __name__ == "__main__":
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
    # Select the data number here (0 to 18)
    #file_name = sys.argv[1]
    data_number = int(sys.argv[1]) - 1
    #input_data_file = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{file_name}"
    input_data_file = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{data_list[(data_number)]}.dat"
    # print(data_number)
    print("Water Network:", data_list[(data_number)])
    if data_number==5:
        optimizer = WaterNetworkOptimizer("newyork_model.mod", input_data_file, data_number, data_list)
    elif data_number == 6:
        optimizer = WaterNetworkOptimizer("blacksburg_model.mod", input_data_file, data_number, data_list)
    else:
        optimizer = WaterNetworkOptimizer("wdnmodel.mod", input_data_file, data_number, data_list)
    # optimizer = WaterNetworkOptimizer(sys.argv[1], sys.argv[2], sys.argv[3])
    optimizer.run()
