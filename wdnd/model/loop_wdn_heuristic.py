import networkx as nx
from amplpy import AMPL
import matplotlib.pyplot as plt
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

    def create_digraph(self):
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = self.ampl.getSet('arcs').to_list()
        self.network_graph = nx.DiGraph()
        self.network_graph.add_nodes_from(nodes_list)
        self.network_graph.add_edges_from(edges_list)
        print(nodes_list)
        print(edges_list)

    def display_results(self):
        """Display relevant results from the optimization."""
        self.ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>1} l[i,j,k];")
        self.ampl.eval("display {(i,j) in arcs}: q[i,j];")
        self.ampl.eval("display h;")
        self.ampl.eval("display solve_result;")
        self.total_cost = self.ampl.get_objective("total_cost").value()
        print("Objective:", self.total_cost)
    
    def plot_graph(self, super_source_out_arc=None, current_cost = None, iteration = 1, edge_weights= None, h = None, D = None, arc=(0,0), l ={}, C = {}):           
        pos = node_position[self.data_number]
        # # Update node positions
        # pos = node_coordinate_d1
        # pos = nx.spectral_layout(self.network_graph)
        cost = {}
        for (i,j) in self.ampl.getSet('arcs'):
            cost[i,j] = sum(l[i,j,k] * C[k] for k in self.pipes)
        plt.figure(figsize=(10, 8))
        # plt.figure(figsize=(15, 11))
        cmap = plt.cm.plasma
        network_graph = self.generate_random_acyclic_from_solution(edge_weights)
        nx.draw_networkx_nodes(network_graph, pos, node_color='lightblue',edgecolors="black", node_size=300,linewidths=0.5, label='Regular Nodes') 
        indegree_2_or_more = [node for node, indeg in network_graph.in_degree() if indeg >= 2]
        if indegree_2_or_more:
            nx.draw_networkx_nodes(network_graph, pos, nodelist=indegree_2_or_more, node_color='orange',edgecolors="orange", node_size=300, label='Nodes with In-Degree ≥ 2')
        # nx.draw_networkx_nodes(self.network_graph, pos, nodelist=list(self.source), node_color='cornflowerblue',edgecolors="black", node_size=300,linewidths=0.5, label='Source node')
        # if not self.visited_nodes:
        #     nx.draw_networkx_nodes(self.network_graph, pos, nodelist=indegree_2_or_more, node_color='orange',edgecolors="black", node_size=300,linewidths=0.5, label='Nodes with In-Degree ≥ 2')
        # if self.sorted_nodes:
        #     nx.draw_networkx_nodes(self.network_graph, pos, nodelist=self.sorted_nodes, node_color='orange',edgecolors="black", node_size=300,linewidths=0.5, label='Visited nodes')
        # if self.visited_nodes:
        #     visited_nodes = [item for item in  self.visited_nodes if item in self.indegree_2_or_more]
        #     nx.draw_networkx_nodes(self.network_graph, pos, nodelist=visited_nodes, node_color='violet',edgecolors="black", node_size=300,linewidths=0.5, label='Sorted nodes')
        nx.draw_networkx_labels(network_graph, pos, font_size=10)
        nx.draw_networkx_edges(network_graph, pos, arrowstyle="->", arrowsize=12, edge_color='black', label='Regular Arcs', arrows=True) # arrows=False
        if super_source_out_arc:
            nx.draw_networkx_edges(network_graph, pos, edgelist=super_source_out_arc,arrowstyle="->", arrowsize=12, edge_color='red', width=1, label='Fix arc direction')
        # if best_arc:
        #     nx.draw_networkx_edges(self.network_graph, pos, edgelist=[best_arc],arrowstyle="->", arrowsize=12, edge_color='magenta', width=1, label = 'Best arc')
        # Annotate node demands
        if h:
            for node, (x, y) in pos.items():
                demand = h.get(node, 0)  # Get the head for the node, default to 0 if not in dictionary
                # plt.text(x, y + 120, f"{demand:.2f}", fontsize=8, color='blue', ha='center')  # Annotate demand below the node
                # plt.text(x, y + 200, f"{demand:.2f}", fontsize=8, color='blue', ha='center')  # Annotate demand below the node
                # plt.text(x, y -100 , f"{demand:.2f}", fontsize=8, color='magenta', ha='center')  # Annotate demand below the node
                # plt.text(mid_x, mid_y + 50 , f"{weight:.2f}", fontsize=8, color='green', ha='center')  # Annotate weight on edge
                # plt.text(mid_x, mid_y + 100, f"{value:.2f}", fontsize=8, color='green', ha='center')  # Annotate weight on edge
        if D:
            for node, (x, y) in pos.items():
                demand = D.get(node, 0)  # Get the demand for the node, default to 0 if not in dictionary
                # plt.text(x-250, y - 400, f"{demand*1000:.2f}", fontsize=10, color='black', ha='center')  # Annotate demand below the node
                plt.text(x-80, y - 130 , f"{demand:.2f}", fontsize=10, color='magenta', ha='center')  # Annotate demand below the node
        if edge_weights:
            for (u, v), weight in edge_weights.items():
                # if self.network_graph.has_edge(u, v):
                mid_x = (pos[u][0] + pos[v][0]) / 2  # Midpoint x-coordinate
                mid_y = (pos[u][1] + pos[v][1]) / 2  # Midpoint y-coordinate
                # plt.text(mid_x+270, mid_y+180 , f"{weight:.2f}", fontsize=10, color='black', ha='center')  # Annotate weight on edge
                plt.text(mid_x+150, mid_y + 80 , f"{weight*1000:.2f}", fontsize=10, color='black', ha='center')  # Annotate weight on edge
        # if cost:
        #     for (u, v), value in cost.items():
        #         # if self.network_graph.has_edge(u, v):
        #         mid_x = (pos[u][0] + pos[v][0]) / 2  # Midpoint x-coordinate
        #         mid_y = (pos[u][1] + pos[v][1]) / 2  # Midpoint y-coordinate
                # plt.text(mid_x, mid_y-300 , f"{round(value)}", fontsize=8, color='purple', ha='center')  # Annotate weight on edge
                # plt.text(mid_x, mid_y + 100, f"{value:.2f}", fontsize=8, color='green', ha='center')  # Annotate weight on edge
        pipe_dia_arc = {}
        for (i,j) in self.arcs:
            list_=[]
            for k in self.pipes:
                if l[i,j,k]>= 1e-5:
                    list_.append(k)
            pipe_dia_arc[i,j] = list_ 
        if pipe_dia_arc:
            for (u, v), weight in pipe_dia_arc.items():
                # if self.network_graph.has_edge(u, v):
                mid_x = (pos[u][0] + pos[v][0]) / 2  # Midpoint x-coordinate
                mid_y = (pos[u][1] + pos[v][1]) / 2  # Midpoint y-coordinate
                plt.text(mid_x - 100, mid_y - 100 , f"{weight}", fontsize=10, color='purple', ha='center')  # Annotate weight on edge
                # plt.text(mid_x-250, mid_y+120 , f"{weight}", fontsize=11, color='#006400', ha='center')  # Annotate weight on edge
        
        regular_node_patch = mpatches.Patch(color='lightblue', label='Regular Nodes')
        indegree_node_patch = mpatches.Patch(color='orange', label='Nodes with In-Degree ≥ 2')
        regular_edge_line = mlines.Line2D([], [], color='black', label='Regular Arcs')
        super_source_edge_line = mlines.Line2D([], [], color='red', label='Fix arc direction')
        # best_edge_line = mlines.Line2D([], [], color='magenta', label='Best Arc')
        # plt.legend(handles=[regular_node_patch, indegree_node_patch, regular_edge_line, super_source_edge_line], loc='lower right')
        
        cost = round(current_cost)
        # res = f"{cost:,}"
        plt.title(f"Total cost: {self.format_indian_number(cost)}")
        # (u,v) = arc
        plt.savefig(f"/home/nitishdumoliya/waterNetwork/wdnd/figure/newfigure/d{self.data_number}_iteration_{iteration}.png")
        plt.box(False)
        plt.show()
    
    def cycle_basis(self):
        root = self.ampl.getSet('Source').to_list()[0]
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = self.ampl.getSet('arcs').to_list() 
        uwg = nx.Graph()
        uwg.add_nodes_from(nodes_list)
        uwg.add_edges_from(edges_list)
        # print("Edges in the undirected graph:", edges_list)
        print("cycle basis for given water network: ",nx.cycle_basis(uwg, root))
        
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
 
    def generate_random_acyclic_graph(self):
        uwg = nx.Graph()
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = self.ampl.getSet('arcs').to_list()
        uwg.add_nodes_from(nodes_list)
        uwg.add_edges_from(edges_list)
        print("Edges in the undirected graph:", edges_list) 
        # Generate a random spanning tree using Wilson's algorithm
        random_tree = nx.random_spanning_tree(uwg) 
        # Retrieve the root from the AMPL source set
        root_l = self.ampl.getSet('Source').to_list()
        root = root_l[0]
        print("Root node:", root)
        # Ensure the root is present in the random tree
        if root not in random_tree.nodes:
            raise ValueError("The specified root must be a node in the graph.")
        # Create a directed graph from the random tree starting from the specified root
        self.network_graph = nx.DiGraph()
        visited = set()
        def dfs(node):
            visited.add(node)
            for neighbor in random_tree.neighbors(node):
                if neighbor not in visited:
                    self.network_graph.add_edge(node, neighbor) 
                    dfs(neighbor)
        # Start DFS from the specified root
        dfs(root)
        # Draw the initial directed tree
        plt.figure(figsize=(15, 10))
        plt.subplot(121)
        nx.draw_spectral(self.network_graph, with_labels=True, node_color='lightgreen', font_weight='bold', arrows=True)
        plt.title("Directed Spanning Tree")
        # Add remaining edges from the original graph and check for cycles
        for u, v in uwg.edges():
            if not self.network_graph.has_edge(u, v):  
                self.network_graph.add_edge(u, v)  
                if not nx.is_directed_acyclic_graph(self.network_graph):  
                    self.network_graph.remove_edge(u, v)  
                    self.network_graph.add_edge(v, u)  
        # Draw the final directed graph after adding remaining edges
        plt.subplot(122)
        nx.draw_spectral(self.network_graph, with_labels=True, node_color='lightgreen', font_weight='bold', arrows=True)
        plt.title("Acyclic Directed Graph")
        plt.show()

    def update_model(self):
        # print("Fix the arcs direction using the acyclic network\n")
        edges_list = [(arc[0],arc[1]) for arc in self.ampl.getSet('arcs')]
        for edge in self.network_graph.edges:
            i, j = edge
            if edge in edges_list:
                self.ampl.eval(f"s.t. flow_direction{i}_{j}: q[{i},{j}] >=0;")
                # self.ampl.eval(f"s.t. head_bound_left{i}_{j}: E[{j}]+P[{j}] <= h[{j}];")
                # self.ampl.eval(f"s.t. head_bound_right{i}_{j}: E[{j}] + P[{j}] <= h[{i}];")
            else:
                self.ampl.eval(f"s.t. flow_direction{i}_{j}: q[{j},{i}] <=0;")
                # self.ampl.eval(f"s.t. head_bound_left{i}_{j}: E[{i}]+P[{i}] <= h[{i}];")
                # self.ampl.eval(f"s.t. head_bound_right{i}_{j}: E[{i}] + P[{i}] <= h[{j}];")  

    def is_valid_edge(self, source, target):
        """Check if adding the directed edge (source -> target) maintains acyclicity."""
        self.network_graph.add_edge(source, target)  # Temporarily add the edge
        is_dag = nx.is_directed_acyclic_graph(self.network_graph)  # Check for acyclicity
        self.network_graph.remove_edge(source, target)  # Remove the edge after checking
        return is_dag  # Return True if it maintains acyclicity     
    
    def check_incoming_arcs(self):
        root = list(self.ampl.getSet('Source'))[0]
        # Iterate over all nodes in the graph
        for node in self.network_graph.nodes():
            # Skip the root node
            if node == root:
                continue
            # Check if the in-degree of the node is at least 1
            if self.network_graph.in_degree(node) < 1:
                # print(f"Node {node} does not have any incoming arcs.")
                return False
        # print("All nodes except the root have at least one incoming arc.")
        return True
            
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

    def update_initial_points(self,l_solution, q_solution, h_solution):
        for (i, j, k), val in l_solution.items():
            self.ampl.eval(f'let l[{i},{j},{k}] := {val};')
        for (i, j), val in q_solution.items():
            self.ampl.eval(f'let q[{i},{j}] := {val};')
            if self.data_number ==5:
                self.ampl.eval(f'let q1[{i},{j}] := {self.q1[i,j]};')
                self.ampl.eval(f'let q2[{i},{j}] := {self.q2[i,j]};')
        for i, val in h_solution.items():
            self.ampl.eval(f'let h[{i}] := {val};')

    def update_initial_points1(self,l_solution, q_solution, h_solution, inarc):
        for (i, j, k), val in l_solution.items():
            self.ampl.eval(f'let l[{i},{j},{k}] := {val};')        
        edge_list_network = self.network_graph.edges
        for i, val in h_solution.items():
            self.ampl.eval(f'let h[{i}] := {val};') 
        for (i, j), val in q_solution.items():
            edge = (i, j)
            if edge in edge_list_network: 
                if (i,j) not in inarc:
                    # print(f"self.ampl.eval(let q[{i},{j}] := {val};)")
                    self.ampl.eval(f"let q[{i},{j}] := {val} ;")
                else:
                    # print(f"self.ampl.eval(let q[{i},{j}] := {-val};)")
                    self.ampl.eval(f"let q[{i},{j}] := {val} ;")
            else:
                if (j,i) not in inarc:
                    # print(f"self.ampl.eval(let q[{i},{j}] := {val};)")
                    self.ampl.eval(f"let q[{i},{j}] := {val};")
                else:
                    # print(f"self.ampl.eval(let q[{j},{i}] := {-val};)")
                    self.ampl.eval(f"let q[{j},{i}] := {val} ;") 
        # current_duals = {}
        # for con_name, val in self.ampl.get_constraints():
        #     dual_values = val.get_values()
        #     current_duals[con_name] = dual_values

        # Initialize dual values for all constraints
        # for con_name, dual_values in all_duals.items():
            # if con_name in current_duals:
                # Initialize dual values for each constraint
                # self.ampl.get_constraint(con_name).set_values(dual_values)
            # else:
            #     print(f"Skipping initialization for constraint: {con_name} (not in current model)")

    def update_initial_points_with_perturbation1(self, l_solution, q_solution, h_solution,all_duals, inarc, delta=0.1):
        edge_list_network = self.network_graph.edges
        L = self.ampl.getParameter('L').getValues().to_dict()
        # Perturb l values
        for (i, j, k), val in l_solution.items():
            if (i,j) not in inarc:
                if val>= 1e-5:
                    perturbation = random.gauss(0, 1)
                    new_val = val + perturbation
                    # if val >= 1e-5:
                    self.ampl.eval(f'let l[{i},{j},{k}] := {new_val};')
                else:
                    self.ampl.eval(f'let l[{i},{j},{k}] := {0};')
            else:
                if val>= 1e-5:
                    perturbation = random.gauss(0, 1)
                    new_val = val + perturbation
                    self.ampl.eval(f'let l[{i},{j},{k}] := {new_val};')
                else:
                    self.ampl.eval(f'let l[{i},{j},{k}] := {0};')            
        # Perturb h values
        for i, val in h_solution.items():
            perturbation = random.gauss(0, 1)
            new_val = val + perturbation
            self.ampl.eval(f'let h[{i}] := {new_val};')
        # Modify q values based on heuristic
        for (i, j), val in q_solution.items():
            edge = (i, j)
            perturbation = random.gauss(0, 1)
            if edge in edge_list_network:
                if (i, j) not in inarc:
                    self.ampl.eval(f"let q[{i},{j}] := {val + perturbation};")
                else:
                    self.ampl.eval(f"let q[{i},{j}] := {(val + perturbation)};")
            else:
                if (j, i) not in inarc:
                    self.ampl.eval(f"let q[{i},{j}] := {val + perturbation};")
                else:
                    self.ampl.eval(f"let q[{j},{i}] := {(val + perturbation)};") 
        current_duals = {}
        for con_name, val in self.ampl.get_constraints():
            dual_values = val.get_values()
            current_duals[con_name] = dual_values
        # Initialize dual values for all constraints
        for con_name, dual_values in all_duals.items():
            if con_name in current_duals:
                # Initialize dual values for each constraint
                self.ampl.get_constraint(con_name).set_values(dual_values)
    
    def update_initial_points_with_perturbation(self, ampl, l_solution, q_solution, h_solution):
        delta = 20
        # Perturb l values
        for (i, j, k), val in l_solution.items():
            perturbation = random.gauss(0, delta)
            new_val = val + perturbation
            ampl.eval(f'let l[{i},{j},{k}] := {new_val};')
            # if val>= 1e-5:
            #     perturbation = random.gauss(0, delta)
            #     new_val = val + perturbation
            #     ampl.eval(f'let l[{i},{j},{k}] := {new_val};')
            # else:
            #     ampl.eval(f'let l[{i},{j},{k}] := {0};')
        # Perturb h values
        for i, val in h_solution.items():
            perturbation = random.gauss(0, delta)
            new_val = val + perturbation
            ampl.eval(f'let h[{i}] := {new_val};')
        # Modify q values based on heuristic
        for (i, j), val in q_solution.items():
            edge = (i, j)
            perturbation = random.gauss(0, delta)
            if val >= 0:
                ampl.eval(f"let q[{i},{j}] := {val + perturbation};")
            else:
                ampl.eval(f"let q[{i},{j}] := {val + perturbation};") 
        # current_duals = {}
        # for con_name, val in ampl.get_constraints():
        #     dual_values = val.get_values()
        #     current_duals[con_name] = dual_values
        #
        # # Initialize dual values for all constraints
        # for con_name, dual_values in all_duals.items():
        #     if con_name in current_duals:
        #         # Initialize dual values for each constraint
        #         ampl.get_constraint(con_name).set_values(dual_values)
 
    def multistart(self, inarc, current_cost, best_acyclic_flow, improved, super_source_out_arc, iteration):
        improved = False
        max_l = max(self.ampl.getParameter('L').to_dict().values())
        max_q = self.ampl.getParameter('D').getValues().toDict()
        # print(max_q[1])
        # for i, value in max_q.items():
        #     print(value) 
        source = self.ampl.getSet('Source').to_list()
        E = self.ampl.getParameter('E').getValues().toDict()
        P = self.ampl.getParameter('P').getValues().toDict() 
        # Define the number of starts for multistart heuristic
        num_starts = 10 
        # Set a random seed for reproducibility
        random.seed(num_starts) 
        # Loop for multistart heuristic
        for start in range(num_starts): 
            for (i,j) in self.arcs:
                for k in self.pipes:
                    value = random.uniform(0, max_l)  
                    self.ampl.eval(f' let l[{i},{j},{k}] := {self.l[i,j,k]};')
            for (i,j) in self.network_graph.edges:
                if (i,j) not in inarc:
                    if (i,j) in self.arcs:
                        value = random.uniform(0, -max_q[1])
                        self.ampl.eval(f'let q[{i},{j}] := {value};')
                    else:
                        value = random.uniform(max_q[1], 0)
                        self.ampl.eval(f'let q[{j},{i}] := {value};')
                else:
                    if (i,j) in self.arcs:
                        value = random.uniform(0, -max_q[1])
                        self.ampl.eval(f'let q[{i},{j}] := {value};')
                    else:
                        value = random.uniform(max_q[1], 0)
                        self.ampl.eval(f'let q[{j},{i}] := {value};')
            # for (i,j) in self.ampl.get_set("arcs").to_list():
            #     value = random.uniform(max_q[1], -max_q[1])  
            #     self.ampl.eval(f'let q[{i},{j}] := {value};')
            for i in self.nodes:
                value = random.uniform(E[i]+P[i], E[source[0]])  
                self.ampl.eval(f'let h[{i}] := {self.h[i]};')
            self.ampl.set_option("solver", "ipopt")
            self.ampl.set_option("ipopt_options", "outlev = 0 expect_infeasible_problem = yes bound_push = 0.01 bound_frac = 0.001 warm_start_init_point = yes halt_on_ampl_error = yes")
            # ampl.option[""]
            self.ampl.solve()
            if self.ampl.solve_result == 'solved':
                cost = self.ampl.get_objective("total_cost").value()
                print(cost)
                # Update the best solution if the current cost is lower
                if cost < current_cost:
                    improved = True
                    current_cost = cost
                    self.generate_random_acyclic_from_solution()
                    best_acyclic_flow = self.network_graph.copy()
                    l = self.ampl.getVariable('l').getValues().to_dict()
                    q = self.ampl.getVariable('q').getValues().to_dict()
                    h = self.ampl.getVariable('h').getValues().to_dict()
                    D = self.ampl.getParameter('D').getValues().to_dict()
                    self.plot_graph(super_source_out_arc, current_cost, iteration, q, h, D)
                    # print(best_acyclic_flow, improved, current_cost, l, q, h)
                # else:
                #     pass
        if improved:
            return best_acyclic_flow, improved, current_cost, l, q, h
        else:
            return self.best_acyclic_flow, improved, current_cost, self.l, self.q, self.h

    def acyclic_arcs(self):
        network_graph = self.best_acyclic_flow
        indegree_2_or_more = [node for node, indeg in network_graph.in_degree() if indeg >= 2]
        acyclic_arc = set()
        for node in indegree_2_or_more:
            # print("Node:", node,"in_degree:", self.network_graph.in_degree(node))
            for edge in list(network_graph.in_edges(node)):
                (u, v) = edge
                if (u,v) not in self.super_source_out_arc :
                    network_graph.remove_edge(u,v)
                    network_graph.add_edge(v,u)
                    acy_check = nx.is_directed_acyclic_graph(network_graph)
                    in_arc_check = self.check_incoming_arcs()
                    # print("Acyclic", acy_check and in_arc_check)
                    if acy_check and in_arc_check:
                        acyclic_arc.add((u,v))
                        # if (u,v) in self.arcs:
                        #     acyclic_arc.add((u,v))
                        # else:
                        #     acyclic_arc.add((v,u))
                    network_graph.remove_edge(v, u)
                    network_graph.add_edge(u, v)
        return acyclic_arc
    
    def constraint_violations(self, q_values, h_values, l_values, epsilon, solver):
        total_absolute_constraint_violation = 0
        total_relative_constraint_violation = 0 
        con1_gap = {}
        if self.data_number==5:
            q1 = self.q1
            q2 = self.q2
            for i in self.nodes:
                if i not in self.source:
                    con1_rhs = self.D[i]
                    incoming_flow = sum(q1[j, i] + q2[j,i] for j in self.nodes if (j, i) in self.arcs)
                    outgoing_flow = sum(q1[i, j] + q2[i,j] for j in self.nodes if (i, j) in self.arcs)
                    con1_lhs = incoming_flow - outgoing_flow
                    con1_violation = con1_lhs - con1_rhs
                    con1_gap[f"{i}"] = con1_violation
                    total_absolute_constraint_violation += abs(con1_violation)
        else:
            for i in self.nodes:
                if i not in self.source:
                    con1_rhs = self.D[i]
                    incoming_flow = sum(q_values[j, i] for j in self.nodes if (j, i) in self.arcs)
                    outgoing_flow = sum(q_values[i, j] for j in self.nodes if (i, j) in self.arcs)
                    con1_lhs = incoming_flow - outgoing_flow
                    con1_violation = con1_lhs - con1_rhs
                    con1_gap[f"{i}"] = con1_violation
                    total_absolute_constraint_violation += abs(con1_violation)
        #print("con1_gap:", con1_gap) 
        con2_original_gap = {}
        con2_approx_gap = {}
        absolute_violations = {}
        relative_violations = {}
        con2_absolute_constraint_violation = 0
        con2_relative_constraint_violation = 0
        con2_original_violation = 0
        con2_approx_violation = 0

        if self.data_number==5:
            #q1 = self.ampl.get_variable('q1').get_values().to_dict()
            #q2 = self.ampl.get_variable('q2').get_values().to_dict()
            self.exdiam = self.ampl.getParameter('exdiam').to_dict()
            for (i, j) in q1.keys():
                # Original constraint value
                lhs = 2*(h_values[i] - h_values[j])
                alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[i,j] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                original_rhs = q1[i, j] * (abs(q1[i, j])) ** 0.852 * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + q2[i, j] * (abs(q2[i, j])) ** 0.852 * alpha_rhs  
                #original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs
                # Approximated constraint value
                approx_rhs = (q1[i, j]**3 * ((q1[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q1[i,j]**2 + 0.426*epsilon[i,j]**2)) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j]**3 * ((q2[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q2[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                # approx_rhs = (q1[i, j] * ((q1[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j] * ((q2[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

                # approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs
                con2_original_gap[f"{i},{j}"] = lhs - original_rhs
                con2_original_violation += abs(lhs - original_rhs) 
                con2_approx_gap[f"{i},{j}"] = lhs - approx_rhs 
                total_absolute_constraint_violation += abs(lhs - approx_rhs)    
                con2_approx_violation += abs(lhs - approx_rhs) 

                # Compute absolute violation
                absolute_violation =  original_rhs - approx_rhs
                absolute_violations[f"{i},{j}"] = absolute_violation
                con2_absolute_constraint_violation += abs(absolute_violation)

                # Compute relative violation between original_rhs and approx_rhs
                relative_violation = (original_rhs - approx_rhs) / (original_rhs+1e-14)
                relative_violations[f"{i},{j}"] = relative_violation
                con2_relative_constraint_violation += abs(relative_violation)

        elif self.data_number==6:
            self.fixarcs = self.ampl.getSet('fixarcs')
            self.fix_r = self.ampl.getParameter('fix_r').to_dict()
            self.exdiam = self.ampl.getParameter('fixdiam').to_dict()
            for (i, j) in self.arcs:
                if (i,j) not in self.fixarcs:
                    # Original constraint value
                    lhs = h_values[i] - h_values[j]
                    alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                    #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                    original_rhs = q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs  
                    #original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs 
                    # Approximated constraint value
                    approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                    # approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

                    #approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs
                    con2_original_gap[f"{i},{j}"] = lhs - original_rhs
                    con2_approx_gap[f"{i},{j}"] = lhs - approx_rhs

                    total_absolute_constraint_violation += abs(lhs - approx_rhs)    
                    con2_original_violation += abs(lhs - original_rhs) 
                    con2_approx_violation += abs(lhs - approx_rhs) 

                    # Compute absolute violation
                    absolute_violation =  original_rhs - approx_rhs
                    absolute_violations[f"{i},{j}"] = absolute_violation
                    con2_absolute_constraint_violation += abs(absolute_violation)

                    # Compute relative violation between original_rhs and approx_rhs
                    relative_violation = (original_rhs - approx_rhs) / (original_rhs+1e-14)
                    relative_violations[f"{i},{j}"] = relative_violation
                    con2_relative_constraint_violation += abs(relative_violation)

            #print("con2_gap:", con2_gap)
            for (i, j) in self.fixarcs:
                # Original constraint value
                lhs = h_values[i] - h_values[j]
                alpha_rhs = 10.67 * self.L[i, j] / ((self.fix_r[i,j] ** 1.852) * ((self.exdiam[i,j]) ** 4.87))
                #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                original_rhs = q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * 10.67 * self.L[i,j]/(self.fix_r[i,j]**1.852 * self.exdiam[i,j]**4.87) 
                #original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs

                # Approximated constraint value
                approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * 10.67 * self.L[i,j]/(self.fix_r[i,j]**1.852 * self.exdiam[i,j]**4.87) 
                # approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * 10.67 * self.L[i,j]/(self.fix_r[i,j]**1.852 * self.exdiam[i,j]**4.87) 

                #approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs

                #con2_original_violation =  lhs - original_rhs
                con2_original_gap[f"{i},{j}"] = lhs - original_rhs
                con2_original_violation += abs(lhs - original_rhs) 

                #con2_approx_violation =  lhs - approx_rhs
                con2_approx_gap[f"{i},{j}"] = lhs - approx_rhs

                total_absolute_constraint_violation += abs(lhs - approx_rhs)    
                con2_approx_violation += abs(lhs - approx_rhs) 

                # Compute absolute violation
                absolute_violation =  original_rhs - approx_rhs
                absolute_violations[f"{i},{j}"] = absolute_violation
                con2_absolute_constraint_violation += abs(absolute_violation)

                # Compute relative violation between original_rhs and approx_rhs
                relative_violation = (original_rhs - approx_rhs) / (original_rhs+1e-14)
                relative_violations[f"{i},{j}"] = relative_violation
                con2_relative_constraint_violation += abs(relative_violation)
        else:
            for (i, j) in q_values.keys():
                # Original constraint value
                lhs = h_values[i] - h_values[j]
                alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                #alpha_rhs = sum(10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * ((self.d[k]) ** 4.87)) for k in self.pipes)
                original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs
                
                # Approximated constraint value
                approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q_values[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                # approx_rhs = (q_values[i, j] * ((q_values[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

                #approx_rhs = (q_values[i, j]**3 * ((q_values[i, j]**2 + 1e-12) ** 0.426)/(q_values[i,j]**2 + 0.426*1e-12))*alpha_rhs

                con2_original_gap[f"{i},{j}"] = lhs - original_rhs
                con2_original_violation += abs(lhs - original_rhs) 
                
                con2_approx_gap[f"{i},{j}"] = lhs - approx_rhs
                
                total_absolute_constraint_violation += abs(lhs - approx_rhs)    
                con2_approx_violation += abs(lhs - approx_rhs) 

                 # Compute absolute violation
                absolute_violation =  original_rhs - approx_rhs
                absolute_violations[f"{i},{j}"] = absolute_violation
                con2_absolute_constraint_violation += abs(absolute_violation)

                # Compute relative violation between original_rhs and approx_rhs
                relative_violation = (original_rhs - approx_rhs) / (original_rhs + 1e-14)
                relative_violations[f"{i},{j}"] = relative_violation
                con2_relative_constraint_violation += abs(relative_violation)
           
        #print("con2_gap:", con2_gap)
        
        con3_gap = {}
        if self.data_number==6:
            self.fixarcs = self.ampl.getSet('fixarcs')
            for (i,j) in self.arcs:
                if (i,j) not in self.fixarcs:
                    con3_rhs = self.L[i,j]
                    con3_lhs = sum(l_values[i,j,k] for k in self.pipes) 
                    con3_violation = con3_lhs - con3_rhs
                    con3_gap[f"{i},{j}"] = con3_violation 
                    total_absolute_constraint_violation += abs(con3_violation)
        else:
            for (i,j) in self.arcs:
                con3_rhs = self.L[i,j]
                con3_lhs = sum(l_values[i,j,k] for k in self.pipes) 
                con3_violation = con3_lhs - con3_rhs
                con3_gap[f"{i},{j}"] = con3_violation 
                total_absolute_constraint_violation += abs(con3_violation)
        #print("con3_gap:", con3_gap)

        con4_gap = {}
        for (i,j) in self.arcs:
            for k in self.pipes:
                #con4_rhs = self.L[i,j]
                #con4_lhs = l_values[i,j,k]
                con4_violation = max(0,l_values[i,j,k]-self.L[i,j])
                con4_gap[f"{i},{j},{k}"] = con4_violation 
                total_absolute_constraint_violation += abs(con4_violation)
        #print("con4_gap:", con4_gap)
        
        con5_gap = {}
        for j in self.source:
            con5_rhs = self.E[j]
            con5_lhs = h_values[j]
            con5_violation = con5_lhs - con5_rhs
            con5_gap[f"{j}"] = con5_violation 
            total_absolute_constraint_violation += abs(con5_violation)
        #print("con5_gap:", con5_gap)

        con6_gap = {}
        for j in self.nodes:
            if j not in self.source:
                #con6_rhs = self.E[j] + self.P[j]
                #con6_lhs = h_values[j]
                con6_violation = max(0, -h_values[j] + self.E[j] + self.P[j])
                con6_gap[f"{j}"] = con6_violation 
                total_absolute_constraint_violation += abs(con6_violation)
        #print("con6_gap:", con6_gap)

        # print("*******************************************************************************\n")
        # print("Constraints violation:")

        table_data = []
        for constraint, vio in con1_gap.items():
               table_data.append([constraint, f"{con1_gap[constraint]:.8f}"])
        for constraint, vio in con2_approx_gap.items():
               table_data.append([constraint, f"{con2_approx_gap[constraint]:.8f}"])
        for constraint, vio in con3_gap.items():
               table_data.append([constraint, f"{con3_gap[constraint]:.8f}"])
        for constraint, vio in con4_gap.items():
               table_data.append([constraint, f"{con4_gap[constraint]:.8f}"])
        for constraint, vio in con5_gap.items():
               table_data.append([constraint, f"{con5_gap[constraint]:.8f}"])
        for constraint, vio in con6_gap.items():
               table_data.append([constraint, f"{con6_gap[constraint]:.8f}"])

        #headers = ["Constraint ID", "Violation"]
        #print(tabulate(table_data, headers=headers, tablefmt="grid"))
        print("Sum of constraints violation:", total_absolute_constraint_violation)
        #print("*******************************************************************************\n")
        #table_data = []
        #for constraint, vio in con2_original_gap.items():
        #       table_data.append([constraint, f"{con2_original_gap[constraint]:.8f}",  f"{con2_approx_gap[constraint]:.8f}"])

        # print("*******************************************************************************\n")
        #print("Constraint 2 violations:\n")
        #headers = ["Constraint ID", "Original Con Violation", "Approx Con Violation"]
        #print(tabulate(table_data, headers=headers, tablefmt="grid")) 
        #print("\nSum of violation of original con2:", con2_original_violation) 
        #print("Sum of violation of approx con2:", con2_approx_violation)
        table_data = []
        for constraint, vio in relative_violations.items():
            i_str, j_str = constraint.split(',')
            i, j = int(i_str), int(j_str)
            table_data.append([constraint, q_values[i,j], f"{con2_original_gap[constraint]:.8f}",  f"{con2_approx_gap[constraint]:.8f}", f"{absolute_violations[constraint]:.8f}", f"{relative_violations[constraint]:.8f}"])
        # print("*******************************************************************************\n")
        # print("Absolute and relative violations between original and approximation constraint 2:\n")
        # headers = ["Constraint ID", "flow value", "Original Con Violation", "Approx Con Violation", "Absolute Violation", "Relative Violation(in %)"]
        # print(tabulate(table_data, headers=headers, tablefmt="grid"))
        print("Sum of violation of original headloss constraint:", con2_original_violation) 
        print("Sum of violation of approx headloss constraint:", con2_approx_violation)
        print("Con2 sum of absolute violation between original function and approximate function:", con2_absolute_constraint_violation)
        print("Con2 sum of relative violation between original function and approximate function:", con2_relative_constraint_violation)
        # Print total violations
        #print("\nTotal absolute constraint violation:", total_absolute_constraint_violation)
        #print("Total relative constraint violation:", total_relative_constraint_violation)

        # print("*******************************************************************************\n")

    def head_increase(self):
        improved = False
        min_head_node = []
        for j, val in self.h.items():
            print(f"{j}:",self.E[j] + self.P[j], val)
            if round(val,1) == self.E[j] + self.P[j]:
                min_head_node.append(j)
        print("\n*********************************************************************************************")
        print("Iteration :",self.dia_red_iteration + self.iteration -1, "\n")

        self.all_duals = {}
        for con_name, val in self.ampl.get_constraints():
            self.all_duals[con_name] = val.getValues()
        sorted_nodes = []
        dual_dict = self.all_duals["con7"].to_dict()
        print(dual_dict)
        sorted_duals = dict(sorted(dual_dict.items(), key=lambda kv: abs(kv[1]), reverse=True))
        sorted_nodes = list(sorted_duals.keys())
        sorted_nodes = [node for node in sorted_nodes if node not in self.visited_nodes]
        # print("sorted_arcs:", self.visited_arc)
        # print("\nsorted_arcs:", sorted_arcs)
        # if self.data_number == 6:
        #     sorted_nodes = [arc for arc in sorted_arcs if arc not in self.fixarcs]
        # sorted_nodes = [node for node in sorted_nodes if node not in min_head_node]
        print("sorted_arcs:", sorted_nodes)

        print("----------------------------------------------------------------------------------------")
        print(f"{'Node':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")

        for i in sorted_nodes:
            self.visited_nodes.append(i)
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
                # if (x,y) != (i,j):
                # if val >= 1e-3 and val <=1e-3:
                ampl.eval(f'let l[{x},{y},{k}] := {val};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number ==5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            for x, val in self.h.items():
                ampl.eval(f'let h[{x}] := {val};') 
            # self.update_initial_points_with_perturbation(ampl, self.l, self.q, self.h) 
            if i in min_head_node:
                ampl.eval(f"s.t. headloss{i}: h[{i}] >= {self.h[i]} + 1e-4;")
            else:
                ampl.eval(f"s.t. headloss{i}: h[{i}] <= {self.h[i]} - 1e-4;")

            ampl.option['solver'] = "ipopt" 
            # ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"
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
            total_cost = ampl.getObjective("total_cost").value()
            if ampl.solve_result == "solved":
                # self.plot_graph(self.super_source_out_arc, total_cost, 0, q, h, self.D, (0,0), l, self.C)
                if total_cost < self.current_cost:
                    print(f"{str((i)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(total_cost)):<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{ampl.solve_result:<14}{'Yes':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    self.current_cost = total_cost
                    improved = True
                    self.network_graph = self.generate_random_acyclic_from_solution(q1)
                    self.best_acyclic_flow = self.network_graph.copy()
                    self.indegree_2_or_more = [node for node, indeg in self.best_acyclic_flow.in_degree() if indeg >= 2]
                    # print("indegree_2_or_more:", self.indegree_2_or_more)
                    # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
                    #best_arc = (v,u)
                    self.l = l1 
                    self.q = q1
                    self.h = h1 
                    # ampl.eval("display l;")
                    # ampl.eval("display l.rc;")
                    if self.data_number==5:
                        self.q1 = ampl.getVariable('q1').getValues().to_dict()
                        self.q2 = ampl.getVariable('q2').getValues().to_dict()
                    #self.sorted_nodes = sorted(self.indegree_2_or_more, key=lambda node: self.D[node], reverse=True)
                    #print("\nvisited_nodes:", self.visited_nodes)
                    #if self.visited_nodes:
                    #    self.sorted_nodes = [item for item in self.sorted_nodes if item not in self.visited_nodes]
                    #print("sorted_nodes", self.sorted_nodes) 
                    # self.fix_arc_set = list(set(self.super_source_out_arc) | set(self.fix_arc_set))
                else:
                    print(f"{str((i)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(total_cost)):<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{ampl.solve_result:<14}{'No':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
            else:
                print(f"{str((i)):<10}" 
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{self.format_indian_number(round(total_cost)):<14}" 
                    f"{(str(round(solve_time, 2)) + 's'):<12}"
                    f"{ampl.solve_result:<14}{'No':<10}"
                    f"{round(time.time() - self.start_time, 2)}s ")
            if improved:
                self.dia_red_iteration = self.dia_red_iteration + 1
                self.head_increase()
                break
    def max_flow(self):
        improved = False

        print("\n*********************************************************************************************")
        print("Iteration :", self.max_flow_iteration + self.headloss_increase_iteration + self.dia_red_iteration + self.iteration -1, "\n")

        self.all_duals = {}
        for con_name, val in self.ampl.get_constraints():
            self.all_duals[con_name] = val.getValues()
        sorted_arcs = []
        dual_dict = self.all_duals["con2"].to_dict()
        sorted_duals = dict(sorted(dual_dict.items(), key=lambda kv: abs(kv[1]), reverse=True))
        sorted_arcs = list(sorted_duals.keys())
        sorted_arcs = [arc for arc in sorted_arcs if arc not in self.visited_arc]
        # print("sorted_arcs:", self.visited_arc)
        # print("\nsorted_arcs:", sorted_arcs)
        if self.data_number == 6:
            sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fixarcs]
        sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fix_arc_set]
        print("sorted_arcs:", sorted_arcs)

        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")

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
                # if (x,y) != (i,j):
                # if val >= 1e-3 and val <=1e-3:
                ampl.eval(f'let l[{x},{y},{k}] := {val};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number ==5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            for x, val in self.h.items():
                ampl.eval(f'let h[{x}] := {val};') 
            # self.update_initial_points_with_perturbation(ampl, self.l, self.q, self.h) 
 
            # ampl.eval(f"""subject to con3_l_:sum {{(i,j) in arcs}} sum {{k in pipes}} C[k] * l[i,j,k] <= {self.current_cost};""")
            if self.q[i,j]>= 0:
                # print(self.h[i] - self.h[j])
                ampl.eval(f"s.t. max_flow{i}_{j}: q[{i},{j}]<={self.q[i,j]};")
            else:
                # print(self.h[i] - self.h[j])
                ampl.eval(f"s.t. max_flow{i}_{j}: q[{i},{j}]>={self.q[i,j]};")

            ampl.option['solver'] = "ipopt" 
            # ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"
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
                    improved = True
                    self.network_graph = self.generate_random_acyclic_from_solution(q1)
                    self.best_acyclic_flow = self.network_graph.copy()
                    self.indegree_2_or_more = [node for node, indeg in self.best_acyclic_flow.in_degree() if indeg >= 2]
                    # print("indegree_2_or_more:", self.indegree_2_or_more)
                    # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
                    #best_arc = (v,u)
                    self.l = l1 
                    self.q = q1
                    self.h = h1 
                    # ampl.eval("display l;")
                    # ampl.eval("display l.rc;")
                    if self.data_number==5:
                        self.q1 = ampl.getVariable('q1').getValues().to_dict()
                        self.q2 = ampl.getVariable('q2').getValues().to_dict()
                    #self.sorted_nodes = sorted(self.indegree_2_or_more, key=lambda node: self.D[node], reverse=True)
                    #print("\nvisited_nodes:", self.visited_nodes)
                    #if self.visited_nodes:
                    #    self.sorted_nodes = [item for item in self.sorted_nodes if item not in self.visited_nodes]
                    #print("sorted_nodes", self.sorted_nodes) 
                    self.fix_arc_set = list(set(self.super_source_out_arc) | set(self.fix_arc_set))
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
            if improved:
                self.dia_red_iteration = self.dia_red_iteration + 1
                self.max_flow()
                break
 
    def headloss_increase(self):
        improved = False
        # print("\n*********************************************************************************************")
        print("Iteration :", self.headloss_increase_iteration + self.iteration)
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
        print("sorted_arcs:", sorted_arcs)
        # sorted_arcs = [sorted_arcs[0]] 
        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")

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
                ampl.eval(f'let l[{x},{y},{k}] := {val - self.tol};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number ==5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            for x, val in self.h.items():
                ampl.eval(f'let h[{x}] := {val};') 
            
            current_duals = {}
            for con_name, val in ampl.get_constraints():
                dual_values = val.get_values()
                current_duals[con_name] = dual_values

            # Initialize dual values for all constraints
            for con_name, dual_values in self.all_duals.items():
                if con_name in current_duals:
                    ampl.get_constraint(con_name).set_values(dual_values)
 
            # self.update_initial_points_with_perturbation(ampl, self.l, self.q, self.h) 
            
            # ampl.eval(f"""subject to con3_l_:sum {{(i,j) in arcs}} sum {{k in pipes}} C[k] * l[i,j,k] <= {self.current_cost};""")
            # if self.h[i] - self.h[j]>= 0:
            #     # print(self.h[i] - self.h[j])
            #     ampl.eval(f"s.t. headloss1{i}_{j}: 0 <= h[{i}]- h[{j}]<={self.h[i] - self.h[j]} - 1e-3;")
            # else:
            #     # print(self.h[i] - self.h[j])
            #     ampl.eval(f"s.t. headloss2{i}_{j}: {self.h[i] - self.h[j]} + 1e-3 <= h[{i}]- h[{j}] <= 0;")
            ampl.eval(f"s.t. headloss1{i}_{j}: abs(h[{i}]- h[{j}])>=abs({self.h[i] - self.h[j]}) + 1e-2;")
            # ampl.eval(f"subject to con3_l{i}_{j}: sum {{k in pipes}} C[k]*l[{i},{j},k] <= {sum(self.C[k] * self.l[i,j,k] for k in self.pipes)} - 1e-2;")

            ampl.option['solver'] = "ipopt" 
            # ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"
            ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes mu_init = {self.mu_init}"
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
                    self.ampl = ampl
                    self.current_cost = total_cost
                    improved = True
                    self.network_graph = self.generate_random_acyclic_from_solution(q1)
                    self.best_acyclic_flow = self.network_graph.copy()
                    self.indegree_2_or_more = [node for node, indeg in self.best_acyclic_flow.in_degree() if indeg >= 2]
                    # print("indegree_2_or_more:", self.indegree_2_or_more)
                    # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
                    #best_arc = (v,u)
                    self.l = l1 
                    self.q = q1
                    self.h = h1 
                    # ampl.eval("display l;")
                    # ampl.eval("display l.rc;")
                    if self.data_number==5:
                        self.q1 = ampl.getVariable('q1').getValues().to_dict()
                        self.q2 = ampl.getVariable('q2').getValues().to_dict()
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
            if improved:
                self.headloss_increase_iteration = self.headloss_increase_iteration + 1
                self.headloss_increase()
                break

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
        print("Iteration :",self.dia_red_iteration + self.headloss_increase_iteration + self.iteration)

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
        print("sorted_arcs:", sorted_arcs)
        
        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")

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
                ampl.eval(f'let l[{x},{y},{k}] := {val - self.tol};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number ==5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            for x, val in self.h.items():
                ampl.eval(f'let h[{x}] := {val};') 
            
            current_duals = {}
            for con_name, val in ampl.get_constraints():
                dual_values = val.get_values()
                current_duals[con_name] = dual_values

            # Initialize dual values for all constraints
            for con_name, dual_values in self.all_duals.items():
                if con_name in current_duals:
                    # Initialize dual values for each constraint
                    ampl.get_constraint(con_name).set_values(dual_values) 
 
            # self.update_initial_points_with_perturbation(ampl, self.l, self.q, self.h)
            # ampl.eval(f"""subject to con3_l_:sum {{(i,j) in arcs}} sum {{k in pipes}} C[k] * l[i,j,k] <= {self.current_cost};""")
            # ampl.eval(f"subject to con3_{i}_{j}: sum{{k in pipes: k <=  {arc_max_dia[i,j]-1}}} l[{i},{j},k] = L[{i},{j}];")
            # ampl.eval(f"subject to con3_l{i}_{j}_{arc_max_dia[i,j]}: l[{i},{j},{arc_max_dia[i,j]}] = 0;")
            # ampl.eval(f"subject to con3_l{i}_{j}: sum {{k in pipes}} C[k]*l[{i},{j},k] <= {sum(self.C[k] * self.l[i,j,k] for k in self.pipes)}*(1-0.2);")
            # ampl.eval(f"subject to con3_l{i}_{j}: sum {{k in pipes}} C[k]*l[{i},{j},k] <= {self.C[arc_max_dia[i,j]] * self.L[i,j]};")
            for k in self.pipes:
                if k>=arc_max_dia[i,j]:
                    ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = 0;")
            # if self.h[i] - self.h[j]>= 0:
            #     ampl.eval(f"s.t. headloss1{i}_{j}: h[{i}]- h[{j}]>={self.h[i] - self.h[j]} + 1e-0;")
            # else:
            #     ampl.eval(f"s.t. headloss2{i}_{j}: h[{i}]- h[{j}]<={self.h[i] - self.h[j]} - 1e-0;")

            ampl.option['solver'] = "ipopt" 

            ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes mu_init = {self.mu_init}"
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
            if improved:
                self.dia_red_iteration = self.dia_red_iteration + 1
                self.diameter_reduction()
                break

    def iterate_acyclic_flows(self):
        # self.network_graph = self.best_acyclic_flow.copy()
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        # self.plot_graph(self.super_source_out_arc, self.current_cost, 0, self.q, self.h, self.D, (0,0), self.l, self.C)

        # print("\n*********************************************************************************************")
        print("Iteration :",self.iteration)
        improved = False
        self.all_duals = {}
        for con_name, val in self.ampl.get_constraints():
            # Get dual values for each constraint
            dual_values = val.getValues()
            self.all_duals[con_name] = dual_values

        sorted_all_arcs = []
        for node in  self.indegree_2_or_more:
            # print("\nNode:",node)
            for (u,v) in self.network_graph.in_edges(node):
                if (u,v) in self.arcs:
                    sorted_all_arcs.append((u,v))
                else:
                    sorted_all_arcs.append((v,u))
        sorted_all_arcs = [arc for arc in sorted_all_arcs if arc not in self.fix_arc_set]
        sorted_all_arcs = [arc for arc in sorted_all_arcs if arc not in self.visited_arc_reverse]

        for con_name, dual_values in self.all_duals.items():
            if con_name == "con2":
                # print(dual_values, "\n")
                dual_dict = dual_values.to_dict()
                dual_dict = {idx: val for idx, val in dual_dict.items() if idx in sorted_all_arcs}
                # max_index, max_dual = max(dual_dict.items(), key=lambda kv: abs(kv[1]))
                # print(f"  Max abs dual = {max_dual} at index {max_index}")
        sorted_by_abs_dual = dict(sorted(dual_dict.items(), key=lambda kv: abs(kv[1]), reverse=True))
        print("sorted_arcs: ",list(sorted_by_abs_dual.keys()))

        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")

        for edge in sorted_by_abs_dual.keys():
            self.visited_arc_reverse.append(edge)
            (u,v) = edge
            # self.load_model()
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
                ampl.eval(f'let l[{x},{y},{k}] := {val - self.tol};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number ==5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            for x, val in self.h.items():
                ampl.eval(f'let h[{x}] := {val};') 
            
            current_duals = {}
            for con_name, val in ampl.get_constraints():
                dual_values = val.get_values()
                current_duals[con_name] = dual_values

            # Initialize dual values for all constraints
            for con_name, dual_values in self.all_duals.items():
                if con_name in current_duals:
                    # Initialize dual values for each constraint
                    ampl.get_constraint(con_name).set_values(dual_values)               
            
            #self.ampl.eval(f"set inarc := {{{inarc_set}}};")
            #self.ampl.eval(f"set indegree_node := {{{set(self.indegree_2_or_more)}}};")
            fix_arc_set = self.fix_leaf_arc_flow()
            # self.update_initial_points(self.l, self.q, self.h)
            # self.update_initial_points_with_perturbation(self.ampl, self.l, self.q, self.h) 
            if self.q[u,v] >= 0:
                ampl.eval(f"s.t. flow_direction1{u}_{v}: q[{u}, {v}]<=0;")
            else:
                ampl.eval(f"s.t. flow_direction1{u}_{v}: q[{u}, {v}]>=0;")
            # with self.suppress_output():
            ampl.option["solver"] = "ipopt"
            ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = no halt_on_ampl_error = yes mu_init = {self.mu_init}")   #max_iter = 1000
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
                l = ampl.getVariable('l').getValues().to_dict()
                q = ampl.getVariable('q').getValues().to_dict()
                h = ampl.getVariable('h').getValues().to_dict() 

                if self.total_cost < self.current_cost:
                    # self.visited_nodes.append(node)
                    #self.visited_arc.append((v,u))
                    print(f"{str((u, v)):<10}"
                      f"{self.format_indian_number(round(self.current_cost)):<14}"
                      f"{self.format_indian_number(round(self.total_cost)):<14}"
                      f"{(str(round(self.ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                      f"{self.solve_result:<14}{'Yes':<10}"
                      f"{round(time.time() - self.start_time, 2)}s")
                    #print("\n")
                    # self.plot_graph(self.super_source_out_arc, self.total_cost, 0, q, h, self.D, (0,0), l, self.C)
                    self.current_cost = self.total_cost
                    improved = True
                    self.ampl = ampl
                    self.network_graph = self.generate_random_acyclic_from_solution(q)
                    # self.best_acyclic_flow = self.network_graph.copy()
                    # self.indegree_2_or_more = [node for node, indeg in self.best_acyclic_flow.in_degree() if indeg >= 2]
                    # print("indegree_2_or_more:", self.indegree_2_or_more)
                    self.l = l 
                    self.q = q
                    self.h = h 
                    if self.data_number==5:
                        self.q1 = ampl.getVariable('q1').getValues().to_dict()
                        self.q2 = ampl.getVariable('q2').getValues().to_dict()
                    self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set)
                    print("----------------------------------------------------------------------------------------")
                else:
                    print(f"{str((u, v)):<10}"
                      f"{self.format_indian_number(round(self.current_cost)):<14}"
                      f"{self.format_indian_number(round(self.total_cost)):<14}"
                      f"{(str(round(ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                      f"{self.solve_result:<14}{'No':<10}"
                      f"{round(time.time() - self.start_time, 2)}s")
            else:
                # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
                print(f"{str((u, v)):<10}"
                  f"{self.format_indian_number(round(self.current_cost)):<14}"
                  f"{self.format_indian_number(round(self.total_cost)):<14}"
                  f"{(str(round(ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                  f"{self.solve_result:<14}{'No':<10}"
                  f"{round(time.time() - self.start_time, 2)}s")
                #print("\n")
            # print("----------------------------------------------------------------------------------------")
            if improved:
                self.iteration += 1
                self.iterate_acyclic_flows()
                break
        # print("----------------------------------------------------------------------------------------")

    # Function to suppress output
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

    def initial_solution(self):
        ampl = AMPL()
        ampl.reset()
        #if self.data_number==5:
        ampl.read("lp_rel_wdnmodel1.mod")
        ampl.read_data(self.data_file)
        ampl.option['solver'] = "highs" 
        ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = no bound_relax_factor=0 warm_start_init_point = no halt_on_ampl_error = yes "

        #with self.suppress_output():
        ampl.solve()
        q_lp = ampl.getVariable('q').getValues().to_dict()
        l_lp = ampl.getVariable('l').getValues().to_dict()
        h_lp = ampl.getVariable('h').getValues().to_dict()
        print("total cost from lp rel:", ampl.get_objective("total_cost").value()) 
        return l_lp, q_lp, h_lp

    def solve(self):
        self.ampl.option["solver"] = "ipopt"
        self.ampl.set_option("ipopt_options", f"outlev = 0 tol = 1e-9 bound_relax_factor=0  bound_push = {self.bound_push} bound_frac = {self.bound_frac} halt_on_ampl_error = yes halt_on_ampl_error = yes warm_start_init_point = no expect_infeasible_problem = no mu_init = {self.mu_init} warm_start_same_structure = no")   #max_iter = 1000
        # self.ampl.set_option("ipopt_options", f"outlev = 0 tol = 1e-9 bound_relax_factor=0  bound_push = 0.01 bound_frac = 0.01 halt_on_ampl_error = yes halt_on_ampl_error = yes warm_start_init_point = no expect_infeasible_problem = no")   #max_iter = 1000
        # self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = yes halt_on_ampl_error = yes "
        #self.ampl.set_option("ipopt_options", f"outlev = 0  bound_relax_factor=0 warm_start_init_point = no halt_on_ampl_error = yes")   #max_iter = 1000
        #self.ampl.set_option("ipopt_options", f"outlev = 0 warm_start_init_point = no ")   #max_iter = 1000
        self.ampl.option["presolve_eps"] = "6.82e-14"
        self.ampl.option['presolve'] = 1
        #min_demand = self.ampl.getParameter('D_min').getValues().to_list()[0]
        #max_demand = self.ampl.getParameter('D_max').getValues().to_list()[0]
        #max_flow = self.ampl.getParameter('Q_max').getValues().to_list()[0]
        #print("min_demand:", min_demand)
        #print("max_demand:", max_demand)
        #print("max_flow:", max_flow)
        #d_max = self.ampl.getParameter('d_max').getValues().to_list()[0]
        #d_min = self.ampl.getParameter('d_min').getValues().to_list()[0]
        #self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = {epsilon};")
        self.ampl.solve()
        self.solve_result = self.ampl.solve_result
        self.total_cost = self.ampl.get_objective("total_cost").value()
        #self.ampl.eval("display q;")
        # print("Objective:", self.total_cost)
        # print("solve_result: ",self.solve_result)
        solve_time = self.ampl.get_value('_solve_elapsed_time')
        self.solver_time += solve_time
        self.number_of_nlp += 1

    def solve1(self):
        with self.suppress_output():
            # """Solve the optimization problem."""
            ampl.option["solver"] = "ipopt"
            #print("bound_push:", self.bound_push)
            #print("bound_frac:", self.bound_frac)
            #self.bound_push , self.bound_frac = (0.001, 0.000001)
            ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = no halt_on_ampl_error = yes")   #max_iter = 1000
            # self.ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  tol = 1e-9 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = yes halt_on_ampl_error = yes")   #max_iter = 1000
            #self.ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  bound_relax_factor=0 warm_start_init_point = yes halt_on_ampl_error = yes")   #max_iter = 1000
            # self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = yes halt_on_ampl_error = yes "
            ampl.option["presolve_eps"] = "6.82e-14"
            ampl.option['presolve'] = 1
            
            #min_demand = self.ampl.getParameter('D_min').getValues().to_list()[0]
            #max_demand = self.ampl.getParameter('D_max').getValues().to_list()[0]
            #max_flow = self.ampl.getParameter('Q_max').getValues().to_list()[0]

            #print("min_demand:", min_demand)
            #print("max_demand:", max_demand)
            #print("max_flow:", max_flow)
            #d_max = self.ampl.getParameter('d_max').getValues().to_list()[0]
            #d_min = self.ampl.getParameter('d_min').getValues().to_list()[0]

            ampl.solve()
            #self.ampl.eval("display q;")
            self.solve_result = ampl.solve_result
            self.total_cost = ampl.get_objective("total_cost").value()
        #self.ampl.eval("display q;")
        # print("Objective:", self.total_cost)
        # print("solve_result: ",self.solve_result)
        # print("eps:", epsilon)
        solve_time = ampl.get_value('_solve_elapsed_time')
        self.solver_time += solve_time
        self.number_of_nlp += 1

    def run(self):
        """Main function to run the Heuristic Approach."""
        self.start_time = time.time()
        self.bound_push , self.bound_frac = (0.1, 0.1)
        self.mu_init = 0.1
        # print("NLP solve using:  smooth approximation 1, Epsilon selection using absolute error\n")
        # print("NLP solve using: smooth approximation 1, epsilon selection using relative error\n")
        # print("NLP solve using: smooth approximation 2, epsilon selection using absolute error\n")
        print("NLP Model: Smooth Approximatate Model 2, Epsilon Selection using Relative Error")
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
        # self.ampl.eval("display l;")
        # self.ampl.eval("display q;")
        # self.ampl.eval("display h;")
        # self.ampl.eval("display l.rc;")
        # self.ampl.eval("display l;")
        # self.ampl.eval("display l.rc;")
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
        self.visited_arc_reverse = []
        # self.plot_graph(fix_arc_set, self.total_cost, 0, self.q, self.h, self.D, (0,0), self.l, self.C)
        print("---------------------------Reverse Arc Direction Approach------------------------------------")
        self.iteration = 1
        self.iterate_acyclic_flows() 
        # print("\n----------------------------Diameter Reduction Approach--------------------------------------")
        # self.dia_red_iteration = 1
        # self.visited_arc = []
        # self.diameter_reduction()
        # print("\n--------------------------------Head increase Approach---------------------------------------")
        # self.dia_red_iteration = 1
        # self.visited_nodes = []
        # self.head_increase()
        print("\n-----------------------------Headloss increase Approach--------------------------------------")
        self.headloss_increase_iteration = 1
        self.visited_arc = []
        self.headloss_increase()
        print("\n----------------------------Diameter Reduction Approach--------------------------------------")
        self.dia_red_iteration = 1
        self.visited_arc = []
        self.diameter_reduction()
        # print("\n---------------------------Reverse Arc Direction Approach------------------------------------")
        # self.visited_arc_reverse = []
        # self.iteration = 1
        # self.iterate_acyclic_flows() 
        # print("\n-----------------------------------Max Flow Approach-----------------------------------------")
        # self.max_flow_iteration = 1
        # self.visited_arc = []
        # self.max_flow() 
        print("\n************************************Final Best Results*****************************************")
        print("Water Network:", self.data_list[self.data_number])
        # self.eps = self.ampl.get_variable('eps').get_values().to_dict()
        self.eps = self.ampl.getParameter('eps').getValues().to_dict()
        print(f"Final best objective: {self.current_cost}")
        #self.ampl.eval("display q;")
        print("Number of nlp problem solved:", self.number_of_nlp)
        print("Total number of iteration:", self.iteration + self.dia_red_iteration)
        self.constraint_violations(self.q, self.h, self.l, self.eps, "ipopt")
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
