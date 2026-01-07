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
import json
import plotly.graph_objects as go

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
    
    def export_solution(self, iteration=1):
        solution = {
            "iteration": iteration,
            "objective": self.current_cost,
            "q": {f"{i},{j}": v for (i,j), v in self.q.items()},
            # "q": self.q,
            "h": self.h,
            "l": {f"{i},{j},{k}": self.l[i,j,k] 
                  for (i,j) in self.arcs for k in self.pipes if self.l[i,j,k] > 1e-3},
            "D": self.D,
            "L": {f"{i},{j}": v for (i,j), v in self.L.items()},
            "C": self.C,
            "hmin": {f"{i}": self.E[i] if i == list(self.source)[0] else self.E[i] + self.P[i] for i in self.nodes},
            # "pipe_diameters": self.l,
            "nodes": list(self.nodes),
            "arcs": [(i,j) for (i,j) in self.arcs]
        }
    
        with open(f"../figure/json_file/solution_{self.data_number}.json", "w") as f:
            json.dump(solution, f, indent=2)

    def plot_interactive_wdn(self, solution_file, node_pos):
        """
        Interactive WDN plot (server-safe HTML):
        - Green = positive flow, Red = negative flow
        - Arrows show direction (boundary to boundary)
        - Hover on arc: flow + pipe diameters & lengths
        - Hover on node: head, min head, demand
        """
        # --------------------------------------------------
        # Load solution
        # --------------------------------------------------
        with open(solution_file) as f:
            sol = json.load(f)
        def parse_key(k):
            return tuple(map(int, k.replace("(", "").replace(")", "").split(",")))
        q = {parse_key(k): v for k, v in sol["q"].items()}
        L = {parse_key(k): v for k, v in sol.get("L", {}).items()}
        # Pipe info: l[i,j,k]
        pipe_info = {}
        for k_str, val in sol.get("l", {}).items():
            if val < 1e-6:
                continue
            i, j, d = map(int, k_str.replace("(", "").replace(")", "").split(","))
            pipe_info.setdefault((i, j), []).append((d, val))
        h = {int(k): v for k, v in sol["h"].items()}
        hmin = {int(k): v for k, v in sol["hmin"].items()}
        D = {int(k): v for k, v in sol["D"].items()}
        
        # --------------------------------------------------
        # Compute node radius in DATA coordinates
        # --------------------------------------------------
        node_marker_size=22
        fig_width=1800
        fig_height=1800
        xs = [p[0] for p in node_pos.values()]
        ys = [p[1] for p in node_pos.values()]

        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        
        data_width  = xmax - xmin
        data_height = ymax - ymin
        BASE_SIZE = 1800 
        if data_width >= data_height:
            fig_width  = BASE_SIZE
            fig_height = int(BASE_SIZE * data_height / data_width)
        else:
            fig_height = BASE_SIZE
            fig_width  = int(BASE_SIZE * data_width / data_height)

        data_range = max(xs) - min(xs)
        node_radius = (node_marker_size / 2) * (data_width / fig_width)

        # --------------------------------------------------
        # Edge containers
        # --------------------------------------------------
        edge_x_pos, edge_y_pos, edge_text_pos = [], [], []
        edge_x_neg, edge_y_neg, edge_text_neg = [], [], []
        click_x, click_y, click_text = [], [], []
        
        arrows = []
        # --------------------------------------------------
        # EDGES + ARROWS
        # --------------------------------------------------
        for (i, j), flow in q.items():
            if i not in node_pos or j not in node_pos:
                continue
            x0, y0 = node_pos[i]
            x1, y1 = node_pos[j]
            if flow >= 0:
                # Direction: i -> j
                xs, ys = x0, y0
                xe, ye = x1, y1
            else:
                # Direction: j -> i
                xs, ys = x1, y1
                xe, ye = x0, y0
            
            dx, dy = xe - xs, ye - ys
            dist = math.hypot(dx, dy)
            if dist < 1e-8:
                continue
            
            # Shorten to node boundaries (ALWAYS correct now)
            sx = xs + node_radius * dx / dist
            sy = ys + node_radius * dy / dist
            ex = xe - node_radius * dx / dist
            ey = ye - node_radius * dy / dist

            mx = 0.5 * (sx + ex)
            my = 0.5 * (sy + ey)

            # Pipe hover text
            pipes = pipe_info.get((i, j), [])
            if pipes:
                pipe_txt = "<br>".join([f"  {d} : {l:.2f}" for d, l in pipes])
            else:
                pipe_txt = "No pipe selected"
            hover_text = (
                f"<b>Arc {i} → {j}</b><br>"
                f"Flow: {flow:.5f}<br>"
                f"Length: {L.get((i,j),0):.2f}<br>"
                f"<b>Pipes:</b><br>{pipe_txt}"
            )
            click_x.append(mx)
            click_y.append(my)
            click_text.append(hover_text)   # same hover text as arc

            if flow >= 0:
                edge_x_pos += [sx, ex, None]
                edge_y_pos += [sy, ey, None]
                # edge_text_pos.append(hover_text)
                arrows.append(dict(
                ax=sx, ay=sy,
                x=ex, y=ey,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=3,
                arrowsize=1,
                arrowwidth=2,
                arrowcolor="green"
            ))
            else:
                edge_x_neg += [sx, ex, None]
                edge_y_neg += [sy, ey, None]
                # edge_text_neg.append(hover_text)
                arrows.append(dict(
                ax=sx, ay=sy,
                x=ex, y=ey,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=3,
                arrowsize=1,
                arrowwidth=2,
                arrowcolor="red"
            ))
        edge_pos = go.Scatter(
            x=edge_x_pos, y=edge_y_pos,
            mode="lines",
            line=dict(width=2, color="green"),
            hoverinfo="text",
            # hovertext=edge_text_pos,
            name="Positive flow"
        )
        edge_neg = go.Scatter(
            x=edge_x_neg, y=edge_y_neg,
            mode="lines",
            line=dict(width=2, color="red"),
            hoverinfo="text",
            # hovertext=edge_text_neg,
            name="Negative flow"
        )
        arc_click_trace = go.Scatter(
        x=click_x,
        y=click_y,
        mode="markers",
        marker=dict(size=0,color="lightgreen",symbol="circle",opacity=0.0),
        hoverinfo="text",
        hovertext=click_text,
        name="Arc info points",
        showlegend=False,
        )

        # --------------------------------------------------
        # NODES
        # --------------------------------------------------
        node_x, node_y, node_text, node_labels = [], [], [], []
        node_colors = []
        for n, (x, y) in node_pos.items():
            node_x.append(x)
            node_y.append(y)
            node_labels.append(str(n))
            node_text.append(
                f"<b>Node {n}</b><br>"
                f"Head: {h.get(n,0):.2f}<br>"
                f"Min Head: {hmin.get(n,0):.2f}<br>"
                f"Demand: {D.get(n,0):.5f}"
            )
            if n in list(self.source):
                node_colors.append("royalblue")
            else:
                node_colors.append("skyblue")

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode="markers+text",
            text=node_labels,
            textposition="middle center",
            hoverinfo="text",
            hovertext=node_text,
            marker=dict(size=22,color=node_colors,line=dict(width=1.5, color="black")),
            textfont=dict(size=12, color="black"),
            name="Nodes"
        )
        # --------------------------------------------------
        # FIGURE
        # --------------------------------------------------
        fig = go.Figure([edge_pos, edge_neg, arc_click_trace, node_trace])

        fig.update_layout(
            title=f"WDN Network = {data_list[(self.data_number)]}, Total Cost = {sol['objective']:.2f}",
            annotations=arrows,
            # width=fig_width,
            # height=fig_height,
            hovermode="closest",
            showlegend=True,
            xaxis=dict(visible=True,zeroline=False,    # constrain="domain"   # 🔴 KEY LINE
            ),
            yaxis=dict(visible=True,zeroline=False,scaleanchor="x",scaleratio=1),
            margin=dict(l=20, r=20, t=50, b=20),
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)"
        )

        # fig.update_layout(
        #     title=f"WDN Local Solution (Objective = {sol['objective']:.2f})",
        #     annotations=arrows,
        #     hovermode="closest",
        #     showlegend=True,
        #     xaxis=dict(visible=False),
        #     yaxis=dict(visible=False),
        #     margin=dict(l=20, r=20, t=50, b=20)
        # )
        # --------------------------------------------------
        # SAVE (SERVER SAFE)
        # --------------------------------------------------
        fig.write_html(f"../figure/json_file/wdn_interactive_solution{self.data_number+1}.html", 
                       auto_open=False)
        print("✔ Interactive WDN plot saved (HTML)")

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
        print("Iteration :",self.head_increase_iter + self.dia_red_iteration + self.iteration -1, "\n")

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
                ampl.eval(f"s.t. headloss{i}: h[{i}] >= {self.h[i]} + 1e-2;")
            else:
                ampl.eval(f"s.t. headloss{i}: h[{i}] <= {self.h[i]} - 1e-2;")

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
        print("Iteration :", self.max_flow_iteration + self.headloss_increase_iteration + self.dia_red_iteration + self.iteration, "\n")

        # self.sen_score = {}
        # for (i,j) in self.arcs:
        #     self.sen_score[i,j] = -1.852 * (self.h[i] - self.h[j])/(np.abs(self.q[i,j])**2.852)
        #     # print((i,j), self.sen_score[i,j])
        # # print("sen_score:", self.sen_score)
        # edge = min(self.sen_score, key=self.sen_score.get) 
        # print("minimum sen_score:",edge, self.sen_score[edge[0], edge[1]])

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
                ampl.eval(f"s.t. max_flow{i}_{j}: q[{i},{j}]<={self.q[i,j]} - 1e-4;")
            else:
                # print(self.h[i] - self.h[j])
                ampl.eval(f"s.t. max_flow{i}_{j}: q[{i},{j}]>={self.q[i,j]} + 1e-4;")
            # ampl.eval(f"s.t. max_flow{i}_{j}: abs(q[{i},{j}]-{self.q[i,j]}) >= 1e-3;")
            # ampl.eval(f"s.t. max_head: sum{{i in nodes}} h[i] <= {sum(self.h[x] for x in self.nodes)}-1e-2;")

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

    def iterate_two_arc_reversals(self):
        """
        Try reversing two incoming arcs (at the same node) at a time for nodes with indegree >= 3.
        Uses duals (con2) to rank candidate pairs by |dual1| + |dual2|.
        Warm-starts the AMPL model with previous l,q,h and duals. Ensures no directed cycle is created.
        If an improved feasible solution is found, updates self.* and recurses (one improvement per call).
        """
        import itertools
        import networkx as nx
        import time

        # gather candidate nodes (indegree >= 3)
        self.indegree_3_or_more = [n for n, deg in self.network_graph.in_degree() if deg >= 3]
        print("indegree_3_or_more:", self.indegree_3_or_more)
        print("Iteration :", self.iteration)

        # collect duals from the current stored AMPL (self.ampl) or recompute if needed
        # Build self.all_duals if not present or stale
        self.all_duals = {}
        try:
            # If you have a prior AMPL object in self.ampl, use it
            source_ampl = getattr(self, "ampl", None)
            if source_ampl is not None:
                for con_name, val in source_ampl.get_constraints():
                    # attempt both variants used historically
                    try:
                        dual_vals = val.getValues()
                    except Exception:
                        dual_vals = val.get_values()
                    self.all_duals[con_name] = dual_vals
        except Exception as e:
            # fallback: empty duals
            print("Warning: failed reading duals from self.ampl:", e)

        # We only use con2 (as in your earlier code) to build dual dict on arcs
        dual_dict = {}
        if "con2" in self.all_duals:
            try:
                # convert to python dict: prefer .to_dict() if available
                dv = self.all_duals["con2"]
                try:
                    dual_py = dv.to_dict()
                except Exception:
                    # if it's an AMPLValue with getValues or python mapping
                    try:
                        dual_py = dict(dv)
                    except Exception:
                        # last resort: iterate
                        dual_py = {}
                        for k in dv:
                            dual_py[k] = dv[k]
                dual_dict = {idx: val for idx, val in dual_py.items()}
            except Exception as e:
                print("Warning: couldn't parse con2 duals:", e)

        # Build incoming-arc lists for nodes with indegree >=3
        node_incoming = {}
        for node in self.indegree_3_or_more:
            arcs = list(self.network_graph.in_edges(node))
            # normalize orientation to (u,v) as stored in self.arcs if possible
            norm_arcs = []
            for (a, b) in arcs:
                if (a, b) in self.arcs:
                    norm_arcs.append((a, b))
                elif (b, a) in self.arcs:
                    norm_arcs.append((b, a))
                else:
                    # keep as-is
                    norm_arcs.append((a, b))
            # remove arcs that are fixed or already visited as pairs with the same orientation
            norm_arcs = [e for e in norm_arcs if e not in self.fix_arc_set]
            if norm_arcs:
                node_incoming[node] = norm_arcs

        # generate candidate pairs (two incoming arcs of the same node)
        candidate_pairs = []
        for node, arcs in node_incoming.items():
            if len(arcs) < 2:
                continue
            for e1, e2 in itertools.combinations(arcs, 2):
                # do not attempt same arc twice
                if e1 == e2:
                    continue
                # skip if pair already tried (either order)
                pair_key = tuple(sorted([e1, e2]))
                if hasattr(self, "visited_arc_reverse_pairs") and pair_key in self.visited_arc_reverse_pairs:
                    continue
                # compute score = |dual(e1)| + |dual(e2)|
                d1 = abs(dual_dict.get(e1, 0)) if dual_dict else 0
                d2 = abs(dual_dict.get(e2, 0)) if dual_dict else 0
                score = d1 + d2
                candidate_pairs.append((score, node, e1, e2))

        # sort candidate pairs by descending score
        candidate_pairs.sort(reverse=True, key=lambda x: x[0])
        print(f"Found {len(candidate_pairs)} candidate pairs (two-arc reversal)")

        # ensure visited pairs container exists
        if not hasattr(self, "visited_arc_reverse_pairs"):
            self.visited_arc_reverse_pairs = set()

        improved = False

        # header for printing results (same format as your existing output)
        print("----------------------------------------------------------------------------------------")
        print(f"{'ArcPair':<24}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")

        # iterate candidate pairs
        for score, node, e1, e2 in candidate_pairs:
            # mark pair as visited
            pair_key = tuple(sorted([e1, e2]))
            self.visited_arc_reverse_pairs.add(pair_key)

            # quick cycle check: simulate reversing both edges in a copy of graph and check DAG
            Gcopy = self.network_graph.copy()
            # remove original directed edges (take care if they don't exist in that orientation)
            try:
                if Gcopy.has_edge(*e1):
                    Gcopy.remove_edge(*e1)
                else:
                    # maybe stored reversed — try to remove reversed orientation
                    if Gcopy.has_edge(e1[1], e1[0]):
                        Gcopy.remove_edge(e1[1], e1[0])
            except Exception:
                pass
            try:
                if Gcopy.has_edge(*e2):
                    Gcopy.remove_edge(*e2)
                else:
                    if Gcopy.has_edge(e2[1], e2[0]):
                        Gcopy.remove_edge(e2[1], e2[0])
            except Exception:
                pass

            # add reversed versions
            Gcopy.add_edge(e1[1], e1[0])
            Gcopy.add_edge(e2[1], e2[0])

            # if reversing both creates a directed cycle, skip
            if not nx.is_directed_acyclic_graph(Gcopy):
                # optionally you can try to detect if only trivial cycles are made, but skip for safety
                # print("Skipping pair (would create cycle):", e1, e2)
                continue

            # Build a fresh AMPL instance for the test solve (same structure as your earlier code)
            ampl = AMPL()
            ampl.reset()
            if self.data_number == 5:
                ampl.read("newyork_model.mod")
            elif self.data_number == 6:
                ampl.read("blacksburg_model.mod")
            else:
                ampl.read("wdnmodel.mod")
            ampl.read_data(self.data_file)

            # set variables warm-start (l, q, h) - use existing self.l, self.q, self.h
            # we subtract small tol from l (like you did) to keep feasibility margin
            for (x, y, k), val in self.l.items():
                ampl.eval(f'let l[{x},{y},{k}] := {val - self.tol};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number == 5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            for x, val in self.h.items():
                ampl.eval(f'let h[{x}] := {val};')

            # initialize dual values in the new AMPL with self.all_duals if present
            current_duals = {}
            for con_name, val in ampl.get_constraints():
                try:
                    current_duals[con_name] = val.get_values()
                except Exception:
                    try:
                        current_duals[con_name] = val.getValues()
                    except Exception:
                        current_duals[con_name] = None

            for con_name, dual_vals in self.all_duals.items():
                if con_name in current_duals and dual_vals is not None:
                    try:
                        ampl.get_constraint(con_name).set_values(dual_vals)
                    except Exception:
                        # try alternate method if names differ; ignore failures
                        try:
                            ampl.get_constraint(con_name).setValues(dual_vals)
                        except Exception:
                            pass

            # add two flow-direction constraints for e1 and e2 (flip sign requirement)
            def add_direction_constraint(edge):
                (u, v) = edge
                try:
                    # if current q >= 0 then reversing means we require q[u,v] <= 0 (push to opposite sign)
                    if self.q.get((u, v), 0) >= 0:
                        ampl.eval(f"s.t. flow_dir_{u}_{v}: q[{u},{v}] <= 0;")
                    else:
                        ampl.eval(f"s.t. flow_dir_{u}_{v}: q[{u},{v}] >= 0;")
                except Exception:
                    # if q stored reversed, attempt with reversed indexing
                    if self.q.get((v, u), 0) >= 0:
                        ampl.eval(f"s.t. flow_dir_{v}_{u}: q[{v},{u}] <= 0;")
                    else:
                        ampl.eval(f"s.t. flow_dir_{v}_{u}: q[{v},{u}] >= 0;")

            add_direction_constraint(e1)
            add_direction_constraint(e2)

            # solver options (same as your function)
            ampl.option["solver"] = "ipopt"
            ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = no halt_on_ampl_error = yes mu_init = {self.mu_init}")
            ampl.option["presolve_eps"] = "6.82e-14"
            ampl.option['presolve'] = 1

            # solve quietly
            with self.suppress_output():
                try:
                    ampl.solve()
                except Exception as solve_err:
                    # if AMPL raised, treat as no improvement
                    self.solve_result = getattr(ampl, "solve_result", str(solve_err))
                    self.total_cost = None
                    solve_time = ampl.get_value('_solve_elapsed_time') if 'ampl' in locals() else 0.0
                    self.solver_time += solve_time
                    self.number_of_nlp += 1
                    print(f"{str((e1, e2)):<24}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{'---':<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{str(self.solve_result):<14}{'No':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    continue

            # record result and timing
            self.solve_result = ampl.solve_result
            try:
                self.total_cost = ampl.get_objective("total_cost").value()
            except Exception:
                # if objective not available, skip
                self.total_cost = None
            try:
                solve_time = ampl.get_value('_solve_elapsed_time')
            except Exception:
                solve_time = 0.0
            self.solver_time += solve_time
            self.number_of_nlp += 1

            # handle solved case
            if self.solve_result == "solved" and self.total_cost is not None:
                # get variable values
                try:
                    lnew = ampl.getVariable('l').getValues().to_dict()
                except Exception:
                    lnew = self.l.copy()
                try:
                    qnew = ampl.getVariable('q').getValues().to_dict()
                except Exception:
                    qnew = self.q.copy()
                try:
                    hnew = ampl.getVariable('h').getValues().to_dict()
                except Exception:
                    hnew = self.h.copy()

                if self.total_cost < self.current_cost:
                    # success: update everything and change graph orientation for both edges
                    print(f"{str((e1, e2)):<24}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{self.solve_result:<14}{'Yes':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")

                    # update current solution
                    self.current_cost = self.total_cost
                    improved = True
                    # keep the AMPL instance that solved
                    self.ampl = ampl

                    # update graph orientation: remove original edges and add reversed ones
                    # careful to use correct orientation existing in self.network_graph
                    if self.network_graph.has_edge(*e1):
                        self.network_graph.remove_edge(*e1)
                    elif self.network_graph.has_edge(e1[1], e1[0]):
                        self.network_graph.remove_edge(e1[1], e1[0])
                    self.network_graph.add_edge(e1[1], e1[0])

                    if self.network_graph.has_edge(*e2):
                        self.network_graph.remove_edge(*e2)
                    elif self.network_graph.has_edge(e2[1], e2[0]):
                        self.network_graph.remove_edge(e2[1], e2[0])
                    self.network_graph.add_edge(e2[1], e2[0])

                    # update stored solution vectors
                    self.l = lnew
                    self.q = qnew
                    self.h = hnew
                    # copy to backup fields as you did earlier
                    self.l1 = self.l
                    self.q1 = self.q
                    self.h1 = self.h
                    if self.data_number == 5:
                        self.q1 = ampl.getVariable('q1').getValues().to_dict()
                        self.q2 = ampl.getVariable('q2').getValues().to_dict()

                    # update fix_arc_set to include super_source_out_arc and any leaf-fixed arcs
                    fix_arc_set = self.fix_leaf_arc_flow()
                    self.fix_arc_set = list(set(self.super_source_out_arc) | set(fix_arc_set))

                    print("----------------------------------------------------------------------------------------")
                else:
                    # not improved
                    print(f"{str((e1, e2)):<24}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{self.solve_result:<14}{'No':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
            else:
                # not solved/infeasible
                print(f"{str((e1, e2)):<24}"
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{'---' if self.total_cost is None else self.format_indian_number(round(self.total_cost)):<14}"
                    f"{(str(round(solve_time, 2)) + 's'):<12}"
                    f"{self.solve_result:<14}{'No':<10}"
                    f"{round(time.time() - self.start_time, 2)}s")

            # if we improved, recurse to continue search from new graph/solution
            if improved:
                self.iteration += 1
                # recursively continue searching from the improved graph
                self.iterate_two_arc_reversals()
                break

        # finished trying all pairs (or broke after improvement)
        if not improved:
            print("No two-arc reversal improved the solution in this iteration.")

    def iterate_branching_on_q(self, eps=0, max_branch_per_iter=None):
        """
        Continuous branching heuristic on q_{i,j}.
        For each candidate arc (u,v) (ordered by absolute dual), solve two subproblems:
          A) q[u,v] <= q0 - eps
          B) q[u,v] >= q0 + eps
        If either produces an improved solution, accept it and move to the next variable.
        """
        print("Iteration :", self.cont_branching_iteration)
        improved_global = False
        # collect duals from current ampl (original model)
        self.all_duals = {}
        for con_name, val in self.ampl.get_constraints():
            dual_values = val.getValues()
            self.all_duals[con_name] = dual_values

        # build candidate arcs (same as your original logic)
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        sorted_all_arcs = []
        for node in self.indegree_2_or_more:
            for (u, v) in self.network_graph.in_edges(node):
                if (u, v) in self.arcs:
                    sorted_all_arcs.append((u, v))
                else:
                    sorted_all_arcs.append((v, u))

        # filter-out fixed / visited arcs
        sorted_all_arcs = [arc for arc in sorted_all_arcs if arc not in self.fix_arc_set]
        sorted_all_arcs = [arc for arc in sorted_all_arcs if arc not in self.visited_arc_reverse]

        # pick duals from con2 (as in your code) and sort by abs dual
        dual_dict = {}
        for con_name, dual_values in self.all_duals.items():
            if con_name == "con2":
                d = dual_values.to_dict()
                dual_dict = {idx: val for idx, val in d.items() if idx in sorted_all_arcs}
                break
        sorted_by_abs_dual = dict(sorted(dual_dict.items(), key=lambda kv: abs(kv[1]), reverse=True))
        cand_arcs = list(sorted_by_abs_dual.keys())

        if max_branch_per_iter:
            cand_arcs = cand_arcs[:max_branch_per_iter]

        print("sorted_arcs: ", cand_arcs)
        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")
        print("----------------------------------------------------------------------------------------")

        # For each candidate arc, attempt two branches
        for edge in cand_arcs:
            # mark as visited to avoid re-trying same arc next time
            self.visited_arc_reverse.append(edge)
            u, v = edge
            q0 = self.q.get((u, v), 0.0)  # current q value for this arc

            # prepare two branch constraints (A then B)
            branches = [
                ("le", q0 - eps, f"branch_le_{u}_{v}"),
                ("ge", q0 + eps, f"branch_ge_{u}_{v}")
            ]

            # try each branch; if improved, accept and move to next variable
            improved_for_edge = False
            for sense, rhs, cname in branches:
                # load fresh AMPL model and data
                ampl = AMPL()
                ampl.reset()
                if self.data_number == 5:
                    ampl.read("newyork_model.mod")
                elif self.data_number == 6:
                    ampl.read("blacksburg_model.mod")
                else:
                    ampl.read("wdnmodel.mod")
                ampl.read_data(self.data_file)

                # set initial continuous variables (l, q, h) from current solution (use slight perturbation for l as you did)
                for (x, y, k), val in self.l.items():
                    ampl.eval(f'let l[{x},{y},{k}] := {val - self.tol};')
                for (x, y), val in self.q.items():
                    ampl.eval(f'let q[{x},{y}] := {val};')
                    if self.data_number == 5:
                        ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                        ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
                for x, val in self.h.items():
                    ampl.eval(f'let h[{x}] := {val};')

                # initialize duals from parent
                # current_duals = {}
                # for con_name, val in ampl.get_constraints():
                #     current_duals[con_name] = val.getValues()
                # for con_name, dual_values in self.all_duals.items():
                #     if con_name in current_duals:
                #         ampl.get_constraint(con_name).set_values(dual_values)

                # add branch constraint for q[u,v]
                # ensure unique constraint name, and add as <= or >= accordingly
                if sense == "le":
                    ampl.eval(f"s.t. {cname}: q[{u},{v}] <= {rhs};")
                else:
                    ampl.eval(f"s.t. {cname}: q[{u},{v}] >= {rhs};")

                # add any leaf/fix arc constraints you maintain
                fix_arc_set = self.fix_leaf_arc_flow()

                # configure solver & warm-start
                ampl.option["solver"] = "ipopt"
                ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes")
                ampl.option["presolve_eps"] = "6.82e-14"
                ampl.option['presolve'] = 1

                # solve quietly and record results
                with self.suppress_output():
                    try:
                        ampl.solve()
                    except Exception as e:
                        # solver failure (e.g., infeasible or crash)
                        self.solve_result = getattr(ampl, "solve_result", "error")
                        self.total_cost = float('inf')
                        solve_time = ampl.get_value('_solve_elapsed_time') if ampl.get_value('_solve_elapsed_time') is not None else 0.0
                    else:
                        self.solve_result = ampl.solve_result
                        # if objective is unavailable for infeasible solves, guard it
                        try:
                            self.total_cost = ampl.get_objective("total_cost").value()
                        except:
                            self.total_cost = float('inf')
                        solve_time = ampl.get_value('_solve_elapsed_time')

                self.solver_time += solve_time if solve_time is not None else 0.0
                self.number_of_nlp += 1

                # print row
                new_cost_str = self.format_indian_number(round(self.total_cost)) if self.total_cost != float('inf') else "inf"
                print(f"{str((u, v)):<10}"
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{new_cost_str:<14}"
                    f"{(str(round(solve_time, 2)) + 's') if solve_time is not None else '':<12}"
                    f"{self.solve_result:<14}{'Yes' if self.total_cost < self.current_cost else 'No':<10}"
                    f"{round(time.time() - self.start_time, 2)}s")

                # if solved and improved, accept solution and update state
                if self.solve_result == "solved" and self.total_cost < self.current_cost:
                    # pull new variables
                    try:
                        l_new = ampl.getVariable('l').getValues().to_dict()
                        q_new = ampl.getVariable('q').getValues().to_dict()
                        h_new = ampl.getVariable('h').getValues().to_dict()
                    except Exception:
                        # fallback — if variable names differ or missing, treat as no-improvement
                        l_new, q_new, h_new = None, None, None

                    if q_new is not None:
                        # accept new solution
                        self.current_cost = self.total_cost
                        improved_global = True
                        improved_for_edge = True
                        self.cont_branching_iteration += 1

                        # update ampl instance and solution containers
                        self.ampl = ampl
                        self.l = l_new
                        self.q = q_new
                        self.h = h_new
                        if self.data_number == 5:
                            self.q1 = ampl.getVariable('q1').getValues().to_dict()
                            self.q2 = ampl.getVariable('q2').getValues().to_dict()

                        # update network graph to follow new acyclic solution
                        self.network_graph = self.generate_random_acyclic_from_solution(self.q)

                        # update fix arc set (include super source arcs and leaf-fixed arcs)
                        self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set)

                        print("Accepted improvement. Moving to next variable.")
                        print("----------------------------------------------------------------------------------------")

                        # break branch loop and move to next variable
                        break

                # else: branch gave no improvement; continue to next branch
                # cleanup ampl instance (let Python GC handle it); continue trying other branch

            # after trying both branches for this edge, continue to next edge
            # (if improved_for_edge True we already updated state and will proceed with next candidate arc)
            # Optionally, stop early if you want only one improvement per outer call:
            if improved_global:
                self.iterate_branching_on_q(eps=eps, max_branch_per_iter=max_branch_per_iter)
                break

            # after all candidates tried
            # if improved_global:
            # continue further iterations if desired; here I follow original pattern of recursion
            # but to avoid deep recursion we call this method iteratively
            # (you can change to recursive call if you prefer)
            # call next iteration to continue branching from new solution
            # NOTE: to avoid infinite loops you might want to limit iterations externally.
            # self.iterate_branching_on_q(eps=eps, max_branch_per_iter=max_branch_per_iter)
            # else:
            print("No improvement found in this pass.")

    def flow_change_in_cycle(self):
        import networkx as nx
        from amplpy import AMPL
        import time
        print("Iteration:",self.flow_change_in_cycle_iteration)
        improved = False
        start_time_total = time.time()
        # ----------------------------
        # 1. Build undirected graph for cycle detection
        # ----------------------------
        G_undir = nx.Graph()
        G_undir.add_nodes_from(self.network_graph.nodes())
        G_undir.add_edges_from(self.network_graph.edges())
        basic_cycles = nx.cycle_basis(G_undir)
        print("Basic Cycles:", basic_cycles)

        if self.data_number==5:
            delta = 0
        else:
            delta = 0
        # Ensure visited_cycles exists
        if not hasattr(self, "visited_cycles"):
            self.visited_cycles = set()
        #
        # Step 2: Normalize & filter out visited cycles
        normalized = [tuple(sorted(c)) for c in basic_cycles]
        remaining_cycles = [c for c, nc in zip(basic_cycles, normalized) if nc not in self.visited_cycles]
        print("All cycles found:", len(basic_cycles))
        print("Visited cycles:", len(self.visited_cycles))
        print("Remaining new cycles:", len(remaining_cycles))


        print("-" * 90)
        print(f"{'Cycle':<14}{'C_Best_Sol':<14}{'New_Sol':<14}"
            f"{'Solve_Time':<12}{'Result':<14}{'Improved':<10}{'Time':<12}")
        print("-" * 90)

        # ----------------------------
        # 2. Iterate all simple fundamental cycles
        # ----------------------------
        # for cycle in basic_cycles:
        for idx, cycle in enumerate(remaining_cycles, start=1):
            norm_cycle = tuple(sorted(cycle))
            # print(f"\nProcessing New Cycle-{idx}: {norm_cycle}")
            # Store visited cycle
            self.visited_cycles.add(norm_cycle)
            # ---- Normalize cycle to avoid duplicates ----
            # norm_cycle = tuple(sorted(cycle))
            #
            # # ---- Skip if already used ----
            # if norm_cycle in self.visited_cycles:
            #     print(f"Skipping Cycle-{idx} (already used)")
            #     continue
            # Build arc sequence of the cycle
            arcs = [(cycle[k], cycle[(k+1) % len(cycle)]) for k in range(len(cycle))]
            cycle_name = str(arcs)
            t0 = time.time()
            # ----------------------------
            # 3. Load AMPL fresh model
            # ----------------------------
            ampl = AMPL()
            ampl.reset()
            if self.data_number == 5:
                ampl.read("newyork_model.mod")
            elif self.data_number == 6:
                ampl.read("blacksburg_model.mod")
            else:
                ampl.read("wdnmodel.mod")
            ampl.read_data(self.data_file)
            # Set current q as starting point

            # for (x, y, k), val in self.l.items():
            #     ampl.eval(f'let l[{x},{y},{k}] := {val};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number ==5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            # for x, val in self.h.items():
            #     ampl.eval(f'let h[{x}] := {val};') 

            # for (i, j), val in self.q.items():
            # if (i,j) in arcs:
            #     continue
            # else:
            # ampl.eval(f"let q[{i},{j}] := {val};")
            # ----------------------------
            # 4. Apply flow perturbation on cycle
            # ----------------------------
            for (i, j) in arcs:
                if (i, j) in self.arcs:
                    if self.data_number == 5:
                        ampl.eval(f"let q1[{i},{j}] := {self.q1[i,j] - delta};")
                        ampl.eval(f"let q2[{i},{j}] := {self.q2[i,j] - delta};")
                    else:
                        ampl.eval(f"let q[{i},{j}] := {self.q[i,j] - delta};")
                else:
                    if self.data_number == 5:
                        # reverse arc direction
                        ampl.eval(f"let q1[{j},{i}] := {self.q1[j,i] + delta};")
                        ampl.eval(f"let q2[{j},{i}] := {self.q2[j,i] + delta};")
                    else:
                        ampl.eval(f"let q[{j},{i}] := {self.q[j,i] + delta};")

            # for (i, j) in arcs:
            #     if (i, j) in self.arcs:
            #         ampl.eval(f"s.t. cycle_flow_{i}_{j}: q[{i},{j}] >= {self.q[i,j] - delta};")
            #     else:
            #         # reverse arc direction
            #         ampl.eval(f"s.t. cycle_flow_{j}_{i}: q[{j},{i}] <= {self.q[j,i] + delta};")

            # ----------------------------
            # 5. Solve
            # ----------------------------
            ampl.option["solver"] = "ipopt"
            # ampl.set_option("ipopt_options", "outlev=0 tol=1e-9 warm_start_init_point=yes")
            ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"
            # ampl.option["presolve_eps"]= "7.19e-10"
            try:
                with self.suppress_output():
                    ampl.solve()
                solve_time = ampl.get_value('_solve_elapsed_time')
                self.solver_time += solve_time
                self.number_of_nlp += 1
            except:
                solve_time = time.time() - t0
                print(f"{idx:<14}"
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{'-':<14}"
                    f"{str(round(solve_time, 2))+'s':<12}"
                    f"{'failed':<14}{'No':<10}"
                    f"{round(time.time() - start_time_total, 2)}s")
                continue
            solve_time = time.time() - t0
            result = ampl.solve_result
            new_cost = ampl.get_objective("total_cost").value()
            # ----------------------------
            # 6. Improvement check
            # ----------------------------
            if result == "solved" and new_cost < self.current_cost:
                print(f"{idx:<14}"
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{self.format_indian_number(round(new_cost)):<14}"
                    f"{str(round(solve_time, 2))+'s':<12}"
                    f"{result:<14}{'Yes':<10}"
                    f"{round(time.time() - start_time_total, 2)}s")
                # Update internal state
                self.current_cost = new_cost
                self.l = ampl.getVariable('l').getValues().to_dict()
                self.q = ampl.getVariable('q').getValues().to_dict()
                if self.data_number == 5:
                    self.q1 = ampl.getVariable('q1').getValues().to_dict()
                    self.q2 = ampl.getVariable('q2').getValues().to_dict()
                self.h = ampl.getVariable('h').getValues().to_dict()
                # Update graph based on new flow
                self.network_graph = self.generate_random_acyclic_from_solution(self.q)
                self.indegree_2_or_more = [n for n, indeg in self.network_graph.in_degree() if indeg >= 2]
                self.best_acyclic_flow = self.network_graph.copy()
                # Get updated h and l
                self.h = ampl.getVariable('h').getValues().to_dict()
                self.l = ampl.getVariable('l').getValues().to_dict()
                improved = True
                # self.headloss_increase_iteration += 1
                # self.headloss_increase()
                # break
                print("-" * 90)
            else:
                print(f"{idx:<14}"
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{self.format_indian_number(round(new_cost)):<14}"
                    f"{str(round(solve_time, 2))+'s':<12}"
                    f"{result:<14}{'No':<10}"
                    f"{round(time.time() - start_time_total, 2)}s")
            if improved:
                self.flow_change_in_cycle_iteration += 1
                self.flow_change_in_cycle()
                break

    def headloss_increase(self):
        improved = False
        # print("\n*********************************************************************************************")
        print("Iteration :", self.headloss_increase_iteration)

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
        print("sorted_arcs:", sorted_arcs)
        # sorted_arcs = [sorted_arcs[0]] 
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
            #     ampl.eval(f"s.t. headloss1{i}_{j}: h[{i}]- h[{j}]<={self.h[i] - self.h[j]} - 1e-3;")
            # else:
            #     # print(self.h[i] - self.h[j])
            #     ampl.eval(f"s.t. headloss2{i}_{j}: {self.h[i] - self.h[j]} + 1e-3 <= h[{i}]- h[{j}];")
            ampl.eval(f"s.t. headloss1{i}_{j}: abs(h[{i}]- h[{j}])>=abs({self.h[i] - self.h[j]}) + 1e-2;")
            # ampl.eval(f"s.t. max_head: sum{{i in nodes}} h[i] <= {sum(self.h[x] for x in self.nodes)}-1e-2;")
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
            self.total_cost = ampl.getObjective("total_cost").value()
            if ampl.solve_result == "solved":
                # self.plot_graph(self.super_source_out_arc, total_cost, 0, q, h, self.D, (0,0), l, self.C)
                if self.total_cost < self.current_cost:
                    print(f"{str((i,j)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{ampl.solve_result:<14}{'Yes':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    self.ampl = ampl
                    self.current_cost = self.total_cost
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
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(solve_time, 2)) + 's'):<12}"
                        f"{ampl.solve_result:<14}{'No':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
            else:
                print(f"{str((i,j)):<10}" 
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{self.format_indian_number(round(self.total_cost)):<14}" 
                    f"{(str(round(solve_time, 2)) + 's'):<12}"
                    f"{ampl.solve_result:<14}{'No':<10}"
                    f"{round(time.time() - self.start_time, 2)}s ")

            if improved:
                self.headloss_increase_iteration = self.headloss_increase_iteration + 1
                self.headloss_increase()
                break


    def diameter_reduction1(self):
        """
        Iteratively try to reduce diameters on promising arcs.
        For each candidate arc (sorted by |dual_con2|), build a reduced NLP:
          - allow flow change only on arcs in fundamental loops that include (i,j)
          - fix q on all other arcs to current solution
          - disable larger diameters on (i,j)
        Solve the reduced NLP; accept first improving solution and recurse.
        """
        improved = False
        arc_max_dia = {}

        # ---------------------------------------------------------
        #  Build arc_max_dia: max commercial diameter used on each arc
        # ---------------------------------------------------------
        if getattr(self, "data_number", None) == 6:
            try:
                self.fixarcs = self.ampl.getSet('fixarcs')
            except Exception:
                self.fixarcs = set()

            for (i, j, d), val in self.l.items():
                if (i, j) not in self.fixarcs and (j, i) not in self.fixarcs:
                    if val > 1e-3:
                        arc_max_dia.setdefault((i, j), d)
                        arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)
        else:
            for (i, j, d), val in self.l.items():
                if val > 1e-3:
                    arc_max_dia.setdefault((i, j), d)
                    arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)

        print("Iteration :", getattr(self, "dia_red_iteration", 0))

        # ---------------------------------------------------------
        # Compute fundamental loops using cycle_basis
        # ---------------------------------------------------------
        G_undir = nx.Graph()
        G_undir.add_nodes_from(self.network_graph.nodes())
        G_undir.add_edges_from(self.network_graph.edges())

        basic_cycles = nx.cycle_basis(G_undir)
        # print("Basic Cycles:", basic_cycles)

        self.loops = []
        self.loop_sign = []

        for cycle in basic_cycles:
            loop_arcs = []
            signs = {}

            for k in range(len(cycle)):
                u = cycle[k]
                v = cycle[(k + 1) % len(cycle)]

                if (u, v) in self.arcs:
                    loop_arcs.append((u, v))
                    signs[(u, v)] = +1
                elif (v, u) in self.arcs:
                    loop_arcs.append((v, u))
                    signs[(v, u)] = -1
                else:
                    raise ValueError(f"Arc ({u},{v}) or ({v},{u}) not found in directed graph")

            self.loops.append(sorted(loop_arcs))
            self.loop_sign.append(signs)

        print(f"Found {len(self.loops)} fundamental loops using cycle_basis().") 

        source = list(self.source)[0]
        # Build BFS spanning tree from source
        T = nx.bfs_tree(G_undir, source)

        # Build path_sign dictionary: key = (j,i,k)
        self.paths = {}
        path_sign = {}
        for j in self.nodes:
            if j == source:
                continue

            try:
                path_nodes = nx.shortest_path(T, source, j)
            except:
                raise RuntimeError(f"No path found from source {source} to {j}")

            # convert nodes → arcs in the path
            path_arcs = [(path_nodes[t], path_nodes[t+1])
                for t in range(len(path_nodes)-1)]

            path = []
            for (u,v) in path_arcs:
                if (u,v) in self.arcs:
                    path.append((u,v))
                else:
                    path.append((v,u))
            self.paths[j] = path

            # assign signs
            for (i,k) in self.arcs:
                if (i,k) in path_arcs:
                    path_sign[(j,i,k)] = 1
                elif (k,i) in path_arcs:
                    path_sign[(j,i,k)] = -1
                else:
                    path_sign[(j,i,k)] = 0

        # ---------------------------------------------------------
        # Read duals, sort arcs by |dual|
        # ---------------------------------------------------------
        self.all_duals = {}

        try:
            for con_name, val in self.ampl.get_constraints():
                self.all_duals[con_name] = val.getValues()
        except Exception:
            self.all_duals = {}

        sorted_arcs = []
        try:
            dual_dict = self.all_duals["con2"].to_dict()
            sorted_duals = dict(sorted(dual_dict.items(),
                                       key=lambda kv: abs(kv[1]),
                                       reverse=True))
            sorted_arcs = list(sorted_duals.keys())
        except Exception:
            sorted_arcs = list(self.arcs)

        # Filter invalid arcs
        sorted_arcs = [arc for arc in sorted_arcs
            if arc not in getattr(self, "visited_arc", [])]

        if getattr(self, "data_number", None) == 6:
            sorted_arcs = [arc for arc in sorted_arcs
                if arc not in getattr(self, "fixarcs", set())]

        sorted_arcs = [arc for arc in sorted_arcs
            if arc not in getattr(self, "fix_arc_set", set())]

        sorted_arcs = [arc for arc in sorted_arcs
            if arc in arc_max_dia and arc_max_dia[arc] != 1]

        print("sorted_arcs:", sorted_arcs)

        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"
            f"{'Solve_Time':<12}{'Solve_Result':<14}"
            f"{'Improved':<10}{'Time':<12}")
        print("----------------------------------------------------------------------------------------")

        if not hasattr(self, "loops") or not hasattr(self, "loop_sign"):
            raise RuntimeError("Fundamental loops not available")

        # ---------------------------------------------------------
        # Try each arc for diameter reduction
        # ---------------------------------------------------------
        for (i, j) in sorted_arcs:

            self.visited_arc.append((i, j))

            # Find all loop arcs affecting this arc
            affected_arcs = set()
            for loop_arcs in self.loops:
                if (i, j) in loop_arcs or (j, i) in loop_arcs:
                    affected_arcs.update(loop_arcs)

            if not affected_arcs:
                affected_arcs = {(i, j)}

            # ------------------------------------------------------
            # Build reduced AMPL model
            # ------------------------------------------------------
            ampl = AMPL()
            ampl.reset()

            if getattr(self, "data_number", None) == 5:
                ampl.read("newyork_model.mod")
            elif getattr(self, "data_number", None) == 6:
                ampl.read("blacksburg_model.mod")
            else:
                ampl.read("loop_wdnmodel.mod")

            ampl.read_data(self.data_file)

            # Warm-start q, l, h
            for (x, y, k), val in self.l.items():
                v = float(val) if val is not None else 0.0
                try:
                    ampl.eval(f"let l[{x},{y},{k}] := {v};")
                except Exception:
                    pass

            for (x, y), val in self.q.items():
                try:
                    ampl.eval(f"let q[{x},{y}] := {float(val)};")
                except Exception:
                    pass

            for nd, val in self.h.items():
                try:
                    ampl.eval(f"let h[{nd}] := {float(val)};")
                except Exception:
                    pass

            # Fix flows outside affected arcs
            # for (u, v) in self.arcs:
            #     if (u, v) not in affected_arcs:
            #         q_fix = float(self.q.get((u, v), 0.0))
            #         try:
            #             ampl.eval(f"subject to fix_q_{u}_{v}: q[{u},{v}] = {q_fix};")
            #             for k in self.pipes:
            #                 ampl.eval(f"subject to fix_l_{u}_{v}_{k}: l[{u},{v},{k}] = {self.l[u,v,k]};")
            #         except Exception:
            #             pass

            # Disable larger diameters on (i,j)
            max_d = arc_max_dia.get((i, j), None)
            if max_d is None:
                max_d = max(self.pipes) if self.pipes else None

            if max_d is not None:
                for k in self.pipes:
                    if k >= max_d:
                        try:
                            ampl.eval(f"subject to disable_dia_{i}_{j}_{k}: "
                                f"l[{i},{j},{k}] = 0;")
                        except Exception:
                            pass

            # number of loops
            nL = len(self.loops)
            ampl.eval(f"set Loops := 1..{nL};")

            # declare LOOP set and sign param
            # ampl.eval("set LOOP{Loops} within arcs;")
            ampl.eval("param sign{arcs, Loops} default 0;")
            
            # populate LOOP and sign
            for l_idx, arc_list in enumerate(self.loops, start=1):
                # print(l_idx, arc_list)
                # 1. assign arcs of this loop
                ampl_tuples = [tuple((u, v)) for (u, v) in arc_list]
                # print(ampl_tuples)
                ampl.eval(f'set LOOP{l_idx} within {{i in nodes, j in nodes: i != j}};')
                ampl.set[f'LOOP{l_idx}'] = ampl_tuples
                # loops_set = ampl.get_set(f"LOOP{l_idx}")
                # print(loops_set)
                # for (i,j) in loops_set:
                #     print(i,j)

                for (u, v), s in self.loop_sign[l_idx - 1].items():
                    ampl.param["sign"][u, v, l_idx] = s

                # ---------------------------------------------------------------
                # ADD CYCLE BASIS LOOP HEAD LOSS CONSTRAINTS IN AMPL
                # ---------------------------------------------------------------
                ampl.eval(
                    f"s.t. headloss{l_idx}: "
                        f"sum{{(i,j) in LOOP{l_idx}}} sign[i,j,{l_idx}] * q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426 * eps[i,j]^2) * (sum{{k in pipes}} 10.67 * l[i,j,k] / (R[k]^1.852 * d[k]^4.87)) = 0;"
                )
                # c = ampl.getConstraint(f"headloss{l_idx}")
                # print(c)
                # ampl.eval(f"expand con3;")
                # ampl.eval(f"expand headloss{l_idx};")

            ampl.eval("param path_sign{nodes, arcs} default 0;")
            ps = ampl.param["path_sign"]

            for (u,v,k), val in path_sign.items():
                ps[u, v, k] = val
            # print(list(self.paths.values()))
            for key, arc_list in self.paths.items():
                # print(idx, arc_list)
                ampl_tuples = [tuple((u,v)) for (u,v) in arc_list]
                ampl.eval(f'set PATH{key} within {{i in nodes, j in nodes: i != j}};')
                ampl.set[f'PATH{key}'] = ampl_tuples

                # ---------------------------------------------------------------
                # ADD HEAD BOUNDS CONSTRAINTS IN AMPL
                # ---------------------------------------------------------------
                ampl.eval(f"""
                subject to head_eqn{key}:
                    h[{key}] = 
                        E[{source}] 
                        - sum{{(i,k) in PATH{key}}}
                             path_sign[{key},i,k] *
                             ( q[i,k]^3 * (q[i,k]^2 + eps[i,k]^2)^0.426 )
                             / ( q[i,k]^2 + 0.426*eps[i,k]^2 )
                             * ( sum{{p in pipes}} 10.67*l[i,k,p]/(R[p]^1.852 * d[p]^4.87) );
                """)

            ampl.option['solver'] = "ipopt" 

            ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"
            #ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no  bound_relax_factor=0 warm_start_init_point = yes halt_on_ampl_error = yes"
            #with self.suppress_output():
            ampl.option["presolve_eps"]= "2.82e-9"

            with self.suppress_output():
                ampl.solve()
            solve_time = ampl.get_value('_solve_elapsed_time')
            self.solver_time += solve_time
            self.number_of_nlp += 1

            l1 = ampl.getVariable('l').getValues().to_dict()
            q1 = ampl.getVariable('q').getValues().to_dict()
            h1 = ampl.getVariable('h').getValues().to_dict()
            # ampl.eval("display l;")
            # ampl.eval("display q;")
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
        
        # sorted_arcs = [sorted_arcs[0]]
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
            # for (u,v) in self.arcs:
            #     if (u,v) != (i,j):
            #         # ampl.eval(f"s.t. arc_dia{u}_{v}: sum{{k in pipes: k <=  {arc_max_dia[u,v]+1}}} l[{u},{v},k] = L[{u},{v}];")
            #         for k in self.pipes:
            #             if k>=arc_max_dia[u,v]+2:
            #                 ampl.eval(f"subject to con3__{u}_{v}_{k}: l[{u},{v},{k}] = 0;")
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
            ampl.option["presolve_eps"]= "7.19e-13"

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


    def solve_nlp_model(self, arc_max_dia, i,j):
        ampl = AMPL()
        ampl.reset()
        if self.data_number==5:
            ampl.read("newyork_model.mod")
        elif self.data_number==6:
            ampl.read("blacksburg_model.mod")
        else:
            if self.lp:
                ampl.read("lp_model2.mod")
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
 
        for k in self.pipes:
            if k>=arc_max_dia[i,j]:
                ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = 0;")
        
        if self.lp:
            for (u,v) in self.arcs:
                # ampl.eval(f"s.t. arc_flow{u}_{v}: q[{u}, {v}] = {self.q[u,v]};")
                ampl.param["q"][u, v] = self.q[u,v]

        if self.lp:
            ampl.option['solver'] = "cplexamp" 
        else:
            ampl.option['solver'] = "ipopt" 
        ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"
        #ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no  bound_relax_factor=0 warm_start_init_point = yes halt_on_ampl_error = yes"
        #with self.suppress_output():
        ampl.option["presolve_eps"]= "1.52e-09"

        with self.suppress_output():
            ampl.solve()
        # solve_time = ampl.get_value('_solve_elapsed_time')
        # self.solver_time += solve_time
        # self.number_of_nlp += 1
        return ampl

    def arc_dia_red(self):
        improved = False
        arc_max_dia = {}
        if self.data_number == 6:
            self.fixarcs = self.ampl.getSet('fixarcs')
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
        
        self.all_duals = {}
        for con_name, val in self.ampl.get_constraints():
            self.all_duals[con_name] = val.getValues()
        sorted_arcs = []
        dual_dict = self.all_duals["con2"].to_dict()
        sorted_duals = dict(sorted(dual_dict.items(), key=lambda kv: abs(kv[1]), reverse=True))
        sorted_arcs = list(sorted_duals.keys())
        sorted_arcs = [arc for arc in sorted_arcs if arc not in self.visited_arc]
        if self.data_number == 6:
            sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fixarcs]
        sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fix_arc_set]
        sorted_arcs = [arc for arc in sorted_arcs if arc_max_dia[arc[0], arc[1]] != 1]
        
        # print("sorted_arcs:", sorted_arcs)
        
        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")
        print("----------------------------------------------------------------------------------------")
        # arc_list = {}
        for (i,j) in sorted_arcs:
            self.visited_arc.append((i,j))
            self.lp = True
            ampl = self.solve_nlp_model(arc_max_dia, i, j)
            # l1 = ampl.getVariable('l').getValues().to_dict()
            # q1 = ampl.getVariable('q').getValues().to_dict()
            # h1 = ampl.getVariable('h').getValues().to_dict()
            total_cost = ampl.getObjective("total_cost").value()
 
            if ampl.solve_result == "solved":
                if abs(self.current_cost - total_cost) / max(self.current_cost, total_cost) < 0.05:
                    # arc_list[i,j] = total_cost
                    self.lp = False
                    ampl = self.solve_nlp_model(arc_max_dia, i, j)
                    solve_time = ampl.get_value('_solve_elapsed_time')
                    self.solver_time += solve_time
                    self.number_of_nlp += 1

                    l1 = ampl.getVariable('l').getValues().to_dict()
                    q1 = ampl.getVariable('q').getValues().to_dict()
                    h1 = ampl.getVariable('h').getValues().to_dict()
                    total_cost = ampl.getObjective("total_cost").value()
                    if ampl.solve_result == "solved":
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
                            self.l = l1 
                            self.q = q1
                            self.h = h1 
                            if self.data_number==5:
                                self.q1 = ampl.getVariable('q1').getValues().to_dict()
                                self.q2 = ampl.getVariable('q2').getValues().to_dict()
                            print("----------------------------------------------------------------------------------------")
                        else:
                            print(f"{str((i,j)):<10}"
                                f"{self.format_indian_number(round(self.current_cost)):<14}"
                                f"{self.format_indian_number(round(total_cost)):<14}"
                                f"{(str(round(solve_time, 2)) + 's'):<12}"
                                f"{ampl.solve_result:<14}{'No':<10}"
                                f"{round(time.time() - self.start_time, 2)}s")
                    else:
                        print(f"{str((i,j)):<10}" 
                            f"{self.format_indian_number(round(self.current_cost)):<14}"
                            f"{self.format_indian_number(round(total_cost)):<14}" 
                            f"{(str(round(solve_time, 2)) + 's'):<12}"
                            f"{ampl.solve_result:<14}{'No':<10}"
                            f"{round(time.time() - self.start_time, 2)}s ")

            if improved:
                self.dia_red_iteration = self.dia_red_iteration + 1
                self.arc_dia_red()
                break



    def diameter_reduction2(self):
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
                ampl.read("reduced_wdnmodel.mod")
            ampl.read_data(self.data_file)

            # for (x, y, k), val in self.l.items():
            #     ampl.eval(f'let l[{x},{y},{k}] := {val};')
            # for (x, y), val in self.q.items():
            #     ampl.eval(f'let q[{x},{y}] := {val};')
            #     if self.data_number ==5:
            #         ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
            #         ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            # for x, val in self.h.items():
            #     ampl.eval(f'let h[{x}] := {val};') 

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
            
            # for k in self.pipes:
            #     if k>=arc_max_dia[i,j]:
            #         ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = 0;")
            
            # ampl.eval(f"s.t. flow_value{i}_{j}: q[{i}, {j}]>=abs({self.q[i,j]});")
            # if self.h[i] - self.h[j]>= 0:
            #     ampl.eval(f"s.t. headloss1{i}_{j}: h[{i}]- h[{j}]>={self.h[i] - self.h[j]} + 1e-0;")
            # else:
            #     ampl.eval(f"s.t. headloss2{i}_{j}: h[{i}]- h[{j}]<={self.h[i] - self.h[j]} - 1e-0;")
            # for key, val in self.l.items():
            #     u,v,k = key
            #     if (key[0], key[1]) != (i,j):
            #         ampl.eval(f"subject to con3__{u}_{v}_{k}: l[{u},{v},{k}] = {self.l[u,v,k]};")

            for k in self.pipes:
                if k==arc_max_dia[i,j]-1:
                    ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = {self.L[i,j]};")
                else:
                    ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = 0;")

            ampl.option['solver'] = "ipopt" 

            ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = no halt_on_ampl_error = yes"
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
            ampl.eval("display q;")
            ampl.eval("display h;")
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
            # break

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
        # remaining_arcs = [c for c, nc in sorted_by_abs_dual if nc not in self.visited_arc_reverse]
        # print("All cycles found:", len(sorted_by_abs_dual))
        # print("Visited cycles:", len(self.visited_arc_reverse))
        # print("Remaining new cycles:", len(remaining_arcs))


        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")
        print("----------------------------------------------------------------------------------------")

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

            # for (i,j) in self.arcs:
            #     if (i,j)!=(u,v):
            #         if self.q[i,j]>=0:
            #             ampl.eval(f"s.t. flow_dir{i}_{j}: q[{i},{j}]>=0;")
            #         else:
            #             ampl.eval(f"s.t. flow_dir{i}_{j}: q[{i},{j}]<=0;")

            #self.ampl.eval(f"set inarc := {{{inarc_set}}};")
            #self.ampl.eval(f"set indegree_node := {{{set(self.indegree_2_or_more)}}};")
            fix_arc_set = self.fix_leaf_arc_flow()
            # self.update_initial_points(self.l, self.q, self.h)
            # self.update_initial_points_with_perturbation(self.ampl, self.l, self.q, self.h) 
            if self.q[u,v] >= 0:
                ampl.eval(f"s.t. flow_direction1{u}_{v}: q[{u}, {v}]<=0;")
            else:
                ampl.eval(f"s.t. flow_direction1{u}_{v}: q[{u}, {v}]>=0;")

            for (i,j) in self.visited_arc_reverse:
                if self.q[i,j]>=0:
                    ampl.eval(f"s.t. flow_direction{i}_{j}: q[{i}, {j}]>=0;")
                else:
                    ampl.eval(f"s.t. flow_direction{i}_{j}: q[{i}, {j}]<=0;")

            # with self.suppress_output():
            ampl.option["solver"] = "ipopt"
            ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes max_iter = 200")   #max_iter = 1000
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
                    # print(self.q)
                    # ampl.eval("display q;")
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
                    # ampl.eval("display q;")
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
            # headloss = 0
            # for (i,j) in self.arcs:
            #     headloss = headloss + np.abs(h[i] - h[j])
            # print("headloss:", headloss)

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
        self.ampl.set_option("ipopt_options", f"outlev = 3 tol = 1e-9 bound_relax_factor=0  bound_push = {self.bound_push} bound_frac = {self.bound_frac} halt_on_ampl_error = yes warm_start_init_point = no expect_infeasible_problem = no")   #max_iter = 1000
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
            self.ampl.option["solver"] = "ipopt"
            #print("bound_push:", self.bound_push)
            #print("bound_frac:", self.bound_frac)
            #self.bound_push , self.bound_frac = (0.001, 0.000001)
            self.ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = no halt_on_ampl_error = yes")   #max_iter = 1000
            # self.ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  tol = 1e-9 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = yes halt_on_ampl_error = yes")   #max_iter = 1000
            #self.ampl.set_option("ipopt_options", f"outlev = 0 expect_infeasible_problem = no  bound_relax_factor=0 warm_start_init_point = yes halt_on_ampl_error = yes")   #max_iter = 1000
            # self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = yes halt_on_ampl_error = yes "
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

            self.ampl.solve()
            #self.ampl.eval("display q;")
            self.solve_result = self.ampl.solve_result
            self.total_cost = self.ampl.get_objective("total_cost").value()
        #self.ampl.eval("display q;")
        # print("Objective:", self.total_cost)
        # print("solve_result: ",self.solve_result)
        # print("eps:", epsilon)
        solve_time = self.ampl.get_value('_solve_elapsed_time')
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
        # self.ampl.eval("display l;")
        # self.ampl.eval("display q;")
        # self.ampl.eval("display h;")
        # self.ampl.eval("display {i in nodes} E[i] + P[i];")
        # self.ampl.eval("display {i in nodes} h[i] - E[i] - P[i];")
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
        self.export_solution()
        self.plot_interactive_wdn(f"../figure/json_file/solution_{self.data_number}.json",node_position[self.data_number])

        print("\n--------------------------Continuous Variable Branching Approach-------------------------------")
        self.cont_branching_iteration = 1
        self.visited_arc_reverse = []
        # self.iterate_branching_on_q()

        print("\n-----------------------------Flow change in cycle Approach---------------------------------")
        self.flow_change_in_cycle_iteration = self.cont_branching_iteration + 1
        # self.flow_change_in_cycle()

        print("---------------------------Reverse Arc Direction Approach------------------------------------")
        self.iteration = 1
        self.visited_arc_reverse = []
        # self.iterate_acyclic_flows() 

        # print("\n-----------------------------------Flow change in cycle Approach---------------------------")
        # self.flow_change_in_cycle_iteration = self.iteration + 1
        # self.flow_change_in_cycle()

        # print("\n----------------------------Diameter Reduction Approach--------------------------------------")
        # self.dia_red_iteration = 1
        # self.visited_arc = []
        # self.diameter_reduction()

        # print("\n--------------------------------Head increase Approach---------------------------------------")
        # self.dia_red_iteration = 1
        # self.visited_nodes = []
        # self.head_increase()

        print("\n-----------------------------Headloss increase Approach------------------------------------")
        self.headloss_increase_iteration = self.flow_change_in_cycle_iteration + 1
        self.visited_arc = []
        # self.headloss_increase()

        print("\n----------------------------Diameter Reduction Approach------------------------------------")
        self.dia_red_iteration = self.headloss_increase_iteration + 1
        self.visited_arc = []
        # self.diameter_reduction()
        # self.arc_dia_red()

        # print("\n--------------------------Continuous Variable Branching Approach---------------------------")
        # self.cont_branching_iteration = self.dia_red_iteration + 1
        # self.visited_arc_reverse = []
        # self.iterate_branching_on_q()

        # print("\n---------------------------Reverse 2 Arc Direction Approach------------------------------------")
        # self.visited_arc_reverse = []
        # self.iteration = 1
        # self.iterate_acyclic_flows() 
        # self.iterate_two_arc_reversals()

        # print("\n-----------------------------------Max Flow Approach-----------------------------------------")
        self.max_flow_iteration = 1
        self.visited_arc = []
        # self.max_flow() 
        # print("\n-----------------------------------Head increase Approach-----------------------------------------")
        self.head_increase_iter = 1
        self.visited_nodes = []
        # self.head_increase()
        # print("\n-----------------------------------Flow change in cycle Approach-----------------------------------------")
        # self.flow_change_in_cycle_iteration = 1
        # self.flow_change_in_cycle()
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
