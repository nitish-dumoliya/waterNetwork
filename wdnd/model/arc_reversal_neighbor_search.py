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
import dash
from dash import dcc, html, Input, Output
import plotly.io as pio

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
        self.max_iter = 2000
        self.bound_push , self.bound_frac = (0.1, 0.1)
        self.mu_init = 1e-2        

        # Trust-region parameters
        self.eta_min = 0.01
        self.eta = self.eta_min 
        self.alpha_min = 0.12
        self.alpha = self.alpha_min
        self.Delta_min = 0.0001
        self.Delta = self.Delta_min 
        self.Delta_h = 0.001
        self.eta_max = 1.00
        self.Delta_max = 10.00   # adjust based on flow scale
        self.alpha_shrink = 0.1
        self.alpha_expand = 1.1 
        self.max_tr_failures = 50
        self.tr_failure_count = 0
        self.Terminate = False
        self.success_streak = 0
        self.fail_streak = 0


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
        self.Q_max = self.ampl.getParameter('Q_max').to_list()[0]
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
    
    def export_solution(self, iteration):
        solution = {
            "iteration": iteration,
            "objective": self.current_cost,
            "solve_time": self.elapsed_time, 
            "q": {f"{i},{j}": float(v) for (i,j), v in self.q.items()},
            # "q": self.q,
            "h": self.h,
            "l": {f"{i},{j},{k}": self.l[i,j,k] 
                  for (i,j) in self.arcs for k in self.pipes if self.l[i,j,k] > 1e-3},
            "D": self.D,
            "Diameters": {f"{i}": self.d[i] for i in self.pipes},
            "L": {f"{i},{j}": v for (i,j), v in self.L.items()},
            "C": self.C,
            "hmin": {f"{i}": self.E[i] if i == list(self.source)[0] else self.E[i] + self.P[i] for i in self.nodes},
            # "pipe_diameters": self.l,
            "nodes": list(self.nodes),
            "arcs": [(i,j) for (i,j) in self.arcs],
            "source": list(self.source)
        }
    
        with open(f"../figure/json_file/d{self.data_number + 1}/solution_{iteration}.json", "w") as f:
            json.dump(solution, f, indent=2)

    def export_solution_iteration(self, iteration):
        solution = {
            "iteration": iteration,
            "objective": self.current_cost,
            "solve_time": time.time() - self.start_time, 
            "q": {f"{i},{j}": float(v) for (i,j), v in self.q.items()},
            # "q": self.q,
            "h": self.h,
            "l": {f"{i},{j},{k}": self.l[i,j,k] 
                  for (i,j) in self.arcs for k in self.pipes if self.l[i,j,k] > 1e-3},
            "D": self.D,
            "Diameters": {f"{i}": self.d[i] for i in self.pipes},
            "L": {f"{i},{j}": v for (i,j), v in self.L.items()},
            "C": self.C,
            "hmin": {f"{i}": self.E[i] if i == list(self.source)[0] else self.E[i] + self.P[i] for i in self.nodes},
            # "pipe_diameters": self.l,
            "nodes": list(self.nodes),
            "arcs": [(i,j) for (i,j) in self.arcs],
            "source": list(self.source)
        }
    
        with open(f"../figure/json_file/d{self.data_number+1}/solution_{self.iteration}.json", "w") as f:
            json.dump(solution, f, indent=2)

    def build_plot(self,iteration, solution_file, node_pos, data_number, heuristic_approach, node_head_diff, arc_diff, edge):
        # --------------------------------------------------
        # Load solution
        # --------------------------------------------------
        with open(solution_file) as f:
            sol = json.load(f)
        def parse_key(k):
            return tuple(map(int, k.replace("(", "").replace(")", "").split(",")))
        q = {parse_key(k): v for k, v in sol["q"].items()}
        L = {parse_key(k): v for k, v in sol.get("L", {}).items()}
        # Pipe info l[i,j,d]
        pipe_info = {}
        for k_str, val in sol.get("l", {}).items():
            if val < 1e-6:
                continue
            i, j, d = map(int, k_str.replace("(", "").replace(")", "").split(","))
            pipe_info.setdefault((i, j), []).append((d, val))
        h    = {int(k): v for k, v in sol["h"].items()}
        hmin = {int(k): v for k, v in sol["hmin"].items()}
        D    = {int(k): v for k, v in sol["D"].items()}
        # source = sol["source"]
        source = sol.get("source", [])
        diameters = {int(k): v for k, v in sol["Diameters"].items()}
        # --------------------------------------------------
        # Figure scaling + node radius
        # --------------------------------------------------
        node_marker_size = 22
        BASE = 1800
        xs = [p[0] for p in node_pos.values()]
        ys = [p[1] for p in node_pos.values()]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        data_w = xmax - xmin
        data_h = ymax - ymin
        if data_w >= data_h:
            fig_w = BASE
            fig_h = int(BASE * data_h / data_w)
        else:
            fig_h = BASE
            fig_w = int(BASE * data_w / data_h)
        node_radius = (node_marker_size / 2) * (data_w / fig_w)
        # --------------------------------------------------
        # FIXED DIAMETER ORDER FROM DICTIONARY (DETERMINISTIC)
        # --------------------------------------------------
        # Use dictionary keys as the dataset order
        DATASET_DIAMETERS = list(diameters.keys())
        # --------------------------------------------------
        # VISUAL PALETTES (must cover all diameters)
        # --------------------------------------------------
        thickness_levels = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        color_palette = [
             "#1f77b4", "#ff7f0e", "#2ca02c", "#008080",   # ← replaced red
             "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
             "#bcbd22", "#17becf", "#393b79", "#637939",
             "#8c6d31", "#843c39", "#7b4173", "#3182bd",
             "#31a354", "#756bb1", "#636363", "#e6550d",
             "#969696", "#6baed6"
        ]

        assert len(DATASET_DIAMETERS) <= len(color_palette)
        assert len(DATASET_DIAMETERS) <= len(thickness_levels)
        # --------------------------------------------------
        # DIAMETER → INDEX (STABLE)
        # --------------------------------------------------
        DIAM_INDEX = {d: i for i, d in enumerate(DATASET_DIAMETERS)}
        # --------------------------------------------------
        # BUILD STYLE MAP (INDEPENDENT OF SOLUTION)
        # --------------------------------------------------
        all_diams = {d for pipes in pipe_info.values() for d, _ in pipes}
        diam_to_style = {}
        for d in all_diams:
            idx = DIAM_INDEX[d]
            diam_to_style[d] = {
                "width": thickness_levels[idx],
                "color": color_palette[idx]
            }
        # --------------------------------------------------
        # Containers
        # --------------------------------------------------
        edge_groups = {}  # d → lines
        arrows = []
        click_x, click_y, click_text = [], [], []
        flow_tx, flow_ty, flow_txt = [], [], []
        pipe_tx, pipe_ty, pipe_txt = [], [], []
        shapes = []

        traces = []
        # --------------------------------------------------
        # EDGES (arc hover + multi-diameter coloring)
        # --------------------------------------------------
        for (i, j), flow in q.items(): 
            if i not in node_pos or j not in node_pos:
                continue
            x0, y0 = node_pos[i]
            x1, y1 = node_pos[j]
            xs, ys, xe, ye = (x0, y0, x1, y1) if flow >= 0 else (x1, y1, x0, y0)
            dx, dy = xe - xs, ye - ys
            dist = math.hypot(dx, dy)
            if dist < 1e-9:
                continue
            sx = xs + node_radius * dx / dist
            sy = ys + node_radius * dy / dist
            ex = xe - node_radius * dx / dist
            ey = ye - node_radius * dy / dist
            mx, my = (sx + ex) / 2, (sy + ey) / 2
            # --------------------------------------------------
            # Hover info (UNCHANGED – arc level)
            # --------------------------------------------------
            pipes = pipe_info.get((i, j), [])
            pipe_label = ", ".join([f"D{d}:{l:.1f}" for d, l in pipes])
            hover = (
                f"<b>Arc {i} → {j}</b><br>"
                f"Flow: {flow:.6f}<br>"
                f"Length: {L.get((i,j),0):.2f}<br>"
                f"<b>Pipes:</b><br>{pipe_label}"
            )
            click_x.append(mx)
            click_y.append(my)
            click_text.append(hover)
            flow_tx.append(mx)
            flow_ty.append(my)
            flow_txt.append(f"{flow:.5f}")
            pipe_tx.append(mx)
            pipe_ty.append(my - 0.02 * data_h)
            pipe_txt.append(pipe_label)
            # --------------------------------------------------
            # Multi-diameter segmentation
            # --------------------------------------------------
            total_len = sum(l for _, l in pipes)
            if total_len <= 0:
                continue
            cur_len = 0.0
            # for d, l in pipes:
            for d, l in sorted(pipes, key=lambda x: x[0], reverse=True):
                frac_s = cur_len / total_len
                frac_e = (cur_len + l) / total_len
                seg_sx = sx + frac_s * dx
                seg_sy = sy + frac_s * dy
                seg_ex = sx + frac_e * dx
                seg_ey = sy + frac_e * dy
                style = diam_to_style[d]
                edge_groups.setdefault(d, {"x": [], "y": []})
                edge_groups[d]["x"] += [seg_sx, seg_ex, None]
                edge_groups[d]["y"] += [seg_sy, seg_ey, None]
                # traces.append(go.Scatter(
                # x=[seg_sx, seg_ex],
                # y=[seg_sy, seg_ey],
                # mode="lines",
                # line=dict(
                #     color=style["color"],
                #     width=style["width"],
                #     # cap="round"
                # ),
                # hoverinfo="skip",
                # showlegend=False   # legend added separately
                # ))
                cur_len += l
            # --------------------------------------------------
            # Arrow (single arrow for whole arc)
            # --------------------------------------------------
            arrow_color = "red" if (i, j) == edge else "black"
            d_arrow = max(d for d, _ in pipes)
            style = diam_to_style[d_arrow]
            arrows.append(dict(
                ax=sx, ay=sy, x=ex, y=ey,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=3,
                arrowsize=2,
                arrowwidth=1,
                arrowcolor= arrow_color
                # arrowcolor= style["color"]
            ))
        # --------------------------------------------------
        # CREATE TRACES (FIXED LEGEND ORDER)
        # --------------------------------------------------
        # traces = []
        for d in DATASET_DIAMETERS:
            if d not in edge_groups:
                continue
            data = edge_groups[d]
            style = diam_to_style[d]
            traces.append(go.Scatter(
                x=data["x"],
                y=data["y"],
                mode="lines",
                line=dict(
                    width=style["width"],
                    color=style["color"]
                ),
                name=f"Diameter D{d}: {diameters[d]} m",
                hoverinfo="skip"
            ))

        traces.append(go.Scatter(
            x=click_x, y=click_y,
            mode="markers",
            marker=dict(size=0, opacity=0),
            hoverinfo="text",
            hovertext=click_text,
            showlegend=False
        ))

        # color the revered arc or diameter reduction arc
        i,j = edge[0], edge[1]
        if i in node_pos and j in node_pos:
            x0, y0 = node_pos[i]
            x1, y1 = node_pos[j]

            traces.append(go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(
                    color="red",
                    width=5
                ),
                hoverinfo="skip",
                showlegend=False
            ))

        # --------------------------------------------------
        # Nodes
        # --------------------------------------------------
        nx_src, ny_src, label_src, text_src = [], [], [], []
        nx_dem, ny_dem, label_dem, text_dem = [], [], [], []
        hx_src, hy_src, htxt_src = [], [], []
        hx_dem, hy_dem, htxt_dem = [], [], []

        hx_min, hy_min, hmin_txt = [], [], []
        hx_demval, hy_demval, dem_txt = [], [], []

        # -------------------------------
        # Populate
        # -------------------------------
        for n, (x, y) in node_pos.items():
            hover_text = (
                f"<b>Node {n}</b><br>"
                f"Head: {h.get(n,0):.2f}<br>"
                f"Min Head: {hmin.get(n,0):.2f}<br>"
                f"Demand: {D.get(n,0):.5f}"
            )
            if n in source:
                nx_src.append(x); ny_src.append(y)
                label_src.append(str(n))
                text_src.append(hover_text)
                hx_src.append(x)
                hy_src.append(y + 0.03 * data_h)
                htxt_src.append(f"{h.get(n,0):.2f}")
            else:
                nx_dem.append(x); ny_dem.append(y)
                label_dem.append(str(n))
                text_dem.append(hover_text)
                hx_dem.append(x)
                hy_dem.append(y + 0.03 * data_h)
                htxt_dem.append(f"{h.get(n,0):.2f}")
            # ---------------------------------
            # Min head (below node)
            # ---------------------------------
            hx_min.append(x - 0.03 * data_w)
            hy_min.append(y - 0.02 * data_h)
            hmin_txt.append(f"{hmin.get(n,0):.2f}")

            # ---------------------------------
            # Demand (below-left)
            # ---------------------------------
            hx_demval.append(x)
            hy_demval.append(y - 0.035 * data_h)
            dem_txt.append(f"{D.get(n,0):.4f}")
        # -------------------------------
        # Source nodes trace
        # -------------------------------
        traces.append(go.Scatter(
            x=nx_src, y=ny_src,
            mode="markers+text",
            text=label_src,
            textposition="middle center",
            hoverinfo="text",
            hovertext=text_src,
            marker=dict(
                size=node_marker_size + 4,
                color="royalblue",
                symbol="circle",
                line=dict(width=2, color="black")
            ),
            name="Source Nodes"
        ))
        # -------------------------------
        # Demand nodes trace
        # -------------------------------
        traces.append(go.Scatter(
            x=nx_dem, y=ny_dem,
            mode="markers+text",
            text=label_dem,
            textposition="middle center",
            hoverinfo="text",
            hovertext=text_dem,
            marker=dict(
                size=node_marker_size,
                color="skyblue",
                line=dict(width=2, color="black")
            ),
            name="Demand Nodes"
        ))
        # --------------------------------------------------
        # Min head requirement (below node)
        # --------------------------------------------------
        traces.append(go.Scatter(
            x=hx_min, y=hy_min,
            mode="text",
            text=hmin_txt,
            textfont=dict(size=12, color="darkred"),
            name="Min Head Requirement",
            hoverinfo="skip",
            visible="legendonly"
        ))
        # --------------------------------------------------
        # Node demand (below-left)
        # --------------------------------------------------
        traces.append(go.Scatter(
            x=hx_demval, y=hy_demval,
            mode="text",
            text=dem_txt,
            textfont=dict(size=12, color="darkgreen"),
            name="Node demand",
            hoverinfo="skip",
            visible="legendonly"
        ))
        # --------------------------------------------------
        # Toggleable text layers
        # --------------------------------------------------
        traces += [
            go.Scatter(x=flow_tx, y=flow_ty, mode="text",
                       text=flow_txt, name="Flow Values",
                       # visible="legendonly"
                       ),
            go.Scatter(x=pipe_tx, y=pipe_ty, mode="text",
                       text=pipe_txt, name="Pipe Info",
                       # visible="legendonly"
                       ),
            go.Scatter(x=hx_src, y=hy_src, mode="text",
                       text=htxt_src, name="Source Node Head",
                       # visible="legendonly"
                       ),
            go.Scatter(x=hx_dem, y=hy_dem, mode="text",
                       text=htxt_dem, name="Demand Node Head",
                       # visible="legendonly"
                       )
        ] 

        # -------------------------------
        # Difference arcs
        # -------------------------------
        mx, my, htxt_arc = [], [], [] 
        for (i, j), info in arc_diff.items():
            if i not in node_pos or j not in node_pos:
                continue
            x0, y0 = node_pos[i]
            x1, y1 = node_pos[j] 
            # Midpoint
            xm = 0.5 * (x0 + x1)
            ym = 0.5 * (y0 + y1)
            txt = f"<b>Arc ({i} → {j})</b><br>"
            if "flow" in info:
                f = info["flow"]
                txt += (
                    f"Flow change<br>"
                    f"q₁ = {f['q_sol1']}<br>"
                    f"q₂ = {f['q_sol2']}<br>"
                    f"Δq = {f['delta']}<br>"
                )
            if "pipes" in info:
                txt += (
                    f"Pipes changed<br>"
                    f"sol1: {info['pipes']['sol1']}<br>"
                    f"sol2: {info['pipes']['sol2']}"
                )
            mx.append(xm)
            my.append(ym)
            htxt_arc.append(txt)
        traces.append(
            go.Scatter(
                x=mx,
                y=my,
                mode="markers",
                marker=dict(
                    size=12,
                    color="rgba(0,0,0,0)"   # 👈 fully invisible
                ),
                hoverinfo="text",
                hovertext=htxt_arc,
                showlegend=False,
                legendgroup="diff",
                visible="legendonly"
            )
        )

        arc_set = set(arc_diff.keys())
        ax_diff, ay_diff = [], []
        for (i, j) in arc_set:
            if i not in node_pos or j not in node_pos:
                continue 
            x0, y0 = node_pos[i]
            x1, y1 = node_pos[j] 
            ax_diff += [x0, x1, None]
            ay_diff += [y0, y1, None] 
        traces.append(go.Scatter(
            x=ax_diff, y=ay_diff,
            mode="lines",
            line=dict(
                width=3,
                color="red",
                dash="dash"
            ),
            hoverinfo="text",
            hovertext=["Arc changed"] * len(ax_diff),
            name="Affected arcs",
            legendgroup="diff",
            visible="legendonly"   # 👈 OFF by default
        ))

        # -------------------------------
        # Difference nodes
        # -------------------------------
        node_set = set(node_head_diff.keys())
        dx, dy, dlabel, dtext = [], [], [], []
        nocdx, nocdy, nocdlabel, nocdtext = [], [], [], []
        nocsx, nocsy, nocslabel, nocstext  = [], [], [], []
        # for n, info in node_head_diff.items():
        for n, (x, y) in node_pos.items():
            if n in node_set:
                # x, y = node_pos[n]
                dx.append(x)
                dy.append(y)
                dlabel.append(str(n))
                # dtext.append(f"Node {n}<br>Head changed")
                info = node_head_diff[n]
                dtext.append(f"<b>Node {n}</b><br>"
                             f"Head changed<br>"
                             f"h₁ = {info['h_sol1']}<br>"
                             f"h₂ = {info['h_sol2']}<br>"
                             f"Δh = {info['delta']}")
            else:
                if n in source:
                    nocsx.append(x)
                    nocsy.append(y)
                    nocslabel.append(str(n)) 
                    nocstext.append(f"<b>Source Node {n}</b><br>"
                                    f"Head fixed<br>"
                                    f"h = {h.get(n,0):.2f}<br>"
                                    )

                else:
                    nocdx.append(x)
                    nocdy.append(y)
                    nocdlabel.append(str(n))
                    nocdtext.append(f"<b>Node {n}</b><br>"
                             f"Head Unchanged<br>"
                             # f"h₁ = {info['h_sol1']}<br>"
                             f"h₂ = {h.get(n,0):.2f}<br>"
                             # f"Δh = {info['delta']}"
                                    )
        traces.append(go.Scatter(
            x=dx, y=dy,
            mode="markers + text",
            text=dlabel,
            textposition='middle center',
            marker=dict(
                size=node_marker_size,
                color="skyblue",
                # color="rgba(135, 206, 235, 0.0)",   # transparent red
                symbol="circle",
                line=dict(width=2, color="red")
            ),
            hoverinfo="text",
            hovertext=dtext,
            name="Affected nodes",
            legendgroup="diff",
            legendgrouptitle=dict(text="Local solution difference"),
            visible="legendonly"   # 👈 OFF by default
        )) 
        traces.append(go.Scatter(
            x=nocdx, y=nocdy,
            mode="markers + text",
            text=nocdlabel,
            textposition='middle center',
            marker=dict(
                size=node_marker_size,
                color="skyblue",
                # color="rgba(135, 206, 235, 0.0)",   # transparent red
                symbol="circle",
                line=dict(width=2, color="black")
            ),
            hoverinfo="text",
            hovertext=nocdtext,
            name="Unaffected nodes",
            legendgroup="diff",
            # legendgrouptitle=dict(text="Local solution difference"),
            visible="legendonly"   # 👈 OFF by default
        ))
        traces.append(go.Scatter(
            x=nocsx, y=nocsy,
            mode="markers + text",
            text=nocslabel,
            textposition='middle center',
            marker=dict(
                size=node_marker_size + 4,
                color="royalblue",
                # color="rgba(135, 206, 235, 0.0)",   # transparent red
                symbol="circle",
                line=dict(width=2, color="black")
            ),
            hoverinfo="text",
            hovertext=nocstext,
            # name="Unaffected nodes",
            legendgroup="diff",
            # legendgrouptitle=dict(text="Local solution difference"),
            # visible="legendonly"   # 👈 OFF by default
            showlegend=False
        ))

        # --------------------------------------------------
        # Figure
        # --------------------------------------------------
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
        fig = go.Figure(traces)
        # fig.update_traces(line=dict(simplify=False))
        fig.update_layout(
            # title=f"WDN = {data_list[data_number]}, Total Cost = {sol['objective']:.2f}",
            title=dict(
                text=f"Network: {data_list[data_number]} | Heuristic Approach: {heuristic_approach} | Iteration: {iteration} | Arc: {edge} | Objective: {sol['objective']:.2f} | Time: {sol['solve_time']:.2f} s",
                font=dict(
                    size=24,        # Increase font size here
                    color='black',  # Optional: title color
                    family="Arial"  # Optional: font family
                ),
                # x=0.5,             # Optional: center the title
                # xanchor='center'
            ),
            # title = f"Network: {data_list[data_number]} | Heuristic Approach: {heuristic_approach} | Iteration: {self.iteration} | Arc: {edge} | Objective: {sol['objective']:.2f} | Time: {sol['solve_time']:.2f} s",
            annotations=arrows,
            dragmode="pan",
            hovermode="closest",
            # hovermode="x unified",
            autosize=False,
            width=1900,
            height=1100,
            xaxis=dict(visible=True,showgrid=True,zeroline=True, fixedrange=False, showline=True,linecolor='black',mirror=True,ticks='outside',),
            yaxis=dict(visible=True,showgrid=True,zeroline=True, scaleanchor="x", fixedrange=False, showline=True,linecolor='black',mirror=True,ticks='outside',),
            # xaxis=dict(
            #     scaleanchor="y",
            #     scaleratio=1,
            #     constrain="domain"
            # ),
            # yaxis=dict(
            #     scaleanchor="x",
            #     scaleratio=1,
            #     constrain="domain"
            # ),
            plot_bgcolor="white",
            paper_bgcolor="white",
            legend=dict(
                title=dict(
                    text="Network Elements",
                    font=dict(
                        size=20,          # Bigger title
                        family="Arial",
                        color="black"
                    )
                ),
                font=dict(
                    size=16,              # Legend item font size
                    family="Arial",
                    color="black"
                ),
                orientation="v",
                x=1.02,
                y=1,
                xanchor="left",
                yanchor="top",
                bgcolor="rgba(255,255,255,0.85)",   # Light background
                bordercolor="rgba(0,0,0,0.3)",
                borderwidth=2,
                itemsizing="constant",
                itemwidth=40
            )

            # legend=dict(
            #     title=dict(text="Network Elements"),
            #     orientation="v",
            #     x=1.02,
            #     y=1,
            #     bordercolor="rgba(0,0,0,0.15)",
            #     borderwidth=1
            # ),
            # shapes=shapes
        )
        fig.write_html(
            f"../figure/json_file/d{data_number+1}/wdn_interactive_solution{iteration}.html",
            auto_open=False,
            config={"scrollZoom": True, 
                    }
        )
        pio.write_image(fig, f"../figure/json_file/d{data_number+1}/wdn_interactive_solution{self.iteration}.pdf", format="pdf")
        # network_info = f"Network: {data_list[data_number]} | Objective: {sol['objective']:.2f} | Time: {sol['solve_time']:.2f} s"
        # print("✔ Diameter-colored & thickness-separated WDN plot saved")
        return fig

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

    def fix_leaf_arc_flow1(self):
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

    def compare_two_local_solutions(self,
        l1, q1, h1,
        l2, q2, h2,
        tol=1e-5
    ):
        """
        Compare two locally optimal solutions (dictionary-based)
        Returns:
            node_head_diff : dict
            arc_diff       : dict
        """
        # --------------------------------------------------
        # 1. Node head differences
        # --------------------------------------------------
        node_head_diff = {}
    
        for i in h1:
            if abs(h1.get(i, 0.0) - h2.get(i, 0.0)) > tol:
                node_head_diff[i] = {
                    "h_sol1": round(h1[i], 4),
                    "h_sol2": round(h2[i], 4),
                    "delta": round(h2[i] - h1[i], 4)
                }
        # --------------------------------------------------
        # 2. Arc differences
        # --------------------------------------------------
        arc_diff = {}
        arcs = set(q1.keys()) | set(q2.keys())
        for (i, j) in arcs:
            arc_changed = False
            info = {}
            # -------------------------
            # Flow change
            # -------------------------
            q1_ij = q1.get((i, j), 0.0)
            q2_ij = q2.get((i, j), 0.0)
            if abs(q1_ij - q2_ij) > tol:
                arc_changed = True
                info["flow"] = {
                    "q_sol1": round(q1_ij, 4),
                    "q_sol2": round(q2_ij, 4),
                    "delta": round(q2_ij - q1_ij, 4)
                }
            # -------------------------
            # Pipe length / diameter structure change
            # -------------------------
            pipes1 = {
                k: round(l1[(i, j, k)], 4)
                for (ii, jj, k) in l1
                if ii == i and jj == j and l1[(ii, jj, k)] > tol
            }
            pipes2 = {
                k: round(l2[(i, j, k)], 4)
                for (ii, jj, k) in l2
                if ii == i and jj == j and l2[(ii, jj, k)] > tol
            }
            if pipes1 != pipes2:
                arc_changed = True
                info["pipes"] = {
                    "sol1": pipes1,
                    "sol2": pipes2
                }
            if arc_changed:
                arc_diff[(i, j)] = info
        return node_head_diff, arc_diff

    def reduced_nlp_model(self, iteration, l_cvx, q_cvx, h_cvx, l, q, h, z_star, sorted_arcs, l_points, q_points, dual_dict):
        arc_max_dia = {}
        if self.data_number == 6:
            self.fixarcs = self.ampl.getSet('fixarcs')
            #print("fixarcs:",self.fixarcs)
            for (i, j, d), val in l_cvx.items():
                if (i,j) not in self.fixarcs or (j,i) not in self.fixarcs:
                    if val > 1e-3:
                        if (i, j) not in arc_max_dia:
                            arc_max_dia[(i, j)] = d
                        else:
                            arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)
        else:
            for (i, j, d), val in l_cvx.items():
                if val > 1e-3:
                    if (i, j) not in arc_max_dia:
                        arc_max_dia[(i, j)] = d
                    else:
                        arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)

        ampl_cvxnlp = AMPL()
        ampl_cvxnlp.reset()
        ampl_cvxnlp.read("wdn.mod")
        if self.data_number == 5:
            ampl_cvxnlp.eval("""param exdiam{arcs};
                                param excost{arcs};
                                param d_max = max{i in pipes} d[i];
                                param d_min = min{i in pipes} d[i];

                                param R{arcs};		   # Roughness of each commercial pipe
                                param R_min = min{(i,j) in arcs} R[i,j];
                                param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
                                param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);
                                subject to con1{j in nodes diff Source}:
                                    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
                                ;
                                subject to con2{(i,j) in arcs}: 
                                    h[i] - h[j]  = (q1[i,j])^3 *((((q1[i,j])^2 + eps[i,j]^2)^0.426) /((q1[i,j])^2 + 0.426*eps[i,j]^2)) *omega * L[i,j] / ( (R[i,j]^1.852) * (exdiam[i,j])^4.87) ;
                                subject to con2_{(i,j) in arcs}: 
                                    h[i] - h[j]  = (q2[i,j])^3 *((((q2[i,j])^2 + eps[i,j]^2)^0.426) /((q2[i,j])^2 + 0.426*eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k]/(R[i,j]^1.852 * d[k]^4.87)) ;
                                subject to con9{(i,j) in arcs}: q[i,j] = q1[i,j] + q2[i,j];

                                subject to con3{(i,j) in arcs}: 
                                    sum{k in pipes} l[i,j,k] = L[i,j]
                                ;
                                subject to con4{(i,j) in arcs , k in pipes}: 
                                    l[i,j,k] <= L[i,j]
                                ;
                                """)
        elif self.data_number == 6:
            ampl_cvxnlp.eval(""" 
                                set fixarcs within {i in nodes, j in nodes: i != j};
                                param d_min; # minimum diameter (should match the minimum discrete diameter)
                                param d_max; # maximum diameter (should match the maximum discrete diameter)
                                param fix_r{fixarcs};
                                param fix_c{fixarcs};
                                param fixdiam{fixarcs};
                                param R{pipes};		   # Roughness of each commercial pipe
                                param R_min = min{k in pipes} R[k];
                                param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
                                param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);
                                subject to con2{(i,j) in arcs diff fixarcs}: 
                                   h[i] - h[j]  = (q[i,j])^3 *((((q[i,j])^2 + eps[i,j]^2)^0.426) /((q[i,j])^2 + 0.426*eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k]/(R[k]^1.852 * d[k]^4.87));
                                subject to con2_{(i,j) in fixarcs}:
                                   h[i] - h[j]  = (q[i,j])^3 *((((q[i,j])^2 + eps[i,j]^2)^0.426) /((q[i,j])^2 + 0.426*eps[i,j]^2)) * omega * L[i,j] / (fix_r[i,j]^1.852 * fixdiam[i,j]^4.87);
                                subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];
                                subject to con4{(i,j) in arcs diff fixarcs, k in pipes}: l[i,j,k] <= L[i,j];
                                """)
        else:
            ampl_cvxnlp.eval("""param R{pipes};		   # Roughness of each commercial pipe
                                param d_max = max{i in pipes} d[i];
                                param d_min = min{i in pipes} d[i];

                                param R_min = min{k in pipes} R[k];
                                param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
                                param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);
                                subject to con1{j in nodes diff Source}:
                                    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
                                ;
                                subject to con2{(i,j) in arcs}: 
                                    h[i] - h[j]  =  (q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426*eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
                                subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];
                                subject to con4{(i,j) in arcs , k in pipes}: 
                                    l[i,j,k] <= L[i,j]
                                ;
                                """)  

        ampl_cvxnlp.read_data(self.data_file)

        for (x, y, k), val in l.items():
            ampl_cvxnlp.eval(f'let l[{x},{y},{k}] := {val};')
        for (x, y), val in q.items():
            ampl_cvxnlp.eval(f'let q[{x},{y}] := {val};')
            if self.data_number ==5:
                ampl_cvxnlp.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                ampl_cvxnlp.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
        for x, val in h.items():
            ampl_cvxnlp.eval(f'let h[{x}] := {val};') 
        
        # for (u), val in dual_dict.items():
        #     ampl_cvxnlp.param["lambda"][u] = val + 1e+1

        # current_duals = {}
        # for con_name, val in ampl_cvxnlp.get_constraints():
        #     dual_values = val.get_values()
        #     current_duals[con_name] = dual_values
        #
        # for con_name, dual_values in self.all_duals.items():
        #     if con_name in current_duals:
        #         ampl_cvxnlp.get_constraint(con_name).set_values(dual_values) 

        ampl_cvxnlp.eval("""param l_ref {arcs, pipes};""")
        ampl_cvxnlp.eval("""param q_ref {arcs};""")
        ampl_cvxnlp.eval("""param h_ref {nodes};""")


        for (u,v,k), val in l_cvx.items():
            ampl_cvxnlp.param["l_ref"][u, v, k] = val

        for (u,v), val in q_cvx.items():
            ampl_cvxnlp.param["q_ref"][u, v] = val
        for (u), val in h_cvx.items():
            ampl_cvxnlp.param["h_ref"][u] = val

        for (u,v) in self.arcs:
        # for (u,v) in sorted_arcs[:len(sorted_arcs)-10]:
            if self.data_number == 5:
                if self.q[u,v]>=0:
                    ampl_cvxnlp.eval(f"s.t. flow_neigh_{u}_{v}: {q[u,v] - self.Delta} <= q[{u}, {v}] <= {q[u,v] + self.Delta};")
                    ampl_cvxnlp.eval(f"s.t. flow_neigh_q1{u}_{v}: {self.q1[u,v] - self.Delta} <= q1[{u}, {v}] <= {self.q1[u,v] + self.Delta};")
                    ampl_cvxnlp.eval(f"s.t. flow_neigh_q2{u}_{v}: {self.q2[u,v] - self.Delta} <= q2[{u}, {v}] <= {self.q2[u,v] + self.Delta};")
                else:
                    ampl_cvxnlp.eval(f"s.t. flow_neigh_{u}_{v}:{q[u,v] - self.Delta} <= q[{u}, {v}] <= {q[u,v] + self.Delta};")
                    ampl_cvxnlp.eval(f"s.t. flow_neigh_q1{u}_{v}:{self.q1[u,v] - self.Delta} <= q1[{u}, {v}] <= {self.q1[u,v] + self.Delta};")
                    ampl_cvxnlp.eval(f"s.t. flow_neigh_q2{u}_{v}:{self.q2[u,v] - self.Delta} <= q2[{u}, {v}] <= {self.q2[u,v] + self.Delta};")
            else:
                if self.q[u,v]>=0:
                    ampl_cvxnlp.eval(f"s.t. flow_neigh_{u}_{v}: {q[u,v] - self.Delta} <= q[{u}, {v}] <= {q[u,v] + self.Delta};")
                else:
                    ampl_cvxnlp.eval(f"s.t. flow_neigh_{u}_{v}:{q[u,v] - self.Delta} <= q[{u}, {v}] <= {q[u,v] + self.Delta};")

        # for (u,v) in self.arcs:
        #     for k in self.pipes:
        #         ampl_cvxnlp.eval(f"""subject to TR_l_lower_{u}_{v}_{k}: l[{u},{v},{k}] >= {l[u,v,k] - self.eta};""")
        #         ampl_cvxnlp.eval(f"""subject to TR_l_upper_{u}_{v}_{k}: l[{u},{v},{k}] <= {l[u,v,k] + self.eta};""")

        ampl_cvxnlp.option['solver'] = "ipopt" 

        ampl_cvxnlp.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"
        ampl_cvxnlp.option["presolve_eps"]= "1.86e-7"
        return ampl_cvxnlp

    def solve_lp(self, q):
        ampl_lp = AMPL()
        ampl_lp.reset()

        ampl_lp.read("lp_model2.mod")
        # ampl.read(self.model_file)
        ampl_lp.read_data(self.data_file)
        # for (x, y, k), val in self.l.items():
        #     ampl.eval(f'let l[{x},{y},{k}] := {val};')
        # for (x, y), val in self.q.items():
        #     ampl.eval(f'let q[{x},{y}] := {val};')
                # if self.data_number ==5:
                #     ampl_cvxnlp.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                #     ampl_cvxnlp.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
        # for idx, val in self.h.items():
        #     ampl.eval(f'let h[{idx}] := {val};')

        for (u,v), val in q.items():
            ampl_lp.param["q"][u, v] = val

        ampl_lp.option['solver'] = "cplexamp"
        ampl_lp.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = yes bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = no halt_on_ampl_error = yes"
        ampl_lp.option["presolve_eps"]= "7.19e-7"

        with self.suppress_output():
            ampl_lp.solve()
        return ampl_lp

    def local_solution_improvement_heuristic(self):
        # print("Iteration :",self.iteration)
        abs_flows = sorted(
            abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-3
        )
        m = len(abs_flows)
        # median_flow = abs_flows[m // 2]   # floor(m/2)
        median_flow = abs_flows[m//2] 
        print(f"\nIteration {self.iteration} | Alpha={self.alpha}| median_flow= {median_flow} | Delta={self.Delta:.6f}" )
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

        self.all_duals = {}
        for con_name, con in self.ampl.get_constraints():
            self.all_duals[con_name] = con.getValues()
        
        # Remove fixed and already-visited arcs
        sorted_all_arcs = [
            arc for arc in self.arcs
            if arc not in self.fix_arc_set
        ]
        
        dual_dict = {}
        for con_name, dual_values in self.all_duals.items():
            if con_name == "con7":
                tmp = dual_values.to_dict()
                dual_dict = {
                    node: val for node, val in tmp.items()
                }
                break   # only con2 is needed

        # if self.iteration == 1:
      
        cvx_nlp = self.reduced_nlp_model(self.iteration, self.l_star, self.q_star, self.h_star, self.l, self.q, self.h, self.z_star, sorted_all_arcs, self.l_points, self.q_points, dual_dict)

        with self.suppress_output():
            cvx_nlp.solve()
        all_duals_reduced_model = {}
        for con_name, con in cvx_nlp.get_constraints():
            all_duals_reduced_model[con_name] = con.getValues()

        if cvx_nlp.solve_result == "solved":
            l_star = cvx_nlp.getVariable("l").getValues().to_dict()
            # self.L_ref.append((l_star))
            self.z_star = sum(self.C[k] * l_star[u, v, k] for (u, v) in self.arcs for k in self.pipes)
            # z = cvx_nlp.getObjective("total_energy_cost").value()
            solve_time = cvx_nlp.get_value('_solve_elapsed_time')
            self.l_star = cvx_nlp.getVariable('l').getValues().to_dict()
            self.q_star = cvx_nlp.getVariable('q').getValues().to_dict()
            self.h_star = cvx_nlp.getVariable('h').getValues().to_dict()
            if self.data_number==5:
                self.q1_star = cvx_nlp.getVariable('q1').getValues().to_dict()
                self.q2_star = cvx_nlp.getVariable('q2').getValues().to_dict()

            # print("-------------------NCVX NLP l Var--------------------:")
            # self.ampl.eval("display {(i,j) in arcs, k in pipes: l[i,j,k]>=1e-3}: l[i,j,k];")
            # self.ampl.eval("display {(i,j) in arcs: abs(q[i,j])>=0}: q[i,j];")
            # self.ampl.eval("display h;")
            # print("----------------Reduced NLP l Var------------------------:")
            # cvx_nlp.eval("display {(i,j) in arcs, k in pipes: l[i,j,k]>=1e-3}: l[i,j,k];")
            # cvx_nlp.eval("display {(i,j) in arcs: abs(q[i,j])>=0}: q[i,j];")
            # cvx_nlp.eval("display h;")
            # cvx_nlp.eval("display {i in nodes diff Source} max(0, E[i]+P[i] - h[i]);")

            # cvx_nlp.eval("display l;")
            self.l_points.append(self.l_star)
            self.q_points.append(self.q_star)

            # ampl_lp = self.solve_lp(q)
            # l = ampl_lp.getVariable('l').getValues().to_dict() 
            # q = ampl_lp.getVariable('q').getValues().to_dict() 
            # h = ampl_lp.getVariable('h').getValues().to_dict() 
            
            self.all_duals = {}
            for con_name, con in cvx_nlp.get_constraints():
                self.all_duals[con_name] = con.getValues()

            dual_dict = {}
            for con_name, dual_values in self.all_duals.items():
                if con_name == "flow_balance":
                    tmp = dual_values.to_dict()
                    dual_dict = {
                        arc: val for arc, val in tmp.items()
                    }
                    break   # only con2 is needed
            # print("dual_values of flow_balance:", dual_dict)

            print("z_star:", self.z_star, "solve_time:", solve_time)
            
            # arc_min_dia = {}
            # if self.data_number == 6:
            #     self.fixarcs = set(self.ampl.getSet('fixarcs'))
            #     for (i, j, d), val in self.l_star.items():
            #         # skip fixed arcs (both directions fixed)
            #         if (i, j) in self.fixarcs or (j, i) in self.fixarcs:
            #             continue
            #         if val > 1e-3:
            #             if (i, j) not in arc_min_dia:
            #                 arc_min_dia[(i, j)] = d
            #             else:
            #                 arc_min_dia[(i, j)] = min(arc_min_dia[(i, j)], d)
            # else:
            #     for (i, j, d), val in self.l_star.items():
            #         if val > 1e-3:
            #             if (i, j) not in arc_min_dia:
            #                 arc_min_dia[(i, j)] = d
            #             else:
            #                 arc_min_dia[(i, j)] = min(arc_min_dia[(i, j)], d)

            if self.z_star < self.current_cost:
                ampl_feas = AMPL()
                ampl_feas.reset()
                if self.data_number==5:
                    ampl_feas.read("newyork_model.mod")
                elif self.data_number==6:
                    ampl_feas.read("blacksburg_model.mod")
                else:
                    ampl_feas.read("wdnmodel_feasibility.mod")
                ampl_feas.read_data(self.data_file)

                for (x, y, k), val in self.l_star.items():
                    ampl_feas.eval(f'let l[{x},{y},{k}] := {val};')
                for (x, y), val in self.q_star.items():
                    ampl_feas.eval(f'let q[{x},{y}] := {val};')
                    if self.data_number ==5:
                        ampl_feas.eval(f'let q1[{x},{y}] := {self.q1_star[x,y]};')
                        ampl_feas.eval(f'let q2[{x},{y}] := {self.q2_star[x,y]};')
                for x, val in self.h_star.items():
                    ampl_feas.eval(f'let h[{x}] := {val};') 
                # ampl_feas.eval(f"""subject to con3_l_:sum {{(i,j) in arcs}} sum {{k in pipes}} C[k] * l[i,j,k] <= {self.current_cost};""")
                ampl_feas.eval("#minimize total_cost : sum{i in nodes diff Source} (-h[i] + E[i] + P[i])^2;")
                ampl_feas.option["solver"] = "ipopt"
                ampl_feas.set_option("ipopt_options", f"outlev = 0 max_iter = {self.max_iter} mu_init = {self.mu_init} expect_infeasible_problem = no  tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes")   #max_iter = 1000
                ampl_feas.option["presolve_eps"] = "6.82e-14"
                ampl_feas.option['presolve'] = 1
                with self.suppress_output():
                    ampl_feas.solve()
                solve_result = ampl_feas.solve_result
                solve_time = ampl_feas.get_value('_solve_elapsed_time')
                print("feasibility model| solve_result:",solve_result, "solve time:", solve_time)
                self.solver_time += solve_time
                self.number_of_nlp += 1
                self.l_star = ampl_feas.getVariable('l').getValues().to_dict()
                self.q_star = ampl_feas.getVariable('q').getValues().to_dict()
                self.h_star = ampl_feas.getVariable('h').getValues().to_dict() 

                ampl = AMPL()
                ampl.reset()
                if self.data_number==5:
                    ampl.read("newyork_model.mod")
                elif self.data_number==6:
                    ampl.read("blacksburg_model.mod")
                else:
                    ampl.read("wdnmodel.mod")
                ampl.read_data(self.data_file)

                for (x, y, k), val in self.l_star.items():
                    ampl.eval(f'let l[{x},{y},{k}] := {val};')
                for (x, y), val in self.q_star.items():
                    ampl.eval(f'let q[{x},{y}] := {val};')
                    if self.data_number ==5:
                        ampl.eval(f'let q1[{x},{y}] := {self.q1_star[x,y]};')
                        ampl.eval(f'let q2[{x},{y}] := {self.q2_star[x,y]};')
                for x, val in self.h_star.items():
                    ampl.eval(f'let h[{x}] := {val};') 
                    # ampl.eval(f'let h[{x}] := {self.E[x]+self.P[x] + 1e-2};') 

                # current_duals = {}
                # for con_name, val in ampl.get_constraints():
                #     dual_values = val.get_values()
                #     current_duals[con_name] = dual_values
                #
                # for con_name, dual_values in all_duals_reduced_model.items():
                #     if con_name in current_duals:
                #         # if con_name=="con7":
                #         ampl.get_constraint(con_name).set_values(dual_values) 

                if self.data_number==6:
                    # arcs = [(u,v) for (u,v) in self.arcs if (u,v) not in self.fixarcs]
                    ampl.eval("subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
                    # for (i,j) in arcs:
                    #     ampl.eval(f"s.t. arc_dia{i}_{j}: sum{{k in pipes: k >=  {arc_min_dia[i,j]}}} l[{i},{j},k] = L[{i},{j}];")
                else:
                    ampl.eval("subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
                    # for (i,j) in self.arcs:
                    #     ampl.eval(f"s.t. arc_dia{i}_{j}: sum{{k in pipes: k >=  {arc_min_dia[i,j]}}} l[{i},{j},k] = L[{i},{j}];")
 
                # with self.suppress_output():
                ampl.option["solver"] = "ipopt"
                ampl.set_option("ipopt_options", f"outlev = 0 max_iter = {self.max_iter} mu_init = {self.mu_init} expect_infeasible_problem = no  tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes")   #max_iter = 1000
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

                    if self.total_cost < self.current_cost - 1e-2:
                        # node_head_diff, arc_diff = self.compare_two_local_solutions(self.l, self.q, self.h, l, q, h)
                        # print(node_head_diff)
                        # node_set = set(node_head_diff.keys())
                        # arc_set = set(arc_diff.keys())
                        # no_change_arcs = [arc  for arc in self.arcs if arc not in arc_diff.keys()]
                        # print("no_change_arcs: ", no_change_arcs)
                        print(
                            f"{self.format_indian_number(round(self.current_cost)):<14}"
                            f"{self.format_indian_number(round(self.total_cost)):<14}"
                            f"{(str(round(self.ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                            f"{self.solve_result:<14}{'Yes':<10}"
                            f"{round(time.time() - self.start_time, 2)}s")
                        self.current_cost = self.total_cost
                        improved = True
                        self.ampl = ampl
                        # self.network_graph = self.generate_random_acyclic_from_solution(q)
                        self.l = l 
                        self.q = q
                        self.h = h 
                        # print(self.q)
                        # ampl.eval("display q;")
                        if self.data_number==5:
                            self.q1 = ampl.getVariable('q1').getValues().to_dict()
                            self.q2 = ampl.getVariable('q2').getValues().to_dict()
                        print("----------------------------------------------------------------------------------------")
                    else: 
                        print(
                            f"{self.format_indian_number(round(self.current_cost)):<14}"
                            f"{self.format_indian_number(round(self.total_cost)):<14}"
                            f"{(str(round(ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                            f"{self.solve_result:<14}{'No':<10}"
                            f"{round(time.time() - self.start_time, 2)}s")
                else:
                    print(#f"{str((i, j)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                        f"{self.solve_result:<14}{'No':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")

            abs_flows = sorted(
                abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4
            )
            m = len(abs_flows)
            # median_flow = abs_flows[m // 2]   # floor(m/2)
            median_flow = abs_flows[m//2]   # floor(m/2)
            # print("median_flow:", median_flow)
            if improved:
                # self.export_solution_iteration(self.iteration)
                # json_file = f"/home/nitishdumoliya/waterNetwork/wdnd/figure/json_file/d{self.data_number+1}/solution_{self.iteration}.json"
                # node_pos = node_position[self.data_number]
                # heuristic_approach = "Arc Reversal"
                # self.build_plot(self.iteration, json_file, node_pos, self.data_number, heuristic_approach, node_head_diff, arc_diff, edge)

                self.iteration += 1
                self.alpha = self.alpha_shrink * self.alpha
                # self.alpha = self.alpha_min
                self.Delta = self.alpha * median_flow
                self.Delta_h = self.Delta_h
                print(self.Delta)
                self.tr_failure_count = 0
                self.Terminate = False
                self.fail_streak = 0
                self.local_solution_improvement_heuristic()
            else:
                # Terminate = False
                self.alpha = self.alpha_expand * self.alpha
                self.Delta = self.alpha * median_flow
                self.Delta_h = self.Delta_h
                print(self.Delta)
                self.fail_streak += 1
                self.iteration += 1

                # Terminate2 = any(self.l[i,j,k] + self.eta > self.L[i,j]
                #                  for (i,j) in self.arcs for k in self.pipes
                # )
                if self.do_arc_reversal:
                    if self.fail_streak >= 3:
                        self.fail_streak = 0
                        # self.visited_arc_reverse = []
                        self.iterate_acyclic_flows()
                        return
                    else:
                        self.local_solution_improvement_heuristic()
                else:
                    if self.Terminate or self.fail_streak>=self.total_run:
                    # if self.Terminate:
                        print(self.fail_streak)
                        print("Trust-region exhausted → declaring local optimum.")
                        return
                    else:
                        self.Terminate = all(
                        abs(self.q[i, j]) + (self.Delta) > self.Q_max
                        for (i, j) in sorted_all_arcs
                        )
                        self.local_solution_improvement_heuristic()

    def iterate_acyclic_flows(self):
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

        # self.network_graph = self.best_acyclic_flow.copy()
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        # self.plot_graph(self.super_source_out_arc, self.current_cost, 0, self.q, self.h, self.D, (0,0), self.l, self.C)
        # print("\n*********************************************************************************************")
        print("Iteration :",self.iteration)
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
        sorted_all_arcs = [
            arc for arc in sorted_all_arcs
            if arc not in self.fix_arc_set
        ]
        
        sorted_arcs = [
            arc for arc in sorted_all_arcs
            if arc not in self.visited_arc_reverse
        ]
        
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
        for (i, j), dual_val in dual_dict.items():
            self.sen_score[(i, j)] = -dual_val * np.abs(self.h[i] - self.h[j])

        # print("sen_score:", self.sen_score)

        # --------------------------------------------------
        # Rank arcs using sensitivity score (absolute value)
        # --------------------------------------------------
        sorted_arcs = [
            arc for arc, _ in sorted(
                self.sen_score.items(),
                key=lambda kv: abs(kv[1]),
                reverse=True
            )
        ]
        # sorted_arcs = [
        #     arc for arc, _ in sorted(
        #         dual_dict.items(),
        #         key=lambda kv: abs(kv[1]),
        #         reverse=True
        #     )
        # ]       
        # print("sorted_arcs:", sorted_arcs)

        # remaining_arcs = [c for c, nc in sorted_by_abs_dual if nc not in self.visited_arc_reverse]
        # print("All cycles found:", len(sorted_by_abs_dual))
        # print("Visited cycles:", len(self.visited_arc_reverse))
        # print("Remaining new cycles:", len(remaining_arcs))

        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")
        print("----------------------------------------------------------------------------------------")

        # for edge in sorted_arcs[:min(20, len(sorted_arcs))]:
        for edge in sorted_arcs:
            self.visited_arc_reverse.append(edge)
            (i,j) = edge
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

            # for (i,j) in self.arcs:
            #     if (i,j)!=(u,v):
            #         if self.q[i,j]>=0:
            #             ampl.eval(f"s.t. flow_dir{i}_{j}: q[{i},{j}]>=0;")
            #         else:
            #             ampl.eval(f"s.t. flow_dir{i}_{j}: q[{i},{j}]<=0;")

            #self.ampl.eval(f"set inarc := {{{inarc_set}}};")
            #self.ampl.eval(f"set indegree_node := {{{set(self.indegree_2_or_more)}}};")
            # fix_arc_set = self.fix_leaf_arc_flow(ampl)

            # for (u,v) in fix_arc_set:
            #     ampl.eval(f"s.t. fix_arc_flow_{u}_{v}: q[{u}, {v}] = {self.q[u,v]};")
            # self.update_initial_points(self.l, self.q, self.h)
            # self.update_initial_points_with_perturbation(self.ampl, self.l, self.q, self.h) 

            if self.q[i,j] >= 0:
                ampl.eval(f"s.t. flow_direction1{i}_{j}: q[{i}, {j}]<=0;")
            else:
                ampl.eval(f"s.t. flow_direction1{i}_{j}: q[{i}, {j}]>=0;")

            if self.data_number==6:
                # arcs = [(u,v) for (u,v) in self.arcs if (u,v) not in self.fixarcs]
                ampl.eval("subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
                # for (i,j) in arcs:
                #     ampl.eval(f"s.t. arc_dia{i}_{j}: sum{{k in pipes: k >=  {arc_min_dia[i,j]}}} l[{i},{j},k] = L[{i},{j}];")
            else:
                ampl.eval("subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];")

            # with self.suppress_output():
            ampl.option["solver"] = "ipopt"
            ampl.set_option("ipopt_options", f"outlev = 0 max_iter = {self.max_iter} mu_init = {self.mu_init} expect_infeasible_problem = no  tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes")   #max_iter = 1000
            # ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes mu_strategy = adaptive recalc_y = no"
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
                    node_head_diff, arc_diff = self.compare_two_local_solutions(self.l, self.q, self.h, l, q, h)
                    # print(node_head_diff)
                    # node_set = set(node_head_diff.keys())
                    # arc_set = set(arc_diff.keys())
                    # no_change_arcs = [arc  for arc in self.arcs if arc not in arc_diff.keys()]
                    # print("no_change_arcs: ", no_change_arcs)
                    print(f"{str((i, j)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(self.ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                        f"{self.solve_result:<14}{'Yes':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    #print("\n")
                    # self.plot_graph(self.super_source_out_arc, self.total_cost, 0, q, h, self.D, (0,0), l, self.C)
                    self.current_cost = self.total_cost
                    improved = True
                    self.is_improved_in_arc_reversal = True
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
                    # self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set)
                    print("----------------------------------------------------------------------------------------")
                    print("-------------------Adaptive Two-Level Neighborhood Search (ATLNS)----------------------------")
                    # self.iteration = self.iteration + 1
                    self.l_points = []
                    self.q_points = []
                    self.z_star = 0
                    self.l_star = self.l
                    self.q_star = self.q
                    self.h_star = self.h
                    self.alpha = self.alpha_min
                    abs_flows = sorted(abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4)
                    m = len(abs_flows)
                    # self.Delta = self.alpha_min*abs_flows[m // 2]           
                    self.Delta = self.alpha*abs_flows[m//2]     
                    self.local_solution_improvement_heuristic()

                else: 
                    print(f"{str((i, j)):<10}"
                        f"{self.format_indian_number(round(self.current_cost)):<14}"
                        f"{self.format_indian_number(round(self.total_cost)):<14}"
                        f"{(str(round(ampl.get_value('_solve_elapsed_time'), 2)) + 's'):<12}"
                        f"{self.solve_result:<14}{'No':<10}"
                        f"{round(time.time() - self.start_time, 2)}s")
                    # ampl.eval("display q;")
            else:
                # self.plot_graph(self.super_source_out_arc, total_cost, 0, q1, h1, self.D, (0,0), l1, self.C)
                print(f"{str((i, j)):<10}"
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
                # self.export_solution_iteration(self.iteration)
                # json_file = f"/home/nitishdumoliya/waterNetwork/wdnd/figure/json_file/d{self.data_number+1}/solution_{self.iteration}.json"
                # node_pos = node_position[self.data_number]
                # heuristic_approach = "Arc Reversal"
                # self.build_plot(self.iteration, json_file, node_pos, self.data_number, heuristic_approach, node_head_diff, arc_diff, edge)

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

    def solve(self):
        self.ampl.option["solver"] = "ipopt"
        self.ampl.set_option("ipopt_options", f"outlev = 5 tol = 1e-9 bound_relax_factor=0  bound_push = {self.bound_push} bound_frac = {self.bound_frac} halt_on_ampl_error = yes warm_start_init_point = no expect_infeasible_problem = no")   #max_iter = 1000
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

    def run(self):
        """Main function to run the Heuristic Approach."""
        self.start_time = time.time()

        # print("NLP solve using:  smooth approximation 1, Epsilon selection using absolute error\n")
        # print("NLP solve using: smooth approximation 1, epsilon selection using relative error\n")
        # print("NLP solve using: smooth approximation 2, epsilon selection using absolute error\n")
        print("NLP Model: Smooth Approximatate WDN Design Model 2, Epsilon Selection using Relative Error")
        print("NLP Solver: Ipopt")
        print("************************Initial Solution of Approximate WDN Design Model**********************")
        self.load_model()
        fix_arc_set = self.fix_leaf_arc_flow(self.ampl)
        # for (u,v) in fix_arc_set:
        #     self.ampl.eval(f"s.t. fix_arc_flow_{u}_{v}: q[{u}, {v}] = {self.q[u,v]};")
        print("fix_arc_set:",fix_arc_set)
        self.super_source_out_arc = self.fix_arc_set()
        print("super_source_out_arc:", self.super_source_out_arc, "\n")

        if self.data_number==6:
            self.ampl.eval("subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
        else:
            self.ampl.eval("subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];")

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
        l_initial = self.l
        q_initial = self.q
        h_initial = self.h
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
        self.iteration = 0
        self.elapsed_time = time.time() - self.start_time
        # [os.remove(f) for f in os.listdir(f"../figure/json_file/d{self.data_number + 1}/") if f.startswith("wdn_interactive")]
        # self.export_solution(0)
        # # fig = self.(f"../figure/json_file/solution_{self.data_number}.json",node_position[self.data_number])
        # json_file = f"/home/nitishdumoliya/waterNetwork/wdnd/figure/json_file/d{self.data_number+1}/solution_{self.iteration}.json"
        # node_pos = node_position[self.data_number]
        # heuristic_approach = "Original Model"
        # self.build_plot(self.iteration, json_file, node_pos, self.data_number, heuristic_approach, {}, {}, (0,0))

        print("---------------------------Reverse Arc Direction Approach------------------------------------")
        self.iteration = self.iteration + 1
        self.visited_arc_reverse = []
        # self.iterate_acyclic_flows_new() 
        self.do_arc_reversal = True
        self.is_improved_in_arc_reversal = False
        self.iterate_acyclic_flows() 

        if self.is_improved_in_arc_reversal == False:
            print("-------------------Adaptive Two-Level Neighborhood Search (ATLNS)----------------------------")
            self.iteration = self.iteration + 1
            self.l_points = []
            self.q_points = []
            self.z_star = 0
            self.l_star = self.l
            self.q_star = self.q
            self.h_star = self.h
            self.do_arc_reversal = False
            self.alpha = self.alpha_min
            abs_flows = sorted(abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4)
            m = len(abs_flows)
            # self.Delta = self.alpha_min*abs_flows[m // 2]           
            self.Delta = self.alpha*abs_flows[m//2]
            self.alpha_min = 0.0001
            self.alpha_shrink = 0.1
            self.alpha_expand = 1.4 
            self.total_run = 5
            self.local_solution_improvement_heuristic()

        # print("local points:", self.local_points)
        # self.export_solution_iteration(self.iteration)
        # json_file = f"/home/nitishdumoliya/waterNetwork/wdnd/figure/json_file/d{self.data_number + 1}/solution_{self.iteration}.json"
        # node_head_diff, arc_diff = self.compare_two_local_solutions(l_initial, q_initial, h_initial, self.l, self.q, self.h)
        #
        # node_pos = node_position[self.data_number]
        # heuristic_approach = "Reverse Arc"
        # self.build_plot(2, json_file, node_pos, self.data_number, heuristic_approach, node_head_diff, arc_diff, (0,0))

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
        print("Total number of iteration:", self.iteration)
        # self.constraint_violations(self.q, self.h, self.l, self.eps, "ipopt")
        self.elapsed_time = time.time() - self.start_time
        solver_time = self.solver_time
        print(f"Solver_time: {solver_time:.2f} seconds")
        # print(f"Heuristic elapsed time: {elapsed_time:.2f} seconds = {elapsed_time/60:.2f} minutes.\n")
        print(f"Heuristic elapsed time: {self.elapsed_time:.2f} seconds")
        # self.export_solution(0)

        print("***********************************************************************************************\n")
        # self.launch_dashboard(fig)

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
    # --------------------------------------------------
