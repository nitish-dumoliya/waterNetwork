"""
Water Network Design Optimizer
================================
Heuristic approach for optimizing water distribution network (WDN) pipe design
using NLP (IPOPT), arc reversal, and local improvement strategies.

Dependencies:
    amplpy, networkx, numpy, scipy, matplotlib, plotly, dash, sklearn,
    optuna, pyswarm, jax, tabulate
"""

# ── Standard Library ──────────────────────────────────────────────────────────
import contextlib
import copy
import json
import math
import os
import random
import sys
import time
import warnings

warnings.filterwarnings("ignore")

# ── Third-Party ───────────────────────────────────────────────────────────────
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

from amplpy import AMPL
from collections import Counter
from networkx.algorithms import cycle_basis
from sklearn.decomposition import PCA
from scipy.integrate import solve_ivp
from scipy.linalg import svd, eigh
from numpy.linalg import eig
from tabulate import tabulate

import jax
import jax.numpy as jnp
from jax import jacfwd, jacrev

# ── Local ─────────────────────────────────────────────────────────────────────
from network_layout import node_position


# ══════════════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════════════

DATA_LIST = [
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
    "d14_balerma",
]

# IPOPT / solver defaults
DEFAULT_TOL = 1e-9
DEFAULT_MAX_ITER = 2000
DEFAULT_BOUND_PUSH = 1e-1
DEFAULT_BOUND_FRAC = 1e-1
DEFAULT_MU_INIT = 1e-2
PRESOLVE_EPS_MAIN = "6.82e-14"
PRESOLVE_EPS_LP = "7.19e-7"


# ══════════════════════════════════════════════════════════════════════════════
# WaterNetworkOptimizer
# ══════════════════════════════════════════════════════════════════════════════

class WaterNetworkOptimizer:
    """
    Adaptive Variable Neighborhood Search (AVNS) heuristic for water
    distribution network (WDN) pipe-sizing optimization.

    The optimizer uses a smooth NLP approximation (IPOPT) of the Hazen-Williams
    head-loss equations, combined with:
      - Local solution improvement (trust-region-style NLP re-solves)
      - Arc-reversal perturbation (changes flow direction on selected arcs)
      - Diameter-reduction perturbation (forces smaller diameters on arcs)

    Parameters
    ----------
    model_file : str
        Path to the AMPL model file (.mod).
    data_file : str
        Path to the AMPL data file (.dat).
    data_number : int
        Zero-based index into DATA_LIST identifying which network is being solved.
    data_list : list[str]
        Names of all available networks (used for labelling outputs).
    """

    # ── Construction ──────────────────────────────────────────────────────────

    def __init__(self, model_file: str, data_file: str, data_number: int, data_list: list):
        self.ampl = AMPL()
        self.model_file = model_file
        self.data_file = data_file
        self.data_number = data_number
        self.data_list = data_list

        # Solution state
        self.total_cost = None
        self.network_graph = None
        self.solve_result = None
        self.best_acyclic_flow = None
        self.current_cost = None

        # Solver tracking
        self.solver_time = 0
        self.number_of_nlp = 0
        self.start_time = None

        # IPOPT options
        self.tol = DEFAULT_TOL
        self.max_iter = DEFAULT_MAX_ITER
        self.bound_push = DEFAULT_BOUND_PUSH
        self.bound_frac = DEFAULT_BOUND_FRAC
        self.mu_init = DEFAULT_MU_INIT

        # Trust-region / neighbourhood parameters
        self.alpha_shrink = 0.1
        self.alpha_q_min = 0.2
        self.alpha_expand = 1.2
        self.alpha_q = self.alpha_q_min

        self.eta_l_min = 0.2
        self.eta_h_min = 0.2
        self.eta_l_expend = 1.5
        self.eta_h_expend = 1.5
        self.eta_l = self.eta_l_min
        self.eta_h = self.eta_h_min

        # Termination / iteration tracking
        self.total_run = 10
        self.max_tr_failures = 50
        self.tr_failure_count = 0
        self.Terminate = False
        self.success_streak = 0
        self.fail_streak = 0
        self.k_neigh = 1

        # Hazen-Williams constants
        self.delta = 0.1
        self.p = 1.852
        self.omega = 10.67

    # ── Model Loading ─────────────────────────────────────────────────────────

    def load_model(self):
        """Load the AMPL model and data files and cache commonly-used sets/parameters."""
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)

        self.nodes = self.ampl.getSet("nodes")
        self.source = self.ampl.getSet("Source")
        self.arcs = self.ampl.getSet("arcs")
        self.pipes = self.ampl.getSet("pipes")

        if self.data_number == 6:
            self.fixarcs = self.ampl.getSet("fixarcs")

        self.L = self.ampl.getParameter("L").to_dict()
        self.D = self.ampl.getParameter("D").to_dict()
        self.C = self.ampl.getParameter("C").to_dict()
        self.P = self.ampl.getParameter("P").to_dict()
        self.R = self.ampl.getParameter("R").to_dict()
        self.E = self.ampl.getParameter("E").to_dict()
        self.d = self.ampl.getParameter("d").to_dict()
        self.Q_max = self.ampl.getParameter("Q_max").to_list()[0]
        self.eps = self.ampl.getParameter("eps").to_dict()
 
        self.alpha = self.ampl.get_parameter("alpha").get_values().to_dict()
        self.R_min = min(self.R[k] for k in self.pipes)
        self.R_max = max(self.R[k] for k in self.pipes)
        self.d_min = min(self.d[k] for k in self.pipes)
        self.d_max = max(self.d[k] for k in self.pipes)

    # ── Graph Utilities ───────────────────────────────────────────────────────

    def create_digraph(self):
        """Build a NetworkX DiGraph from the AMPL node/arc sets."""
        nodes_list = list(self.ampl.getSet("nodes"))
        edges_list = self.ampl.getSet("arcs").to_list()

        self.network_graph = nx.DiGraph()
        self.network_graph.add_nodes_from(nodes_list)
        self.network_graph.add_edges_from(edges_list)

        print(nodes_list)
        print(edges_list)

    def generate_random_acyclic_from_solution(self, q: dict) -> nx.DiGraph:
        """
        Build a directed acyclic graph whose edge directions are determined
        by the sign of the flow variables *q*.

        Parameters
        ----------
        q : dict
            Flow values keyed by arc tuple (i, j).

        Returns
        -------
        nx.DiGraph
        """
        self.network_graph = nx.DiGraph()
        self.network_graph.add_nodes_from(self.nodes)

        for (i, j) in self.arcs:
            if q[i, j] >= 0:
                self.network_graph.add_edge(i, j)
            else:
                self.network_graph.add_edge(j, i)

        return self.network_graph

    def generate_random_acyclic_graph(self):
        """
        Generate a random acyclic directed graph by first building a random
        spanning tree (Wilson's algorithm) and then orienting additional arcs
        to maintain acyclicity.
        """
        uwg = nx.Graph()
        nodes_list = list(self.ampl.getSet("nodes"))
        edges_list = self.ampl.getSet("arcs").to_list()
        uwg.add_nodes_from(nodes_list)
        uwg.add_edges_from(edges_list)

        random_tree = nx.random_spanning_tree(uwg)
        root_l = self.ampl.getSet("Source").to_list()
        root = root_l[0]

        if root not in random_tree.nodes:
            raise ValueError("The specified root must be a node in the graph.")

        self.network_graph = nx.DiGraph()
        visited = set()

        def dfs(node):
            visited.add(node)
            for neighbor in random_tree.neighbors(node):
                if neighbor not in visited:
                    self.network_graph.add_edge(node, neighbor)
                    dfs(neighbor)

        dfs(root)

        # Visualise the spanning tree
        plt.figure(figsize=(15, 10))
        plt.subplot(121)
        nx.draw_spectral(
            self.network_graph, with_labels=True,
            node_color="lightgreen", font_weight="bold", arrows=True,
        )
        plt.title("Directed Spanning Tree")

        # Add remaining edges while preserving acyclicity
        for u, v in uwg.edges():
            if not self.network_graph.has_edge(u, v):
                self.network_graph.add_edge(u, v)
                if not nx.is_directed_acyclic_graph(self.network_graph):
                    self.network_graph.remove_edge(u, v)
                    self.network_graph.add_edge(v, u)

        plt.subplot(122)
        nx.draw_spectral(
            self.network_graph, with_labels=True,
            node_color="lightgreen", font_weight="bold", arrows=True,
        )
        plt.title("Acyclic Directed Graph")
        plt.show()

    def cycle_basis(self):
        """Print the cycle basis of the undirected network graph."""
        root = self.ampl.getSet("Source").to_list()[0]
        nodes_list = list(self.ampl.getSet("nodes"))
        edges_list = self.ampl.getSet("arcs").to_list()

        uwg = nx.Graph()
        uwg.add_nodes_from(nodes_list)
        uwg.add_edges_from(edges_list)

        print("cycle basis for given water network:", nx.cycle_basis(uwg, root))

    # ── Feasibility Checks ────────────────────────────────────────────────────

    def is_valid_edge(self, source: int, target: int) -> bool:
        """Return True if adding edge (source → target) keeps the graph acyclic."""
        self.network_graph.add_edge(source, target)
        is_dag = nx.is_directed_acyclic_graph(self.network_graph)
        self.network_graph.remove_edge(source, target)
        return is_dag

    def check_incoming_arcs(self) -> bool:
        """
        Return True iff every non-source node has at least one incoming arc.
        """
        root = list(self.ampl.getSet("Source"))[0]
        for node in self.network_graph.nodes():
            if node != root and self.network_graph.in_degree(node) < 1:
                return False
        return True

    # ── Presolve Utilities ────────────────────────────────────────────────────

    def _is_cycle(
        self,
        graph: nx.Graph,
        start_node: int,
        end_node: int,
        visited_copy: dict,
        parent: int,
    ) -> bool:
        """Recursive DFS helper to detect whether *end_node* lies on a cycle."""
        visited_copy[start_node] = True
        for neighbor in graph.neighbors(start_node):
            if not visited_copy[neighbor]:
                if self._is_cycle(graph, neighbor, end_node, visited_copy, start_node):
                    return True
            elif parent != neighbor and end_node == neighbor:
                return True
        return False

    def _presolve(
        self,
        graph: nx.Graph,
        node: int,
        visited: dict,
        parent: int,
        set_arc: list,
    ) -> list:
        """
        DFS traversal that records which arcs have a fixed direction
        (those that lie on a cycle from *node* back to itself).
        """
        visited_copy = visited.copy()
        in_cycle = self._is_cycle(graph, node, node, visited_copy, parent)
        visited[node] = True

        if in_cycle:
            for neighbor in graph.neighbors(node):
                if parent != neighbor:
                    set_arc.append((node, neighbor))
        else:
            for neighbor in graph.neighbors(node):
                if parent != neighbor:
                    set_arc.append((node, neighbor))
                    self._presolve(graph, neighbor, visited, node, set_arc)

        return set_arc

    def fix_arc_set(self) -> list:
        """
        Return the list of arcs with a topologically fixed flow direction
        (arcs that emanate from the source through a tree traversal).
        """
        graph = nx.Graph()
        arc_set = self.ampl.getSet("arcs").to_list()
        graph.add_edges_from(arc_set)

        visited = {node: False for node in graph.nodes()}
        source = self.ampl.getSet("Source").to_list()[0]
        set_arc = []
        return self._presolve(graph, source, visited, -1, set_arc)

    def fix_leaf_arc_flow(self, ampl: AMPL) -> set:
        """
        Iteratively fix flow values for leaf arcs (degree-1 nodes) by propagating
        demand balances inward from the network periphery.

        Returns
        -------
        set
            Set of arcs (i, j) whose flows have been fixed.
        """
        graph = nx.Graph()
        arc_set = ampl.getSet("arcs").to_list() if hasattr(ampl, "getSet") else self.ampl.getSet("arcs").to_list()
        graph.add_edges_from(arc_set)

        D = self.ampl.getParameter("D").getValues().to_dict()
        source = self.ampl.getSet("Source").to_list()[0]
        fixed_arcs: set = set()
        Q_max = self.ampl.getParameter("Q_max").getValues().to_list()[0]
        D[source] = -Q_max

        while True:
            leaf_nodes = [n for n in graph.nodes if graph.degree[n] == 1]
            if not leaf_nodes:
                break

            for leaf in leaf_nodes:
                neighbor = next(graph.neighbors(leaf))

                if (neighbor, leaf) in arc_set:
                    edge = (neighbor, leaf)
                    if edge not in fixed_arcs:
                        if leaf == source:
                            D[neighbor] += D[leaf]
                            source = neighbor
                        else:
                            D[neighbor] += D[leaf]
                        fixed_arcs.add(edge)
                else:
                    edge = (leaf, neighbor)
                    if edge not in fixed_arcs:
                        flow_value = -D[leaf]
                        if leaf == source:
                            D[neighbor] -= flow_value
                            source = neighbor
                        elif neighbor == source:
                            D[neighbor] -= D[leaf]
                        else:
                            D[neighbor] += -flow_value
                        fixed_arcs.add(edge)

                graph.remove_node(leaf)

        return fixed_arcs

    # fix_leaf_arc_flow1 is an identical variant kept for internal use
    fix_leaf_arc_flow1 = fix_leaf_arc_flow

    # ── Acyclic Arc Utilities ─────────────────────────────────────────────────

    def acyclic_arcs(self) -> set:
        """
        For nodes with in-degree ≥ 2, try reversing each incoming arc and
        check whether the resulting graph is still acyclic with valid in-degrees.

        Returns
        -------
        set
            Arcs that can be reversed while maintaining acyclicity.
        """
        network_graph = self.best_acyclic_flow
        indegree_2_or_more = [
            node for node, indeg in network_graph.in_degree() if indeg >= 2
        ]
        acyclic_arc = set()

        for node in indegree_2_or_more:
            for (u, v) in list(network_graph.in_edges(node)):
                if (u, v) not in self.super_source_out_arc:
                    network_graph.remove_edge(u, v)
                    network_graph.add_edge(v, u)

                    if nx.is_directed_acyclic_graph(network_graph) and self.check_incoming_arcs():
                        acyclic_arc.add((u, v))

                    network_graph.remove_edge(v, u)
                    network_graph.add_edge(u, v)

        return acyclic_arc

    # ── Formatting ────────────────────────────────────────────────────────────

    @staticmethod
    def format_indian_number(num: int) -> str:
        """
        Format an integer using the Indian numbering convention
        (e.g., 1,23,45,678).
        """
        s = str(num)
        if len(s) <= 3:
            return s
        last_three = s[-3:]
        remaining = s[:-3]
        groups = [
            remaining[max(i - 2, 0): i]
            for i in range(len(remaining), 0, -2)
        ][::-1]
        return ",".join(groups) + "," + last_three

    # ── Output Suppression ────────────────────────────────────────────────────

    @contextlib.contextmanager
    def suppress_output(self):
        """Context manager that silences stdout and stderr."""
        with open(os.devnull, "w") as devnull:
            old_stdout, old_stderr = sys.stdout, sys.stderr
            sys.stdout = sys.stderr = devnull
            try:
                yield
            finally:
                sys.stdout, sys.stderr = old_stdout, old_stderr

    # ── Solution I/O ──────────────────────────────────────────────────────────

    def export_solution(self, iteration: int):
        """Serialise the current best solution to a JSON file."""
        solution = self._build_solution_dict(iteration, time_elapsed=self.elapsed_time)
        path = f"../figure/json_file/d{self.data_number + 1}/solution_{iteration}.json"
        with open(path, "w") as f:
            json.dump(solution, f, indent=2)

    def export_solution_iteration(self, iteration: int):
        """Serialise the current solution (timestamped from start) to a JSON file."""
        solution = self._build_solution_dict(
            iteration, time_elapsed=time.time() - self.start_time
        )
        path = f"../figure/json_file/d{self.data_number + 1}/solution_{self.iteration}.json"
        with open(path, "w") as f:
            json.dump(solution, f, indent=2)

    def _build_solution_dict(self, iteration: int, time_elapsed: float) -> dict:
        """Return a JSON-serialisable dict of the current solution."""
        return {
            "iteration": iteration,
            "objective": self.current_cost,
            "solve_time": time_elapsed,
            "q": {f"{i},{j}": float(v) for (i, j), v in self.q.items()},
            "h": self.h,
            "l": {
                f"{i},{j},{k}": self.l[i, j, k]
                for (i, j) in self.arcs
                for k in self.pipes
                if self.l[i, j, k] > 1e-3
            },
            "D": self.D,
            "Diameters": {f"{i}": self.d[i] for i in self.pipes},
            "L": {f"{i},{j}": v for (i, j), v in self.L.items()},
            "C": self.C,
            "hmin": {
                f"{i}": self.E[i] if i == list(self.source)[0] else self.E[i] + self.P[i]
                for i in self.nodes
            },
            "nodes": list(self.nodes),
            "arcs": [(i, j) for (i, j) in self.arcs],
            "source": list(self.source),
        }

    # ── Visualisation ─────────────────────────────────────────────────────────

    def display_results(self):
        """Print key AMPL decision variables and objective to stdout."""
        self.ampl.eval("display {(i,j) in arcs, k in pipes:l[i,j,k]>1} l[i,j,k];")
        self.ampl.eval("display {(i,j) in arcs}: q[i,j];")
        self.ampl.eval("display h;")
        self.ampl.eval("display solve_result;")
        self.total_cost = self.ampl.get_objective("total_cost").value()
        print("Objective:", self.total_cost)

    def plot_graph(
        self,
        super_source_out_arc=None,
        current_cost=None,
        iteration=1,
        edge_weights=None,
        h=None,
        D=None,
        arc=(0, 0),
        l=None,
        C=None,
    ):
        """Render the network graph using Matplotlib."""
        l = l or {}
        C = C or {}
        pos = node_position[self.data_number]

        # Per-arc pipe cost
        cost = {
            (i, j): sum(l[i, j, k] * C[k] for k in self.pipes)
            for (i, j) in self.ampl.getSet("arcs")
        }

        plt.figure(figsize=(10, 8))
        cmap = plt.cm.plasma
        network_graph = self.generate_random_acyclic_from_solution(edge_weights)

        # Nodes
        nx.draw_networkx_nodes(
            network_graph, pos, node_color="lightblue",
            edgecolors="black", node_size=300, linewidths=0.5,
            label="Regular Nodes",
        )
        indegree_2_or_more = [
            node for node, indeg in network_graph.in_degree() if indeg >= 2
        ]
        if indegree_2_or_more:
            nx.draw_networkx_nodes(
                network_graph, pos, nodelist=indegree_2_or_more,
                node_color="orange", edgecolors="orange", node_size=300,
                label="Nodes with In-Degree ≥ 2",
            )

        nx.draw_networkx_labels(network_graph, pos, font_size=10)
        nx.draw_networkx_edges(
            network_graph, pos, arrowstyle="->", arrowsize=12,
            edge_color="black", label="Regular Arcs", arrows=True,
        )
        if super_source_out_arc:
            nx.draw_networkx_edges(
                network_graph, pos, edgelist=super_source_out_arc,
                arrowstyle="->", arrowsize=12, edge_color="red",
                width=1, label="Fix arc direction",
            )

        # Node demand annotations
        if D:
            for node, (x, y) in pos.items():
                demand = D.get(node, 0)
                plt.text(x - 80, y - 130, f"{demand:.2f}", fontsize=10,
                         color="magenta", ha="center")

        # Edge flow annotations
        if edge_weights:
            for (u, v), weight in edge_weights.items():
                mid_x = (pos[u][0] + pos[v][0]) / 2
                mid_y = (pos[u][1] + pos[v][1]) / 2
                plt.text(mid_x + 150, mid_y + 80, f"{weight * 1000:.2f}",
                         fontsize=10, color="black", ha="center")

        # Pipe diameter annotations
        pipe_dia_arc = {
            (i, j): [k for k in self.pipes if l.get((i, j, k), 0) >= 1e-5]
            for (i, j) in self.arcs
        }
        for (u, v), diameters in pipe_dia_arc.items():
            mid_x = (pos[u][0] + pos[v][0]) / 2
            mid_y = (pos[u][1] + pos[v][1]) / 2
            plt.text(mid_x - 100, mid_y - 100, f"{diameters}",
                     fontsize=10, color="purple", ha="center")

        plt.title(f"Total cost: {self.format_indian_number(round(current_cost))}")
        plt.savefig(
            f"/home/nitishdumoliya/waterNetwork/wdnd/figure/newfigure/"
            f"d{self.data_number}_iteration_{iteration}.png"
        )
        plt.box(False)
        plt.show()

    # ── Constraint Violation Reporting ────────────────────────────────────────

    def constraint_violations(
        self,
        q_values: dict,
        h_values: dict,
        l_values: dict,
        epsilon: dict,
        solver: str,
    ):
        """
        Compute and print a detailed constraint violation report for the
        current NLP solution.

        Checks:
          con1 – Flow balance at demand nodes
          con2 – Hazen-Williams head-loss (original vs. smooth approximation)
          con3 – Pipe-length equality
          con4 – Individual pipe-length upper bounds
          con5 – Source head equality
          con6 – Minimum pressure at demand nodes
        """
        total_abs_viol = 0.0

        # ── con1: Flow balance ────────────────────────────────────────────────
        con1_gap = {}
        flow_vars = (
            {**{k: self.q1[k] + self.q2[k] for k in self.q1}, **{}}
            if self.data_number == 5
            else q_values
        )

        for i in self.nodes:
            if i not in self.source:
                rhs = self.D[i]
                inflow = sum(q_values[j, i] for j in self.nodes if (j, i) in self.arcs)
                outflow = sum(q_values[i, j] for j in self.nodes if (i, j) in self.arcs)
                viol = inflow - outflow - rhs
                con1_gap[str(i)] = viol
                total_abs_viol += abs(viol)

        # ── con2: Head-loss ───────────────────────────────────────────────────
        con2_original_gap: dict = {}
        con2_approx_gap: dict = {}
        absolute_violations: dict = {}
        relative_violations: dict = {}

        con2_original_viol = con2_approx_viol = 0.0
        con2_abs_diff = con2_rel_diff = 0.0

        arc_iterator = (
            self.arcs
            if self.data_number not in (5, 6)
            else (
                [(i, j) for (i, j) in self.arcs if (i, j) not in self.fixarcs]
                if self.data_number == 6
                else list(q_values.keys())
            )
        )

        for (i, j) in arc_iterator:
            lhs = h_values[i] - h_values[j]

            if self.data_number == 6:
                alpha_rhs = sum(
                    10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * (self.d[k] ** 4.87))
                    for k in self.pipes
                )
            elif self.data_number == 5:
                alpha_rhs = sum(
                    10.67 * l_values[i, j, k] / ((self.R[i, j] ** 1.852) * (self.d[k] ** 4.87))
                    for k in self.pipes
                )
            else:
                alpha_rhs = sum(
                    10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * (self.d[k] ** 4.87))
                    for k in self.pipes
                )

            q_ij = q_values[i, j]
            original_rhs = q_ij * abs(q_ij) ** 0.852 * alpha_rhs
            eps_ij = epsilon[i, j]
            approx_rhs = (
                q_ij ** 3
                * ((q_ij ** 2 + eps_ij ** 2) ** 0.426)
                / (q_ij ** 2 + 0.426 * eps_ij ** 2)
            ) * alpha_rhs

            key = f"{i},{j}"
            con2_original_gap[key] = lhs - original_rhs
            con2_approx_gap[key] = lhs - approx_rhs
            total_abs_viol += abs(lhs - approx_rhs)
            con2_original_viol += abs(lhs - original_rhs)
            con2_approx_viol += abs(lhs - approx_rhs)

            abs_diff = original_rhs - approx_rhs
            absolute_violations[key] = abs_diff
            con2_abs_diff += abs(abs_diff)

            rel_diff = abs_diff / (original_rhs + 1e-14)
            relative_violations[key] = rel_diff
            con2_rel_diff += abs(rel_diff)

        # ── con3: Pipe-length sum ─────────────────────────────────────────────
        con3_gap = {}
        for (i, j) in self.arcs:
            if self.data_number == 6 and (i, j) in self.fixarcs:
                continue
            viol = sum(l_values[i, j, k] for k in self.pipes) - self.L[i, j]
            con3_gap[f"{i},{j}"] = viol
            total_abs_viol += abs(viol)

        # ── con4: Individual pipe-length upper bound ──────────────────────────
        con4_gap = {}
        for (i, j) in self.arcs:
            for k in self.pipes:
                viol = max(0, l_values[i, j, k] - self.L[i, j])
                con4_gap[f"{i},{j},{k}"] = viol
                total_abs_viol += abs(viol)

        # ── con5: Source head ─────────────────────────────────────────────────
        con5_gap = {}
        for j in self.source:
            viol = h_values[j] - self.E[j]
            con5_gap[str(j)] = viol
            total_abs_viol += abs(viol)

        # ── con6: Minimum pressure ────────────────────────────────────────────
        con6_gap = {}
        for j in self.nodes:
            if j not in self.source:
                viol = max(0, -h_values[j] + self.E[j] + self.P[j])
                con6_gap[str(j)] = viol
                total_abs_viol += abs(viol)

        # ── Print summary ─────────────────────────────────────────────────────
        print("Sum of constraints violation:", total_abs_viol)
        print("Sum of violation of original headloss constraint:", con2_original_viol)
        print("Sum of violation of approx headloss constraint:", con2_approx_viol)
        print("Con2 sum of absolute diff (original vs approx):", con2_abs_diff)
        print("Con2 sum of relative diff (original vs approx):", con2_rel_diff)

    # ── Solution Comparison ───────────────────────────────────────────────────

    def compare_two_local_solutions(
        self,
        l1: dict, q1: dict, h1: dict,
        l2: dict, q2: dict, h2: dict,
        tol: float = 1e-5,
    ) -> tuple:
        """
        Compare two locally optimal solutions.

        Returns
        -------
        node_head_diff : dict
            Nodes where head values differ beyond *tol*.
        arc_diff : dict
            Arcs where flow or pipe composition differ beyond *tol*.
        """
        # Node head differences
        node_head_diff = {
            i: {
                "h_sol1": round(h1[i], 4),
                "h_sol2": round(h2[i], 4),
                "delta": round(h2[i] - h1[i], 4),
            }
            for i in h1
            if abs(h1.get(i, 0.0) - h2.get(i, 0.0)) > tol
        }

        # Arc differences
        arc_diff = {}
        for (i, j) in set(q1) | set(q2):
            info = {}
            q1_ij, q2_ij = q1.get((i, j), 0.0), q2.get((i, j), 0.0)

            if abs(q1_ij - q2_ij) > tol:
                info["flow"] = {
                    "q_sol1": round(q1_ij, 4),
                    "q_sol2": round(q2_ij, 4),
                    "delta": round(q2_ij - q1_ij, 4),
                }

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
                info["pipes"] = {"sol1": pipes1, "sol2": pipes2}

            if info:
                arc_diff[(i, j)] = info

        return node_head_diff, arc_diff

    # ── Helper: dict addition ─────────────────────────────────────────────────

    @staticmethod
    def add_dict(d1: dict, d2: dict) -> dict:
        """Return a new dict that is the element-wise sum of *d1* and *d2*."""
        result = d1.copy()
        for k, v in d2.items():
            result[k] = result.get(k, 0.0) + v
        return result

    # ── Helper: arc max diameter ──────────────────────────────────────────────

    def _arc_max_diameter(self) -> dict:
        """Return a dict mapping each arc (i,j) → largest active pipe diameter."""
        arc_max_dia = {}
        fixarcs = set()
        if self.data_number == 6:
            fixarcs = set(self.ampl.getSet("fixarcs"))

        for (i, j, d), val in self.l.items():
            if self.data_number == 6 and ((i, j) in fixarcs or (j, i) in fixarcs):
                continue
            if val > 1e-3:
                if (i, j) not in arc_max_dia:
                    arc_max_dia[(i, j)] = d
                else:
                    arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)

        return arc_max_dia

    # ── AMPL Sub-model Builders ───────────────────────────────────────────────

    def _build_base_ampl(self) -> AMPL:
        """
        Create and return a fresh AMPL instance loaded with the appropriate
        model and data files for the current network.
        """
        ampl = AMPL()
        ampl.reset()

        model_map = {
            5: "newyork_model.mod",
            6: "blacksburg_model.mod",
        }
        ampl.read(model_map.get(self.data_number, "wdnmodel.mod"))
        ampl.read_data(self.data_file)
        return ampl

    def _set_ipopt_options(self, ampl: AMPL, outlev: int = 0, warm_start: bool = True):
        """Apply standard IPOPT options to *ampl*."""
        ws = "yes" if warm_start else "no"
        ampl.option["solver"] = "ipopt"
        ampl.set_option(
            "ipopt_options",
            f"outlev = {outlev} max_iter = {self.max_iter} mu_init = {self.mu_init} "
            f"expect_infeasible_problem = no tol = {self.tol} "
            f"bound_push = {self.bound_push} bound_frac = {self.bound_frac} "
            f"warm_start_init_point = {ws} halt_on_ampl_error = yes",
        )
        ampl.option["presolve_eps"] = PRESOLVE_EPS_MAIN
        ampl.option["presolve"] = 1

    def _initialise_ampl_vars(self, ampl: AMPL):
        """Set AMPL variable initial values from the current solution."""
        for (x, y, k), val in self.l.items():
            ampl.eval(f"let l[{x},{y},{k}] := {val};")
        for (x, y), val in self.q.items():
            ampl.eval(f"let q[{x},{y}] := {val};")
            if self.data_number == 5:
                ampl.eval(f"let q1[{x},{y}] := {self.q1[x, y]};")
                ampl.eval(f"let q2[{x},{y}] := {self.q2[x, y]};")
        for x, val in self.h.items():
            ampl.eval(f"let h[{x}] := {val};")

    def _add_pipe_length_constraint(self, ampl: AMPL):
        """Add the pipe-length equality constraint (con3) to *ampl*."""
        if self.data_number == 6:
            ampl.eval(
                "subject to con3{(i,j) in arcs diff fixarcs}: "
                "sum{k in pipes} l[i,j,k] = L[i,j];"
            )
        else:
            ampl.eval(
                "subject to con3{(i,j) in arcs}: "
                "sum{k in pipes} l[i,j,k] = L[i,j];"
            )

    def _transfer_duals(self, source_ampl: AMPL, target_ampl: AMPL):
        """Copy dual variable values from *source_ampl* to *target_ampl*."""
        current_duals = {
            name: con.get_values()
            for name, con in target_ampl.get_constraints()
        }
        for name, dual_values in self.all_duals.items():
            if name in current_duals:
                target_ampl.get_constraint(name).set_values(dual_values)

    # ── NLP Sub-problem: QP Approximation (Linearised Model) ─────────────────

    def qp_approximation(
        self,
        iteration: int,
        l: dict,
        q: dict,
        h: dict,
        l_: dict,
        q_: dict,
        h_: dict,
        z_star: float,
        sorted_arcs: list,
        l_points: list,
        q_points: list,
        dual_dict: dict,
        edge,
    ) -> AMPL:
        """
        Build a linearised (QP) model around the current solution point.

        The QP introduces delta variables (d_l, d_q, d_h) and minimises the
        first-order change in pipe cost subject to linearised flow-balance,
        head-loss, and bound constraints.

        Returns
        -------
        AMPL
            Populated AMPL instance ready to be solved with CPLEX.
        """
        ampl_qp = AMPL()
        ampl_qp.reset()
        ampl_qp.read("wdn.mod")

        # ── Network-specific model extensions ─────────────────────────────────
        if self.data_number == 5:
            ampl_qp.eval("""
                param exdiam{arcs};
                param excost{arcs};
                param d_max = max{i in pipes} d[i];
                param d_min = min{i in pipes} d[i];
                param R{arcs};
                param R_min = min{(i,j) in arcs} R[i,j];
                param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
                param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);
                subject to con1{j in nodes diff Source}:
                    sum{i in nodes : (i,j) in arcs} q[i,j]
                    - sum{i in nodes : (j,i) in arcs} q[j,i] = D[j];
                subject to con2{(i,j) in arcs}:
                    h[i] - h[j] = (q1[i,j])^3
                    * (((q1[i,j])^2 + eps[i,j]^2)^0.426
                    / ((q1[i,j])^2 + 0.426*eps[i,j]^2))
                    * omega * L[i,j] / ((R[i,j]^1.852) * (exdiam[i,j])^4.87);
                subject to con2_{(i,j) in arcs}:
                    h[i] - h[j] = (q2[i,j])^3
                    * (((q2[i,j])^2 + eps[i,j]^2)^0.426
                    / ((q2[i,j])^2 + 0.426*eps[i,j]^2))
                    * sum{k in pipes}(omega * l[i,j,k]/(R[i,j]^1.852 * d[k]^4.87));
                subject to con9{(i,j) in arcs}: q[i,j] = q1[i,j] + q2[i,j];
                subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];
                subject to con4{(i,j) in arcs, k in pipes}: l[i,j,k] <= L[i,j];
            """)
        elif self.data_number == 6:
            ampl_qp.eval("""
                set fixarcs within {i in nodes, j in nodes: i != j};
                param d_min;
                param d_max;
                param fix_r{fixarcs};
                param fix_c{fixarcs};
                param fixdiam{fixarcs};
                param R{pipes};
                param R_min = min{k in pipes} R[k];
                param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
                param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);
                subject to con2{(i,j) in arcs diff fixarcs}:
                    h[i] - h[j] = (q[i,j])^3
                    * (((q[i,j])^2 + eps[i,j]^2)^0.426
                    / ((q[i,j])^2 + 0.426*eps[i,j]^2))
                    * sum{k in pipes}(omega * l[i,j,k]/(R[k]^1.852 * d[k]^4.87));
                subject to con2_{(i,j) in fixarcs}:
                    h[i] - h[j] = (q[i,j])^3
                    * (((q[i,j])^2 + eps[i,j]^2)^0.426
                    / ((q[i,j])^2 + 0.426*eps[i,j]^2))
                    * omega * L[i,j] / (fix_r[i,j]^1.852 * fixdiam[i,j]^4.87);
                subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];
                subject to con4{(i,j) in arcs diff fixarcs, k in pipes}: l[i,j,k] <= L[i,j];
            """)
        else:
            ampl_qp.eval("""
                param R{pipes};
                param d_max = max{i in pipes} d[i];
                param d_min = min{i in pipes} d[i];
                param R_min = min{k in pipes} R[k];
                param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
                param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);
            """)

        ampl_qp.read_data(self.data_file)

        # ── Reference parameters ──────────────────────────────────────────────
        ampl_qp.eval("param l_ref {arcs, pipes};")
        ampl_qp.eval("param q_ref {arcs};")
        ampl_qp.eval("param h_ref {nodes};")
        ampl_qp.eval("param x_ref {arcs};")
        ampl_qp.eval("param hessian {arcs};")

        for (u, v, k), val in l.items():
            ampl_qp.param["l_ref"][u, v, k] = val
        for (u, v), val in q.items():
            ampl_qp.param["q_ref"][u, v] = val
        for u, val in h.items():
            ampl_qp.param["h_ref"][u] = val
        for (u, v) in self.arcs:
            ampl_qp.param["x_ref"][u, v] = dual_dict[u, v]

        # Approximate Hessian diagonal for quadratic term
        for (u, v) in self.arcs:
            val = q[u, v]
            if val == 0:
                hess = 0.0
            elif val > 0:
                hess = 1.577904 / val ** 0.148
            else:
                hess = -1.577904 / (-val) ** 0.148
            ampl_qp.param["hessian"][u, v] = hess

        # ── Objective: minimise pipe-cost change (via epigraph variable beta) ─
        ampl_qp.eval("minimize total_cost : beta;")
        ampl_qp.eval(
            "s.t. con_obj: "
            "sum{(i,j) in arcs} sum{k in pipes} d_l[i,j,k]*C[k] - beta <= 0;"
        )

        # ── Linearised flow balance ───────────────────────────────────────────
        ampl_qp.eval("""
            s.t. con1{j in nodes diff Source}:
                sum{i in nodes : (i,j) in arcs} q_ref[i,j]
                - sum{i in nodes : (j,i) in arcs} q_ref[j,i] - D[j]
                + sum{i in nodes : (i,j) in arcs} d_q[i,j]
                - sum{i in nodes : (j,i) in arcs} d_q[j,i] = 0;
        """)

        # ── Linearised head-loss per arc ──────────────────────────────────────
        for (u, v) in self.arcs:
            A_ref = sum(
                10.67 * l[u, v, k] / (self.R[k] ** 1.852 * self.d[k] ** 4.87)
                for k in self.pipes
            )
            g2 = h[u] - h[v] - q[u, v] * abs(q[u, v]) ** 0.852 * A_ref
            q_ref_uv = q[u, v]
            ampl_qp.eval(f"""
                s.t. con2{u}_{v}:
                    {g2}
                    - {A_ref} * 1.852 * abs({q_ref_uv})^0.852 * (d_q[{u},{v}])
                    - {q_ref_uv} * abs({q_ref_uv})^0.852
                      * sum{{k in pipes}} (omega / ((R[k]^1.852) * (d[k])^4.87))
                      * (d_l[{u},{v},k])
                    + d_h[{u}] - d_h[{v}] = 0;
            """)

        # ── Pipe-length equality (linearised) ─────────────────────────────────
        ampl_qp.eval(
            "s.t. con3{(i,j) in arcs}: "
            "sum{k in pipes} l_ref[i,j,k] - L[i,j] "
            "+ sum{k in pipes} d_l[i,j,k] = 0;"
        )

        # ── Head bounds ───────────────────────────────────────────────────────
        ampl_qp.eval("s.t. con4{i in Source}: d_h[i] = 0;")
        ampl_qp.eval(
            "s.t. con5{j in nodes diff Source}: "
            "-h_ref[j] + E[j] + P[j] - d_h[j] <= 0;"
        )

        # ── Delta bounds ──────────────────────────────────────────────────────
        for (u, v) in self.arcs:
            for k in self.pipes:
                ampl_qp.eval(
                    f"s.t. con6{u}_{v}_{k}: "
                    f"{-l[u, v, k]} <= d_l[{u},{v},{k}] <= L[{u},{v}] - {l[u, v, k]};"
                )
            ampl_qp.eval(
                f"s.t. con7{u}_{v}: "
                f"-Q_max - {q[u, v]} <= d_q[{u},{v}] <= Q_max - {q[u, v]};"
            )

        ampl_qp.option["solver"] = "cplex"
        ampl_qp.option["ipopt_options"] = (
            f"outlev = 0 expect_infeasible_problem = yes "
            f"bound_relax_factor = 0 tol = 1e-9 "
            f"bound_push = {self.bound_push} bound_frac = {self.bound_frac} "
            f"warm_start_init_point = yes halt_on_ampl_error = yes max_iter = 500"
        )
        ampl_qp.option["presolve_eps"] = "1.86e-7"

        return ampl_qp

    # ── NLP Sub-problem: Reduced (Perturbed) NLP ──────────────────────────────

    def reduced_nlp_model(
        self,
        iteration: int,
        l_cvx: dict,
        q_cvx: dict,
        h_cvx: dict,
        l: dict,
        q: dict,
        h: dict,
        z_star: float,
        sorted_arcs: list,
        l_points: list,
        q_points: list,
        dual_dict: dict,
    ) -> AMPL:
        """
        Build a perturbed NLP starting from the current solution with a
        random perturbation scaled by the neighbourhood radius parameters
        (eta_l, eta_h, alpha/Delta).

        Returns
        -------
        AMPL
            Populated AMPL instance ready to be solved with IPOPT.
        """
        ampl_nlp = self._build_base_ampl()

        # Trust-region radii
        delta_rho = 0.15
        rho_curr = min(delta_rho * self.k_neigh, 1.0)
        rho_prev = max(rho_curr - delta_rho, 0.0)

        # Perturb pipe lengths
        for (u, v) in self.arcs:
            for k in self.pipes:
                r = random.uniform(rho_prev, rho_curr)
                cand = l[u, v, k] + r * self.eta_l * self.L[u, v]
                ampl_nlp.eval(f"let l[{u},{v},{k}] := {cand};")

        # Perturb flows
        for (u, v) in self.arcs:
            r = random.uniform(rho_prev, rho_curr)
            if self.data_number == 5:
                ampl_nlp.eval(f"let q1[{u},{v}] := {self.q1[u, v] + r * self.Delta};")
                ampl_nlp.eval(f"let q2[{u},{v}] := {self.q2[u, v] + r * self.Delta};")
            else:
                ampl_nlp.eval(f"let q[{u},{v}] := {q[u, v] + r * self.Delta};")

        # Perturb heads
        for u in self.nodes:
            if u not in list(self.source):
                r = random.uniform(rho_prev, rho_curr)
                cand = h[u] + r * self.eta_h * abs(h[u])
                ampl_nlp.eval(f"let h[{u}] := {cand};")

        # Transfer dual values to warm-start the NLP
        current_duals = {
            name: con.get_values()
            for name, con in ampl_nlp.get_constraints()
        }
        for name, dual_values in self.all_duals.items():
            if name in current_duals:
                ampl_nlp.get_constraint(name).set_values(dual_values)

        # self._add_pipe_length_constraint(ampl_nlp)

        ampl_nlp.option["solver"] = "ipopt"
        ampl_nlp.option["ipopt_options"] = (
            f"outlev = 0 "
            f"bound_relax_factor = 0 "
            f"bound_push = {self.bound_push} bound_frac = {self.bound_frac} "
            f"warm_start_init_point = yes "
        )
        ampl_nlp.option["presolve_eps"] = "1.86e-7"

        return ampl_nlp

    # ── Main Solve ────────────────────────────────────────────────────────────

    def solve(self):
        """Solve the initial NLP relaxation."""
        self.ampl.option["solver"] = "ipopt"
        self.ampl.set_option(
            "ipopt_options",
            f"outlev = 1 bound_relax_factor = 0 "
            f"bound_push = {self.bound_push} bound_frac = {self.bound_frac} "
            f"warm_start_init_point = no",
        )
        self.ampl.option["presolve_eps"] = PRESOLVE_EPS_MAIN
        self.ampl.option["presolve"] = 1
        self.ampl.solve()

        self.solve_result = self.ampl.solve_result
        self.total_cost = self.ampl.get_objective("total_cost").value()
        solve_time = self.ampl.get_value("_solve_elapsed_time")
        self.solver_time += solve_time
        self.number_of_nlp += 1

    # ── Heuristic: Local Solution Improvement (Perturbed NLP) ─────────────────

    def local_solution_improvement_heuristic_new(self):
        """
        Adaptive local search using perturbed NLP re-solves (reduced NLP).

        At each call the trust-region radius is adjusted:
        - Improvement found  → shrink radius, reset fail streak, recurse.
        - No improvement     → expand radius, increment fail streak.
        - Termination        → return when flow range is exhausted or
                               fail streak exceeds *total_run*.
        """
        abs_flows = sorted(
            abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4
        )
        m = len(abs_flows)
        median_flow = abs_flows[m // 2]
        improved = False

        # Cache dual variables
        self.all_duals = {
            name: con.getValues()
            for name, con in self.ampl.get_constraints()
        }

        sorted_all_arcs = [a for a in self.arcs if a not in self.fix_arc_set]

        dual_dict = {}
        for con_name, dual_values in self.all_duals.items():
            if con_name == "con7":
                tmp = dual_values.to_dict()
                dual_dict = {node: val for node, val in tmp.items()}
                break

        cvx_nlp = self.reduced_nlp_model(
            self.iteration, self.l_star, self.q_star, self.h_star,
            self.l, self.q, self.h, self.z_star,
            sorted_all_arcs, self.l_points, self.q_points, dual_dict,
        )

        with self.suppress_output():
            cvx_nlp.solve()

        self.solver_time += cvx_nlp.get_value("_solve_elapsed_time")
        self.number_of_nlp += 1

        if cvx_nlp.solve_result != "solved":
            print("solve_result:", cvx_nlp.solve_result)
            return

        self.z_star = cvx_nlp.getObjective("total_cost").value()
        self.l_star = cvx_nlp.getVariable("l").getValues().to_dict()
        self.q_star = cvx_nlp.getVariable("q").getValues().to_dict()
        self.h_star = cvx_nlp.getVariable("h").getValues().to_dict()

        if self.data_number == 5:
            self.q1_star = cvx_nlp.getVariable("q1").getValues().to_dict()
            self.q2_star = cvx_nlp.getVariable("q2").getValues().to_dict()

        self.Z_original[self.number_of_nlp] = self.z_star

        if self.z_star < self.current_cost - 1e-4:
            # ── Improvement accepted ──────────────────────────────────────────
            self._log_nlp_row(
                nlp_num=self.number_of_nlp,
                old_cost=self.current_cost,
                new_cost=self.z_star,
                solve_time=self.ampl.get_value("_solve_elapsed_time"),
                solve_result=self.solve_result,
                improved=True,
            )
            self.current_cost = self.z_star
            self.Z_best[self.number_of_nlp] = self.current_cost
            improved = True
            self.ampl = cvx_nlp
            self.l = cvx_nlp.getVariable("l").getValues().to_dict()
            self.q = cvx_nlp.getVariable("q").getValues().to_dict()
            self.h = cvx_nlp.getVariable("h").getValues().to_dict()
            self.network_graph = self.generate_random_acyclic_from_solution(self.q)

            if self.data_number == 5:
                self.q1 = cvx_nlp.getVariable("q1").getValues().to_dict()
                self.q2 = cvx_nlp.getVariable("q2").getValues().to_dict()
            print("---" * 28)
        else:
            self._log_nlp_row(
                nlp_num=self.number_of_nlp,
                old_cost=self.current_cost,
                new_cost=self.z_star,
                solve_time=cvx_nlp.get_value("_solve_elapsed_time"),
                solve_result=self.solve_result,
                improved=False,
            )

        # Update median flow for next iteration
        abs_flows = sorted(
            abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4
        )
        m = len(abs_flows)
        median_flow = abs_flows[m // 2]

        if improved:
            self.local_iteration += 1
            self.eta_l = self.eta_l_min
            self.eta_h = self.eta_h_min
            self.alpha_q = self.alpha_q_min
            self.Delta = self.alpha_q * median_flow
            self.tr_failure_count = 0
            self.k_neigh = 1
            self.Terminate = False
            self.fail_streak = 0
            self._print_iteration_header(self.local_iteration)
            self.local_solution_improvement_heuristic_new()
        else:
            self.eta_l *= self.eta_l_expend
            self.eta_h *= self.eta_h_expend
            self.alpha_q *= self.alpha_expand
            self.Delta = self.alpha_q * median_flow
            self.fail_streak += 1
            self.k_neigh += 1

            self.Terminate = all(
                abs(self.q[i, j]) + self.Delta > self.Q_max
                for (i, j) in sorted_all_arcs
            )

            should_exit = self.Terminate or self.fail_streak >= self.total_run
            if self.do_local_improvement:
                if should_exit:
                    self.fail_streak = 0
                    self.local_iteration += 1
                    print("---" * 28)
                    return
                else:
                    self.local_solution_improvement_heuristic_new()
            else:
                if should_exit:
                    print("---" * 28)
                    print("Local Search Exits.")
                    return
                else:
                    self.local_solution_improvement_heuristic_new()

    # ── Heuristic: Local Solution Improvement (QP → NLP) ─────────────────────

    def local_solution_improvement_heuristic(self):
        """
        Two-phase local improvement:
          1. Solve the linearised QP to get a descent direction.
          2. Follow the direction with a full NLP re-solve.

        Trust-region parameters are adapted based on success/failure.
        """
        abs_flows = sorted(
            abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4
        )
        m = len(abs_flows)
        median_flow = abs_flows[m // 2]
        improved = False

        sorted_all_arcs = [a for a in self.arcs if a not in self.fix_arc_set]

        cvx_nlp = self.qp_approximation(
            self.iteration, self.l, self.q, self.h,
            self.l, self.q, self.h, self.z_star,
            sorted_all_arcs, self.l_points, self.q_points,
            self.dual_dict, None,
        )

        with self.suppress_output():
            cvx_nlp.solve()

        # Refresh duals and dual dict after QP
        self.all_duals = {
            name: con.get_values()
            for name, con in cvx_nlp.get_constraints()
        }
        self.dual_dict = {
            (u, v): self.all_duals[f"con2{u}_{v}"].to_list()[0]
            for (u, v) in self.arcs
        }

        self.solver_time += cvx_nlp.get_value("_solve_elapsed_time")
        self.number_of_nlp += 1
        self.z_star = cvx_nlp.getObjective("total_cost").value()

        if cvx_nlp.solve_result != "solved":
            self._handle_qp_failure(cvx_nlp, sorted_all_arcs, median_flow)
            return

        # Update starred variables (accumulated delta)
        self.l_star = self.add_dict(self.l_star, cvx_nlp.getVariable("d_l").getValues().to_dict())
        self.q_star = self.add_dict(self.q_star, cvx_nlp.getVariable("d_q").getValues().to_dict())
        self.h_star = self.add_dict(self.h_star, cvx_nlp.getVariable("d_h").getValues().to_dict())

        if self.z_star >= self.current_cost - 1e-4:
            self._log_qp_row(cvx_nlp, improved=False)
            self._handle_qp_no_improvement(cvx_nlp, sorted_all_arcs, median_flow)
            return

        self._log_qp_row(cvx_nlp, improved=True)
        self.Z_red[self.local_iteration] = self.z_star

        # Phase 2: NLP re-solve from QP solution
        ampl = self._build_base_ampl()
        for (x, y, k), val in self.l_star.items():
            ampl.eval(f"let l[{x},{y},{k}] := {val};")
        for (x, y), val in self.q_star.items():
            ampl.eval(f"let q[{x},{y}] := {val};")
            if self.data_number == 5:
                ampl.eval(f"let q1[{x},{y}] := {self.q1_star[x, y]};")
                ampl.eval(f"let q2[{x},{y}] := {self.q2_star[x, y]};")
        for x, val in self.h_star.items():
            ampl.eval(f"let h[{x}] := {val};")

        self._add_pipe_length_constraint(ampl)
        self._set_ipopt_options(ampl, outlev=0, warm_start=True)

        with self.suppress_output():
            ampl.solve()

        self.solve_result = ampl.solve_result
        self.total_cost = ampl.get_objective("total_cost").value()
        self.solver_time += ampl.get_value("_solve_elapsed_time")
        self.number_of_nlp += 1
        self.Z_original[self.local_iteration] = self.total_cost

        if self.solve_result == "solved" and self.total_cost < self.current_cost:
            l = ampl.getVariable("l").getValues().to_dict()
            q = ampl.getVariable("q").getValues().to_dict()
            h = ampl.getVariable("h").getValues().to_dict()

            self._log_nlp_row(
                self.number_of_nlp, self.current_cost, self.total_cost,
                self.ampl.get_value("_solve_elapsed_time"),
                self.solve_result, True, prefix="NLP",
            )
            self.current_cost = self.total_cost
            self.Z_best[self.local_iteration] = self.current_cost
            improved = True
            self.ampl = ampl
            self.network_graph = self.generate_random_acyclic_from_solution(q)
            self.l, self.q, self.h = l, q, h

            if self.data_number == 5:
                self.q1 = ampl.getVariable("q1").getValues().to_dict()
                self.q2 = ampl.getVariable("q2").getValues().to_dict()
            print("---" * 28)
        else:
            self._log_nlp_row(
                self.number_of_nlp, self.current_cost, self.total_cost,
                ampl.get_value("_solve_elapsed_time"),
                self.solve_result, False, prefix="NLP",
            )

        # Recompute median
        abs_flows = sorted(
            abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4
        )
        m = len(abs_flows)
        median_flow = abs_flows[m // 2]

        if improved:
            self.local_iteration += 1
            self.eta_l = self.eta_l_min
            self.eta_h = self.eta_h_min
            self.alpha_q = self.alpha_min
            self.Delta = self.alpha_q * median_flow
            self.tr_failure_count = 0
            self.Terminate = False
            self.fail_streak = 0
            self._print_iteration_header(self.local_iteration)
            self.local_solution_improvement_heuristic()
        else:
            self.eta_l *= self.eta_l_expend
            self.eta_h *= self.eta_h_expend
            self.alpha_q *= self.alpha_expand
            self.Delta = self.alpha_q * median_flow
            self.fail_streak += 1
            self.local_iteration += 1

            self.Terminate = all(
                abs(self.q[i, j]) + self.Delta > self.Q_max
                for (i, j) in sorted_all_arcs
            )
            should_exit = self.Terminate or self.fail_streak >= self.total_run
            if self.do_local_improvement:
                if should_exit:
                    self.fail_streak = 0
                    self.iteration += 1
                    print("---" * 28)
                    return
                else:
                    self.local_solution_improvement_heuristic()
            else:
                if should_exit:
                    print("---" * 28)
                    print("Local Search Exits.")
                    return
                else:
                    self.local_solution_improvement_heuristic()

    # ================================================================
    # INITIALIZE NLP MODEL ONLY ONCE
    # ================================================================
    def initialize_local_search_model(self):
        self.local_ampl = AMPL()
        model_map = {
            5: "newyork_model.mod",
            6: "blacksburg_model.mod",
        }
        self.local_ampl.read(
            model_map.get(self.data_number, "wdnmodel.mod")
        )
        self.local_ampl.read_data(self.data_file)
        self.local_ampl.option["solver"] = "ipopt"
        self.local_ampl.option["ipopt_options"] = (
            "outlev=0 "
            "bound_relax_factor=0 "
            "warm_start_init_point=yes "
            # f"bound_push = {self.bound_push} bound_frac = {self.bound_frac} "
            # "warm_start_bound_push=1e-9 "
            # "warm_start_mult_bound_push=1e-9 "
        )
        self.local_ampl.option["presolve_eps"] = "1.86e-7"
        # ------------------------------------------------------------
        # perturbation parameters
        # ------------------------------------------------------------
        self.local_ampl.eval("""
            param q_start{arcs};
            param h_start{nodes};
            param l_start{arcs,pipes};
            param perturb_q{arcs} default 0;
            param perturb_h{nodes} default 0;
            param perturb_l{arcs,pipes} default 0;
        """)
        # ------------------------------------------------------------
        # optional perturbation equations
        # ------------------------------------------------------------
    
        # self.local_ampl.eval("""
        #
        #     subject to q_warm{(i,j) in arcs}:
        #
        #         q[i,j]
        #         =
        #         q_start[i,j] + perturb_q[i,j];
        #
        #     subject to h_warm{i in nodes diff Source}:
        #
        #         h[i]
        #         =
        #         h_start[i] + perturb_h[i];
        #
        # """)

    def local_solution_improvement_heuristic_fast(self):
        ampl = self.local_ampl
        abs_flows = sorted(
            abs(self.q[i, j])
            for (i, j) in self.arcs
            if abs(self.q[i, j]) > 1e-4
        )
        median_flow = abs_flows[len(abs_flows)//2]
        improved = False
        # ============================================================
        # CACHE DUALS
        # ============================================================
        self.all_duals = {
            name: con.getValues()
            for name, con in self.ampl.get_constraints()
        }
        sorted_all_arcs = [
            a for a in self.arcs
            if a not in self.fix_arc_set
        ]
        # ============================================================
        # TRUST REGION
        # ============================================================
        delta_rho = 0.15
        rho_curr = min(
            delta_rho * self.k_neigh,
            1.0
        )
        rho_prev = max(
            rho_curr - delta_rho,
            0.0
        )
        # ============================================================
        # UPDATE PARAMETERS ONLY
        # ============================================================
        for (u, v) in self.arcs:
            # --------------------------------------------------------
            # q warm start
            # --------------------------------------------------------
            ampl.param["q_start"][u, v] = self.q[u, v]
            r = random.uniform(rho_prev, rho_curr)
            ampl.param["perturb_q"][u, v] = (
                r * self.Delta
            )
            # --------------------------------------------------------
            # variable warm start
            # --------------------------------------------------------
            ampl.var["q"][u, v] = (
                self.q[u, v]
                + ampl.param["perturb_q"][u, v]
            )
        # ============================================================
        # HEADS
        # ============================================================
        for u in self.nodes:
            ampl.param["h_start"][u] = self.h[u]
            if u not in list(self.source):
                r = random.uniform(rho_prev, rho_curr)
                ampl.param["perturb_h"][u] = (
                    r
                    * self.eta_h
                    * abs(self.h[u])
                )
                ampl.var["h"][u] = (
                    self.h[u]
                    + ampl.param["perturb_h"][u]
                )
            else:
                ampl.param["perturb_h"][u] = 0
                ampl.var["h"][u] = self.h[u]
        # ============================================================
        # PIPE LENGTHS
        # ============================================================
        for (u, v) in self.arcs:
            for k in self.pipes:
                ampl.param["l_start"][u, v, k] = (
                    self.l[u, v, k]
                )
                r = random.uniform(rho_prev, rho_curr)
                cand = (
                    self.l[u, v, k]
                    + r * self.eta_l * self.L[u, v]
                )
                ampl.var["l"][u, v, k] = cand
        # ============================================================
        # DUAL WARM START
        # ============================================================
        current_duals = {
            name: con.get_values()
            for name, con in ampl.get_constraints()
        }
        for name, dual_values in self.all_duals.items():
            if name in current_duals:
                ampl.get_constraint(name).set_values(
                    dual_values
                )
        # ============================================================
        # SOLVE
        # ============================================================
        with self.suppress_output():
            ampl.solve()
        self.number_of_nlp += 1
        self.solver_time += ampl.get_value(
            "_solve_elapsed_time"
        )
        if ampl.solve_result != "solved":
            print("solve_result:", ampl.solve_result)
            return
        # ============================================================
        # RETRIEVE SOLUTION
        # ============================================================
        self.z_star = ampl.getObjective(
            "total_cost"
        ).value()
        self.l_star = ampl.getVariable(
            "l"
        ).getValues().to_dict()
        self.q_star = ampl.getVariable(
            "q"
        ).getValues().to_dict()
        self.h_star = ampl.getVariable(
            "h"
        ).getValues().to_dict()
        self.Z_original[
            self.number_of_nlp
        ] = self.z_star
        # ============================================================
        # IMPROVEMENT
        # ============================================================
        if self.z_star < self.current_cost - 1e-4:
            self._log_nlp_row(
                nlp_num=self.number_of_nlp,
                old_cost=self.current_cost,
                new_cost=self.z_star,
                solve_time=ampl.get_value("_solve_elapsed_time"),
                solve_result=ampl.solve_result,
                improved=True,
            )
            self.current_cost = self.z_star
            self.Z_best[
                self.number_of_nlp
            ] = self.current_cost
            improved = True
            self.ampl = ampl
            self.l = self.l_star.copy()
            self.q = self.q_star.copy()
            self.h = self.h_star.copy()
            self.network_graph = (
                self.generate_random_acyclic_from_solution(
                    self.q
                )
            )
            print("---" * 28)
        else:
            self._log_nlp_row(
                nlp_num=self.number_of_nlp,
                old_cost=self.current_cost,
                new_cost=self.z_star,
                solve_time=ampl.get_value("_solve_elapsed_time"),
                solve_result=ampl.solve_result,
                improved=False,
            )
        # ============================================================
        # UPDATE TRUST REGION
        # ============================================================
        abs_flows = sorted(
            abs(self.q[i, j])
            for (i, j) in self.arcs
            if abs(self.q[i, j]) > 1e-4
        )
        median_flow = abs_flows[len(abs_flows)//2]
        if improved:
            self.local_iteration += 1
            self.eta_l = self.eta_l_min
            self.eta_h = self.eta_h_min
            self.alpha_q = self.alpha_q_min
            self.Delta = self.alpha_q * median_flow
            self.fail_streak = 0
            self.k_neigh = 1
            self.Terminate = False
            self._print_iteration_header(
                self.local_iteration
            )
            self.local_solution_improvement_heuristic_fast()
        else:
            self.eta_l *= self.eta_l_expend
            self.eta_h *= self.eta_h_expend
            self.alpha_q *= self.alpha_expand
            self.Delta = self.alpha_q * median_flow
            self.fail_streak += 1
            self.k_neigh += 1
            self.Terminate = all(
                abs(self.q[i, j]) + self.Delta > self.Q_max
                for (i, j) in sorted_all_arcs
            )
            should_exit = (
                self.Terminate
                or self.fail_streak >= self.total_run
            )
            if should_exit:
                print("---" * 28)
                print("Local Search Exits.")
                return
            else:
                self.local_solution_improvement_heuristic_fast()

    # ── Heuristic: Arc Reversal ────────────────────────────────────────────────
    # ─────────────────────────────────────────────────────────────
    # INITIALIZE ONCE
    # ─────────────────────────────────────────────────────────────
    def initialize_arc_reversal_model(self):
        self.arc_ampl = AMPL()
        self.arc_ampl.read(self.model_file)
        self.arc_ampl.read_data(self.data_file)
        self.arc_ampl.option["solver"] = "ipopt"
        self.arc_ampl.option["ipopt_options"] = (
            "outlev=0 "
            "warm_start_init_point=yes "
            "bound_relax_factor=0 "
            # f"bound_push = {self.bound_push} bound_frac = {self.bound_frac} "
            # "warm_start_bound_push=1e-9 "
            # "warm_start_mult_bound_push=1e-9 "
        )
        # ---------------------------------------------------------
        # dynamic direction parameters
        # ---------------------------------------------------------
        self.arc_ampl.eval("""
            param dir_active{arcs} binary default 0;
            param dir_rhs{arcs} default 0;
        """)
        # ---------------------------------------------------------
        # dynamic direction constraints
        # ---------------------------------------------------------
        self.arc_ampl.eval("""
            subject to flow_direction{(i,j) in arcs}:
                dir_rhs[i,j] * q[i,j]
                >=
                -Q_max * (1 - dir_active[i,j]);
        """)
        self.base_ampl = self.arc_ampl

    def iterate_acyclic_flows_fast(self):
        ampl = self.base_ampl
        self.indegree_2_or_more = [
            node
            for node, indeg
            in self.network_graph.in_degree()
            if indeg >= 2
        ]
        print("Iteration:", self.iteration)
        improved = False
        # =========================================================
        # CACHE DUALS
        # =========================================================
        self.all_duals = {
            name: con.getValues()
            for name, con in self.ampl.get_constraints()
        }
        # =========================================================
        # CANDIDATE ARCS
        # =========================================================
        sorted_all_arcs = []
        for node in self.indegree_2_or_more:
            for (u, v) in self.network_graph.in_edges(node):
                arc = (
                    (u, v)
                    if (u, v) in self.arcs
                    else (v, u)
                )
                sorted_all_arcs.append(arc)
        sorted_all_arcs = [
            a for a in sorted_all_arcs
            if a not in self.fix_arc_set
        ]
        sorted_arcs = [
            a for a in sorted_all_arcs
            if a not in self.visited_arc_reverse
        ]
        # =========================================================
        # DUAL SENSITIVITY
        # =========================================================
        dual_dict = {}
        for con_name, dual_values in self.all_duals.items():
            if con_name == "con2":
                tmp = dual_values.to_dict()
                dual_dict = {
                    arc: val
                    for arc, val in tmp.items()
                    if arc in sorted_arcs
                }
                break
        self.sen_score = {}
        for (i, j), dual_val in dual_dict.items():
            resistance_term = sum(
                10.67
                / (
                    float(self.R[k]) ** 1.852
                    * self.d[k] ** 4.87
                )
                for k in self.pipes
            )
            self.sen_score[(i, j)] = (
                -dual_val
                * (
                    2
                    + 1.852
                    * abs(self.q[i, j]) ** 0.852
                    * resistance_term
                    + resistance_term
                    * abs(self.q[i, j]) ** 1.852
                )
            )
        sorted_arcs = [
            arc
            for arc, _
            in sorted(
                self.sen_score.items(),
                key=lambda kv: abs(kv[1]),
                reverse=True
            )
        ]
        print("sorted_arcs:", sorted_arcs)
        self._print_arc_reversal_header()
        # =========================================================
        # MAIN LOOP
        # =========================================================
        while sorted_arcs:
            (i, j) = sorted_arcs.pop(0)
            self.visited_arc_reverse.append((i, j))
            # =====================================================
            # RESET ALL DIRECTION FLAGS
            # =====================================================
            for (u,v) in self.arcs:
                ampl.param["dir_active"][u,v] = 0
                ampl.param["dir_rhs"][u,v] = 1
            # =====================================================
            # ACTIVATE CURRENT REVERSAL
            # =====================================================
            # activate one arc
            ampl.param["dir_active"][i,j] = 1
            
            if self.q[i,j] >= 0:
                ampl.param["dir_rhs"][i,j] = -1
            else:
                ampl.param["dir_rhs"][i,j] = 1
            
            # =====================================================
            # WARM START
            # =====================================================
            self._initialise_ampl_vars(ampl)
            self._transfer_duals(self.ampl, ampl)

            # =====================================================
            # SOLVE
            # =====================================================
            with self.suppress_output():
                ampl.solve()
            self.solve_result = ampl.solve_result
            self.total_cost = ampl.get_objective(
                "total_cost"
            ).value()
            self.solver_time += ampl.get_value(
                "_solve_elapsed_time"
            )
            self.number_of_nlp += 1
            # =====================================================
            # RETRIEVE SOLUTION
            # =====================================================
            if self.solve_result == "solved":
                l_new = ampl.getVariable(
                    "l"
                ).getValues().to_dict()
                q_new = ampl.getVariable(
                    "q"
                ).getValues().to_dict()
                h_new = ampl.getVariable(
                    "h"
                ).getValues().to_dict()
                self.Z_original[
                    self.number_of_nlp
                ] = self.z_star
                # =================================================
                # IMPROVEMENT
                # =================================================
                if self.total_cost < self.current_cost:
                    self._log_arc_row(
                        self.number_of_nlp,
                        (i, j),
                        self.current_cost,
                        self.total_cost,
                        ampl.get_value("_solve_elapsed_time"),
                        self.solve_result,
                        True
                    )
                    self.current_cost = self.total_cost
                    improved = True
                    self.is_improved_in_arc_reversal = True
                    # ---------------------------------------------
                    # update incumbent
                    # ---------------------------------------------
                    self.l = l_new
                    self.q = q_new
                    self.h = h_new
                    self.ampl = ampl
                    self.network_graph = (
                        self.generate_random_acyclic_from_solution(q_new)
                    )
                    self.Z_best[
                        self.number_of_nlp
                    ] = self.current_cost
                    # ---------------------------------------------
                    # remove simultaneously reversed arcs
                    # ---------------------------------------------
                    for (x, y) in list(sorted_arcs):
                        if (
                            q_new[x, y]
                            * self.q.get((x, y), 0)
                            < 0
                        ):
                            if (x, y) not in self.reversed_arcs:
                                self.reversed_arcs.append((x, y))
                                if (x, y) in sorted_arcs:
                                    sorted_arcs.remove((x, y))
                    print("---" * 28)
                    self._start_local_improvement_after_reversal()
                else:
                    self._log_arc_row(
                        self.number_of_nlp,
                        (i, j),
                        self.current_cost,
                        self.total_cost,
                        ampl.get_value("_solve_elapsed_time"),
                        self.solve_result,
                        False
                    )
            else:
                self._log_arc_row(
                    self.number_of_nlp,
                    (i, j),
                    self.current_cost,
                    self.total_cost,
                    ampl.get_value("_solve_elapsed_time"),
                    self.solve_result,
                    False
                )
            # =====================================================
            # RESTART SEARCH
            # =====================================================
            if improved:
                self.iteration = self.local_iteration
                self.iterate_acyclic_flows_fast()
                break

    def iterate_acyclic_flows(self):
        """
        Arc-reversal heuristic: for each candidate arc (ordered by sensitivity
        score), force the flow direction to reverse and re-solve the full NLP.
        If a lower-cost solution is found, switch to it and restart.
        """
        self.indegree_2_or_more = [
            node for node, indeg in self.network_graph.in_degree() if indeg >= 2
        ]
        print("Iteration:", self.iteration)
        improved = False

        # Cache duals
        self.all_duals = {
            name: con.getValues()
            for name, con in self.ampl.get_constraints()
        }

        # Candidate arcs: incoming arcs to high-indegree nodes
        sorted_all_arcs = []
        for node in self.indegree_2_or_more:
            for (u, v) in self.network_graph.in_edges(node):
                arc = (u, v) if (u, v) in self.arcs else (v, u)
                sorted_all_arcs.append(arc)

        sorted_all_arcs = [a for a in sorted_all_arcs if a not in self.fix_arc_set]
        sorted_arcs = [a for a in sorted_all_arcs if a not in self.visited_arc_reverse]

        # Sensitivity scores (con2 duals)
        dual_dict = {}
        for con_name, dual_values in self.all_duals.items():
            if con_name == "con2":
                tmp = dual_values.to_dict()
                dual_dict = {arc: val for arc, val in tmp.items() if arc in sorted_arcs}
                break

        self.sen_score = {}
        for (i, j), dual_val in dual_dict.items():
            if self.data_number == 5:
                self.sen_score[(i, j)] = -dual_val * (
                    2
                    + 1.852 * abs(self.q[i, j]) ** 0.852
                    * sum(10.67 * self.l[i, j, k] / (float(self.R[i, j]) ** 1.852 * self.d[k] ** 4.87) for k in self.pipes)
                    + abs(self.q[i, j]) ** 1.852
                    * sum(10.67 / (float(self.R[i, j]) ** 1.852 * self.d[k] ** 4.87) for k in self.pipes)
                )
            else:
                self.sen_score[(i, j)] = -dual_val * (
                    2
                    + 1.852 * abs(self.q[i, j]) ** 0.852
                    * sum(10.67 * self.l[i, j, k] / (float(self.R[k]) ** 1.852 * self.d[k] ** 4.87) for k in self.pipes)
                    + sum(10.67 / (float(self.R[k]) ** 1.852 * self.d[k] ** 4.87) for k in self.pipes)
                    * abs(self.q[i, j]) ** 1.852
                )

        sorted_arcs = [
            arc for arc, _ in sorted(
                self.sen_score.items(), key=lambda kv: abs(kv[1]), reverse=True
            )
        ]
        print("sorted_arcs:", sorted_arcs)

        self._print_arc_reversal_header()

        while sorted_arcs:
            (i, j) = sorted_arcs.pop(0)
            self.visited_arc_reverse.append((i, j))

            ampl = self._build_base_ampl()
            self._initialise_ampl_vars(ampl)
            self._transfer_duals(self.ampl, ampl)

            # Force flow reversal
            if self.q[i, j] >= 0:
                ampl.eval(f"s.t. flow_direction1{i}_{j}: q[{i},{j}] <= 0;")
            else:
                ampl.eval(f"s.t. flow_direction1{i}_{j}: q[{i},{j}] >= 0;")

            # self._add_pipe_length_constraint(ampl)
            self._set_ipopt_options(ampl, outlev=0, warm_start=True)

            with self.suppress_output():
                ampl.solve()

            self.solve_result = ampl.solve_result
            self.total_cost = ampl.get_objective("total_cost").value()
            self.solver_time += ampl.get_value("_solve_elapsed_time")
            self.number_of_nlp += 1

            if self.solve_result == "solved":
                l = ampl.getVariable("l").getValues().to_dict()
                q = ampl.getVariable("q").getValues().to_dict()
                h = ampl.getVariable("h").getValues().to_dict()
                self.Z_original[self.number_of_nlp] = self.z_star

                if self.total_cost < self.current_cost:
                    self._log_arc_row(self.number_of_nlp, (i, j), self.current_cost,
                                      self.total_cost, self.ampl.get_value("_solve_elapsed_time"),
                                      self.solve_result, True)
                    self.current_cost = self.total_cost
                    improved = True
                    self.is_improved_in_arc_reversal = True
                    self.ampl = ampl
                    self.l, self.q, self.h = l, q, h
                    self.network_graph = self.generate_random_acyclic_from_solution(q)
                    self.Z_best[self.number_of_nlp] = self.current_cost

                    if self.data_number == 5:
                        self.q1 = ampl.getVariable("q1").getValues().to_dict()
                        self.q2 = ampl.getVariable("q2").getValues().to_dict()

                    # Remove any arcs whose direction also changed
                    for (x, y) in list(sorted_arcs):
                        if q[x, y] * self.q.get((x, y), 0) < 0:
                            if (x, y) not in self.reversed_arcs:
                                self.reversed_arcs.append((x, y))
                                if (x, y) in sorted_arcs:
                                    sorted_arcs.remove((x, y))

                    print("---" * 28)
                    self._start_local_improvement_after_reversal()
                else:
                    self._log_arc_row(self.number_of_nlp, (i, j), self.current_cost,
                                      self.total_cost, ampl.get_value("_solve_elapsed_time"),
                                      self.solve_result, False)
                    for (x, y) in list(sorted_arcs):
                        if q[x, y] * self.q.get((x, y), 0) < 0:
                            if (x, y) not in self.reversed_arcs:
                                self.reversed_arcs.append((x, y))
                                if (x, y) in sorted_arcs:
                                    sorted_arcs.remove((x, y))
            else:
                self._log_arc_row(self.number_of_nlp, (i, j), self.current_cost,
                                  self.total_cost, ampl.get_value("_solve_elapsed_time"),
                                  self.solve_result, False)

            if improved:
                self.iteration = self.local_iteration
                self.iterate_acyclic_flows()
                break

    def _start_local_improvement_after_reversal(self):
        """Reset neighbourhood parameters and launch local improvement after an arc reversal."""
        print("---" * 28)
        print("------- Adaptive Variable Neighborhood Search Heuristic --------")
        self.l_points = []
        self.q_points = []
        self.z_star = 0
        self.l_star = self.l
        self.q_star = self.q
        self.h_star = self.h
        self.alpha_q = self.alpha_q_min
        self.eta_l = self.eta_l_min
        self.eta_h = self.eta_h_min

        abs_flows = sorted(abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4)
        m = len(abs_flows)
        self.Delta = self.alpha_q * abs_flows[m // 2]

        self.local_iteration = self.iteration + 1
        self.local_improvement = False
        self.fail_streak = 0
        self._print_iteration_header(self.local_iteration)
        # self.local_solution_improvement_heuristic_new()
        self.initialize_local_search_model()
        self.local_solution_improvement_heuristic_fast()

    # ── Heuristic: Diameter Reduction ────────────────────────────────────────

    def diameter_reduction(self):
        """
        Diameter-reduction heuristic: for each arc (ordered by sensitivity),
        try forcing only smaller diameters and check whether cost decreases.
        """
        arc_max_dia = self._arc_max_diameter()
        print("Iteration:", self.iteration)
        improved = False

        self.all_duals = {
            name: con.getValues()
            for name, con in self.ampl.get_constraints()
        }

        sorted_all_arcs = [
            a for a in self.arcs
            if a not in self.fix_arc_set and a not in self.visited_arc
        ]
        if self.data_number == 6:
            fixarcs = set(self.ampl.getSet("fixarcs"))
            sorted_all_arcs = [a for a in sorted_all_arcs if a not in fixarcs]

        # Only consider arcs that currently use a diameter other than the smallest
        sorted_arcs = [a for a in sorted_all_arcs if arc_max_dia.get(a) != 1]

        # Rank by dual magnitude
        dual_dict = {}
        for con_name, dual_values in self.all_duals.items():
            if con_name == "con2":
                tmp = dual_values.to_dict()
                dual_dict = {arc: val for arc, val in tmp.items() if arc in sorted_arcs}
                break

        sorted_arcs = sorted(dual_dict, key=lambda a: abs(dual_dict[a]), reverse=True)

        print("---" * 28)
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"
              f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")
        print("---" * 28)

        for (i, j) in sorted_arcs:
            self.visited_arc.append((i, j))

            ampl = self._build_base_ampl()
            self._initialise_ampl_vars(ampl)
            self._add_pipe_length_constraint(ampl)

            # Force all diameters ≥ current maximum to zero
            for k in self.pipes:
                if k >= arc_max_dia[i, j]:
                    ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = 0;")
            ampl.eval(
                f"subject to con3_{i}_{j}: "
                f"sum{{k in pipes: k <= {arc_max_dia[i, j] - 1}}} l[{i},{j},k] = L[{i},{j}];"
            )

            ampl.option["solver"] = "ipopt"
            ampl.option["ipopt_options"] = (
                f"outlev = 0 max_iter = {self.max_iter} mu_init = {self.mu_init} "
                f"expect_infeasible_problem = no bound_relax_factor = 0 tol = {self.tol} "
                f"bound_push = {self.bound_push} bound_frac = {self.bound_frac} "
                f"warm_start_init_point = yes halt_on_ampl_error = yes"
            )
            ampl.option["presolve_eps"] = "7.19e-13"

            with self.suppress_output():
                ampl.solve()

            solve_time = ampl.get_value("_solve_elapsed_time")
            self.solver_time += solve_time
            self.number_of_nlp += 1

            total_cost = ampl.getObjective("total_cost").value()
            l1 = ampl.getVariable("l").getValues().to_dict()
            q1 = ampl.getVariable("q").getValues().to_dict()
            h1 = ampl.getVariable("h").getValues().to_dict()

            status = ampl.solve_result
            if status == "solved" and total_cost < self.current_cost:
                print(
                    f"{str((i, j)):<10}"
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{self.format_indian_number(round(total_cost)):<14}"
                    f"{str(round(solve_time, 2)) + 's':<12}"
                    f"{status:<14}{'Yes':<10}"
                    f"{round(time.time() - self.start_time, 2)}s"
                )
                self.current_cost = total_cost
                self.ampl = ampl
                improved = True
                self.is_improved_in_diameter_reduction = True
                self.network_graph = self.generate_random_acyclic_from_solution(q1)
                self.best_acyclic_flow = self.network_graph.copy()
                self.l, self.q, self.h = l1, q1, h1

                if self.data_number == 5:
                    self.q1 = ampl.getVariable("q1").getValues().to_dict()
                    self.q2 = ampl.getVariable("q2").getValues().to_dict()

                print("---" * 28)
                print("--- Adaptive Two-Level Neighborhood Search (ATLNS) ---")
                self.l_star, self.q_star, self.h_star = self.l, self.q, self.h
                self.alpha_q = self.alpha_q_min

                abs_flows = sorted(abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4)
                m = len(abs_flows)
                self.Delta = self.alpha_q * abs_flows[m // 2]
            else:
                print(
                    f"{str((i, j)):<10}"
                    f"{self.format_indian_number(round(self.current_cost)):<14}"
                    f"{self.format_indian_number(round(total_cost)):<14}"
                    f"{str(round(solve_time, 2)) + 's':<12}"
                    f"{status:<14}{'No':<10}"
                    f"{round(time.time() - self.start_time, 2)}s"
                )

            if improved:
                self.iteration += 1
                self.diameter_reduction()
                break

    # ── Logging Helpers ───────────────────────────────────────────────────────

    def _log_nlp_row(
        self,
        nlp_num: int,
        old_cost: float,
        new_cost: float,
        solve_time: float,
        solve_result: str,
        improved: bool,
        prefix: str = "",
    ):
        tag = f"{prefix}{nlp_num}"
        print(
            f"{tag:<8}"
            f"{self.format_indian_number(round(old_cost)):<14}"
            f"{self.format_indian_number(round(new_cost)):<14}"
            f"{str(round(solve_time, 2)) + 's':<12}"
            f"{solve_result:<14}"
            f"{'Yes' if improved else 'No':<10}"
            f"{round(time.time() - self.start_time, 2)}s"
        )

    def _log_qp_row(self, cvx_nlp: AMPL, improved: bool):
        self._log_nlp_row(
            self.number_of_nlp,
            self.current_cost,
            self.z_star,
            cvx_nlp.get_value("_solve_elapsed_time"),
            cvx_nlp.solve_result,
            improved,
            prefix="QP",
        )

    def _log_arc_row(
        self,
        nlp_num: int,
        arc: tuple,
        old_cost: float,
        new_cost: float,
        solve_time: float,
        solve_result: str,
        improved: bool,
    ):
        print(
            f"{nlp_num:<5}"
            f"{str(arc):<10}"
            f"{self.format_indian_number(round(old_cost)):<14}"
            f"{self.format_indian_number(round(new_cost)):<14}"
            f"{str(round(solve_time, 2)) + 's':<12}"
            f"{solve_result:<14}"
            f"{'Yes' if improved else 'No':<10}"
            f"{round(time.time() - self.start_time, 2)}s"
        )

    def _handle_qp_failure(self, cvx_nlp: AMPL, sorted_all_arcs: list, median_flow: float):
        self._log_qp_row(cvx_nlp, improved=False)
        self.eta_l *= self.eta_l_expend
        self.eta_h *= self.eta_h_expend
        self.alpha_q *= self.alpha_expand
        self.Delta = self.alpha_q * median_flow
        self.fail_streak += 1
        self.local_iteration += 1

        if self.fail_streak >= self.total_run:
            self.fail_streak = 0
            self.iteration += 1
            print("---" * 28)
        else:
            self.local_solution_improvement_heuristic()

    def _handle_qp_no_improvement(self, cvx_nlp: AMPL, sorted_all_arcs: list, median_flow: float):
        self._log_qp_row(cvx_nlp, improved=False)
        self.eta_l *= self.eta_l_expend
        self.eta_h *= self.eta_h_expend
        self.alpha_q *= self.alpha_expand
        self.Delta = self.alpha_q * median_flow
        self.fail_streak += 1
        self.local_iteration += 1

        self.Terminate = all(
            abs(self.q[i, j]) + self.Delta > self.Q_max
            for (i, j) in sorted_all_arcs
        )
        if self.fail_streak >= self.total_run:
            self.fail_streak = 0
            self.iteration += 1
            print("---" * 28)
        else:
            self.local_solution_improvement_heuristic()

    @staticmethod
    def _print_iteration_header(iteration: int):
        print(f"Iteration {iteration}")
        print("---" * 28)
        print(
            f"{'NLP':<5}{'C_Best_Sol':<14}{'New_Sol':<14}"
            f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}"
        )
        print("---" * 28)

    @staticmethod
    def _print_arc_reversal_header():
        print("---" * 28)
        print(
            f"{'NLP':<5}{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"
            f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}"
        )
        print("---" * 28)

    # ── AMPL Model Update ─────────────────────────────────────────────────────

    def update_model(self):
        """
        Fix arc flow directions in the AMPL model based on the current
        acyclic directed graph.
        """
        edges_list = [(arc[0], arc[1]) for arc in self.ampl.getSet("arcs")]
        for i, j in self.network_graph.edges:
            if (i, j) in edges_list:
                self.ampl.eval(f"s.t. flow_direction{i}_{j}: q[{i},{j}] >= 0;")
            else:
                self.ampl.eval(f"s.t. flow_direction{i}_{j}: q[{j},{i}] <= 0;")

    # ── Plotly Interactive Plot ───────────────────────────────────────────────

    def build_plot(
        self,
        iteration: int,
        solution_file: str,
        node_pos: dict,
        data_number: int,
        heuristic_approach: str,
        node_head_diff: dict,
        arc_diff: dict,
        edge: tuple,
    ):
        """
        Construct an interactive Plotly figure for a given solution JSON,
        export it as HTML and PDF, and return the figure object.
        """
        # ── Load solution ─────────────────────────────────────────────────────
        with open(solution_file) as f:
            sol = json.load(f)

        def parse_key(k):
            return tuple(map(int, k.replace("(", "").replace(")", "").split(",")))

        q = {parse_key(k): v for k, v in sol["q"].items()}
        L = {parse_key(k): v for k, v in sol.get("L", {}).items()}

        pipe_info = {}
        for k_str, val in sol.get("l", {}).items():
            if val < 1e-6:
                continue
            i, j, d = map(int, k_str.replace("(", "").replace(")", "").split(","))
            pipe_info.setdefault((i, j), []).append((d, val))

        h = {int(k): v for k, v in sol["h"].items()}
        hmin = {int(k): v for k, v in sol["hmin"].items()}
        D = {int(k): v for k, v in sol["D"].items()}
        source = sol.get("source", [])
        diameters = {int(k): v for k, v in sol["Diameters"].items()}

        # ── Figure scaling ────────────────────────────────────────────────────
        node_marker_size = 22
        BASE = 1800
        xs_pos = [p[0] for p in node_pos.values()]
        ys_pos = [p[1] for p in node_pos.values()]
        data_w = max(xs_pos) - min(xs_pos)
        data_h = max(ys_pos) - min(ys_pos)
        node_radius = (node_marker_size / 2) * (data_w / BASE)

        # ── Diameter styling ──────────────────────────────────────────────────
        DATASET_DIAMETERS = list(diameters.keys())
        thickness_levels = list(range(3, 25))
        color_palette = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#008080",
            "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
            "#bcbd22", "#17becf", "#393b79", "#637939",
            "#8c6d31", "#843c39", "#7b4173", "#3182bd",
            "#31a354", "#756bb1", "#636363", "#e6550d",
            "#969696", "#6baed6",
        ]
        DIAM_INDEX = {d: i for i, d in enumerate(DATASET_DIAMETERS)}

        all_diams = {d for pipes in pipe_info.values() for d, _ in pipes}
        diam_to_style = {
            d: {"width": thickness_levels[DIAM_INDEX[d]], "color": color_palette[DIAM_INDEX[d]]}
            for d in all_diams
        }

        # ── Build edge traces ─────────────────────────────────────────────────
        edge_groups: dict = {}
        arrows = []
        click_x, click_y, click_text = [], [], []
        flow_tx, flow_ty, flow_txt = [], [], []
        pipe_tx, pipe_ty, pipe_txt = [], [], []
        traces = []

        for (i, j), flow in q.items():
            if i not in node_pos or j not in node_pos:
                continue
            x0, y0 = node_pos[i]
            x1, y1 = node_pos[j]
            xs_e, ys_e, xe_e, ye_e = (x0, y0, x1, y1) if flow >= 0 else (x1, y1, x0, y0)
            dx, dy = xe_e - xs_e, ye_e - ys_e
            dist = math.hypot(dx, dy)
            if dist < 1e-9:
                continue

            sx = xs_e + node_radius * dx / dist
            sy = ys_e + node_radius * dy / dist
            ex = xe_e - node_radius * dx / dist
            ey = ye_e - node_radius * dy / dist
            mx, my = (sx + ex) / 2, (sy + ey) / 2

            pipes = pipe_info.get((i, j), [])
            pipe_label = ", ".join(f"D{d}:{lv:.1f}" for d, lv in pipes)
            hover = (
                f"<b>Arc {i} → {j}</b><br>"
                f"Flow: {flow:.6f}<br>"
                f"Length: {L.get((i, j), 0):.2f}<br>"
                f"<b>Pipes:</b><br>{pipe_label}"
            )
            click_x.append(mx); click_y.append(my); click_text.append(hover)
            flow_tx.append(mx); flow_ty.append(my); flow_txt.append(f"{flow:.5f}")
            pipe_tx.append(mx); pipe_ty.append(my - 0.02 * data_h)
            pipe_txt.append(pipe_label)

            # Multi-diameter segmentation
            total_len = sum(lv for _, lv in pipes)
            if total_len > 0:
                cur_len = 0.0
                for d, lv in sorted(pipes, key=lambda x: x[0], reverse=True):
                    frac_s = cur_len / total_len
                    frac_e = (cur_len + lv) / total_len
                    seg_sx = sx + frac_s * dx
                    seg_sy = sy + frac_s * dy
                    seg_ex = sx + frac_e * dx
                    seg_ey = sy + frac_e * dy
                    style = diam_to_style[d]
                    edge_groups.setdefault(d, {"x": [], "y": []})
                    edge_groups[d]["x"] += [seg_sx, seg_ex, None]
                    edge_groups[d]["y"] += [seg_sy, seg_ey, None]
                    cur_len += lv

            # Arrow
            arrow_color = "red" if (i, j) == edge else "black"
            arrows.append(dict(
                ax=sx, ay=sy, x=ex, y=ey,
                xref="x", yref="y", axref="x", ayref="y",
                showarrow=True, arrowhead=3, arrowsize=2,
                arrowwidth=1, arrowcolor=arrow_color,
            ))

        # ── Diameter traces ───────────────────────────────────────────────────
        for d in DATASET_DIAMETERS:
            if d not in edge_groups:
                continue
            data = edge_groups[d]
            style = diam_to_style[d]
            traces.append(go.Scatter(
                x=data["x"], y=data["y"], mode="lines",
                line=dict(width=style["width"], color=style["color"]),
                name=f"Diameter D{d}: {diameters[d]} m",
                hoverinfo="skip",
            ))

        traces.append(go.Scatter(
            x=click_x, y=click_y, mode="markers",
            marker=dict(size=0, opacity=0),
            hoverinfo="text", hovertext=click_text, showlegend=False,
        ))

        # Highlight the modified edge
        i_e, j_e = edge
        if i_e in node_pos and j_e in node_pos:
            x0, y0 = node_pos[i_e]
            x1, y1 = node_pos[j_e]
            traces.append(go.Scatter(
                x=[x0, x1], y=[y0, y1], mode="lines",
                line=dict(color="red", width=5),
                hoverinfo="skip", showlegend=False,
            ))

        # ── Node traces ───────────────────────────────────────────────────────
        nx_src, ny_src, label_src, text_src = [], [], [], []
        nx_dem, ny_dem, label_dem, text_dem = [], [], [], []
        hx_src, hy_src, htxt_src = [], [], []
        hx_dem, hy_dem, htxt_dem = [], [], []
        hx_min, hy_min, hmin_txt = [], [], []
        hx_demval, hy_demval, dem_txt = [], [], []

        for n, (x, y) in node_pos.items():
            hover_text = (
                f"<b>Node {n}</b><br>"
                f"Head: {h.get(n, 0):.2f}<br>"
                f"Min Head: {hmin.get(n, 0):.2f}<br>"
                f"Demand: {D.get(n, 0):.5f}"
            )
            if n in source:
                nx_src.append(x); ny_src.append(y)
                label_src.append(str(n)); text_src.append(hover_text)
                hx_src.append(x); hy_src.append(y + 0.03 * data_h)
                htxt_src.append(f"{h.get(n, 0):.2f}")
            else:
                nx_dem.append(x); ny_dem.append(y)
                label_dem.append(str(n)); text_dem.append(hover_text)
                hx_dem.append(x); hy_dem.append(y + 0.03 * data_h)
                htxt_dem.append(f"{h.get(n, 0):.2f}")

            hx_min.append(x - 0.03 * data_w)
            hy_min.append(y - 0.02 * data_h)
            hmin_txt.append(f"{hmin.get(n, 0):.2f}")
            hx_demval.append(x)
            hy_demval.append(y - 0.035 * data_h)
            dem_txt.append(f"{D.get(n, 0):.4f}")

        traces += [
            go.Scatter(
                x=nx_src, y=ny_src, mode="markers+text",
                text=label_src, textposition="middle center",
                hoverinfo="text", hovertext=text_src,
                marker=dict(size=node_marker_size + 4, color="royalblue",
                            symbol="circle", line=dict(width=2, color="black")),
                name="Source Nodes",
            ),
            go.Scatter(
                x=nx_dem, y=ny_dem, mode="markers+text",
                text=label_dem, textposition="middle center",
                hoverinfo="text", hovertext=text_dem,
                marker=dict(size=node_marker_size, color="skyblue",
                            line=dict(width=2, color="black")),
                name="Demand Nodes",
            ),
            go.Scatter(
                x=hx_min, y=hy_min, mode="text", text=hmin_txt,
                textfont=dict(size=12, color="darkred"),
                name="Min Head Requirement", hoverinfo="skip",
                visible="legendonly",
            ),
            go.Scatter(
                x=hx_demval, y=hy_demval, mode="text", text=dem_txt,
                textfont=dict(size=12, color="darkgreen"),
                name="Node demand", hoverinfo="skip",
                visible="legendonly",
            ),
            go.Scatter(x=flow_tx, y=flow_ty, mode="text",
                       text=flow_txt, name="Flow Values"),
            go.Scatter(x=pipe_tx, y=pipe_ty, mode="text",
                       text=pipe_txt, name="Pipe Info"),
            go.Scatter(x=hx_src, y=hy_src, mode="text",
                       text=htxt_src, name="Source Node Head"),
            go.Scatter(x=hx_dem, y=hy_dem, mode="text",
                       text=htxt_dem, name="Demand Node Head"),
        ]

        # ── Arc-difference overlay ────────────────────────────────────────────
        mx_diff, my_diff, htxt_arc = [], [], []
        for (i, j), info in arc_diff.items():
            if i not in node_pos or j not in node_pos:
                continue
            x0, y0 = node_pos[i]
            x1, y1 = node_pos[j]
            xm, ym = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
            txt = f"<b>Arc ({i} → {j})</b><br>"
            if "flow" in info:
                f_info = info["flow"]
                txt += (f"Flow change<br>q₁ = {f_info['q_sol1']}<br>"
                        f"q₂ = {f_info['q_sol2']}<br>Δq = {f_info['delta']}<br>")
            if "pipes" in info:
                txt += (f"Pipes changed<br>sol1: {info['pipes']['sol1']}<br>"
                        f"sol2: {info['pipes']['sol2']}")
            mx_diff.append(xm); my_diff.append(ym); htxt_arc.append(txt)

        ax_diff, ay_diff = [], []
        for (i, j) in arc_diff:
            if i not in node_pos or j not in node_pos:
                continue
            x0, y0 = node_pos[i]; x1, y1 = node_pos[j]
            ax_diff += [x0, x1, None]; ay_diff += [y0, y1, None]

        traces += [
            go.Scatter(
                x=mx_diff, y=my_diff, mode="markers",
                marker=dict(size=12, color="rgba(0,0,0,0)"),
                hoverinfo="text", hovertext=htxt_arc,
                showlegend=False, legendgroup="diff", visible="legendonly",
            ),
            go.Scatter(
                x=ax_diff, y=ay_diff, mode="lines",
                line=dict(width=3, color="red", dash="dash"),
                name="Affected arcs", legendgroup="diff", visible="legendonly",
            ),
        ]

        # ── Node-difference overlay ───────────────────────────────────────────
        node_set = set(node_head_diff.keys())
        dx_diff, dy_diff, dlabel, dtext = [], [], [], []
        nocdx, nocdy, nocdlabel, nocdtext = [], [], [], []
        nocsx, nocsy, nocslabel, nocstext = [], [], [], []

        for n, (x, y) in node_pos.items():
            if n in node_set:
                info = node_head_diff[n]
                dx_diff.append(x); dy_diff.append(y); dlabel.append(str(n))
                dtext.append(
                    f"<b>Node {n}</b><br>Head changed<br>"
                    f"h₁ = {info['h_sol1']}<br>h₂ = {info['h_sol2']}<br>"
                    f"Δh = {info['delta']}"
                )
            elif n in source:
                nocsx.append(x); nocsy.append(y); nocslabel.append(str(n))
                nocstext.append(
                    f"<b>Source Node {n}</b><br>Head fixed<br>"
                    f"h = {h.get(n, 0):.2f}<br>"
                )
            else:
                nocdx.append(x); nocdy.append(y); nocdlabel.append(str(n))
                nocdtext.append(
                    f"<b>Node {n}</b><br>Head Unchanged<br>h₂ = {h.get(n, 0):.2f}<br>"
                )

        traces += [
            go.Scatter(
                x=dx_diff, y=dy_diff, mode="markers + text",
                text=dlabel, textposition="middle center",
                marker=dict(size=node_marker_size, color="skyblue",
                            symbol="circle", line=dict(width=2, color="red")),
                hoverinfo="text", hovertext=dtext,
                name="Affected nodes", legendgroup="diff",
                legendgrouptitle=dict(text="Local solution difference"),
                visible="legendonly",
            ),
            go.Scatter(
                x=nocdx, y=nocdy, mode="markers + text",
                text=nocdlabel, textposition="middle center",
                marker=dict(size=node_marker_size, color="skyblue",
                            symbol="circle", line=dict(width=2, color="black")),
                hoverinfo="text", hovertext=nocdtext,
                name="Unaffected nodes", legendgroup="diff",
                visible="legendonly",
            ),
            go.Scatter(
                x=nocsx, y=nocsy, mode="markers + text",
                text=nocslabel, textposition="middle center",
                marker=dict(size=node_marker_size + 4, color="royalblue",
                            symbol="circle", line=dict(width=2, color="black")),
                hoverinfo="text", hovertext=nocstext,
                legendgroup="diff", showlegend=False,
            ),
        ]

        # ── Layout ────────────────────────────────────────────────────────────
        title_text = (
            f"Network: {DATA_LIST[data_number]} | "
            f"Heuristic Approach: {heuristic_approach} | "
            f"Iteration: {iteration} | Arc: {edge} | "
            f"Objective: {sol['objective']:.2f} | "
            f"Time: {sol['solve_time']:.2f} s"
        )
        fig = go.Figure(traces)
        fig.update_layout(
            title=dict(
                text=title_text,
                font=dict(size=24, color="black", family="Arial"),
            ),
            annotations=arrows,
            dragmode="pan",
            hovermode="closest",
            autosize=False,
            width=1900,
            height=1100,
            xaxis=dict(
                visible=True, showgrid=True, zeroline=True, fixedrange=False,
                showline=True, linecolor="black", mirror=True, ticks="outside",
            ),
            yaxis=dict(
                visible=True, showgrid=True, zeroline=True, scaleanchor="x",
                fixedrange=False, showline=True, linecolor="black",
                mirror=True, ticks="outside",
            ),
            plot_bgcolor="white",
            paper_bgcolor="white",
            legend=dict(
                title=dict(
                    text="Network Elements",
                    font=dict(size=20, family="Arial", color="black"),
                ),
                font=dict(size=16, family="Arial", color="black"),
                orientation="v", x=1.02, y=1,
                xanchor="left", yanchor="top",
                bgcolor="rgba(255,255,255,0.85)",
                bordercolor="rgba(0,0,0,0.3)", borderwidth=2,
                itemsizing="constant", itemwidth=40,
            ),
        )

        out_dir = f"../figure/json_file/d{data_number + 1}"
        fig.write_html(
            f"{out_dir}/wdn_interactive_solution{iteration}.html",
            auto_open=False,
            config={"scrollZoom": True},
        )
        pio.write_image(
            fig,
            f"{out_dir}/wdn_interactive_solution{self.iteration}.pdf",
            format="pdf",
        )
        return fig

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

        with self.suppress_output():
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

        # if self.data_number == 6:
        #     ampl.eval("subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
        # else:
        #     ampl.eval("subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];")

        ampl.option['solver'] = 'ipopt'
        ampl.option["ipopt_options"] = (
            "outlev=0 "
            "bound_relax_factor=0 "
            # f"bound_push={self.bound_push} "
            # f"bound_frac={self.bound_frac} "
            "warm_start_init_point=yes "
        )

        with self.suppress_output():
            ampl.solve()

        solve_result = ampl.get_value("solve_result")
        q_sol = ampl.get_variable('q').get_values().to_dict()
        h_sol = ampl.get_variable('h').get_values().to_dict()
        l_sol = ampl.get_variable('l').get_values().to_dict()
        cost = ampl.getObjective("total_cost").value()
        solve_time = ampl.get_value('_solve_elapsed_time')

        return ampl, solve_result, q_sol, h_sol, l_sol, cost, solve_time

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
        # heuristic_start = time.time() - self.start_time
    
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
                "outlev=0 bound_relax_factor=0 warm_start_init_point=yes warm_start_bound_push=1e-9 warm_start_mult_bound_push=1e-9"
            )
            ampl.option["presolve_eps"] = "6.82e-14"
            ampl.option["presolve"]     = 1

            with self.suppress_output():
                ampl.solve()
            self.all_duals = {}
            for con_name, con in ampl.get_constraints():
                self.all_duals[con_name] = con.getValues()

            solve_time = ampl.get_value("_solve_elapsed_time")

            q = ampl.getVariable("q").getValues().to_dict()
            h = ampl.getVariable("h").getValues().to_dict()
            y = ampl.getVariable("y").getValues().to_dict()
            z = ampl.getVariable("z").getValues().to_dict()
            print("\n-------------------------------- Solving Recover Model --------------------------")

            rec_ampl, rec_result, l_trial, recovered_cost, t_rec = self.solve_recover_model1(y)
            print(rec_result)
            print(recovered_cost)

            print("\n-------------------- Solving Original Model with Initialize ---------------------")
            orig_ampl, orig_result, q_new, h_new, l_new, final_cost, t_orig = self.solve_original_with_init(l_trial, q, h)

            print(orig_result, final_cost)

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
            ("Total arcs",        str(len(list(self.arcs)))),
        ], title="Configuration")

        # print_segment_table(self.seg_index, title="Initial segment assignments")

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
                # cumul_time          = time.time() - heuristic_start
                cumul_time          = time.time() - self.start_time

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

        # if solve_log:
        #     print_solve_table(solve_log, title="Complete solve history")

        # print_segment_table(best_solution["seg"], title="Best solution — segment assignments")

        return (
            best_solution["q"],
            best_solution["h"],
            best_solution["y"],
            best_solution["z"],
        )

    # ============================================================
    # 🔷 ORIGINAL → REDUCED
    # ============================================================
    def map_original_to_reduced(self, l, q, h):
        # ============================================================
        # Compute alpha values
        # ============================================================
        alpha = {
            k: self.omega / (self.R[k]**1.852 * self.d[k]**4.87)
            for k in self.pipes
        }
        # ============================================================
        # Segment information
        # ============================================================
        NP = len(list(self.pipes))
        # Segments between consecutive pipe types
        segs = list(range(1, NP))
        slope = {}
        intercept = {}
        # ============================================================
        # Compute piecewise linear coefficients
        # ============================================================
        for s in segs:
            alpha1 = alpha[s]
            alpha2 = alpha[s + 1]
            C1 = self.C[s]
            C2 = self.C[s + 1]
            # ---- slope ----
            slope[s] = (C1 - C2) / (alpha1 - alpha2)
            # ---- intercept ----
            intercept[s] = (
                alpha1 * C2 - alpha2 * C1
            ) / (alpha1 - alpha2)

        # print("slope:", slope)
        # print("intercept:", intercept)
        # ============================================================
        # Compute reduced variables
        # ============================================================
        y = {}
        z = {}
        for (i, j) in self.arcs:
            # ---- compute y ----
            y[(i, j)] = sum(
                alpha[k] * l[(i, j, k)] / self.L[i, j]
                for k in self.pipes
            )
            # ---- compute z (epigraph value) ----
            s = self.find_seg(y[(i,j)])[0]
            z[(i, j)] = (
                (
                    slope[s] * y[(i, j)]
                    + intercept[s]
                ) * self.L[i, j]
            )
        return y, z

    
    def find_seg(self, yv):
        alpha_vals = sorted(set(self.alpha.values()), reverse=True)
    
        NP = len(alpha_vals)

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
    
    def true_cost(self,y_sol, seg_sol):
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

    def solve_exact_reduced_model(self):
        # ============================================================
        # Initialize AMPL model
        # ============================================================
        ampl = AMPL()
        ampl.read("exact_reduced_wdn.mod")
        ampl.read_data(self.data_file)
        # ============================================================
        # Optional warm-start initialization
        # ============================================================
        for (i, j), val in self.q.items():
            ampl.var["q"][i, j] = val
        for i, val in self.h.items():
            ampl.var["h"][i] = val
        for (i, j) in self.arcs:
            ampl.var["y"][i, j] = self.y[i,j] 
            ampl.var["z"][i, j] = self.z[i,j] 

        # ============================================================
        # DUAL WARM START
        # ============================================================
        # current_duals = {
        #     name: con.get_values()
        #     for name, con in ampl.get_constraints()
        # }
        # for name, dual_values in self.all_duals.items():
        #     if name in current_duals:
        #         ampl.get_constraint(name).set_values(
        #             dual_values
        #         )
        # ============================================================
        # Solver selection
        # ============================================================
        ampl.option["solver"] = "ipopt"
        # Alternative solvers
        # ampl.option["solver"] = "baron"
        # ============================================================
        # IPOPT options
        # ============================================================
        ampl.option["ipopt_options"] = (
            "outlev=1 "
            "bound_relax_factor=0 "
            # "tol=1e-6 "
            # f"bound_push={self.bound_push} "
            # f"bound_frac={self.bound_frac} "
            "warm_start_init_point=yes "
            # "halt_on_ampl_error=yes "
            "warm_start_bound_push=1e-9 "
            "warm_start_mult_bound_push=1e-9 "
        )
        # ============================================================
        # BARON options
        # ============================================================
        ampl.option["baron_options"] = (
            "maxtime=3600 "
            "outlev=2 "
            "barstats "
            "version "
            "objbound"
        )
        # ============================================================
        # Gurobi options
        # ============================================================
        ampl.option["gurobi_options"] = (
            "outlev 1 "
            "presolve 1 "
            "timelimit 3600"
        )
        # ============================================================
        # SCIP options
        # ============================================================
        ampl.option["scip_options"] = (
            "outlev 1 "
            "timelimit 3600 "
            "lim:gap=1e-9 "
            "chk:feastol=1e-5 "
            "chk:feastolrel=0"
        )
        # Alternative Gurobi settings
        # ampl.option["gurobi_options"] = (
        #     "outlev 1 "
        #     "presolve 1 "
        #     "timelimit 3600 "
        #     "warmstart=0 "
        #     "barconvtol=1e-9 "
        #     "feastol=1e-5 "
        #     "chk:epsrel=0 "
        #     "mipgap=1e-9 "
        #     "NumericFocus=1"
        # )
        # ampl.eval("option presolve_eps 6.98e-11;")
        # ============================================================
        # Solve optimization model
        # ============================================================
        ampl.solve()
        # ============================================================
        # Retrieve solver status
        # ============================================================
        solve_result = ampl.get_value("solve_result")
        # ============================================================
        # Retrieve solution variables
        # ============================================================
        q_sol = (ampl.get_variable("q").get_values().to_dict())
        y_sol = (ampl.get_variable("y").get_values().to_dict())
        h_sol = (ampl.get_variable("h").get_values().to_dict())
        z_sol = (ampl.get_variable("z").get_values().to_dict())
        
        # ============================================================
        # Retrieve objective value and solve time
        # ============================================================
        cost = ampl.getObjective("total_cost").value()
        solve_time = ampl.get_value("_solve_elapsed_time")
        # ============================================================
        # Retrieve parameter values
        # ============================================================
        # self.alpha = (
        #     ampl.get_parameter("alpha")
        #     .get_values()
        #     .to_dict()
        # )
        # ============================================================
        # Optional recovered flow computation
        # ============================================================
        # lambda_vals = (
        #     ampl.getVariable('lambda')
        #     .getValues()
        #     .toDict()
        # )
        # q_p = (
        #     ampl.getParameter('q_p')
        #     .getValues()
        #     .toDict()
        # )
        # B = (
        #     ampl.getParameter('B')
        #     .getValues()
        #     .toDict()
        # )
        # cycles = ampl.getSet('cycles').to_list()
        # q_sol = {}
        # for (i, j) in self.arcs:
        #     val = q_p.get((i, j), 0.0)
        #     for c in cycles:
        #         val += (
        #             B.get((i, j, c), 0.0)
        #             * lambda_vals.get(c, 0.0)
        #         )
        #     q_sol[(i, j)] = val
        # print("Recovered flows q:", q_sol)
        # ============================================================
        # Optional debugging
        # ============================================================
        # ampl.eval("display z;")
        # ampl.eval("display x;")
        # obj = 0
        # for (i, j) in self.arcs:
        #     obj += z_sol[i, j]**0.5
        # print("Objective:", obj)
        # for (i, j) in self.arcs:
        #     for k in self.pipes:
        #         print(
        #             f"l[{i},{j},{k}]:",
        #             (
        #                 (y_sol[i, j]
        #                 - alpha[k+1] * self.L[i, j])
        #                 / alpha[k]
        #             ) - alpha[k+1]
        #         )
        # ============================================================
        # Return results
        # ============================================================
        return (
            ampl,
            solve_result,
            z_sol,
            q_sol,
            h_sol,
            y_sol,
            cost,
            solve_time,
        )


    # ── Main Entry Point ──────────────────────────────────────────────────────

    def run(self):
        """
        Execute the full heuristic pipeline:

        1. Load model and solve the initial smooth NLP relaxation.
        2. Run adaptive local improvement (perturbed NLP re-solves).
        3. Apply arc-reversal perturbations.
        4. Print final statistics.
        """
        self.start_time = time.time()
        self.Z_red: dict = {}
        self.Z_original: dict = {}
        self.Z_best: dict = {}

        # ── Step 1: Initial NLP ───────────────────────────────────────────────
        print("NLP Model: Smooth Approximate WDN Design Model 2, Epsilon via Relative Error")
        print("NLP Solver: IPOPT")
        print("=" * 80)
        print("Initial Solution of Approximate WDN Design Model")
        print("=" * 80)

        self.load_model()
        fix_arc_set_result = self.fix_leaf_arc_flow(self.ampl)
        print("fix_arc_set:", fix_arc_set_result)

        self.super_source_out_arc = self.fix_arc_set()
        print("super_source_out_arc:", self.super_source_out_arc, "\n")


        self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set_result)        
        self.sorted_arcs = [
            arc for arc in self.arcs
            if arc not in self.fix_arc_set
        ]

        # if self.data_number == 6:
        #     self.ampl.eval(
        #         "subject to con3{(i,j) in arcs diff fixarcs}: "
        #         "sum{k in pipes} l[i,j,k] - L[i,j] = 0;"
        #     )
        # else:
        #     self.ampl.eval(
        #         "subject to con3{(i,j) in arcs}: "
        #         "sum{k in pipes} l[i,j,k] - L[i,j] = 0;"
        #     )

        self.solve()
        print(f"Objective:    {self.total_cost}")
        print(f"Solve result: {self.solve_result}")
        print(f"Solve time:   {self.ampl.get_value('_solve_elapsed_time'):.2f} s\n")

        if self.solve_result != "solved":
            print("IPOPT did not solve the initial problem optimally. Exiting.")
            sys.exit()

        self.current_cost = self.total_cost
        self.l = self.ampl.getVariable("l").getValues().to_dict()
        self.q = self.ampl.getVariable("q").getValues().to_dict()
        self.h = self.ampl.getVariable("h").getValues().to_dict()

        if self.data_number == 5:
            self.q1 = self.ampl.getVariable("q1").getValues().to_dict()
            self.q2 = self.ampl.getVariable("q2").getValues().to_dict()

        # Cache initial duals and con2 dual dict
        self.all_duals = {
            name: con.getValues()
            for name, con in self.ampl.get_constraints()
        }
        self.dual_dict = {}
        for con_name, dual_values in self.all_duals.items():
            if con_name == "con2":
                self.dual_dict = {node: val for node, val in dual_values.to_dict().items()}
                break

        # ── Step 2: Build acyclic graph and initialise neighbourhood ──────────
        print("=" * 80)
        print("Improve the Initial Solution")
        print("=" * 80)

        self.network_graph = self.generate_random_acyclic_from_solution(self.q)
        self.indegree_2_or_more = [
            node for node, indeg in self.network_graph.in_degree() if indeg >= 2
        ]
        self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set_result)
        self.best_acyclic_flow = self.network_graph.copy()
        self.visited_nodes: list = []
        self.sorted_nodes: list = []

        self.iteration = 0
        self.Z_original[self.iteration] = self.current_cost
        self.Z_red[self.iteration] = None
        self.elapsed_time = time.time() - self.start_time

        # Neighbourhood parameters
        self.iteration = 1
        self.l_points: list = []
        self.q_points: list = []
        self.z_star = 0.0
        self.l_star = self.l
        self.q_star = self.q
        self.h_star = self.h

        abs_flows = sorted(abs(self.q[i, j]) for (i, j) in self.arcs if abs(self.q[i, j]) > 1e-4)
        m = len(abs_flows)
        self.Delta = self.alpha_q * abs_flows[m // 2]

        self.local_iteration = 1
        self.do_local_improvement = False
        self.local_improvement = False
        self.do_arc_reversal = True
        self.total_run = 3

        self._print_iteration_header(self.local_iteration)
        # self.local_solution_improvement_heuristic_new()
        self.initialize_local_search_model()
        self.local_solution_improvement_heuristic_fast()
        # ── Step 3: Arc-Reversal Heuristic ────────────────────────────────────
        print("---" * 28)
        print("Reverse Arc Direction Approach")
        print("---" * 28)

        self.iteration = self.local_iteration + 1
        self.visited_arc_reverse: list = []
        self.do_arc_reversal = True
        self.do_local_improvement = True
        self.is_improved_in_arc_reversal = False
        self.do_diameter_reduction = False
        self.reversed_arcs: list = []
        # self.iterate_acyclic_flows()
        self.initialize_arc_reversal_model()
        self.iterate_acyclic_flows_fast()

        print("\n-------------- Piecewise Linear Water Distribution Network Design Model----------------")
        print("\n-------------------------------- Solving Exact Reduced Model --------------------------")

        self.y, self.z = self.map_original_to_reduced(self.l, self.q, self.h)
        print(sum(self.z[i,j] for (i,j) in self.arcs))

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
        
        self.segs = list(self.ampl.getSet('segs'))
        self.alpha = self.ampl.get_parameter("alpha").get_values().to_dict()
        self.slope = self.ampl.getParameter('slope').to_dict()
        self.intercept = self.ampl.getParameter('intercept').to_dict()
        # print("slope:", self.slope)
        # print("intercept:", self.intercept)
        self.alpha_min = min(self.alpha[k] for k in self.segs)
        self.alpha_max = max(self.alpha[k] for k in self.segs)

        self.y_lb = {}
        self.y_ub = {}
        for (i, j) in self.arcs:
            self.y_lb[(i, j)] = self.omega * self.L[(i, j)] / (self.R_max**1.852 * self.d_max**4.87)
            self.y_ub[(i, j)] = self.omega * self.L[(i, j)] / (self.R_min**1.852 * self.d_min**4.87)
 
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
            self.seg_index[(u, v)] = self.find_seg(y_val)

        print(self.seg_index)

        cost = self.true_cost(self.y, self.seg_index)
        print(cost)
        print("\n-------------------------------- Active Segment Index Reduction Based Heuristic --------------------------")
        best_global = {
            "q": q_sol.copy(), "h": h_sol.copy(),
            "y": y_sol.copy(), "z": z_sol.copy(),
            "seg": self.seg_index.copy(), "cost": self.best_cost
        }

        self.iteration = 1
        best_global["q"],best_global["h"],best_global["y"],best_global["z"] = self.sequential_segment_cut_heuristic()

        #
        # ── Convergence Plot ──────────────────────────────────────────────────
        plt.figure()
        iters_org = list(self.Z_original.keys())
        objs_org = list(self.Z_original.values())
        plt.scatter(iters_org, objs_org, marker="^", s=30, label="Original NLP objective")

        iters_best = list(self.Z_best.keys())
        objs_best = [self.Z_best[i] for i in iters_best]
        plt.scatter(iters_best, objs_best, color="red", marker="*", s=80, label="Best objective found")

        for i, z in zip(iters_best, objs_best):
            plt.annotate(
                str(i), (i, z),
                textcoords="offset points", xytext=(0, 6),
                ha="center", fontsize=8,
            )

        plt.xlabel("Iteration")
        plt.ylabel("Objective value")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.close()

        # ── Final Report ──────────────────────────────────────────────────────
        print("\n" + "=" * 80)
        print("Final Best Results")
        print("=" * 80)
        print(f"Water Network:            {self.data_list[self.data_number]}")
        print(f"Final best objective:     {self.current_cost}")
        print(f"NLP problems solved:      {self.number_of_nlp}")
        print(f"Total iterations:         {self.iteration}")
        self.elapsed_time = time.time() - self.start_time
        print(f"Solver time:              {self.solver_time:.2f} s")
        print(f"Total elapsed time:       {self.elapsed_time:.2f} s")
        print("=" * 80)

# ══════════════════════════════════════════════════════════════════════════════
# Entry Point
# ══════════════════════════════════════════════════════════════════════════════

def _select_model_file(data_number: int) -> str:
    """Return the appropriate AMPL model filename for the given network index."""
    model_map = {5: "newyork_model.mod", 6: "blacksburg_model.mod"}
    return model_map.get(data_number, "wdnmodel.mod")


if __name__ == "__main__":
    data_number = int(sys.argv[1]) - 1
    network_name = DATA_LIST[data_number]
    data_file = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{network_name}.dat"
    model_file = _select_model_file(data_number)

    print(f"Water Network: {network_name}")

    optimizer = WaterNetworkOptimizer(model_file, data_file, data_number, DATA_LIST)
    optimizer.run()
