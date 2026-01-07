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
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

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

        # Regression-guided selection data & hyperparams
        self.training_X = []    # list of feature vectors (np arrays)
        self.training_y = []    # observed improvements (old_cost - new_cost)
        self.src_model = None   # dict: { 'beta_std', 'model', 'scaler' }
        self.min_samples_for_model = 20
        self.retrain_every = 10
        self._samples_since_retrain = 0

        self.bound_push = 0.1
        self.bound_frac = 0.1

    # ----------------- small utilities -----------------
    @contextlib.contextmanager
    def suppress_output(self):
        # Lightweight no-op context manager to match prior usage
        try:
            yield
        finally:
            pass
    # ----------------- IO / model loading -----------------
    def load_model(self):
        """Load the model and data."""
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)
        self.nodes = self.ampl.getSet('nodes').to_list()
        self.source = self.ampl.getSet('Source').to_list()
        self.arcs = self.ampl.getSet('arcs').to_list()
        self.pipes = self.ampl.getSet('pipes').to_list()
        if self.data_number == 6:
            self.fixarcs = set(self.ampl.getSet('fixarcs').to_list())
        # Parameters to dicts
        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
        # optional
        try:
            self.eps = self.ampl.getParameter('eps').getValues().to_dict()
        except Exception:
            self.eps = {}
        self.delta = 0.1
        self.p = 1.852
        self.omega = 10.67

    # ----------------- graph helpers -----------------
    def generate_random_acyclic_from_solution(self, q):
        # Build directed graph using sign of flow in q
        G = nx.DiGraph()
        G.add_nodes_from(self.nodes)
        for (i, j) in self.arcs:
            # q may be a dict-like; ensure safe indexing
            v = q.get((i, j), 0.0) if isinstance(q, dict) else q[i, j]
            if v >= 0:
                G.add_edge(i, j)
            else:
                G.add_edge(j, i)
        # Remove cycles if any by simple break (keep it simple)
        try:
            while True:
                cycle = nx.find_cycle(G, orientation='original')
                if not cycle:
                    break
                # remove one edge from cycle (last)
                u, v, _ = cycle[-1]
                if G.has_edge(u, v):
                    G.remove_edge(u, v)
        except Exception:
            pass
        return G

    def format_indian_number(self,num):
        num_str = str(num)
        if len(num_str) <= 3:
            return num_str
        else:
            last_three = num_str[-3:]
            remaining = num_str[:-3]
            remaining = ','.join([remaining[max(i - 2, 0):i] for i in range(len(remaining), 0, -2)][::-1])
            return remaining + ',' + last_three

    def fix_leaf_arc_flow(self):
        graph = nx.Graph()
        arc_set = self.ampl.getSet('arcs').to_list()
        graph.add_edges_from(arc_set)
        D = self.ampl.getParameter('D').getValues().to_dict()
        source = self.ampl.getSet('Source').to_list()[0]
        fixed_arcs = set()
        try:
            Qmax = self.ampl.getParameter('Q_max').getValues().to_list()[0]
        except Exception:
            Qmax = 0
        D[source] = -Qmax

        while True:
            leaf_nodes = [node for node in graph.nodes if graph.degree[node] == 1]
            if not leaf_nodes:
                break
            for leaf in leaf_nodes:
                neighbor = next(graph.neighbors(leaf))
                if (neighbor, leaf) in arc_set:
                    edge = (neighbor, leaf)
                    if edge not in fixed_arcs:
                        if leaf == source:
                            flow_value = D[leaf]
                            D[neighbor] = (D[neighbor] + flow_value)
                            source = neighbor
                        else:
                            flow_value = D[leaf]
                            D[neighbor] = D[neighbor] + flow_value
                        fixed_arcs.add(edge)
                    graph.remove_node(leaf)
                else:
                    edge = (leaf, neighbor)
                    if edge not in fixed_arcs:
                        if leaf == source:
                            flow_value = -D[leaf]
                            D[neighbor] = D[neighbor] - flow_value
                            source = neighbor
                        elif neighbor == source:
                            flow_value = -D[leaf]
                            D[neighbor] = D[neighbor] - D[leaf]
                        else:
                            flow_value = -D[leaf]
                            D[neighbor] += -flow_value
                        fixed_arcs.add(edge)
                    graph.remove_node(leaf)
        return fixed_arcs

    def is_cycle(self, graph, start_node, end_node, visited_copy, parent):
        visited_copy[start_node] = True
        for neighbor in graph.neighbors(start_node):
            if not visited_copy[neighbor]:
                isCycle = self.is_cycle(graph, neighbor, end_node, visited_copy, start_node)
                if isCycle:
                    return True
            else:
                if parent != neighbor:
                    if end_node == neighbor:
                        return True
        return False

    def presolve(self, graph, node, visited, parent, set_arc):
        visited_copy = visited.copy()
        isCycle = self.is_cycle(graph, node, node, visited_copy, parent)
        visited[node] = True
        if isCycle:
            for neighbor in graph.neighbors(node):
                if parent != neighbor:
                    set_arc.append((node, neighbor))
            return set_arc
        else:
            for neighbor in graph.neighbors(node):
                if parent != neighbor:
                    set_arc.append((node, neighbor))
                    self.presolve(graph, neighbor, visited, node, set_arc)
        return set_arc

    def fix_arc_set(self):
        graph = nx.Graph()
        arc_set = self.ampl.getSet('arcs').to_list()
        graph.add_edges_from(arc_set)
        visited = {node: False for node in graph.nodes()}
        source = self.ampl.getSet('Source').to_list()[0]
        set_arc = []
        set_ = self.presolve(graph, source, visited, -1, set_arc)
        return set_


    # place imports near top of your file if not already present
    # from sklearn.preprocessing import StandardScaler
    # from sklearn.linear_model import Ridge
    # import numpy as np
    
    # -----------------------------
    # 1) Helper: find current pipe index for an arc
    # -----------------------------
    def _current_pipe_index_for_arc(self, i, j, tol=1e-6):
        """Return the pipe index(s) k where l[i,j,k] > tol (selected in solution).
           If none, return k with largest l fraction (fallback)."""
        ks = [k for (u, v, k), val in self.l.items() if (u, v) == (i, j) and val > tol]
        if ks:
            # if multiple, return the one with largest length
            ks_sorted = sorted(ks, key=lambda k: self.l.get((i, j, k), 0.0), reverse=True)
            return ks_sorted[0]
        # fallback: choose the pipe with max d value index used in problem (not ideal)
        try:
            return max(self.pipes)
        except Exception:
            return list(self.pipes)[-1]
    
    # -----------------------------
    # 2) Unit cost derivative approximation (discrete commercial pipes)
    # -----------------------------
    def _approx_cost_derivative_wrt_d(self, i, j):
        """
        Approximate dC/dd for arc (i,j) by looking at the currently selected pipe
        and the next smaller commercial pipe.
        Returns derivative and delta_d used for finite diff (positive).
        """
        k_cur = self._current_pipe_index_for_arc(i, j)
        # find index of k_cur in sorted pipes by diameter value
        pipe_list = sorted(self.pipes, key=lambda k: self.d[k], reverse=True)  # largest -> smallest
        # find position of current pipe (choose pipe with matching index if present)
        # We must match by diameter value; so find k in pipe_list with same d as k_cur
        try:
            d_cur = float(self.d[k_cur])
        except Exception:
            d_cur = float(self.d[pipe_list[0]])
    
        # find all pipe indices sorted by diameter descending -> we want next smaller diameter
        sorted_by_d = sorted(self.pipes, key=lambda k: float(self.d[k]), reverse=True)
        # find the index of the pipe whose diameter is (approximately) d_cur
        pos = None
        for idx, k in enumerate(sorted_by_d):
            if abs(float(self.d[k]) - d_cur) < 1e-9:
                pos = idx
                break
        if pos is None:
            # fallback: consider smallest pipe as next if cannot find
            pos = 0
    
        # next smaller: move to next position (pos+1) if exists else None
        if pos + 1 < len(sorted_by_d):
            k_next = sorted_by_d[pos + 1]
            d_next = float(self.d[k_next])
            cost_cur = float(self.C[k_cur])
            cost_next = float(self.C[k_next])
            # derivative approx (cost_next - cost_cur) / (d_next - d_cur)
            if abs(d_next - d_cur) > 1e-12:
                dC_dd = (cost_next - cost_cur) / (d_next - d_cur)
                delta_d = (d_cur - d_next)  # positive amount diameter would decrease when switching to next smaller
            else:
                dC_dd = 0.0
                delta_d = 0.0
        else:
            # no smaller pipe exists; derivative ~ 0
            dC_dd = 0.0
            delta_d = 0.0
    
        return dC_dd, delta_d, k_cur
    
    # -----------------------------
    # 3) Compute F(q, eps) factor used in con2
    # -----------------------------
    def _F_of_q_eps(self, q_val, eps_val):
        """Compute the q-dependent multiplier F(q,eps) used in your smooth con2."""
        # handle zeros robustly
        q = float(q_val)
        eps = float(eps_val) if eps_val is not None else 0.0
        # Use same formula as in AMPL: F = q^3 * (q^2 + eps^2)^0.426 / (q^2 + 0.426*eps^2)
        denom = (q*q + 0.426 * eps * eps)
        if abs(denom) < 1e-16:
            return 0.0
        return (q**3) * ((q*q + eps*eps)**0.426) / denom
    
    # -----------------------------
    # 4) Compute sensitivity targets for all arcs
    # -----------------------------
    def compute_sensitivity_targets(self):
        """
        For each arc (i,j) compute approximate Δcost if diameter reduced to next smaller commercial pipe.
        Returns dict targets[(i,j)] = approx_delta_cost (positive means cost increase? See sign convention).
        We'll compute target as old_cost - new_cost (positive when improvement), so
           target = - (dC/dd + lambda * dh/dd) * delta_d
        because reducing d by delta_d yields cost change ≈ (dC/dd + lambda * dh/dd) * (-delta_d)
        but we prefer target = old - new = -Δcost (so positive = improvement).
        """
        targets = {}
        # get duals for con2 (headloss) from your stored all_duals or query AMPL
        try:
            # prefer stored all_duals if available
            dual_dict = self.all_duals.get("con2").to_dict() if hasattr(self, 'all_duals') and "con2" in self.all_duals else self.ampl.get_constraint('con2').getValues().to_dict()
        except Exception:
            # fallback: try to read constraints fresh
            dual_dict = {}
            for con_name, val in self.ampl.get_constraints():
                if con_name == "con2":
                    dual_dict = val.getValues().to_dict()
                    break
    
        for (i, j) in self.arcs:
            # flows, heads
            q_val = float(self.q.get((i, j), 0.0))
            hdiff = float(self.h.get(i, 0.0) - self.h.get(j, 0.0))
            eps_val = float(self.eps.get((i, j), 0.0)) if hasattr(self, 'eps') and (i,j) in getattr(self,'eps',{}) else 0.0
    
            # compute F = q-dependent factor
            F = self._F_of_q_eps(q_val, eps_val)
    
            # compute S(i,j) constants: for each pipe in arc, const = omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87)
            # We'll need l_i_j_k for the currently selected pipe index
            # find current pipe index k_cur
            try:
                k_cur = self._current_pipe_index_for_arc(i, j)
            except Exception:
                k_cur = None
    
            # approximate dC/dd and delta_d based on switching from current to next smaller pipe
            dC_dd, delta_d, k_cur = self._approx_cost_derivative_wrt_d(i, j)
    
            # if no delta (cannot reduce) then target zero
            if delta_d <= 0:
                targets[(i, j)] = 0.0
                continue
    
            # get l used for k_cur (length of pipe k_cur on this arc)
            l_len = self.l.get((i, j, k_cur), 0.0)
    
            # compute ∂S/∂d at k_cur: -4.87 * omega * l / (R^1.852 * d^(5.87))
            d_k = float(self.d[k_cur])
            R_k = float(self.R[k_cur])
            if d_k <= 0:
                targets[(i, j)] = 0.0
                continue
    
            const = (self.omega * l_len) / ( (R_k ** 1.852) * (d_k ** 4.87) )  # this is term in S
            # ∂S/∂d_k = -4.87 * (omega * l / (R^1.852 * d^(5.87)))
            dS_dd = -4.87 * (self.omega * l_len) / ( (R_k ** 1.852) * (d_k ** 5.87) )
    
            # ∂h/∂d_k = F * ∂S/∂d_k
            dh_dd = F * dS_dd
    
            # lambda for headloss constraint for this arc (if missing, use 0)
            lambda_hl = float(dual_dict.get((i, j), 0.0)) if dual_dict is not None else 0.0
    
            # Total derivative of cost w.r.t diameter for arc approx:
            # dC/dd + lambda_hl * dh/dd
            total_derivative = dC_dd + lambda_hl * dh_dd
    
            # Delta cost for decreasing diameter by delta_d: Δcost ≈ total_derivative * (-delta_d)
            # So old_cost - new_cost = - Δcost = - total_derivative * (-delta_d) = total_derivative * delta_d
            target = - total_derivative * delta_d
            # sign: positive target means we expect old_cost - new_cost > 0 (improvement)
            targets[(i, j)] = target
    
        return targets
    
    # -----------------------------
    # 5) Convert targets into X,y training arrays
    # -----------------------------
    def build_training_data_from_targets(self, targets):
        """
        Build X,y arrays for regression from the targets dict and current solution.
        Choose features you find predictive — here's a recommended set.
        """
        X = []
        y = []
        for (i, j), targ in targets.items():
            q_val = float(self.q.get((i, j), 0.0))
            headloss = abs(float(self.h.get(i, 0.0) - self.h.get(j, 0.0)))
            Lij = float(self.L.get((i, j), 0.0))
            # current pipe diameter
            try:
                k_cur = self._current_pipe_index_for_arc(i, j)
                d_cur = float(self.d[k_cur])
            except Exception:
                d_cur = float(self.d[min(self.pipes)])
            # dual for con2:
            try:
                dual_val = float(self.all_duals["con2"].to_dict().get((i, j), 0.0))
            except Exception:
                dual_val = 0.0
    
            # feature vector
            features = [
                abs(q_val),
                headloss,
                Lij,
                d_cur,
                dual_val,
                q_val * headloss  # interaction
            ]
            X.append(features)
            y.append(targ)
        return np.array(X, dtype=float), np.array(y, dtype=float)
    
    # -----------------------------
    # 6) Train standardized regression and return beta_std, model, scaler
    # -----------------------------
    def fit_regression_and_get_standardized_beta(self, X, y, alpha=1e-3):
        """
        Fit Ridge on standardized X and standardized y and return standardized coefs.
        """
        # filter degenerate case
        if X.size == 0 or y.size == 0:
            return None, None, None
    
        scaler = StandardScaler()
        Xs = scaler.fit_transform(X)
        # standardize y
        y_mean = y.mean()
        y_std = y.std() if y.std() > 0 else 1.0
        ys = (y - y_mean) / y_std
    
        model = Ridge(alpha=alpha)
        model.fit(Xs, ys)
        beta_std = model.coef_.copy()
    
        # Return beta_std (length = n_features), model, scaler (X scaler only)
        return beta_std, model, scaler

    # ----------------- Regression helpers -----------------
    def compute_features_for_arc(self, arc):
        """
        Build feature vector for arc (i,j) using Method-1 features:
          - |q[i,j]| : absolute flow magnitude
          - headloss : |h[i] - h[j]|
          - L_ij : total arc length
          - d_cur : current diameter index (max of l[i,j,k] > 0)
          - dual_val : dual of headloss constraint (con2)
          - interaction : q_val * headloss
        """
        i, j = arc
    
        # absolute flow
        q_val = abs(self.q.get((i, j), 0.0))
    
        # head difference
        head_diff = self.h.get(i, 0.0) - self.h.get(j, 0.0)
        headloss = abs(head_diff)
    
        # arc length
        L_ij = self.L.get((i, j), 0.0)
    
        # current diameter index
        ks = [k for (u, v, k), val in self.l.items() if (u, v) == (i, j) and val > 1e-9]
        d_cur = max(ks) if ks else 1.0
    
        # dual of headloss constraint
        dual_val = 0.0
        try:
            dual_dict = self.all_duals.get("con2", {}).to_dict()
            dual_val = float(dual_dict.get((i, j), 0.0))
        except Exception:
            dual_val = 0.0
    
        # interaction term
        interaction = q_val * headloss
    
        # final feature vector
        features = np.array([q_val, headloss, L_ij, float(d_cur), dual_val, interaction], dtype=float)
    
        return features

    # ----------------- Regression helpers -----------------
    def compute_features_for_arc1(self, arc):
        """
        Build feature vector for arc (i,j).
        Features:
          - dual of headloss constraint (con2)
          - |q[i,j]| flow magnitude
          - abs(head diff) |h[i] - h[j]|
          - sensitivity proxy: (h_i - h_j) / |q|^2.852 (if q nonzero)
          - cost contribution of arc (sum C[k]*l[i,j,k])
          - current max diameter index (arc_max)
          - reducible (binary) (arc_max > 1)
          - visited (binary)
        """
        i, j = arc
        # duals
        dual_val = 0.0
        try:
            dual_dict = self.all_duals["con2"].to_dict()
            dual_val = float(dual_dict.get((i, j), 0.0))
        except Exception:
            dual_val = 0.0
        # flows and heads
        q_val = abs(self.q.get((i, j), 0.0))
        headdiff = self.h.get(i, 0.0) - self.h.get(j, 0.0)
        abs_headdiff = abs(headdiff)

        # arc_max_dia
        arc_max = 1
        ks = [k for (u, v, k), val in self.l.items() if (u, v) == (i, j) and val > 1e-9]
        if ks:
            arc_max = max(ks)

        sens = 0.0
        if q_val > 1e-9:
            sens = (headdiff) / (q_val ** 2.852)

        cost_contrib = 0.0
        for k in self.pipes:
            cost_contrib += self.C.get(k, 0.0) * self.l.get((i, j, k), 0.0)

        reducible = 1.0 if arc_max > 1 else 0.0
        visited = 1.0 if (i, j) in getattr(self, 'visited_arc', []) else 0.0

        features = np.array([
            dual_val,
            q_val,
            abs_headdiff,
            sens,
            cost_contrib,
            float(arc_max),
            reducible,
            visited
        ], dtype=float)

        return features

    def fit_standardized_regression(self, X, y, alpha=1e-3):
        """
        Fit Ridge on standardized X and standardized y; return beta_std, model, scaler.
        beta_std is model.coef_ (for standardized X/y).
        """
        X = np.array(X, dtype=float)
        y = np.array(y, dtype=float)
        scaler = StandardScaler()
        X_std = scaler.fit_transform(X)

        # standardize y
        if y.std() == 0:
            y_std = (y - y.mean())
        else:
            y_std = (y - y.mean()) / y.std()

        model = Ridge(alpha=alpha)
        model.fit(X_std, y_std)

        beta_std = model.coef_.copy()  # standardized coefficients
        return beta_std, model, scaler

    def select_arc_with_src(self, candidate_arcs):
        """
        Predict improvement for candidate_arcs using current src_model.
        If not enough training data, fallback to dual-based ranking.
        Returns single arc (i,j) chosen.
        """
        if not candidate_arcs:
            return None

        # 1) compute sensitivity-based targets
        targets = self.compute_sensitivity_targets()
        
        # 2) build training data
        X_train, y_train = self.build_training_data_from_targets(targets)
        
        # Optionally print a quick diagnostic
        print("Built training data for regression: samples =", X_train.shape[0])
        
        # 3) fit regression
        beta_std, model, scaler = self.fit_regression_and_get_standardized_beta(X_train, y_train)
        self.src_model = {'beta_std': beta_std, 'model': model, 'scaler': scaler, 'X_train': X_train, 'y_train': y_train}
        
        # 4) you can inspect coefficients:
        if beta_std is not None:
            print("Standardized betas (features: |q|, headloss, L, d_cur, dual, q*headloss):")
            print(np.round(beta_std, 4))


        # fallback if insufficient data
        # if len(self.training_y) < self.min_samples_for_model:
        #     print("True")
        #     # dual dict may not exist; build from all_duals
        #     try:
        #         duals = self.all_duals["con2"].to_dict()
        #         arcs_sorted = sorted(candidate_arcs, key=lambda a: abs(duals.get(a, 0.0)), reverse=True)
        #         return arcs_sorted[0]
        #     except Exception:
        #         return candidate_arcs[0]
        #
        # # ensure model exists
        # if self.src_model is None:
        #     try:
        #         beta_std, model, scaler = self.fit_standardized_regression(self.training_X, self.training_y)
        #         self.src_model = {'beta_std': beta_std, 'model': model, 'scaler': scaler}
        #         self._samples_since_retrain = 0
        #     except Exception:
        #         # fallback
        #         try:
        #             duals = self.all_duals["con2"].to_dict()
        #             arcs_sorted = sorted(candidate_arcs, key=lambda a: abs(duals.get(a, 0.0)), reverse=True)
        #             return arcs_sorted[0]
        #         except Exception:
        #             return candidate_arcs[0]

        # predict improvement for each candidate arc
        best_arc = []
        best_pred = -1e99
        model = self.src_model['model']
        scaler = self.src_model['scaler']

        for arc in candidate_arcs:
            X_feat = self.compute_features_for_arc(arc).reshape(1, -1)
            Xs = scaler.transform(X_feat)
            # model predicts on standardized y; convert to original-scale improvement
            try:
                pred_std = model.predict(Xs)
                # We trained y as (old_cost - new_cost), so positive pred_std indicates improvement expected.
                # We don't have original y_std var here (since target std used in fit), but relative ordering is enough.
                pred = pred_std
            except Exception:
                pred = 0.0

            # small penalty if arc not reducible
            if X_feat[0, -2] == 0.0:
                pred -= 1e-6

            if pred > best_pred:
                best_pred = pred
                best_arc.append(arc)
            # best_arc = pred
        print("best_arc:", best_arc)
        # sorted_arcs = [arc for _, arc in sorted(zip(best_arc, candidate_arcs), key=lambda t: t[0], reverse=True)]
        # if model gives all zeros, fall back to dual
        if best_arc is None:
            try:
                duals = self.all_duals["con2"].to_dict()
                arcs_sorted = sorted(candidate_arcs, key=lambda a: abs(duals.get(a, 0.0)), reverse=True)
                return arcs_sorted[0]
            except Exception:
                return candidate_arcs[0]

        return best_arc[:5]

    # ----------------- diameter reduction (integrated with regression) -----------------
    def diameter_reduction(self):
        improved = False
        arc_max_dia = {}
        # build arc_max_dia mapping
        if self.data_number == 6:
            try:
                self.fixarcs = set(self.ampl.getSet('fixarcs').to_list())
            except Exception:
                self.fixarcs = set()
            for (i, j, d), val in self.l.items():
                if ((i, j) not in self.fixarcs) and ((j, i) not in self.fixarcs):
                    if val > 1e-3:
                        arc_max_dia[(i, j)] = max(arc_max_dia.get((i, j), d), d)
        else:
            for (i, j, d), val in self.l.items():
                if val > 1e-3:
                    arc_max_dia[(i, j)] = max(arc_max_dia.get((i, j), d), d)

        print("Iteration :", getattr(self, 'dia_red_iteration', 0))

        # collect duals
        self.all_duals = {}
        for con_name, val in self.ampl.get_constraints():
            self.all_duals[con_name] = val.getValues()

        dual_dict = self.all_duals.get("con2")
        if dual_dict is None:
            dual_dict = {}
        else:
            dual_dict = dual_dict.to_dict()
        # initial sorting by |dual|
        sorted_arcs = sorted(list(dual_dict.keys()), key=lambda kv: abs(dual_dict.get(kv, 0.0)), reverse=True)
        # filter visited etc.
        sorted_arcs = [arc for arc in sorted_arcs if arc not in getattr(self, 'visited_arc', [])]
        if self.data_number == 6 and hasattr(self, 'fixarcs'):
            sorted_arcs = [arc for arc in sorted_arcs if arc not in self.fixarcs]
        sorted_arcs = [arc for arc in sorted_arcs if arc not in getattr(self, 'fix_arc_set', [])]
        sorted_arcs = [arc for arc in sorted_arcs if arc_max_dia.get((arc[0], arc[1]), 1) != 1]

        if not sorted_arcs:
            print("No candidate arcs for diameter reduction.")
            return
        print("sorted_arcs:", sorted_arcs)
        # Use regression model to select the single best arc to try
        selected_arc = self.select_arc_with_src(sorted_arcs)
        print("selected_arcs:", selected_arc)
        if selected_arc is None:
            selected_arc = sorted_arcs[0]

        sorted_arcs = selected_arc
        print("sorted_arcs (selected):", sorted_arcs)

        print("----------------------------------------------------------------------------------------")
        print(f"{'Arc':<10}{'C_Best_Sol':<14}{'New_Sol':<14}"f"{'Solve_Time':<12}{'Solve_Result':<14}{'Improved':<10}{'Time':<12}")
        print("----------------------------------------------------------------------------------------")

        for (i, j) in sorted_arcs:
            # mark visited for now so we don't retry immediately
            if not hasattr(self, 'visited_arc'):
                self.visited_arc = []
            self.visited_arc.append((i, j))

            ampl = AMPL()
            ampl.reset()
            # read appropriate model file
            if self.data_number == 5:
                ampl.read("newyork_model.mod")
            elif self.data_number == 6:
                ampl.read("blacksburg_model.mod")
            else:
                ampl.read("wdnmodel.mod")
            ampl.read_data(self.data_file)

            # set initial points from current solution
            for (x, y, k), val in self.l.items():
                ampl.eval(f'let l[{x},{y},{k}] := {val};')
            for (x, y), val in self.q.items():
                ampl.eval(f'let q[{x},{y}] := {val};')
                if self.data_number == 5:
                    ampl.eval(f'let q1[{x},{y}] := {self.q1[x,y]};')
                    ampl.eval(f'let q2[{x},{y}] := {self.q2[x,y]};')
            for x, val in self.h.items():
                ampl.eval(f'let h[{x}] := {val};')

            # apply diameter reduction constraint: set larger diameters to zero
            for k in self.pipes:
                if k >= arc_max_dia[(i, j)]:
                    ampl.eval(f"subject to con3__{i}_{j}_{k}: l[{i},{j},{k}] = 0;")

            ampl.option['solver'] = "ipopt"
            ampl.option["ipopt_options"] = f"outlev = 0 expect_infeasible_problem = no bound_relax_factor = 0 tol = 1e-9 bound_push = {self.bound_push} bound_frac = {self.bound_frac} warm_start_init_point = yes halt_on_ampl_error = yes"

            with self.suppress_output():
                ampl.solve()

            # timing and bookkeeping
            solve_time = 0.0
            try:
                solve_time = ampl.get_value('_solve_elapsed_time')
            except Exception:
                solve_time = 0.0
            self.solver_time += solve_time
            self.number_of_nlp += 1

            # read solution
            try:
                l1 = ampl.getVariable('l').getValues().to_dict()
                q1 = ampl.getVariable('q').getValues().to_dict()
                h1 = ampl.getVariable('h').getValues().to_dict()
            except Exception:
                l1, q1, h1 = {}, {}, {}

            try:
                total_cost = ampl.getObjective("total_cost").value()
            except Exception:
                total_cost = None

            # Record training sample: features for the arc and observed improvement
            try:
                X_feat = self.compute_features_for_arc((i, j))
                # if total_cost unavailable, treat as no improvement
                if total_cost is None:
                    y_obs = 0.0
                else:
                    y_obs = float(self.current_cost - total_cost)  # positive => improvement
                self.training_X.append(X_feat)
                self.training_y.append(y_obs)
                self._samples_since_retrain += 1
            except Exception as e:
                # in case of any failure, skip recording this sample
                # but do not stop execution
                print("Warning: failed to record training sample:", e)

            # periodically retrain model
            if self._samples_since_retrain >= self.retrain_every and len(self.training_y) >= self.min_samples_for_model:
                try:
                    beta_std, model, scaler = self.fit_standardized_regression(self.training_X, np.array(self.training_y))
                    self.src_model = {'beta_std': beta_std, 'model': model, 'scaler': scaler}
                    self._samples_since_retrain = 0
                    # print learned betas for debugging
                    print("SRC betas:", np.round(beta_std, 4))
                except Exception as e:
                    print("Warning: failed to retrain SRC model:", e)

            # process solve result
            solved_flag = False
            try:
                solved_flag = (ampl.solve_result == "solved")
            except Exception:
                # fallback if attribute not found
                solved_flag = (getattr(ampl, 'solve_result', 'solved') == "solved")

            if solved_flag and (total_cost is not None):
                if total_cost < self.current_cost:
                    # improvement accepted
                    print(f"{str((i,j)):<10}"
                          f"{self.format_indian_number(round(self.current_cost)):<14}"
                          f"{self.format_indian_number(round(total_cost)):<14}"
                          f"{(str(round(solve_time, 2)) + 's'):<12}"
                          f"{ampl.solve_result if hasattr(ampl, 'solve_result') else 'solved':<14}{'Yes':<10}"
                          f"{round(time.time() - self.start_time, 2)}s")
                    # update current best
                    self.current_cost = total_cost
                    self.ampl = ampl
                    improved = True
                    self.network_graph = self.generate_random_acyclic_from_solution(q1)
                    self.best_acyclic_flow = self.network_graph.copy()
                    self.indegree_2_or_more = [node for node, indeg in self.best_acyclic_flow.in_degree() if indeg >= 2]
                    self.l = l1
                    self.q = q1
                    self.h = h1
                    if self.data_number == 5:
                        try:
                            self.q1 = ampl.getVariable('q1').getValues().to_dict()
                            self.q2 = ampl.getVariable('q2').getValues().to_dict()
                        except Exception:
                            pass
                    # compute head violations (debug)
                    headvio = 0
                    for u in self.nodes:
                        if self.h.get(u, 0.0) >= self.E.get(u, 0.0) + self.P.get(u, 0.0):
                            pass
                        else:
                            headvio += (self.E.get(u, 0.0) + self.P.get(u, 0.0) - self.h.get(u, 0.0))
                    print("headvio:", headvio)
                    self.fix_arc_set = list(set(self.super_source_out_arc) | set(getattr(self, 'fix_arc_set', [])))
                    print("----------------------------------------------------------------------------------------")
                else:
                    print(f"{str((i,j)):<10}"
                          f"{self.format_indian_number(round(self.current_cost)):<14}"
                          f"{self.format_indian_number(round(total_cost)):<14}"
                          f"{(str(round(solve_time, 2)) + 's'):<12}"
                          f"{ampl.solve_result if hasattr(ampl, 'solve_result') else 'not-solved':<14}{'No':<10}"
                          f"{round(time.time() - self.start_time, 2)}s")
            else:
                print(f"{str((i,j)):<10}"
                      f"{self.format_indian_number(round(self.current_cost)):<14}"
                      f"{self.format_indian_number(round(total_cost)) if total_cost is not None else 'N/A':<14}"
                      f"{(str(round(solve_time, 2)) + 's'):<12}"
                      f"{ampl.solve_result if hasattr(ampl, 'solve_result') else 'error':<14}{'No':<10}"
                      f"{round(time.time() - self.start_time, 2)}s ")

            if improved:
                self.dia_red_iteration = getattr(self, 'dia_red_iteration', 0) + 1
                # recursively try further reductions (matches your original behavior)
                self.diameter_reduction()
                break

    def solve(self):
        self.ampl.option["solver"] = "ipopt"
        self.ampl.set_option("ipopt_options", f"outlev = 0 tol = 1e-9 bound_relax_factor=0  bound_push = {self.bound_push} bound_frac = {self.bound_frac} halt_on_ampl_error = yes warm_start_init_point = no expect_infeasible_problem = no")   #max_iter = 1000
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


    # ----------------- main run -----------------
    def run(self):
        """Main function to run the Heuristic Approach."""
        self.start_time = time.time()
        self.bound_push , self.bound_frac = (0.1, 0.1)
        self.mu_init = 0.1
        print("NLP Model: Smooth Approximatate WDN Design Model 2, Epsilon Selection using Relative Error")
        print("NLP Solver: Ipopt")
        print("************************Initial Solution of Approximate WDN Design Model**********************")
        self.load_model()
        fix_arc_set = self.fix_leaf_arc_flow()
        print("fix_arc_set:",fix_arc_set)
        self.super_source_out_arc = self.fix_arc_set()
        print("super_source_out_arc:", self.super_source_out_arc, "\n")

        # initial solve (use self.ampl from load_model)
        self.solve()
        print("Objective: ",self.total_cost)
        print("Solve_result: ",self.solve_result)
        try:
            print("Solve_time:", self.ampl.get_value('_solve_elapsed_time'),"\n")
        except Exception:
            print("\n")
        if self.solve_result != 'solved':
            print("IPOPT did not solve the initial problem optimally. Exiting Heuristic.")
            sys.exit()
        self.current_cost = self.total_cost
        try:
            self.l = self.ampl.getVariable('l').getValues().to_dict()
            self.q = self.ampl.getVariable('q').getValues().to_dict()
            self.h = self.ampl.getVariable('h').getValues().to_dict()
        except Exception:
            self.l, self.q, self.h = {}, {}, {}
        # self.ampl.eval("display l;")
        # self.ampl.eval("display q;")
        # self.ampl.eval("display h;")
        # self.ampl.eval("display {i in nodes} E[i] + P[i];")
        # self.ampl.eval("display {i in nodes} h[i] - E[i] - P[i];")
        if self.data_number==5:
            try:
                self.q1 = self.ampl.getVariable('q1').getValues().to_dict()
                self.q2 = self.ampl.getVariable('q2').getValues().to_dict()
            except Exception:
                pass
        print("*****************************Improve the Initial Solution*************************************\n")
        self.super_source_out_arc = self.fix_arc_set()
        self.network_graph = self.generate_random_acyclic_from_solution(self.q)
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set)
        self.best_acyclic_flow = self.network_graph.copy()
        self.visited_nodes = []
        self.sorted_nodes = []

        print("\n----------------------------Diameter Reduction Approach------------------------------------")
        self.dia_red_iteration = self.headloss_increase_iteration + 1 if hasattr(self, 'headloss_increase_iteration') else 1
        self.visited_arc = []
        self.diameter_reduction()

        print("\n************************************Final Best Results*****************************************")
        print("Water Network:", self.data_list[self.data_number])
        try:
            self.eps = self.ampl.getParameter('eps').getValues().to_dict()
        except Exception:
            self.eps = {}
        print(f"Final best objective: {self.current_cost}")
        print("Number of nlp problem solved:", self.number_of_nlp)
        print("Total number of iteration:", getattr(self, 'iteration', 0) + getattr(self, 'headloss_increase_iteration', 0) + getattr(self, 'dia_red_iteration', 0))
        elapsed_time = time.time() - self.start_time
        solver_time = self.solver_time
        print(f"Solver_time: {solver_time:.2f} seconds")
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

    data_number = int(sys.argv[1]) - 1
    input_data_file = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{data_list[(data_number)]}.dat"
    print("Water Network:", data_list[(data_number)])
    if data_number==5:
        optimizer = WaterNetworkOptimizer("newyork_model.mod", input_data_file, data_number, data_list)
    elif data_number == 6:
        optimizer = WaterNetworkOptimizer("blacksburg_model.mod", input_data_file, data_number, data_list)
    else:
        optimizer = WaterNetworkOptimizer("wdnmodel.mod", input_data_file, data_number, data_list)
    optimizer.run()

