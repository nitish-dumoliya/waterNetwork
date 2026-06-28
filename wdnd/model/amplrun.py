import time
import sys
from tabulate import tabulate
from amplpy import AMPL
import contextlib
import os
import numpy as np
import math
import subprocess, os, tempfile

from collections import defaultdict

class WaterNetworkSolver:
    def __init__(self, solver_name, data_file, data_number):
        # self.model_file = model_file
        self.data_file = data_file
        self.data_number = data_number

        self.solver_name = solver_name
        self.ampl = AMPL()
        # To store solutions
        self.q_init = {}
        self.h_init = {}
        self.l_init = {}
        self.eps_init = {}
        self.q1_init = {}
        self.q2_init = {}
        self.eps = {}

    def read_model_and_data(self):

        if self.data_number == 3:
            self.model_file = "trn_model.mod"
        elif self.data_number == 4:
            self.model_file = "newyork_model.mod"
        elif self.data_number == 7:
            self.model_file = "blacksburg_model.mod"
        elif self.data_number == 9:
            self.model_file = "bakryun_model.mod"
        elif self.data_number == 13:
            self.model_file = "farhadgerd_model.mod"
        elif self.data_number==24:
            self.model_file = "kl_model.mod"

        # elif self.data_number == 22:
        #     self.model_file = "marchi_rural_model.mod"
        else:
            self.model_file = "wdnmodel.mod"

        print("WDN Model Name:", self.model_file)

        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)

        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')

        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()

        if self.data_number==3:
            self.fixarcs = self.ampl.getSet("fixed_arcs")
            self.unfixed_arcs = self.ampl.getSet("unfixed_arcs")
            self.parallel_arcs = self.ampl.getSet("parallel_arcs")
            self.exroughness = self.ampl.getParameter("exroughness").to_dict()
            self.exdiam = self.ampl.getParameter("exdiam").to_dict()


        if self.data_number in [7, 24]:
            self.fixarcs = self.ampl.getSet("fixarcs")
 
        if self.data_number==9:
            self.fixarcs = self.ampl.getSet("fixarcs")
            self.unfixed_arcs = self.ampl.getSet("unfixed_arcs")
            self.parallel_arcs = self.ampl.getSet("parallel_arcs")

        if self.data_number==13:
            self.parallel_arcs = self.ampl.getSet("parallel_arcs")

        if self.data_number==24:
            self.fix_L = self.ampl.getParameter("fix_L").to_dict()
            print(self.fix_L)

        # if self.data_number==24:
        #     fix_arcs = {(477, 433):661.194101403913, (1476, 1477):1571.93425823395 }
        #     for kay, len in fix_arcs.items():
        #         (i,j) = (key[0],key[1])
        #         self.ampl.eval(f"s.t. con_headloss_{i}_{j}:h[{i}] - h[{j}]  = (q[{i},{j}])^3 *((((q[{i},{j}])^2 + eps[{i},{j}]^2)^0.426) /((q[{i},{j}])^2 + 0.426*eps[{i},{j}]^2)) * sum{{k in pipes}}(omega * l[{i},{j},k]/(R[k]^1.852 * d[k]^4.87));")

        # if self.data_number==15:
        #     non_fix_arcs = [(99,1), (1,2), (2,3), (3,4), (4,6), (6,7), (7,11), (11, 16), (16, 17)]
        #     for (i,j) in non_fix_arcs:
        #         # self.ampl.eval(f"s.t. con_headloss_{i}_{j}:h[{i}] - h[{j}]  = (q[{i},{j}])^3 *((((q[{i},{j}])^2 + eps[{i},{j}]^2)^0.426) /((q[{i},{j}])^2 + 0.426*eps[{i},{j}]^2)) * sum{{k in pipes}}(omega * l[{i},{j},k]/(R[k]^1.852 * d[k]^4.87));")
        #         self.ampl.eval(
        #             f"s.t. con_headloss_{i}_{j}: "
        #             f"h[{i}] - h[{j}] = "
        #             f"(q[{i},{j}])^3 * "
        #             f"((((q[{i},{j}])^2 + eps[{i},{j}]^2)^0.426) / "
        #             f"((q[{i},{j}])^2 + 0.426*eps[{i},{j}]^2)) * "
        #             f"sum{{k in pipes}}(omega * l[{i},{j},k]/(R[k]^1.852 * d[k]^4.87));"
        #             )
        #         self.ampl.eval(f"""
        #             subject to con3__{i}_{j}: sum{{k in pipes}} l[{i},{j},k] = L[{i},{j}];
        #             subject to con4_{i}_{j}{{k in pipes}}: l[{i},{j},k] <= L[{i},{j}];
        #         """)

        # n_vars = self.ampl.getValue("_nvars")
        # n_cons = self.ampl.getValue("_ncons")
        # n_nl_cons = self.ampl.getValue("_snl")
        # n_nl_obj  = self.ampl.getValue("_nlobj")

        # print("Variables:", n_vars)
        # print("Constraints:", n_cons)
        # print("Nonlinear constraints:", n_nl_cons)
        # print("Nonlinear objective:", n_nl_obj)


    def compute_adaptive_eps(self, demand):
        """
        Adaptive epsilon calculation based on the minimum demand value.
        """
        epsilon = 1e-14 + 1e-10*(demand)**2
        #if max_demand <= 0.1:
        #    epsilon =  1e-13
        #elif max_demand >= 0.1 and max_demand <= 0.2:
        #epsilon = 1e-14
        #elif max_demand >= 0.2:
        #    epsilon = 1e-6

        return epsilon

    def constraint_violations_(self, q_values, h_values, l_values, epsilon, solver):
        """
        Compute and report all constraint violations for the current solution.
    
        Constraints checked:
            con1: Flow balance at non-source nodes
            con2: Head-loss (original vs. approximated)
            con3: Pipe length summation per arc
            con4: Individual pipe segment length upper bound
            con5: Source node head equality
            con6: Minimum pressure head at demand nodes
        """
        total_absolute_constraint_violation = 0.0

        # ------------------------------------------------------------------
        # Constraint 1: Flow balance at non-source nodes
        # ------------------------------------------------------------------
        con1_gap = {}

        if self.data_number == 5:
            q1 = self.ampl.get_variable('q1').get_values().to_dict()
            q2 = self.ampl.get_variable('q2').get_values().to_dict()
    
            for i in self.nodes:
                if i in self.source:
                    continue
                incoming = sum(q1[j, i] + q2[j, i] for j in self.nodes if (j, i) in self.arcs)
                outgoing = sum(q1[i, j] + q2[i, j] for j in self.nodes if (i, j) in self.arcs)
                violation = (incoming - outgoing) - self.D[i]
                con1_gap[str(i)] = violation
                total_absolute_constraint_violation += abs(violation)
        else:
            for i in self.nodes:
                if i in self.source:
                    continue
                incoming = sum(q_values[j, i] for j in self.nodes if (j, i) in self.arcs)
                outgoing = sum(q_values[i, j] for j in self.nodes if (i, j) in self.arcs)
                violation = (incoming - outgoing) - self.D[i]
                con1_gap[str(i)] = violation
                total_absolute_constraint_violation += abs(violation)
    
        # ------------------------------------------------------------------
        # Constraint 2: Head-loss constraint (original vs. approximated)
        # ------------------------------------------------------------------
        con2_original_gap = {}
        con2_approx_gap = {}
        absolute_violations = {}
        relative_violations = {}
    
        con2_original_violation = 0.0
        con2_approx_violation = 0.0
        con2_absolute_constraint_violation = 0.0
        con2_relative_constraint_violation = 0.0
    
        def _alpha_rhs(i, j):
            """Weighted head-loss resistance coefficient for arc (i, j)."""
            return sum(
                10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * (self.d[k] ** 4.87))
                for k in self.pipes
            )
    
        def _original_rhs(q, alpha):
            return q * (abs(q) ** 0.852) * alpha
    
        def _approx_rhs(q, eps, alpha):
            return (q**3 * ((q**2 + eps**2) ** 0.426) / (q**2 + 0.426 * eps**2)) * alpha
    
        def _record_con2(key, lhs, orig_rhs, apprx_rhs):
            """Update all con2 tracking dictionaries and running totals."""

            nonlocal con2_original_violation, con2_approx_violation
            nonlocal con2_absolute_constraint_violation, con2_relative_constraint_violation
            nonlocal total_absolute_constraint_violation

            con2_original_gap[key] = lhs - orig_rhs
            con2_approx_gap[key]   = lhs - apprx_rhs
    
            abs_viol = orig_rhs - apprx_rhs
            rel_viol = abs_viol / (orig_rhs + 1e-14)
    
            absolute_violations[key] = abs_viol
            relative_violations[key] = rel_viol
    
            con2_original_violation            += abs(lhs - orig_rhs)
            con2_approx_violation              += abs(lhs - apprx_rhs)
            con2_absolute_constraint_violation += abs(abs_viol)
            con2_relative_constraint_violation += abs(rel_viol)
            total_absolute_constraint_violation += abs(lhs - apprx_rhs)
    
        if self.data_number == 5:
            self.exdiam = self.ampl.getParameter('exdiam').to_dict()
    
            for (i, j) in q1.keys():
                lhs = 2 * (h_values[i] - h_values[j])
                alpha = sum(
                    10.67 * l_values[i, j, k] / ((self.R[i, j] ** 1.852) * (self.d[k] ** 4.87))
                    for k in self.pipes
                )
                hl_coeff = 10.67 * self.L[i, j] / (self.R[i, j] ** 1.852 * self.exdiam[i, j] ** 4.87)
    
                orig_rhs = (
                    q1[i, j] * (abs(q1[i, j]) ** 0.852) * hl_coeff
                    + q2[i, j] * (abs(q2[i, j]) ** 0.852) * alpha
                )
                apprx_rhs = (
                    _approx_rhs(q1[i, j], epsilon[i, j], hl_coeff)
                    + _approx_rhs(q2[i, j], epsilon[i, j], alpha)
                )
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
        elif self.data_number == 6:
            self.fixarcs  = self.ampl.getSet('fixarcs')
            self.fix_r    = self.ampl.getParameter('fix_r').to_dict()
            self.exdiam   = self.ampl.getParameter('fixdiam').to_dict()
    
            # Variable arcs
            for (i, j) in self.arcs:
                if (i, j) in self.fixarcs:
                    continue
                lhs      = h_values[i] - h_values[j]
                alpha    = _alpha_rhs(i, j)
                orig_rhs = _original_rhs(q_values[i, j], alpha)
                apprx_rhs = _approx_rhs(q_values[i, j], epsilon[i, j], alpha)
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
            # Fixed-diameter arcs
            for (i, j) in self.fixarcs:
                lhs = h_values[i] - h_values[j]
                hl_coeff = 10.67 * self.L[i, j] / (self.fix_r[i, j] ** 1.852 * self.exdiam[i, j] ** 4.87)
                orig_rhs  = _original_rhs(q_values[i, j], hl_coeff)
                apprx_rhs = _approx_rhs(q_values[i, j], epsilon[i, j], hl_coeff)
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
        else:
            for (i, j) in q_values:
                lhs       = h_values[i] - h_values[j]
                alpha     = _alpha_rhs(i, j)
                orig_rhs  = _original_rhs(q_values[i, j], alpha)
                apprx_rhs = _approx_rhs(q_values[i, j], epsilon[i, j], alpha)
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
        # ------------------------------------------------------------------
        # Constraint 3: Pipe lengths must sum to arc length
        # ------------------------------------------------------------------
        con3_gap = {}
        arcs_to_check = (
            [(i, j) for (i, j) in self.arcs if (i, j) not in self.fixarcs]
            if self.data_number == 6
            else list(self.arcs)
        )
        for (i, j) in arcs_to_check:
            violation = sum(l_values[i, j, k] for k in self.pipes) - self.L[i, j]
            con3_gap[f"{i},{j}"] = violation
            total_absolute_constraint_violation += abs(violation)
    
        # ------------------------------------------------------------------
        # Constraint 4: Individual pipe segment ≤ arc length
        # ------------------------------------------------------------------
        con4_gap = {}
        for (i, j) in self.arcs:
            for k in self.pipes:
                violation = max(0.0, l_values[i, j, k] - self.L[i, j])
                con4_gap[f"{i},{j},{k}"] = violation
                total_absolute_constraint_violation += violation
    
        # ------------------------------------------------------------------
        # Constraint 5: Source node head equality
        # ------------------------------------------------------------------
        con5_gap = {}
        for j in self.source:
            violation = h_values[j] - self.E[j]
            con5_gap[str(j)] = violation
            total_absolute_constraint_violation += abs(violation)
    
        # ------------------------------------------------------------------
        # Constraint 6: Minimum pressure head at demand nodes
        # ------------------------------------------------------------------
        con6_gap = {}
        for j in self.nodes:
            if j in self.source:
                continue
            violation = max(0.0, self.E[j] + self.P[j] - h_values[j])
            con6_gap[str(j)] = violation
            total_absolute_constraint_violation += violation
    
        # ------------------------------------------------------------------
        # Reporting
        # ------------------------------------------------------------------
        separator = "*" * 79
    
        # print(f"\n{separator}")
        # print("CONSTRAINT VIOLATIONS\n")
        #
        # table_data = (
        #     [(k, f"{v:.8f}") for k, v in con1_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con2_approx_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con3_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con4_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con5_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con6_gap.items()]
        # )
        # print(tabulate(table_data, headers=["Constraint ID", "Violation"], tablefmt="grid"))
        print(f"\n{separator}")
        print(f"\nTotal absolute constraint violation: {total_absolute_constraint_violation:.2e}")
    
        print("HEAD-LOSS CONSTRAINT: ORIGINAL vs. APPROXIMATION\n")
    
        table_data = []
        for key, rel_viol in relative_violations.items():
            i_str, j_str = key.split(',')
            i, j = int(i_str), int(j_str)
            table_data.append([
                key,
                q_values[i, j],
                f"{con2_original_gap[key]:.2e}",
                f"{con2_approx_gap[key]:.2e}",
                f"{absolute_violations[key]:.2e}",
                f"{rel_viol:.2e}",
            ])
    
        headers = [
            "Arc (i,j)",
            "Flow",
            "Original Violation",
            "Approx Violation",
            "Absolute Diff",
            "Relative Diff",
        ]
        # print(tabulate(table_data, headers=headers, tablefmt="grid"))
        print(f"\nSum of original head-loss violation : {con2_original_violation:.2e}")
        print(f"Sum of approx  head-loss violation  : {con2_approx_violation:.2e}")
        print(f"Sum of absolute diff (orig vs approx): {con2_absolute_constraint_violation:.2e}")
        print(f"Sum of relative diff (orig vs approx): {con2_relative_constraint_violation:.2e}")
        print(f"\n{separator}\n")

    def constraint_violations(self, q_values, h_values, l_values, epsilon):
        """
        Compute and report all constraint violations for the current solution.
    
        Constraints checked:
            con1 : Flow balance at non-source nodes
            con2 : Head-loss (original vs. approximated)
            con3 : Pipe length summation per arc
            con4 : Individual pipe segment length upper bound
            con5 : Source node head equality
            con6 : Minimum pressure head at demand nodes
    
        Supported data_number values: 3, 4, 7, 9, 24
        """
        total_absolute_constraint_violation = 0.0
    
        # ── Build arc-neighbour index once (avoids O(N²) membership tests) ────
        # incoming[j] = list of i where (i,j) in arcs
        # outgoing[i] = list of j where (i,j) in arcs
        arc_set = set(self.arcs)
        incoming = defaultdict(list)
        outgoing = defaultdict(list)
        for (i, j) in arc_set:
            incoming[j].append(i)
            outgoing[i].append(j)
    
        # ══════════════════════════════════════════════════════════════════════
        # Constraint 1 – Flow balance at non-source nodes
        # ══════════════════════════════════════════════════════════════════════
        con1_gap = {}
    
        if self.data_number == 4:
            q1 = self.q1
            q2 = self.q2
    
            # Accumulate net-flow per node in one pass over arcs
            net_flow = defaultdict(float)  # net_flow[i] = inflow - outflow at i
            for (i, j) in q1:
                f = q1[i, j] + q2[i, j]
                net_flow[j] += f
                net_flow[i] -= f
    
            for i in self.nodes:
                if i in self.source:
                    continue
                violation = net_flow[i] - self.D[i]
                con1_gap[str(i)] = violation
                total_absolute_constraint_violation += abs(violation)
    
        else:
            # Works for data_number 3, 7, 9, 24 (all use scalar q_values)
            net_flow = defaultdict(float)
            for (i, j), f in q_values.items():
                net_flow[j] += f
                net_flow[i] -= f
    
            for i in self.nodes:
                if i in self.source:
                    continue
                violation = net_flow[i] - self.D[i]
                con1_gap[str(i)] = violation
                total_absolute_constraint_violation += abs(violation)
    
        # ══════════════════════════════════════════════════════════════════════
        # Constraint 2 – Head-loss (original vs. approximated)
        # ══════════════════════════════════════════════════════════════════════
        con2_original_gap  = {}
        con2_approx_gap    = {}
        absolute_violations = {}
        relative_violations = {}
    
        con2_original_violation           = 0.0
        con2_approx_violation             = 0.0
        con2_absolute_constraint_violation = 0.0
        con2_relative_constraint_violation = 0.0
    
        # ── Precompute alpha for every (i,j,k) triplet once ──────────────────
        # alpha_cache[(i,j)] = sum_k  omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87)
        OMEGA = 10.67  # Hazen-Williams SI constant (same as param omega in AMPL)
    
        # Precompute per-pipe denominator  R[k]^1.852 * d[k]^4.87
        if self.data_number==4:
            # pipe_denom = self.alpha
            pipe_denom = {k: (100 ** 1.852) * (self.d[k] ** 4.87) for k in self.pipes}
        elif self.data_number==13:
            pipe_denom = {k: (130 ** 1.852) * (self.d[k] ** 4.87) for k in self.pipes}
        else:
            pipe_denom = {k: (self.R[k] ** 1.852) * (self.d[k] ** 4.87) for k in self.pipes}
    
        def _alpha_cache_build(arcs_iter):
            cache = {}
            for (i, j) in arcs_iter:
                cache[i, j] = sum(
                    OMEGA * l_values[i, j, k] / pipe_denom[k]
                    for k in self.pipes
                )
            return cache
    
        def _hl_original(q, alpha):
            """True Hazen-Williams head-loss: q |q|^0.852 * alpha"""
            return q * (abs(q) ** 0.852) * alpha
    
        def _hl_approx(q, eps, alpha):
            """Smooth approximation of head-loss."""
            q2e2 = q * q + eps * eps
            return (q ** 3) * (q2e2 ** 0.426) / (q * q + 0.426 * eps * eps) * alpha
    
        def _record_con2(key, lhs, orig_rhs, apprx_rhs):
            nonlocal con2_original_violation, con2_approx_violation
            nonlocal con2_absolute_constraint_violation, con2_relative_constraint_violation
            nonlocal total_absolute_constraint_violation
    
            con2_original_gap[key] = lhs - orig_rhs
            con2_approx_gap[key]   = lhs - apprx_rhs
    
            abs_viol = orig_rhs - apprx_rhs
            rel_viol = abs_viol / (orig_rhs + 1e-14)
    
            absolute_violations[key] = abs_viol
            relative_violations[key] = rel_viol
            # print(key, abs_viol) 
            con2_original_violation            += abs(lhs - orig_rhs)
            con2_approx_violation              += abs(lhs - apprx_rhs)
            con2_absolute_constraint_violation += abs(abs_viol)
            con2_relative_constraint_violation += abs(rel_viol)
            total_absolute_constraint_violation += abs(lhs - apprx_rhs)
    
        # ── data_number 4 ──────────────────────────────────────────────────────
        if self.data_number == 4:
            self.exdiam = self.ampl.getParameter('exdiam').to_dict()
    
            # alpha uses R[i,j] (arc-level roughness) — different from pipe R[k]
            for (i, j) in q1:
                lhs = 2.0 * (h_values[i] - h_values[j])
    
                alpha = sum(
                    OMEGA * l_values[i, j, k] / ((self.R[i, j] ** 1.852) * (self.d[k] ** 4.87))
                    for k in self.pipes
                )
                hl_coeff = (OMEGA * self.L[i, j]
                            / (self.R[i, j] ** 1.852 * self.exdiam[i, j] ** 4.87))
    
                orig_rhs  = (_hl_original(q1[i, j], hl_coeff)
                             + _hl_original(q2[i, j], alpha))
                apprx_rhs = (_hl_approx(q1[i, j], epsilon[i, j], hl_coeff)
                             + _hl_approx(q2[i, j], epsilon[i, j], alpha))
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
        # ── data_number 7 ──────────────────────────────────────────────────────
        elif self.data_number == 7:
            self.fixarcs = self.ampl.getSet('fixarcs')
            self.fix_r   = self.ampl.getParameter('fix_r').to_dict()
            self.exdiam  = self.ampl.getParameter('fixdiam').to_dict()
    
            var_arcs   = [(i, j) for (i, j) in self.arcs if (i, j) not in self.fixarcs]
            alpha_map  = _alpha_cache_build(var_arcs)
    
            for (i, j) in var_arcs:
                lhs       = h_values[i] - h_values[j]
                alpha     = alpha_map[i, j]
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(q_values[i, j], alpha),
                             _hl_approx(q_values[i, j], epsilon[i, j], alpha))
    
            for (i, j) in self.fixarcs:
                lhs      = h_values[i] - h_values[j]
                hl_coeff = (OMEGA * self.L[i, j]
                            / (self.fix_r[i, j] ** 1.852 * self.exdiam[i, j] ** 4.87))
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(q_values[i, j], hl_coeff),
                             _hl_approx(q_values[i, j], epsilon[i, j], hl_coeff))
    
        # ── data_number 3 ──────────────────────────────────────────────────────
        # Standard network: all arcs are variable-diameter, no fixed arcs.
        # Uses q_values (single flow per arc), l_values for pipe lengths.
        elif self.data_number == 3:
            var_arcs   = [(i, j) for (i, j) in self.unfixed_arcs]
            alpha_map = _alpha_cache_build(var_arcs)
            for (i, j) in var_arcs:
                lhs       = h_values[i] - h_values[j]
                alpha     = alpha_map[i, j]
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(q_values[i, j], alpha),
                             _hl_approx(q_values[i, j], epsilon[i, j], alpha))
            var_arcs   = [(i, j) for (i, j) in self.parallel_arcs]
            alpha_map = _alpha_cache_build(var_arcs)
            for (i, j) in var_arcs:
                lhs       = h_values[i] - h_values[j]
                alpha     = alpha_map[i, j]
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(self.q2[i, j], alpha),
                             _hl_approx(self.q2[i, j], epsilon[i, j], alpha))
                alpha1     = OMEGA * self.L[i,j] / ( (self.exroughness[i,j]**1.852) * (self.exdiam[i,j])**4.87)
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(self.q1[i, j], alpha1),
                             _hl_approx(self.q1[i, j], epsilon[i, j], alpha1))

            var_arcs   = [(i, j) for (i, j) in self.fixarcs if (i,j) not in self.parallel_arcs]
            # alpha_map = _alpha_cache_build(var_arcs)
            for (i, j) in var_arcs:
                lhs       = h_values[i] - h_values[j]
                alpha1     = OMEGA * self.L[i,j] / ( (self.exroughness[i,j]**1.852) * (self.exdiam[i,j])**4.87)
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(q_values[i, j], alpha1),
                             _hl_approx(q_values[i, j], epsilon[i, j], alpha1))


        # ── data_number 9 ──────────────────────────────────────────────────────
        # Mixed network: unfixed_arcs (variable diameter) + parallel_arcs
        # (existing pipe + new pipe in parallel), rest are fixed-diameter.
        # Mirrors the AMPL model in document 3 (set unfixed_arcs / parallel_arcs).
        elif self.data_number == 9:
            # Retrieve sets and parameters needed
            unfixed_arcs  = set(self.ampl.getSet('unfixed_arcs'))
            parallel_arcs = set(self.ampl.getSet('parallel_arcs'))
    
            # q1_9 = self.ampl.get_variable('q1').get_values().to_dict()  # flow on new pipe
            # q2_9 = self.ampl.get_variable('q2').get_values().to_dict()  # flow on existing pipe
            # l1_9 = {(i, j, k): self.ampl.get_variable('l1').get_values().to_dict().get((i, j, k), 0.0)
            #         for (i, j) in parallel_arcs for k in self.pipes}
            
            q1_9 = self.q1
            q2_9 = self.q2
            l1_9 = self.l1

            fix_r_9   = self.ampl.getParameter('fix_r').to_dict()
            fixdiam_9 = self.ampl.getParameter('fixdiam').to_dict()
    
            # Pre-cache alpha for unfixed and parallel arcs separately
            alpha_unfixed  = _alpha_cache_build(unfixed_arcs)
    
            # alpha for parallel arcs uses l1 (new pipe lengths)
            alpha_parallel = {}
            for (i, j) in parallel_arcs:
                alpha_parallel[i, j] = sum(
                    OMEGA * l1_9[i, j, k] / pipe_denom[k]
                    for k in self.pipes
                )
    
            # Unfixed arcs: single flow q[i,j], variable-diameter new pipe
            for (i, j) in unfixed_arcs:
                lhs   = h_values[i] - h_values[j]
                alpha = alpha_unfixed[i, j]
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(q_values[i, j], alpha),
                             _hl_approx(q_values[i, j], epsilon[i, j], alpha))

            # Parallel arcs: q1 on new pipe, q2 on existing (fixed-diameter) pipe
            for (i, j) in parallel_arcs:
                lhs = h_values[i] - h_values[j]
    
                alpha_new = alpha_parallel[i, j]
                hl_exist  = (OMEGA * self.L[i, j]
                             / (fix_r_9[i, j] ** 1.852 * fixdiam_9[i, j] ** 4.87))
    
                # con3: h_i - h_j = HL(q1, new pipe)
                _record_con2(f"{i},{j}_new", lhs,
                             _hl_original(q1_9[i, j], alpha_new),
                             _hl_approx(q1_9[i, j], epsilon[i, j], alpha_new))
    
                # con4: h_i - h_j = HL(q2, existing pipe)
                _record_con2(f"{i},{j}_exist", lhs,
                             _hl_original(q2_9[i, j], hl_exist),
                             _hl_approx(q2_9[i, j], epsilon[i, j], hl_exist))
    
            # Fixed arcs (fixarcs minus parallel_arcs): single flow on existing pipe
            fixarcs_9 = set(self.ampl.getSet('fixarcs'))
            fix_r_all = self.ampl.getParameter('fix_r').to_dict()
            exdiam_9  = self.ampl.getParameter('fixdiam').to_dict()
            for (i, j) in (fixarcs_9 - parallel_arcs):
                lhs      = h_values[i] - h_values[j]
                hl_coeff = (OMEGA * self.L[i, j]
                            / (fix_r_all[i, j] ** 1.852 * exdiam_9[i, j] ** 4.87))
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(q_values[i, j], hl_coeff),
                             _hl_approx(q_values[i, j], epsilon[i, j], hl_coeff))
 
        # ── data_number 4 ──────────────────────────────────────────────────────
        # arcs diff parallel_arcs : single new pipe  l[i,j,k],  flow q[i,j]
        # parallel_arcs           : two new pipes
        #                             pipe-1  l1[i,j,k], flow q1[i,j]
        #                             pipe-2  l2[i,j,k], flow q2[i,j]
        # R is arc-indexed: self.R[i,j]
        elif self.data_number == 13:
            parallel_arcs_4 = set(self.ampl.getSet('parallel_arcs'))
            non_par_4       = arc_set - parallel_arcs_4
 
            q1_4 = self.q1
            q2_4 = self.q2
            l1_4 = self.l1
            l2_4 = self.l2
 
            # pipe_denom_arc: per (arc, pipe) because R is arc-indexed
            def _alpha_arc(i, j, lmap):
                r_pow = self.R[i, j] ** 1.852
                return sum(OMEGA * lmap[i, j, k] / (r_pow * self.d[k] ** 4.87)
                           for k in self.pipes)
 
            # con2: single new pipe, flow q
            for (i, j) in non_par_4:
                lhs   = h_values[i] - h_values[j]
                alpha = _alpha_arc(i, j, l_values)
                _record_con2(f"{i},{j}",
                             lhs,
                             _hl_original(q_values[i, j], alpha),
                             _hl_approx(q_values[i, j], epsilon[i, j], alpha))
 
            # con3: parallel pipe-1, flow q1
            for (i, j) in parallel_arcs_4:
                lhs   = h_values[i] - h_values[j]
                alpha = _alpha_arc(i, j, l1_4)
                _record_con2(f"{i},{j}_pipe1",
                             lhs,
                             _hl_original(q1_4[i, j], alpha),
                             _hl_approx(q1_4[i, j], epsilon[i, j], alpha))
 
            # con4: parallel pipe-2, flow q2
            for (i, j) in parallel_arcs_4:
                lhs   = h_values[i] - h_values[j]
                alpha = _alpha_arc(i, j, l2_4)
                _record_con2(f"{i},{j}_pipe2",
                             lhs,
                             _hl_original(q2_4[i, j], alpha),
                             _hl_approx(q2_4[i, j], epsilon[i, j], alpha))

        # ── data_number 24 ─────────────────────────────────────────────────────
        # Mixed network: arcs (variable-diameter) + fixarcs (existing pipes,
        # with parallel new pipe in l1).  Mirrors the AMPL model in document 4.
        elif self.data_number == 24:
            fixarcs_24 = set(self.ampl.getSet('fixarcs'))
            # fix_r_24   = self.ampl.getParameter('fix_r').to_dict()
            # fixdiam_24 = self.ampl.getParameter('fixdiam').to_dict()
            fix_L_24   = self.ampl.getParameter('fix_L').to_dict()
    
            q1_24 = self.q1  # flow on new parallel pipe
            q2_24 = self.q2  # flow on existing pipe
            l1_24_raw = self.l1
    
            # alpha for regular arcs (variable-diameter, l[i,j,k])
            alpha_var = _alpha_cache_build(
                [(i, j) for (i, j) in self.arcs]
            )
 
            # alpha for new parallel pipe on fixarcs (l1[i,j,k])
            alpha_fix_new = {}
            for (i, j) in fixarcs_24:
                alpha_fix_new[i, j] = sum(
                    OMEGA * l1_24_raw[i, j, k] / pipe_denom[k]
                    for k in self.pipes
                )
    
            # Regular (non-fixed) arcs: con2 in the AMPL model
            for (i, j) in self.arcs:
                if (i, j) in fixarcs_24:
                    continue
                lhs   = h_values[i] - h_values[j]
                alpha = alpha_var[i, j]
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(q_values[i, j], alpha),
                             _hl_approx(q_values[i, j], epsilon[i, j], alpha))
    
            # Fixed arcs: con2_ (new parallel pipe, q1) + con2__ (existing, q2)
            for (i, j) in fixarcs_24:
                lhs = h_values[i] - h_values[j]
    
                alpha_new = alpha_var[i, j]
                # hl_exist  = (OMEGA * fix_L_24[i, j]
                #              / (fix_r_24[i, j] ** 1.852 * fixdiam_24[i, j] ** 4.87))
    
                _record_con2(f"{i},{j}_new", lhs,
                             _hl_original(q1_24[i, j], alpha_new),
                             _hl_approx(q1_24[i, j], epsilon[i, j], alpha_new))
    
                _record_con2(f"{i},{j}_exist", lhs,
                             _hl_original(q2_24[i, j], alpha_new),
                             _hl_approx(q2_24[i, j], epsilon[i, j],alpha_new))
    
        # ── Generic fallback (same logic as original 'else') ──────────────────
        else:
            alpha_map = _alpha_cache_build(q_values)
            for (i, j) in q_values:
                lhs   = h_values[i] - h_values[j]
                alpha = alpha_map[i, j]
                _record_con2(f"{i},{j}", lhs,
                             _hl_original(q_values[i, j], alpha),
                             _hl_approx(q_values[i, j], epsilon[i, j], alpha))

        # ══════════════════════════════════════════════════════════════════════
        # Constraint 3 – Pipe lengths must sum to arc length
        # ══════════════════════════════════════════════════════════════════════
        con3_gap = {}

        if self.data_number == 7:
            arcs_to_check = [(i, j) for (i, j) in self.arcs if (i, j) not in self.fixarcs]
        elif self.data_number==3:
            arcs_to_check = self.unfixed_arcs.to_list() + self.parallel_arcs.to_list()
        elif self.data_number==9:
            arcs_to_check = self.unfixed_arcs
            arc_len_sum = defaultdict(float)
            for (i,j,k), val in self.l1.items():
                arc_len_sum[i,j] += val

            for (i, j) in self.parallel_arcs:
                violation = arc_len_sum[i, j] - self.L[i, j]
                con3_gap[f"{i},{j}"] = violation
                total_absolute_constraint_violation += abs(violation)
        elif self.data_number==13:
            arcs_to_check = list(set(self.arcs) - set(self.parallel_arcs))

            arc_len_sum = defaultdict(float)
            for (i,j,k), val in self.l1.items():
                arc_len_sum[i,j] += val

            for (i, j) in self.parallel_arcs:
                violation = arc_len_sum[i, j] - self.L[i, j]
                con3_gap[f"{i},{j}"] = violation
                total_absolute_constraint_violation += abs(violation)
            
            arc_len_sum = defaultdict(float)
            for (i,j,k), val in self.l2.items():
                arc_len_sum[i,j] += val

            for (i, j) in self.parallel_arcs:
                violation = arc_len_sum[i, j] - self.L[i, j]
                con3_gap[f"{i},{j}"] = violation
                total_absolute_constraint_violation += abs(violation)

        else:
            arcs_to_check = list(self.arcs)

        # Accumulate per-arc totals in one pass over l_values keys
        arc_len_sum = defaultdict(float)
        for (i, j, k), val in l_values.items():
            if (i, j) in set(arcs_to_check):
                arc_len_sum[i, j] += val

        for (i, j) in arcs_to_check:
            violation = arc_len_sum[i, j] - self.L[i, j]
            con3_gap[f"{i},{j}"] = violation
            total_absolute_constraint_violation += abs(violation)

        # ══════════════════════════════════════════════════════════════════════
        # Constraint 4 – Individual pipe segment ≤ arc length
        # ══════════════════════════════════════════════════════════════════════
        con4_gap = {}
        for (i, j, k), val in l_values.items():
            if (i, j) not in arc_set:
                continue
            violation = max(0.0, val - self.L[i, j])
            if violation:
                con4_gap[f"{i},{j},{k}"] = violation
                total_absolute_constraint_violation += violation

        # ══════════════════════════════════════════════════════════════════════
        # Constraint 5 – Source node head equality
        # ══════════════════════════════════════════════════════════════════════
        con5_gap = {}
        for j in self.source:
            violation = h_values[j] - self.E[j]
            con5_gap[str(j)] = violation
            total_absolute_constraint_violation += abs(violation)

        # ══════════════════════════════════════════════════════════════════════
        # Constraint 6 – Minimum pressure head at demand nodes
        # ══════════════════════════════════════════════════════════════════════
        con6_gap = {}
        for j in self.nodes:
            if j in self.source:
                continue
            violation = max(0.0, self.E[j] + self.P[j] - h_values[j])
            if violation:
                con6_gap[str(j)] = violation
                total_absolute_constraint_violation += violation
    
        # ══════════════════════════════════════════════════════════════════════
        # Reporting
        # ══════════════════════════════════════════════════════════════════════
        separator = "*" * 79
    
        print(f"\nTotal absolute constraint violation: {total_absolute_constraint_violation:.2e}")
        print(f"\n{separator}")
        print("HEAD-LOSS CONSTRAINT: ORIGINAL vs. APPROXIMATION\n")
        print(f"\nSum of original head-loss violation : {con2_original_violation:.2e}")
        print(f"Sum of approx  head-loss violation  : {con2_approx_violation:.2e}")
        print(f"Sum of absolute diff (orig vs approx): {con2_absolute_constraint_violation:.2e}")
        print(f"Sum of relative diff (orig vs approx): {con2_relative_constraint_violation:.2e}")
        print(f"\n{separator}\n")


    def constraint_violations2(self, q_values, h_values, l_values, epsilon, solver):
        """
        Compute and report all constraint violations for the current solution.
    
        Constraints checked:
            con1: Flow balance at non-source nodes
            con2: Head-loss (original vs. approximated)
            con3: Pipe length summation per arc
            con4: Individual pipe segment length upper bound
            con5: Source node head equality
            con6: Minimum pressure head at demand nodes
        """
        total_absolute_constraint_violation = 0.0
    
        # ------------------------------------------------------------------
        # Constraint 1: Flow balance at non-source nodes
        # ------------------------------------------------------------------
        con1_gap = {}
    
        if self.data_number == 5:
            q1 = self.ampl.get_variable('q1').get_values().to_dict()
            q2 = self.ampl.get_variable('q2').get_values().to_dict()
    
            for i in self.nodes:
                if i in self.source:
                    continue
                incoming = sum(q1[j, i] + q2[j, i] for j in self.nodes if (j, i) in self.arcs)
                outgoing = sum(q1[i, j] + q2[i, j] for j in self.nodes if (i, j) in self.arcs)
                violation = (incoming - outgoing) - self.D[i]
                con1_gap[str(i)] = violation
                total_absolute_constraint_violation += abs(violation)
        else:
            for i in self.nodes:
                if i in self.source:
                    continue
                incoming = sum(q_values[j, i] for j in self.nodes if (j, i) in self.arcs)
                outgoing = sum(q_values[i, j] for j in self.nodes if (i, j) in self.arcs)
                violation = (incoming - outgoing) - self.D[i]
                con1_gap[str(i)] = violation
                total_absolute_constraint_violation += abs(violation)
    
        # ------------------------------------------------------------------
        # Constraint 2: Head-loss constraint (original vs. approximated)
        # ------------------------------------------------------------------
        con2_original_gap = {}
        con2_approx_gap = {}
        absolute_violations = {}
        relative_violations = {}
    
        con2_original_violation = 0.0
        con2_approx_violation = 0.0
        con2_absolute_constraint_violation = 0.0
        con2_relative_constraint_violation = 0.0
    
        def _alpha_rhs(i, j):
            """Weighted head-loss resistance coefficient for arc (i, j)."""
            return sum(
                10.67 * l_values[i, j, k] / ((self.R[k] ** 1.852) * (self.d[k] ** 4.87))
                for k in self.pipes
            )
    
        def _original_rhs(q, alpha):
            return q * (abs(q) ** 0.852) * alpha
    
        def _approx_rhs(q, eps, alpha):
            return (q**3 * ((q**2 + eps**2) ** 0.426) / (q**2 + 0.426 * eps**2)) * alpha
    
        def _record_con2(key, lhs, orig_rhs, apprx_rhs):
            """Update all con2 tracking dictionaries and running totals."""

            nonlocal con2_original_violation, con2_approx_violation
            nonlocal con2_absolute_constraint_violation, con2_relative_constraint_violation
            nonlocal total_absolute_constraint_violation

            con2_original_gap[key] = lhs - orig_rhs
            con2_approx_gap[key]   = lhs - apprx_rhs
    
            abs_viol = orig_rhs - apprx_rhs
            rel_viol = abs_viol / (orig_rhs + 1e-14)
    
            absolute_violations[key] = abs_viol
            relative_violations[key] = rel_viol
    
            con2_original_violation            += abs(lhs - orig_rhs)
            con2_approx_violation              += abs(lhs - apprx_rhs)
            con2_absolute_constraint_violation += abs(abs_viol)
            con2_relative_constraint_violation += abs(rel_viol)
            total_absolute_constraint_violation += abs(lhs - apprx_rhs)
    
        if self.data_number == 5:
            self.exdiam = self.ampl.getParameter('exdiam').to_dict()
    
            for (i, j) in q1.keys():
                lhs = 2 * (h_values[i] - h_values[j])
                alpha = sum(
                    10.67 * l_values[i, j, k] / ((self.R[i, j] ** 1.852) * (self.d[k] ** 4.87))
                    for k in self.pipes
                )
                hl_coeff = 10.67 * self.L[i, j] / (self.R[i, j] ** 1.852 * self.exdiam[i, j] ** 4.87)
    
                orig_rhs = (
                    q1[i, j] * (abs(q1[i, j]) ** 0.852) * hl_coeff
                    + q2[i, j] * (abs(q2[i, j]) ** 0.852) * alpha
                )
                apprx_rhs = (
                    _approx_rhs(q1[i, j], epsilon[i, j], hl_coeff)
                    + _approx_rhs(q2[i, j], epsilon[i, j], alpha)
                )
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
        elif self.data_number == 6:
            self.fixarcs  = self.ampl.getSet('fixarcs')
            self.fix_r    = self.ampl.getParameter('fix_r').to_dict()
            self.exdiam   = self.ampl.getParameter('fixdiam').to_dict()
    
            # Variable arcs
            for (i, j) in self.arcs:
                if (i, j) in self.fixarcs:
                    continue
                lhs      = h_values[i] - h_values[j]
                alpha    = _alpha_rhs(i, j)
                orig_rhs = _original_rhs(q_values[i, j], alpha)
                apprx_rhs = _approx_rhs(q_values[i, j], epsilon[i, j], alpha)
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
            # Fixed-diameter arcs
            for (i, j) in self.fixarcs:
                lhs = h_values[i] - h_values[j]
                hl_coeff = 10.67 * self.L[i, j] / (self.fix_r[i, j] ** 1.852 * self.exdiam[i, j] ** 4.87)
                orig_rhs  = _original_rhs(q_values[i, j], hl_coeff)
                apprx_rhs = _approx_rhs(q_values[i, j], epsilon[i, j], hl_coeff)
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
        else:
            for (i, j) in q_values:
                lhs       = h_values[i] - h_values[j]
                alpha     = _alpha_rhs(i, j)
                orig_rhs  = _original_rhs(q_values[i, j], alpha)
                apprx_rhs = _approx_rhs(q_values[i, j], epsilon[i, j], alpha)
                _record_con2(f"{i},{j}", lhs, orig_rhs, apprx_rhs)
    
        # ------------------------------------------------------------------
        # Constraint 3: Pipe lengths must sum to arc length
        # ------------------------------------------------------------------
        con3_gap = {}
        arcs_to_check = (
            [(i, j) for (i, j) in self.arcs if (i, j) not in self.fixarcs]
            if self.data_number == 6
            else list(self.arcs)
        )
        for (i, j) in arcs_to_check:
            violation = sum(l_values[i, j, k] for k in self.pipes) - self.L[i, j]
            con3_gap[f"{i},{j}"] = violation
            total_absolute_constraint_violation += abs(violation)
    
        # ------------------------------------------------------------------
        # Constraint 4: Individual pipe segment ≤ arc length
        # ------------------------------------------------------------------
        con4_gap = {}
        for (i, j) in self.arcs:
            for k in self.pipes:
                violation = max(0.0, l_values[i, j, k] - self.L[i, j])
                con4_gap[f"{i},{j},{k}"] = violation
                total_absolute_constraint_violation += violation
    
        # ------------------------------------------------------------------
        # Constraint 5: Source node head equality
        # ------------------------------------------------------------------
        con5_gap = {}
        for j in self.source:
            violation = h_values[j] - self.E[j]
            con5_gap[str(j)] = violation
            total_absolute_constraint_violation += abs(violation)
    
        # ------------------------------------------------------------------
        # Constraint 6: Minimum pressure head at demand nodes
        # ------------------------------------------------------------------
        con6_gap = {}
        for j in self.nodes:
            if j in self.source:
                continue
            violation = max(0.0, self.E[j] + self.P[j] - h_values[j])
            con6_gap[str(j)] = violation
            total_absolute_constraint_violation += violation
    
        # ------------------------------------------------------------------
        # Reporting
        # ------------------------------------------------------------------
        separator = "*" * 79
    
        # print(f"\n{separator}")
        # print("CONSTRAINT VIOLATIONS\n")
        #
        # table_data = (
        #     [(k, f"{v:.8f}") for k, v in con1_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con2_approx_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con3_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con4_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con5_gap.items()]
        #     + [(k, f"{v:.8f}") for k, v in con6_gap.items()]
        # )
        # print(tabulate(table_data, headers=["Constraint ID", "Violation"], tablefmt="grid"))
        print(f"\nTotal absolute constraint violation: {total_absolute_constraint_violation:.2e}")
    
        print(f"\n{separator}")
        print("HEAD-LOSS CONSTRAINT: ORIGINAL vs. APPROXIMATION\n")
    
        table_data = []
        for key, rel_viol in relative_violations.items():
            i_str, j_str = key.split(',')
            i, j = int(i_str), int(j_str)
            table_data.append([
                key,
                q_values[i, j],
                f"{con2_original_gap[key]:.2e}",
                f"{con2_approx_gap[key]:.2e}",
                f"{absolute_violations[key]:.2e}",
                f"{rel_viol:.2e}",
            ])
    
        headers = [
            "Arc (i,j)",
            "Flow",
            "Original Violation",
            "Approx Violation",
            "Absolute Diff",
            "Relative Diff",
        ]
        print(tabulate(table_data, headers=headers, tablefmt="grid"))
        print(f"\nSum of original head-loss violation : {con2_original_violation:.2e}")
        print(f"Sum of approx  head-loss violation  : {con2_approx_violation:.2e}")
        print(f"Sum of absolute diff (orig vs approx): {con2_absolute_constraint_violation:.2e}")
        print(f"Sum of relative diff (orig vs approx): {con2_relative_constraint_violation:.2e}")
        print(f"\n{separator}\n")

    def constraint_violations2(self, q_values, h_values, l_values, epsilon, solver):
        total_absolute_constraint_violation = 0
        total_relative_constraint_violation = 0

        con1_gap = {}
        if self.data_number==5:
            q1 = self.ampl.get_variable('q1').get_values().to_dict()
            q2 = self.ampl.get_variable('q2').get_values().to_dict()
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
                # original_rhs =  q_values[i, j] * (abs(q_values[i, j])) ** 0.852 * alpha_rhs
                # Approximated constraint value
                approx_rhs = (q1[i, j]**3 * ((q1[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q1[i,j]**2 + 0.426*epsilon[i,j]**2)) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j]**3 * ((q2[i, j]**2 + epsilon[i,j]**2) ** 0.426)/(q2[i,j]**2 + 0.426*epsilon[i,j]**2)) * alpha_rhs
                # approx_rhs = (q1[i, j] * ((q1[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * 10.67 * self.L[i,j]/(self.R[i,j]**1.852 * self.exdiam[i,j]**4.87) + (q2[i, j] * ((q2[i, j]**2 + epsilon[i,j]**2) ** 0.426)) * alpha_rhs

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
                # approx_rhs = self.a[i,j]*q_values[i, j]*abs(q_values[i,j]) + self.b[i,j]*q_values[i,j]

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

        print("*******************************************************************************\n")
        print("Constraints violation:\n")

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
        print("\nSum of constraints violation:", total_absolute_constraint_violation)

        #print("*******************************************************************************\n")
        #table_data = []
        #for constraint, vio in con2_original_gap.items():
        #       table_data.append([constraint, f"{con2_original_gap[constraint]:.8f}",  f"{con2_approx_gap[constraint]:.8f}"])

        print("*******************************************************************************\n")
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

        print("*******************************************************************************\n")
        print("Absolute and relative violations between original and approximation constraint 2:\n")
        headers = ["Constraint ID", "flow value", "Original Con Violation", "Approx Con Violation", "Absolute Violation", "Relative Violation"]
        print(tabulate(table_data, headers=headers, tablefmt="grid"))
        print("\nSum of violation of original headloss constraint:", con2_original_violation) 
        print("Sum of violation of approx headloss constraint:", con2_approx_violation)
        print("\nCon2 sum of absolute violation between original function and approximate function:", con2_absolute_constraint_violation)
        print("Con2 sum of relative violation between original function and approximate function:", con2_relative_constraint_violation)

        # Print total violations
        #print("\nTotal absolute constraint violation:", total_absolute_constraint_violation)
        #print("Total relative constraint violation:", total_relative_constraint_violation)

        print("*******************************************************************************\n")


    def kkt_condition_violation(self, l, q, h, eps, lam, x, y, u, v, w, tol=1e-6):
        #  Stationary Condition Violation
        violations = {}
        # dL/dl = 0
        for (i,j) in self.arcs:
            for k in self.pipes:
                lhs = self.C[k] - 10.67*q[i,j]*abs(q[i,j])**0.852 * x[i,j] / (self.R[k]**1.852 * self.d[k]**4.87) + y[i,j] - w[i,j,k]
                # lhs = self.C[k] - (10.67 * q[i,j]**3 * (q[i,j]**2 + eps[i,j]**2)**0.426 / (q[i,j]**2 + 0.426 * eps[i,j]**2)) * x[i,j] / (self.R[k]**1.852 * self.d[k]**4.87) + y[i,j] - w[i,j,k]
                if abs(lhs) > tol:
                    violations.setdefault("con1", []).append(((i,j,k), lhs))
        # dL/dq = 0
        for (i,j) in self.arcs:
            lhs = lam[j] - lam[i] - 1.852 * x[i,j] * abs(q[i,j])**0.852 * sum(10.67*l[i,j,k]/(self.R[k]**1.852 * self.d[k]**4.87) for k in self.pipes)
            if abs(lhs) > tol:
                violations.setdefault("con2", []).append(((i,j), lhs))
        # dL/dh_s = 0
        for s in self.source:
            lhs = -sum(x[s,j] for j in self.nodes if (s,j) in self.arcs) + u[s]
            if abs(lhs) > tol:
                violations.setdefault("con4", []).append((s, lhs))

        # dL/dh_j = 0
        for j in self.nodes:
            if j not in self.source:
                lhs = -sum(x[i,j] for i in self.nodes if (i,j) in self.arcs) + sum(x[j,i] for i in self.nodes if (j,i) in self.arcs) - v[j]
                if abs(lhs) > tol:
                    violations.setdefault("con3", []).append((j, lhs))

        # Dual Feasibility violation
        # con5
        for j in self.nodes: 
            if j not in self.source:
                if v[j] <= -tol:  # should be >= 0
                    violations.setdefault("con5", []).append((j, v[j]))

        # Complementary Slackness Violations
        # con6
        for j in self.nodes:
            if j not in self.source:
                lhs = (-h[j] + self.E[j] + self.P[j]) * v[j]
                if abs(lhs) > tol:
                    violations.setdefault("con6", []).append((j, lhs))

        # Primal Feasibility constraints violation
        for (i,j) in self.arcs:
            lhs = h[i] - h[j] - q[i,j] * abs(q[i,j])**0.852 * sum(10.67*l[i,j,k]/(self.R[k]**1.852 * self.d[k]**4.87) for k in self.pipes)
            if abs(lhs) > tol:
                violations.setdefault("con2", []).append(((i,j), lhs))

        return violations



    def solve_original_model_without_init(self):
        print(f"\n-------------------------------- Solving with {self.solver_name} --------------------------")
        self.ampl.option['solver'] = sys.argv[1]
        self.ampl.option["ipopt_options"] = "outlev = 5 bound_relax_factor=0  bound_push = 0.1 bound_frac = 0.1 warm_start_init_point = no halt_on_ampl_error = yes"

        self.ampl.option["bonmin_options"] = "bonmin.bb_log_level 0 bonmin.nlp_log_level 0 bonmin.num_resolve_at_root = 10 expect_infeasible_problem = no bonmin.time_limit = 600 option_file_name = ipopt.opt print_user_options = yes bonmin.nlp_log_at_root = 5 option_file_name = ipopt.opt linear_solver = ma57 mu_strategy adaptive"
        # self.ampl.option["bonmin_options"] = "bonmin.bb_log_level = 5 outlev = 4 option_file_name = ipopt.opt"

        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 2 nodemethod = 2 concurrentmethod = 3 nonconvex = 2  warmstart = 1 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 600 warmstart = 0 method = 1  mipgapabs = 1e-6 mipgap = 1e-9 barconvtol = 1e-9 sol:chk:feastol = 1e-5 sol:chk:feastolrel = 1e-9 NumericFocus = 1 tech:optionfile = gurobiOpt.prm" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        self.ampl.option["gurobi_options"] = "outlev 1 outlev_mp = 1 presolve 1 aggregate = 1 timelimit 600 alg:numericfocus = 1 obbt = 0 pre:scale = 1 method = 2 nodemethod = 1   nonconvex = 2 mipfocus = 1 nlpheur = 1 varbranch 0  mipgapabs = 1e-6 mipgap = 1e-6 alg:feastol = 1e-6 pre:feastol = 1e-6 pre:feastolrel = 1e-9 chk:feastol = 1e-6 chk:feastolrel = 1e-9  mip:heurfrac = 0.05" 
        # self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 4 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 feastol = 1e-6" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 0 timelimit 3600 NumericFocus = 1" # iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 basis = 1 premiqcpform = 2 preqlin = 2"# intfeastol = 1e-5 feastol = 1e-6 chk:epsrel = 1e-6 checkinfeas chk:inttol = 1e-5 scale = 3 aggregate = 1 intfocus = 1  BarHomogeneous = 1  startnodelimit = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10

        # self.ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 version objbound wantsol = 2 iisfind = 4 threads = 8 epsr = 1e-9" # lsolver = conopt
        self.ampl.option["baron_options"]= "optfile = optfile version objbound wantsol = 2 outlev = 2 barstats prloc=1 trace=1 " # lsolver = conopt
        #self.ampl.option["baron_options"]= "optfile = optfile" # lsolver = conopt
        self.ampl.option["scip_options"] = "param:read = scip.set" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        # self.ampl.option["scip_options"] = "outlev 1 timelimit 3600 heu:settings = 0 method = p lim:absgap=1e-6 lim:gap = 1e-9 chk:feastol = 1e-6 chk:feastolrel=1e-9 param:read = scip.set" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        # self.ampl.option["scip_options"] = "outlev  1 "
        self.ampl.option["knitro_options"] = "maxtime_real = 600 outlev = 4 opttol_abs=1e-6 opttol = 1e-6 feastol_abs = 1.0e-6 feastol = 1.0e-9  ms_enable = 1 ms_maxsolves = 10"
        #self.ampl.option["conopt_options"]= "outlev = 4"
        # self.ampl.option["presolve"] = "1"
        self.ampl.option["presolve_eps"] = "8.53e-15" 

        print(f"{sys.argv[1]} solver outputs:")
        print("---------------------\n")
        self.ampl.solve()
        # self.ampl.eval('write gsol "warmstart.sol";')

        total_cost = self.ampl.getObjective("total_cost").value()
        # print("total_cost:", total_cost, "\n")

        print("*******************************************************************************\n")
        self.q = self.ampl.get_variable('q').get_values().to_dict()
        self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.l = self.ampl.get_variable('l').get_values().to_dict()
        self.eps = self.ampl.getParameter('eps').get_values().to_dict()

        # print(self.l)
        # self.ampl.eval("display l;")
        # self.ampl.eval("display {(i,j) in arcs, k in pipes: l[i,j,k]>1e-6}: l[i,j,k];")

        if self.data_number in [3,4]:
            self.q1 = self.ampl.get_variable('q1').get_values().to_dict()
            self.q2 = self.ampl.get_variable('q2').get_values().to_dict()
        elif self.data_number in [9,24]:
            self.q1 = self.ampl.get_variable('q1').get_values().to_dict()
            self.q2 = self.ampl.get_variable('q2').get_values().to_dict()
            self.l1 = self.ampl.get_variable('l1').get_values().to_dict()

        if self.data_number==13:
            self.q1 = self.ampl.get_variable('q1').get_values().to_dict()
            self.q2 = self.ampl.get_variable('q2').get_values().to_dict()
            self.l1 = self.ampl.get_variable('l1').get_values().to_dict()
            self.l2 = self.ampl.get_variable('l2').get_values().to_dict()

        self.constraint_violations(self.q, self.h, self.l, self.eps)

        #ampl.eval("display con1.body;")
        #ampl.eval("display con2.body;")
        #ampl.eval("display con3.body;")
        #ampl.eval("display con4.body;")
        #ampl.eval("display con5.body;")
        #ampl.eval("display con6_.body;")
        #ampl.eval("display con7.body;")
        #ampl.eval("display con8.body;")

        solve_time = self.ampl.get_value('_solve_elapsed_time')
        total_cost = self.ampl.getObjective("total_cost").value()

        print(f"Total cost using {sys.argv[1]}:", total_cost)
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")

        # Extract solutions
        self.q_init = self.ampl.get_variable('q').get_values().to_dict()
        self.h_init = self.ampl.get_variable('h').get_values().to_dict()
        self.l_init = self.ampl.get_variable('l').get_values().to_dict()
        # self.eps_init = self.ampl.get_variable('eps').get_values().to_dict()
        if self.data_number ==4:
            self.q1_init = self.ampl.get_variable('q1').get_values().to_dict()
            self.q2_init = self.ampl.get_variable('q2').get_values().to_dict()
        print("*******************************************************************************\n")

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

    def second_solve(self,model_file, solver_name):
        print(f"\n-------------------------------- Solving with {self.solver_name} --------------------------")

        # self.ampl.reset()
        #self.read_model_and_data()
        self.ampl.reset()
        self.ampl.read(model_file)
        self.ampl.read_data(self.data_file)

        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')

        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
        d_min = self.ampl.getParameter('d_min').getValues().to_list()
        r_min = min(self.ampl.getParameter('R').getValues().to_list())
        #epsilon = max(ampl.getParameter('eps').getValues().to_list())
        MaxK = {}
        for (i,j) in self.arcs:
            MaxK[i,j] = 10.67*self.L[i,j]/((r_min[0]**1.852) * (d_min[0]**4.87)) 
            # print(MaxK[i,j])
        eps = {}
        # MaxK = self.ampl.getParameter('MaxK').getValues().to_dict()
        for (i,j) in self.arcs:
            val = ((10**(-1))/(0.36061*MaxK[i,j]))**(1/1.852)
            # val = ((10**(-1))/(0.07508 *MaxK[i,j]))**(1/1.852)
            order = math.floor(math.log10(abs(val)))   # find exponent
            coeff = round(val / (10 ** order), 1)
            eps[i, j] = coeff * (10 ** order)
            #eps[i,j] = 0.0048 * (10**(-2))

        # print("epsilon:", eps)
        # self.ampl.param["eps"] = eps

        # q1 = [0.005053, 0.003739, 0.066567, 0.066567, 0.066567, 0.68658, 0.0011742, 0.0004075, 0.00040753, 0.0004075, 0.0370986, 0.005988, 0.004891, 0.029482]
        # q1 = [0.105772, 0.078264, 1.393482, 1.393482, 1.393482, 14.37264, 0.02458, 0.0085311, 0.0085311, 0.0085311, 0.776606, 0.125357, 0.102378, 0.61715]
        # q1 = [0.271782, 0.201101, 3.58056, 3.58056, 3.58056, 36.93059, 0.063159, 0.02192076, 0.0219208, 0.0219208, 1.99549, 0.3221, 0.263062, 1.58578]
        # q1 = [0.36625, 0.271001, 4.82512, 4.82512, 4.82512, 49.76716, 0.085113, 0.02954012, 0.0295401, 0.0295401, 2.6891, 0.43406, 0.354499, 2.13698]
        # self.ampl.param["q1"] = q1[data_number]
        # print(q1[data_number])

        # Set initial values
        #q_var = self.ampl.get_variable('q')
        #h_var = self.ampl.get_variable('h')
        #l_var = self.ampl.get_variable('l')
        #eps_var = self.ampl.get_variable('eps')

        #for idx in self.q_init:
        #   q_var[idx].set_value(self.q_init[idx])
        #for idx in self.h_init:
        #   h_var[idx].set_value(self.h_init[idx])
        #for idx in self.l_init:
        #   l_var[idx].set_value(self.l_init[idx])
        #for idx in self.eps_init:
        #   eps_var[idx].set_value(self.eps_init[idx])

        # for (i, j, k), val in self.l_init.items():
        #    self.ampl.eval(f'let l[{i},{j},{k}] := {val};')
        # for (i, j), val in self.q_init.items():
        #    self.ampl.eval(f'let q[{i},{j}] := {val};')
        #    if self.data_number ==5:
        #        self.ampl.eval(f'let q1[{i},{j}] := {self.q1_init[i,j]};')
        #        self.ampl.eval(f'let q2[{i},{j}] := {self.q2_init[i,j]};')
        # for i, val in self.h_init.items():
        #    self.ampl.eval(f'let h[{i}] := {val};')

        # for (i, j), val in self.eps_init.items():
        #     self.ampl.eval(f'let eps[{i},{j}] := {val};')

        #for (i, j, k), val in self.l_init.items():
        #    self.ampl.var["l"][i,j,k].fix(val)
        #for (i, j), val in self.q_init.items():
        #    self.ampl.var["q"][i,j].fix(val)
        #    if self.data_number ==5:
        #        self.ampl.var["q1"][i,j].fix(self.q1[i,j])
        #        self.ampl.var["q2"][i,j].fix(self.q2[i,j])
        #for i, val in self.h_init.items():
        #    self.ampl.var["h"][i].fix(val)

        # Change solver and solve
        self.ampl.option['solver'] = solver_name

        self.ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"

        # self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = no bound_relax_factor=0 tol = 1e-9 bound_push = 0.1 bound_frac = 0.1 warm_start_init_point = no halt_on_ampl_error = yes "

        #ampl.set_option("ipopt_options", "outlev = 0 expect_infeasible_problem = yes bound_push = 0.001 bound_frac = 0.001 nlp_scaling_method = gradient-based  warm_start_init_point = yes halt_on_ampl_error = yes warm_start_bound_push=1e-9 warm_start_mult_bound_push=1e-9")   #max_iter = 1000
        self.ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 2 warm_start_init_point = no bonmin.num_resolve_at_root = 10 tol = 1e-9 expect_infeasible_problem = no bound_relax_factor = 0 bound_push = 0.1 bound_frac = 0.1 bonmin.time_limit = 600 option_file_name = ipopt.opt print_user_options = yes outlev = 1 bonmin.nlp_log_at_root = 5 linear_solver = ma57"
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 2 nodemethod = 2 concurrentmethod = 3 nonconvex = 2  warmstart = 1 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 600 warmstart = 0 method = 1  mipgapabs = 1e-6 mipgap = 1e-9 barconvtol = 1e-9 sol:chk:feastol = 1e-5 sol:chk:feastolrel = 1e-9 NumericFocus = 1 tech:optionfile = gurobiOpt.prm" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        self.ampl.option["gurobi_options"] = "outlev 1 outlev_mp = 1 presolve 1 aggregate = 1 timelimit 600 alg:numericfocus = 1 obbt = 0 pre:scale = 1 method = 2 nodemethod = 1   nonconvex = 2 mipfocus = 1 nlpheur = 1 varbranch 0  mipgapabs = 1e-6 mipgap = 1e-9 alg:feastol = 1e-6 pre:feastol = 1e-6 pre:feastolrel = 1e-9 chk:feastol = 1e-6 chk:feastolrel = 1e-9  mip:heurfrac = 0.05" 
        # self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 4 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 feastol = 1e-6" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 0 timelimit 3600 NumericFocus = 1" # iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 basis = 1 premiqcpform = 2 preqlin = 2"# intfeastol = 1e-5 feastol = 1e-6 chk:epsrel = 1e-6 checkinfeas chk:inttol = 1e-5 scale = 3 aggregate = 1 intfocus = 1  BarHomogeneous = 1  startnodelimit = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10

        # self.ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 version objbound wantsol = 2 iisfind = 4 threads = 8 epsr = 1e-9" # lsolver = conopt
        self.ampl.option["baron_options"]= "optfile = optfile version objbound wantsol = 2 outlev = 2 barstats" # lsolver = conopt
        #self.ampl.option["baron_options"]= "optfile = optfile" # lsolver = conopt
        self.ampl.option["scip_options"] = "param:read = scip.set" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        # self.ampl.option["scip_options"] = "outlev 1 timelimit 3600 heu:settings = 0 method = p lim:absgap=1e-6 lim:gap = 1e-9 chk:feastol = 1e-6 chk:feastolrel=1e-9 param:read = scip.set" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        # self.ampl.option["scip_options"] = "outlev  1 "
        self.ampl.option["knitro_options"] = "maxtime_real = 600 outlev = 4 opttol_abs=1e-6 opttol = 1e-9 feastol_abs = 1.0e-6 feastol = 1.0e-9  ms_enable = 1 ms_maxsolves = 10"
        #self.ampl.option["conopt_options"]= "outlev = 4"
        self.ampl.option["presolve"] = "1"
        self.ampl.option["presolve_eps"] = "8.53e-15" 
        # self.ampl.option["hsllib"]= "/usr/local/lib/libma57.so"
        self.ampl.option["hsllib"]= "/home/nitishdumoliya/local/lib/x86_64-linux-gnu/libcoinhsl.so"
        # self.ampl.option["presolve_eps"] = "3.75e-14"
        #print(f"{self.solver_name} solver outputs:\n")

        #min_demand = self.ampl.getParameter('D_min').getValues().to_list()[0]
        #max_demand = self.ampl.getParameter('D_max').getValues().to_list()[0]
        #max_flow = self.ampl.getParameter('Q_max').getValues().to_list()[0]

        #print("min_demand:", min_demand)
        #print("max_demand:", max_demand)
        #print("max_flow:", max_flow)
        #d_max = self.ampl.getParameter('d_max').getValues().to_list()[0]
        #d_min = self.ampl.getParameter('d_min').getValues().to_list()[0]
        #max_L = max(self.L[i,j] for (i,j) in self.arcs)
        #R_min = min(self.R[k] for k in self.pipes)
        #MaxK = 10.67*max_L/((R_min**1.852) * (d_min**4.87))

        #epsilon = ((10**(-6))/(0.07508*MaxK))**(1/0.926)
        #epsilon = (10**(-6)/(0.04001571*MaxK))**(1/1.852)

        #eps = self.ampl.getParameter('eps').get_values().to_dict()
        #print("eps:",eps) 

        #epsilon = 1e-3
        #epsilon = self.compute_adaptive_eps(min_demand/1000)

        #print("eps:", epsilon,"\n")


        #eps = ampl.getParameter('eps').to_list()
        #for (i,j) in eps.items():
        #eps[i,j].setValue(epsilon)
        #self.ampl.eval(f"param MaxK{{arcs}};")
        #self.ampl.eval(f"subject to MaxK{{(i,j) in arcs}}: MaxK[i,j] = 10.67 * L[i,j]/({R_min**1.852}*d_min^4.87);")
        #self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = ((10^(-6))/(0.07508*MaxK[i,j]))^(1/0.0926);")
        #self.ampl.eval(f"subject to eps_selection{{(i,j) in arcs}}: eps[i,j] = 1e-6 + q[i,j]^2;")

        #self.ampl.eval(f"subject to flow_1_2: q[1,2] >= 0;")
        #self.ampl.eval(f"subject to flow_2_3: q[2,3] >= 0;")
        #self.ampl.eval(f"subject to flow_2_4: q[2,4] >= 0;")
        #self.ampl.eval(f"subject to flow_3_5: q[3,5] >= 0;")
        #self.ampl.eval(f"subject to flow_4_5: q[4,5] >= 0;")
        #self.ampl.eval(f"subject to flow_4_6: q[4,6] >= 0;")
        #self.ampl.eval(f"subject to flow_6_7: q[6,7] <= 0;")
        #self.ampl.eval(f"subject to flow_7_5: q[5,7] <= 0;")
        #self.ampl.eval(f"subject to flow_7_5: h[7] - h[5] <= 0;")
        #self.ampl.eval(f"subject to flow_7: abs(q[6,7]) - abs(q[5,7]) = D[7];")
        #self.ampl.eval(f"subject to flow_5: abs(q[5,7]) + abs(q[3,5]) + abs(q[4,5]) = D[5];")
        # Define MaxK as a computed parameter
        #self.ampl.eval(f"""param MaxK{{(i,j) in arcs}} := 10.67 * L[i,j] / ({R_min}^1.852 * {d_min}^4.87);""")

        # Now use MaxK in your constraint
        #self.ampl.eval("""subject to Epsilon_Selection{(i,j) in arcs}:eps[i,j] = (1e-6 / (0.07508 * MaxK[i,j]))^(1 / 0.0926);""")
        #self.ampl.eval("display {i in 1.._ncons} (_conname[i]);")
        #self.ampl.eval("expand ;")
        # self.ampl.eval(f"solution {warmstart.sol};")
        self.ampl.solve()

        # Check constraint violations
        self.q = self.ampl.get_variable('q').get_values().to_dict()
        self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.l = self.ampl.get_variable('l').get_values().to_dict()
        self.eps = self.ampl.getParameter('eps').get_values().to_dict()
        # self.eps = self.ampl.get_variable('eps').get_values().to_dict()
        # self.a = self.ampl.get_variable('a').get_values().to_dict()
        # self.b = self.ampl.get_variable('b').get_values().to_dict()

        # self.ampl.eval("display eps;")
        # self.ampl.eval("display l;")
        # self.ampl.eval("display q;")
        # self.ampl.eval("display h;")
        # self.ampl.eval("display q1;")
        # self.ampl.eval("display q2;")
        # self.ampl.eval("display eps;")
        # print("eps: ", self.eps[next(iter(self.eps))])
        # for (i,j) in self.arcs:
        #     if np.abs(self.q[i,j]) <=1e-3:
        #         print(f"q[{i},{j}]:",self.q[i,j])

        # current_duals = {}
        # v = {}
        # for con_name, val in self.ampl.get_constraints():
        #     dual_values = val.get_values()
        #     current_duals[con_name] = dual_values
        #     if con_name == "con1":
        #         lam = dual_values.to_dict()
        #         # print(lam)
        #     if con_name == "con2":
        #         x = dual_values.to_dict()
        #         print(x)
        #     if con_name == "con3":
        #         y = dual_values.to_dict()
        #         print(y)
        #     if con_name == "con5":
        #         w = dual_values.to_dict()
        #     if con_name == "con6":
        #         u = dual_values.to_dict() 
        #     if con_name == "con7":
        #         v = dual_values.to_dict()
        #         # for j in val.keys():
        #         #     v[j] = val[j]
        #         # print(v)

        # tol = 1e-6
        # self.ampl.eval("display con3;")
        # self.ampl.eval("display con3.dual;")

        def dualmodelsolve(file):
            ampl = AMPL()
            ampl.reset()
            # ampl.read("lagDualWdnModel.mod")
            ampl.read(file)
            ampl.read_data(self.data_file)

            for (i, j, k), value in self.l.items():
                ampl.param['l'][i, j, k] = value
            for (i, j), value in self.q.items():
                ampl.param['q'][i, j] = value
            for i, value in self.h.items():
                ampl.param['h'][i] = value
            for (i,j), value in self.eps.items():
                ampl.param['eps'][i,j] = value

            ampl.option['solver'] = "cplex"
            ampl.option["presolve_eps"] = "7.583e-13"
            ampl.solve()
            lam = ampl.get_variable('lam').get_values().to_dict()
            x = ampl.get_variable('x').get_values().to_dict()
            y = ampl.get_variable('y').get_values().to_dict()
            u = ampl.get_variable('u').get_values().to_dict()
            v = ampl.get_variable('v').get_values().to_dict()
            w = ampl.get_variable('w').get_values().to_dict()
            # s = ampl.get_variable('s').get_values().to_dict()

            # print("lam: ", lam)
            # print("x: ", x)
            # print("y: ", y)
            # print("u: ", u)
            # print("v: ", v)
            # print("w: ", w)
            # print("s: ", s)
            return lam, x, y, u, v, w
        # lam, x, y, u, v, w = dualmodelsolve("approximatekkt.mod")
        # dualmodelsolve("originalkkt.mod")

        # self.constraint_violations(self.q, self.h, self.l, self.eps, self.solver_name)
        # violations = self.kkt_condition_violation(self.l, self.q, self.h, self.eps, lam, x, y, u, v, w)
        # print("violations: ", violations)
        solve_time = self.ampl.get_value('_solve_elapsed_time')
        self.total_cost = self.ampl.getObjective("total_cost").value()

        print(f"Total cost using {self.solver_name}:", self.total_cost)
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")

        # Extract solutions
        #q_sol = self.ampl.get_variable('q').get_values().to_dict()
        #h_sol = self.ampl.get_variable('h').get_values().to_dict()
        #l_sol = self.ampl.get_variable('l').get_values().to_dict()

        #self.ampl.eval("display q;")
        #self.ampl.eval("display q2;")
        #self.ampl.eval("display l;")
        #self.ampl.eval("display {(i,j) in arcs}: h[i] - h[j];")
        #self.ampl.eval("display {(i,j) in arcs, k in pipes: l[i,j,k]>1e-6}: l[i,j,k];")
        #self.ampl.eval("display {(i,j) in arcs} h[i] - h[j] - q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));")

        #self.ampl.eval("display {(i,j) in arcs} h[i] - h[j]  - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));")

        #self.ampl.eval("display {(i,j) in arcs} sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));")

        #self.ampl.eval("display sum {(i,j) in arcs} (q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)));")

        #self.ampl.eval("display sum {(i,j) in arcs} (q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) - (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)))/(q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)));")

        print("*******************************************************************************\n")


    def solve_content_model(self):
        self.ampl.reset()
        self.ampl.read("wdn_content_model.mod")
        self.ampl.read_data(self.data_file)

        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')

        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()

        for cor, val in self.l.items():
            #ampl.eval(f"s.t. l_val{cor[0]}_{cor[1]}_{cor[2]}: l[{cor[0]},{cor[1]},{cor[2]}] = {val};")
            self.ampl.param['l'][cor[0], cor[1], cor[2]] = val


        # Set initial values
        #q_var = self.ampl.get_variable('q')
        #h_var = self.ampl.get_variable('h')
        #l_var = self.ampl.get_variable('l')

        #for idx in self.q_init:
        #    q_var[idx].set_value(self.q_init[idx])
        #for idx in self.h_init:
        #    h_var[idx].set_value(self.h_init[idx])
        #for idx in self.l_init:
        #    l_var[idx].set_value(self.l_init[idx])

        # Change solver and solve
        self.ampl.option['solver'] = self.solver_name
        # self.ampl.option["hsllib"]= "/usr/local/lib/libma57.so"        

        self.ampl.option["mmultistart_options"] = "--presolve 1 --log_level 3 --eval_within_bnds 1 --nlp_engine IPOPT"
        self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = no halt_on_ampl_error = yes hsllib = /usr/local/lib/libma57.so"
        self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 warmstart = 1 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0 mipgap = 1e-9 NumericFocus = 1" 
        self.ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 barstats version objbound" # lsolver = conopt
        self.ampl.option["scip_options"] = "outlev  1 timelimit 3600 lim:gap = 1e-9 chk:feastol = 1e-5 chk:feastolrel=0 " #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        self.ampl.option["knitro_options"]= "maxtime_real = 3600 outlev = 4 threads=8 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10"
        #self.ampl.option["conopt_options"]= "outlev = 4"
        self.ampl.option["presolve"] = "1"
        self.ampl.option["presolve_eps"] = "8.53e-15"
        self.ampl.solve()
        #self.q = self.ampl.get_variable('q').get_values().to_dict()
        #self.h = self.ampl.get_variable('h').get_values().to_dict()
        #self.l = self.ampl.get_variable('l').get_values().to_dict()
        #self.eps = self.ampl.getParameter('eps').get_values().to_dict()
        #self.ampl.eval("display eps;")
        #self.ampl.eval("display q;")
        #self.ampl.eval("display q1;")
        #self.ampl.eval("display q2;")
        #self.ampl.eval("display eps;")

        #self.constraint_violations(self.q, self.h, self.l, self.eps, self.solver_name)

        solve_time = self.ampl.get_value('_solve_elapsed_time')
        self.total_cost = self.ampl.getObjective("total_cost").value()

        print(f"Total cost using {self.solver_name}:", self.total_cost)
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")

        self.ampl.eval("display sum{(i,j) in arcs, k in pipes} l[i,j,k]*C[k];")
        # Extract solutions
        #q_sol = self.ampl.get_variable('q').get_values().to_dict()
        #h_sol = self.ampl.get_variable('h').get_values().to_dict()
        #l_sol = self.ampl.get_variable('l').get_values().to_dict()

        self.ampl.eval("display q;")
        self.ampl.eval("display con1.dual;")
        self.ampl.eval("display E[1] - sum{k in pipes} (10.67*l[1,2,k]*(q[1,2]^1.852)/(R[k]^1.852 * d[k]^4.87)) ;")
        #self.ampl.eval("display q2;")
        #self.ampl.eval("display l;")
        #self.ampl.eval("display {(i,j) in arcs}: h[i] - h[j];")
        #self.ampl.eval("display {(i,j) in arcs, k in pipes: l[i,j,k]>1e-6}: l[i,j,k];")

        # Extract duals (node potentials)
        con1 = self.ampl.getConstraint("con1")
        duals = {}

        print("\nOriginal duals (up to constant):")
        for n in self.nodes:
            duals[n] = con1[n].dual()
            print(f"  {n}: {duals[n]:.6g}")

        source = list(self.source)[0]
        H_source = self.E[source]

        lambda_source = duals[source]
        print(lambda_source)

        shift = H_source - lambda_source

        print("\nShifted duals (physical heads):")
        for n in self.nodes:
            head = duals[n] + shift
            print(f"  {n}: {head:.6g}")

        print(f"\nAt source node '{source}': shifted dual = {duals[source] + shift} (should equal H_source = {H_source})")

    def solve_reduced_model1(self):
        import os, subprocess, tempfile, time
    
        print(f"\n-------------------------------- Solving Reduced Model using {self.solver_name} --------------------------")
    
        if self.data_number == 5:
            self.model_file = "newyork_epigraph_model.mod"
        elif self.data_number == 6:
            self.model_file = "blacksburg_epigraph_model.mod"
        else:
            self.model_file = "exact_reduced_wdn.mod"
        print("WDN Model Name:", self.model_file)
    
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)
    
        self.ampl.option['solver'] = self.solver_name
        self.ampl.option["presolve"] = "0"
        self.ampl.option["presolve_eps"] = "8.53e-15"
    
        self.ampl.option["ipopt_options"] = "outlev=5 expect_infeasible_problem=no bound_relax_factor=0 bound_push=0.01 bound_frac=0.01 warm_start_init_point=no halt_on_ampl_error=yes"
        self.ampl.option["gurobi_options"] = "outlev=1 outlev_mp=1 presolve=1 aggregate=1 timelimit=600 alg:numericfocus=1 obbt=0 pre:scale=1 method=2 nodemethod=1 nonconvex=2 mipfocus=1 nlpheur=1 varbranch=0 mipgapabs=1e-6 mipgap=1e-6 alg:feastol=1e-6 pre:feastol=1e-6 pre:feastolrel=1e-9 chk:feastol=1e-6 chk:feastolrel=1e-9 mip:heurfrac=0.05"
        self.ampl.option["baron_options"] = "optfile=optfile version objbound wantsol=2 outlev=2 barstats"
        self.ampl.option["scip_options"] = "param:read=scip.set"
        self.ampl.option["knitro_options"] = "maxtime_real=600 outlev=4 opttol_abs=1e-6 opttol=1e-6 feastol_abs=1e-6 feastol=1e-9 ms_enable=1 ms_maxsolves=10"
 
        self.ampl.option["bonmin_options"] = "bonmin.bb_log_level 0 bonmin.nlp_log_level 0 bonmin.num_resolve_at_root = 10 expect_infeasible_problem = no bonmin.time_limit = 600 option_file_name = ipopt.opt print_user_options = yes bonmin.nlp_log_at_root = 5 option_file_name = ipopt.opt linear_solver = ma57 mu_strategy adaptive"
        # ── BONMIN: subprocess path ──────────────────────────────────────────────
        if self.solver_name == 'bonmin':
            tmpdir  = tempfile.mkdtemp()
            stub    = os.path.join(tmpdir, 'bonmin_prob')
            nl_file  = stub + '.nl'
            sol_file = stub + '.sol'
            opt_file = os.path.join(tmpdir, 'bonmin.opt')
    
            # write .nl file
            self.ampl.eval(f"option presolve 0; write g{stub};")
            if not os.path.exists(nl_file):
                raise RuntimeError(f"NL file not written: {nl_file}")
            print(f"NL file written ({os.path.getsize(nl_file)} bytes)")
    
            # write bonmin options file
            with open(opt_file, 'w') as f:
                f.write(
                    "bonmin.bb_log_level 1\n"
                    "bonmin.nlp_log_level 2\n"
                    "bonmin.num_resolve_at_root 10\n"
                    "bonmin.time_limit 600\n"
                    "bonmin.nlp_log_at_root 5\n"
                    "tol 1e-9\n"
                    "bound_relax_factor 0\n"
                    "bound_push 0.1\n"
                    "bound_frac 0.1\n"
                    "warm_start_init_point no\n"
                    "expect_infeasible_problem no\n"
                    "linear_solver mumps\n"
                )
    
            bonmin_exe = os.path.expanduser('~/bonmin_new2/bin/bonmin')
            env = os.environ.copy()
            env['LD_LIBRARY_PATH'] = (
                os.path.expanduser('~/bonmin_new2/lib') + ':' +
                os.path.expanduser('~/local/lib') + ':' +
                env.get('LD_LIBRARY_PATH', '')
            )
    
            print("Running bonmin...")
            t0 = time.time()
            proc = subprocess.run(
                [bonmin_exe, stub],
                capture_output=True, text=True,
                cwd=tmpdir, env=env
            )
            solve_time = time.time() - t0
            print(proc.stdout[-3000:])
            if proc.stderr:
                print("STDERR:", proc.stderr[-300:])
    
            if not os.path.exists(sol_file):
                raise RuntimeError(f"No .sol file produced. Return code: {proc.returncode}")
    
            print(f"Sol file: {os.path.getsize(sol_file)} bytes")
    
            # debug: print raw sol file
            with open(sol_file) as f:
                print("=== RAW SOL ===\n", f.read()[:1000])
    
            # read solution back into AMPL
            self.ampl.eval(f"solution '{sol_file}';")
            self.ampl.eval("display total_cost;")
    
            self.objective = self.ampl.getObjective("total_cost").value()
            print(f"Total cost using bonmin: {self.objective}")
            print(f"bonmin solve time: {solve_time:.2f} seconds")
    
        # ── ALL OTHER SOLVERS ────────────────────────────────────────────────────
        else:
            self.ampl.solve()
            solve_time = self.ampl.get_value('_solve_elapsed_time')
            self.objective = self.ampl.getObjective("total_cost").value()
            print(f"Total cost using {self.solver_name}: {self.objective}")
            print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")
    
        # ── EXTRACT SOLUTION VARIABLES (all solvers) ─────────────────────────────
        self.q   = self.ampl.get_variable('q').get_values().to_dict()
        self.h   = self.ampl.get_variable('h').get_values().to_dict()
        self.y   = self.ampl.get_variable('y').get_values().to_dict()
    
        self.nodes  = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs   = self.ampl.getSet('arcs')
        self.pipes  = self.ampl.getSet('pipes')
    
        self.L   = self.ampl.getParameter('L').to_dict()
        self.D   = self.ampl.getParameter('D').to_dict()
        self.C   = self.ampl.getParameter('C').to_dict()
        self.P   = self.ampl.getParameter('P').to_dict()
        self.R   = self.ampl.getParameter('R').to_dict()
        self.E   = self.ampl.getParameter('E').to_dict()
        self.d   = self.ampl.getParameter('d').to_dict()
        self.eps = self.ampl.getParameter('eps').to_dict()


    def solve_reduced_model(self):
        print(f"\n-------------------------------- Solving Reduced Model using {self.solver_name} --------------------------")

        if self.data_number==3:
            self.model_file = "trn_epigraph_model.mod"
        elif self.data_number == 4:
            self.model_file = "newyork_epigraph_model.mod"
        elif self.data_number == 7:
            self.model_file = "blacksburg_epigraph_model.mod"
        elif self.data_number == 9:
            self.model_file = "bakryun_epigraph_model.mod"
        elif self.data_number == 13:
            self.model_file = "farhadgerd_epigraph_model.mod"
        else:
            self.model_file = "exact_reduced_wdn.mod"
        print("WDN Model Name:", self.model_file)

        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)

        # if self.data_number==13:
        #     for (i, j), val in y_input.items():
        #         ampl.param["y"][i, j] = val
        #
        #     for (i, j), val in self.y1.items():
        #         ampl.param["y1"][i, j] = val
        #     for (i, j), val in self.y2.items():
        #         ampl.param["y2"][i, j] = val
        # else:
        #     for (i, j), val in y_input.items():
        #         ampl.param["y"][i, j] = val

        self.ampl.option['solver'] = self.solver_name

        self.ampl.option["ipopt_options"] = "outlev = 5 expect_infeasible_problem = no bound_relax_factor=0 bound_push = 0.1 bound_frac = 0.1 warm_start_init_point = no halt_on_ampl_error = yes "

        # self.ampl.option["bonmin_options"] = "bonmin.bb_log_level 0 bonmin.nlp_log_level 0 warm_start_init_point = no bonmin.num_resolve_at_root = 10 expect_infeasible_problem = yes bound_relax_factor = 0 bound_push = 0.01 bound_frac = 0.01 bonmin.time_limit = 600 print_user_options = yes outlev = 1 bonmin.nlp_log_at_root = 5 linear_solver = ma57 outlev = 0 print_level = 0"
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 2 nodemethod = 2 concurrentmethod = 3 nonconvex = 2  warmstart = 1 barconvtol = 1e-9 feastol = 1e-5 chk:epsrel = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 600 warmstart = 0 method = 1  mipgapabs = 1e-6 mipgap = 1e-9 barconvtol = 1e-9 sol:chk:feastol = 1e-5 sol:chk:feastolrel = 1e-9 NumericFocus = 1 tech:optionfile = gurobiOpt.prm" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        self.ampl.option["gurobi_options"] = "outlev 1 outlev_mp = 1 presolve 1 aggregate = 1 timelimit 600 alg:numericfocus = 1 obbt = 0 pre:scale = 1 method = 2 nodemethod = 1   nonconvex = 2 mipfocus = 1 nlpheur = 1 varbranch 0  mipgapabs = 1e-6 mipgap = 1e-6 alg:feastol = 1e-6 pre:feastol = 1e-6 pre:feastolrel = 1e-9 chk:feastol = 1e-6 chk:feastolrel = 1e-9  mip:heurfrac = 0.05" 
        # self.ampl.option["gurobi_options"] = "outlev 1 presolve 1 timelimit 3600 iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 4 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 feastol = 1e-6" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10
        #self.ampl.option["gurobi_options"] = "outlev 1 presolve 0 timelimit 3600 NumericFocus = 1" # iis = 1 iismethod = 0 iisforce = 1 NumericFocus = 1 socp = 2 method = 3 nodemethod = 1 concurrentmethod = 3 nonconvex = 2 varbranch = 0 obbt = 1 warmstart = 1 basis = 1 premiqcpform = 2 preqlin = 2"# intfeastol = 1e-5 feastol = 1e-6 chk:epsrel = 1e-6 checkinfeas chk:inttol = 1e-5 scale = 3 aggregate = 1 intfocus = 1  BarHomogeneous = 1  startnodelimit = 0" #lim:time=10 concurrentmip 8 pool_jobs 0 Threads=1 basis = 1 mipstart = 3 feastol=1e-9 mipfocus = 1 fixmodel = 1 PumpPasses = 10

        # self.ampl.option["baron_options"]= "maxtime = 3600  outlev = 2 version objbound wantsol = 2 iisfind = 4 threads = 8 epsr = 1e-9" # lsolver = conopt
        self.ampl.option["baron_options"]= "optfile = optfile version objbound wantsol = 2 outlev = 2 barstats" # lsolver = conopt
        # self.ampl.option["baron_options"]=  "lsolver = snopt"
        self.ampl.option["scip_options"] = "param:read = scip.set" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        # self.ampl.option["scip_options"] = "outlev 1 timelimit 3600 heu:settings = 0 method = p lim:absgap=1e-6 lim:gap = 1e-9 chk:feastol = 1e-6 chk:feastolrel=1e-9 param:read = scip.set" #cvt/pre/all = 0 pre:maxrounds 1 pre:settings 3 cvt:pre:all 0
        # self.ampl.option["scip_options"] = "outlev  1 "
        self.ampl.option["knitro_options"] = "maxtime_real = 600 outlev = 4 opttol_abs=1e-6 opttol = 1e-6 feastol_abs = 1.0e-6 feastol = 1.0e-9  ms_enable = 1 ms_maxsolves = 10"

        self.ampl.option["bonmin_options"] = "bonmin.bb_log_level 2 bonmin.nlp_log_level 1 bonmin.num_resolve_at_root = 20 expect_infeasible_problem = no bonmin.time_limit = 600 option_file_name = ipopt.opt print_user_options = yes bonmin.nlp_log_at_root = 2 option_file_name = ipopt.opt linear_solver = ma57 mu_strategy adaptive"
        #self.ampl.option["conopt_options"]= "outlev = 4"
        self.ampl.option["presolve"] = "0"
        self.ampl.option["presolve_eps"] = "8.53e-15" 

        # with self.suppress_output():
        self.ampl.solve()

        self.q = self.ampl.get_variable('q').get_values().to_dict()
        self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.y = self.ampl.get_variable('y').get_values().to_dict()
        
        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')

        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
        self.eps = self.ampl.getParameter('eps').to_dict()


        solve_time = self.ampl.get_value('_solve_elapsed_time')
        self.objective = self.ampl.getObjective("total_cost").value()

        print(f"Total cost using {self.solver_name}:", self.objective)
        # self.ampl.eval("display sum{(i,j) in arcs, k in pipes} l[i,j,k]*C[k];")
        print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")
        

    def solve_recover_model(self):
        print(f"\n-------------------------------- Solving Recover Model using cplex --------------------------")
        if self.data_number==3:
            ampl.read("trn_recover_model.mod")
        elif self.data_number==4:
            ampl.read("newyork_recover_model.mod")
        elif self.data_number==7:
            ampl.read("blacksburg_recover_model.mod")
        elif self.data_number==9:
            ampl.read("bakryun_recover_model.mod")
        elif self.data_number==13:
            ampl.read("farhadgerd_recover_model.mod")
        else:
            self.model_file = "recover_wdnmodel.mod"

        ampl = AMPL()
        ampl.reset()
        ampl.read(self.model_file)
        ampl.read_data(self.data_file)
       
        for (u, v), val in self.y.items():
            ampl.param["y"][u, v] = val
    
        #ampl.option['solver'] = "/home/nitishdumoliya/build/bin/ipopt"
        ampl.option['solver'] = "cplex"

        ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = no bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = no halt_on_ampl_error = yes "

        ampl.option["presolve"] = "0"
        ampl.option["presolve_eps"] = "8.53e-15" 
        #with self.suppress_output():
        ampl.solve()

        # self.q = self.ampl.get_variable('q').get_values().to_dict()
        # self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.l = ampl.get_variable('l').get_values().to_dict()

        # for (u,v) in self.arcs:
        #     for k in self.pipes:
        #         if self.l[u,v,k]>=0.0001:
        #             print(f"l[{u},{v},{k}]:", self.l[u,v,k])

        solve_time = ampl.get_value('_solve_elapsed_time')
        self.total_cost = ampl.getObjective("total_cost").value()

        self.constraint_violations(self.q, self.h, self.l, self.eps, self.solver_name)

        print(f"Total Cost:", self.total_cost)
        # self.ampl.eval("display sum{(i,j) in arcs, k in pipes} l[i,j,k]*C[k];")
        print(f"cplex solve time: {solve_time:.2f} seconds")

        # return self.total_cost, solve_time

    def solve_original_model_with_init(self):
        # print(f"\n-------------------------------- Solving Recover Model using {self.solver_name} --------------------------")

        ampl = AMPL()
        ampl.reset()
        ampl.read("wdnmodel.mod")
        ampl.read_data(self.data_file)
       
        # for (u, v), val in self.y.items():
        #     ampl.param["y"][u, v] = val
 
        for (i, j, k), val in self.l.items():
           ampl.eval(f'let l[{i},{j},{k}] := {val};')
        for (i, j), val in self.q.items():
           ampl.eval(f'let q[{i},{j}] := {val};')
           if self.data_number ==5:
               ampl.eval(f'let q1[{i},{j}] := {self.q1[i,j]};')
               ampl.eval(f'let q2[{i},{j}] := {self.q2[i,j]};')
        for i, val in self.h.items():
           ampl.eval(f'let h[{i}] := {val};')
        
        if self.data_number == 6:
            ampl.eval("subject to con3{(i,j) in arcs diff fixarcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
        else:
            ampl.eval("subject to con3{(i,j) in arcs}: sum{k in pipes} l[i,j,k] = L[i,j];")
    
        # ampl.option['solver'] = "ipopt"
        ampl.option['solver'] = "~/build/bin/ipopt"

        ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = no bound_relax_factor=0 tol = 1e-9 bound_push = 0.1 bound_frac = 0.1 warm_start_init_point = yes halt_on_ampl_error = yes "

        ampl.option["presolve"] = "1"
        ampl.option["presolve_eps"] = "8.53e-15" 
        with self.suppress_output():
            ampl.solve()

        self.q = self.ampl.get_variable('q').get_values().to_dict()
        self.h = self.ampl.get_variable('h').get_values().to_dict()
        self.l = ampl.get_variable('l').get_values().to_dict()

        # for (u,v) in self.arcs:
        #     for k in self.pipes:
        #         if self.l[u,v,k]>=0.0001:
        #             print(f"l[{u},{v},{k}]:", self.l[u,v,k])

        solve_time = ampl.get_value('_solve_elapsed_time')
        self.total_cost = ampl.getObjective("total_cost").value()

        # print(f"Total Cost:", self.total_cost)
        # self.ampl.eval("display sum{(i,j) in arcs, k in pipes} l[i,j,k]*C[k];")
        # print(f"{self.solver_name} solve time: {solve_time:.2f} seconds")

        return self.total_cost, solve_time

    def plot_value_function_for_arc(self, arc, n_points=300):
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

    def solve_reduced_recover_original(self):
        self.K = {arc: min(self.C[k]*self.R[k]**1.852 * self.d[k]**4.87/10.67 for k in self.pipes) for arc in self.arcs}

        total_iteration = 10
        iter = 1

        print(f"{'Iter':<6}{'Reduced Obj':<20}{'Reduced Time':<18}{'Recover Obj':<20}{'Recover Time':<18}{'Orig Obj':<20}{'Orig Time':<18}")
        print("-" * 95)

        while iter <= total_iteration:
            # Solve reduced model
            reduced_obj, reduced_time = self.solve_reduced_model()
            # Solve recover model
            recover_obj, recover_time = self.solve_recover_model()
            # Solve original model
            original_model_obj, original_model_time = self.solve_original_model_with_init()
            print(
                f"{iter:<6}"
                f"{reduced_obj:<20.8f}"
                f"{reduced_time:<18.4f}"
                f"{recover_obj:<20.8f}"
                f"{recover_time:<18.4f}"
                f"{original_model_obj:<20.8f}"
                f"{original_model_time:<18.4f}"
            )
            for (u, v), val in self.y.items():
                # print(f"y[{u},{v}]:",self.y[u,v])
                self.K[u, v] = sum(self.C[k]*self.l[u,v,k] for k in self.pipes)/(sum(10.67*self.l[u,v,k]/(self.R[k]**1.852 * self.d[k]**4.87) for k in self.pipes))
                # print(f"K[{u},{v}]:",self.K[u,v], "\n")

            iter += 1

    def run(self):
        self.read_model_and_data()
        self.solve_original_model_without_init()

        # self.solve_reduced_model()
        # self.solve_recover_model()

if __name__ == "__main__":
    data_list = [
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
        "d25_large"
    ]

    # data_list1 = [
    #     "d1_bessa",
    #     "d2_shamir",
    #     "d3_hanoi",
    #     "d4_double_hanoi",
    #     "d5_triple_hanoi",
    #     "d6_newyork",
    #     "d7_blacksburg",
    #     "d8_fossolo_iron",
    #     "d9_fossolo_poly_0",
    #     "d10_fossolo_poly_1",
    #     "d11_kadu",
    #     "d12_pescara",
    #     "d13_modena",
    #     "d14_balerma",
    #     "d15_goyang",
    #     "d16_bakryun",
    #     "d17_farhadgerd",
    #     "d18_ramnagar",
    #     "d20_zj_network",
    #     "d21_trn",
    #     "d22_kl",
    #     "d23_marchi_rural",
    #     "d25_small",
    #     "d26_large",
    #     "d24_cgn"
    # ]

    # Select the data number here (0 to 18)
    data_number = int(sys.argv[2]) - 1
    data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{data_list[(data_number)]}.dat"
    #data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{sys.argv[3]}"
    #data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{sys.argv[3]}"
    #data = f"/home/nitishdumoliya/waterNetwork/data/minlplib_data/{data_list1[(data_number)]}.dat"
    #print("Water Network:", f"{data_list[(data_number)]}.dat")
    print("Water Network:", f"{sys.argv[2]}")

    # print("Results of actual Hazen--Williams headloss constraint\n")
    # print("Results smooth approximation 1 of Hazen--Williams headloss constraint using absolute error\n")
    # print("Results smooth approximation 1 of Hazen--Williams headloss constraint using relative error\n")
    # print("Results smooth approximation 2 of Hazen--Williams headloss constraint using absolute error\n")
    # print("Results smooth approximation 2 of Hazen--Williams headloss constraint using relative error\n")

    print("***********************************************************************************************")

    solver = sys.argv[1] 

    solver_instance = WaterNetworkSolver(solver, data, data_number+1)
    solver_instance.run()
