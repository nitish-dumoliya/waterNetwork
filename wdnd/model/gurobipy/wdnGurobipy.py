import gurobipy as gp
# from gurobipy import GRB
from gurobipy import GRB, quicksum, Model, abs_, nlfunc
import gurobipy_pandas as gppd
gppd.set_interactive()
import numpy as np
import sys
def build_model(data):

    nodes  = data.nodes
    pipes  = data.pipes
    arcs   = data.arcs
    Source = data.Source

    E = data.E
    D = data.D
    P = data.P
    L = data.L
    d = data.d
    C = data.C
    R = data.R

    omega = data.omega

    Q_max = sum(abs(D[i]) for i in nodes if i not in Source)

    m = gp.Model("WDN")

    # =========================
    # Variables
    # =========================
    q = m.addVars(arcs, vtype=GRB.CONTINUOUS, lb=-Q_max, ub=Q_max, name="q")
    h = m.addVars(nodes, vtype=GRB.CONTINUOUS,name="h")
    l = m.addVars(arcs, pipes, vtype=GRB.CONTINUOUS, lb=0.0, name="l")
    phi = m.addVars(arcs, name="phi")
    q_abs = m.addVars(arcs,vtype=GRB.CONTINUOUS, name="q_abs")
    
    for a in arcs:
        m.addConstr(q_abs[a] == abs_(q[a]), name="absconstr")
    # =========================
    # Objective
    # =========================
    m.setObjective(
        gp.quicksum(l[i, j, k] * C[k]
                    for (i, j) in arcs
                    for k in pipes),
        GRB.MINIMIZE
    )

    # =========================
    # Flow balance
    # =========================
    for j in nodes:
        if j in Source:
            continue
        m.addConstr(
            gp.quicksum(q[i, j] for i in nodes if (i, j) in arcs)
            - gp.quicksum(q[j, i] for i in nodes if (j, i) in arcs)
            == D[j],
            name=f"flow[{j}]"
        )

    # =========================
    # Pipe length constraints
    # =========================
    for (i, j) in arcs:
        m.addConstr(
            gp.quicksum(l[i, j, k] for k in pipes) == L[i, j],
            name=f"length_sum[{i},{j}]"
        )
        for k in pipes:
            m.addConstr(
                l[i, j, k] <= L[i, j],
                name=f"length_bound[{i},{j},{k}]"
            )

    # =========================
    # Head constraints
    # =========================
    for i in Source:
        m.addConstr(h[i] == E[i], name=f"source_head[{i}]")

    for i in nodes:
        if i not in Source:
            m.addConstr(
                h[i] >= E[i] + P[i],
                name=f"min_pressure[{i}]"
            )

    # =========================
    # PWL Hazen–Williams
    # =========================
    # xs = np.linspace(-Q_max, Q_max, 50)
    # ys = [x * abs(x)**0.852 for x in xs]
    #
    # for (i, j) in arcs:
    #     m.addGenConstrPWL(q[i, j], phi[i, j], xs.tolist(), ys,
    #                       name=f"PWL_HW[{i},{j}]")
    #
    #     K_ij = gp.quicksum(
    #         omega * l[i, j, k] / (R[k]**1.852 * d[k]**4.87)
    #         for k in pipes
    #     )
    #
    #     m.addConstr(
    #         h[i] - h[j] == phi[i, j] * K_ij,
    #         name=f"headloss[{i},{j}]"
    #     )

    # =========================
    # Original Hazen–Williams
    # =========================

    for (i, j) in arcs:
        K_ij = gp.quicksum(
            omega * l[i, j, k] / (R[k]**1.852 * d[k]**4.87)
            for k in pipes
        )
        constrs[i,j] = gppd.add_constrs(m, h[i] - h[j], GRB.EQUAL, (q[i, j] * q_abs[i, j]**0.852 * K_ij).apply(nlfunc.log))
        # ORIGINAL CONSTRAINT (NOT SUPPORTED)
        # m.addConstr(
        #     h[i] - h[j]
        #     ==
        #     (q[i, j] * q_abs[i, j]**0.852 * K_ij).apply(nlfunc.log),
        #     name=f"HW_original_{i}_{j}"
        # )
    return m

if __name__ == "__main__":
    from data import d1
    # from model import build_model
    data = d1
    m = build_model(data)

    m.Params.NonConvex = 2
    m.optimize()

    print("\nOptimal cost:", m.ObjVal)

