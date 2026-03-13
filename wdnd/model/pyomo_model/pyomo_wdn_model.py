from pyomo.environ import *
from pyomo.opt import SolverFactory
import sys

# =========================
# Model
# =========================
model = AbstractModel()

# =========================
# Sets
# =========================
model.nodes  = Set()
model.pipes  = Set()
model.arcs   = Set(dimen=2)
model.Source = Set()

# =========================
# Parameters
# =========================
model.L    = Param(model.arcs)
model.E    = Param(model.nodes)
model.P    = Param(model.nodes)
model.pmax = Param(model.nodes)
model.vmax = Param(model.arcs)
model.D    = Param(model.nodes)

model.d = Param(model.pipes)
model.C = Param(model.pipes)
model.R = Param(model.pipes)

model.omega = Param(initialize=10.67)
model.p     = Param(initialize=1.852)

# =========================
# Scalar Params (FIXED)
# =========================

def Q_max_rule(m):
    return sum(abs(m.D[i]) for i in m.nodes if i not in m.Source)

model.Q_max = Param(initialize=Q_max_rule)

def d_min_rule(m):
    return min(m.d[k] for k in m.pipes)

model.d_min = Param(initialize=d_min_rule)

def R_min_rule(m):
    return min(m.R[k] for k in m.pipes)

model.R_min = Param(initialize=R_min_rule)

# =========================
# Arc-based Params
# =========================
def MaxK_rule(m, i, j):
    return m.omega * m.L[i, j] / (m.R_min**1.852 * m.d_min**4.87)

model.MaxK = Param(model.arcs, rule=MaxK_rule)

# =========================
# Expressions
# =========================
def eps_rule(m, i, j):
    return 0.0535 * (1e-3 / m.MaxK[i, j])**0.54

model.eps = Expression(model.arcs, rule=eps_rule)

# =========================
# Variables
# =========================
model.l = Var(model.arcs, model.pipes, domain=NonNegativeReals)
model.q = Var(model.arcs)
model.h = Var(model.nodes)

# =========================
# Objective
# =========================
def total_cost_rule(m):
    return sum(
        m.l[i, j, k] * m.C[k]
        for (i, j) in m.arcs
        for k in m.pipes
    )

model.total_cost = Objective(rule=total_cost_rule, sense=minimize)

# =========================
# Constraints
# =========================

# Flow balance
def flow_balance_rule(m, j):
    if j in m.Source:
        return Constraint.Skip
    return (
        sum(m.q[i, j] for i in m.nodes if (i, j) in m.arcs)
        - sum(m.q[j, i] for i in m.nodes if (j, i) in m.arcs)
        == m.D[j]
    )

model.con1 = Constraint(model.nodes, rule=flow_balance_rule)

# Original Hazen–Williams
def hw_original_rule(model, i, j):
    return (
        model.h[i] - model.h[j]
        ==
        model.q[i, j] * abs(model.q[i, j])**0.852
        * sum(
            model.omega * model.l[i, j, k]
            / (model.R[k]**1.852 * model.d[k]**4.87)
            for k in model.pipes
        )
    )

model.con2 = Constraint(model.arcs, rule=hw_original_rule)

# Smooth Approximation 1
def hw_smooth1_rule(model, i, j):
    return (
        model.h[i] - model.h[j]
        ==
        (
            model.q[i, j]**3
            * (model.q[i, j]**2 + model.eps[i, j]**2)**0.426
            / (model.q[i, j]**2 + 0.426 * model.eps[i, j]**2)
        )
        * sum(
            model.omega * model.l[i, j, k]
            / (model.R[k]**1.852 * model.d[k]**4.87)
            for k in model.pipes
        )
    )

# model.con2_approx1 = Constraint(model.arcs, rule=hw_smooth1_rule)

# Smooth Approximation 2
def hw_smooth2_rule(model, i, j):
    return (
        model.h[i] - model.h[j]
        ==
        model.q[i, j]
        * (model.q[i, j]**2 + model.eps[i, j]**2)**0.426
        * sum(
            model.omega * model.l[i, j, k]
            / (model.R[k]**1.852 * model.d[k]**4.87)
            for k in model.pipes
        )
    )

# model.con2 = Constraint(model.arcs, rule=hw_smooth2_rule)

# Pipe length consistency
def length_sum_rule(m, i, j):
    return sum(m.l[i, j, k] for k in m.pipes) == m.L[i, j]

model.con3 = Constraint(model.arcs, rule=length_sum_rule)

def length_bound_rule(m, i, j, k):
    return m.l[i, j, k] <= m.L[i, j]

model.con4 = Constraint(model.arcs, model.pipes, rule=length_bound_rule)

# Head constraints
def source_head_rule(m, i):
    return m.h[i] == m.E[i]

model.con6 = Constraint(model.Source, rule=source_head_rule)

def min_pressure_rule(m, i):
    if i in m.Source:
        return Constraint.Skip
    return m.h[i] >= m.E[i] + m.P[i]

model.con7 = Constraint(model.nodes, rule=min_pressure_rule)

# Flow bounds
def flow_bounds_rule(m, i, j):
    return (-m.Q_max, m.q[i, j], m.Q_max)

model.con8 = Constraint(model.arcs, rule=flow_bounds_rule)

# =========================
# Solve
# =========================
if __name__ == "__main__":
    data_file = sys.argv[1] if len(sys.argv) > 1 else "data.dat"
    instance = model.create_instance(data_file)
    solver = SolverFactory("knitro")
    # solver = SolverFactory("gurobi_persistent")
    # solver.set_instance(instance)
    # solver.options["outlev"] = 3
    # solver.options['NonConvex'] = 2
    # solver.options['FuncNonlinear'] = 1
    # solver.options['linear_solver'] = 'mumps'  # key option
    results = solver.solve(instance, tee=True)

    if (results.solver.status == SolverStatus.ok and
        results.solver.termination_condition == TerminationCondition.optimal):
        print("Total cost =", value(instance.total_cost))
        # instance.l.display()
        # instance.q.display()
        # instance.h.display() 
    else:
        print("Solver failed.")
        print("Status:", results.solver.status)
        print("Termination:", results.solver.termination_condition)
