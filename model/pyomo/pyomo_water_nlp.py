from pyomo.environ import *
#from pyomo.environ import smooth_abs

# Create model
model = ConcreteModel()

# ******************** SETS ******************** #
model.nodes = Set(initialize=[1, 2, 3, 4, 5, 6, 7])
model.pipes = Set(initialize=range(1, 15))
model.arcs = Set(initialize=[(1, 2), (2, 4), (2, 3), (3, 5), (4, 5), (4, 6), (7, 5), (6, 7)])
model.Source = Set(initialize=[1])

# ******************** PARAMETERS ******************** #
L_data = {(1,2):1000, (2,4):1000, (2,3):1000, (3,5):1000, (4,5):1000, (4,6):1000, (7,5):1000, (6,7):1000}
model.L = Param(model.arcs, initialize=L_data)

E_data = {1:210, 2:150, 3:160, 4:155, 5:150, 6:165, 7:160}
model.E = Param(model.nodes, initialize=E_data)

P_data = {1:0, 2:30, 3:30, 4:30, 5:30, 6:30, 7:30}
model.P = Param(model.nodes, initialize=P_data)

D_data = {1:-311.1087, 2:27.7777, 3:27.777, 4:33.333, 5:75, 6:91.666, 7:55.555}
model.D = Param(model.nodes, initialize=D_data)

d_data = {i: val for i, val in enumerate([25.4, 50.8, 76.2, 101.6, 152.4, 203.2, 254, 304.8, 355.6, 406.4, 457.2, 508, 558.8, 609.6], 1)}
model.d = Param(model.pipes, initialize=d_data)

C_data = {i: val for i, val in enumerate([2, 5, 8, 11, 16, 23, 32, 50, 60, 90, 130, 170, 300, 550], 1)}
model.C = Param(model.pipes, initialize=C_data)

R_data = {i: 130 for i in model.pipes}
model.R = Param(model.pipes, initialize=R_data)

model.omega = Param(initialize=10.67)
model.delta = Param(initialize=0.1)
model.p = Param(initialize=1.852)
model.eps = Param(initialize=1e-10)

Q_max = sum(D_data[k] for k in model.nodes if k not in model.Source)
model.Q_max = Param(initialize=Q_max)

# ******************** VARIABLES ******************** #
model.l = Var(model.arcs, model.pipes, within=NonNegativeReals)
model.q = Var(model.arcs, within=Reals)
model.h = Var(model.nodes, within=Reals)

# ******************** OBJECTIVE ******************** #
def objective_rule(model):
    return sum(model.l[i, j, k] * model.C[k] for (i, j) in model.arcs for k in model.pipes)
model.total_cost = Objective(rule=objective_rule, sense=minimize)

# ******************** CONSTRAINTS ******************** #
# Flow conservation constraint
def con1_rule(model, j):
    return sum(model.q[i, j] for i in model.nodes if (i, j) in model.arcs) - \
           sum(model.q[j, i] for i in model.nodes if (j, i) in model.arcs) == model.D[j]
model.con1 = Constraint(model.nodes, rule=con1_rule)

def smooth_sign(q, eps=1e-6):
    """Smooth approximation of sign(q) for Pyomo."""
    return q / ((q**2 + eps)**0.5)

# Original (Head-loss constraint)Hazen-Williams constraint 
def con2(model, i, j):
    return (model.h[i] - model.h[j] == 
            (model.q[i, j]) * abs(model.q[i, j]) ** 0.852 * (0.001 ** 1.852) *
            sum(model.omega * model.l[i, j, k] / ((model.R[k] ** 1.852) * (model.d[k] / 1000) ** 4.87) 
                for k in model.pipes))

# First order approximation of Head-loss constraint
def con2_approx1(model, i, j):
    term1 = 0.00000278 * model.q[i, j] * (abs(model.q[i, j]) + 1000 * model.eps) ** 0.852
    term2 = 0.002368316 * model.eps * model.q[i, j] / (abs(model.q[i, j]) + 1000 * model.eps) ** 0.148
    sum_term = sum(model.omega * model.l[i, j, k] / ((model.R[k] ** 1.852) * (model.d[k] / 1000) ** 4.87) 
                   for k in model.pipes)
    return model.h[i] - model.h[j] == (term1 - term2) * sum_term

# First order approximation of Head-loss constraint
def con2_approx2(model, i, j):
    return (model.h[i] - model.h[j] == 
            ((model.q[i, j] * (abs(model.q[i, j]) + 148 * model.eps)) / 
             ((abs(model.q[i, j]) + 1000 * model.eps) ** 0.148)) * (0.001 ** 1.852) *
            sum(model.omega * model.l[i, j, k] / ((model.R[k] ** 1.852) * (model.d[k] / 1000) ** 4.87) 
                for k in model.pipes))

# Second order approximation of Head-loss constraint
def con2_approx3(model, i, j):
    expr = (0.001 ** 1.852) * model.q[i, j] * (abs(model.q[i, j]) + 1000 * model.eps) ** 0.852
    expr -= (0.002368316 * model.eps * model.q[i, j] / (abs(model.q[i, j]) + 1000 * model.eps) ** 0.148)
    expr += ((0.175255362 * (model.eps) ** 2) * model.q[i, j] / ((abs(model.q[i, j]) + 1000 * model.eps) ** 1.148))
    return model.h[i] - model.h[j] == expr * sum(model.omega * model.l[i, j, k] / ((model.R[k] ** 1.852) * (model.d[k] / 1000) ** 4.87) for k in model.pipes)

model.con2 = Constraint(model.arcs, rule=con2_approx2)  #con2_approx3

model.con3 = Constraint(model.arcs, rule=lambda model, i, j: sum(model.l[i, j, k] for k in model.pipes) == model.L[i, j])

model.con4 = Constraint(model.arcs, model.pipes, rule=lambda model, i, j, k: model.l[i, j, k] <= model.L[i, j])

model.con5 = Constraint(model.Source, rule=lambda model, i: model.h[i] == model.E[i])

model.con6 = Constraint(model.nodes - model.Source, rule=lambda model, i: model.h[i] >= model.E[i] + model.P[i])

model.con7 = Constraint(model.arcs, rule=lambda model, i, j: model.q[i, j] >= -model.Q_max)
model.con8 = Constraint(model.arcs, rule=lambda model, i, j: model.q[i, j] <= model.Q_max)

# ******************** SOLVE MODEL ******************** #
def solve_model():
    solver = SolverFactory('ipopt')
    results = solver.solve(model, tee=True)
    model.display()

    print(f"Objective Value: {model.total_cost()}")  

if __name__ == "__main__":
    solve_model()
