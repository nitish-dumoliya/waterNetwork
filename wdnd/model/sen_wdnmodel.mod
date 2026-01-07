reset;

# ======================= sIPOPT SUFFIXES =========================
suffix sens_state_0, IN;
suffix sens_state_1, IN;
suffix sens_state_value_1, IN;
suffix sens_sol_state_1, OUT;
suffix sens_init_constr, IN;

option presolve 0;

# ======================= SETS ====================================
set nodes;
set pipes;
set arcs within {i in nodes, j in nodes: i != j};
set Source;

# ======================= PARAMETERS ==============================
param L{arcs};
param E{nodes};
param P{nodes};
param pmax{nodes};
param D{nodes};
param C{pipes};
param R{pipes};

param omega := 10.67;
param p := 1.852;

param Q_max = sum{k in nodes diff Source} D[k];

# ---------- NOMINAL DIAMETER (IMPORTANT) -------------------------
param d{pipes};
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);

param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];
param d_max = max{i in pipes} d[i];
param d_min = min{i in pipes} d[i];

param R_min = min{k in pipes} R[k];

param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
#param eps{(i,j) in arcs} := (1e-5 / (0.07508 * MaxK[i,j]))^(1/1.852);
#param eps{(i,j) in arcs} := 0.0953*(1e-2/MaxK[i,j])^(0.54);
param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);

# ======================= VARIABLES ================================
var l{arcs,pipes} >= 0;
var q{arcs};
var h{nodes};

# 👉 Diameter is now a VARIABLE (not parameter!)
var nominal_d{pipes}>=d_min, <= d_max;

# ======================= OBJECTIVE ================================
minimize total_cost:
    sum{(i,j) in arcs, k in pipes} l[i,j,k] * C[k];

# ======================= CONSTRAINTS ==============================

subject to flow_balance{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs} q[i,j] - sum{i in nodes : (j,i) in arcs} q[j,i] = D[j];

# Hazen–Williams (smooth version)
subject to headloss{(i,j) in arcs}:
    h[i] - h[j]  =  (q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426*eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * nominal_d[k]^4.87));

subject to length_cons{(i,j) in arcs}:
    sum{k in pipes} l[i,j,k] = L[i,j];

subject to length_bounds{(i,j) in arcs , k in pipes}: 
    l[i,j,k] <= L[i,j]
;
subject to source_head{i in Source}:
    h[i] = E[i];

subject to min_pressure{i in nodes diff Source}:
    h[i] >= E[i] + P[i];

subject to flow_bounds{(i,j) in arcs}:
    -Q_max <= q[i,j] <= Q_max;

# ======================= INITIAL CONSTRAINTS ======================
# These are REQUIRED for sIPOPT

subject to const_d{k in pipes}:
    nominal_d[k] = d[k];

#*******************************************************************************************#
