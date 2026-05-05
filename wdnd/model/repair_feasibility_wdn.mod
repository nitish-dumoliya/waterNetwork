# ============================================================
# EXACT REDUCED WDN DESIGN MODEL
# ============================================================

# =========================
# SETS
# =========================
set nodes;
set pipes ordered;
set arcs within {i in nodes, j in nodes: i != j};
set Source;

param NP integer > 1 := card(pipes);
set segs := 1..NP-1;

# Map ordered positions to actual pipe labels
#param pipe_id{1..NP} symbolic;

# =========================
# PARAMETERS
# =========================
param L{arcs};
param E{nodes};
param P{nodes};
param pmax{nodes};
param D{nodes};

param d{pipes};
param C{pipes};
param R{pipes};

param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);

param omega := 10.67;
param p := 1.852;

# Resistance coefficient alpha_k
param alpha{k in pipes} :=
    omega / (R[k]^1.852 * d[k]^4.87);

param alpha_min := min{k in pipes} alpha[k];
param alpha_max := max{k in pipes} alpha[k];

param R_min = min{k in pipes} R[k];
param d_min = min{k in pipes} d[k];
param c_min = min{k in pipes} C[k];
param c_max = max{k in pipes} C[k];
param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];

param MaxK{(i,j) in arcs} :=
    omega * L[i,j] / (R_min^1.852 * d_min^4.87);

param eps{(i,j) in arcs} :=
    0.0535 * (1e-3 / MaxK[i,j])^(0.54);

# Segment slope
param slope{s in segs} :=
    ( C[s] - C[s+1] ) /
    ( alpha[s] - alpha[s+1] );

# Segment intercept
param intercept{s in segs} :=
    (alpha[s]*C[s+1] - alpha[s+1] * C[s]) / ( alpha[s] - alpha[s+1] );

# =========================
# VARIABLES
# =========================
var q{arcs};
var h{nodes};
var y{arcs}>=alpha_min,<=alpha_max;
var z{arcs};
#var delta_flow_balance{nodes};
var delta_headloss{arcs};
var delta_y{arcs}>=0;
# =========================
# OBJECTIVE
# =========================

#minimize total_infeasibility:sum{i in nodes diff Source} delta_flow_balance[i]^2 + sum{(i,j) in arcs} delta_headloss[i,j]^2 ;
#minimize total_infeasibility: sum{(i,j) in arcs} (delta_headloss[i,j]^2 + delta_y[i,j]^2);
param rho >= 0 default 1;

#minimize infeas_obj:
#    rho * (
#        sum{(i,j) in arcs} (delta_headloss[i,j]^2)
#    );
# =========================
# CONSTRAINTS
# =========================

subject to flow_balance{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs} q[i,j] - sum{i in nodes : (j,i) in arcs} q[j,i] - D[j] = 0;

subject to con2{(i,j) in arcs}:
#    h[i] - h[j] - (q[i,j] * abs(q[i,j])^0.852) * y[i,j] * L[i,j] = 0;
    h[i] - h[j] - ( q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426 * eps[i,j]^2)) * y[i,j] * L[i,j] = delta_headloss[i,j];

subject to source_head{i in Source}:
    h[i] - E[i] = 0;

subject to min_pressure{i in nodes diff Source}:
    h[i] >= E[i] + P[i];

subject to exact_cost{(i,j) in arcs, s in segs}:
    z[i,j] >= L[i,j]*(slope[s] * y[i,j] + intercept[s]);

#subject to z_bounds_l{(i,j) in arcs}:
#    c_min*L[i,j] <= z[i,j] ;

#subject to z_bounds_r{(i,j) in arcs}:
#    z[i,j] <= c_max*L[i,j];

#subject to z_upper{(i,j) in arcs}:
#    z[i,j] <= L[i,j]*(c_max + ((c_min - c_max)/(alpha_max - alpha_min))*(y[i,j] - alpha_min));

subject to con8{(i,j) in arcs}: 
   -Q_max <= q[i,j] <= Q_max
;

