#*******************************************SETS******************************************#
set nodes;                             # Nodes
set pipes;                             # Commercial pipe types
set arcs within {i in nodes, j in nodes: i != j};
set Source;                            # Source nodes

#***************************************PARAMETERS***************************************#
param L{arcs};                         # Arc length
param E{nodes};                        # Elevation
param P{nodes};                        # Minimum pressure
param D{nodes};                        # Demand

param d{pipes};                        # Pipe diameter
param C{pipes};                        # Cost per unit length
param R{pipes};                        # Roughness
param pmax{nodes};		   # Maximum pressure required at each node
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);

param omega := 10.67;                  # Hazen–Williams constant
param p := 1.852;                      # HW exponent

param Q_max = sum{k in nodes diff Source} D[k];

#param eta = 0.4;     # relative trust region for l
#param Delta = 0.2;   # absolute trust region for q
param y_ref{arcs} default 0;
param s_ref{arcs} default 0;

#-------------------------------- Energy coefficient ------------------------------------#
# alpha_k = ω / (2 (p+1) R^p d^4.87)
param alpha{k in pipes} :=
    omega / (2 * (p + 1)) / (R[k]^p * d[k]^4.87);

#****************************************VARIABLES****************************************#
var l{arcs, pipes} >= 0;               # Length assigned to pipe type
var q{arcs};                      # Flow (direction fixed for convexity)
var h{nodes};                          # Hydraulic head

var z{arcs};                      # Flow (direction fixed for convexity)
#***************************************OBJECTIVE****************************************#
# Convex cost + convex energy dissipation
minimize total_energy_cost: (sum{(i,j) in arcs} sum{k in pipes} C[k] * l[i,j,k]) + (sum {(i,j) in arcs} z[i,j]) - sum{k in Source} E[k]*(sum{i in nodes: (k,i) in arcs} q[k,i]);

#**************************************CONSTRAINTS***************************************#

#----------------------------- Flow conservation ----------------------------------------#
subject to flow_balance{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs} q[i,j] - sum{i in nodes : (j,i) in arcs} q[j,i] = D[j];

#----------------------------- Length partition -----------------------------------------#
#subject to length_sum{(i,j) in arcs}:
#    sum{k in pipes} l[i,j,k] = L[i,j];

subject to length_bound{(i,j) in arcs, k in pipes}:
    l[i,j,k] <= L[i,j];

#subject to head_lb{(i,j) in arcs}:
#    abs(h[i] - h[j]) >= z[i,j];
##----------------------------- Head reference --------------------------------------------#
#subject to source_head{i in Source}:
#    h[i] = E[i];
#
#subject to min_pressure{i in nodes diff Source}:
#    h[i] >= E[i] + P[i];

#----------------------------- Flow bound ------------------------------------------#
#subject to flow_cap{(i,j) in arcs}:
#    -Q_max <= q[i,j] <= Q_max;

