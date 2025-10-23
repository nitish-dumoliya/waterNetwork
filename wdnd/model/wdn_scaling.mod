#*******************************************SETS******************************************#
set nodes;
set pipes;
set arcs within {i in nodes, j in nodes: i != j};
set Source;

#****************************************PARAMETERS***************************************#
param L{arcs};
param E{nodes};
param P{nodes};
param pmax{nodes};
param D{nodes};
param d{pipes};
param C{pipes};
param R{pipes};
param omega := 10.67;
param p := 1.852;
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);
param Q_max = sum{k in nodes diff Source} D[k];
param D_max = max{i in nodes diff Source} D[i];
param D_min = min{i in nodes diff Source} D[i];
param R_min = min{k in pipes} R[k];
param d_min = min{i in pipes} d[i];
param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);

param L_max = max{(i,j) in arcs} L[i,j];
# Scaling parameter
param s_q := Q_max;
param scale := 1000;
#param eps{(i,j) in arcs} := (1e-6 / (0.07508 * MaxK[i,j]))^(1 / 1.852);
param eps{(i,j) in arcs} := (1e-4 / (0.36061 * MaxK[i,j]))^(0.54);
#param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);
#param eps{(i,j) in arcs} := 1.7235*(1e-3/MaxK[i,j])^(0.54);
#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0;
#var qs{(i,j) in arcs} >= -1, <= 1;
var h{nodes};
var q{arcs}>=-s_q,<=s_q;
var x{arcs};
var y{arcs};
#****************************************OBJECTIVE****************************************#
minimize total_cost:
    sum{(i,j) in arcs} sum{k in pipes} l[i,j,k] * C[k];

#****************************************CONSTRAINTS**************************************#
subject to con1_scaled{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs} q[i,j] - sum{i in nodes : (j,i) in arcs} q[j,i] = D[j];

#subject to con2_scaled{(i,j) in arcs}:
#    (h[i] - h[j]) = (q[i,j]^3 * (((q[i,j]/eps[i,j])^2 + 1)^0.426) / (eps[i,j]^1.148 * ((q[i,j]/eps[i,j])^2 + 0.426 * 1))) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
subject to con2_scaled{(i,j) in arcs}:
    #(h[i] - h[j]) = ((q[i,j])^3 * (((q[i,j])^2 + (eps[i,j])^2)^0.426) / ((q[i,j])^2 + 0.426 * eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #(h[i] - h[j])*(scale)^1.852 = ((q[i,j]*scale)^3 * ((q[i,j]*scale)^2 + (eps[i,j]*scale)^2)^0.426 / ((q[i,j]*scale)^2 + 0.426*(eps[i,j]*scale)^2 )) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #(h[i] - h[j])*(scale)^1.852 = ((x[i,j])^3 * ((x[i,j])^2 + (y[i,j])^2)^0.426 / ((x[i,j])^2 + 0.426*(y[i,j])^2 )) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #(h[i] - h[j]) =eps[i,j]^1.852 * ((q[i,j]/eps[i,j])^3 * ((q[i,j]/eps[i,j])^2 + (1)^2)^0.426 / ((q[i,j]/eps[i,j])^2 + 0.426*(1)^2 )) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    (h[i] - h[j])*(scale^0.852)  =  (q[i,j] * ((scale*q[i,j])^2 + (scale*eps[i,j])^2)^0.426) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #(h[i] - h[j])/(eps[i,j]^0.148)  =  ((q[i,j]/eps[i,j]) * ((q[i,j]/eps[i,j])^2 + (1)^2)^0.426) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));

subject to con3{(i,j) in arcs}:
    sum{k in pipes} l[i,j,k] = L[i,j];

subject to con4{(i,j) in arcs, k in pipes}:
    l[i,j,k] <= L[i,j];

subject to con5{i in Source}:
    h[i] = E[i];

subject to con6{i in nodes diff Source}:
    h[i] >= (E[i] + P[i]);

#subject to con7{(i,j) in arcs}:
#    x[i,j] = q[i,j]*scale;
#subject to con8{(i,j) in arcs}:
#    y[i,j] = eps[i,j]*scale;



