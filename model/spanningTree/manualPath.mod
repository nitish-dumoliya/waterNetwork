#**************************************SETS******************************************#
set nodes;              # Set of nodes/vertexes
set pipes;              # Set of commercial pipes available
set arcs within {i in nodes, j in nodes: i != j};           # Set of arcs/links/edges
set Source;             # Source node ID
set infSet := {1, 28, 24, 23};

#***********************************PARAMETERS***************************************#
param L{arcs};          # Total length of each arc/link
param E{nodes};         # Elevation of each node
param P{nodes};         # Minimum pressure required at each node
param D{nodes};         # Demand of each node
param d{pipes};         # Diameter of each commercial pipe
param C{pipes};         # Cost per unit length of each commercial pipe
param R{pipes};          # Roughness of each commercial pipe
param omega := 10.68;   # SI Unit Constant for Hazen Williams Equation
param delta := 0.1;
param p:= 1.852;
param a := (15*(delta)^(p-1))/8 + ((p-1)*p*delta^(p-1))/8 - 7*p*(delta^(p-1))/8;
param b := (-5*(delta)^(p-3))/4 - ((p-1)*p*delta^(p-3))/4 + 5*p*(delta^(p-3))/4; 
param c := (3*(delta)^(p-5))/8 + ((p-1)*p*delta^(p-5))/8 - 3*p*(delta^(p-5))/8;
param vmax{arcs} default (sum {k in nodes diff Source} D[k]/1000)/((3.14/4)*(d[1]/1000)^2);

var v{nodes};

#let D[1] := sum{i in nodes diff Source} D[i]*v[i];

#***********************************VARIABLES****************************************#
var q{arcs};            # Flow variable
var h{nodes}>=0;           # Head
var X{arcs, pipes} binary;
var x{arcs} binary;
#var dia{arcs};          # Diameter of each pipe
#var R_{arcs};
#***********************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"

minimize total_cost : 
    sum{(i,j) in arcs} L[i,j]*(sum{k in pipes} C[k]*X[i,j,k])*x[i,j];

#**********************************CONSTRAINTS***************************************#

s.t. binaryCons {(i,j) in arcs}:
    sum {k in pipes} X[i,j,k]  = 1;

#s.t. diameterCons {(i,j) in arcs}:
#    sum {k in pipes} d[k]*X[i,j,k] = dia[i,j];

#subject to Roughness_cons{(i,j) in arcs}:
#    sum {k in pipes} R[k]*X[i,j,k] = R_[i,j];

s.t. con1{j in nodes diff Source}: 
    sum{i in nodes : (i,j) in arcs}q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]*v[j];

s.t. con2{(i,j) in arcs}: 
    (h[i] - h[j])*x[i,j] = (0.001^1.852)*(q[i,j] * abs(q[i,j])^0.852) * omega * L[i,j]*(sum{k in pipes} (X[i,j,k] / (R[k]^1.852 * (d[k]/1000)^4.87)));

# s.t. con3{(i,j) in arcs}: 
#    (if -0.1<=q[i,j]<=0.1  then 
#                                 (0.001^1.852)*(c*(q[i,j]^5) + b*(q[i,j]^3) + a*q[i,j])* omega * L[i,j] / ( ((sum{k in pipes}R[k]*X[i,j,k])^1.852) * ((sum{k in pipes}d[k]*X[i,j,k])/1000)^4.87)
#   else 
#         (q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) *omega * L[i,j] / ( ((sum{k in pipes}R[k]*X[i,j,k])^1.852) * ((sum{k in pipes}d[k]*X[i,j,k])/1000)^4.87)) = h[i] - h[j] ;

s.t. con3 {k in Source}: h[k] = E[k];

s.t. con4{i in nodes diff Source}: 
    (h[i] - E[i] - P[i])*v[i] >=0;

s.t. flow_bound1 {(i,j) in arcs}:
    (-sum{k in nodes diff Source} D[k])*x[i,j] <= q[i,j] ;

s.t. flow_bound2 {(i,j) in arcs}:
    q[i,j] <= (sum {k in nodes diff Source} D[k])*x[i,j];

#************************************************************************************#
