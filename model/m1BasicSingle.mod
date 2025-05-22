#**************************************SETS******************************************#
set nodes;              # Set of nodes/vertexes
set pipes;              # Set of commercial pipes available
set arcs within {i in nodes, j in nodes: i != j};           # Set of arcs/links/edges

#***********************************PARAMETERS***************************************#
param L{arcs};          # Total length of each arc/link
param E{nodes};         # Elevation of each node
param P{nodes};         # Minimum pressure required at each node
param D{nodes};         # Demand of each node
param d{pipes};         # Diameter of each commercial pipe
param C{pipes};         # Cost per unit length of each commercial pipe
param R{pipes};          # Roughness of each commercial pipe
set Source;             # Source node ID

param omega := 10.67;   # SI Unit Constant for Hazen Williams Equation
param delta := 0.1;
param p:= 1.852;
param a := (15*(delta)^(p-1))/8 + ((p-1)*p*delta^(p-1))/8 - 7*p*(delta^(p-1))/8;
param b := (-5*(delta)^(p-3))/4 - ((p-1)*p*delta^(p-3))/4 + 5*p*(delta^(p-3))/4; 
param c := (3*(delta)^(p-5))/8 + ((p-1)*p*delta^(p-5))/8 - 3*p*(delta^(p-5))/8;
param vmax{arcs} default (sum {k in nodes diff Source} D[k]/1000)/((3.14/4)*(d[1]/1000)^2);
param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];

#***********************************VARIABLES****************************************#
var q{arcs};            # Flow variable
var h{nodes};           # Head
var X{arcs, pipes} binary;

#***********************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"

minimize total_cost : 
    sum{(i,j) in arcs} L[i,j]*(sum{k in pipes} C[k]*X[i,j,k]);

#**********************************CONSTRAINTS***************************************#

s.t. binaryCons {(i,j) in arcs}:
    sum {k in pipes} X[i,j,k] = 1;

s.t. con1{j in nodes diff Source}: 
    sum{i in nodes : (i,j) in arcs}q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j];

#s.t. con2{(i,j) in arcs}: 
#    h[i] - h[j] = (0.00000277971326776)*(q[i,j] * abs(q[i,j])^0.852) * omega * L[i,j]*(sum{k in pipes} (X[i,j,k] / (R[k]^1.852 * (d[k]/1000)^4.87)));

s.t. con2{(i,j) in arcs}: 
     h[i] - h[j]  = q[i,j]^3 *(((q[i,j]^2 + 0.001)^0.426) /(q[i,j]^2+0.00426)) * (0.00000277971326776) * omega * L[i,j] * sum{k in pipes} (X[i,j,k]/( (R[k]^1.852) * (d[k]/1000)^4.87));

# s.t. con2{(i,j) in arcs}: 
#    (if -0.1<=q[i,j]<=0.1  then 
#                                 (0.001^1.852)*(c*(q[i,j]^5) + b*(q[i,j]^3) + a*q[i,j])* omega * L[i,j] / ( ((sum{k in pipes}R[k]*X[i,j,k])^1.852) * ((sum{k in pipes}d[k]*X[i,j,k])/1000)^4.87)
#   else 
#         (q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) *omega * L[i,j] / ( ((sum{k in pipes}R[k]*X[i,j,k])^1.852) * ((sum{k in pipes}d[k]*X[i,j,k])/1000)^4.87)) = h[i] - h[j] ;

s.t. con3 {k in Source}: h[k] = E[k];

s.t. con4{i in nodes diff Source}: 
    h[i] >= E[i] + P[i];

subject to con5{(i,j) in arcs}: -Q_max <= q[i,j];
subject to con6{(i,j) in arcs}: q[i,j] <= Q_max;
#************************************************************************************#
