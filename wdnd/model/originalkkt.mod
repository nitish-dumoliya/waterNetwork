#*******************************************SETS******************************************#
set nodes;			   # Set of nodes/vertices
set pipes;			   # Set of commercial pipes available
set arcs within {i in nodes, j in nodes: i != j};	# Set of arcs/links/edges
set Source;		       # Set of source nodes 

#****************************************PARAMETERS***************************************#
param L{arcs };		   # Total length of each arc/link
param E{nodes};		   # Elevation of each node
param P{nodes};		   # Minimum pressure required at each node
param pmax{nodes};		   # Maximum pressure required at each node
param D{nodes};		   # Demand of each node
param d{pipes};		   # Diameter of each commercial pipe
param C{pipes};		   # Cost per unit length of each commercial pipe
param R{pipes};		   # Roughness of each commercial pipe
param omega := 10.67;  # SI Unit Constant for Hazen Williams Equation
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);
param p:= 1.852;
param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];
param d_max = max{i in pipes} d[i];
param d_min = min{i in pipes} d[i];
param delta := 0.01;
param R_min = min{k in pipes} R[k];
param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
#param eps{(i,j) in arcs} := 0.0535*(1e-4/MaxK[i,j])^(0.54);
param eps{arcs};
#****************************************VARIABLES****************************************#
param l{arcs,pipes};	# Length of each commercial pipe for each arc/link
param q{arcs};	            # Flow variable
param h{nodes};	            # Head

var lam{nodes};
var x{arcs};
var y{arcs};
var u{Source};
var v{nodes diff Source}>=0;
var w{arcs, pipes}>=0;
#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : 0;	

#****************************************CONSTRAINTS**************************************#

#Stationay conditions
subject to con1{(i,j) in arcs , k in pipes}:
    C[k] - 10.67*q[i,j]*abs(q[i,j])^0.852 * x[i,j] / (R[k]^1.852 * d[k]^4.87) + y[i,j]  - w[i,j,k]= 0
    #C[k] - (10.67 * q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426*eps[i,j]^2) * x[i,j] / (R[k]^1.852 * d[k]^4.87)) + y[i,j] = 0
;

subject to con2{(i,j) in arcs}:
    lam[j] - lam[i] - 1.852 * x[i,j] * abs(q[i,j])^0.852 * (sum{k in pipes} 10.67 * l[i,j,k] / (R[k]^1.852 * d[k]^4.87)) = 0
    #lam[j] - lam[i] - x[i,j] * ( (q[i,j]^2 + 0.426*eps[i,j]^2) * ( 3*q[i,j]^2*(q[i,j]^2 + eps[i,j]^2)^0.426 + 0.852*q[i,j]^4*(q[i,j]^2 + eps[i,j]^2)^(-0.574) ) - 2*q[i,j]^4*(q[i,j]^2 + eps[i,j]^2)^0.426 ) / (q[i,j]^2 + 0.426*eps[i,j]^2)^2 * (sum{k in pipes} 10.67 * l[i,j,k] / (R[k]^1.852 * d[k]^4.87)) = 0
;

subject to con3{s in Source}:
    - (sum{j in nodes: (s,j) in arcs} x[s,j]) + u[s] = 0
;

subject to con4{j in nodes diff Source}:
    - sum{i in nodes: (i,j) in arcs} x[i,j] + sum{i in nodes: (j,i) in arcs} x[j,i] - v[j] = 0
;

# Primal Feasibility conditions
#subject to con5{j in nodes diff Source}:
#    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
#;

#subject to con6{(i,j) in arcs}: 
#    #h[i] - h[j]  =  (q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426*eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
#    h[i] - h[j]  =  q[i,j]*abs(q[i,j])^0.852 * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));

#subject to con7{(i,j) in arcs}: 
#    sum{k in pipes} l[i,j,k] = L[i,j]
#;

#subject to con8{(i,j) in arcs , k in pipes}: 
#    l[i,j,k] <= L[i,j]
#;

#subject to con9{i in Source}: 
#    h[i] = E[i]
#;

#subject to con10{i in nodes diff Source}: h[i] >= (E[i] + P[i]) ;

# Dual Feasibility conditions
#subject to con11{j in nodes diff Source}:
#    v[j] >= 0
#;

# Complementary sleckness conditions
subject to con12{j in nodes diff Source}:
    (-h[j] + E[j] + P[j]) * v[j] = 0
;

#*******************************************************************************************#
