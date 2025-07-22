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
param eps{(i,j) in arcs} := (1e-6 / (0.07508 * MaxK[i,j]))^(1 / 0.926);
#param eps{(i,j) in arcs} := (1e-3 / (0.04001571 * MaxK[i,j]))^(1 / 1.852);
#param eps{(i,j) in arcs} := (1e-4 * R_min^1.852 * d_min^4.87 / (0.07508 * 10.67 * L[i,j]))^(1 / 0.926);

#****************************************VARIABLES****************************************#
var l{arcs,pipes} ;	# Length of each commercial pipe for each arc/link
var q{arcs};	            # Flow variable
#var h{nodes};	            # Head
#var eps{arcs};
#var x{arcs, pipes}>=0,<=1;

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
#minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k] + sum{(i,j) in arcs} sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87))  * (abs(q[i,j])^(2.852))/2.852 ; #+ sum{u in Source} E[u] * (sum{(i,u) in arcs} q[i,u] - sum{(j,u) in arcs} q[j,u]);

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;
#subject to con2{i in Source}: 
#    h[i] = E[i]
#;

subject to con3{(i,j) in arcs}: 
    sum{k in pipes} l[i,j,k] = L[i,j]
;

subject to con4{(i,j) in arcs , k in pipes}: 
    0 <= l[i,j,k] <= L[i,j]
;
#subject to con7{(i,j) in arcs}: -Q_max <= q[i,j];
#subject to con8{(i,j) in arcs}: q[i,j] <= Q_max;

#*******************************************************************************************#
