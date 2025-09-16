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
param R{arcs};		   # Roughness of each commercial pipe
param omega := 10.67;  # SI Unit Constant for Hazen Williams Equation
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);
param delta := 0.1;
param exdiam{arcs};
param excost{arcs};
param p:= 1.852;
param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];
param d_min = min{i in pipes} d[i];
param d_max = max{i in pipes} d[i];

param R_min = min{(i,j) in arcs} R[i,j];

param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);

param eps{(i,j) in arcs} := (1e-6 / (0.07508 * MaxK[i,j]))^(1 / 0.926);


#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs},>=-Q_max,<=Q_max;	            # Flow variable
var q1{arcs},>=-Q_max,<=Q_max;	            # Flow variable
var q2{arcs},>=-Q_max,<=Q_max;	            # Flow variable
#var h{nodes};	            # Head
#var eps{arcs}>=0;

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
#minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k] ;	
minimize total_cost : sum{(i,j) in arcs} ((sum{k in pipes} (omega * l[i,j,k] / (R[i,j]^1.852) * d[k]^4.87))  * (abs(q1[i,j])^(2.852))/2.852) + sum{(i,j) in arcs} (omega * L[i,j] / ( (R[i,j]^1.852) * (exdiam[i,j]^4.87))*((abs(q2[i,j])^2.852 / 2.852))) ; # + sum{u in Source} sum{(i,j) in arcs} E[u]* q[i,j];  #sum{u in Source} E[u] * (sum{(i,u) in arcs} q[i,u] - sum{(j,u) in arcs} q[j,u]);

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs }(q1[i,j] + q2[i,j]) -  sum{i in nodes : (j,i) in arcs}(q1[j,i] + q2[j,i]) =  D[j]
;
subject to con3{(i,j) in arcs}: 
    sum{k in pipes} l[i,j,k] = L[i,j]
;

subject to con4{(i,j) in arcs , k in pipes}: 
    l[i,j,k] <= L[i,j]
;

#subject to con5{i in Source}: 
#    h[i] = E[i]
#;

#subject to con6{i in nodes diff Source}: h[i] >= (E[i] + P[i]) ;

#subject to con7{(i,j) in arcs}: -Q_max <= q[i,j];
#subject to con8{(i,j) in arcs}: q[i,j] <= Q_max;

subject to con9{(i,j) in arcs}: q[i,j] = q1[i,j] + q2[i,j];

#subject to con10{(i,j) in arcs}: q1[i,j]*q2[i,j] >= 0;
#subject to con10{(i,j) in arcs}: q1[i,j]*(q2[i,j]^2 + eps[i,j])^0.5 = q2[i,j]*(q1[i,j]^2 + eps[i,j])^0.5;
#*******************************************************************************************#
