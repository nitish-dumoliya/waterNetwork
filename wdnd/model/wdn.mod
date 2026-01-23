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
param omega := 10.67;  # SI Unit Constant for Hazen Williams Equation
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);
param p:= 1.852;
param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];
#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0;	# Length of each commercial pipe for each arc/link
var q{arcs};	            # Flow variable
var h{nodes};	            # Head
#var y{arcs}>=1e-4;
var q1{arcs};	            # Flow variable
var q2{arcs};	            # Flow variable

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	

#****************************************CONSTRAINTS**************************************#
subject to con5{i in Source}: 
    h[i] = E[i]
;

#subject to con6{i in nodes diff Source}: h[i] >= (E[i] + P[i]) ;
#subject to con6_{i in nodes diff Source}: h[i] <= max{j in Source} E[j] ;
subject to con7_{i in nodes diff Source}: h[i] >= 0 ;
#subject to con7{(i,j) in arcs}: -Q_max <= q[i,j];
#subject to con8{(i,j) in arcs}: q[i,j] <= Q_max;

#*******************************************************************************************#
