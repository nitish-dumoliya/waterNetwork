#*******************************************SETS******************************************#
set nodes;			   # Set of nodes/vertices
set pipes;			   # Set of commercial pipes available
set arcs within {i in nodes, j in nodes: i != j};	# Set of arcs/links/edges
set Source;		       # Set of source nodes 

#****************************************PARAMETERS***************************************#
param L{arcs };		   # Total length of each arc/link
param E{nodes};		   # Elevation of each node
param P{nodes};		   # Minimum pressure required at each node
param D{nodes};		   # Demand of each node
param C{pipes};		   # Cost per unit length of each commercial pipe
param R{pipes};		   # Roughness of each commercial pipe
param omega := 10.68;  # SI Unit Constant for Hazen Williams Equation
param vmax{arcs};

#****************************************VARIABLES****************************************#
var q{arcs};	            # Flow variable
var h{nodes};	            # Head
var d{arcs};		   # Diameter of each arcs

#****************************************OBJECTIVE****************************************#
#minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	
minimize total_cost : 0;	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes}: 
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;
subject to con2{(i,j) in arcs}: 
    h[i] - h[j] = (q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * (omega * L[i,j] / ( (130^1.852) * (d[i,j])^4.87))
;

subject to con5{i in Source}: 
    h[i] = E[i]
;
subject to con6{i in nodes diff Source}: 
    h[i] >= E[i] + P[i]
;
subject to con7{(i,j) in arcs}:
    -sum{k in nodes diff Source} D[k] <= q[i,j]
;
subject to con8{(i,j) in arcs}:
    q[i,j] <= sum{k in nodes diff Source} D[k]
;


subject to con9{(i,j) in arcs}:
    0.0254 <= d[i,j] <= 0.6096;
#*******************************************************************************************#
