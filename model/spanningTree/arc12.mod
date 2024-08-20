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
param d{pipes};		   # Diameter of each commercial pipe
param C{pipes};		   # Cost per unit length of each commercial pipe
param R{pipes};		   # Roughness of each commercial pipe
param omega := 10.68;  # SI Unit Constant for Hazen Williams Equation
param vmax{arcs} default (sum {k in nodes diff Source} D[k]/1000)/((3.14/4)*(d[1]/1000)^2);

param Qmax:= sum{k in nodes diff Source} D[k];
#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{k in pipes}l[1,2,k]*C[k];	

#****************************************CONSTRAINTS**************************************#

subject to arc12_head: 
    sum{k in pipes } (omega * l[1,2,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) <= (E[1]-(E[2]+P[2]))/((Qmax * abs(Qmax)^0.852) * (0.001^1.852))
;


subject to lengthCon: 
    sum{k in pipes} l[1,2,k] = L[1,2]
;
