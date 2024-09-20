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

#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;    # Length of each commercial pipe for each arc/link
var q1{arcs}>=0;	            # Flow direction to original arc
var q2{arcs}>=0;	            # Opposite flow direction to original arc
var h{nodes};			    # Head
var z{arcs} binary;

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes} l[i,j,k]*C[k];	

#****************************************CONSTRAINTS**************************************#
s.t. con1{j in nodes }: sum{i in nodes : (i,j) in arcs }(q1[i,j] - q2[i,j]) -  sum{i in nodes : (j,i) in arcs}(q1[j,i] - q2[j,i]) =  D[j];

s.t. con2{(i,j) in arcs}: h[i] - h[j] = (q1[i,j]^1.852- q2[i,j]^1.852) * (0.001^1.852) * sum{k in pipes } omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87);

s.t. con3{(i,j) in arcs }: sum{k in pipes} l[i,j,k] = L[i,j];

s.t. con4{(i,j) in arcs , k in pipes}: l[i,j,k] <= L[i,j];

s.t. con5{i in Source}: h[i] = E[i];

s.t. con6{i in nodes diff Source}: h[i] >= E[i] + P[i];

s.t. con8{(i,j) in arcs }: q1[i,j] <= (sum{k in nodes diff Source} D[k])*z[i,j];

s.t. con9{(i,j) in arcs }: q2[i,j] <= (sum{k in nodes diff Source} D[k])*(1-z[i,j])  ;

# s.t. q12: q[1,2] = sum{k in nodes diff Source} D[k] ;

#*****************************************************************************************
