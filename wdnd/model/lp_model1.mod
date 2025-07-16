### SETS ###
set nodes;			### Set of nodes/vertexes
set pipes;			### Set of commercial pipes available
set arcs within {i in nodes, j in nodes: i != j};	### Set of arcs/links/edges
set   Source;		### Source node ID

### PARAMETERS ###
param L{arcs};		### Total length of each arc/link
param E{nodes};		### Elevation of each node
param P{nodes};		### Minimum pressure required at each node
param D{nodes};		### Demand of each node
param d{pipes};		### Diameter of each commercial pipe
param C{pipes};		### Cost per unit length of each commercial pipe
param R{pipes};		### Roughness of each commercial pipe
param vmax{arcs} default (sum {k in nodes diff Source} D[k]/1000)/((3.14/4)*(d[1]/1000)^2);

param omega := 10.68;			### SI Unit Constant for Hazen Williams Equation

param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];
param d_max = max{i in pipes} d[i];
param d_min = min{i in pipes} d[i];
param pmax{nodes};		   # Maximum pressure required at each node

### VARIABLES ###
#var l_lp{arcs,pipes} >= 0 ;			### Length of each commercial pipe for each arc/link
var q_lp{arcs},>=-Q_max,<=Q_max;	### Flow variable
#var h_lp{nodes};					### Head

### OBJECTIVE ###
minimize total_cost : 0;	### Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"

### CONSTRAINTS ###

s.t. con1{j in nodes diff Source}: sum{i in nodes : (i,j) in arcs}q_lp[i,j] -  sum{i in nodes : (j,i) in arcs}q_lp[j,i] =  D[j];


#s.t. con1{(i,j) in arcs}: (h_lp[i] - h_lp[j])= ((q_lp[i,j] * abs(q_lp[i,j])^0.852) * (0.001^1.852) ) * sum{k in pipes } omega * l_lp[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87);

#s.t. con2{(i,j) in arcs}: sum{k in pipes} l_lp[i,j,k] = L[i,j];

#s.t. con3{i in Source}: h_lp[i] = E[i];

#s.t. con4{i in nodes diff Source}: h_lp[i] >= E[i] + P[i];

#s.t. con5{(i,j) in arcs, k in pipes}: l_lp[i,j,k] <= L[i,j];

# for d.dat dataset 
########################### Iteration 1  #####################
