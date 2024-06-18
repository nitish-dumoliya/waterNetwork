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

### VARIABLES ###
var l_lp{arcs,pipes} >= 0 ;			### Length of each commercial pipe for each arc/link
param q_lp{arcs};	### Flow variable
var h_lp{nodes};					### Head

### OBJECTIVE ###
minimize total_cost : sum{(i,j) in arcs } sum{k in pipes}l_lp[i,j,k]*C[k];	### Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"

### CONSTRAINTS ###
s.t. con3{(i,j) in arcs}: h_lp[i] - h_lp[j] = (q_lp[i,j] * abs(q_lp[i,j])^0.852) * (0.001^1.852) * sum{k in pipes } omega * l_lp[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87);

s.t. con4{(i,j) in arcs}: sum{k in pipes} l_lp[i,j,k] = L[i,j];

s.t. con5{i in Source}: h_lp[i] = E[i];

s.t. con2{i in nodes diff Source}: h_lp[i] >= E[i] + P[i];

s.t. bound1{(i,j) in arcs, k in pipes}: l_lp[i,j,k] <= L[i,j];

# for d.dat dataset 
########################### Iteration 1  #####################
