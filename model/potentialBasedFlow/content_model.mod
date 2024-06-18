#******************************************************************************************

### SETS ###
set nodes;			### Set of nodes/vertexes
set pipes;			### Set of commercial pipes available
set arcs within {i in nodes, j in nodes: i != j};	### Set of arcs/links/edges

#******************************************************************************************

### PARAMETERS ###
param L{arcs};	        ### Total length of each arc/link
param E{nodes};		### Elevation of each node
param P{nodes};		### Minimum pressure required at each node
param D{nodes};		### Demand of each node
param d{pipes};		### Diameter of each commercial pipe
param C{pipes};		### Cost per unit length of each commercial pipe
param R{pipes};		### Roughness of each commercial pipe
param vmax{arcs};
set   Source;		   ### Source node ID
param omega := 10.68;			### SI Unit Constant for Hazen Williams Equation

#******************************************************************************************

### VARIABLES ###
var l{arcs,pipes} >= 0 ;	### Length of each commercial pipe for each arc/link
var q{arcs};	### Flow variable

#******************************************************************************************

### OBJECTIVE ###
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)  * (abs(q[i,j])^(2.852))/2.852);

### CONSTRAINTS ###
s.t. con1{j in nodes}: sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j];

s.t. con4{(i,j) in arcs }: sum{k in pipes} l[i,j,k] = L[i,j];

#******************************************************************************************
