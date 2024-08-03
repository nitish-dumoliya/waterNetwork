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

param q{arcs};
### VARIABLES ###
var lambda{arcs} ;			### Length of each commercial pipe for each arc/link
var beta{arcs};	### Flow variable
var gamma{nodes};					### Head

### OBJECTIVE ###
maximize total_cost : sum{(i,j) in arcs } L[i,j]*beta[i,j] + E[1]*gamma[1] - sum{j in nodes diff Source} gamma[j]*(E[j] + P[j]);	### Total dual cost 
### CONSTRAINTS ###
s.t. con1{(i,j) in arcs, k in pipes}: -(q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * omega * lambda[i,j] / ( (R[k]^1.852) * (d[k]/1000)^4.87) + beta[i,j] <= C[k];

s.t. con2{j in nodes diff Source}: sum{(j,i) in arcs} lambda[j,i] - sum{(i,j) in arcs} lambda[i,j] = gamma[j];

s.t. con3{s in Source}: sum{(s,j) in arcs}lambda[s,j] + gamma[s] = 0;

s.t. con4{j in nodes diff Source}: gamma[j] <=0;
##########################################################################
