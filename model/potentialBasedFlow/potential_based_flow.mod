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
var h{nodes};			### Head
var z1{arcs} binary;
var z2{arcs} binary;

#******************************************************************************************

### OBJECTIVE ###
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes} l[i,j,k]*C[k];	### Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
# minimize total_cost : sum{(i,j) in arcs} sum{k in pipes} l[i,j,k]*(z1[i,j] + z2[i,j])*C[k];	### Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"

### Variable bounds ###
s.t. bound1{(i,j) in arcs , k in pipes}: l[i,j,k] <= L[i,j];

### CONSTRAINTS ###
s.t. con1{j in nodes diff Source }: sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j];

s.t. con3{(i,j) in arcs}: h[i] - h[j] = ((z1[i,j]-z2[i,j])* abs(q[i,j])^1.852) * (0.001^1.852) * sum{k in pipes } omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87);

s.t. con4{(i,j) in arcs }: sum{k in pipes} l[i,j,k] = L[i,j];

s.t. con5{i in Source}: h[i] = E[i];

s.t. con2{i in nodes diff Source}: h[i] >= E[i] + P[i];

s.t. con6{(i,j) in arcs }: z1[i,j] + z2[i,j] = 1;

s.t. con8{(i,j) in arcs }: q[i,j] <= (sum{k in nodes diff Source} D[k])*z1[i,j];

s.t. con9{(i,j) in arcs }: -(sum{k in nodes diff Source} D[k])*z2[i,j] <= q[i,j] ;

s.t. q12: q[1,2] = sum{k in nodes diff Source} D[k] ;

#******************************************************************************************

# data ../Data/d1_Sample_input_cycle_twoloop.dat;
# data ../Data/d2_Sample_input_cycle_hanoi.dat;

# Solve the optimization problem
# option solver baron;
# option ipopt_options "outlev = 3";
# option baron_options 'maxtime = -1 lsolmsg  outlev = 1';
# option knitro_options 'outlev = 4 ms_enable = 1 ms_maxsolves = 5 mip_multistart = 1';

# option solver "/home/nitishdumoliya/Nitish/minotaur/build/bin/mmultistart" ; 
# option mmultistart_options "--presolve 1 --log_level 6 --eval_within_bnds 1 " ;
# option solver cplexamp;
# option presolve 1;
# option presolve_eps 0.001;
# solve;

# display total_cost;
# display l;
# display q;
# display h;
# display z1;
# display z2;

# var avg_dia;
# for {(i,j) in arcs} {
#     # let avg_dia := (sum{k in pipes} (d[k]*l[i,j,k])/L[i,j]);

#     let avg_dia := (L[i,j] / sum{k in pipes} l[i,j,k]/(d[k]**4.87) )**(1/4.87);

#     printf "Velocity in arc (%s %s) : %g\n", i,j,  q[i,j]/((3.14/4)* avg_dia**2);
# }

#******************************************************************************************

