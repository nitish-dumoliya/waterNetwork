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
param R{pipes};		   # Roughness of each commercial pipe
param omega := 10.67;  # SI Unit Constant for Hazen Williams Equation
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);
param p:= 1.852;
param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];
param d_max = max{i in pipes} d[i];
param d_min = min{i in pipes} d[i];

param delta := 0.01;

param R_min = min{k in pipes} R[k];

param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
param q1;
param q2:=Q_max;

param eps{(i,j) in arcs} := 0.0153*(1e-2/MaxK[i,j])^(0.54);
#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs}>=-Q_max, <=Q_max;	            # Flow variable
var h{nodes};	            # Head

var a{arcs};
var b{arcs};
var B1{arcs};
var B2{arcs};
var B3{arcs};
var B4{arcs};
var B5{arcs};

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;

# Smooth-Approximation of Hazen-Williams Constraint
subject to con2{(i,j) in arcs}: 
    h[i] - h[j] = a[i,j]*q[i,j]*abs(q[i,j]) + b[i,j]*q[i,j]
; 

subject to var_a{(i,j) in arcs}:
    a[i,j] = ((sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)))*B4[i,j] - b[i,j] * B3[i,j]) / B1[i,j];

subject to var_b{(i,j) in arcs}:
    b[i,j] = (sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)))*((B5[i,j] * B1[i,j]) - (B4[i,j] * B3[i,j])) / ((B2[i,j]*B1[i,j]) - (B3[i,j]^2));

subject to var_A{(i,j) in arcs}:
    B3[i,j] = (q2^0.296 - q1^0.296) / 0.296;

subject to var_B{(i,j) in arcs}:
    B1[i,j] = (q2^1.296 - q1^1.296) / 1.296;

subject to var_C{(i,j) in arcs}:
    B4[i,j] = (q2^1.148 - q1^1.148) / 1.148;

subject to var_D{(i,j) in arcs}:
    B2[i,j] = (q2^(-0.704) - q1^(-0.704)) / (-0.704);

subject to var_E{(i,j) in arcs}:
    B5[i,j] = (q2^0.148 - q1^0.148) / 0.148;

subject to con3{(i,j) in arcs}: 
    sum{k in pipes} l[i,j,k] = L[i,j]
;
subject to con4{(i,j) in arcs , k in pipes}: 
    l[i,j,k] <= L[i,j]
;
subject to con5{(i,j) in arcs , k in pipes}: 
    l[i,j,k] >= 0
;
subject to con6{i in Source}: 
    h[i] = E[i]
;
subject to con7{i in nodes diff Source}: h[i] >= (E[i] + P[i]) ;
#*******************************************************************************************#
