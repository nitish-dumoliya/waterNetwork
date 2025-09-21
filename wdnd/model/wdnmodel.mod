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
param a := (15/8)*delta**(p-1) + (1/8)*(p-1)*p*delta**(p-1) - (7/8)*p*delta**(p-1);
param b := - (5/4)*delta**(p-3) - (1/4)*(p-1)*p*delta**(p-3) - (5/4)*p*delta**(p-3);
param c := (3/8)*delta**(p-5) + (1/8)*(p-1)*p*delta**(p-5) - (3/8)*p*delta**(p-5);

param R_min = min{k in pipes} R[k];

param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);

#param eps{(i,j) in arcs} := (1e-6 / (0.07508 * MaxK[i,j]))^(1 / 0.926);
#param eps{(i,j) in arcs} := 4.047*(1e-4)^(1/1.852)*1e-4;
param eps{(i,j) in arcs} := (0.0535/(MaxK[i,j])^(0.54)) * (1e-1)^(0.54);
#param eps{(i,j) in arcs} := (1e-3 / (0.04001571 * MaxK[i,j]))^(1 / 1.852);
#param eps{(i,j) in arcs} := (1e-4 * R_min^1.852 * d_min^4.87 / (0.07508 * 10.67 * L[i,j]))^(1 / 0.926);
#param eps{(i,j) in arcs} := (1/0.58023)*(1.54e-8 / (MaxK[i,j]*0.0000100002395709))^(1 / 1.852);
#param eps{(i,j) in arcs} := 0.3074/(MaxK[i,j]^0.54) *(0.01^0.54);

#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs}>=-Q_max, <=Q_max;	            # Flow variable
var h{nodes};	            # Head
#var eps{arcs}>=1e-11, <=1;

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;

# hazen-Williams Constraint 
#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = q[i,j]*abs(q[i,j])^0.852 * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));

#subject to epsilon_upper1{(i,j) in arcs}:
    #eps[i,j] * (1e-2 + sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87)))^(0.54)  <= 0.0535 * (abs(h[i]-h[j]) + 1e-4)^(0.54);
    #eps[i,j]  <= 0.0535 / (1e-6 + sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87)))^(0.54) * (abs(h[i]-h[j]) + 1e-6)^(0.54);
    #eps[i,j] * (sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87)))^(0.54)  <= 0.0535 * (abs(h[i]-h[j]) + 1e-2)^(0.54);
    #eps[i,j] = 0.0535*(abs(q[i,j])+1e-4);
#    eps[i,j] <= (0.0535 / (MaxK[i,j]^0.54) ) * (1e-3)^(0.54);

# Smooth-Approximation of Hazen-Williams Constraint
subject to con2{(i,j) in arcs}: 
    h[i] - h[j]  =  (q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426*eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));

#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j] = (q[i,j] * ((q[i,j]^2 + eps[i,j])^0.426) * ( (q[i,j]^2 + 0.713 * eps[i,j]) / (q[i,j]^2 + 0.287 * eps[i,j]) )) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));

#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = (q[i,j]*abs(q[i,j])) *(((abs(q[i,j]) + eps[i,j])^0.852) /(abs(q[i,j]) + 0.852*eps[i,j]))  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));

#subject to con2{(i,j) in arcs }: 
#    (if -delta<=q[i,j]<=delta  then  
#        (q[i,j])^3 *((((q[i,j])^2 + eps[i,j])^0.426) /((q[i,j])^2 + 0.426*eps[i,j]))  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87))
#    else 
#		(q[i,j] * abs(q[i,j])^0.852) * sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)) = h[i] - h[j]  
#;

# Bragalli approximation
#subject to con2{(i,j) in arcs }: 
#    (if -delta<=q[i,j]<=delta  then 
#        (c*(q[i,j]^5) + b*(q[i,j]^3) + a*q[i,j])*(sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)) 
#    else 
#		(q[i,j] * abs(q[i,j])^0.852) * sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)) = h[i] - h[j]  
#;

subject to con3{(i,j) in arcs}: 
    sum{k in pipes} l[i,j,k] = L[i,j]
;

subject to con4{(i,j) in arcs , k in pipes}: 
    l[i,j,k] <= L[i,j]
;

subject to con5{i in Source}: 
    h[i] = E[i]
;

subject to con6{i in nodes diff Source}: h[i] >= (E[i] + P[i]) ;
#subject to con6_{i in nodes diff Source}: h[i] <= sum{j in Source} E[j] ;

#subject to con7{(i,j) in arcs}: -Q_max <= q[i,j];
#subject to con8{(i,j) in arcs}: q[i,j] <= Q_max;

#*******************************************************************************************#
