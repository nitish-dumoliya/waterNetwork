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
#param eps{arcs};

#param eps{(i,j) in arcs} := (1e-5 / (0.07508 * MaxK[i,j]))^(1/1.852);
#param eps{(i,j) in arcs} := 0.0953*(1e-2/MaxK[i,j])^(0.54);
param eps{(i,j) in arcs} := 0.0535*(1e-3/MaxK[i,j])^(0.54);

#param eps{(i,j) in arcs} := (1e-6 / (1.267 * MaxK[i,j]))^(1 / 1.852);

#param eps{(i,j) in arcs} := (1e-6 / (0.36061 * MaxK[i,j]))^(0.54);
#param eps{(i,j) in arcs} := 0.0153*(1e-2/MaxK[i,j])^(0.54);

#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs};	            # Flow variable
var h{nodes};	            # Head
#var eps{arcs}>=1e-12, <=1;
#var x{arcs};
#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	
#minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k] - sum {(i,j) in arcs} abs(h[i] - h[j]);	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;

# hazen-Williams Constraint 
#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = q[i,j]*abs(q[i,j])^0.852 * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));

# Smooth-Approximation of Hazen-Williams Constraint
subject to con2{(i,j) in arcs}: 
     #h[i] - h[j]  = q[i,j]*abs(q[i,j])^0.852 * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));
    h[i] - h[j]  =  (q[i,j]^3 * (q[i,j]^2 + eps[i,j]^2)^0.426 / (q[i,j]^2 + 0.426*eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #(h[i] - h[j])  =  q[i,j] * ((q[i,j]^2 + eps[i,j]^2))^0.426 * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #h[i] - h[j]  =  (q[i,j] * (q[i,j]^2 + 0.574 * eps[i,j]^2) / (q[i,j]^2 + eps[i,j]^2)^0.574) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #h[i] - h[j]  =  (q[i,j] * (q[i,j]^2 + eps[i,j]^2)^0.426 - (0.426 * eps[i,j]^2 / (q[i,j]^2 + eps[i,j]^2)^0.574)) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #(h[i] - h[j])*(q[i,j]^2 + eps[i,j]^2)^0.574  =  (q[i,j] * (q[i,j]^2 + eps[i,j]^2)) * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));
    #(h[i] - h[j])  =  q[i,j] * (x[i,j])^0.426 * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87));

#subject to var_x{(i,j) in arcs}: 
#    x[i,j] = q[i,j]^2 + eps[i,j]^2
#;
#subject to x_bounds_left{(i,j) in arcs}: 
#    eps[i,j]^2 <= x[i,j] <= Q_max^2 + eps[i,j]^2
#;

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

#subject to con3{(i,j) in arcs}: 
#    sum{k in pipes} l[i,j,k] = L[i,j]
#;
subject to con4{(i,j) in arcs , k in pipes}: 
    l[i,j,k] <= L[i,j]
;
#subject to con5{(i,j) in arcs , k in pipes}: 
#    l[i,j,k] >= 0
#;
subject to con6{i in Source}: 
    h[i] = E[i]
;
subject to con7{i in nodes diff Source}: h[i] >= (E[i] + P[i]) ;
#subject to con6_{i in nodes diff Source}: h[i] <= sum{j in Source} E[j] ;
subject to con8{(i,j) in arcs}: 
   -Q_max <= q[i,j] <= Q_max
;
#*******************************************************************************************#
