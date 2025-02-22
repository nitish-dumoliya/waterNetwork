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
param omega := 10.67;  # SI Unit Constant for Hazen Williams Equation
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);
param delta := 0.1;
param p:= 1.852;
param a := (15*(delta)^(p-1))/8 + ((p-1)*p*delta^(p-1))/8 - 7*p*(delta^(p-1))/8;
param b := (-5*(delta)^(p-3))/4 - ((p-1)*p*delta^(p-3))/4 + 5*p*(delta^(p-3))/4; 
param c := (3*(delta)^(p-5))/8 + ((p-1)*p*delta^(p-5))/8 - 3*p*(delta^(p-5))/8;
param Q_max = sum{k in nodes diff Source} D[k]*15850.3;

param eps default 1e-6;  # Small smoothing parameter

#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs};	            # Flow variable
var h{nodes};	            # Head
var t{nodes}>=0;
#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k]/3.28084;	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]*15850.3
;

#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = q[i,j]*abs(q[i,j])^0.852 * sum{k in pipes} (10.44 * l[i,j,k] / ( ((R[k])^1.852) * (d[k])^4.87))
#;

#subject to con2{(i,j) in arcs }:
#    (if -delta<=q[i,j]<=delta  then
#        ((q[i,j] * abs(q[i,j])*(abs(q[i,j])+eps)^0.852) / (abs(q[i,j]) + 0.852*eps)) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87))
#    else 
#		(q[i,j] * abs(q[i,j])^0.852) * sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)) = h[i] - h[j]  
#;


#subject to con2{(i,j) in arcs }: 
#    (if -delta<=q[i,j]<=delta  then 
#        (c*(q[i,j]^5) + b*(q[i,j]^3) + a*q[i,j])*(sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)) 
#    else 
#		(q[i,j] * abs(q[i,j])^0.852) * sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87)) = h[i] - h[j]  
#;

#f_first_order_approx1
#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = (q[i,j]*(abs(q[i,j]) + eps)^0.852 - 0.852*eps*q[i,j]/(abs(q[i,j])+eps)^0.148) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]*39.3701)^4.87));

#f_first_order_approx2
#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = ((q[i,j]*(abs(q[i,j])+0.148*eps)) / (abs(q[i,j]) + eps)^0.148) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));

#f_first_order_approx3
subject to con2{(i,j) in arcs}:
     h[i] - h[j]  = ((q[i,j] * abs(q[i,j])*(abs(q[i,j])+eps)^0.852) / (abs(q[i,j]) + 0.852*eps)) * sum{k in pipes} (10.44 * l[i,j,k] / ( ((R[k])^1.852) * (d[k])^4.87));

#f_second_order_approx
#subject to con2{(i,j) in arcs}:
#     h[i] - h[j]  = (q[i,j]*(abs(q[i,j]) + eps)^0.852 - 0.852*eps*q[i,j]/(abs(q[i,j])+eps)^0.148 + (0.063048*(eps)^2) * q[i,j]/(abs(q[i,j])+eps)^1.148) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));

subject to con3{(i,j) in arcs}: 
    sum{k in pipes} l[i,j,k] = L[i,j]*3.28084
;
subject to con4{(i,j) in arcs , k in pipes}: 
    l[i,j,k] <= L[i,j]*3.28084
;

subject to con5{i in Source}: 
    h[i] = E[i]*3.28084
;

subject to con6_{i in nodes diff Source}: h[i] >= (E[i] + P[i])*3.28084 ;

subject to con7{(i,j) in arcs}: -Q_max <= q[i,j];
subject to con8{(i,j) in arcs}: q[i,j] <= Q_max;

#*******************************************************************************************#
