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
param delta := 0.1;
param p:= 1.852;
param a := (15*(delta)^(p-1))/8 + ((p-1)*p*delta^(p-1))/8 - 7*p*(delta^(p-1))/8;
param b := (-5*(delta)^(p-3))/4 - ((p-1)*p*delta^(p-3))/4 + 5*p*(delta^(p-3))/4; 
param c := (3*(delta)^(p-5))/8 + ((p-1)*p*delta^(p-5))/8 - 3*p*(delta^(p-5))/8;

#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs};	            # Flow variable
var h{nodes};	            # Head
var u{nodes}>=0;            # Lagrangian Multipliers
var x{arcs};                # Lagrangian Multipliers
var v{arcs,pipes};          # Lagrangian Multipliers
#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
# minimize total_cost : (sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k]) + sum{j in nodes diff Source} u[j]*(E[j]+P[j]-h[j]) + sum{(i,j) in arcs , k in pipes} v[i,j,k]*(l[i,j,k]-L[i,j]);
# minimize total_cost : (sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k]) + sum{(i,j) in arcs , k in pipes} v[i,j,k]*(l[i,j,k]-L[i,j]);	
# 
# minimize total_cost : (sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k]) + sum{j in nodes diff Source} u[j]*(E[j]+P[j]-h[j]) ;	
minimize total_cost : (sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k]) + sum{(i,j) in arcs} x[i,j]*(h[i] - h[j] - (q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * sum{k in pipes } omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));	
# minimize total_cost : (sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k]) + sum{(i,j) in arcs} v[i,j]*(sum{k in pipes}l[i,j,k]-L[i,j]);	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes}: 
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;

# subject to con2{(i,j) in arcs}:
#     h[i] - h[j] = (q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * sum{k in pipes } omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)
# ;
# subject to con3{(i,j) in arcs }: 
#                        (if -0.1<=q[i,j]<=0.1  then 
#                            (0.001^1.852)*(c*(q[i,j]^5) + b*(q[i,j]^3) + a*q[i,j])*(sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) 
#                        else 
# 						    (q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) = h[i] - h[j] 
#;

subject to con4{(i,j) in arcs}: 
    sum{k in pipes} l[i,j,k] = L[i,j]
;
subject to con5{(i,j) in arcs , k in pipes}: 
    l[i,j,k] <= L[i,j]
;
subject to con6{i in Source}: 
    h[i] = E[i]
;
# subject to con7{i in nodes diff Source}: 
#     h[i] >= E[i] + P[i]
# ;
subject to con8{(i,j) in arcs}:
    -sum{k in nodes diff Source} D[k] <= q[i,j]
;
subject to con9{(i,j) in arcs}:
    q[i,j] <= sum{k in nodes diff Source} D[k]
;
#*******************************************************************************************#
