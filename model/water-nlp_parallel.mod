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
param vmax{arcs} default (sum {k in nodes diff Source} D[k]/1000)/((3.14/4)*(d[1]/1000)^2);
param delta := 0.1;
param p:= 1.852;
param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];
param d_max = max{i in pipes} d[i]/1000;
param d_min = min{i in pipes} d[i]/1000;


#param eps default 1e-4;  # Small smoothing parameter

#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs};	            # Flow variable
var h{nodes};	            # Head
var t{nodes}>=0;

var eps{arcs}>=0;   # Small smoothing variable
var q1{arcs}>=0 ;	### Flow variable
var q2{arcs}>=0 ;	### Flow variable

var z{arcs}, binary;

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;

#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = q[i,j]*abs(q[i,j])^0.852 * (0.001^1.852) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));

#s.t. con2{(i,j) in arcs}: (h[i] - h[j]) = ((q1[i,j])^1.852 -(q2[i,j])^1.852) * (1000^3.018) * sum{k in pipes } omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87);

subject to con2{(i,j) in arcs}: 
     h[i] - h[j]  = (((q1[i,j])^2 *(((q1[i,j] + 1e+3 * eps[i,j] )^0.852) /(q1[i,j] + 0.852 * 1e+3 *eps[i,j]))) - ((q2[i,j])^2 *(((q2[i,j] + 1e+3 * eps[i,j] )^0.852) /(q2[i,j] + 0.852 * 1e+3 *eps[i,j]))) ) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));

#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = (q[i,j])^3 *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) /((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) * (0.001^1.852)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87));


#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = ((q[i,j])^3/((q[i,j])^2 + 0.426 * 1e+6 *eps[i,j])) *((((q[i,j])^2 + 1e+6 * eps[i,j] )^0.426) ) * (1000^3.018)  * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));


subject to con3{(i,j) in arcs}: 
    sum{k in pipes} l[i,j,k] = L[i,j]
;

subject to con4{(i,j) in arcs , k in pipes}: 
    l[i,j,k] <= L[i,j]
;

subject to con5{i in Source}: 
    h[i] = E[i]
;

subject to con6{i in nodes diff Source}: h[i] >= E[i] + P[i] ;
#subject to con6_{i in nodes diff Source}: h[i] <= sum{j in Source} E[j] ;

subject to con7{(i,j) in arcs}: -Q_max <= q[i,j];
subject to con8{(i,j) in arcs}: q[i,j] <= Q_max;

s.t. con9{(i,j) in arcs }:  q1[i,j] <= Q_max*z[i,j];
# s.t. con9{(i,j) in arcs }:  q[i,j] >=0;
s.t. con10{(i,j) in arcs }:  q2[i,j] <= Q_max*(1-z[i,j]);

s.t. con11{(i,j) in arcs}: q[i,j] = q1[i,j]-q2[i,j];



#subject to eps_selection{(i,j) in arcs}: 0.07508*(sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87))) * eps[i,j]^0.926 = 1e-6;
subject to eps_selection{(i,j) in arcs}: 0.04001571*((omega/( (R[1]^1.852) * (d_min)^4.87))) * eps[i,j]^1.852 = 1e-6;
#subject to eps_selection{(i,j) in arcs}: eps[i,j] = (1e-6/(0.07508*(sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)))))^(1/0.926);

#subject to eps_selection{(i,j) in arcs}: eps[i,j] <= 1e-8;
#*******************************************************************************************#
