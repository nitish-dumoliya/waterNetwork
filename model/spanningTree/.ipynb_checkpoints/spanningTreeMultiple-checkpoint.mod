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

param Qmax:= sum{k in nodes diff Source} D[k];
#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs} ;	            # Flow variable
var h{nodes};	            # Head
var x{arcs}, binary;	            #
# var z{arcs, pipes}, binary;

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
minimize total_cost : sum{(i,j) in arcs} (sum{k in pipes}l[i,j,k]*C[k])*x[i,j];	
#minimize total_cost : sum{i in nodes} h[i];	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}: 
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;

# subject to demand_con1 {(i,j) in arcs }:
#     -(Qmax - D[j])<= q[i,j] <= (Qmax-D[i])
# ;

subject to con2{(i,j) in arcs}: 
    (h[i] - h[j])*x[i,j] = (q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * sum{k in pipes } omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)
;

# subject to con3{(i,j) in arcs}: 
#     (if -0.1<=q[i,j]<=0.1  then 
#         (0.001^1.852)*(q[i,j])*(sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) 
#      else  
#  		(q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) = (h[i] - h[j])*x[i,j] 
# ;

#subject to con3{(i,j) in arcs}: 
#    (if -0.1<=q[i,j]<=0.1  then 
#        (0.001^1.852)*(c*(q[i,j]^5) + b*(q[i,j]^3) + a*q[i,j])*(sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) 
#     else  
# 		(q[i,j] * abs(q[i,j])^0.852) * (0.001^1.852) * sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) = x[i,j]*(h[i] - h[j]) 
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

subject to con7{i in nodes diff Source}: 
    h[i] >= E[i] + P[i]
;

subject to flow_con{(i,j) in arcs }:
     q[i,j] <= Qmax*x[i,j]
;

subject to flow_con1{(i,j) in arcs }:
    -Qmax*(x[i,j]) <= q[i,j] 
;

# subject to con8{(i,j) in arcs, k in pipes}: 
#     (q[i,j]/1000)/((3.14/4)*(d[k]/1000)^2) <= vmax[i,j] + 1000*(1-z[i,j,k])
# ;
# subject to con9{(i,j) in arcs, k in pipes}: 
#     z[i,j,k] <= L[i,j] * l[i,j,k]
# ;
# subject to con10{(i,j) in arcs, k in pipes}: 
#     z[i,j,k] >= l[i,j,k] / L[i,j]
# ;

# subject to con9{(i,j) in arcs, k in pipes}: 
#     l[i,j,k] <= L[i,j] * z[i,j,k]
# ;

#******************************Solve the optimization problem*****************************#
# data ../Data/d1_Sample_input_cycle_twoloop.dat;
# data ../Data/d2_Sample_input_cycle_hanoi.dat;
# data ../Data/d3_Sample_input_double_hanoi.dat;
# data ../Data/d4_Sample_input_triple_hanoi.dat;
# data ../Data/d5_Taichung_input.dat;
# data ../Data/d6_HG_SP_1_4.dat;
# data ../Data/d7_HG_SP_2_3.dat;
# data ../Data/d8_HG_SP_3_4.dat;
# data ../Data/d9_HG_SP_4_2.dat;
# data ../Data/d10_HG_SP_5_5.dat;
# data ../Data/d11_HG_SP_6_3.dat;
# data ../Data/d12.dat;
# data ../Data/d13.dat;
# data ../Data/d14_NewYork.dat;
# data ../Data/d15_foss_poly_0.dat;
# data ../Data/d16_foss_iron.dat;
# data ../Data/d17_foss_poly_1.dat;
# data ../Data/d19_modena.dat;

# option solver "ipopt";
# option solver "baron";
# option knitro_options "outlev = 1 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 30 ms_maxtime_real = 10 mip_multistart=1 maxtime_real =60";
# option baron_options 'maxtime = -1  outlev = 1 lsolver=conopt  barstats objbound  epsr = 0.05  threads = 12  prloc = 1 prfreq=1000 prtime 10';
# option solver "/home/nitishdumoliya/minotaur/build/bin/mmultistart" ; 
# option mmultistart_options "--presolve 1 --log_level 3 --eval_within_bnds 1 " ;
# solve;

#*************************************Display the results********************************#
# display total_cost;
# display l;
# display q;
# display h;
# display x;

# for {(i,j) in arcs, k in pipes : l[i,j,k]>1} {
#     printf "Velocity in arc (%s %s %s) : %g\n", i,j,k,  (q[i,j]/1000)/((3.14/4)* (d[k]/1000)**2);
# }

