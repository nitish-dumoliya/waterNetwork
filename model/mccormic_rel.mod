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
param vmax{arcs};

param Qmax:= sum {k in nodes diff Source} D[k];
param max_length := max{(i,j) in arcs} L[i,j];

param M := ( (Qmax/1000)^1.852) *  omega * max_length / ( (R[1]^1.852) * (d[1]/1000)^4.87);
#****************************************VARIABLES****************************************#
var l{arcs,pipes} >= 0 ;	# Length of each commercial pipe for each arc/link
var q{arcs} ;	            # Flow variable
var h{nodes};	            # Head
var x{arcs, pipes};	            #
var z{arcs};
var y{arcs};

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
maximize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}: 
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;

subject to con2{(i,j) in arcs}: 
    h[i] - h[j] = (0.001^1.852)*sum{k in pipes } omega * x[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)
;

subject to con_x1 {(i,j) in arcs, k in pipes}:
    x[i,j,k] >= -(Qmax^1.852)*l[i,j,k]
;
subject to con_x4 {(i,j) in arcs, k in pipes}:
    x[i,j,k] >= (Qmax^1.852)*l[i,j,k] + L[i,j]*z[i,j] - L[i,j]*(Qmax)^1.852
;
subject to con_x2 {(i,j) in arcs, k in pipes}:
    x[i,j,k] <= (Qmax^1.852)*l[i,j,k]
;
subject to con_x3 {(i,j) in arcs, k in pipes}:
    x[i,j,k] <= L[i,j]*z[i,j] - (Qmax^1.852)*l[i,j,k] + L[i,j]*(Qmax^1.852)
;

subject to con_z1 {(i,j) in arcs}:
    z[i,j] <= Qmax*y[i,j]
;
subject to con_z2 {(i,j) in arcs}:
    z[i,j] <= -Qmax*y[i,j] + ((Qmax)^0.852)*q[i,j] + Qmax*(Qmax)^0.852
;
subject to con_z3 {(i,j) in arcs}:
    z[i,j] >= -Qmax*y[i,j]
;
subject to con_z4 {(i,j) in arcs}:
    z[i,j] >= ((Qmax)^0.852)*q[i,j] + Qmax*y[i,j] - Qmax*(Qmax)^0.852
;

subject to con_y1 {(i,j) in arcs}:
    y[i,j] <= (Qmax)^0.852
;
subject to con_y2 {(i,j) in arcs}:
    y[i,j] >= ((Qmax)^0.852)*q[i,j]/Qmax
;
subject to con_y3 {(i,j) in arcs}:
    y[i,j] >= ((Qmax)^0.852 )*q[i,j]/(-Qmax)
;

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
subject to con81{(i,j) in arcs}: 
    q[i,j]<= Qmax
;
subject to con82{(i,j) in arcs}: 
    q[i,j]>= -Qmax
;

#******************************Solve the optimization problem*****************************#
# data ../data/d1_Sample_input_cycle_twoloop.dat;
# data ../data/d2_Sample_input_cycle_hanoi.dat;
# data ../Data/d3_Sample_input_double_hanoi.dat;
# data ../Data/d4_Sample_input_triple_hanoi.dat;
# data ../Data/d5_Taichung_input.dat;
# data ../Data/d6_HG_SP_1_4.dat;
# data ../Data/d7_HG_SP_2_3.dat;
# data ../data/d8_HG_SP_3_4.dat;
# data ../Data/d9_HG_SP_4_2.dat;
# data ../Data/d10_HG_SP_5_5.dat;
# data ../Data/d11_HG_SP_6_3.dat;
# data ../data/d12.dat;
# data ../data/d13.dat;
# data ../data/d14_NewYork.dat;
# data ../data/ata/d15_foss_poly_0.dat;
# data ../data/d16_foss_iron.dat;
# data ../data/d17_foss_poly_1.dat;
# data ../data/d19_modena.dat;

# option solver "/home/nitishdumoliya/dist/bin/ipopt";

# option solver "cplexamp";

# option knitro_options "outlev = 1 threads=12 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10 ms_maxtime_real = 10 mip_multistart=1 maxtime_real =60 mip_numthreads=1";

# option baron_options 'maxtime = -1  lsolmsg  outlev = 1';
# option solver "/home/nitishdumoliya/minotaur/build/bin/mqg" ; 
# option mmultistart_options "--presolve 1 --log_level 6  --eval_within_bnds 1 " ;
# option mbnb_options "--presolve 0 --nlp_engine IPOPT --log_level 6 --eval_within_bnds 1 " ;
# solve;

#*************************************Display the results********************************#
# display total_cost;
# display l;
# display q;
# display h;
# display z;
# display x;
# display con1.body, con1.lb, con1.ub;

# for {(i,j) in arcs, k in pipes : l[i,j,k]>1} {
#     printf "Velocity in arc (%s %s %s) : %g\n", i,j,k,  (q[i,j]/1000)/((3.14/4)* (d[k]/1000)**2);
# }
