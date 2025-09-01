#-----------------------------------------SETS--------------------------------------------#
set nodes;
set pipes;
set arcs within {i in nodes, j in nodes: i != j};
set Source;

#--------------------------------------PARAMETERS-----------------------------------------#
param L{arcs };
param E{nodes};
param P{nodes};
param pmax{nodes};
param D{nodes};
param d{pipes}; 
param C{pipes};
param R{pipes};
param omega := 10.67;
param vmax{arcs} default (sum {k in nodes diff Source} D[k])/((3.14/4)*(d[1])^2);
param delta := 0.1;
param p:= 1.852;
param Q_max = sum{k in nodes diff Source} D[k];
param D_min = min{i in nodes diff Source} D[i];
param D_max = max{i in nodes diff Source} D[i];
param d_max = max{i in pipes} d[i];
param d_min = min{i in pipes} d[i];

param R_min = min{k in pipes} R[k];
param eps{(i,j) in arcs} := (5.35*1e-5)^2;
param MaxK{(i,j) in arcs} := omega * L[i,j] / (R_min^1.852 * d_min^4.87);
#---------------------------------------VARIABLES----------------------------------------#
var l{arcs,pipes} >= 0 ;
var q{arcs}>=-Q_max, <=Q_max;
var h{nodes};

var q1{arcs}>=0,<=Q_max ;	
var q2{arcs}>=0,<=Q_max ;
#var hl{arcs}>=0 ;
var hl1{arcs}>=0 ;
var hl2{arcs}>=0 ;
var z{arcs}<=1,>=-1;
#var z2{arcs} binary ;

#---------------------------------------OBJECTIVE----------------------------------------#
minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	

#--------------------------------------CONSTRAINTS---------------------------------------#
subject to con1{j in nodes diff Source}:
     sum{i in nodes : (i,j) in arcs }(q1[i,j]-q2[i,j]) -  sum{i in nodes : (j,i) in arcs}(q1[j,i]-q2[j,i]) =  D[j]
;
#subject to con2_positive{(i,j) in arcs}: 
#     hl1[i,j]  = q1[i,j]^1.852 * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87))
#;
#subject to con2_negative{(i,j) in arcs}: 
#     hl2[i,j]  = q2[i,j]^1.852 * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87))
#;
subject to con2_positive{(i,j) in arcs}: 
     hl1[i,j]  = (q1[i,j])^3 * ((((q1[i,j])^2 + eps[i,j])^0.426) /((q1[i,j])^2 + 0.426*eps[i,j])) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87))
;
subject to con2_negative{(i,j) in arcs}: 
     hl2[i,j]  = (q2[i,j])^3 * ((((q2[i,j])^2 + eps[i,j])^0.426) /((q2[i,j])^2 + 0.426*eps[i,j])) * sum{k in pipes} (omega * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87))
;

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

subject to con7{(i,j) in arcs}: q1[i,j] <= Q_max*(1+z[i,j]);
subject to con8{(i,j) in arcs}: q2[i,j] <= Q_max*(1-z[i,j]);

subject to con9{(i,j) in arcs}: q[i,j] = q1[i,j]-q2[i,j];
#subject to con10{(i,j) in arcs}: z[i,j]^2 + (1-z[i,j])^2 <= 1;
#subject to con10{(i,j) in arcs}: (1+z[i,j])*(1-z[i,j]) = 0;
#subject to con10{(i,j) in arcs}: abs(z[i,j]) = 1;

subject to con11{(i,j) in arcs}: h[i] - h[j] = (hl1[i,j] - hl2[i,j]) ;
subject to con12{(i,j) in arcs}: hl1[i,j] <= (Q_max^1.852) * MaxK[i,j] * (1+z[i,j]) ;
subject to con13{(i,j) in arcs}: hl2[i,j] <= (Q_max^1.852) * MaxK[i,j] * (1-z[i,j]) ;
#-----------------------------------------------------------------------------------------#
