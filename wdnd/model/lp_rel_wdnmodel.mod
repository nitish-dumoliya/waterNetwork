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

param eps{(i,j) in arcs} := (1e-6 / (0.07508 * MaxK[i,j]))^(1 / 0.926);
#param eps{(i,j) in arcs} := (1e-3 / (0.04001571 * MaxK[i,j]))^(1 / 1.852);
#param eps{(i,j) in arcs} := (1e-4 * R_min^1.852 * d_min^4.87 / (0.07508 * 10.67 * L[i,j]))^(1 / 0.926);

param n_exp{(i,j) in arcs} := 
    (6 + log10(0.07508 * MaxK[i,j])) / 0.926;

#param eps{(i,j) in arcs} := 10^(-ceil(n_exp[i,j]));
#****************************************VARIABLES****************************************#
param l{arcs,pipes};	# Length of each commercial pipe for each arc/link
var q{arcs};	            # Flow variable
var h{nodes};	            # Head

var w{arcs, pipes};
var x{arcs};
var y{arcs};


#var eps{arcs};
#var x{arcs, pipes}>=0,<=1;
#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
#minimize total_cost : sum{(i,j) in arcs} sum{k in pipes}l[i,j,k]*C[k];	
minimize total_cost : 0;	

#****************************************CONSTRAINTS**************************************#
subject to con1{j in nodes diff Source}:
    sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] =  D[j]
;

# hazen-Williams Constraint 
#subject to con2{(i,j) in arcs}: 
#     h[i] - h[j]  = sum{k in pipes} (omega * w[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));

subject to con2{(i,j) in arcs}: 
     h[i] - h[j]  = sum{k in pipes} (omega * x[i,j] * l[i,j,k] / ( (R[k]^1.852) * (d[k])^4.87));


#subject to con3{(i,j) in arcs}: 
#    sum{k in pipes} l[i,j,k] = L[i,j]
#;

# 4a
#subject to con_4a {(i,j) in arcs, k in pipes}:
#    w[i,j,k] >= - (Q_max^1.852) * l[i,j,k];

# 4b
#subject to con_4b {(i,j) in arcs, k in pipes}:
#    w[i,j,k] >= (Q_max^1.852) * l[i,j,k] + x[i,j] * L[i,j] - (Q_max^1.852) * L[i,j];

# 4c
#subject to con_4c {(i,j) in arcs, k in pipes}:
#    w[i,j,k] <= (Q_max^1.852) * l[i,j,k];

# 4d
#subject to con_4d {(i,j) in arcs, k in pipes}:
#    w[i,j,k] <= x[i,j] * L[i,j] - (Q_max^1.852) * l[i,j,k] + (Q_max^1.852) * L[i,j];


# 5a
subject to con_5a {(i,j) in arcs}:
    x[i,j] >= -Q_max * y[i,j];

# 5b
subject to con_5b {(i,j) in arcs}:
    x[i,j] >= Q_max * y[i,j] + q[i,j] * (Q_max^0.852) - (Q_max^1.852);

# 5c
subject to con_5c {(i,j) in arcs}:
    x[i,j] <= Q_max * y[i,j];

# 5d
subject to con_5d {(i,j) in arcs}:
    x[i,j] <= q[i,j] * (Q_max^0.852);


# 6a
subject to con_6a {(i,j) in arcs}:
    y[i,j] >= ((Q_max^0.852) / Q_max) * q[i,j];

# 6b
subject to con_6b {(i,j) in arcs}:
    y[i,j] >= - ((Q_max^0.852) / Q_max) * q[i,j];

# 6c
subject to con_6c {(i,j) in arcs}:
    y[i,j] <= (Q_max^0.852);

#subject to con7{(i,j) in arcs , k in pipes}: 
#    l[i,j,k] <= L[i,j]
#;

subject to con8{i in Source}: 
    h[i] = E[i]
;

subject to con9{i in nodes diff Source}: h[i] >= (E[i] + P[i]) ;

subject to con10a{(i,j) in arcs}: -Q_max <= q[i,j];
subject to con10b{(i,j) in arcs}: q[i,j] <= Q_max;

#*******************************************************************************************#
