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
#param eps{(i,j) in arcs} := 0.0535*(1e-4/MaxK[i,j])^(0.54);
param eps{arcs};
#****************************************VARIABLES****************************************#
param l{arcs,pipes};	# Length of each commercial pipe for each arc/link
param q{arcs};	            # Flow variable
param h{nodes};	            # Head

var lam{nodes};
var x{arcs};
var y{arcs};
var u{Source};
var v{nodes diff Source}>=0;
var w{arcs, pipes}>=0;
var s{arcs, pipes}>=0;

#****************************************OBJECTIVE****************************************#
# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"
maximize total_cost : sum{(i,j) in arcs} sum{k in pipes} l[i,j,k]*C[k] + 
		      sum{(i,j) in arcs} q[i,j]*(lam[j]- lam[i]) - sum{j in nodes} lam[j]*D[j] + 
		      #sum{j in nodes diff Source} lam[j] * (sum{i in nodes : (i,j) in arcs }q[i,j] -  sum{i in nodes : (j,i) in arcs}q[j,i] -  D[j]) + 
		      sum{(i,j) in arcs} x[i,j] * (h[i] - h[j] - q[i,j]*abs(q[i,j])^0.852 * sum{k in pipes}(omega * l[i,j,k] / (R[k]^1.852 * d[k]^4.87))) + 
		      sum{(i,j) in arcs} y[i,j] * (sum{k in pipes} l[i,j,k] - L[i,j]) + 
		      sum{i in Source} u[i]*(h[i] - E[i]) + 
		      sum{j in nodes diff Source} v[j]*(-h[j]+E[j]+P[j])+ 
		      sum{(i,j) in arcs} sum{k in pipes} (-l[i,j,k]*w[i,j,k] + (l[i,j,k]-L[i,j])*s[i,j,k]) 
		      ;
#*******************************************************************************************#
