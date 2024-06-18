### SETS ###
set nodes;			### Set of nodes/vertexes
set pipes;			### Set of commercial pipes available
set arcs within {i in nodes, j in nodes: i != j};	### Set of arcs/links/edges

### PARAMETERS ###
param L{arcs};		   ### Total length of each arc/link
param E{nodes};		### Elevation of each node
param P{nodes};		### Minimum pressure required at each node
param D{nodes};		### Demand of each node
param d{pipes};		### Diameter of each commercial pipe
param C{pipes};		### Cost per unit length of each commercial pipe
param R{pipes};		### Roughness of each commercial pipe
# param vmax{arcs};
set   Source;		   ### Source node ID
param omega := 10.68;			### SI Unit Constant for Hazen Williams Equation
param delta := 0.1;
param p:= 1.852;
param a := (15*(delta)^(p-1))/8 + ((p-1)*p*delta^(p-1))/8 - 7*p*(delta^(p-1))/8;
param b := (-5*(delta)^(p-3))/4 - ((p-1)*p*delta^(p-3))/4 + 5*p*(delta^(p-3))/4; 
param c := (3*(delta)^(p-5))/8 + ((p-1)*p*delta^(p-5))/8 - 3*p*(delta^(p-5))/8;

### VARIABLES ###
var l{arcs,pipes} >= 0 ;			### Length of each commercial pipe for each arc/link
var q{arcs} ;	### Flow variable
var q1{arcs}>=0 ;	### Flow variable
var q2{arcs}>=0 ;	### Flow variable
var h{nodes};					### Head
# var x{arcs}, binary;
var z{arcs}, binary;
### OBJECTIVE ###
minimize total_cost : sum{(i,j) in arcs  } sum{k in pipes} l[i,j,k]*C[k]*x[i,j];	### Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"

### Variable bounds ###
s.t. bound1{(i,j) in arcs , k in pipes}: l[i,j,k] <= L[i,j];

### CONSTRAINTS ###

s.t. con1{j in nodes }: sum{i in nodes : (i,j) in arcs }(q1[i,j]-q2[i,j]) -  sum{i in nodes : (j,i) in arcs}(q1[j,i]-q2[j,i]) =  D[j];

s.t. con3{(i,j) in arcs}: (h[i] - h[j])*x[i,j] = ((q1[i,j])^1.852 -(q2[i,j])^1.852) * (0.001^1.852) * sum{k in pipes } omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87);

# s.t. con3{(i,j) in arcs }: (if -0.1<=q1[i,j]-q2[i,j]<=0.1  then 
#                               (0.001^1.852)*(c*(q1[i,j]-q2[i,j])^5 + b*(q1[i,j]-q2[i,j])^3 + a*(q1[i,j]-q2[i,j]))*(sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) 
#                            else 
# 						      ((q1[i,j])^1.852 -(q2[i,j])^1.852) * (0.001^1.852) * sum{k in pipes} omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]/1000)^4.87)) = (h[i] - h[j])*x[i,j] ;

s.t. con4{(i,j) in arcs }: sum{k in pipes} l[i,j,k] = L[i,j]*x[i,j];

s.t. con5{i in Source}: h[i] = E[i];

s.t. con2{i in nodes diff Source}: h[i] >=  E[i] + P[i];

# s.t. con6: sum {(i,j) in arcs} x[i,j] = 31;

s.t. con8{(i,j) in arcs }:  q1[i,j] <= -(sum{k in Source} D[k])*z[i,j];
# s.t. con9{(i,j) in arcs }:  q[i,j] >=0;
s.t. con9{(i,j) in arcs }:  q2[i,j] <= -(sum{k in Source} D[k])*(1-z[i,j]);

s.t. flow_con{(i,j) in arcs}: q[i,j] = q1[i,j]-q2[i,j];

# ## Two loop subtour
# s.t. sub_tour_con1: x[2,3]+x[2,4]+x[4,5]+x[3,5] <= 3; 
# s.t. sub_tour_con2: x[4,5]+x[4,6]+x[6,7]+x[7,5] <= 3; 
# s.t. sub_tour_con3: x[2,3]+x[3,5]+x[7,5]+x[6,7]+x[4,6]+x[2,4] <= 5; 

# s.t. x_con1: x[6,7] = 0;
# s.t. x_con2: x[3,5] = 0;
# s.t. z_con1: x[7,5] = 0;

# s.t. degree_con{j in nodes diff Source}: sum{(i,j) in arcs}x[i,j] = 1;
# s.t. degree_con{j in nodes }: sum{(i,j) in arcs}x[i,j] + sum{(j,i) in arcs}x[j,i]  >= 1;

# s.t. sub_tour_con1: x[2,3]+x[2,4]+x[4,5]+x[3,5] = 3; 
# s.t. sub_tour_con2: x[4,6]+x[6,7]+x[7,5] = 2; 

# s.t. sub_tour_con1: x[2,3]+x[2,4]+x[3,5] = 3; 
# s.t. sub_tour_con2: x[4,5]+x[4,6]+x[6,7]+x[7,5] = 2; 

### Hanoi subtour
s.t. sub_tour_con1: x[23,28]+x[23,24]+x[24,25]+x[25,32] + x[32,31]+x[30,31]+x[29,30]+x[28,29] <= 7;
s.t. sub_tour_con2: x[3,20]+x[20,23]+x[23,24]+x[24,25] + x[26,25]+x[27,26]+x[16,27]+x[17,16] + x[18,17] + x[19,18] +x[3,19] <= 10;
s.t. sub_tour_con3: x[3,4]+x[4,5]+x[5,6]+x[6,7] + x[7,8]+x[8,9]+x[9,10]+ x[10,14] + x[14,15] + x[15,16] +x[17,16] + x[18,17] + x[19,18] +x[3,19] <= 13;



### New York Network subtour

# s.t. sub_tour_con1: x[1,2] + x[2,3] + x[3,4] + x[4,5] + x[5,6] + x[6,7] + x[7,8] + x[8,9] + x[11,9] + x[12,11] +  x[13,12] + x[14,13] + x[15,14] + x[1,15] <= 13;
# s.t. sub_tour_con2: x[11,9] + x[11,20] + x[20,16] + x[9,16] <= 3;
