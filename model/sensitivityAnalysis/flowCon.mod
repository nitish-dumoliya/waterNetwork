# SETS 
set nodes;			# Set of nodes/vertexes
set pipes;			# Set of commercial pipes available
set arcs within {i in nodes, j in nodes: i != j};	# Set of arcs/links/edges
set Source;		# Source node ID

# PARAMETERS 
param L{arcs};		# Total length of each arc/link
param E{nodes};		# Elevation of each node
param P{nodes};		# Minimum pressure required at each node
param D{nodes};		# Demand of each node
param d{pipes};		# Diameter of each commercial pipe
param C{pipes};		# Cost per unit length of each commercial pipe
param R{pipes};		# Roughness of each commercial pipe
param vmax{arcs} default (sum {k in nodes diff Source} D[k]/1000)/((3.14/4)*(d[1]/1000)^2);
param q_lp{arcs};	# Fixed Flow from nlp solution 

# VARIABLES 
var q_del{arcs};	# Flow variable

# OBJECTIVE 
minimize total_cost: 0;	# Total cost as a sum of "length of the commercial pipe * cost per unit length of the commercial pipe"

# CONSTRAINTS 
subject to con1{j in nodes}: sum{i in nodes: (i,j) in arcs} (q_lp[i,j]+q_del[i,j]) - sum{i in nodes: (j,i) in arcs} (q_lp[j,i]+q_del[j,i]) = D[j];

###################################################################
