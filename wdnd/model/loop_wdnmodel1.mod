# ----------------- SETS -----------------
set nodes;
set pipes;
set arcs within {i in nodes, j in nodes: i != j};
set Source;                # should contain exactly one source node for this model

# Topology sets produced externally (use the NetworkX script below to create):
set tree_arcs within arcs;  # spanning tree arcs (undirected orientation must match arcs order)
set chords within arcs;     # arcs not in the spanning tree
set loops;                  # one loop id per chord

# loop_sign[loop,(i,j)] = -1,0,1
param loop_sign{loops, arcs} integer, default 0;

# path_sign[node,(i,j)] = -1,0,1
# path_sign gives the contribution of arc (i,j) to the path from source to node:
# +1 if arc appears with same orientation as path, -1 if reversed, 0 if not on path.
param path_sign{nodes, arcs} integer, default 0;

# ----------------- PARAMETERS -----------------
param L{arcs};              # total length of each arc/link (input)
param E{nodes};             # elevation / source head reference values
param P{nodes};             # minimum pressure requirement (head - elevation >= P)
param pmax{nodes};          # not used in this model; keep if needed
param D{nodes};             # demand (positive outflow) at each node
param d{pipes};
param C{pipes};
param R{pipes};
param omega := 10.67;
param p := 1.852;
param Q_max default sum{k in nodes diff Source} D[k];

# eps must be > 0 to keep smooth denominator safe
param eps{arcs} > 0;

# ----------------- VARIABLES -----------------
var l{arcs, pipes} >= 0;    # length of chosen commercial pipe in each arc
var q{arcs};                # signed flow on arcs (no heads variable)

# ----------------- AUXILIARY PARAM/EXPRESSIONS -----------------
# Kcoeff per arc (resistance coefficient aggregated over chosen pipes)
param Kcoeff{(i,j) in arcs} :=
    sum{k in pipes} ( omega * l[i,j,k] / ( (R[k]^1.852) * (d[k]^4.87) ) );

# Smooth Hazen-Williams head-loss Phi for arc (i,j) as an expression in q:
# Phi_ij(q) = Kcoeff[i,j] * ( q^3 * (q^2 + eps^2)^0.426 / ( q^2 + 0.426*eps^2 ) )
# We'll use an AMPL indexed expression via 'let' style (AMPL supports expressions in constraints directly).

# ----------------- OBJECTIVE -----------------
minimize total_cost:
    sum{(i,j) in arcs} sum{k in pipes} l[i,j,k] * C[k];

# ----------------- CONSTRAINTS -----------------
# Mass balance at demand nodes (sources excluded or include but D[source]=...? Here exclude)
subject to mass_balance{j in nodes diff Source}:
    sum{i in nodes: (i,j) in arcs} q[i,j] - sum{i in nodes: (j,i) in arcs} q[j,i] = D[j];

# Length partition constraints (same as before)
subject to total_length{(i,j) in arcs}:
    sum{k in pipes} l[i,j,k] = L[i,j];

subject to length_upper{(i,j) in arcs, k in pipes}:
    l[i,j,k] <= L[i,j];

subject to length_nonneg{(i,j) in arcs, k in pipes}:
    l[i,j,k] >= 0;

# Flow bounds
subject to q_bounds{(i,j) in arcs}:
    -Q_max <= q[i,j] <= Q_max;

# ----------------- LOOP ENERGY EQUATIONS -----------------
# For each loop ℓ: sum_{arcs} loop_sign[ℓ,(i,j)] * Phi_ij(q[i,j]) = 0
subject to loop_energy{ℓ in loops}:
    sum{(i,j) in arcs} loop_sign[ℓ,(i,j)]
        * ( (q[i,j]^3) * ((q[i,j]^2 + eps[i,j]^2)^0.426)
            / ( q[i,j]^2 + 0.426 * eps[i,j]^2 ) ) * Kcoeff[i,j]
    = 0;

# ----------------- NODE HEAD COMPUTATION & MINIMUM PRESSURE -----------------
# For each node j, compute its head from the source head and the tree-path head-losses:
# H_j = E[source] - sum_{(i,k) in arcs} path_sign[j,(i,k)] * Phi_{i,k}(q[i,k])
# Enforce H_j >= E[j] + P[j]  (minimum required head is elevation + pressure head P[j])

# We'll treat Source as a single-element set; pick the unique source s:
param s default 0;
# user must set s to the source node id in the data or ensure Source has the chosen source
# For safety, allow using the only element of Source if s==0: (some AMPL data handling needed)
# Here we assume s is provided in the .dat (set Source has one element and you set param s := that; )

subject to min_pressure{j in nodes diff {s}}:
    ( E[s]
      - sum{(i,k) in arcs} path_sign[j,(i,k)]
            * ( (q[i,k]^3) * ((q[i,k]^2 + eps[i,k]^2)^0.426)
                / ( q[i,k]^2 + 0.426 * eps[i,k]^2 ) ) * Kcoeff[i,k]
    ) >= ( E[j] + P[j] );

# If you want to also ensure source head equals elevation E[s], that is implicit since we use E[s] as reference.
# End model

