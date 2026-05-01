#***************************SETS***************************#
set pipes;

#************************PARAMETERS************************#
param L;                     # total length for this arc
param y;                     # fixed value (y_ij)

param C{pipes};              # cost
param alpha{pipes};          # alpha_k

#************************VARIABLES************************#
var l{pipes} >= 0;

#************************OBJECTIVE************************#
minimize total_cost:
    sum{k in pipes} C[k] * l[k];

#************************CONSTRAINTS***********************#
subject to FlowConstraint:
    sum{k in pipes} alpha[k] * l[k] = y;

subject to LengthConstraint:
    sum{k in pipes} l[k] = L;
