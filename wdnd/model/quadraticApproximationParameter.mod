# -------------------------------
# Parameters
# -------------------------------
param q2 := 0.42043;
param eps_tol:=1e-1;
param alpha := 199592.673693674;

# -------------------------------
# Variables
# -------------------------------
var a;
var b;
var q1>=0;
var A;
var B;
var C;
var D;
var E;
var q_eps;

# -------------------------------
# Objective Function
# -------------------------------
minimize obj: 0;

# -------------------------------
# Constraints
# -------------------------------
subject to root:
    ((a * q_eps^2 + b * q_eps) - alpha * (q_eps^1.852)) / (alpha * q_eps^1.852) - eps_tol = 0;

subject to var_q_eps:
    q_eps = 0.85 * b / (0.15 * a);

subject to var_a:
    a = (C - b * A) / B;

subject to var_b:
    b = ((A * C) / (D * B) - (E / D)) / (1 + ((A^2) / (D * B)));

subject to var_A:
    A = (q2^0.3 - q1^0.3) / (0.3 * alpha^2);

subject to var_B:
    B = (q2^1.3 - q1^1.3) / (1.3 * alpha^2);

subject to var_C:
    C = (q2^1.15 - q1^1.15) / (1.15 * alpha);

subject to var_D:
    D = (q2^(-0.7) - q1^(-0.7)) / (0.7 * alpha^2);

subject to var_E:
    E = (q2^0.15 - q1^0.15) / (0.15 * alpha);

# -------------------------------
# Solver Option and Output
# -------------------------------
option solver scip;

solve;

display q1, obj;
