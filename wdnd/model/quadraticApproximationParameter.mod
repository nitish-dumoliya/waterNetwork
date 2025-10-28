# -------------------------------
# Parameters
# -------------------------------
param q2 := 0.42043;
param eps_tol:=1e-1;
param L:=1000;
param R := 120;
param d := 0.1;
param alpha := (10.65*L)/(R^1.852 * d^4.871);

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
minimize obj: q1;

# -------------------------------
# Constraints
# -------------------------------
subject to root:
    ((a * q_eps^2 + b * q_eps) - alpha * (q_eps^1.85) * (1-eps_tol)) = 0;

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
    D = (q2^(-0.7) - q1^(-0.7)) / (-0.7 * alpha^2);

subject to var_E:
    E = (q2^0.15 - q1^0.15) / (0.15 * alpha);

# -------------------------------
# Solver Option and Output
# -------------------------------
option solver baron;

option baron_options "outlev = 1";

solve;

display q1, obj;
