# -------------------------------
# Parameters
# -------------------------------
param q2 := 0.42043;
param eps_tol:=1e-1;
param L:=1000;
param R := 120;
param d := 0.1;
param alpha := (10.67*L)/(R^1.852 * d^4.871);

# -------------------------------
# Variables
# -------------------------------
var a;
var b;
var q1>=0;
var B1;
var B2;
var B3;
var B4;
var B5;
var q_eps>=0;

# -------------------------------
# Objective Function
# -------------------------------
minimize obj: q1;

# -------------------------------
# Constraints
# -------------------------------
subject to root:
    (a * q_eps^2 + b * q_eps) - alpha * (q_eps^1.852) * (1-eps_tol) = 0;

subject to var_q_eps:
    q_eps = 0.852 * b / (0.148 * a);

subject to var_a:
    a = (alpha*B4 - b * B3) / B1;

subject to var_b:
    b = alpha*((B5*B1) - (B4*B3)) / (B2*B1 - B3^2);

subject to var_A:
    B1 = (q2^1.296 - q1^1.296) / (1.296);

subject to var_B:
    B2 = (q2^(-0.704) - q1^(-0.704)) / (-0.704);

subject to var_C:
    B3 = (q2^0.296 - q1^0.296) / (0.296);

subject to var_D:
    B4 = (q2^(1.148) - q1^(1.148)) / (1.148);

subject to var_E:
    B5 = (q2^0.148 - q1^0.148) / (0.148);

# -------------------------------
# Solver Option and Output
# -------------------------------
option solver scip;

option baron_options "outlev = 1";

solve;

display q1, obj;
