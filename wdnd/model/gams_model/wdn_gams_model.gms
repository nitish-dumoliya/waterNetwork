* ============================================================
* Water Distribution Network Design Model
* Hazen-Williams with Smooth Approximations
* ============================================================

* ---- Set Declarations ----
set nodes
set Source(nodes)
set pipes 
;

alias(nodes,N,M);

set arcs(N,M)
;

Parameter elevation(nodes)
          D(nodes)
          P(nodes)
          L(N,M)
          vmax(N,M)
          R(pipes)
          dia(pipes)
          C(pipes)
	  epsi(N,M)
;

Variables h(nodes)
          q(N,M)
          total_cost
;

Positive Variable length(N,M,pipes);

Equations obj
          Balance(nodes)
          HydraulicHead_Flow(N,M)
	  length_limit(N,M)
	  elevation_constraint(Source)
;

Scalar omega / 10.67 / ;
Scalar epsilon / 1e-6 / ;
* ============================================================
* Equation Definitions
* ============================================================
obj..
    total_cost =E= sum(arcs, sum(pipes, length(arcs, pipes) * C(pipes)));

Balance(nodes)$(not Source(nodes))..
    D(nodes) =E=
        - sum(M$arcs(nodes,M), q(nodes,M))
        + sum(N$arcs(N,nodes), q(N,nodes));

* ---- Head-loss: uncomment one formulation ----

* (1) Original Hazen-Williams via signPower (DNLP):
*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E=
*        signPower(q(N,M), 1.852)
*        * sum(pipes,
*            omega * length(N,M,pipes)
*            / ( R(pipes)**1.852 * dia(pipes)**4.87 ));

* (2) Original HW explicit form (DNLP):
*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E=
*        q(N,M) * abs(q(N,M))**0.852
*        * sum(pipes,
*            omega * length(N,M,pipes)
*            / ( R(pipes)**1.852 * dia(pipes)**4.87 ));

* (3) Smooth Approximation 1: f1 = K*q*(q^2+eps^2)^0.426 (NLP):
*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E=
*        q(N,M) * ( q(N,M)**2 + epsi(N,M)**2 )**0.426
*        * sum(pipes,
*            omega * length(N,M,pipes)
*            / ( R(pipes)**1.852 * dia(pipes)**4.87 ));

* (4) Smooth Approximation 2: f2 = K*q^3*(q^2+eps^2)^0.426/(q^2+0.426*eps^2) (NLP):
HydraulicHead_Flow(N,M)$(arcs(N,M))..
    h(N) - h(M) =E=
        q(N,M)**3
        * ( q(N,M)**2 + epsi(N,M)**2 )**0.426
        / ( q(N,M)**2 + 0.426 * epsi(N,M)**2 )
        * sum(pipes,
            omega * length(N,M,pipes)
            / ( R(pipes)**1.852 * dia(pipes)**4.87 ));

*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E=
*        ( q(N,M)**3 * ( q(N,M)**2 + eps(N,M)**2 )**0.426
*          / ( q(N,M)**2 + 0.426 * eps(N,M)**2 ) )
*        * sum(pipes,
*            omega * length(N,M,pipes)
*            / ( R(pipes)**1.852 * dia(pipes)**4.87 ));

length_limit(N,M)$(arcs(N,M))..
    sum(pipes, length(N,M,pipes)) =E= L(N,M);

elevation_constraint(Source)..
    h(Source) =E= elevation(Source);

* ============================================================
* Load Data
* ============================================================
*$include d2_small.gmsdat
$include d1_bessa.gmsdat

* ============================================================
* Derived Parameters (computed after data is loaded)
* ============================================================
Parameter R_min;
Parameter R_max;
Parameter d_min;
R_min = smin(pipes, R(pipes));
R_max = smax(pipes, R(pipes));
d_min = smin(pipes, dia(pipes));

Parameter MaxK(N,M);
MaxK(N,M)$(arcs(N,M)) = omega * L(N,M) / ( R_min**1.852 * d_min**4.87 );

* ---- Epsilon selection: uncomment one ----
*Parameter epsi(N,M);

* Relative error for smooth approximation 1 (active):
*epsi(N,M)$(arcs(N,M)) = 0.0953 * ( 1e-2 / MaxK(N,M) )**0.54 ;

* Relative error for smooth approximation 2:
*epsi(N,M)$(arcs(N,M)) = 0.0535 * ( 1e-2 / MaxK(N,M) )**0.54 ;

* Absolute error for smooth approximation 1:
*epsi(N,M)$(arcs(N,M)) = ( 1e-5 / (0.36061 * MaxK(N,M)) )**0.54;

* Absolute error for smooth approximation 1:
epsi(N,M)$(arcs(N,M)) = ( 1e-2 / (0.07508 * MaxK(N,M)) )**0.54;

Parameter Q_max;
Q_max = sum(nodes$(not Source(nodes)), D(nodes));

* ============================================================
* Variable Bounds
* ============================================================
h.lo(nodes)          = P(nodes) + elevation(nodes);
q.up(arcs(N,M))      =  Q_max;
q.lo(arcs(N,M))      = -Q_max;
*q.l(arcs(N,M))       =  0.01;
length.lo(N,M,pipes) =  0;
length.up(N,M,pipes) =  L(N,M);

* Better initial point
q.l(arcs(N,M))  = Q_max / card(arcs);
h.l(nodes)      = elevation(nodes) + P(nodes) + 10;
length.l(N,M,pipes) = L(N,M) / card(pipes);

* ============================================================
* Solve
* ============================================================
model gms_water_nlp / all /;
Option solver   = baron;
Option solprint = off;
option subsystems;

* Use DNLP for formulations 1, 2, 5 (non-smooth)
* Use NLP  for formulations 3, 4  (smooth approximations)
gms_water_nlp.optfile = 1;

solve gms_water_nlp using NLP min total_cost;

* ============================================================
* Display Results
* ============================================================
display total_cost.l;
display q.l;
display h.l;
display length.l;
