option subsystems;

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
	  MaxK(N,M)
	  epsilon(N,M)
;

Scalar omega / 10.67 / ;
*Scalar epsilon / 1e-6 / ;


Variables q(N,M)
          total_cost
;

Positive Variable length(N,M,pipes);
Positive Variable h(N);

Equations obj
          Balance(nodes)
          HydraulicHead_Flow(N,M)
	  length_limit(N,M)
	  elevation_constraint(Source)
;

obj.. total_cost =E= sum(arcs, sum(pipes, length(arcs,pipes)*C(pipes)));

Balance(nodes)$(not Source(nodes))..
  D(nodes) =E=  - SUM(M$arcs(nodes,M), q(nodes,M)) + SUM(N$arcs(N,nodes), q(N,nodes));

*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E= signPower(q(N,M),1.852) * sum(pipes, (omega * length(N,M,pipes)) / ((R(pipes)**1.852) * ((dia(pipes))**4.87) ));

*------------------------------------------
* Headloss constraint for each arc
*------------------------------------------
HydraulicHead_Flow(N,M)$(arcs(N,M))..
    h(N) - h(M) =E= (q(N,M)**3 * (q(N,M)**2 + epsilon(N,M)**2)**0.426 / (q(N,M)**2 + 0.426 * epsilon(N,M)**2)) *  sum(pipes, (omega * length(N,M,pipes)) / ((R(pipes)**1.852) * ((dia(pipes))**4.87) ));
*h(N) - h(M) =E= (q(N,M)*(q(N,M)**2 + epsilon(N,M)**2)**0.426) *  sum(pipes, (omega * length(N,M,pipes)) / ((R(pipes)**1.852) * ((dia(pipes))**4.87) ));


*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M)
*    =E=
*    q(N,M) * (q(N,M)**2 + epsilon(N,M)**2)** 0.426
*    * sum(pipes,
*        (omega * length(N,M,pipes))
*        / (R(pipes)**1.852 * dia(pipes)**4.87)
*      );

*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E= sign(q(N,M)) * ((abs(q(N,M)) )**1.852)*(0.001**1.852) * sum(pipes, (omega * length(N,M,pipes)) / ((R(pipes)**1.852) * ((dia(pipes)/1000)**4.87) ));

*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E= ((q(N,M) * abs(q(N,M)) * ((abs(q(N,M)) + epsilon[N,M])**0.852)) / (abs(q(N,M)) + 0.852*epsilon[N,M])) * sum(pipes, (omega * length(N,M,pipes)) / ((R(pipes)**1.852) * ((dia(pipes))**4.87)) );

length_limit(N,M)$(arcs(N,M)) ..
    sum(pipes, length(N,M,pipes)) =E= L(N,M);

elevation_constraint(Source)..
    h(Source) =E= elevation(Source);

$include d1.gmsdat

Parameter Q_max;
Q_max = sum(nodes$(not Source(nodes)), D(nodes));

Scalar d_min;
d_min = smin(pipes, dia(pipes));
Scalar d_max;
d_max = smax(pipes, dia(pipes));

*---------------------------
* Minimum pipe roughness
*---------------------------
Scalar R_min;
R_min = smin(pipes, R(pipes));

*---------------------------
* Minimum pipe diameter
*---------------------------

*---------------------------
* Maximum K for each arc
*---------------------------
Parameter MaxK(N,M);
MaxK(N,M)$(arcs(N,M)) = omega * L(N,M) / (R_min**1.852 * d_min**4.87);

*---------------------------
* Epsilon for smoothing in headloss
*---------------------------
Parameter epsilon(N,M);
epsilon(N,M)$(arcs(N,M)) = 0.0535 * (1e-3 / MaxK(N,M))**0.54;
*epsilon(N,M)$(arcs(N,M)) = 1e-1;

h.lo(nodes) = P(nodes) + elevation(nodes);
q.up(arcs(N,M)) =  Q_max;
q.lo(arcs(N,M)) =  -Q_max;
length.up(N,M,pipes) =  L(N,M);
length.lo(N,M,pipes) =  0;

model gms_water_nlp / all /;


Option solver = knitro;

*gms_water_nlp.optfile = 1;

Option solprint = off ;


solve gms_water_nlp using NLP min total_cost;

display total_cost.l;
display q.l;
display h.l;
display length.l;
