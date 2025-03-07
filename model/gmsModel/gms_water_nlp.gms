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

obj.. total_cost =E= sum(arcs, sum(pipes, length(arcs,pipes)*C(pipes)));

Balance(nodes)$(not Source(nodes))..
  D(nodes) =E=  - SUM(M$arcs(nodes,M), q(nodes,M)) + SUM(N$arcs(N,nodes), q(N,nodes));

HydraulicHead_Flow(N,M)$(arcs(N,M))..
    h(N) - h(M) =E= signPower(q(N,M),1.852)*(0.001**1.852) * sum(pipes, (omega * length(N,M,pipes)) / ((R(pipes)**1.852) * ((dia(pipes)/1000)**4.87) ));

*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E= sign(q(N,M)) * ((abs(q(N,M)) )**1.852)*(0.001**1.852) * sum(pipes, (omega * length(N,M,pipes)) / ((R(pipes)**1.852) * ((dia(pipes)/1000)**4.87) ));

*HydraulicHead_Flow(N,M)$(arcs(N,M))..
*    h(N) - h(M) =E= ((q(N,M) * abs(q(N,M)) * ((abs(q(N,M)) + 1000*epsilon)**0.852)) / (abs(q(N,M)) + 852*epsilon))*(0.001**1.852) * sum(pipes, (omega * length(N,M,pipes)) / ((R(pipes)**1.852) * ((dia(pipes)/1000)**4.87)) );

length_limit(N,M)$(arcs(N,M)) ..
    sum(pipes, length(N,M,pipes)) =E= L(N,M);

elevation_constraint(Source)..
    h(Source) =E= elevation(Source);

$include d1.gmsdat

Parameter Q_max;
Q_max = sum(nodes$(not Source(nodes)), D(nodes));

*dmin = smin(pipe, dia(pipe));
*dmax = smax(pipe, dia(pipe));

h.lo(nodes) = P(nodes) + elevation(nodes);
q.up(arcs(N,M)) =  Q_max;
q.lo(arcs(N,M)) =  -Q_max;
length.up(N,M,pipes) =  L(N,M);
length.lo(N,M,pipes) =  0;

model gms_water_nlp / all /;

Option solver = IPOPT;

Option solprint = off ;

solve gms_water_nlp using DNLP min total_cost;

display total_cost.l;
display q.l;
display h.l;
display length.l;
