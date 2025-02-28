$ONTROL IPOPT

SETS
    nodes /1*7/
    pipes /1*14/
    arcs(nodes,nodes) /1.2, 2.4, 2.3, 3.5, 4.5, 4.6, 7.5, 6.7/
    Source(nodes) /1/;

PARAMETERS
    L(arcs) /1.2 1000, 2.4 1000, 2.3 1000, 3.5 1000, 4.5 1000, 4.6 1000, 7.5 1000, 6.7 1000/
    vmax(arcs) /1.2 2, 2.4 2, 2.3 2, 3.5 2, 4.5 2, 4.6 2, 7.5 2, 6.7 2/
    d(pipes) /1 25.4, 2 50.8, 3 76.2, 4 101.6, 5 152.4, 6 203.2, 7 254, 8 304.8, 9 355.6, 10 406.4, 11 457.2, 12 508, 13 558.8, 14 609.6/
    E(nodes) /1 210, 2 150, 3 160, 4 155, 5 150, 6 165, 7 160/
    P(nodes) /1 0, 2 30, 3 30, 4 30, 5 30, 6 30, 7 30/
    D(nodes) /1 -311.1087, 2 27.7777, 3 27.777, 4 33.333, 5 75, 6 91.666, 7 55.555/
    C(pipes) /1 2, 2 5, 3 8, 4 11, 5 16, 6 23, 7 32, 8 50, 9 60, 10 90, 11 130, 12 170, 13 300, 14 550/
    R(pipes) /1 130, 2 130, 3 130, 4 130, 5 130, 6 130, 7 130, 8 130, 9 130, 10 130, 11 130, 12 130, 13 130, 14 130/;

VARIABLES
    l(arcs,pipes)
    q(arcs)
    h(nodes)
    t(nodes);

POSITIVE VARIABLES l, q, t;
EQUATIONS
    con1(nodes),
    con2(arcs),
    con3(arcs),
    con4(arcs,pipes),
    con5(Source),
    con6(nodes),
    con7(arcs),
    con8(arcs),
    objective;

objective.. sum((i,j), sum(k, l(i,j,k) * C(k))) =E= z;

con1(j).. sum((i,j), q(i,j)) - sum((j,i), q(j,i)) =E= D(j);

con2(i,j).. h(i) - h(j) =E= ((0.001^1.852) * q(i,j) * (abs(q(i,j)) + 1000 * eps)^0.852 - (0.002368316 * eps * q(i,j) / (abs(q(i,j)) + 1000 * eps)^0.148) + ((0.175255362 * (eps)^2) * q(i,j) / ((abs(q(i,j)) + 1000 * eps)^1.148))) * sum(k, (10.67 * l(i,j,k) / ((R(k)^1.852) * (d(k)/1000)^4.87)));

con3(i,j).. sum(k, l(i,j,k)) =E= L(i,j);

con4(i,j,k).. l(i,j,k) =L= L(i,j);

con5(i).. h(i) =E= E(i);

con6(i).. h(i) =G= E(i) + P(i);

con7(i,j).. q(i,j) =G= -Q_max;

con8(i,j).. q(i,j) =L= Q_max;

MODEL water_network /all/;
SOLVE water_network USING NLP MINIMIZING z;

DISPLAY z.l;

