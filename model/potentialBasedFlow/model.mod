
param samples > 0;
param demand{1..samples} >= 0;
param cost > 0;
param recover < cost;
param retail >= 0;
param minrev;
param maxrev;
param alpha >= 0, <= 1;
param beta >= 0, <= 1;

var nu >= minrev, <= maxrev;
var excess{1..samples} >= 0, <= maxrev - minrev;
var order >= 0;
var sales{i in 1..samples} >= 0, <= demand[i];
var discount{1..samples} >= 0;
var profit{1..samples} >= minrev, <= maxrev;

var cvar = nu + 1 / ((1 - alpha) * samples) * sum{i in 1..samples} excess[i];
var average_profit = (1/samples) * sum{i in 1..samples} profit[i];

maximize prof:
    - beta * cvar + (1-beta) * average_profit;

s.t. c1 {i in 1..samples}: profit[i] == -cost * order + retail * sales[i] + recover * discount[i];
s.t. c2 {i in 1..samples}: sales[i] + discount[i] == order;
s.t. c3 {i in 1..samples}: -profit[i] - nu <= excess[i];
