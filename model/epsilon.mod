var epsilon>=0;
var x;

maximize relative_error:   abs(x^1.852-(x^3*(x^2 + epsilon)^0.426)/(x^2 + 0.426*epsilon)) ;
#maximize constant: 0 ;

#subject to relative_error_con: 
#    abs(1-(x^1.148*(x^2 + epsilon)^0.426)/(x^2 + 0.426*epsilon)) = 1e-6
#    (1.852 * x^6 + 2.64095 * epsilon * x^4 + 1.278 * epsilon^2 * x^2 - 1.852 * (x^2)^0.426 * (x^2 +epsilon)^0.574 * (x^2 + 0.426 * epsilon)^2) = 0
#;

subject to epsilon_bound:
    epsilon = 1e-6 
;

subject to x_bound:
    1e-6 <= x <= 1
;

option solver baron;
option ipopt_options "outlev = 1";
option mmultistart_options "--presolve 1,--log_level 6,--nlp_engine IPOPT,--eval_within_bnds 1";
#option solver baron;
option knitro_options "outlev = 0 ms_enable = 1 ms_maxsolves = 20 mip_multistart=1";
# option knitro_options "outlev =1 ";
option baron_options  "maxtime = 3600  outlev = 1";
# option conopt_options  "outlev = 1 ";

# option presolve 1;
option presolve_eps 1.71e-14;

solve ;
display epsilon;
display x;

display abs(1-(x^1.148*(x^2 + epsilon)^0.426)/(x^2 + 0.426*epsilon));  

