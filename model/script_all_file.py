import sys
from tabulate import tabulate
from amplpy import AMPL
ampl = AMPL()



data_list = [
    "d1_Sample_input_cycle_twoloop",
    "d2_Sample_input_cycle_hanoi",
    "d3_Sample_input_double_hanoi",
    "d4_Sample_input_triple_hanoi",
    "d5_Taichung_input",
    "d6_HG_SP_1_4",
    "d7_HG_SP_2_3",
    "d8_HG_SP_3_4",
    "d9_HG_SP_4_2",
    "d10_HG_SP_5_5",
    "d11_HG_SP_6_3",
    "d12",
    "d13",
    "d15_foss_poly_0",
    "d16_foss_iron",
    "d17_foss_poly_1",
    "d18_pescara",
    "d19_modena"
]

for i in range(len(data_list)):
    print("\nWater Network:", data_list[i], "\n")
    data = f"/home/nitishdumoliya/waterNetwork/data/{data_list[i]}.dat"
    ampl.reset()
    ampl.read("water-nlp2.mod")
    ampl.read_data(data)
    ampl.option["solver"] = "ipopt"
    ampl.set_option("ipopt_options", "outlev = 0  expect_infeasible_problem = yes tol = 1e-9 bound_relax_factor=0  bound_push = 0.01 bound_frac = 0.01 nlp_scaling_method = none")   #max_iter = 1000

    ampl.option["presolve_eps"] = "8.53e-15"

    
    min_demand = ampl.getParameter('D_min').getValues().to_list()[0]
    max_demand = ampl.getParameter('D_max').getValues().to_list()[0]
    max_flow = ampl.getParameter('Q_max').getValues().to_list()[0]
    d_max = ampl.getParameter('d_max').getValues().to_list()[0]
    d_min = ampl.getParameter('d_min').getValues().to_list()[0]
    
    print("min_demand:", min_demand/1000)
    print("max_demand:", max_demand/1000)
    print("max_flow:", max_flow/1000)

    nodes = ampl.getSet('nodes')
    source = ampl.getSet('Source')
    arcs = ampl.getSet('arcs')
    pipes = ampl.getSet('pipes')
    
    L = ampl.getParameter('L').to_dict()
    D = ampl.getParameter('D').to_dict()
    C = ampl.getParameter('C').to_dict()
    P = ampl.getParameter('P').to_dict()
    R = ampl.getParameter('R').to_dict()
    E = ampl.getParameter('E').to_dict()
    d = ampl.getParameter('d').to_dict()
    
    max_L = max(L[i,j] for (i,j) in arcs)
    R_min = min(R[k] for k in pipes)


    MaxK = 10.67 / ((R_min ** 1.852) * ((d_min) ** 4.87))

    print("K_rhs:", MaxK)
        
    epsilon = (10**(-6)/(0.07508*MaxK))**(1/0.926)
    
    print("epsilon:", epsilon)

    ampl.close 
