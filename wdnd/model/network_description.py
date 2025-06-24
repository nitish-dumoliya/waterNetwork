import sys
import csv
from amplpy import AMPL

ampl = AMPL()

data_list = [
    "d1_bessa",
    "d2_shamir",
    "d3_hanoi",
    "d4_double_hanoi",
    "d5_triple_hanoi",
    "d6_newyork",
    "d7_blacksburg",
    "d8_fossolo_iron",
    "d9_fossolo_poly_0",
    "d10_fossolo_poly_1",
    "d11_kadu",
    "d12_pescara",
    "d13_modena",
    "d14_balerma"
]

# Prepare CSV file
csv_filename = "network_summary.csv"
csv_file = open(csv_filename, mode='w', newline='')
csv_writer = csv.writer(csv_file)

# Write header row
csv_writer.writerow([
    "Network", "Min Length", "Max Length", "Total Length",
    "Min Demand", "Max Demand", "Max Flow", "Min Dia", "Max Dia", "Epsilon"
])

# Loop over network instances
for i in range(len(data_list)):
    print("\nWater Network:", data_list[i], "\n")
    data = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{data_list[i]}.dat"
    ampl.reset()

    # Load appropriate model
    if i == 5:
        ampl.read("newyork_model.mod")
    elif i == 6:
        ampl.read("blacksburg_model.mod")
    else:
        ampl.read("wdnmodel.mod")

    # Load data
    ampl.read_data(data)

    # Solver settings
    ampl.option["solver"] = "ipopt"
    ampl.set_option("ipopt_options",
        "outlev=0 expect_infeasible_problem=yes tol=1e-9 bound_relax_factor=0 "
        "bound_push=0.01 bound_frac=0.01 nlp_scaling_method=none"
    )
    ampl.option["presolve_eps"] = "8.53e-15"

    # Length statistics
    L_values = [item[2] for item in ampl.getParameter('L').getValues().to_list()]
    min_length = min(L_values)
    max_length = max(L_values)
    total_length = sum(L_values)

    # Demand, flow, diameter limits
    min_demand = ampl.getParameter('D_min').getValues().to_list()[0]
    max_demand = ampl.getParameter('D_max').getValues().to_list()[0]
    max_flow = ampl.getParameter('Q_max').getValues().to_list()[0]
    d_max = ampl.getParameter('d_max').getValues().to_list()[0]
    d_min = ampl.getParameter('d_min').getValues().to_list()[0]
    #L = ampl.getParameter('L').getValues().to_dict()
    
    L = [item[2] for item in ampl.getParameter('L').getValues().to_list()]
    R = [item[1] for item in ampl.getParameter('R').getValues().to_list()]
    R_min = min(R)
    L_max = max(L)
    print(R_max)
    MaxK = 10.67*L_max/((R_min**1.852) * (d_min**4.87))
 
    epsilon = ((10**(-6))/(0.07508*MaxK))**(1/0.926)
        
    # Write one row to CSV (unrounded values)
    csv_writer.writerow([
        data_list[i],
        min_length,
        max_length,
        total_length,
        min_demand,
        max_demand,
        max_flow,
        d_min,
        d_max,
        epsilon
    ])

# Close CSV file
csv_file.close()
print(f"\n✅ Summary written to: {csv_filename}")

