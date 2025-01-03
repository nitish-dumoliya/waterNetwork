import amplpy
import random

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
        "d14_NewYork",
        "d15_foss_poly_0",
        "d16_foss_iron",
        "d17_foss_poly_1",
        "d18_pescara",
        "d19_modena"
    ]

# Select the data number here (0 to 18)
data_number = 4
input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data_list[data_number]}.dat"
print("Water Network:", data_list[data_number],"\n")

# Initialize AMPL environment
ampl = amplpy.AMPL()

# Load the model into AMPL
ampl.read("../m1Basic.mod")

# Load data into AMPL (assumed to be available in a .dat file)
ampl.read_data(input_data_file)

max_l = max(ampl.getParameter('L').to_dict().values())
max_q = ampl.getParameter('D').getValues().toDict()
# print(max_q[1])
# for i, value in max_q.items():
#     print(value)

source = ampl.getSet('Source').to_list()
E = ampl.getParameter('E').getValues().toDict()
P = ampl.getParameter('P').getValues().toDict()

# Define the number of starts for multistart heuristic
num_starts = 50
best_cost = float('inf')
best_solution = {}

# Set a random seed for reproducibility
random.seed(num_starts)

# Loop for multistart heuristic
for start in range(num_starts):

    for (i,j) in ampl.get_set("arcs").to_list():
        for k in ampl.get_set("pipes").to_list():
            value = random.uniform(0, max_l)  
            ampl.eval(f' let l[{i},{j},{k}] := {value};')

    for (i,j) in ampl.get_set("arcs").to_list():
        value = random.uniform(max_q[1], -max_q[1])  
        ampl.eval(f'let q[{i},{j}] := {value};')

    for i in ampl.get_set("nodes").to_list():
        value = random.uniform(E[i]+P[i], E[source[0]])  
        ampl.eval(f'let h[{i}] := {value};')

    # Solve the optimization problem using IPOPT
    ampl.set_option("solver", "ipopt")
    ampl.set_option("ipopt_options", "outlev = 0")
    # ampl.option[""]
    ampl.solve()

    # Check if the solution is feasible
    if ampl.solve_result == 'solved':
        current_cost = ampl.get_objective("total_cost").value()
        print(current_cost)

        # Update the best solution if the current cost is lower
        if current_cost < best_cost:
            best_cost = current_cost
    # ampl.close()

# # Output the best solution found
print(f"Best Cost: {best_cost}")
print("Best Solution (lengths):", best_solution)

