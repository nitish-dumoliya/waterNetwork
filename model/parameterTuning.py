import sys
import numpy as np
from amplpy import AMPL

# Initialize AMPL
ampl = AMPL()

# Read model and data
ampl.read(sys.argv[1])
ampl.read_data(sys.argv[2])
ampl.option["solver"] = "ipopt"

# Define parameter search space
bound_push_values = np.linspace(1e-6, 0.004 , 6)  # Values: 0, 0.1, ..., 0.5
bound_frac_values = np.linspace(1e-6, 0.02 , 6)

best_objective = float("inf")
best_params = None

# Store results
results = []

for bound_push in bound_push_values:
    for bound_frac in bound_frac_values:
        print(f"Testing bound_push={bound_push}, bound_frac={bound_frac}")

        # Set IPOPT options dynamically
        ipopt_options = f"outlev=0 expect_infeasible_problem=yes bound_push={bound_push} bound_frac={bound_frac} nlp_scaling_method=gradient-based"
        ampl.set_option("ipopt_options", ipopt_options)
        ampl.option["presolve_eps"] = "8.53e-15"

        # Solve model
        ampl.solve()

        # Extract objective value
        total_cost = ampl.getObjective("total_cost").value()
        print("total cost:", total_cost, "best cost:", best_objective)
        
        # Extract constraint violations
        #violation = ampl.getValue("sum {c in con2} abs(con2[c].dual)")

        # Store results
        results.append((bound_push, bound_frac, total_cost))

        # Update best parameters
        if total_cost < best_objective :  # Ensure feasibility
            best_objective = total_cost
            best_params = (bound_push, bound_frac)

# Print best found parameters
print("\nBest Parameters:")
print(f"Bound Push: {best_params[0]}")
print(f"Bound Frac: {best_params[1]}")
print(f"Best Objective Value: {best_objective}")

# Extract solution values
l_init = ampl.getVariable('l').getValues().to_dict()
q_init = ampl.getVariable('q').getValues().to_dict()
h_init = ampl.getVariable('h').getValues().to_dict()

# Close AMPL session
ampl.close()
