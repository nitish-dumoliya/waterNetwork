import optuna
import numpy as np
import networkx as nx
from amplpy import AMPL
import matplotlib.pyplot as plt
import numpy as np
import time
import copy
import sys
import os
import contextlib
from sklearn.decomposition import PCA
import random

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from tabulate import tabulate
import warnings
warnings.filterwarnings("ignore")


class AMPLTuner:
    def __init__(self, input_data_file):
        self.ampl = AMPL()
        self.bound_push = None
        self.bound_frac = None
        self.data = input_data_file

    def load_model(self):
        """ Load or reset the AMPL model """
        self.ampl.reset()
        self.ampl.read("../water-nlp.mod")  # Load your model file
        self.ampl.readData(f"{self.data}")  # Load data if needed

    def solve(self):
        """ Solve the AMPL optimization model """
        self.ampl.solve()

    def objective(self, trial):
        """ Objective function for Bayesian Optimization """
        
        # Sample hyperparameters from Optuna
        bound_push = trial.suggest_float("bound_push", 0.0001, 0.01)
        bound_frac = trial.suggest_float("bound_frac", 0.0001, 0.01)

        print(f"Testing bound_push={bound_push}, bound_frac={bound_frac}")
        self.bound_push = bound_push
        self.bound_frac = bound_frac

        # Load the model
        self.load_model()
        self.ampl.option["solver"] = "ipopt"
        # Set IPOPT options dynamically
        ipopt_options = f"outlev=0 expect_infeasible_problem=yes bound_push={bound_push} bound_frac={bound_frac} nlp_scaling_method=gradient-based"
        self.ampl.set_option("ipopt_options", ipopt_options)

        # Solve the model
        self.solve()

        # Extract objective value
        total_cost = self.ampl.getObjective("total_cost").value()
        print(f"Total cost: {total_cost}")

        return total_cost  # Minimize total cost




        
if __name__ == "__main__":
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
    data_number = int(sys.argv[1]) -1
    input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data_list[(data_number)]}.dat"
    print("Water Network:", data_list[(data_number)],"\n")

 
    # Initialize AMPL tuner
    ampl_tuner = AMPLTuner(input_data_file)
    
    # Create an Optuna study
    study = optuna.create_study(direction="minimize")
    study.optimize(ampl_tuner.objective, n_trials=5)
    
    # Print best parameters
    print("\nBest Found Parameters:")
    print(f"Bound Push: {study.best_params['bound_push']}")
    print(f"Bound Frac: {study.best_params['bound_frac']}")
    print(f"Best Objective Value: {study.best_value}")
    
