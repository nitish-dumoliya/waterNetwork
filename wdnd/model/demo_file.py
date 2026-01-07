from matplotlib.cbook import violin_stats
import networkx as nx
from amplpy import AMPL
import matplotlib.pyplot as plt
from networkx.algorithms import cycle_basis
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
import optuna
import warnings
warnings.filterwarnings("ignore")
from pyswarm import pso
import multiprocessing
from network_layout import node_position
import math
class WaterNetworkOptimizer:
    def __init__(self, model_file, data_file, data_number, data_list):
        self.ampl = AMPL()
        self.model_file = model_file
        self.data_list = data_list
        self.data_file = data_file
        self.data_number = data_number
        self.total_cost = None
        self.network_graph = None
        self.solve_result = None
        self.solver_time = 0
        self.best_acyclic_flow = None
        self.number_of_nlp = 0
        self.tol = 0

    def load_model(self):
        """Load the model and data."""
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file) 
        self.nodes = self.ampl.getSet('nodes')
        self.source = self.ampl.getSet('Source')
        self.arcs = self.ampl.getSet('arcs')
        self.pipes = self.ampl.getSet('pipes')
        if self.data_number == 6:
            self.fixarcs = self.ampl.getSet('fixarcs')
        self.L = self.ampl.getParameter('L').to_dict()
        self.D = self.ampl.getParameter('D').to_dict()
        self.C = self.ampl.getParameter('C').to_dict()
        self.P = self.ampl.getParameter('P').to_dict()
        self.R = self.ampl.getParameter('R').to_dict()
        self.E = self.ampl.getParameter('E').to_dict()
        self.d = self.ampl.getParameter('d').to_dict()
        #self.eps = self.ampl.getParameter('eps').to_dict()
        self.delta = 0.1
        self.p = 1.852
        self.omega = 10.67


    def solve(self):
        self.ampl.option["solver"] = "ipopt"
        # self.ampl.set_option("ipopt_options", f"outlev = 0 tol = 1e-9 bound_relax_factor=0  bound_push = {self.bound_push} bound_frac = {self.bound_frac} halt_on_ampl_error = yes warm_start_init_point = no expect_infeasible_problem = no")   #max_iter = 1000
        # self.ampl.set_option("ipopt_options", f"outlev = 0 tol = 1e-9 bound_relax_factor=0  bound_push = 0.01 bound_frac = 0.01 halt_on_ampl_error = yes halt_on_ampl_error = yes warm_start_init_point = no expect_infeasible_problem = no")   #max_iter = 1000
        # self.ampl.option["ipopt_options"] = "outlev = 0 expect_infeasible_problem = yes bound_relax_factor=0 bound_push = 0.01 bound_frac = 0.01 warm_start_init_point = yes halt_on_ampl_error = yes "
        #self.ampl.set_option("ipopt_options", f"outlev = 0  bound_relax_factor=0 warm_start_init_point = no halt_on_ampl_error = yes")   #max_iter = 1000
        #self.ampl.set_option("ipopt_options", f"outlev = 0 warm_start_init_point = no ")   #max_iter = 1000
        self.ampl.option["presolve_eps"] = "6.82e-14"
        self.ampl.option['presolve'] = 0


        self.ampl.solve()
        self.solve_result = self.ampl.solve_result
        self.total_cost = self.ampl.get_objective("total_cost").value()
        #self.ampl.eval("display q;")
        # print("Objective:", self.total_cost)
        # print("solve_result: ",self.solve_result)
        solve_time = self.ampl.get_value('_solve_elapsed_time')
        self.solver_time += solve_time
        self.number_of_nlp += 1

    def run(self):
        """Main function to run the Heuristic Approach."""
        self.start_time = time.time()
        self.bound_push , self.bound_frac = (0.1, 0.1)
        self.mu_init = 0.1
        print("NLP Model: Smooth Approximatate WDN Design Model 2, Epsilon Selection using Relative Error")
        print("NLP Solver: Ipopt")
        print("************************Initial Solution of Approximate WDN Design Model**********************")
        self.load_model()
        fix_arc_set = self.fix_leaf_arc_flow()
        print("fix_arc_set:",fix_arc_set)
        self.super_source_out_arc = self.fix_arc_set()
        print("super_source_out_arc:", self.super_source_out_arc, "\n")

        self.solve()
        print("Objective: ",self.total_cost)
        print("Solve_result: ",self.solve_result)
        print("Solve_time:", self.ampl.get_value('_solve_elapsed_time'),"\n")
        if self.solve_result != 'solved':
            print("IPOPT did not solve the initial problem optimally. Exiting Heuristic.")
            sys.exit()
        self.current_cost = self.total_cost
        self.l = self.ampl.getVariable('l').getValues().to_dict()
        self.q = self.ampl.getVariable('q').getValues().to_dict()
        self.h = self.ampl.getVariable('h').getValues().to_dict()
        if self.data_number==5:
            self.q1 = self.ampl.getVariable('q1').getValues().to_dict()
            self.q2 = self.ampl.getVariable('q2').getValues().to_dict()
        print("*****************************Improve the Initial Solution*************************************\n")
        self.super_source_out_arc = self.fix_arc_set()
        self.network_graph = self.generate_random_acyclic_from_solution(self.q)
        # print("Fix the flow direction in optimization model and solve the updated model")
        self.indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        self.fix_arc_set = list(set(self.super_source_out_arc) | fix_arc_set)        
        self.best_acyclic_flow = self.network_graph.copy() 
        print("\n--------------------------Global Search Method-------------------------------")
        self.iteration = 1
        # self.global_search()

        print("\n************************************Final Best Results*****************************************")
        print("Water Network:", self.data_list[self.data_number])
        # self.eps = self.ampl.get_variable('eps').get_values().to_dict()
        self.eps = self.ampl.getParameter('eps').getValues().to_dict()
        # cost = 0
        # for (i,j) in self.arcs:
        #     for k in self.pipes:
        #         cost = cost + self.C[k]*self.l[i,j,k]
        # print("cost:", cost)
        print(f"Final best objective: {self.current_cost}")
        #self.ampl.eval("display q;")
        print("Number of nlp problem solved:", self.number_of_nlp)
        print("Total number of iteration:", self.iteration + self.headloss_increase_iteration + self.dia_red_iteration)
        self.constraint_violations(self.q, self.h, self.l, self.eps, "ipopt")
        elapsed_time = time.time() - self.start_time
        solver_time = self.solver_time
        print(f"Solver_time: {solver_time:.2f} seconds")
        # print(f"Heuristic elapsed time: {elapsed_time:.2f} seconds = {elapsed_time/60:.2f} minutes.\n")
        print(f"Heuristic elapsed time: {elapsed_time:.2f} seconds")
        print("***********************************************************************************************\n")

if __name__ == "__main__":
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
    # Select the data number here (0 to 18)
    #file_name = sys.argv[1]
    data_number = int(sys.argv[1]) - 1
    #input_data_file = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{file_name}"
    input_data_file = f"/home/nitishdumoliya/waterNetwork/wdnd/data/{data_list[(data_number)]}.dat"
    # print(data_number)
    print("Water Network:", data_list[(data_number)])
    if data_number==5:
        optimizer = WaterNetworkOptimizer("newyork_model.mod", input_data_file, data_number, data_list)
    elif data_number == 6:
        optimizer = WaterNetworkOptimizer("blacksburg_model.mod", input_data_file, data_number, data_list)
    else:
        optimizer = WaterNetworkOptimizer("wdnmodel.mod", input_data_file, data_number, data_list)
    # optimizer = WaterNetworkOptimizer(sys.argv[1], sys.argv[2], sys.argv[3])
    optimizer.run()
