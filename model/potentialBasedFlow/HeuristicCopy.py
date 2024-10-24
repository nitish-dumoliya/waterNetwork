import networkx as nx
from amplpy import AMPL
import matplotlib.pyplot as plt
import time
import copy
import sys
import os
import contextlib

class WaterNetworkOptimizer:
    def __init__(self, model_file, data_file, data_number):
        self.ampl = AMPL()
        self.model_file = model_file
        self.data_file = data_file
        self.data_number = data_number
        self.total_cost = None
        self.network_graph = None
        self.solve_result = None
        self.solver_time = 0

    def load_model(self):
        """Load the model and data."""
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)
    

    def create_digraph(self):
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = self.ampl.getSet('arcs').to_list()
        self.network_graph = nx.DiGraph()
        self.network_graph.add_nodes_from(nodes_list)
        self.network_graph.add_edges_from(edges_list)
        print(nodes_list)
        print(edges_list)
    
    def display_results(self):
        """Display relevant results from the optimization."""
        self.ampl.eval("display {(i,j) in arcs}: q[i,j];")
        self.ampl.eval("display h;")
        self.ampl.eval("display solve_result;")
        self.total_cost = self.ampl.get_objective("total_cost").value()
        print("Objective:", self.total_cost)
        # self.ampl.eval("display {(i,j) in arcs} h[i] - h[j];")
        # self.ampl.eval("display {i in nodes} h[i] - (E[i] + P[i]);")
        
    def display_best_results(self,ampl):
        ampl.eval("display {(i,j) in arcs}: q[i,j];")
        ampl.eval("display l;")
        ampl.eval("display h;")
        ampl.eval("display solve_result;")
        total_cost = ampl.get_objective("total_cost").value()
        print("Objective:", total_cost)
    
    # def plot_graph(self):
    #     print("Edges of the graph:",self.network_graph.edges())
    #     plt.figure(figsize=(10, 8))
    #     nx.draw_spectral(self.network_graph, with_labels=True)
    #     plt.show()


    def plot_graph(self, super_source_out_arc=None, best_arc=None):
        # Print the edges of the graph
        print("Edges of the graph:", self.network_graph.edges())
        # Identify nodes with indegree >= 2
        indegree_2_or_more = [node for node, indeg in self.network_graph.in_degree() if indeg >= 2]
        # Draw the graph using spectral layout
        plt.figure(figsize=(10, 8))
        pos = nx.spectral_layout(self.network_graph)
        # Draw all nodes in default color (light blue)
        nx.draw_networkx_nodes(self.network_graph, pos, node_color='lightblue', node_size=200)
        # Highlight nodes with indegree >= 2 in a different color (orange)
        if indegree_2_or_more:
            nx.draw_networkx_nodes(self.network_graph, pos, nodelist=indegree_2_or_more, node_color='orange', node_size=200)
        # Draw labels for all nodes
        nx.draw_networkx_labels(self.network_graph, pos)
        # Draw all edges in default color (black)
        nx.draw_networkx_edges(self.network_graph, pos, edge_color='black')
        # If there is an edge to highlight, draw it in a different color (red)
        if super_source_out_arc:
            nx.draw_networkx_edges(self.network_graph, pos, edgelist=super_source_out_arc, edge_color='red', width=1)
        if best_arc:
            nx.draw_networkx_edges(self.network_graph, pos, edgelist=[best_arc], edge_color='magenta', width=1)
        # Show the plot
        plt.show()
        
    def cycle_basis(self):
        root = self.ampl.getSet('Source').to_list()[0]
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = self.ampl.getSet('arcs').to_list()
        
        uwg = nx.Graph()
        # Add nodes and edges to the undirected graph
        uwg.add_nodes_from(nodes_list)
        uwg.add_edges_from(edges_list)
        # print("Edges in the undirected graph:", edges_list)
        print("cycle basis for given water network: ",nx.cycle_basis(uwg, root))
    
    def generate_random_acyclic_from_solution(self):
        print("Generate the acyclic network using ipopt solution")
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = self.ampl.getSet('arcs').to_list()
        
        self.network_graph = nx.DiGraph()
        self.network_graph.add_nodes_from(nodes_list)
        
        q = self.ampl.getVariable('q').getValues().toDict()
        for (i,j) in edges_list:
            if q[i,j] >= 0:
                self.network_graph.add_edge(i,j)
            else:
                self.network_graph.add_edge(j,i)
        
    def generate_random_acyclic_graph(self):
        uwg = nx.Graph()
        
        # Retrieve nodes and edges from AMPL
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = self.ampl.getSet('arcs').to_list()
        
        # Add nodes and edges to the undirected graph
        uwg.add_nodes_from(nodes_list)
        uwg.add_edges_from(edges_list)
        print("Edges in the undirected graph:", edges_list)
        
        # Generate a random spanning tree using Wilson's algorithm
        random_tree = nx.random_spanning_tree(uwg)
        
        # Retrieve the root from the AMPL source set
        root_l = self.ampl.getSet('Source').to_list()
        root = root_l[0]
        print("Root node:", root)

        # Ensure the root is present in the random tree
        if root not in random_tree.nodes:
            raise ValueError("The specified root must be a node in the graph.")

        # Create a directed graph from the random tree starting from the specified root
        self.network_graph = nx.DiGraph()
        visited = set()

        def dfs(node):
            visited.add(node)
            for neighbor in random_tree.neighbors(node):
                if neighbor not in visited:
                    self.network_graph.add_edge(node, neighbor)  # Direct edge from parent to child
                    dfs(neighbor)

        # Start DFS from the specified root
        dfs(root)

        # Draw the initial directed tree
        plt.figure(figsize=(15, 10))
        plt.subplot(121)
        nx.draw_spectral(self.network_graph, with_labels=True, node_color='lightgreen', font_weight='bold', arrows=True)
        plt.title("Directed Spanning Tree")

        # Add remaining edges from the original graph and check for cycles
        for u, v in uwg.edges():
            if not self.network_graph.has_edge(u, v):  # Only consider non-tree edges
                self.network_graph.add_edge(u, v)  # Try adding the edge in the (u -> v) direction
                if not nx.is_directed_acyclic_graph(self.network_graph):  # Check for cycles
                    self.network_graph.remove_edge(u, v)  # Remove the edge creating a cycle
                    self.network_graph.add_edge(v, u)  # Add it in the reverse direction (v -> u)

        # Draw the final directed graph after adding remaining edges
        plt.subplot(122)
        nx.draw_spectral(self.network_graph, with_labels=True, node_color='lightgreen', font_weight='bold', arrows=True)
        plt.title("Acyclic Directed Graph")
        plt.show()
    
    def build_acyclic_graph(self):
        """Build a directed graph from the AMPL model."""
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        edges_list = [(arc[0],arc[1]) for arc in self.ampl.getSet('arcs')]

        # Add directed edges while ensuring acyclicity
        for edge in edges_list:
            # Here you can implement a check to ensure that adding the edge won't create a cycle
            # For example, only add the edge if the target node has a higher index than the source
            source, target = edge
            if self.is_valid_edge(source, target):
                self.network_graph.add_edge(source, target)
            else:
                print("graph is not acyclic")
                print("We change the arc direction")
                self.network_graph.add_edge(target, source)
                is_dag_ = nx.is_directed_acyclic_graph(self.network_graph)  # Check for acyclicity
                if is_dag_:
                    print("graph is acyclic after changing the arc direction")
                else:
                    print("graph is not acyclic after changing the arc direction")

        print("Nodes in the graph:", nodes_list)
        print("Directed edges in the graph:", self.network_graph.edges)

        # return self.network_graph

    def update_model(self):
        # print("Fix the arcs direction using the acyclic network\n")
        edges_list = [(arc[0],arc[1]) for arc in self.ampl.getSet('arcs')]
        for edge in self.network_graph.edges:
            i, j = edge
            if edge in edges_list:
                self.ampl.eval(f"s.t. flow_direction{i}_{j}: q[{i},{j}] >=0;")
                # self.ampl.eval(f"s.t. head_bound_left{i}_{j}: E[{j}]+P[{j}] <= h[{j}];")
                # self.ampl.eval(f"s.t. head_bound_right{i}_{j}: E[{j}] + P[{j}] <= h[{i}];")
            else:
                self.ampl.eval(f"s.t. flow_direction{i}_{j}: q[{j},{i}] <=0;")
                # self.ampl.eval(f"s.t. head_bound_left{i}_{j}: E[{i}]+P[{i}] <= h[{i}];")
                # self.ampl.eval(f"s.t. head_bound_right{i}_{j}: E[{i}] + P[{i}] <= h[{j}];")

    def is_valid_edge(self, source, target):
        """Check if adding the directed edge (source -> target) maintains acyclicity."""
        self.network_graph.add_edge(source, target)  # Temporarily add the edge
        is_dag = nx.is_directed_acyclic_graph(self.network_graph)  # Check for acyclicity
        self.network_graph.remove_edge(source, target)  # Remove the edge after checking
        return is_dag  # Return True if it maintains acyclicity     
    
    def check_incoming_arcs(self):
        root = list(self.ampl.getSet('Source'))[0]
        # Iterate over all nodes in the graph
        for node in self.network_graph.nodes():
            # Skip the root node
            if node == root:
                continue
            # Check if the in-degree of the node is at least 1
            if self.network_graph.in_degree(node) < 1:
                # print(f"Node {node} does not have any incoming arcs.")
                return False
        # print("All nodes except the root have at least one incoming arc.")
        return True
    
    def is_cycle(self, graph, start_node, end_node, visited_copy, parent):
        visited_copy[start_node] = True
        # print(f"Is node {end_node} in cycle?")
        for neighbor in graph.neighbors(start_node):
            # print("visited",neighbor,visited_copy[neighbor])
            if not visited_copy[neighbor]:
                # print("neighbor of node", start_node, "is", neighbor)
                isCycle = self.is_cycle(graph, neighbor, end_node, visited_copy, start_node)
                if isCycle:
                    return True
            else:
                # print("parent:", parent)
                if parent != neighbor:
                    if end_node == neighbor:
                        # print(f"Node {end_node} is in cycle")
                        return True
        return False

    def presolve(self, graph, node, visited, parent, set_arc):
        visited_copy = visited.copy()
        # print(visited_copy)
        isCycle = self.is_cycle(graph, node, node, visited_copy, parent)
        print(f"Is node {node} in cycle?",isCycle)
        visited[node] = True
        if isCycle:
            for neighbor in graph.neighbors(node):
                if parent!=neighbor:
                    set_arc.append((node,neighbor))
                    print("Fix the arc", (node, neighbor))
            return set_arc
        else:
            for neighbor in graph.neighbors(node):
                if parent != neighbor:
                    set_arc.append((node,neighbor))
                    # print(set_arc)
                    print("Fix the arc", (node, neighbor))
                    # print("neighbor:", neighbor)
                    self.presolve(graph, neighbor, visited, node, set_arc)
        return set_arc

    def fix_arc_set(self):
        graph = nx.Graph()
        self.load_model()
        arc_set = self.ampl.getSet('arcs').to_list()
        graph.add_edges_from(arc_set)
        visited = {node: False for node in graph.nodes()}
        source = self.ampl.getSet('Source').to_list()[0]
        set_arc = []
        print("\nPresolve the model for fixing the arc direction")
        set_ = self.presolve(graph, source, visited, -1, set_arc)
        print("fixed arc direction:",set_, "\n") 
        return set_

    def fix_arc(self, data_number):
        fix_arc_list = [[(1,2), (2,3), (2,4)],
                [(1,2), (2,3), (3,4), (3,20), (3, 19)],
                [(1, 2), (2, 33), (2, 3)],
                [(1, 2), (2, 33), (2, 63), (2, 3)],
                [(20, 19), (20, 10)],
                [(1, 30), (1, 6), (1, 22)],
                [(1, 59), (1, 56)],
                [(1, 62), (62, 17), (62, 22), (62, 55)],
                [(1, 90), (1, 25), (1, 125)],
                [(1, 4), (1, 7), (1, 65)],
                [(1, 39), (1, 25), (1, 153)],
                [(1, 31), (31, 38), (31, 28)],
                [(8, 9), (9, 3), (3, 7), (3, 10), (10, 2), (2, 4), (2, 6)],
                [(1, 15), (1, 2)],
                [(37, 1), (1, 36), (1, 17), (1, 31)],
                [(37, 1), (1, 36), (1, 17), (1, 31)],
                [(37, 1), (1, 36), (1, 17), (1, 31)],
                [(15, 14), (14, 24), (24, 17), (24, 25)]
               ]
        return fix_arc_list[data_number]
    
    def iterate_arc(self,best_acyclic_flow, improved, current_cost, l, q, h, best_arc, super_source_out_arc):
        improved = False
        self.network_graph = best_acyclic_flow.copy()
        print("Acyclic network arcs direction: ",self.network_graph.edges())
        # self.plot_graph(super_source_out_arc, best_arc)
        # try:
        #     super_source_out_arc.remove(best_arc)
        # except:
        #     pass
        print("Fixed arc set:",super_source_out_arc)
        BEST_ARC = []
        BEST_ARC.append(best_arc)
        for node in self.network_graph.nodes():
            #Check if the node has an indegree >= 2
            print("Node:", node,"in_degree:", self.network_graph.in_degree(node))
            if self.network_graph.in_degree(node) >= 2:
                for u,v in list(self.network_graph.in_edges(node)):
                    if (u,v) not in super_source_out_arc and (u,v) not in BEST_ARC:
                        self.network_graph.remove_edge(u,v)
                        self.network_graph.add_edge(v,u)
                        acy_check = nx.is_directed_acyclic_graph(self.network_graph)
                        in_arc_check = self.check_incoming_arcs()
                        if acy_check and in_arc_check:
                            self.load_model()
                            self.update_model()
                            self.solve()
                            if self.solve_result == "solved":
                                if self.total_cost < current_cost:
                                    print("Arc", (u,v),"Acyclic:", acy_check and in_arc_check, "Best optimal: ", current_cost, "New optimal: ", self.total_cost, "Solve_result: ", self.solve_result, "Improved: Yes")
                                    current_cost = self.total_cost
                                    improved = True
                                    best_acyclic_flow = self.network_graph.copy()
                                    # try:
                                    #     super_source_out_arc.remove(best_arc)
                                    # except:
                                    #     pass
                                    best_arc = (v,u)
                                    l = self.ampl.getVariable('l').getValues().to_list()
                                    q = self.ampl.getVariable('q').getValues().to_list()
                                    h = self.ampl.getVariable('h').getValues().to_list()
                                    # Restore the original direction of the arc if no improvement or if the graph becomes cyclic
                                    # print(f"Restoring original direction.\n")
                                    self.network_graph.remove_edge(v, u)
                                    self.network_graph.add_edge(u, v)  # Restore the original arc direction
                                else:
                                    # print(f"No improvement with reversed edge ({u}, {v}). Restoring original direction.\n")
                                    # Restore the original direction of the arc if no improvement or if the graph becomes cyclic
                                    print("Arc", (u,v),"Acyclic:", acy_check and in_arc_check, "Best optimal: ", current_cost, "New optimal: ", self.total_cost, "solve_result: ", self.solve_result, "Improved: No")
                                    self.network_graph.remove_edge(v, u)
                                    self.network_graph.add_edge(u, v)  # Restore the original arc direction
                            else:
                                print("Arc", (u,v),"Acyclic:",acy_check and in_arc_check," Best optimal: ", current_cost, "New optimal: ", self.total_cost, "solve_result: ", self.solve_result)
                                self.network_graph.remove_edge(v, u)
                                self.network_graph.add_edge(u, v)  # Restore the original arc direction                       
                        else:
                            print("Arc", (u,v),"Acyclic:",acy_check and in_arc_check)
                            self.network_graph.remove_edge(v, u)
                            self.network_graph.add_edge(u, v)  # Restore the original arc direction                       
                print(" ")
        # print("\n")
        return best_acyclic_flow, improved, current_cost, l, q, h, best_arc
        
    def iterate_acyclic_flows(self):
        """Iterate to find improved acyclic flows by attempting arc reversals while maintaining acyclicity."""
        improved = True  # Flag to keep iterating as long as improvements are found
        
        best_acyclic_flow = self.network_graph.copy()
        
        if self.solve_result == "solved":
            current_cost = self.total_cost  # Keep track of the current cost
            l = self.ampl.getVariable('l').getValues().to_list()
            q = self.ampl.getVariable('q').getValues().to_list()
            h = self.ampl.getVariable('h').getValues().to_list()

        elif self.solve_result != "solved":
            current_cost = 10e+14
            l = None
            q = None
            h = None

        # acyclicity_time = 0
        iteration = 1
        best_arc = None
        super_source_out_arc = self.fix_arc_set()
        # super_source_out_arc.append(best_arc)
        # print("super_source_out_arc:",super_source_out_arc)
        # super_source_out_arc = self.fix_arc(self.data_number)

        while improved:
            print("Iteration :",iteration)
            best_acyclic_flow, improved, current_cost, l, q, h, best_arc = self.iterate_arc(best_acyclic_flow, improved, current_cost, l, q, h, best_arc, super_source_out_arc)
            # if (current_cost == 10e+14) or (self.solve_result !="solved"):
            #     print("\nSearch optimal acyclic network using ipopt solution")
            #     self.load_model()
            #     self.solve()
            #     print("Objective: ",self.total_cost)
            #     print("solve_result: ",self.solve_result)
            #     self.generate_random_acyclic_from_solution()
            #     # self.generate_random_acyclic_graph()
            #     self.update_model()
            #     self.solve()
            #     print("Objective: ",self.total_cost)
            #     print("solve_result: ",self.solve_result)
            #     self.best_ampl = self.ampl
            #     # self.display_results()
            #     self.iterate_acyclic_flows()

            # print("Current best acyclic network:")
            # plt.figure(figsize=(10, 8))
            # nx.draw_spectral(best_acyclic_flow, with_labels=True)
            # plt.show()
            # super_source_out_arc = self.fix_arc()
            # super_source_out_arc.append(best_arc)
            iteration += 1
            # print(f"Current best solution: {current_cost}")
            # print(" ")
        print("*********************************************************Final best results***************************************************************")
        print(f"Final best objective: {current_cost}")
        # print("length of the arcs: ", l, "\n")
        # print("flow in the arcs: ", q, "\n")
        # print("head value at node: ", h, "\n")
        # print("acyclicity_time:", acyclicity_time)
    
    # Function to suppress output
    @contextlib.contextmanager
    def suppress_output(self):
        # Open devnull to suppress the output
        with open(os.devnull, 'w') as devnull:
            # Redirect stdout and stderr
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            sys.stdout = devnull
            sys.stderr = devnull
            try:
                yield
            finally:
                # Restore original stdout and stderr
                sys.stdout = old_stdout
                sys.stderr = old_stderr   

    def solve(self):
        with self.suppress_output():
            """Solve the optimization problem."""
            self.ampl.option["solver"] = "ipopt"
            # self.ampl.option["solver"] = "/home/nitishdumoliya/minotaur/build/bin/mmultistart"
            self.ampl.set_option("ipopt_options", "outlev = 0  print_user_options = no  sb = yes ")
            self.ampl.set_option("knitro_options", "outlev = 0 ms_enable = 1 ms_maxsolves = 30 mip_multistart=1")
            self.ampl.set_option("baron_options","maxtime = 20  outlev = 1 lsolver=knitro firstloc 1 barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")
            self.ampl.option["presolve_eps"] = " 8.53e-15"
            self.ampl.option['presolve'] = 1
            # self.ampl.option['solver_msg'] = 0
            # self.ampl.option['show_stats'] = 0
            self.ampl.solve()
            self.solve_result = self.ampl.solve_result
            # Get the total cost
            self.total_cost = self.ampl.get_objective("total_cost").value()
            # print("Objective:", self.total_cost)
            # self.best_flow = list(self.network_graph.edges)  # Store the best flow configuration
            # print("solve_result: ",self.solve_result)
        # Get and display the solve time
        solve_time = self.ampl.get_value('_solve_elapsed_time')
        # print(f"Solve time: {solve_time} seconds")
        self.solver_time += solve_time
        
    def run(self):
        """Main method to run the optimization process."""
        start_time = time.time()
        self.load_model()
        # self.solve()
        # self.display_results()
        # print("Objective: ",self.total_cost)
        # print("solve_result: ",self.solve_result)
        # self.create_digraph()
        # self.plot_graph()
        # self.generate_random_acyclic_from_solution()
        # print("Acyclic network:",self.network_graph.edges(), "\n")
        self.generate_random_acyclic_graph()
        self.update_model()
        self.solve()
        
        # Get the solver results
        print("Objective: ",self.total_cost)
        print("solve_result: ",self.solve_result)
        # self.display_results()
        self.iterate_acyclic_flows()
        # self.display_best_results(self.best_ampl)
        elapsed_time = time.time() - start_time
    
        print("Solver_time:",self.solver_time)
        print("Heuristic elapsed time:", elapsed_time)

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
    data_number = 1
    input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data_list[data_number]}.dat"
    print("Water Network:", data_list[data_number],"\n")
    optimizer = WaterNetworkOptimizer("../m1Basic.mod", input_data_file, data_number)
    optimizer.run()
