import networkx as nx
import time
from amplpy import AMPL
import matplotlib.pyplot as plt

class WaterNetworkOptimizer:
    def __init__(self, model_file, data_file):
        self.ampl = AMPL()
        self.model_file = model_file
        self.data_file = data_file
        self.total_cost = None
        self.network_graph = None
        self.best_flow = None

    def load_model(self):
        """Load the model and data."""
        self.ampl.reset()
        self.ampl.read(self.model_file)
        self.ampl.read_data(self.data_file)

    def build_graph(self):
        """Build a directed graph from the AMPL model."""
        nodes_list = [i for i in self.ampl.getSet('nodes')]
        # edges_list = self.ampl.getParameter('L').getValues()
        edges_list = self.ampl.getSet('arcs').to_list()
    
        self.network_graph = nx.DiGraph()
        self.network_graph.add_nodes_from(nodes_list)
        self.network_graph.add_edges_from(edges_list)
        print(nodes_list)
        print(edges_list)

    def plot_graph(self):
        nx.draw(self.network_graph, with_labels=True)
        plt.show()
     
    def solve(self):
        """Solve the optimization problem."""
        self.ampl.option["solver"] = "ipopt"
        self.ampl.set_option("ipopt_options", "outlev = 1")
        self.ampl.solve()

        # Get the total cost
        self.total_cost = self.ampl.get_objective("total_cost").value()
        print("Objective:", self.total_cost)

    def display_results(self):
        """Display relevant results from the optimization."""
        self.ampl.eval("display {(i,j) in arcs}: q[i,j];")
        self.ampl.eval("display h;")
        self.ampl.eval("display {(i,j) in arcs} h[i] - h[j];")
        self.ampl.eval("display {i in nodes} h[i] - (E[i] + P[i]);")

    def iterate_aclyclic_flows(self):
        """Iterate to find improved acyclic flows."""
        improved = True
       
        while improved:
            improved = False
            current_cost = self.total_cost
            
            for u, v in self.network_graph.edges():
                # Attempt to reverse the direction of the arc (u, v) to (v, u)
                self.network_graph.remove_edge(u, v)
                self.network_graph.add_edge(v, u, weight=self.network_graph[v][u]['weight'])  # Reversing the arc
                
                # Check if the new graph is still acyclic
                if nx.is_directed_acyclic_graph(self.network_graph):
                    self.ampl.eval("set arcs := " + str(self.network_graph.edges()))
                    self.solve()
                
                    if self.total_cost < current_cost:
                        current_cost = self.total_cost
                        self.best_flow = self.network_graph.edges(data=True)
                        improved = True
                        print("Improved cost found:", current_cost)
                
                # Restore the original direction of the arc
                self.network_graph.remove_edge(v, u)
                self.network_graph.add_edge(u, v, weight=self.network_graph[u][v]['weight'])  # Restore original arc direction

    def run(self):
        """Main method to run the optimization process."""
        start_time = time.time()
        self.load_model()
        self.build_graph()
        self.plot_graph()
        self.solve()
        self.display_results()
        #self.iterate_aclyclic_flows()
        elapsed_time = time.time() - start_time
        print("Elapsed time:", elapsed_time)
        print("Best flow configuration:", self.best_flow)

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
    data_number = 0
    input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data_list[data_number]}.dat"
    optimizer = WaterNetworkOptimizer("../m1Basic.mod", input_data_file)
    optimizer.run()

