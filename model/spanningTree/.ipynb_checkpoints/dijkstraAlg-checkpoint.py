import networkx as nx

from amplpy import AMPL

data_list=[
        # "data1",
        # "twoloop",
        # "d1_Sample_input_cycle_twoloop",
        "d2_Sample_input_cycle_hanoi",
        # "d3_Sample_input_double_hanoi",
        # "d4_Sample_input_triple_hanoi",
        # "d5_Taichung_input",
        # "d6_HG_SP_1_4",
        # "d7_HG_SP_2_3",
        # "d8_HG_SP_3_4",
        # "d9_HG_SP_4_2",
        # "d10_HG_SP_5_5",
        # "d11_HG_SP_6_3",
        # "d12",
        # "d13",
        # "d14_NewYork",
        # "d15_foss_poly_0",
        # "d16_foss_iron",
        # "d17_foss_poly_1",
        # "d18_pescara",
        # "d19_modena"
        ]

# Python program for Dijkstra's single source shortest path algorithm. The program is
# for adjacency matrix representation of the graph

class Graph():
 
    def __init__(self, vertices):
        self.V = vertices
        self.graph = [[0 for column in range(vertices)]
                      for row in range(vertices)]
 
    def printSolution(self, dist):
        print("Vertex \t Distance from Source")
        for node in range(self.V):
            print(node, "\t\t", dist[node])
 
    # A utility function to find the vertex with minimum distance value, from the set of 
    # vertices not yet included in shortest path tree
    def minDistance(self, dist, sptSet):
 
        # Initialize minimum distance for next node
        min = 1e7
 
        # Search not nearest vertex not in the shortest path tree
        for v in range(self.V):
            if dist[v] < min and sptSet[v] == False:
                min = dist[v]
                min_index = v
 
        return min_index
 
    # Function that implements Dijkstra's single sourcen shortest path algorithm for a graph
    # represented using adjacency matrix representation
    def dijkstra(self, src):
 
        dist = [1e7] * self.V
        dist[src] = 0
        sptSet = [False] * self.V
 
        for cout in range(self.V):
 
            # Pick the minimum distance vertex from the set of vertices not yet processed.
            # u is always equal to src in first iteration
            u = self.minDistance(dist, sptSet)
 
            # Put the minimum distance vertex in the
            # shortest path tree
            sptSet[u] = True
 
            # Update dist value of the adjacent vertices of the picked vertex only if the 
            # current distance is greater than new distance and the vertex in not in the 
            # shortest path tree
            for v in range(self.V):
                if (self.graph[u][v] > 0 and
                   sptSet[v] == False and
                   dist[v] > dist[u] + self.graph[u][v]):
                    dist[v] = dist[u] + self.graph[u][v]
 
        self.printSolution(dist)

def Ampl(data_list):
    ampl = AMPL()
    ampl.reset()
    # ampl.read("spanningTreeMultiple.mod")
    ampl.read("spanningTreeSingle.mod")
    # ampl.read("manualPath.mod")
    input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data_list[0]}.dat"
    ampl.read_data(input_data_file)
    return ampl

nodes_list = []

def distanceMatrix(ampl):
    nodeList = ampl.getSet('nodes').to_list()
    arcSet = ampl.getSet('arcs').to_list()
    L = ampl.getParameter('L').getValues().toDict()
    print(arcSet)

    n = len(nodeList)
    
    # Initialize the distance matrix 
    dist = [[float('inf') for column in range(n)] 
                for row in range(n)]  
    print(dist[0][0])
    # Set the diagonal to 0, as the distance from a node to itself is 0
    for i in range(n):
        dist[i][i] = 0
    
    for (i,j) in arcSet:
        dist[i-1][j-1] = L[i,j]
        dist[j-1][i-1] = L[i,j]

    return dist


ampl = Ampl(data_list)
distance_matrix = distanceMatrix(ampl)
# print(distance_matrix)

g= Graph(32)
g.graph = distance_matrix
g.dijkstra(0)
 

'''

TREE_COST = []
LOOP_COST = []
for data in data_list:
    print("Water Network File : ",data)
    print("  ")

    print("#************************ Results Minimum Cost Spanning Tree********************#")    
    print(" ")

    tree_ampl = AMPL()
    tree_ampl.reset()
    # tree_ampl.read("spanningTreeMultiple.mod")
    tree_ampl.read("spanningTreeSingle.mod")
    # tree_ampl.read("manualPath.mod")
    input_data_file = f"/home/nitishdumoliya/waterNetwork/data/{data}.dat"
    tree_ampl.read_data(input_data_file)
    nodes_list = []

    for i in tree_ampl.getSet('nodes'):
        nodes_list.append(i)
    # print("Nodes :",nodes_list)
    edge_set = tree_ampl.getSet('arcs')

    edges_list = tree_ampl.getParameter('L').getValues()
    G = edges_list.toDict()
    # print(G)
    for v in sorted(G.items(), key = lambda item:item[1]):
        print(v)
    # {k: v for k, v in sorted(G.items(), key=lambda item: item[1])}
    uwg = nx.Graph()
    uwg.add_nodes_from(nodes_list)
    # print(uwg.nodes())

    uwg.add_weighted_edges_from(edges_list)
    # print(uwg.edges())

    ##################################################################################################
    
    # Calculate a minimum spanning tree of an undirected weighted graph with the kruskal algorithm
    mst = nx.minimum_spanning_tree(uwg, algorithm='kruskal')
    print(mst)
    # G = mst.edges(data=True)
    # print(G)
    # {k: v for k, v in sorted(G.items(), key=lambda item: item[1])}
    # print(sorted(mst.edges(data=True)))
    print("Minimum spanning tree is",mst.edges() )
    print(" ")
    for (i,j) in edge_set:
        if (i,j) in mst.edges():
            print(f"tree_ampl.eval(s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 1;)")
            tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 1;")
        elif (j,i) in mst.edges():
            tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 1;")
        else:
            print((i,j))
            tree_ampl.eval(f"s.t. fixed_binary_{i}_{j}: x[{i},{j}]= 0;")
            tree_ampl.eval(f"s.t. h_{i}: h[{i}] >= E[{i}] + P[{i}];")
            tree_ampl.eval(f"s.t. h_{j}: h[{j}] >= E[{j}] + P[{j}];")

    ####################################################################################################
    print("======================Solver Results====================")
    tree_ampl.option["solver"] = "cplexamp"
    tree_ampl.option["cplex_options"] = "iisfind=1 mipdisplay = 5"
    # tree_ampl.option["solver"] = "/home/nitishdumoliya/minotaur/build/bin/mmultistart"
    # tree_ampl.set_option("mmultistart_options","--presolve 1,--log_level 6,--eval_within_bnds 1")
    # tree_ampl.option["bonmin_options"] = "bonmin.bb_log_level 5 bonmin.nlp_log_level 0 "
    # tree_ampl.option["ipopt_options"] = " outlev = 0"
    tree_ampl.option["knitro_options"] = "outlev = 0 threads=1 feastol = 1.0e-7 feastol_abs = 1.0e-7 ms_enable = 1 ms_maxsolves = 10 ms_maxtime_real = 30 mip_multistart=1 maxtime_real =60"
    # tree_ampl.option["knitro_options"] = "outlev = 3 ms_enable 1  ms_maxsolves 2 "
    tree_ampl.option["presolve_eps"]="1.64e-12"
    # tree_ampl.set_option("baron_options","maxtime = -1  outlev = 1 lsolver=conopt  barstats deltaterm 1 objbound    threads = 12  prloc = 1 prfreq=1000 prtime 10")

    # tree_ampl.option["presolve"]="1"
    tree_ampl.solve()
    # tree_ampl.eval("show;")
    # tree_ampl.eval("expand;")

    # tree_ampl.eval("display X;")
    # tree_ampl.eval("display x;")
    tree_ampl.eval("display l;")
    #tree_ampl.eval("display {(i,j) in arcs} : sum{k in pipes}l[i,j,k]*C[k] ;")
    # tree_ampl.eval("display q;")
    tree_ampl.eval("display h;")
    # tree_ampl.eval("display {j in 1.._ncons} (_conname[j]);")
    # tree_ampl.eval("display x;")
    # tree_ampl.eval("display sum{i in nodes} h[i];")
    # tree_ampl.eval("display {(i,j) in arcs} h[i]-h[j];")
    # tree_ampl.eval("display sum{(i,j) in arcs} abs(h[i]-h[j]);")
    # tree_ampl.eval("display total_cost;")
    totalcost = tree_ampl.get_objective("total_cost")
    print("Objective for minimum spanning tree is:", totalcost.value())
    TREE_COST.append(totalcost.value())
    print("==========================================================")
    break




print("Tree Objective Value list :", TREE_COST)
print(" ")
print("Loop Objective Value list :", LOOP_COST)

'''
