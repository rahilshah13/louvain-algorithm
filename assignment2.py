import networkx as nx
import community
import math
import time

files = [
    ("2013-coactivation-matrix.txt", 1), 
    ("1993-macaque71.txt", 1), 
    ("1991-cerebral-cortex-fv30.txt", 1), 
    ("2007-macaque47.txt", 1), 
    ("1995-cerebral-cortex-cat.txt", 1), 
    ("2012-interareal-macaque.txt", 2)
]

communities = {}
# for creating graphs 1 - 5
def createGraphWeighted(fn, g):
    with open(fn, "r") as f:
        f.readline()                
        for line in f.readlines():
            d = line.strip().split("\t")
            g.add_edge(d[0], d[1], weight=float(d[2]))

# for creating graph 6
def createGraphWeightedDirected(fn, g):
    with open(fn, "r") as f:
        f.readline()     
        for line in f.readlines():
            d = line.strip().split("\t")
            g.add_edge(d[2], d[3], weight=-math.log(float(d[4])))

######################################
def goesToCommunity(g, e, c):
    return g.nodes[e[0]]["community"] == c or g.nodes[e[1]]["community"] == c

def inCommunity(g, e, c):
    return g.nodes[e[0]]["community"] == c and g.nodes[e[1]]["community"] == c

# dQ = (1/m) * (k_in - ( (sigma_in + k_in) *  k_i / 2m) )
def getDeltaQ(g, node, new_community):
    m = g.size(weight="weight")
    k_i = sum([g[u][v]["weight"] for u, v in g.edges(node)]) 
    k_in = sum([g[u][v]["weight"] for u, v in g.edges(node) if goesToCommunity(g, (u, v), new_community)])
    sigma_in = sum([g[u][v]["weight"] for u, v in g.edges(node) if inCommunity(g, (u, v), new_community)])
    return (1/m) * (k_in - ((sigma_in + k_in) *  k_i / (2*m)))

# comms: c -> {nodes}
def getQ(g, comms):
    m, total = g.size(weight="weight"), 0
    for c in comms:
        for i in comms[c]:
            for j in g.neighbors(i):
                if g.nodes[j]["community"] == g.nodes[i]["community"]:
                    total += (g[i][j]["weight"] - ((g.degree(i, weight="weight")*g.degree(j, weight="weight"))/ (2*m)))
    return total / (2*m)

def collapseGraph(g, comms):
    collapsed, visited = nx.Graph(), set()

    for c in comms:

        if len(comms[c]) == 0: 
            continue

        for i in comms[c]:
            if not g.has_node(i):
                continue
            for j in g.neighbors(i):
                if g.nodes[j]["community"] == g.nodes[i]["community"]:
                    if (i, j) not in visited:
                        visited |= {(i, j), (j, i)}
                        if collapsed.has_edge(c, c):
                            collapsed[c][c]["weight"] += g[i][j]["weight"]
                        else: 
                            collapsed.add_edge(c, c, weight=g[i][j]["weight"])
                else:
                    if collapsed.has_edge(c, j):
                        collapsed[c][j]["weight"] += g[i][j]["weight"]
                    else:
                        collapsed.add_edge(c, g.nodes[j]["community"], weight=g[i][j]["weight"])
                        collapsed.nodes[g.nodes[j]["community"]]["community"] = g.nodes[j]["community"]
        if collapsed.has_node(c):
            collapsed.nodes[c]["community"] = c

    return collapsed

def updateGraphCommunities(OG, comms):
    for c in comms:
        for i in comms[c]:
            OG.nodes[i]["community"] = c
    return OG
'''
OG is the original graph
g is used for the heirarchical decompositon
output format: [(modularity, #_communities), ...] in level order
'''
def louvainMethod(g, minC):
    global communities
    OG, results, old_n_comm = g.copy(), [], len(communities)
    results.append((getQ(OG, communities), old_n_comm))

    while True:
        g_new, comm_new = g.copy(), communities.copy()

        for node in g.nodes():
            dqResults, visited = [], set([g.nodes[node]["community"]])

            for neighbor in g.neighbors(node):
                c = g.nodes[neighbor]["community"]
                if c not in visited:
                    dqResults.append((getDeltaQ(g, node, c), c))
                visited.add(c)

            best = max(dqResults)
            # update node's community
            if best[0] > 0:
                for e in comm_new[g_new.nodes[node]["community"]]:
                    comm_new[best[1]].add(e)
                comm_new[g_new.nodes[node]["community"]].clear()
                g_new.nodes[node]["community"] = best[1]


        # get modularity and collapse graph
        OG = updateGraphCommunities(OG, comm_new)
        results.append((getQ(OG, comm_new), len([c for c in comm_new if len(comm_new[c]) != 0])))
        g, communities = collapseGraph(g_new, comm_new), comm_new

        #print("modularity g: ", results[-1][0], ", num comms: ", results[-1][1])

        if (old_n_comm - results[-1][1]) <= 0 or results[-1][1] <= minC:
            return results

        old_n_comm = results[-1][1]

# main loop
for file, type in files:
    G = nx.Graph()
    communities = dict()

    if type == 1: 
        createGraphWeighted(file, G)
    elif type == 2: 
        createGraphWeightedDirected(file, G)

    # initialize each node to its own community, Q = 0.0
    for node in G.nodes():
        G.nodes[node]["community"] = node
        communities[node] = set([node])   

    # calculate results
    start = time.time()
    results, output = louvainMethod(G, 3), ""
    bestMod, bestL = max(results)[0], 0
    for i, r in enumerate(results):
        if r[0] == bestMod:
            bestL = i
        output += "\tlevel {}: Q={}, n_comms={}\n".format(i, r[0], r[1])

    # package results
    pkg_comms = community.best_partition(G, weight='weight')
    pkg_mod = community.modularity(pkg_comms, G, weight='weight')
    pkg_comms = [pkg_comms[c] for c in pkg_comms if pkg_comms[c] != 0]
    pkg_avg = sum(pkg_comms) / len(pkg_comms)

    # print results
    print('File {} - {} nodes and {} edges\n\
    my results:\n{}\n\
    \tbest result: Q={} at level {}\n\
    \ttime: {} seconds\n\n\
    pkg results:\n\
    \tQ={}, n_comms={}, avg_size={}\
    '.format(file, G.number_of_nodes(), G.number_of_edges(), output, bestMod, bestL, time.time()-start, pkg_mod, len(pkg_comms), pkg_avg))
    print('--------------------------------------')