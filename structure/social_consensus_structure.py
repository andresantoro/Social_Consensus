
#Packages
import numpy as np
import networkx as nx
from matplotlib import pyplot
import json
import random
import sys

#Function that converts a graph in networkx format into a dictionary (nodes, first_neighbours)
def graph_to_dict(G):
    graph=dict()
    for i in G.edges:
        node_i = i[0]
        node_j = i[1]
        #Scanning over the edges and inserting the links in the dictionary
        if node_i not in graph:
            graph[node_i]=[node_j]
        else:
            graph[node_i].append(node_j)
        if node_j not in graph:
            graph[node_j]=[node_i]
        else:
            graph[node_j].append(node_i)        
    return graph

#Function that computes the degree distribution of a graph [input should be a dictionary]
def compute_degree_distribution(graph):
    dgr=np.zeros(max(max(graph.keys()),len(graph.keys())))
    for i in graph.keys():
        dgr[i-1]= len(graph[i])
    pyplot.hist(dgr,density=True)

#Function that prints the graph structure into a file [.json format] "filename"
def print_network_onfile(graph,filename):
    with open(filename, 'w+') as fp:
        json.dump(graph, fp)
        
#Function that loads the graph structure from the file "filename" into a dictionary
#with the format: [key:nodes,values:first_neighbours]
def load_network_fromfile(filename):
    json_data=open(filename).read()
    graph=json.loads(json_data)
    graph={int(k):[int(i) for i in v] for k,v in graph.items()}
    return(graph)

#Function that creates the networkx structure from the edge list dictionary. It returns the graph G as output
def create_graph_fromedgelist(graph_structure):
    G=nx.Graph()
    for key in graph_structure.keys():
        for l in graph_structure[key]:
            G.add_edge(key,l)
    return(G)

#Function that computes the betweeness centrality of a generic graph G.
#It returns a dictionary in the form {key:nodes, values:rankings}
def betweeness_fromG(G):
    btw=nx.betweenness_centrality(G)
    rank_btw=sorted(btw,key=btw.get, reverse=True)
    ranking_btw =dict()
    for node,rank in enumerate(rank_btw):
        ranking_btw[node]=rank
    return(ranking_btw)
    
#Function that assigns the influence (Leaders/Followers) to each node uniformly at random   
def alpha_from_rank(graph_structure,ranking,ratio= 0.5,alphaL=0.9,alphaF=0.1):
    N=len(graph_structure.keys())
    N_lead=int(round(N)*ratio)
    N_followers= N-N_lead
    permutation=np.arange(0,N)
    random.shuffle(permutation)
    lead_follow_dict=dict()
    for i,j in enumerate(permutation):
        if i < N_lead:
            lead_follow_dict[str(j)]=str(alphaL)
        else:
            lead_follow_dict[str(j)]=str(alphaF)
    #print(lead_follow_dict)
    return(lead_follow_dict)

#Function that assigns the influence (Leaders/Followers) to each node proportional to the betweeness centrality   
def alpha_from_rank_betweeness(graph_structure,ranking,ratio= 0.5,alphaL=0.9,alphaF=0.1):
    N=len(graph_structure.keys())
    N_lead=int(round(N)*ratio)
    N_followers= N-N_lead
    lead_follow_dict=dict()
    #first -> ranking, values -> node
    for i,j in ranking.items():
        if i < N_lead:
            lead_follow_dict[str(j)]=str(alphaL)
        else:
            lead_follow_dict[str(j)]=str(alphaF)
    return(lead_follow_dict)

def print_influence_distribution(influence_distr,filename):
    with open(filename, 'w+') as fp:
        json.dump(influence_distr, fp)

def create_power_law_degree_sequence(N,kmin,kmax,exponent):
    degree_vector = np.arange(kmin,kmax+1)
    degree_distribution = degree_vector**(-exponent)
    degree_distribution /= np.sum(degree_distribution)
    degree_distribution = np.concatenate((np.array([0]*int(kmin)), 
        degree_distribution))
    cumulative_degree_distribution = np.array([
        np.sum(degree_distribution[0:i]) 
        for i in range(len(degree_distribution)+1)])
    #generate n expected degree from the dist
    u_vec = np.random.random(N)
    expected_degree_sequence = [np.searchsorted(cumulative_degree_distribution,
    u, side = 'right')-1 for u in u_vec]

    return expected_degree_sequence


#<k> = 2*K/N -> Remember that the network should be fully connected for the case of ER random graphs
N = 100 # nodes
K = 400 # edges
m = 4    # number of stubs for the Barabasi-Albert random graph
kmin = 6 #minimal degree for power-law graph
kmax = np.floor(10*np.sqrt(N)) #maximal degree allowed
exponent = 2.25

#Parameters for the Watts-Strogatz random graph
k_neigh=8
p_WS=0.05


#Parameters for the planted partition graph
gamma = 2
avg_degree=8
N_in_G=50
N_modules=2
p_out= avg_degree/(1.0*N_in_G *(gamma+1))
p_in = p_out*gamma




if int(sys.argv[1]) == 0:
    #Creating an ER random graphs with N number of nodes and K number of edges
    if str(sys.argv[3]) == 'ER':
        G = nx.gnm_random_graph(N, K)

    #Creating a Barabasi-Albert with m stubs on each step
    if str(sys.argv[3]) == 'BA':
        G=nx.barabasi_albert_graph(N,m) 

    #Create random graph with a power-law degree distribution
    # G = nx.expected_degree_graph(create_power_law_degree_sequence(
    #     N,kmin,kmax,exponent), selfloops=False)

    #Create the Watts-Strogatz random graph
    if str(sys.argv[3]) == 'WS':
        G =nx.watts_strogatz_graph(N, k_neigh, p_WS, seed=None)
    #print(len(G.edges()))
    #Create the planted partion random graph // Number of groups,  Number of vertices in each group, prob. of connecting vertices within a group, prob. of connected vertices between groups
    if str(sys.argv[3]) == 'PP':
        G = nx.planted_partition_graph(N_modules,N_in_G, p_in, p_out,seed=95)
    

    #Creating the dictionary structure from the graph G
    graph_structure = graph_to_dict(G)

    #Computing the degree distribution from the dictionary structure of a graph
    #compute_degree_distribution(graph_structure)

    #Printing the network on a file
    string1= 'network_structure_{0}.json'.format(str(sys.argv[3]))
    print_network_onfile(graph_structure,string1)
    ranking=betweeness_fromG(G)
    influence_distr=alpha_from_rank_betweeness(graph_structure,ranking,float(sys.argv[2])/(1.0*N))
    #print(influence_distr)
    string2= 'influence_distribution_{0}.json'.format(str(sys.argv[3]))
    #influence_distr=alpha_from_rank(graph_structure,ranking,float(sys.argv[2]))
    print_influence_distribution(influence_distr,string2)


    #Loading the network as a dictionary structure from a file
else:
    string1= 'network_structure_{0}.json'.format(str(sys.argv[3]))
    graph_structure=load_network_fromfile(string1)
    #np.random([])
    #print(graph_structure)
    G=create_graph_fromedgelist(graph_structure)
    ranking=betweeness_fromG(G)
    print(float(sys.argv[2])/(1.0*100))
    #influence_distr=alpha_from_rank(graph_structure,ranking,float(sys.argv[2])/(1.0*100))
    influence_distr=alpha_from_rank_betweeness(graph_structure,ranking,float(sys.argv[2])/(1.0*N))
    string2= 'influence_distribution_{0}.json'.format(str(sys.argv[3]))
    print_influence_distribution(influence_distr,string2)
#print(graph_structure)