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


#<k> = 2*K/N -> Remember that the network should be fully connected for the case of ER random graphs
N = 500 # nodes
K = 2000 # edges
m = 4    # number of stubs for the Barabasi-Albert random graph

if int(sys.argv[1]) == 0:
    #Creating an ER random graphs with N number of nodes and K number of edges
    G = nx.gnm_random_graph(N, K)

    #Creating a Barabasi-Albert with m stubs on each step
    #G=nx.barabasi_albert_graph(N,m) 

    #Creating the dictionary structure from the graph G
    graph_structure = graph_to_dict(G)

    #Computing the degree distribution from the dictionary structure of a graph
    #compute_degree_distribution(graph_structure)

    #Printing the network on a file
    print_network_onfile(graph_structure,'network_structure_ER.json')
    ranking=betweeness_fromG(G)
    #influence_distr=alpha_from_rank_betweeness(graph_structure,ranking,float(sys.argv[2])/(1.0*N))
    #print(influence_distr)
    influence_distr=alpha_from_rank(graph_structure,ranking,float(sys.argv[2]))
    print_influence_distribution(influence_distr,'influence_distribution.json')


    #Loading the network as a dictionary structure from a file
else:
    graph_structure=load_network_fromfile("network_structure_ER.json")
    #np.random([])
    #print(graph_structure)
    G=create_graph_fromedgelist(graph_structure)
    ranking=betweeness_fromG(G)
    print(float(sys.argv[2])/(1.0*100))
    influence_distr=alpha_from_rank(graph_structure,ranking,float(sys.argv[2])/(1.0*100))
    #influence_distr=alpha_from_rank_betweeness(graph_structure,ranking,float(sys.argv[2])/(1.0*N))
    print_influence_distribution(influence_distr,'influence_distribution.json')
#print(graph_structure)