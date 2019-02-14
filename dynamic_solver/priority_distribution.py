#Packages
import numpy as np
import networkx as nx
from matplotlib import pyplot
import json
import random
import sys

#Function that assigns the influence (Leaders/Followers) to each node uniformly at random
def alpha_uniform(N,ranking,ratio= 0.5,alphaL=0.9,alphaF=0.1):
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
def alpha_from_rank(N,ranking,ratio= 0.5,alphaL=0.9,alphaF=0.1):
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

#Function that prints the priority distribution into a file [.json format] "filename"
def print_priority_onfile(priority_dict,filename):
    with open(filename, 'w+') as fp:
        json.dump(priority_dict, fp)


def create_centralized_priority_list(N, priority_min=5, exponent=0):
    index = np.arange(1,N+1)*1.
    ressource_proportion = index**(-exponent)
    ressource_proportion /= np.sum(ressource_proportion)
    priority_list = (np.ones(N)*priority_min +
                                (N-priority_min-1)*ressource_proportion)
    return priority_list


if __name__ == '__main__':
    if str(sys.argv[2]) == "CENTRALIZED":
        N = 200
        priority_min = 5
        exponent = 1

        if int(sys.argv[1]) == 0:
            #print only one time the distribution
            priority_list = sorted(create_centralized_priority_list(
                N, priority_min,exponent), reverse=True)
            priority_dict = {i: priority_list[i] for i in range(N)}

            #Printing the priority dictionary
            string1= 'priority_dict_{0}.json'.format(str(sys.argv[2]))
            print_priority_onfile(priority_dict,string1)

        #Print influence dist
        ranking = {i:i for i in range(N)} #node are sorted in priority order
        influence_dist=alpha_from_rank(N,ranking,float(sys.argv[3])/(1.0*N))
        string2= 'influence_distribution_{0}.json'.format(str(sys.argv[2]))
        print_influence_distribution(influence_dist,string2)

