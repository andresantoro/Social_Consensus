#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve the consensus dynamics for a given structure and set of parameters.

Author: Guillaume St-Onge
"""

import os
import sys
import argparse
import numpy as np
import json
from scipy.integrate import quad
from scipy.special import erfc
from scipy.stats import kstest

#Try to import fast module, otherwise rely on python
cpp_module_installed = True
try:
    from FastConsensusSolver import QuenchedSolver
    from FastConsensusSolver import AnnealedSolver
except ImportError:
    cpp_module_installed = False

def null_model_distribution(N):
    return lambda x: N**2*0.5*np.sum([(-1)**k*(x*N-k)**(N-1)*np.sign(x*N-k)/(gamma(k+1)*gamma(N-k+1))
        for k in range(0,N+1)])

def asymptotic_null_model_distribution(N):
    sigma = 1/np.sqrt(12*N)
    return lambda x: np.exp(-(x-0.5)**2/(2*sigma**2))/np.sqrt(2*np.pi*sigma**2)

def null_model_cumulative(N):
    return lambda x: quad(null_model_distribution(N), 0, x)[0]

def asymptotic_null_model_cumulative(N):
    return lambda x: 0.5*erfc(np.sqrt(6*N)*(0.5-x))

def democratic_fairness(N, output_sample):
    #use the kolmogorov-smirnov test as a measure of fairness (0 = most fair)
    return np.abs(kstest(output_sample, asymptotic_null_model_cumulative(N))[0])

def democratic_pointwise_k(input_sample,output_sample,k):
		sum_dk=0
		for i in range(0,len(output_sample)):
			sum_dk += abs(input_sample[i]-output_sample[i])**k
		return sum_dk/len(input_sample)

def sigmoidal_influence(eta):
    return lambda alpha_i, alpha_j : 1/(1+np.exp(6*(alpha_i-alpha_j+0.5)))

def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--network_file', type=str,
        help="File of the network")
    parser.add_argument('-i', '--influence_file', type=str,
        help="File for the influence of each node", default = None)
    parser.add_argument('-m', '--model', type=str,
        help="Name for the model used (influence)", default="linear")
    parser.add_argument('-mm', '--mode', type=str,
        help="Mode of the network", choices=['quenched', 'annealed'],
        default='quenched')
    parser.add_argument('-p', '--param', type=float,
        help="Parameter (eta) for the influence model", default=0.5)
    parser.add_argument('-t', '--tol', type=float,
        help="Standard deviation for consensus", default=0.01)
    parser.add_argument('-s', '--sample_size', type=int,
        help="Number of sample", default=100)
    parser.add_argument('--seed', type=int, help='Seed for the RNG', default=42)
    parser.add_argument('--allout', action='store_true',
        help="Output all history vector for each sample")
    parser.add_argument('-o', '--output_file', type=str,
        help="JSON File to output all history", default="allout.json")
    args = parser.parse_args(arguments)

    #verify the installation of cpp module
    if not cpp_module_installed:
        print("The cpp solver module needs to be installed properly")

    else:
        #Import the network
        with open(args.network_file, 'r') as fp:
            network_dict = json.load(fp)
            #convert key items to integer
            network_dict = {int(k):v for k,v in network_dict.items()}

        #Import the influence of each node if specified
        if args.influence_file == None:
            influence_dict = {k : 0.5 for k,v in network_dict.items()}
        else:
            with open(args.influence_file, 'r') as fp:
                influence_dict = json.load(fp)
                #convert key items to integer
                influence_dict = {int(k):float(v) for k,v in influence_dict.items()}

        #initialize output dictionary for allout and vectors
        if args.allout:
            output_dict = dict()
            for i in range(args.sample_size):
                output_dict["initial_state"] = dict()
                output_dict["history_vector"] = dict()
        mean_final_state = np.zeros(args.sample_size)
        consensus_time = np.zeros(args.sample_size)


        #initialize linear influence model
        if args.model=="linear":
            #initialize solver
            if args.mode=="quenched":
                S = QuenchedSolver(network_dict, influence_dict,
                                            args.param, args.seed)
            elif args.mode=="annealed":
                #priority dict is constructed based on the degree of nodes
                priority_dict = {node: len(neighborhood) for node, neighborhood
                                in network_dict.items()}
                S = AnnealedSolver(priority_dict, influence_dict,
                                            args.param, args.seed)

        #initialize k-fairness
        k = 2
        k_fairness=np.zeros(args.sample_size)

        #get a sample of consensus formation
        for i in range(0,args.sample_size):
            S.reach_consensus(args.tol)
            k_fairness[i]=democratic_pointwise_k(S.get_initial_state_vector(),
                                                 S.get_state_vector(),k)
            mean_final_state[i] = S.get_mean()
            consensus_time[i] = S.get_time()
            #output initial state and all variations
            if args.allout == True:
                output_dict["initial_state"][i] = \
                        S.get_initial_state_vector()
                output_dict["history_vector"][i] = \
                        S.get_history_vector()
            S.reset_all()

        #output the data
        if args.allout:
            with open(args.output_file, "w") as wf:
                json.dump(output_dict, wf)
        else:
            #Output average mesures
            print(np.mean(consensus_time), np.std(consensus_time),
                  np.mean(k_fairness), np.std(k_fairness),
                  democratic_fairness(len(network_dict), mean_final_state))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
