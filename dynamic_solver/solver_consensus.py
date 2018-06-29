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
from ConsensusSolver import *
from scipy.special import erfc
from scipy.stats import kstest


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
    parser.add_argument('-p', '--param', type=float, 
        help="Parameter (eta) for the influence model", default=0.5)
    parser.add_argument('-t', '--tol', type=float,
        help="Standard deviation for consensus", default=0.01)
    parser.add_argument('-s', '--sample_size', type=int,
        help="Number of sample", default=100)
    parser.add_argument('--allout', action='store_true', 
        help="Output the result for each consensus")
    args = parser.parse_args(arguments)

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
            influence_dict = {int(k):v for k,v in influence_dict.items()}

    #solve the linear influence model
    if args.model=="linear":
        #initialize solver and result
        S = ConsensusSolver(network_dict, influence_dict, eta=args.param)
        result = np.zeros((args.sample_size,2))
        #get a sample of consensus formation
        for i in range(0,args.sample_size):
            t, x = S.reach_consensus(args.tol)
            result[i][0] = t
            result[i][1] = x
            S.reset_state()

    #solve the sigmoidal influence model
    elif args.model=="sigmoidal":
        #initialize solver and result
        S = ConsensusSolver(network_dict, influence_dict, 
            influence_function=sigmoidal_influence(args.param))
        result = np.zeros((args.sample_size,2))
        #get a sample of consensus formation
        for i in range(0,args.sample_size):
            t, x = S.reach_consensus(args.tol)
            result[i][0] = t
            result[i][1] = x
            S.reset_state()

    if args.allout == True:
        for sample in result:
            print(sample[0], sample[1])
    else:
        #output mean time to consensus and standard deviation of the consensus opinion
        print(np.mean(result[:,0]), democratic_fairness(len(network_dict),result[:,1]))
        



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
