#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solve the consensus dynamics for a given structure and set of parameters

Author: Guillaume St-Onge
"""

import os
import sys
import argparse
import numpy as np
from numpy.random import randint

def sigmoidal_influence(alpha_i, alpha_j):
    pass

def consensus_step(x, influence_model, alpha, network):
    size = len(network)
    #choose speaker-listener
    j = randint(0,size)
    i = network[randint(0,len(network[j]))]
    #update listener state
    x[i] += (x[j]-x[i])*influence_model(alpha[i],alpha[j])

def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--network_file', type=str, 
        help="File of the network")
    parser.add_argument('-i', '--influence_file', type=str, help="")
    parser.add_argument('-m', '--model', type=str, 
        help="Name for the model used")
    parser.add_argument('-d', '--delta', type=float, 
        help="maximal variation allowed per encounter")
    parser.add_argument('-p', '--params', type=float)

    args = parser.parse_args(arguments)





if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))