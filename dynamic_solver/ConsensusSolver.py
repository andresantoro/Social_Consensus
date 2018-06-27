#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solver class for the social consensus dynamics

Author: Guillaume St-Onge
"""

import os
import sys
import argparse
import numpy as np
from numpy.random import randint
from numpy.random import random

def sigmoidal_influence(eta):
    return lambda alpha_i, alpha_j : 1/(1+np.exp(eta*(alpha_i-alpha_j)))

class ConsensusSolver:

    def __init__(self, 
        network_dict, #dictionary that fixes the structure
        influence_list=None, #list/array containing the influence of each node
        eta=None, #inverse temperature for the sigmoidal_influence fonction
        influence_function=None, #can parse any influence function
        x0 = None): #initial state to parse -- otherwise x0 \in [0,1]^N
    
        #initialize solver
        if influence_list == None:
            self.influence_list = random(len(network_dict))
        else:
            self.influence_list = influence_list
        self.network_dict = network_dict
        if x0 == None:
            self.x = random(len(network_dict))
        else:
            self.x = x0
        if influence_function != None:
            self.influence_function = influence_function
        else:
            if eta != None:
                self.influence_function = sigmoidal_influence(eta)
            else:
                self.influence_function = sigmoidal_influence(eta=1)

    def consensus_step(self):
        #choose listener
        i = randint(0,len(self.network_dict))
        j = self.network_dict[1][randint(0,len(self.network_dict[i]))]
        
        #update listener state
        self.x[i] += (self.x[j]-self.x[i])*self.influence_function(
            self.influence_list[i], self.influence_list[j])

    def reset_state(self):
        self.x = random(len(self.network_dict))

    def generate_random_state(self):
        return random(len(self.network_dict))

    def reach_consensus(self, tol):
        t = 0
        while np.std(self.x) > tol:
            self.consensus_step()
            t += 1
        return t, np.mean(self.x)



if __name__ == '__main__':
    #test of the class
    network_dict = {0: [1,2,3], 1: [0,2,4], 2: [0,1], 3: [0,4], 4: [3,1]}
    S = ConsensusSolver(network_dict)
    t, x = S.reach_consensus(0.01)
    print(t,x)
        