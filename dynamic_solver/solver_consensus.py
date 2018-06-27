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
from ConsensusSolver import *


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--network_file', type=str, 
        help="File of the network", default=None)
    parser.add_argument('-i', '--influence_file', type=str, help="")
    parser.add_argument('-m', '--model', type=str, 
        help="Name for the model used (influence)", default="sigmoid")
    parser.add_argument('-p', '--param', type=float)
    parser.add_argument('-s', '--sample_size', type=int)
    args = parser.parse_args(arguments)
    
    edge_list = np.loadtxt(args.edge_list)
    alpha = np.loadtxt(args.influence_file)
    if args.model=="sigmoid":
        



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
