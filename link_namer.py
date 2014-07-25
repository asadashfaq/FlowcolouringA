#! /usr/bin/env/python
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import aurespf.solvers as au
from aurespf.tools import get_q
from aurespf.tools import *
from EUgrid import *
from link_colour_less import track_flows, get_link_direction
from new_linkcolouralgorithm_less import track_link_usage_total

"""
Script to print link names and IDs in the format: ID 'from country' - 'to country'
eg: 0 AT-CH
"""

def link_label(num,N):
    """
    Translate link number into a string like: 'Country1-Country2'.
    """
    label = get_link_direction(num,N)
    return str(label[0].label)+'-'+str(label[1].label)


N = EU_Nodes_usage('linear.npz')
F = np.load('./results/linear-flows.npy')

for link in range(len(F)):
    print link,link_label(link,N)
