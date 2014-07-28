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
Script to print link IDs  and labels.
"""

def link_label(num,N):
    """
    Translate link number into a string like: 'Country1-Country2'.
    """
    label = get_link_direction(num,N)
    return str(label[0].label)+'-'+str(label[1].label)

def link_namer(N=None,F=None):
    """
    Return a list of link labels: 'from-to'
    """
    if N == None:
        N = EU_Nodes_usage('linear.npz')
    if F == None:
        F = np.load('./results/linear-flows.npy')
    lnames = []
    for link in range(len(F)):
        lnames.append(link_label(link,N))
    return lnames

def node_namer(N=None):
    """
    Return a list of node labels eg. 'DK'
    """
    if N == None:
        N = EU_Nodes_usage('linear.npz')
    nnames = []
    for node in N:
        nnames.append(node.label)
    return nnames

def link_dict(N=None,F=None):
    """
    Return a dictionary of node labels and IDs of their directly connected links
    eg. 'AT':[0,1,2,3,4,5]
    """
    if N == None:
        N = EU_Nodes_usage('linear.npz')
    if F == None:
        F = np.load('./results/linear-flows.npy')
    ldict = {}
    for link in range(len(F)):
        label = get_link_direction(link,N)
        for l in label:
            if str(l.label) in ldict:
                ldict[str(l.label)].append(link)
            else:
                ldict.update({str(l.label):[link]})
    return ldict
