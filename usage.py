#! /usr/bin/env/python
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import aurespf.solvers as au
from aurespf.tools import get_q, get_quant_caps
from EUgrid import *
from link_colour_less import track_flows, get_link_direction
from new_linkcolouralgorithm_less import track_link_usage_total
from link_namer import node_namer,link_dict

"""
Script to calculate nodes' usages (stakes) of links using the all new and
improved model to rule all previous models.

In time this script should incorporate three transmission paradigms:
- localised / linear / selfish
- synchronised / square / cooperative
- market

The output should contain solutions for:
- import stakes
- export stakes
- combined stakes

"""

"""
Initialisation
"""
if len(sys.argv)<2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

def bin_maker(F_matrix,q,lapse):
    """
    Create a lot of bins to calculate each nodes average usage of a link. The
    last bin should be placed such that its right side is aligned with the 99%
    quantile of the flow on the link. This value is set as the capacity of the
    link.
    """
    bin_max = np.ceil(q) # last bin ends at 99% quantile
    nbins = 12 # this number is not at all arbitrary
    bin_size = bin_max/nbins
    bin_edges = np.linspace(bin_size,bin_max),nbins) # value at the right side of each bin, last bin is 99% quantile

    H_temp = []
    H = np.zeros((nbins,2)) # [number of events, mean usage]
    for b in range(nbins-1):
        for t in range(lapse):
            if b*bin_size <= F_matrix[t,0] < (b+1)*bin_size:
                H_temp.append(F_matrix[t,1])
        if len(H_temp)>0:
            H[b,0] = len(H_temp)
            H[b,1] = sum(H_temp)/len(H_temp)
        else: # no data in the bin
            H[b,0] = 0
            H[b,1] = 0
        H_temp=[]

    # move events further to the right of the 99% quantile to the last bin and keep their value
    b = nbins-1
    for t in range(lapse):
        if b*bin_size <= F_matrix[t,0]:
            H_temp.append(F_matrix[t,1])
    if len(H_temp)>0:
        H[b,0] = len(H_temp)
        H[b,1] = sum(H_temp)/len(H_temp)
    else: # no data in the bin
        H[b,0] = 0
        H[b,1] = 0
    H_temp=[]
    return H,bin_edges

def bin_prob(bin_id,H):
    """
    Probability of occurence of given flow
    """
    return H[bin_id,0]/sum(H[:,0])

def bin_CDF(bin_id,H):
    """
    Cumulative distribution function to be used on the bins created in the
    bin_maker.
    """
    i = 0
    P = 0
    while i <= bin_id:
        P += H[i,0]
        i += 1
    P = P/sum(H[:,0])
    return P

def node_stake(H,bin_edges):
    """
    Calculate a node's stake in a specific link
    """
    flows = np.append([0],bin_edges)
    nbins = len(bin_edges)
    c = 0 # node stake
    for i in range(nbins-1):
        c += (flows[i+1]-flows[i])/(1-bin_CDF(i,H))
        for j in range(i+1,nbins):
            c += bin_prob(j,H)*H[j,1]
    return c

if 'solve' in task:
    """
    Calculate nodes' stakes and save results to file.
    """
    lapse = 1000 # number of hours to include
    
    # pick one of three transmission paradigms
    N = EU_Nodes_usage('linear.npz')
    F = abs(np.load('./results/linear-flows.npy'))
    
    # Do everything below for both import and export usages
    direction = ['import','export']
    # for loop over direction
    Usages = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy') # change linear, export
    
    # Calculate usages and save to file
    Node_stakes = np.zeros((len(N),len(F))) # empty array for calculated usages
    # Get 99% quantile of link flow to determine bin size
    quantiles = [get_q(abs(F[link,:lapse]),.99) for link in range(len(F))]

    for node in range(len(N)):
        print node+1,'/ '+str(len(N))
        for link in range(len(F)):
            # Stacking and sorting data
            F_vert = np.reshape(F[link,:lapse],(len(F[link,:lapse]),1))
            exp_vert = np.reshape(Usages[link,node,:lapse],(len(Usages[link,node,:lapse]),1))
            F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
            F_matrix[F_matrix[:,0].argsort()]
            
            H,bin_edges = bin_maker(F_matrix,quantiles[link],lapse)
            Node_stakes[node,link] = node_stake(H,bin_edges)
            
    # save results to file for faster and better plotting in usage_plotting.py
    np.save('Node_stakes_linear_export.npy',Node_stakes)
    print 'Saved Node_stakes to Node_stakes_linear_export.npy'
    np.save('quantiles.npy',quantiles)
    print 'Saved 99% quantiles to quantiles.npy'


# plotting mode
# calculate double integral
# collect usages to compare all nodes

# aggregate results from import and export usages and compare all nodes


