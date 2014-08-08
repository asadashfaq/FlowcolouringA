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

def bin_maker(F_matrix,99q):
    """
    Create a lot of bins to calculate each nodes average usage of a link. The
    last bin should be placed such that its right side is aligned with the 99%
    quantile of the flow on the link. This value is set as the capacity of the
    link.
    """
    bin_max = np.ceil(99q) # last bin ends at 99% quantile
    nbins = 12 # this number is not at all arbitrary
    bin_size = bin_max/nbins
    bin_means = np.linspace(.5*bin_size,bin_max-(.5*bin_size),nbins) # mean values of bins

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

    # move events further to the right to the last bin
    b = 11
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

    # save bin sums and number of events


#    if summed:
#        part_sum = np.multiply(bin_means,bin_size)
#        bin_sum = sum(np.multiply(H[:,1],part_sum))
#        return np.array([bin_means,H[:,1]]),bin_sum
#    else:
#        return bin_means,H


def bin_CDF:(flow,bins)
    """
    Cumulative distribution function to be used on the bins created in the
    bin_maker.
    """

    # input a certain flow and a bin series
    # sum the bins until the given flow is reached


lapse = 10 # number of hours to include

# pick one of three transmission paradigms
N = EU_Nodes_usage('linear.npz')
F = abs(np.load('./results/linear-flows.npy'))

# Do everything below for both import and export usages
direction = ['import','export']
# for loop over direction
Usages = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy') # change linear, export

# Calculate usages and save to file
Node_usages = np.zeros((len(N),len(F))) # empty array for calculated usages
export_usage = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy') # change linear
for node in [0]: # range(len(N)):
    print node+1,'/ '+str(len(N))
    for link in [0]: #range(len(F)):
        # Get 99% quantile of link flow and determine bin size
        qq = get_q(abs(F[link,:lapse]),.99) # Get 99% quantile of link flow to determine bin size
        bin_size = .08*qq  # percentage is determined in the _convergence_ mode.
        # Stacking and sorting data
        F_vert = np.reshape(F[link,:lapse],(len(F[link,:lapse]),1))
        exp_vert = np.reshape(export_usage[link,node,:lapse],(len(export_usage[link,node,:lapse]),1))
        F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
        F_matrix[F_matrix[:,0].argsort()]
        
        # Calculate weigthed sums of usages and save for later use
        M,bin_sum = bin_maker(bin_size,F_matrix,summed=1)
        Node_usages[node,link] = bin_sum

# save results to file for faster and better plotting in usage_plotting.py
np.save('Node_usages_linear_export.npy',Node_usages)
print 'Saved N_usages to N_usages.npy'



# plotting mode
# calculate double integral
# collect usages to compare all nodes

# aggregate results from import and export usages and compare all nodes

