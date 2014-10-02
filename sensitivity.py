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
Script to calculate nodes' usage of links in the same way as
usage_old.py and usage.py but for different networks: N = 8, N = 83.

This script incorporates three transmission paradigms referred to as _modes_ below:
- localised / linear / selfish
- synchronised / square / cooperative
- market / RND

The output contain solutions for:
- import usages
- export usages
- combined import/export usages

All ouput files are saved to ./sensitivity/

Call the script using only one of the following command line arguments:
- solve: solve network and save solutions


"""

if len(sys.argv)<2:
    raise Exception('Not enough inputs')
else:
    task = str(sys.argv[1:])

networks = ['regions', 'superRegions']
schemes = ['linear', 'square', 'RND']
lapse = 70128  # number of hours to include
    
def networkSolver(network):
    print('Solving :'+str(network))
    if network == 'regions':
        Nodes = EU_Nodes_regions()
    elif network == 'superRegions':
        Nodes = EU_Nodes_superRegions()

    for scheme in schemes:
        print 'Solving: '+scheme
        N,F = au.solve(Nodes,mode=scheme+' copper',lapse=lapse)
        N.save_nodes(scheme,path='./sensitivity/'+network+'-')
        np.save('./sensitivity/'+network+'-'+scheme+'-flows',F)

def calc_usage(network):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """
    for scheme in schemes:
        if network == 'regions':
            N = EU_Nodes_regions('../sensitivity/regions-'+scheme+'.npz')
            F = np.load('./sensitivity/regions-'+scheme+'-flows.npy')
        elif network == 'superRegions':
            N = EU_Nodes_superRegions('../sensitivity/superRegions-'+scheme+'.npz')
            F = np.load('./sensitivity/superRegions-'+scheme+'-flows.npy')
    
        """
        N2 is a new nodes object containing individual powermixes for import and
        export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
        """
        N2,power_mixes_total = track_flows(N,F,lapse=lapse)
    
        """
        track_link_usage_total tracks each nodes usage of all links. The results
        are saved to files '..._links_ex_...' and '..._links_im_...'.
        """
        boxplot,boxplotlabel = track_link_usage_total(N2,F,mode=scheme,alph='same',lapse=lapse,sensitivity=True)
        return

"""
Solving flows for different export schemes
"""
if 'solve' in task:
    print 'Mode selected: Solving network flows'
    print 'Lapse: '+str(lapse)+' hours'
    
    p = Pool(2)
    p.map(networkSolver,networks)


"""
Calculate powermixes and nodes' usages of links and save results to file.
"""
if 'usage' in task:
    print 'Mode selected: Calculate usages.'
    print 'Lapse: '+str(lapse)+' hours.'
    print 'Running on 3 cores.'
    p = Pool(2)
    p.map(calc_usage,networks)

