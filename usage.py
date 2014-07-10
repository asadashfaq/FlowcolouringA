#! /usr/bin/env/python
from __future__ import division
import sys
import numpy as np
from pylab import plt
import aurespf.solvers as au
from aurespf.tools import *
from EUgrid import *
from link_colour_less import track_flows
from new_linkcolouralgorithm_less import track_link_usage_total

"""
Script to create scatter plots of nodes' usages of links for three export schemes:
- Linear (most localised)
- Square (cooperative)
- Market (RND)

Call the program with only one of the following command line arguments:
- solve : solves the network
- usage : calculates nodes' usages of links
- plot : plot results

Example call:
python usage.py solve

"""

"""
Initialisation
"""
if len(sys.argv)<2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1])

modes = ['linear','square','RND']
lapse = 10

"""
Solving flows for different export schemes
"""
if 'solve' in task:
    alphas = np.ones(30)*.7
    gammas = np.ones(30)
    
    Nodes = EU_Nodes_usage()
    
    for m in modes:
        N,F = au.solve(Nodes,mode=m+' copper verbose',lapse=lapse)
        N.save_nodes(m)
        np.save('./results/'+m+'-flows',F)

"""
Calculate powermixes
"""
if 'usage' in task:
    for m in modes:
        N = EU_Nodes_usage(m+'.npz')
        F = np.load('./results/'+m+'-flows.npy')
        
        """
        N2 is a new nodes object containing individual powermixes for import and
        export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
        """
        N2,power_mixes_total = track_flows(N,F,lapse=lapse)

        """
        track_link_usage_total tracks each nodes usage of all links. The results
        are saved to files '..._links_ex_...' and '..._links_im_...'.
        """
        boxplot,boxplotlabel = track_link_usage_total(N2,F,mode=m,lapse=lapse)

"""
Scatterplots of usages
"""
if 'plot' in task:
    for m in modes:
        # do something



