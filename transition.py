#! /usr/bin/env/python
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import aurespf.solvers as au
from aurespf.tools import get_q, get_quant_caps
#from EUgrid import *
from link_colour_less import track_flows, get_link_direction
from new_linkcolouralgorithm_less import track_link_usage_total
from link_namer import node_namer,link_dict
#from usage import bin_maker, bin_prob, bin_CDF, node_contrib
from europe_plusgrid import *

"""
Script to investigate the transition of a node's usage of the links as
the transmission is constrained

"""

def calc_usage(N,F,lapse,direction,b):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """
    
    """
    N2 is a new nodes object containing individual powermixes for import and
    export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
    """
    N2,power_mixes_total = track_flows(N,F,lapse=lapse)

    """
    track_link_usage_total tracks each nodes usage of all links. The results
    are saved to files '..._links_ex_...' and '..._links_im_...'.
    """
    boxplot,boxplotlabel = track_link_usage_total(N2,F,mode=direction,lapse=lapse,constrained=b)
    return

def caller(direction):
    for b in np.linspace(0.05,1.5,30):
        F = np.load('./ConstrainedFlowData/Europe_aHE_'+str(b)+'q99_DC_'+direction+'_flows.npy')
        calc_usage(N,F,lapse,direction,b)

lapse = 8760 # 280512
N = europe_plus_Nodes(load_filename='../ConstrainedFlowData/Europe_aHE_copper_DC_lin.npz')
directions = ['lin', 'sqr']
p = Pool(2)
p.map(caller,directions)
