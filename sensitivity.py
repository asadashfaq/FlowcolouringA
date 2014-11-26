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
from functions import binMaker, bin_prob, bin_CDF, node_contrib

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
- solve:    solve network and save solutions
- trace:    run flow tracing and save results
- usage:    calculate nodes' usage and save results

usage can be followed by 'length' to include modelling of link lengths.
"""

if len(sys.argv)<2:
    raise Exception('Not enough inputs')
else:
    task = str(sys.argv[1:])

networks = ['superRegions', 'regions']
schemes = ['linear', 'square'] #'RND'
directions = ['import', 'export', 'combined']
lapse = 70128  # number of hours to include

if 'length' in task: length = True # modelling of link lengths
else: length = False

def regionSolver(scheme):
    print('Solving regions')
    Nodes = EU_Nodes_regions()

    N,F = au.solve(Nodes,mode=scheme+' copper',lapse=lapse)
    N.save_nodes(scheme,path='./sensitivity/regions-')
    np.save('./sensitivity/regions-'+scheme+'-flows',F)

def sRegionSolver(scheme):
    print('Solving super regions')
    Nodes = EU_Nodes_superRegions()

    N,F = au.solve(Nodes,mode=scheme+' copper',lapse=lapse)
    N.save_nodes(scheme,path='./sensitivity/superRegions-')
    np.save('./sensitivity/superRegions-'+scheme+'-flows',F)

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
        if network == 'regions':
            admat = './settings/mergedadmat.txt'
        elif network == 'superRegions':
            admat = './settings/superregionadmat.txt'
        N2,power_mixes_total = track_flows(N,F,admat=admat,lapse=lapse)
    
        """
        track_link_usage_total tracks each nodes usage of all links. The results
        are saved to files '..._links_ex_...' and '..._links_im_...'.
        """
        boxplot,boxplotlabel = track_link_usage_total(N2,F,mode=scheme,alph='same',lapse=lapse,sensitivity=True)
    return

def calc_contribution(network, scheme, verbose=None):
    """
    Calculate nodes' contribution for a single mode and save results to file.
    Adapted from _solver_ in usage.py
    Set verbose to True to follow the progress in the command line. Not recommended for multi processing.
    """
    # pick one of three transmission paradigms
    if verbose:
        print('Scheme: '+scheme)
    
    if network == 'superRegions':
        N = EU_Nodes_superRegions('../sensitivity/superRegions-'+scheme+'.npz')
        if length: lengths = 'superRegions'
    elif network == 'regions':
        N = EU_Nodes_regions('../sensitivity/regions-'+scheme+'.npz')
        if length: lengths = 'regions'
    else:
        raise Exception('Wrong network!')
    
    # Get 99% quantile of link flow
    F = abs(np.load('./sensitivity/'+network+'-'+scheme+'-flows.npy'))
    quantiles = [get_q(abs(F[link,:lapse]),.99) for link in range(len(F))]
    
    # Do everything below for both import and export usages unless we want the combined case
    for direction in directions:
        if verbose:
            print('Direction: '+direction)
        if direction == 'combined':
            Usages1 = np.load('./sensitivity/linkcolouring/'+network+'-old_'+scheme+'_copper_link_mix_import_all_alpha=same.npy')
            Usages2 = np.load('./sensitivity/linkcolouring/'+network+'-old_'+scheme+'_copper_link_mix_export_all_alpha=same.npy')
            Usages = (Usages1 + Usages2)*.5
            Usages1,Usages2 = None,None
        else:
            Usages = np.load('./sensitivity/linkcolouring/'+network+'-old_'+scheme+'_copper_link_mix_'+direction+'_all_alpha=same.npy')
        print('Loaded '+scheme+' usages')
        
        # Calculate usages and save to file
        Node_contributions = np.zeros((len(N),len(F))) # empty array for calculated usages
        for node in range(len(N)):
            if verbose:
                print node+1,'/ '+str(len(N))
            for link in range(len(F)):
                # Stacking and sorting data
                F_vert = np.reshape(F[link,:lapse],(len(F[link,:lapse]),1))
                exp_vert = np.reshape(Usages[link,node,:lapse],(len(Usages[link,node,:lapse]),1))
                F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
                F_matrix[F_matrix[:,0].argsort()]
                
                H,bin_edges = binMaker(F_matrix,quantiles[link],lapse,N_bins)
                Node_contributions[node,link] = node_contrib(H, bin_edges, linkID=link, lengths=lengths)

        # save results to file
        if not length:
            np.save('./sensitivity/'+network+'-Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy',Node_contributions)
            print('Saved Node_contributions to ./sensitivity/'+network+'-Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')
        else:
            np.save('./sensitivity/'+network+'-Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy',Node_contributions)
            print('Saved Node_contributions to ./sensitivity/'+network+'-Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
        if direction == 'import':
            np.save('./sensitivity/'+network+'-quantiles_'+scheme+'_'+str(lapse)+'.npy',quantiles)
            print('Saved 99% quantiles to ./sensitivity/'+network+'-quantiles_'+scheme+'_'+str(lapse)+'.npy')

"""
Solving flows for different export schemes
"""
if 'solve' in task:
    print 'Mode selected: Solving network flows'
    print 'Lapse: '+str(lapse)+' hours'
    
    p = Pool(3)
    p.map(regionSolver,schemes)
    p.map(sRegionSolver,schemes)

"""
Calculate powermixes and nodes' usages of links and save results to file.
"""
if 'trace' in task:
    print 'Mode selected: Flowtraving.'
    print 'Lapse: '+str(lapse)+' hours.'
    for network in networks:
        calc_usage(network)

"""
Calculate nodes' contributions and save results to file.
"""
if 'usage' in task:
    print 'Mode selected: Calculate usage'
    N_bins = 90 # This number has been determined from a calculation of convergence. See 'convergence' in _usage.py_.
    for network in networks:
        for scheme in schemes:
            calc_contribution(network,scheme)
