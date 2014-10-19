#! /usr/bin/env/python
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from aurespf.tools import get_q, get_quant_caps
from link_colour_less import track_flows, get_link_direction
from new_linkcolouralgorithm_less import track_link_usage_total
from link_namer import node_namer, link_namer, link_dict
from functions import binMaker, bin_prob, bin_CDF, node_contrib
from europe_plusgrid import *
from EUgrid import EU_Nodes

"""
Script to investigate the transition of a node's usage of the links as
the transmission is constrained

All ouput files are saved to ./constrained/

Call the script using only one of the following command line arguments:
- trace:    run flow tracing and save results
- usage:    calculate node's usage of links
- plot:		create figures of usage as function of transmission constraint
"""

if len(sys.argv)<2:
    raise Exception('Not enough inputs')
else:
    task = str(sys.argv[1:])

def calc_usage(N,F,lapse,scheme,b):
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
    boxplot,boxplotlabel = track_link_usage_total(N2,F,mode=scheme,lapse=lapse,constrained=b)
    return

def calc_contribution(b, verbose=False):
    """
    Calculate nodes' contribution for a single scheme and save results to file.
    Set verbose to True to follow the progress in the command line. Not recommended for multi processing.
    """
    for scheme in schemes:
        if verbose:
            print('scheme: '+scheme) 
        F = abs(np.load('./ConstrainedFlowData/Europe_aHE_'+str(b)+'q99_DC_'+scheme+'_flows.npy'))
        # Get 99% quantile of link flow
        quantiles = [get_q(abs(F[link,:lapse]),.99) for link in range(len(F))]
        
        # Do everything below for both import and export usages unless we want the combined case
        for direction in directions:
            if verbose:
                print('Direction: '+direction)
            if direction == 'combined':
                Usages1 = np.load('./constrained/linkcolouring/one_year/'+scheme+'_link_mix_import_b_'+str(b)+'.npy')
                Usages2 = np.load('./constrained/linkcolouring/one_year/'+scheme+'_link_mix_export_b_'+str(b)+'.npy')
                Usages = Usages1+Usages2
                Usages1,Usages2 = None,None
            else:
                Usages = np.load('./constrained/linkcolouring/one_year/'+scheme+'_link_mix_'+str(direction)+'_b_'+str(b)+'.npy')
            print('Loaded '+scheme+' usages')
            
            # Calculate usages and save to file
            Node_contributions = np.zeros((nNodes,len(F))) # empty array for calculated usages
            for node in range(nNodes):
                if verbose:
                    print node+1,'/ '+str(nNodes)
                for link in range(len(F)):
                    # Stacking and sorting data
                    F_vert = np.reshape(F[link,:lapse],(len(F[link,:lapse]),1))
                    exp_vert = np.reshape(Usages[link,node,:lapse],(len(Usages[link,node,:lapse]),1))
                    F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
                    F_matrix[F_matrix[:,0].argsort()]
                    
                    H,bin_edges = binMaker(F_matrix,quantiles[link],lapse,N_bins)
                    Node_contributions[node,link] = node_contrib(H,bin_edges)
                    
            # save results to file 
            np.save('./constrained/Node_contrib_'+scheme+'_'+direction+'_b_'+str(b)+'.npy',Node_contributions)
            print('Saved Node_contributions to ./constrained/Node_contrib_'+scheme+'_'+direction+'_b_'+str(b)+'.npy')
            if direction == 'import':
                np.save('./constrained/quantiles_'+scheme+'_b_'+str(b)+'.npy',quantiles)
                print('Saved 99% quantiles to ./constrained/quantiles_'+scheme+'_b_'+str(b)+'.npy')

def caller(scheme):
    N = europe_plus_Nodes(load_filename='../ConstrainedFlowData/Europe_aHE_copper_DC_lin.npz')
#    for b in np.linspace(0.05,1.5,30):
    for b in np.linspace(0.05,1.45,15):
        F = np.load('./ConstrainedFlowData/Europe_aHE_'+str(b)+'q99_DC_'+scheme+'_flows.npy')
        calc_usage(N,F,lapse,scheme,b)

def plotter():
    """
    Make a lot of figures for a single mode and put them in
    ./constrained/figures/scheme/
    This function is very ugly because it was written very fast without thinking much
    """
    # Country's usage of all links as function of b
    for node in nodes:
        plt.figure(figsize=(16,7))
        usages = []
        totalUsagesLin = []
        scheme='lin'
        for b in np.linspace(0.05,1.5,30):
            # Load data and results
            F = abs(np.load('./ConstrainedFlowData/Europe_aHE_'+str(b)+'q99_DC_'+scheme+'_flows.npy'))
            quantiles = np.load('./constrained/quantiles_'+scheme+'_b_'+str(b)+'.npy')
            N_usages = np.load('./constrained/Node_contrib_'+scheme+'_combined_b_'+str(b)+'.npy')
            usages.append(.5*N_usages[node]/quantiles)
            totalUsagesLin.append(np.sum(.5*N_usages/quantiles,1))
        usages = np.array(usages).transpose()
        ax1 = plt.subplot(121)
        plt.pcolormesh(usages)
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        plt.ylabel('link number')
        plt.title(names[node]+'\'s usage of links, linear')
        ax1.set_xticks(np.linspace(0,30,16))
        ax1.set_xticklabels(np.linspace(0,1.5,16))
        ax1.set_yticks(np.linspace(.5,49.5,50))
        ax1.set_yticklabels(link_names,ha="right",va="center",fontsize=6)

        usages = []
        totalUsagesSqr = []
        scheme='sqr'
        for b in np.linspace(0.05,1.5,30):
            # Load data and results
            F = abs(np.load('./ConstrainedFlowData/Europe_aHE_'+str(b)+'q99_DC_'+scheme+'_flows.npy'))
            quantiles = np.load('./constrained/quantiles_'+scheme+'_b_'+str(b)+'.npy')
            N_usages = np.load('./constrained/Node_contrib_'+scheme+'_combined_b_'+str(b)+'.npy')
            usages.append(.5*N_usages[node]/quantiles)
            totalUsagesSqr.append(np.sum(.5*N_usages/quantiles,1))
        usages = np.array(usages).transpose()
        ax2 = plt.subplot(122)
        plt.pcolormesh(usages)
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        plt.title(names[node]+'\'s usage of links, square')
        ax2.set_xticks(np.linspace(0,30,16))
        ax2.set_xticklabels(np.linspace(0,1.5,16))
        ax2.set_yticks(np.linspace(.5,49.5,50))
        ax2.set_yticklabels(link_names,ha="right",va="center",fontsize=6)
        plt.savefig('./constrained/figures/node-usage-'+str(node)+'.png', bbox_inches='tight')

    # Figure comparing nodes' total network usage as function of b
    totalUsagesLin = np.array(totalUsagesLin).transpose()
    totalUsagesSqr = np.array(totalUsagesSqr).transpose()
    plt.figure(figsize=(16,7))
    ax3 = plt.subplot(121)
    ax3.set_xticks(np.linspace(0,30,16))
    ax3.set_xticklabels(np.linspace(0,1.5,16))
    ax3.set_yticks(np.linspace(.5,29.5,30))
    ax3.set_yticklabels(names,ha="right",va="center",fontsize=9)
    plt.pcolormesh(totalUsagesLin)
    plt.colorbar()
    plt.xlabel(r'$\beta$')
    plt.title('Total network usage, linear')
    ax4 = plt.subplot(122)
    ax4.set_xticks(np.linspace(0,30,16))
    ax4.set_xticklabels(np.linspace(0,1.5,16))
    ax4.set_yticks(np.linspace(.5,29.5,30))
    ax4.set_yticklabels(names,ha="right",va="center",fontsize=9)
    plt.pcolormesh(totalUsagesSqr)
    plt.colorbar()
    plt.xlabel(r'$\beta$')
    plt.title('Total network usage, square')
    plt.savefig('./constrained/figures/total-network-usage.png', bbox_inches='tight')

    # All countries usage of a single link as function of b
    for link in links:
        plt.figure(figsize=(16,7))
        usages = []
        scheme='lin'
        for b in np.linspace(0.05,1.5,30):
            # Load data and results
            F = abs(np.load('./ConstrainedFlowData/Europe_aHE_'+str(b)+'q99_DC_'+scheme+'_flows.npy'))
            quantiles = np.load('./constrained/quantiles_'+scheme+'_b_'+str(b)+'.npy')
            N_usages = np.load('./constrained/Node_contrib_'+scheme+'_combined_b_'+str(b)+'.npy')
            usages.append(.5*N_usages[:,link]/quantiles[link])
        usages = np.array(usages).transpose()
        ax5 = plt.subplot(121)
        ax5.set_xticks(np.linspace(0,30,16))
        ax5.set_xticklabels(np.linspace(0,1.5,16))
        ax5.set_yticks(np.linspace(.5,29.5,30))
        ax5.set_yticklabels(names,ha="right",va="center",fontsize=9)
        plt.pcolormesh(usages)
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        plt.title('Usage of link '+link_names[link]+', linear')

        usages = []
        scheme='sqr'
        for b in np.linspace(0.05,1.5,30):
            # Load data and results
            F = abs(np.load('./ConstrainedFlowData/Europe_aHE_'+str(b)+'q99_DC_'+scheme+'_flows.npy'))
            quantiles = np.load('./constrained/quantiles_'+scheme+'_b_'+str(b)+'.npy')
            N_usages = np.load('./constrained/Node_contrib_'+scheme+'_combined_b_'+str(b)+'.npy')
            usages.append(.5*N_usages[:,link]/quantiles[link])
        usages = np.array(usages).transpose()
        ax6 = plt.subplot(122)
        ax6.set_xticks(np.linspace(0,30,16))
        ax6.set_xticklabels(np.linspace(0,1.5,16))
        ax6.set_yticks(np.linspace(.5,29.5,30))
        ax6.set_yticklabels(names,ha="right",va="center",fontsize=9)
        plt.pcolormesh(usages)
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        plt.title('Usage of link '+link_names[link]+', square')
        plt.savefig('./constrained/figures/link-usage-'+str(link)+'.png', bbox_inches='tight')

lapse = 8760 # 280512
schemes = ['lin', 'sqr']
directions = ['import', 'export', 'combined']

"""
Calculate powermixes and nodes' usages of links and save results to file.
"""
if 'trace' in task:
    print('Mode selected: flow tracing')
    p = Pool(2)
    p.map(caller,schemes)

"""
Calculate nodes' contributions and save results to file.
"""
if 'usage' in task:
    print('Mode selected: usage calculation')
    N = np.load('./ConstrainedFlowData/Europe_aHE_copper_DC_lin.npz', mmap_mode='r+')
    nNodes = len(N['id'])
    N_bins = 16
    b = np.linspace(0.05,1.5,30)
    p = Pool(8)
    p.map(calc_contribution, b)

"""
Create various plots of usage and save figures to ./figures/
"""
if 'plot' in task:
    N = EU_Nodes()
    F = abs(np.load('./ConstrainedFlowData/Europe_aHE_0.05q99_DC_lin_flows.npy'))
    link_dic = link_dict(N,F)
    nodes = range(len(N))
    links = range(len(F))
    names = [str(N[i].label) for i in range(len(N))]
    link_names = link_namer(N,F)
    N = None
    plotter()
