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
Script to create various plots of nodes' usages of links for three export schemes:
- Linear (most localised)
- Square (cooperative)
- Market (RND)

"""


"""
Initialisation
"""
if len(sys.argv)<2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

modes = ['linear','square','RND']   # Export schemes. 'RND' is also known as 'Market'.


"""
Helper functions
"""
def link_label(num,N):
    """
    Translate link number into a string like: 'Country1-Country2'.
    """
    label = get_link_direction(num,N)
    return str(label[0].label)+'-'+str(label[1].label)


if 'country' in task:
    """
    Figures of how important each country is for the flow on each link
    """
    print "Calculating weighted sum of each country's conditional usage."

    # Load data
    N = EU_Nodes_usage('linear.npz')
    F = abs(np.load('./results/linear-flows.npy'))
    export_usage = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy')
    N_usages = np.zeros((len(N),len(F))) # empty array for weighted sums of conditional usages
    links = range(len(F))
    names = []

    for node in range(len(N)):
        names.append(str(N[node].label))
        print node
        for link in range(len(F)):
            print link
            # Get 99% quantile of link flow and determine bin size
            qq = get_q(abs(F[link,:lapse]),.99) # Get 99% quantile of link flow to determine bin size
            bin_size = .04*qq  # percentage is determined in the _convergence_ mode.
            # Stacking and sorting data
            F_vert = np.reshape(F[link,:lapse],(len(F[link,:lapse]),1))
            exp_vert = np.reshape(export_usage[link,node,:lapse],(len(export_usage[link,node,:lapse]),1))
            F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
            F_matrix.sort(axis=0)
            
            # Calculate weigthed sums of usages and save for later use
            M,bin_sum = bin_maker(bin_size,F_matrix,summed=1)
            N_usages[node,link] = bin_sum
    
    # save results to file for faster plotting
    np.save('N_usages.npy',N_usages)
    print 'saved N_usages to N_usages.npy'
    usage_sums = np.sum(N_usages,0) # total usage of each link

    # plot each nodes usage of each link normed to the total usage of the link
    for node in range(len(N)):
        plt.figure()
        ax = plt.subplot(111)

        plt.bar(links,np.divide(N_usages[node,:],usage_sums),width=1)
        
        node_id = str(N[node].label)
        ax.set_title(node_id+' usage of links normalised to total individual link usage')
        ax.set_xlabel(r'Link')
        ax.set_ylabel(r'$U_{nl}/\sum_n U_{nl}$')
        plt.savefig('./figures/country-importance-'+str(node)+'.png')
        plt.close('all')
    print 'Saved results to: ./figures/country-importance-node.png'

    # plot each nodes usage of each link normed to the country's total network usage
    print "Calculating link importance for country's network usage"
    network_usages = np.sum(N_usages,1) # each country's total network usage
    for node in [0]: # range(len(N)):
        plt.figure()
        ax = plt.subplot(111)

        plt.bar(links,np.divide(N_usages[node,:],network_usages[node]),width=1)
        
        node_id = str(N[node].label)
        ax.set_title(node_id+' usage of links normalised to total network usage')
        ax.set_xlabel(r'Link')
        ax.set_ylabel(r'$U_{nl}/\sum_l U_{nl}$')
        plt.savefig('./figures/link-importance-'+str(node)+'.png')
        plt.close('all')
    print 'Saved results to: ./figures/link-importance-node.png'

    # compare each nodes total usage of the network to the others normalised to the average usage of the network
    avg_usage = sum(network_usages)/len(network_usages)
    plt.figure()
    ax = plt.subplot(111)

    nodes = range(len(N))
    plt.bar(nodes,np.divide(network_usages,avg_usage),width=1)
    plt.plot([0,30],[1,1],'--k',lw=2)

    ax.set_xticks(np.linspace(0,len(N),len(N)+1))
    ax.set_xticklabels(names,rotation=60,ha="right",va="top")
    
    ax.set_title('Comparing network usages for weighted sums of link usages')
    ax.set_xlabel(r'Node')
    ax.set_ylabel(r'Network usage normalised to avg. network usage')
    plt.savefig('./figures/network-usage.png')
    plt.close('all')
    print 'Saved results to: ./figures/network-usage.png'
