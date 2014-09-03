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
Script to calculate nodes' usages of links using the all new and improved model
to rule all previous models.

In time this script should incorporate three transmission paradigms:
- localised / linear / selfish
- synchronised / square / cooperative
- market

The output should contain solutions for:
- import usages
- export usages
- combined usages

"""

"""
Initialisation
"""
if len(sys.argv)<2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

def link_label(num,N):
    """
    Translate link number into a string like: 'Country1-Country2'.
    """
    label = get_link_direction(num,N)
    return str(label[0].label)+'-'+str(label[1].label)

def bin_maker(F_matrix,q,lapse,nbins=None):
    """
    Create a lot of bins to calculate each nodes average usage of a link. The
    last bin should be placed such that its right side is aligned with the 99%
    quantile of the flow on the link. This value is set as the capacity of the
    link.
    """
    bin_max = np.ceil(q) # last bin ends at 99% quantile
    if not nbins:
        nbins = 12 # this number is not at all arbitrary
    bin_size = bin_max/nbins
    bin_edges = np.linspace(bin_size,bin_max,nbins) # value at the right side of each bin, last bin is 99% quantile

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

def node_contrib(H,bin_edges):
    """
    Calculate a node's contribution to a specific links capacity
    """
    flows = np.append([0],bin_edges)
    bin_size = flows[2]-flows[1]
    nbins = len(bin_edges)
    C = 0 # total contribution
    for i in range(nbins-1):
        c1,c2 = 0,0 # partial contributions
        if i == 0:
            c1 += (flows[i+1]-flows[i])
        else:
            c1 += (flows[i+1]-flows[i])/(1-bin_CDF(i,H))
        l = i+1
        while l < nbins:
            c2 += bin_prob(l,H)*H[l,1]/flows[l]
            l += 1
        C += c1*c2
    return C

def convergence(test):
    """
    Calculate usages and save to file
    """
    Node_contributions = np.zeros((len(N),len(F))) # empty array for calculated usages
    for node in range(len(N)):
        for link in range(len(F)):
            # Stacking and sorting data
            F_vert = np.reshape(F[link,:lapse],(len(F[link,:lapse]),1))
            exp_vert = np.reshape(Usages[link,node,:lapse],(len(Usages[link,node,:lapse]),1))
            F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
            F_matrix[F_matrix[:,0].argsort()]
            H,bin_edges = bin_maker(F_matrix,quantiles[link],lapse,nbins=test)
            Node_contributions[node,link] = node_contrib(H,bin_edges)
    # save results to file 
    np.save('./convergence/Node_contrib_'+str(test)+'.npy',Node_contributions)


if 'converge' in task:
    """
    Check for convergence with altering bin size
    """
    lapse = 70128 # number of hours to include
    test_bins = range(10,31,2) # how many bins to test
    print('Geting ready to test convergence')
    print('Lapse: '+str(lapse)+', Number of runs: '+str(len(test_bins)))


    
    # Pick a transmission paradigm
    N = EU_Nodes_usage('linear.npz')
    F = abs(np.load('./results/linear-flows.npy'))
    Usages = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy')
    # Get 99% quantile of link flow
    quantiles = [get_q(abs(F[link,:lapse]),.99) for link in range(len(F))]
    np.save('./convergence/quantiles.npy',quantiles)

    # Calculate usages
    print('Populating workers')
    p = Pool(6)
    p.map(convergence,test_bins)

    # Compare results
    print('Comparing results')
    mins = np.zeros(len(test_bins))
    means = np.zeros(len(test_bins))
    maxs = np.zeros(len(test_bins))
    quantiles = np.load('./convergence/quantiles.npy')
    index = 0
    for i in test_bins:
        N_usages = np.load('./convergence/Node_contrib_'+str(i)+'.npy')
        Total = np.sum(N_usages,0)/quantiles
        mins[index] = Total.min()
        means[index] = Total.mean()
        maxs[index] = Total.max()
        index+=1

    # Plot comparison
    print('Plotting results')
    np.savez('./convergence/results.npz',mins,means,maxs)
    print('Saved results to ./convergence/results.npz')
    plt.figure()
    plt.plot(test_bins,mins,'-r',lw=2)
    plt.plot(test_bins,means,'-k',lw=2)
    plt.plot(test_bins,maxs,'-b',lw=2)
    plt.xlabel('Number of bins')
    plt.ylabel(r'$C/T^{99\%}$')
    plt.title('Convergence check for different bin sizes')
    plt.savefig('./figures/convergence.png')
    print('Saved results to ./figures/convergence.png')

if 'solve' in task:
    """
    Calculate nodes' contributions and save results to file.
    """
    print('Solving')
    lapse = 70128 # number of hours to include
    
    # pick one of three transmission paradigms
    N = EU_Nodes_usage('linear.npz')
    F = abs(np.load('./results/linear-flows.npy'))
    
    # Do everything below for both import and export usages
    direction = ['import','export']
    # for loop over direction
    Usages = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy') # change linear, export
    
    # Calculate usages and save to file
    Node_contributions = np.zeros((len(N),len(F))) # empty array for calculated usages
    # Get 99% quantile of link flow
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
            Node_contributions[node,link] = node_contrib(H,bin_edges)
            
    # save results to file 
    np.save('Node_contrib_linear_export.npy',Node_contributions)
    print('Saved Node_contributions to Node_contrib_linear_export.npy')
    np.save('quantiles.npy',quantiles)
    print('Saved 99% quantiles to quantiles.npy')


if 'plot' in task:
    """
    Create various plots of usage and save figures to ./figures/
    """
    print('Getting ready to plot')

    # Load data and results
    N = EU_Nodes_usage('linear.npz')
    F = abs(np.load('./results/linear-flows.npy'))
    N_usages = np.load('Node_contrib_linear_export.npy')
    quantiles = np.load('quantiles.npy')
    link_dic = link_dict(N,F)
    links = range(len(F))
    nodes = range(len(N))
    names = [] # array for node labels

    # For every node, plot usage of every link, C_n
    print('Plotting node usages')
    for node in nodes:
        node_id = str(N[node].label)
        names.append(node_id)
        plt.figure()
        ax = plt.subplot(111)
        plt.bar(links,np.divide(N_usages[node,:],quantiles),width=1,color="#0000ff")
        # highlight directly connected links
        for link in link_dic[node_id]:
            plt.bar(link,np.divide(N_usages[node,link],quantiles[link]),width=1,color="#ff0000")
        ax.set_title(node_id+' usage of links')
        ax.set_xlabel('Link')
        ax.set_ylabel(r'Usage, $C_n$')
        plt.savefig('./figures/link-usage-'+str(node)+'.png')
    print('Saved figures to ./figures/link-usage-node.png')

    # For every link, plot each node's usage, C_n
    print('Plotting usage of each link')
    node_ids = np.array(range(len(N))).reshape((len(N),1))
    node_mean_load = [n.mean for n in N]
    node_mean_load = np.reshape(node_mean_load,(len(node_mean_load),1))
    for link in links:
        label = link_label(link,N)
        # sort node names for x-axis
        Usages = np.reshape(N_usages[:,link],(len(N_usages),1))
        data = np.hstack([Usages,node_ids,node_mean_load])
        data_sort = data[data[:,2].argsort()]
        names_sort = [names[int(i)] for i in data_sort[:,1]]
        # flip order so largest is first
        names_sort = names_sort[::-1]
        data_sort = data_sort[:,0][::-1]
        
        plt.figure()
        ax = plt.subplot(111)
        plt.bar(nodes,np.divide(data_sort,quantiles[link]),width=1,color="#0000ff")
        ax.set_xticks(np.linspace(1,len(N)+1,len(N)+1))
        ax.set_xticklabels(names_sort,rotation=60,ha="right",va="top")
        ax.set_title('Usage of link '+label)
        ax.set_ylabel(r'Usage, $C_n$')
        plt.savefig('./figures/node-usage-'+str(link)+'.png')
    print('Saved figures to ./figures/node-usage-link.png')

    # Compare node transmission to mean load
    print('Plotting node comparison')
    # sort node names for x-axis
    Total_usage = np.sum(N_usages,1)
    node_mean_load = [n.mean for n in N]
    normed_usage = Total_usage/node_mean_load
    normed_usage = np.reshape(normed_usage,(len(normed_usage),1))
    node_mean_load = np.reshape(node_mean_load,(len(node_mean_load),1))
    data = np.hstack([normed_usage,node_ids,node_mean_load])
    data_sort = data[data[:,2].argsort()]
    # flip order so largest is first
    data_sort = data_sort[:,0][::-1]

    plt.figure()
    ax = plt.subplot(111)
    plt.bar(nodes,data_sort,width=1,color="#0000ff")
    ax.set_xticks(np.linspace(1,len(N)+1,len(N)+1))
    ax.set_xticklabels(names_sort,rotation=60,ha="right",va="top")
    ax.set_title('Comparison of total network usage')
    ax.set_ylabel(r'[MW$_T$/MW$_L$]')
    plt.axis([0,len(N),0,1.05*max(data_sort)])
    plt.savefig('./figures/network-usage.png')
    print('Saved figures to ./figures/network-usage.png')
