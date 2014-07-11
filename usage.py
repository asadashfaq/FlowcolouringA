#! /usr/bin/env/python
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
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
lapse = 365*24


"""
Helper functions
"""
def usage(m):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """
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
    boxplot,boxplotlabel = track_link_usage_total(N2,F,mode=m,alph='same',lapse=lapse)
    return

def scatter_plotter(F,Fmax,usage,direction,mode,lapse):
    """
    Scatter plots of nodes' import/export usages of links saved to ./figures/.
    """
    colors = ['#b30000','#ff0000','#ff4d4d']
    nodes = usage.shape[1]
    links = usage.shape[0]
    for l in range(links):
        plt.figure()
        ax = plt.subplot(111)
        for t in range(lapse):
            linkflow = abs(F[l,t]*np.ones(nodes))
            usages = usage[l,:,t]/Fmax[l]
            """
            Find the three highest contributors and plot in red colors.
            """
            maxs = []
            for k in range(3):
                i = usages.argmax()
                maxs.append(usages[i])
                usages = np.delete(usages,i)
                plt.plot(linkflow[k],maxs[k],'.',color=str(colors[k]))
            """
            Plot other nodes' usages.
            """
            plt.plot(linkflow[3:],usages,'.k')
        ax.set_title(str(mode)+' '+str(direction)+' flows on link #'+str(l+1))
        ax.set_xlabel(r'$F_l(t)$ [MW]')
        ax.set_ylabel(r'$H_{ln}/\max(F_l)$')

        # Shrink x-axis to make room for legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.85, box.height])

        names = ['1st', '2nd', '3rd', 'Rest']
        ax.legend((names),loc='center left', bbox_to_anchor=(1,0.5),title='Contributions')

        plt.savefig('./figures/'+str(mode)+'-'+str(direction)+'-flows-'+str(l+1)+'.png')
    return

def top_plotter(top,N,F,Fmax,usage,direction,mode,lapse):
    """
    Plots of top 5 contributors to each links transmission capacity
    """
    colors = ['#e32636','#008000','#007fff','#fdee00','#ed9121']
    nodes = usage.shape[1]
    links = usage.shape[0]
    for l in range(links):
        plt.figure()
        ax = plt.subplot(111)
        max_usages = np.max(usage[l],1)
        max_users = []
        names = ['Rest']
        for k in range(top):
            i = max_usages.argmax()
            max_users.append(i)
            max_usages = np.delete(max_usages,i)
            names.append(str(N[i].label))
        for t in range(lapse):
            linkflow = abs(F[l,t]*np.ones(nodes))
            usages = usage[l,:,t]/Fmax[l]
            bg, = ax.plot(linkflow,usages,'.',color='#b2beb5',alpha=.3,label='rest')
            for k in range(top):
                ax.plot(linkflow[k],usages[max_users[k]],'.',color=str(colors[k]),label=names[k])
        ax.set_title('Top '+str(top)+' '+str(mode)+' '+str(direction)+' flows on link #'+str(l+1))
        ax.set_xlabel(r'$F_l(t)$ [MW]')
        ax.set_ylabel(r'$H_{ln}/\max(F_l)$')

        # Shrink x-axis to make room for legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.85, box.height])

        ax.legend((names),loc='center left', bbox_to_anchor=(1,0.5),title='Country')
        plt.savefig('./figures/top-'+str(mode)+'-'+str(direction)+'-flows-'+str(l+1)+'.png')
    return

"""
Solving flows for different export schemes
"""""
if 'solve' in task:
    alphas = np.ones(30)*.7
    gammas = np.ones(30)
    
    Nodes = EU_Nodes_usage()
    
    for m in modes:
        N,F = au.solve(Nodes,mode=m+' copper verbose',lapse=lapse)
        N.save_nodes(m)
        np.save('./results/'+m+'-flows',F)


"""
Calculate powermixes and nodes' usages of links and save results to file.
"""
if 'usage' in task:
    p = Pool(3)
    p.map(usage,modes)


"""
Scatterplots of usages
"""
if 'plot' in task:
    top = 5
    for m in modes:
        """
        Load flows and find maximum values for normalisation
        """
        N = EU_Nodes_usage(m+'.npz')
        F = np.load('./results/'+m+'-flows.npy')
        Fmax = np.max(abs(F),1)

        """
        Load usages, do scatter plots and plot of top 5 contributors to a links capacity.
        """
        export_usage = np.load('./linkcolouring/old_'+m+'_copper_link_mix_export_all_alpha=same.npy')
        scatter_plotter(F,Fmax,export_usage,'export',m,lapse)
        top_plotter(top,N,F,Fmax,export_usage,'export',m,lapse)
        export_usage = []

        import_usage = np.load('./linkcolouring/old_'+m+'_copper_link_mix_import_all_alpha=same.npy')
        scatter_plotter(F,Fmax,import_usage,'import',m,lapse)
        top_plotter(top,N,F,Fmax,import_usage,'import',m,lapse)
        import_usage = []

