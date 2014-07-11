#! /usr/bin/env/python
from __future__ import division
import sys
import numpy as np
from pylab import plt
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
lapse = 70128

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
Calculate powermixes and nodes' usages of links and save results to file.
"""
def usage(m):
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

if 'usage' in task:
    p = Pool(3)
    p.map(usage,modes)

"""
Scatterplots of usages
"""
if 'plot' in task:
    for m in modes:
        """
        Load flows and find maximum values for normalisation
        """
        F = np.load('./results/'+m+'-flows.npy')
        Fmax = np.max(F,1)

        """
        Load nodes' usage of links for export and import respectively with
        dimensions export_usage = (links, nodes, lapse).
        """
        export_usage = np.load('./linkcolouring/old_'+m+'_copper_link_mix_export_all_alpha=same.npy')
        import_usage = np.load('./linkcolouring/old_'+m+'_copper_link_mix_import_all_alpha=same.npy')

        """
        Plots of export usage
        """
        nodes = export_usage.shape[1]
        links = export_usage.shape[0]
        for l in range(links):
            plt.figure()
            for t in range(lapse):
                linkflow = abs(F[l,t]*np.ones(nodes))
                usages = export_usage[l,:,t]/Fmax[l]
                """
                Find the three highest contributors and plot in red colors.
                """
                maxs = []
                colors = ['#b30000','#ff0000','#ff4d4d']
                for k in range(3):
                    i = usages.argmax()
                    maxs.append(usages[i])
                    usages = np.delete(usages,i)
                    plt.plot(linkflow[k],maxs[k],'.',color=str(colors[k]))
                """
                Plot other nodes' usages.
                """
                plt.plot(linkflow[3:],usages,'.k')
            plt.title('Export flows on link #'+str(l))
            plt.xlabel(r'$F_l(t)$ [MW]')
            plt.ylabel(r'$H_{ln}/\max(F_l)$')
            plt.savefig('./figures/exportflows-'+str(l)+'.png')

        """
        Plots of import usage
        """
        for l in range(links):
            plt.figure()
            for t in range(lapse):
                linkflow = abs(F[l,t]*np.ones(nodes))
                usages = import_usage[l,:,t]/Fmax[l]
                """
                Find the three highest contributors and plot in red colors.
                """
                maxs = []
                colors = ['#b30000','#ff0000','#ff4d4d']
                for k in range(3):
                    i = usages.argmax()
                    maxs.append(usages[i])
                    usages = np.delete(usages,i)
                    plt.plot(linkflow[k],maxs[k],'.',color=str(colors[k]))
                """
                Plot other nodes' usages.
                """
                plt.plot(linkflow,usages,'.k')
            plt.title('Import flows on link #'+str(l))
            plt.xlabel(r'$F_l(t)$ [MW]')
            plt.ylabel(r'$H_{ln}/\max(F_l)$')
            plt.savefig('./figures/importflows-'+str(l)+'.png')
