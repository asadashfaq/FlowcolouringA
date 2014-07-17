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
Script to create scatter plots of nodes' usages of links for three export schemes:
- Linear (most localised)
- Square (cooperative)
- Market (RND)

Call the program with only one of the following command line arguments:
- solve:                solve the network and save solutions
- usage:                calculate nodes' usages of links
- plot:                 produce figures of usage
- conditional:          calculate conditional average usages
- conditional converge: check convergence of conditional average usages with bin sizes

Example call:
python usage.py solve
"""


"""
Initialisation
"""
if len(sys.argv)<2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

modes = ['linear','square','RND']   # Export schemes. 'RND' is also known as 'Market'.
lapse = 365*24  # number of hours to include
print 'Lapse: '+str(lapse)+' hours'


"""
Helper functions
"""
def link_label(num,N):
    """
    Translate link number into a string like: 'Country1-Country2'.
    """
    label = get_link_direction(num,N)
    return str(label[0].label)+'-'+str(label[1].label)

def calc_usage(m):
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

def scatter_plotter(N,F,Fmax,usage,direction,mode,lapse):
    """
    Scatter plots of nodes' import/export usages of links saved to ./figures/.
    """
    colors = ['#670000','#ff0000','#e67400']
    nodes = usage.shape[1]
    links = usage.shape[0]
    for l in range(links):
        diag = []
        diagflow = []
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
            """
            Plot diagonal and 99 quantile of link flow
            """
            diag.append(sum(usages)+sum(maxs))
            diagflow.append(linkflow[0])
        plt.plot(diagflow,diag,'-k')
        qq = get_q(abs(F[l,:lapse]),.99)
        plt.plot([qq,qq],[0,max(diag)],'-k')

        label = link_label(l,N)
        ax.set_title(str(mode)+' '+str(direction)+' flows on link '+label)
        ax.set_xlabel(r'$F_l(t)$ [MW]')
        ax.set_ylabel(r'$H_{ln}/\max(F_l)$')
        plt.axis([0,max(abs(F[l,:lapse])),0,max(diag)])
    
        # Shrink x-axis to make room for legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
    
        names = ['1st', '2nd', '3rd', 'Rest']
        ax.legend((names),loc='center left', bbox_to_anchor=(1,0.5),title='Contributions')
    
        plt.savefig('./figures/'+str(mode)+'-'+str(direction)+'-flows-'+str(l+1)+'.png')
        plt.close('all') # fixes memory leak.
    return

def top_plotter(top,N,F,Fmax,usage,direction,mode,lapse):
    """
    Plots of top 5 contributors to each links transmission capacity
    """
    colors = ['#e32636','#008000','#007fff','#fdee00','#ed9121']
    nodes = usage.shape[1]
    links = usage.shape[0]
    for l in range(links):
        diag = []
        diagflow = []
        plt.figure()
        ax = plt.subplot(111)
        max_usages = np.max(usage[l,:,:lapse],1)
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
            ax.plot(linkflow,usages,'.',color='#b2beb5',alpha=.3,label='rest')
            for k in range(top):
                ax.plot(linkflow[k],usages[max_users[k]],'.',color=str(colors[k]),label=names[k])
            """
            Plot diagonal and 99 quantile of link flow
            """
            diag.append(sum(usages))
            diagflow.append(linkflow[0])
        plt.plot(diagflow,diag,'-k')
        qq = get_q(abs(F[l,:lapse]),.99)
        plt.plot([qq,qq],[0,max(diag)*1.05],'-k')

        label = link_label(l,N)
        ax.set_title('Top '+str(top)+' '+str(mode)+' '+str(direction)+' flows on link '+label)
        ax.set_xlabel(r'$F_l(t)$ [MW]')
        ax.set_ylabel(r'$H_{ln}/\max(F_l)$')
        plt.axis([0,max(abs(F[l,:lapse])),0,max(diag)])

        # Shrink x-axis to make room for legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.85, box.height])

        ax.legend((names),loc='center left', bbox_to_anchor=(1,0.5),title='Country')
        plt.savefig('./figures/top-'+str(mode)+'-'+str(direction)+'-flows-'+str(l+1)+'.png')
        plt.close('all') # fixes memory leak.
    return

def bin_maker(bin_size,F_matrix,summed=None):
    """
    Calculate the conditional usage as a function of the flow on the link according to bin_size
    """
    bin_max = np.ceil(max(F_matrix[:,0])/bin_size)*bin_size # round up to nearest bin_size
    nbins = bin_max/bin_size # number of bins
    bin_means = np.linspace(.5*bin_size,bin_max-(.5*bin_size),nbins) # mean values of bins

    H_temp = []
    H = np.zeros((nbins,4)) # [#nodes, mean usage, min usage, max usage]
    bin_id = 0
    for t in range(lapse):
        if F_matrix[t,0] < (bin_id+1)*bin_size:
            H_temp.append(F_matrix[t,1])
        else:
            if len(H_temp)>0:
                H[bin_id,0] = len(H_temp)
                H[bin_id,1] = sum(H_temp)/len(H_temp)
                H[bin_id,2] = min(H_temp)
                H[bin_id,3] = max(H_temp)
                bin_id += 1
                H_temp = []
            else: # no data in the bin
                H[bin_id,0] = 0
                H[bin_id,1] = 0
                H[bin_id,2] = 0
                H[bin_id,3] = 0
                bin_id += 1
                H_temp = []
    if summed:
        bin_sum = sum(np.multiply(H[:,1],bin_means))
        return np.array([bin_means,H[:,1]]),bin_sum
    else:
        return bin_means,H


"""
Solving flows for different export schemes
"""
if 'solve' in task:
    print 'Mode selected: Solving network flows'
    alphas = np.ones(30)*.7
    gammas = np.ones(30)
    
    Nodes = EU_Nodes_usage()
    
    for m in modes:
        print 'Solving: '+str(modes[m])
        N,F = au.solve(Nodes,mode=m+' copper verbose',lapse=lapse)
        N.save_nodes(m)
        np.save('./results/'+m+'-flows',F)


"""
Calculate powermixes and nodes' usages of links and save results to file.
"""
if (('usage' in task) and ('conditional' not in task)):
    print 'Mode selected: Calculate usages'
    p = Pool(3)
    p.map(calc_usage,modes)


"""
Scatterplots of usages
"""
if 'plot' in task:
    print 'Mode selected: Plotting figures of usages'
    top = 5
    for m in modes:
        print 'Plotting'+str(modes[m])
        """
        Load flows and find maximum values for normalisation
        """
        N = EU_Nodes_usage(m+'.npz')
        F = np.load('./results/'+m+'-flows.npy')
        Fmax = np.max(np.abs(F[:,:lapse]),1)

        """
        Load usages, do scatter plots and plot of top 5 contributors to a links capacity.
        """
        export_usage = np.load('./linkcolouring/old_'+m+'_copper_link_mix_export_all_alpha=same.npy')
        scatter_plotter(N,F,Fmax,export_usage,'export',m,lapse)
        top_plotter(top,N,F,Fmax,export_usage,'export',m,lapse)
        export_usage = [] # frees roughly 800 MB RAM before loading import usage.

        import_usage = np.load('./linkcolouring/old_'+m+'_copper_link_mix_import_all_alpha=same.npy')
        scatter_plotter(N,F,Fmax,import_usage,'import',m,lapse)
        top_plotter(top,N,F,Fmax,import_usage,'import',m,lapse)
        import_usage = []

"""
Looking at conditional average usages in absolute units
"""
if 'conditional' in task:
    print 'Mode selected: Conditional usage'

    if 'converge' in task:
        """
        Apply a number of bin sizes to investigate convergence
        """
        print 'Checking convergence'

        # Load data
        F = abs(np.load('./results/linear-flows.npy')) # [link, flow], len = lapse
        export_usage = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy') # Dimension: link x node x lapse
    
        # Stacking and sorting data
        F_vert = np.reshape(F[0,:lapse],(len(F[0,:lapse]),1))
        exp_vert = np.reshape(export_usage[0,0,:lapse],(len(export_usage[0,0,:lapse]),1))
        F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
        F_matrix.sort(axis=0)
    
        # Plot data as scatter background
        plt.figure()
        plt.plot(F_vert,exp_vert,'.k',alpha=.3)

        # Plot 99% quantile of flow
        qq = get_q(abs(F[l,:lapse]),.99)
        plt.plot([qq,qq],[0,max(exp_vert)],'-k')
    
        # Get rid of large variables in memory
        export_usage = []
        exp_vert = []
        F = []
        F_vert = []

        # Calculate and plot conditional usages
        run = 0
        bin_sizes = np.linspace(100,1000,10)    # bin sizes to be tested
        colors = ['#ff0000','#ff5500','#ffaa00','#ffff00','#aaff00','#00ff00','#00ffff','#00aaff','#0055ff','#0000ff']
        for bin_size in bin_sizes:
            bin_means,H = bin_maker(bin_size,F_matrix)
            plt.plot(bin_means,H[:,1],'-',color=str(colors[run]))
            run += 1

        """
        Pretty plotting
        """
        plt.title('Convergence check: AT\'s usage of the AT-CH link')
        plt.xlabel(r'$F_l$ [MW]')
        plt.ylabel(r'$H_{ln}$ [MW]')
        plt.savefig('./figures/convergence-'+str(lapse)+'.png')
        print 'Saved results to: ./figures/convergence-'+str(lapse)+'.png'

    if 'usage' in task:
        """
        Calculate conditional usages for all nodes and plot the top 5 contributors for each node
        """
        print 'Calculating conditional usages'

        # Load data
        N = EU_Nodes_usage('linear.npz')
        F = abs(np.load('./results/linear-flows.npy'))
        export_usage = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy')

        bin_size = 400  # determined from the _convergence_ mode.
        top = 5         # how many contributors to plot
        
        # insert loop over links here
        link = 0

        # Plot 99% quantile of flow
        plt.figure()
#        qq = get_q(abs(F[link,:lapse]),.99)
#        plt.plot([qq,qq],[0,max(exp_vert)],'-k')

        Bin_means = []
        Total_H = []
        Bin_sums = np.array([])
        for node in range(len(N)):
            # Stacking and sorting data
            F_vert = np.reshape(F[node,:lapse],(len(F[node,:lapse]),1))
            exp_vert = np.reshape(export_usage[link,node,:lapse],(len(export_usage[link,node,:lapse]),1))
            F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
            F_matrix.sort(axis=0)
    
            # Calculate conditional usages
            M,bin_sum = bin_maker(bin_size,F_matrix,summed=1) # M = [[bin_means],[conditional usages]]
            Bin_means.append(M[0])
            Total_H.append(M[1])
            Bin_sums = np.append(Bin_sums,bin_sum)
            plt.plot(M[0],M[1],'-',color='#545454',alpha=.2)
    
        # Plot top 5 contributors' conditional usages for each link.
        max_users = []
        names = []
        for k in range(top):
            i = Bin_sums.argmax()
            max_users.append(i)
            Bin_sums = np.delete(Bin_sums,i)
            names.append(str(N[i].label))

        run = 0
        colors = ['#ff0000','#ff5500','#ffaa00','#ffff00','#aaff00','#00ff00','#00ffff','#00aaff','#0055ff','#0000ff']
        for c in range(top):
            country = int(c)
            plt.plot(Bin_means[max_users[country]],Total_H[max_users[country]],'-',color=str(colors[run]))
            run += 1
        """
        Pretty plotting
        """
        plt.title('title')
        plt.xlabel(r'$F_l$ [MW]')
        plt.ylabel(r'$H_{ln}$ [MW]')
#        plt.savefig('./figures/convergence-'+str(lapse)+'.png')
#        print 'Saved results to: ./figures/convergence-'+str(lapse)+'.png'
        plt.show()
