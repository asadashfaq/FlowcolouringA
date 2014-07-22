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
- usage:                calculate nodes' usages of each link
- plot:                 produce figures of usage
- conditional converge: check convergence of conditional average usages with bin sizes
- conditional usage:    calculate and produce figures of conditional average usages of each link
- conditional country:  figure of each node's weighted usage of all links

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
        plt.plot(diagflow,diag,'-k',lw=.9)
        qq = get_q(abs(F[l,:lapse]),.99)
        plt.plot([qq,qq],[0,max(diag)],'--k')

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
            max_usages[i] = 0
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
    lapse = 365*24  # number of hours to include
    print 'Lapse: '+str(lapse)+' hours'
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
    print 'Mode selected: Calculate usages.'
    lapse = 365*24  # number of hours to include
    print 'Lapse: '+str(lapse)+' hours.'
    print 'Running on 3 cores.'
    p = Pool(3)
    p.map(calc_usage,modes)


"""
Scatterplots of usages
"""
if 'plot' in task:
    print 'Mode selected: Plotting figures of usages'
    lapse = 365*24  # number of hours to include
    print 'Lapse: '+str(lapse)+' hours'
    top = 5
    for m in modes:
        print 'Plotting '+str(m)
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
    print 'Lapse: 70128'
    lapse = 70128

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
        ax = plt.subplot(111)
        plt.plot(F_vert,exp_vert,'.k',alpha=.2)

        # Plot 99% quantile of flow
        qq = get_q(abs(F[0,:lapse]),.99)
        plt.plot([qq,qq],[0,np.ceil(max(exp_vert)/500)*500],'--k',label=r'$99\%$')
        plt.axis([0,np.ceil(max(F_vert)/1000)*1000,0,np.ceil(max(exp_vert)/500)*500])
    
        # Get rid of large variables in memory
        export_usage = []
        exp_vert = []
        F = []
        F_vert = []

        # Calculate and plot conditional usages
        run = 0
        bin_sizes = np.linspace(100,1000,10)    # bin sizes to be tested
        colors = ['#ff0000','#ff5500','#ffaa00','#ffff00','#aaff00','#00ff00','#00ffff','#00aaff','#0055ff','#0000ff']
        names = ['100','200','300','400','500','600','700','800','900','1000']
        for bin_size in bin_sizes:
            bin_means,H = bin_maker(bin_size,F_matrix)
            plt.plot(bin_means,H[:,1],'-',color=str(colors[run]),label=names[run],lw=1.5)
            run += 1

        ax.set_title('Convergence check:\nAT\'s usage of the AT-CH link over '+str(lapse)+' hours')
        ax.set_xlabel(r'$F_l$ [MW]')
        ax.set_ylabel(r'$H_{ln}$ [MW]')

        # Shrink x-axis to make room for legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
        
        handles, labels = ax.get_legend_handles_labels()
        handles = np.append(handles[1:],handles[0])
        labels = np.append(labels[1:],labels[0])
        ax.legend(handles,labels,loc='center left', bbox_to_anchor=(1,0.5),title='Bin size')
         
        plt.savefig('./figures/convergence-'+str(lapse)+'.pdf')
        print 'Saved result to: ./figures/convergence-'+str(lapse)+'.pdf'

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
        
        for link in range(len(F)):
            plt.figure()
            ax = plt.subplot(111)
            Fmax = np.max(np.abs(F[link,:lapse]))
            Bin_means = []
            Total_H = []
            Bin_sums = np.array([])
            for node in range(len(N)):
                # Stacking and sorting data
                F_vert = np.reshape(F[link,:lapse],(len(F[link,:lapse]),1))
                exp_vert = np.reshape(export_usage[link,node,:lapse],(len(export_usage[link,node,:lapse]),1))
                F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
                F_matrix.sort(axis=0)
        
                # Calculate conditional usages
                M,bin_sum = bin_maker(bin_size,F_matrix,summed=1) # M = [[bin_means],[conditional usages]]
                Bin_means.append(M[0])
                Total_H.append(M[1])
                Bin_sums = np.append(Bin_sums,bin_sum)
                p0 = plt.plot(M[0],M[1],'-',color='#545454',alpha=.2,label='background',lw=1.5)
        
            # Plot 99% quantile of link flow
            qq = get_q(abs(F[link,:lapse]),.99)
            pq = plt.plot([qq,qq],[0,Fmax],'--k',label='99%')
            
            # Plot top 5 contributors' conditional usages for each link.
            max_users = []
            names = []
            for k in range(top):
                i = Bin_sums.argmax()
                max_users.append(i)
                Bin_sums[i] = 0 # = np.delete(Bin_sums,i)
                names.append(str(N[i].label))
    
            colors = ['#ff0000','#ff5500','#ffaa00','#ffff00','#aaff00','#00ff00','#00ffff','#00aaff','#0055ff','#0000ff']
            p1 = plt.plot(Bin_means[max_users[0]],Total_H[max_users[0]],'-',color=str(colors[0]),label=str(names[0]),lw=1.5)
            p2 = plt.plot(Bin_means[max_users[1]],Total_H[max_users[1]],'-',color=str(colors[2]),label=str(names[1]),lw=2)
            p3 = plt.plot(Bin_means[max_users[2]],Total_H[max_users[2]],'-',color=str(colors[3]),label=str(names[2]),lw=2)
            p4 = plt.plot(Bin_means[max_users[3]],Total_H[max_users[3]],'-',color=str(colors[4]),label=str(names[3]),lw=2)
            p5 = plt.plot(Bin_means[max_users[4]],Total_H[max_users[4]],'-',color=str(colors[8]),label=str(names[4]),lw=2)
            
            label = link_label(link,N)
            ax.set_title('Top 5 users of the '+str(label)+' link with bin size '+str(bin_size)+' MW')
            ax.set_xlabel(r'$F_l$ [MW]')
            ax.set_ylabel(r'$\left\langle H_{ln}|F_l\right\rangle$ [MW]')
            plt.axis([0,np.ceil(Fmax/1000)*1000,0,np.ceil(np.max(np.max(Total_H))/500)*500])
    
            # Shrink x-axis to make room for legend
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
            
            handles, labels = ax.get_legend_handles_labels()
            handles = np.append(handles[-5:],handles[-6])
            labels = np.append(labels[-5:],labels[-6])
            ax.legend(handles,labels,loc='center left', bbox_to_anchor=(1,0.5),title='Country')
            
            plt.savefig('./figures/conditional-'+str(lapse)+'-'+str(link)+'.png')
            plt.close('all')
        print 'Saved results to: ./figures/conditional-lapse-link.png'

    if 'country' in task:
        """
        Figures of how important each country is for the flow on each link
        """
        print "Calculating weighted sum of each country's conditional usage."

        # Load data
        N = EU_Nodes_usage('linear.npz')
        F = abs(np.load('./results/linear-flows.npy'))
        export_usage = np.load('./linkcolouring/old_linear_copper_link_mix_export_all_alpha=same.npy')
        bin_size = 400  # determined from the _convergence_ mode.
        N_usages = np.zeros((len(N),len(F))) # empty array for weighted sums of conditional usages
        links = range(len(F))

        for node in range(len(N)):
            for link in range(len(F)):
                # Stacking and sorting data
                F_vert = np.reshape(F[link,:lapse],(len(F[link,:lapse]),1))
                exp_vert = np.reshape(export_usage[link,node,:lapse],(len(export_usage[link,node,:lapse]),1))
                F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
                F_matrix.sort(axis=0)
                
                # Calculate weigthed sums of usages and save for later use
                M,bin_sum = bin_maker(bin_size,F_matrix,summed=1)
                N_usages[node,link] = bin_sum

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
        for node in range(len(N)):
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
        
        node_id = str(N[node].label)
        ax.set_title('Comparing network usages for weighted sums of link usages')
        ax.set_xlabel(r'Node')
        ax.set_ylabel(r'Network usage normalised to avg. network usage')
        plt.savefig('./figures/network-usage.png')
        plt.close('all')
        print 'Saved results to: ./figures/network-usage.png'
