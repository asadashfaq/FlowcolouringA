import sys
import numpy as np
from pylab import plt
import multiprocessing as mp
from aurespf.tools import *
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
import matplotlib as mpl
from EUgrid import EU_Nodes_usage

"""
Script that makes network figures of a country's usage for import, export and
the combination
"""

if len(sys.argv)>1:
    task = str(sys.argv[1:])

colwidth = (3.425)
dcolwidth = (2*3.425+0.236)

def drawnet_usage(N=None,mode='linear',direction='combined'):
    """
    Make network figures of a node's usage of links for both import, export and
    combined. Adapted from drawnet() in aurespf.plotting
    """
    
    if N==None:
        N=EU_Nodes_usage()
    G=nx.Graph()
    nodelist=[]

    # Add nodes and labels to networkx object for plotting
    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    # LF is a list of links
    K,h,LF=AtoKh_old(N)

    for l in LF:
        G.add_edge(l[0], l[1], id=l[2])
    
    # Define position of nodes
    pos={}
    pos['AT']=[0.55,0.45]
    pos['FI']=[.95,1.1]
    pos['NL']=[0.40,0.85]
    pos['BA']=[0.65,0.15]
    pos['FR']=[0.15,0.60]
    pos['NO']=[0.5,1.1]
    pos['BE']=[0.275,0.775]
    pos['GB']=[0.10,1.05]
    pos['PL']=[0.75,0.8]
    pos['BG']=[0.9,0.0]
    pos['GR']=[0.7,0.0]
    pos['PT']=[0.0,0.15]
    pos['CH']=[0.4,0.45]
    pos['HR']=[0.75,0.3]
    pos['RO']=[1.0,0.15]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RS']=[0.85,0.15]
    pos['DE']=[0.45,0.7]
    pos['IE']=[0.0,0.95]
    pos['SE']=[0.75,1.0]
    pos['DK']=[0.5,0.875]
    pos['IT']=[0.4,0.2]
    pos['SI']=[0.55,0.3]
    pos['ES']=[0.15,0.35]
    pos['LU']=[0.325,0.575]
    pos['SK']=[0.90,0.55]
    pos['EE']=[1.0,0.985]
    pos['LV']=[0.975,0.87]
    pos['LT']=[0.925,0.77]
   
    # Define color scale for links
    blueDict = {'red': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(1.0, 1.0, 1.0))}
    
    blueDict2 = {'red': ((0.0, 1.0, 1.0),(.3,0,0),(.5, 0.5, 0.5),(1.0, 1, 1)),
             'green': ((0.0, 1.0, 1.0),(.3,0,0),(.5, 0, 0),(1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(.3,1,1),(.5, 0.5, 0.5),(1.0, 0, 0))}

    cmap = LinearSegmentedColormap('blue',blueDict2,1000)
    
    # Load usages for given mode and direction
    N_usages = np.load('./results/Node_contrib_'+mode+'_'+direction+'_70128.npy')
    quantiles = np.load('./results/quantiles_'+str(mode)+'_70128.npy')

    # Pick a particular node
    for n in N:
        # Calculate colors of links
        N_usages[n.id] = N_usages[n.id]/quantiles
        col = [ (cmap(l)) for l in N_usages[n.id]]
        
        # Create a new figure and plot network below
        fig = plt.figure(dpi=400,figsize=(0.85*dcolwidth,0.85*0.8*dcolwidth))
    
        # color bar in bottom of figure
        ax1 = fig.add_axes([0.05,0.08,0.9,.08])
        
        cbl = mpl.colorbar.ColorbarBase(ax1,cmap,orientation='vhorizontal')
        ax1.set_xlabel(mode+' '+direction+r" usage $C_n/C^{\,99\%}$")
        ax1.xaxis.set_label_position('top') 
        
        ax2 = fig.add_axes([-0.05,0.15,1.1,0.95])
        
        # Set color of nodes, highlight one and draw all 
        node_c = [ "#000000" for node in N]
        node_c[n.id] = "#B30000"
        nx.draw_networkx_nodes(G,pos,node_size=500,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))
    
        # Draw links colored by usage of node n
        edges=[(u,v) for (u,v,d) in G.edges(data=True)]
        edge_id = [d['id'] for (u,v,d) in G.edges(data=True)]
        
        color_sort = []
        for i in range(len(col)):
            color_sort.append(col[edge_id[i]])
        nx.draw_networkx_edges(G,pos,edgelist=edges,width=3.5,edge_color=color_sort,alpha=1)
        
        # Draw country names
        nx.draw_networkx_labels(G,pos,font_size=12,font_color='w',font_family='sans-serif')
        ax2.axis('off') 
    
        # Save figure
        plt.savefig("./figures/network_figures/"+mode+"/"+str(n.id)+'_'+str(direction)+".png")


if 'plot' in task:
    print 'Plotting network figures'
    modes = ['linear' ,'square', 'RND']
    directions = ['import', 'export', 'combined']
    N = EU_Nodes_usage()
    for mode in modes:
        print 'Mode: ',mode
        for direction in directions:
            print 'Direction: ',direction
            drawnet_usage(N,mode,direction)
