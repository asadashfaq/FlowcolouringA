import os
import sys
import numpy as np
from multiprocessing import Pool
import networkx as nx
from pylab import plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from aurespf.tools import *
from EUgrid import EU_Nodes_usage, EU_Nodes_regions, EU_Nodes_superRegions
from link_namer import *
from functions import *

"""
Vector flow tracing on top of the up/down stream flow tracing algorithm.

This script can be called from command line with one of the following inputs:
- trace:                vectorize the flow tracing and save results
- plot network:         load results from above and plot network figures for each color
- plot network day:     same as above but split for daytime and nighttime
- plot network total:   plot total network usage for each color
- plot usage:           plot average usage for each color
- sanity:               check whether the individual colors add to the scalar flow tracing

Results from this script go to the folder: ./results/vector/
Figures go to the folder: ./figures/vector.
"""

"""
Initialisation
"""
if len(sys.argv) < 2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

modes = ['linear', 'square']
directions = ['import', 'export', 'combined']
outPath = './results/vector/'
figPath = './figures/vector/'

# Node indices and names sorted after descending mean load
loadOrder = [18, 4, 7, 22, 24, 20, 8, 5, 2, 6, 1, 15, 0, 10, 14,
             9, 11, 12, 16, 21, 17, 19, 3, 26, 13, 29, 27, 23, 28, 25]

loadNames = np.array(['DE', 'FR', 'GB', 'IT', 'ES', 'SE', 'PL', 'NO', 'NL',
                      'BE', 'FI', 'CZ', 'AT', 'GR', 'RO', 'BG', 'PT', 'CH',
                      'HU', 'DK', 'RS', 'IE', 'BA', 'SK', 'HR', 'LT', 'EE',
                      'SI', 'LV', 'LU'], dtype='|S4')


def usageCalc(F, quantiles, Usages, nodes, links, name):
    """
    Calculate usages and save to file
    """
    Node_contributions = np.zeros((nodes, links))  # empty array for calculated usages
    for node in range(nodes):
        for link in range(links):
            # Stacking and sorting data
            F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
            exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
            F_matrix = np.hstack([F_vert, exp_vert])
            F_matrix[F_matrix[:, 0].argsort()]

            H, bin_edges = binMaker(F_matrix, quantiles[link], lapse)
            Node_contributions[node, link] = node_contrib(H, bin_edges, linkID=link)

    np.save(outPath + 'Node_contrib_' + mode + '_' + direction + '_' + str(name) + '.npy', Node_contributions)
    return


def usageCalcDaily(iF, quantiles, iUsages, nodes, links, name):
    """
    Calculate usages for daytime and nighttime and save to file
    """
    days = lapse / 24
    hours = np.append(range(6, days * 24), range(6))
    splitHours = np.split(hours, days * 2)
    dayTime = np.concatenate(splitHours[0::2])
    nightTime = np.concatenate(splitHours[1::2])

    for d in ['day', 'night']:
        if d == 'day':
            Usages = iUsages[:, :, dayTime]
            F = iF[:, dayTime]
        elif d == 'night':
            Usages = iUsages[:, :, nightTime]
            F = iF[:, nightTime]

        Node_contributions = np.zeros((nodes, links))  # empty array for calculated usages
        for node in range(nodes):
            for link in range(links):
                # Stacking and sorting data
                F_vert = np.reshape(F[link, :lapse / 2], (len(F[link, :lapse / 2]), 1))
                exp_vert = np.reshape(Usages[link, node, :lapse / 2], (len(Usages[link, node, :lapse / 2]), 1))
                F_matrix = np.hstack([F_vert, exp_vert])
                F_matrix[F_matrix[:, 0].argsort()]

                H, bin_edges = binMaker(F_matrix, quantiles[link], lapse / 2)
                Node_contributions[node, link] = node_contrib(H, bin_edges, linkID=link)

        np.save(outPath + 'Node_contrib_' + mode + '_' + direction + '_' + d + '_' + str(name) + '.npy', Node_contributions)
    return


def drawnet_usage(N=None, scheme='linear', direction='combined', color='solar'):
    """
    Make network figures of a node's usage of links for both import, export and
    combined. A figure for each color is created.
    Adapted from drawnet() in aurespf.plotting
    """
    colwidth = (3.425)
    dcolwidth = (2 * 3.425 + 0.236)

    if not N:
        N = EU_Nodes_usage()
    G = nx.Graph()
    nodelist = []

    # Add nodes and labels to networkx object for plotting
    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    # LF is a list of links
    K, h, LF = AtoKh_old(N)

    for l in LF:
        G.add_edge(l[0], l[1], id=l[2])

    # Define position of nodes
    pos = {}
    pos['AT'] = [0.55, 0.45]
    pos['FI'] = [.95, 1.1]
    pos['NL'] = [0.40, 0.85]
    pos['BA'] = [0.65, 0.15]
    pos['FR'] = [0.15, 0.60]
    pos['NO'] = [0.5, 1.1]
    pos['BE'] = [0.275, 0.775]
    pos['GB'] = [0.10, 1.05]
    pos['PL'] = [0.75, 0.8]
    pos['BG'] = [0.9, 0.0]
    pos['GR'] = [0.7, 0.0]
    pos['PT'] = [0.0, 0.15]
    pos['CH'] = [0.4, 0.45]
    pos['HR'] = [0.75, 0.3]
    pos['RO'] = [1.0, 0.15]
    pos['CZ'] = [0.75, 0.60]
    pos['HU'] = [1.0, 0.45]
    pos['RS'] = [0.85, 0.15]
    pos['DE'] = [0.45, 0.7]
    pos['IE'] = [0.0, 0.95]
    pos['SE'] = [0.75, 1.0]
    pos['DK'] = [0.5, 0.875]
    pos['IT'] = [0.4, 0.2]
    pos['SI'] = [0.55, 0.3]
    pos['ES'] = [0.15, 0.35]
    pos['LU'] = [0.325, 0.575]
    pos['SK'] = [0.90, 0.55]
    pos['EE'] = [1.0, 0.985]
    pos['LV'] = [0.975, 0.87]
    pos['LT'] = [0.925, 0.77]

    # Define color scale for links
    blueDict = {'red': ((0.0, 1.0, 1.0), (0.15, 0.0, 0.0), (1.0, 0.0, 0.0)),
                'green': ((0.0, 1.0, 1.0), (0.15, 0.0, 0.0), (1.0, 0.0, 0.0)),
                'blue': ((0.0, 1.0, 1.0), (0.15, 1.0, 1.0), (1.0, 0.3, 0.3))}

    orangeDict = {'red': ((0.0, 1.0, 1.0), (0.1, 1.0, 1.0), (1.0, 0.3, 0.3)),
                  'green': ((0.0, 1.0, 1.0), (0.1, 0.65, .65), (1.0, .1, .1)),
                  'blue': ((0.0, 1.0, 1.0), (0.1, 0.0, 0.0), (1.0, 0.0, 0.0))}

    brownDict = {'red': ((0.0, 1.0, 1.0), (0.05, 0.57, 0.57), (1.0, 0.2, 0.2)),
                 'green': ((0.0, 1.0, 1.0), (0.05, 0.36, .36), (1.0, .12, .12)),
                 'blue': ((0.0, 1.0, 1.0), (0.05, 0.15, 0.15), (1.0, 0.01, 0.01))}

    if color == 'wind':
        cmap = LinearSegmentedColormap('blue', blueDict, 1000)
    elif color == 'solar':
        cmap = LinearSegmentedColormap('orange', orangeDict, 1000)
    else:
        cmap = LinearSegmentedColormap('brown', brownDict, 1000)

    # Load usages for given scheme and direction
    N_usages = np.load(outPath + '/Node_contrib_' + scheme + '_' + direction + '_' + color + '.npy')
    quantiles = np.load('./results/quantiles_' + str(scheme) + '_70128.npy')

    # Pick a particular node
    for n in N:
        # Calculate colors of links
        N_usages[n.id] = N_usages[n.id] / quantiles
        col = [(cmap(l)) for l in N_usages[n.id]]

        # Create a new figure and plot network below
        fig = plt.figure(dpi=400, figsize=(0.85 * dcolwidth, 0.85 * 0.8 * dcolwidth))

        # color bar in bottom of figure
        ax1 = fig.add_axes([0.05, 0.08, 0.9, .08])
        cbl = mpl.colorbar.ColorbarBase(ax1, cmap, orientation='horizontal')

        # Label just above color bar
        if scheme == 'linear':
            xlabel = 'Most localised'
        elif scheme == 'square':
            xlabel = 'Synchronised'
        else:
            xlabel = 'Market'
        # ax1.set_xlabel(xlabel+' '+direction+r" usage $C_n/C^{\,99\%}$")
        ax1.set_xlabel(r'$\mathcal{K}_{ln}/\mathcal{K}^T_l$')
        ax1.xaxis.set_label_position('top')

        ax2 = fig.add_axes([-0.05, 0.15, 1.1, 0.95])

        # Set color of nodes, highlight one and draw all
        node_c = ["#000000" for node in N]
        node_c[n.id] = "#B30000"
        nx.draw_networkx_nodes(G, pos, node_size=500, nodelist=nodelist, node_color=node_c, facecolor=(1, 1, 1))

        # Draw links colored by usage of node n
        edges = [(u, v) for (u, v, d) in G.edges(data=True)]
        edge_id = [d['id'] for (u, v, d) in G.edges(data=True)]

        color_sort = []
        for i in range(len(col)):
            color_sort.append(col[edge_id[i]])
        nx.draw_networkx_edges(G, pos, edgelist=edges, width=3.5, edge_color=color_sort, alpha=1)

        # Draw country names
        nx.draw_networkx_labels(G, pos, font_size=12, font_color='w', font_family='sans-serif')
        ax2.axis('off')

        # Save figure
        plt.savefig(figPath + "network/" + scheme + "/" + str(n.id) + '_' + color + '_' + direction + ".png")
        plt.close()


def drawnet_total(N=None, scheme='linear', direction='combined', color='solar'):
    """
    Make network figures for each color of the total network usage.
    """
    colwidth = (3.425)
    dcolwidth = (2 * 3.425 + 0.236)

    if not N:
        N = EU_Nodes_usage()
    G = nx.Graph()
    nodelist = []

    # Add nodes and labels to networkx object for plotting
    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    # LF is a list of links
    K, h, LF = AtoKh_old(N)

    for l in LF:
        G.add_edge(l[0], l[1], id=l[2])

    # Define position of nodes
    pos = {}
    pos['AT'] = [0.55, 0.45]
    pos['FI'] = [.95, 1.1]
    pos['NL'] = [0.40, 0.85]
    pos['BA'] = [0.65, 0.15]
    pos['FR'] = [0.15, 0.60]
    pos['NO'] = [0.5, 1.1]
    pos['BE'] = [0.275, 0.775]
    pos['GB'] = [0.10, 1.05]
    pos['PL'] = [0.75, 0.8]
    pos['BG'] = [0.9, 0.0]
    pos['GR'] = [0.7, 0.0]
    pos['PT'] = [0.0, 0.15]
    pos['CH'] = [0.4, 0.45]
    pos['HR'] = [0.75, 0.3]
    pos['RO'] = [1.0, 0.15]
    pos['CZ'] = [0.75, 0.60]
    pos['HU'] = [1.0, 0.45]
    pos['RS'] = [0.85, 0.15]
    pos['DE'] = [0.45, 0.7]
    pos['IE'] = [0.0, 0.95]
    pos['SE'] = [0.75, 1.0]
    pos['DK'] = [0.5, 0.875]
    pos['IT'] = [0.4, 0.2]
    pos['SI'] = [0.55, 0.3]
    pos['ES'] = [0.15, 0.35]
    pos['LU'] = [0.325, 0.575]
    pos['SK'] = [0.90, 0.55]
    pos['EE'] = [1.0, 0.985]
    pos['LV'] = [0.975, 0.87]
    pos['LT'] = [0.925, 0.77]

    # Define color scale for links
    blueDict = {'red': ((0.0, 1.0, 1.0), (0.3, 0.0, 0.0), (1.0, 0.0, 0.0)),
                'green': ((0.0, 1.0, 1.0), (0.3, 0.0, 0.0), (1.0, 0.0, 0.0)),
                'blue': ((0.0, 1.0, 1.0), (0.3, 1.0, 1.0), (1.0, 0.3, 0.3))}

    orangeDict = {'red': ((0.0, 1.0, 1.0), (0.3, 1.0, 1.0), (1.0, 0.3, 0.3)),
                  'green': ((0.0, 1.0, 1.0), (0.3, 0.65, .65), (1.0, .1, .1)),
                  'blue': ((0.0, 1.0, 1.0), (0.3, 0.0, 0.0), (1.0, 0.0, 0.0))}

    brownDict = {'red': ((0.0, 1.0, 1.0), (0.3, 0.57, 0.57), (1.0, 0.2, 0.2)),
                 'green': ((0.0, 1.0, 1.0), (0.3, 0.36, .36), (1.0, .12, .12)),
                 'blue': ((0.0, 1.0, 1.0), (0.3, 0.15, 0.15), (1.0, 0.01, 0.01))}

    if color == 'wind':
        cmap = LinearSegmentedColormap('blue', blueDict, 1000)
    elif color == 'solar':
        cmap = LinearSegmentedColormap('orange', orangeDict, 1000)
    else:
        cmap = LinearSegmentedColormap('brown', brownDict, 1000)

    # Load usages for given scheme and direction
    N_usages = np.load(outPath + '/Node_contrib_' + scheme + '_' + direction + '_' + color + '.npy')
    quantiles = np.load('./results/quantiles_' + str(scheme) + '_70128.npy')

    # Calculate colors of links
    linkUsages = np.sum(N_usages, 0) / quantiles
    col = [(cmap(l)) for l in linkUsages]

    # Create a new figure and plot network below
    fig = plt.figure(dpi=400, figsize=(0.85 * dcolwidth, 0.85 * 0.8 * dcolwidth))

    # color bar in bottom of figure
    ax1 = fig.add_axes([0.05, 0.08, 0.9, .08])
    cbl = mpl.colorbar.ColorbarBase(ax1, cmap, orientation='horizontal')

    # Label just above color bar
    if scheme == 'linear':
        xlabel = 'Most localised'
    elif scheme == 'square':
        xlabel = 'Synchronised'
    else:
        xlabel = 'Market'
    # ax1.set_xlabel(xlabel+' '+direction+r" usage $C_n/C^{\,99\%}$")
    ax1.set_xlabel(r'$\mathcal{K}_{l}/\mathcal{K}^T_l$')
    ax1.xaxis.set_label_position('top')

    ax2 = fig.add_axes([-0.05, 0.15, 1.1, 0.95])

    # Set color of nodes and draw all
    node_c = ["#000000" for node in N]
    nx.draw_networkx_nodes(G, pos, node_size=500, nodelist=nodelist, node_color=node_c, facecolor=(1, 1, 1))

    # Draw links colored by usage of node n
    edges = [(u, v) for (u, v, d) in G.edges(data=True)]
    edge_id = [d['id'] for (u, v, d) in G.edges(data=True)]

    color_sort = []
    for i in range(len(col)):
        color_sort.append(col[edge_id[i]])
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=3.5, edge_color=color_sort, alpha=1)

    # Draw country names
    nx.draw_networkx_labels(G, pos, font_size=12, font_color='w', font_family='sans-serif')
    ax2.axis('off')

    # Save figure
    plt.savefig(figPath + "network/" + scheme + "/" + 'total_' + color + '_' + direction + ".png")
    plt.close()


def drawnet_day(N=None, scheme='linear', direction='combined', color='solar'):
    """
    Make network figures for each color of the total network usage during the day and night.
    """
    colwidth = (3.425)
    dcolwidth = (2 * 3.425 + 0.236)

    if not N:
        N = EU_Nodes_usage()
    G = nx.Graph()
    nodelist = []

    # Add nodes and labels to networkx object for plotting
    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    # LF is a list of links
    K, h, LF = AtoKh_old(N)

    for l in LF:
        G.add_edge(l[0], l[1], id=l[2])

    # Define position of nodes
    pos = {}
    pos['AT'] = [0.55, 0.45]
    pos['FI'] = [.95, 1.1]
    pos['NL'] = [0.40, 0.85]
    pos['BA'] = [0.65, 0.15]
    pos['FR'] = [0.15, 0.60]
    pos['NO'] = [0.5, 1.1]
    pos['BE'] = [0.275, 0.775]
    pos['GB'] = [0.10, 1.05]
    pos['PL'] = [0.75, 0.8]
    pos['BG'] = [0.9, 0.0]
    pos['GR'] = [0.7, 0.0]
    pos['PT'] = [0.0, 0.15]
    pos['CH'] = [0.4, 0.45]
    pos['HR'] = [0.75, 0.3]
    pos['RO'] = [1.0, 0.15]
    pos['CZ'] = [0.75, 0.60]
    pos['HU'] = [1.0, 0.45]
    pos['RS'] = [0.85, 0.15]
    pos['DE'] = [0.45, 0.7]
    pos['IE'] = [0.0, 0.95]
    pos['SE'] = [0.75, 1.0]
    pos['DK'] = [0.5, 0.875]
    pos['IT'] = [0.4, 0.2]
    pos['SI'] = [0.55, 0.3]
    pos['ES'] = [0.15, 0.35]
    pos['LU'] = [0.325, 0.575]
    pos['SK'] = [0.90, 0.55]
    pos['EE'] = [1.0, 0.985]
    pos['LV'] = [0.975, 0.87]
    pos['LT'] = [0.925, 0.77]

    # Define color scale for links
    blueDict = {'red': ((0.0, 1.0, 1.0), (0.3, 0.0, 0.0), (1.0, 0.0, 0.0)),
                'green': ((0.0, 1.0, 1.0), (0.3, 0.0, 0.0), (1.0, 0.0, 0.0)),
                'blue': ((0.0, 1.0, 1.0), (0.3, 1.0, 1.0), (1.0, 0.3, 0.3))}

    orangeDict = {'red': ((0.0, 1.0, 1.0), (0.3, 1.0, 1.0), (1.0, 0.3, 0.3)),
                  'green': ((0.0, 1.0, 1.0), (0.3, 0.65, .65), (1.0, .1, .1)),
                  'blue': ((0.0, 1.0, 1.0), (0.3, 0.0, 0.0), (1.0, 0.0, 0.0))}

    brownDict = {'red': ((0.0, 1.0, 1.0), (0.3, 0.57, 0.57), (1.0, 0.2, 0.2)),
                 'green': ((0.0, 1.0, 1.0), (0.3, 0.36, .36), (1.0, .12, .12)),
                 'blue': ((0.0, 1.0, 1.0), (0.3, 0.15, 0.15), (1.0, 0.01, 0.01))}

    if color == 'wind':
        cmap = LinearSegmentedColormap('blue', blueDict, 1000)
    elif color == 'solar':
        cmap = LinearSegmentedColormap('orange', orangeDict, 1000)
    else:
        cmap = LinearSegmentedColormap('brown', brownDict, 1000)

    for time in ['day', 'night']:
        # Load usages for given scheme and direction
        N_usages = np.load(outPath + '/Node_contrib_' + scheme + '_' + direction + '_' + time + '_' + color + '.npy')
        quantiles = np.load('./results/quantiles_' + str(scheme) + '_70128.npy')

        # Calculate colors of links
        linkUsages = np.sum(N_usages, 0) / quantiles
        col = [(cmap(l)) for l in linkUsages]

        # Create a new figure and plot network below
        fig = plt.figure(dpi=400, figsize=(0.85 * dcolwidth, 0.85 * 0.8 * dcolwidth))

        # color bar in bottom of figure
        ax1 = fig.add_axes([0.05, 0.08, 0.9, .08])
        cbl = mpl.colorbar.ColorbarBase(ax1, cmap, orientation='horizontal')

        # Label just above color bar
        if scheme == 'linear':
            xlabel = 'Most localised'
        elif scheme == 'square':
            xlabel = 'Synchronised'
        else:
            xlabel = 'Market'
        # ax1.set_xlabel(xlabel+' '+direction+r" usage $C_n/C^{\,99\%}$")
        ax1.set_xlabel(r'$\mathcal{K}_{l}/\mathcal{K}^T_l$')
        ax1.xaxis.set_label_position('top')

        ax2 = fig.add_axes([-0.05, 0.15, 1.1, 0.95])

        # Set color of nodes and draw all
        node_c = ["#000000" for node in N]
        nx.draw_networkx_nodes(G, pos, node_size=500, nodelist=nodelist, node_color=node_c, facecolor=(1, 1, 1))

        # Draw links colored by usage of node n
        edges = [(u, v) for (u, v, d) in G.edges(data=True)]
        edge_id = [d['id'] for (u, v, d) in G.edges(data=True)]

        color_sort = []
        for i in range(len(col)):
            color_sort.append(col[edge_id[i]])
        nx.draw_networkx_edges(G, pos, edgelist=edges, width=3.5, edge_color=color_sort, alpha=1)

        # Draw country names
        nx.draw_networkx_labels(G, pos, font_size=12, font_color='w', font_family='sans-serif')
        ax2.axis('off')

        # Save figure
        plt.savefig(figPath + "day/" + scheme + "/" + color + '_' + direction + '_' + time + ".png")
        plt.close()


def usagePlotter(direction):
    """
    Scatter plots of nodes' import/export usages of links saved to ./figures/.
    """
    legendNames = ['diagonal', r'$99\%$ quantile', 'avg. usage', 'usage']
    modes = ['linear', 'square']
    modeNames = ['localised', 'synchronised']
    names = ['usageS', 'usageW']
    colors = ['#ffa500', '#0000aa']
    for mode in modes:
        N = EU_Nodes_usage(mode + '.npz')
        F = np.load('./results/' + mode + '-flows.npy')
        Fmax = np.max(np.abs(F), 1)
        nodes = len(N)
        links = F.shape[0]

        usageS = np.load(outPath + mode + '_' + direction + '_' + 'usageS.npy')
        usageW = np.load(outPath + mode + '_' + direction + '_' + 'usageW.npy')
        if mode == 'square':
            usageB = np.load(outPath + mode + '_' + direction + '_' + 'usageB.npy')
            names.append('usageB')
            colors.append('#874a2b')

        for node in xrange(nodes):
            nodeLabel = N[node].label
            nodePath = figPath + 'usage/' + nodeLabel.tostring()
            if not os.path.exists(nodePath):
                os.makedirs(nodePath)
            for link in xrange(links):
                linkLabel = link_label(link, N)
                linkflow = abs(F[link, :])
                qq = get_q(abs(F[link]), .99)

                plt.figure()
                ax = plt.subplot()
                nBins = 90
                totUsage = np.zeros((nBins))
                for i, color in enumerate(names):
                    usages = eval(color)
                    usages = usages[link, node, :] / linkflow
                    F_vert = np.reshape(linkflow, (len(linkflow), 1))
                    exp_vert = np.reshape(usages, (len(usages), 1))
                    F_matrix = np.hstack([F_vert, exp_vert])
                    F_matrix[F_matrix[:, 0].argsort()]
                    H, bin_edges = binMaker(F_matrix, qq, lapse=70128, nbins=nBins)
                    plt.plot(bin_edges / qq, H[:, 1], '-', c=colors[i], lw=2)
                    totUsage += H[:, 1]
                plt.plot(bin_edges / qq, totUsage, '-', c="#aa0000", lw=2)

                plt.axis([0, 1, 0, 1])
                ax.set_xticks(np.linspace(0, 1, 11))
                plt.xlabel(r'$|F_l|/\mathcal{K}_l^T$')
                plt.ylabel(r'$\left\langle H_{ln} \right\rangle /|F_l|$')
                if mode == 'square':
                    modeName = modeNames[1]
                    plt.legend(('solar usage', 'wind usage', 'backup usage', 'total usage'), loc=1)
                else:
                    modeName = modeNames[0]
                    plt.legend(('solar usage', 'wind usage', 'total usage'), loc=1)
                plt.title(nodeLabel.tostring() + ' ' + modeName + ' ' + direction + ' flows on link ' + linkLabel)
                plt.savefig(nodePath + '/' + str(link) + '_' + modeName + '_' + direction + '.png', bbox_inches='tight')
                plt.close()


def link_level_bars(levels, usages, quantiles, scheme, direction, color, nnames, lnames, admat=None):
    """
    Bar plots of nodes' link usage of links at different levels.
    """
    if not admat:
        admat = np.genfromtxt('./settings/eadmat.txt')
    if color == 'solar':
        cmap = 'Oranges'
    elif color == 'wind':
        cmap = 'Blues'
    elif color == 'backup':
        cmap = 'Greys'
    nodes, links = usages.shape
    usageLevels = np.zeros((nodes, levels))
    usageLevelsNorm = np.zeros((nodes, levels))
    for node in range(nodes):
        nl = neighbor_levels(node, levels, admat)
        for lvl in range(levels):
            ll = link_level(nl, lvl, nnames, lnames)
            ll = np.array(ll, dtype='int')
            usageSum = sum(usages[node, ll])
            linkSum = sum(quantiles[ll])
            usageLevels[node, lvl] = usageSum / linkSum
            if lvl == 0:
                usageLevelsNorm[node, lvl] = usageSum
            else:
                usageLevelsNorm[node, lvl] = usageSum / usageLevelsNorm[node, 0]
        usageLevelsNorm[:, 0] = 1

    # plot all nodes
    usages = usageLevels.transpose()
    plt.figure(figsize=(11, 3))
    ax = plt.subplot()
    plt.pcolormesh(usages[:, loadOrder], cmap=cmap)
    plt.colorbar().set_label(label=r'$ \sum_l\, (C_{ln}) / \sum_l\, (\mathcal{K}^T_l)$', size=10)
    ax.set_yticks(np.linspace(.5, levels - .5, levels))
    ax.set_yticklabels(range(1, levels + 1))
    ax.yaxis.set_tick_params(width=0)
    ax.xaxis.set_tick_params(width=0)
    ax.set_xticks(np.linspace(1, nodes, nodes))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10)
    plt.ylabel('Link level')
    plt.savefig(figPath + '/levels/' + str(scheme) + '/' + 'total' + '_' + str(direction) + '_' + color + '.png', bbox_inches='tight')
    plt.close()

    # plot all nodes normalised to usage of first level
    usages = usageLevelsNorm.transpose()
    plt.figure(figsize=(11, 3))
    ax = plt.subplot()
    plt.pcolormesh(usages[:, loadOrder], cmap=cmap)
    plt.colorbar().set_label(label=r'$ \sum_l\, (C_{ln}) / \sum_l\, (\mathcal{K}^T_l)$', size=10)
    ax.set_yticks(np.linspace(.5, levels - .5, levels))
    ax.set_yticklabels(range(1, levels + 1))
    ax.yaxis.set_tick_params(width=0)
    ax.xaxis.set_tick_params(width=0)
    ax.set_xticks(np.linspace(1, nodes, nodes))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10)
    plt.ylabel('Link level')
    plt.savefig(figPath + '/levels/' + str(scheme) + '/' + 'total_norm_cont_' + str(direction) + '_' + color + '.png', bbox_inches='tight')
    plt.close()


def link_level_norm(levels, usages, quantiles, scheme, direction, color, nnames, lnames, admat=None):
    """
    Bar plots of nodes' link usage of links at different levels normed to the
    usage at the first level.
    """
    if not admat:
        admat = np.genfromtxt('./settings/eadmat.txt')
    if color == 'solar':
        cmap = 'Oranges'
    elif color == 'wind':
        cmap = 'Blues'
    elif color == 'backup':
        cmap = 'Greys'
    links, nodes, lapse = usages.shape
    usageLevels = np.zeros((nodes, levels))
    usageLevelsNorm = np.zeros((nodes, levels))
    for node in range(nodes):
        nl = neighbor_levels(node, levels, admat)
        for lvl in range(levels):
            ll = link_level(nl, lvl, nnames, lnames)
            ll = np.array(ll, dtype='int')
            usageSum = sum(sum(usages[ll, node, :]))
            linkSum = sum(quantiles[ll])
            usageLevels[node, lvl] = usageSum / linkSum
            if lvl == 0:
                usageLevelsNorm[node, lvl] = usageSum
            else:
                usageLevelsNorm[node, lvl] = usageSum / usageLevelsNorm[node, 0]
        usageLevelsNorm[:, 0] = 1

    # plot all nodes normalised to usage of first level
    usages = usageLevelsNorm.transpose()
    plt.figure(figsize=(11, 3))
    ax = plt.subplot()
    plt.pcolormesh(usages[:, loadOrder], cmap=cmap)
    plt.colorbar().set_label(label=r'$ \sum_l\, (H_{ln}(t)) / \sum_l\, (\mathcal{K}^T_l)$', size=10)
    ax.set_yticks(np.linspace(.5, levels - .5, levels))
    ax.set_yticklabels(range(1, levels + 1))
    ax.yaxis.set_tick_params(width=0)
    ax.xaxis.set_tick_params(width=0)
    ax.set_xticks(np.linspace(1, nodes, nodes))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10)
    plt.ylabel('Link level')
    plt.savefig(figPath + '/levels/' + str(scheme) + '/' + 'total_norm' + '_' + str(direction) + '_' + color + '.png', bbox_inches='tight')
    plt.close()


def link_level_hour(levels, usages, quantiles, scheme, direction, color, nnames, lnames, admat=None):
    """
    Make a color mesh of a node's average hourly usage of links at different
    levels.
    """
    if not admat:
        admat = np.genfromtxt('./settings/eadmat.txt')
    if color == 'solar':
        cmap = 'Oranges'
    elif color == 'wind':
        cmap = 'Blues'
    elif color == 'backup':
        cmap = 'Greys'
    links, nodes, lapse = usages.shape
    usages = np.reshape(usages, (links, nodes, lapse / 24, 24))
    totalHour = np.zeros((levels, 24))
    totalNormed = np.zeros((levels, 24))
    for node in range(nodes):
        nl = neighbor_levels(node, levels, admat)
        hourSums = np.zeros((levels, 24))
        for lvl in range(levels):
            ll = link_level(nl, lvl, nnames, lnames)
            ll = np.array(ll, dtype='int')
            meanSum = np.sum(np.mean(usages[ll, node], axis=1), axis=0)
            linkSum = sum(quantiles[ll])
            hourSums[lvl] = meanSum / linkSum
        totalHour += hourSums

        plt.figure(figsize=(9, 3))
        ax = plt.subplot()
        plt.pcolormesh(hourSums, cmap=cmap)
        plt.colorbar().set_label(label=r'$ \sum_l\, \left\langle H_{ln}(t) \right\rangle / \sum_l\, (\mathcal{K}^T_l)$', size=10)
        ax.set_yticks(np.linspace(.5, levels - .5, levels))
        ax.set_yticklabels(range(1, levels + 1))
        ax.yaxis.set_tick_params(width=0)
        ax.xaxis.set_tick_params(width=0)
        ax.set_xticks(np.linspace(.5, 23.5, 24))
        ax.set_xticklabels(np.array(np.linspace(1, 24, 24), dtype='int'), ha="center", va="top", fontsize=10)
        plt.ylabel('Link level')
        plt.axis([0, 24, 0, 5])
        plt.title(nnames[node] + ' ' + direction + ' ' + color)
        plt.savefig(figPath + '/hourly/' + str(scheme) + '/' + str(node) + '_' + color + '_' + direction + '.png', bbox_inches='tight')
        plt.close()

        hourSums = hourSums / np.sum(hourSums, axis=1)[:, None]
        totalNormed += hourSums
        plt.figure(figsize=(9, 3))
        ax = plt.subplot()
        plt.pcolormesh(hourSums, cmap=cmap)
        plt.colorbar().set_label(label=r'$ \sum_l\, \left\langle H_{ln}(t) \right\rangle / \sum_l\, (\mathcal{K}^T_l)$', size=10)
        ax.set_yticks(np.linspace(.5, levels - .5, levels))
        ax.set_yticklabels(range(1, levels + 1))
        ax.yaxis.set_tick_params(width=0)
        ax.xaxis.set_tick_params(width=0)
        ax.set_xticks(np.linspace(.5, 23.5, 24))
        ax.set_xticklabels(np.array(np.linspace(1, 24, 24), dtype='int'), ha="center", va="top", fontsize=10)
        plt.ylabel('Link level')
        plt.axis([0, 24, 0, 5])
        plt.title(nnames[node] + ' ' + direction + ' ' + color)
        plt.savefig(figPath + '/hourly/' + str(scheme) + '/normed/' + str(node) + '_' + color + '_' + direction + '.png', bbox_inches='tight')
        plt.close()

    # Plot average hourly usage
    totalHour /= nodes
    plt.figure(figsize=(9, 3))
    ax = plt.subplot()
    plt.pcolormesh(totalHour, cmap=cmap)
    plt.colorbar().set_label(label=r'$ \sum_l\, \left\langle H_{ln}(t) \right\rangle / \sum_l\, (\mathcal{K}^T_l)$', size=10)
    ax.set_yticks(np.linspace(.5, levels - .5, levels))
    ax.set_yticklabels(range(1, levels + 1))
    ax.yaxis.set_tick_params(width=0)
    ax.xaxis.set_tick_params(width=0)
    ax.set_xticks(np.linspace(.5, 23.5, 24))
    ax.set_xticklabels(np.array(np.linspace(1, 24, 24), dtype='int'), ha="center", va="top", fontsize=10)
    plt.ylabel('Link level')
    plt.axis([0, 24, 0, 5])
    plt.title('Average ' + direction + ' ' + color)
    plt.savefig(figPath + '/hourly/' + str(scheme) + '/total_' + color + '_' + direction + '.png', bbox_inches='tight')
    plt.close()

    totalNormed /= nodes
    plt.figure(figsize=(9, 3))
    ax = plt.subplot()
    plt.pcolormesh(totalNormed, cmap=cmap)
    plt.colorbar().set_label(label=r'$ \sum_l\, \left\langle H_{ln}(t) \right\rangle / \sum_l\, (\mathcal{K}^T_l)$', size=10)
    ax.set_yticks(np.linspace(.5, levels - .5, levels))
    ax.set_yticklabels(range(1, levels + 1))
    ax.yaxis.set_tick_params(width=0)
    ax.xaxis.set_tick_params(width=0)
    ax.set_xticks(np.linspace(.5, 23.5, 24))
    ax.set_xticklabels(np.array(np.linspace(1, 24, 24), dtype='int'), ha="center", va="top", fontsize=10)
    plt.ylabel('Link level')
    plt.axis([0, 24, 0, 5])
    plt.title('Average ' + direction + ' ' + color)
    plt.savefig(figPath + '/hourly/' + str(scheme) + '/normed/total_' + color + '_' + direction + '.png', bbox_inches='tight')
    plt.close()

if 'trace' in task:
    print('tracing')
    for mode in modes:
        print(str(mode))
        N = np.load('./results/' + mode + '_pm.npz', mmap_mode='r')
        F = abs(np.load('./results/' + mode + '-flows.npy'))
        quantiles = [get_q(abs(F[link]), .99) for link in range(len(F))]
        nodes = 30
        names = ['solar', 'wind']
        meanLoads = np.reshape(N['mean'], (nodes, 1))
        genS = N['normsolar']
        genW = N['normwind']
        genSum = genS + genW
        genSum[np.where(genSum == 0)] = 1
        if mode == 'square':
            genB = np.divide(N['balancing'], meanLoads)
            genSum += genB
            normGenB = genB / genSum
            names.append('backup')
        normGenS = genS / genSum
        normGenW = genW / genSum

        for direction in directions:
            print(str(direction))
            if not os.path.exists(outPath + mode + '_' + direction + '_' + 'usageS.npy'):
                if direction == 'combined':
                    Usages = np.load('./linkcolouring/old_' + mode + '_copper_link_mix_import_all_alpha=same.npy')
                    Usages += np.load('./linkcolouring/old_' + mode + '_copper_link_mix_export_all_alpha=same.npy')
                    Usages /= 2
                else:
                    Usages = np.load('./linkcolouring/old_' + mode + '_copper_link_mix_' + direction + '_all_alpha=same.npy')

                links, nodes, lapse = Usages.shape
                usageS = np.zeros((links, nodes, lapse))
                usageW = np.zeros((links, nodes, lapse))
                for l in xrange(links):
                    usageS[l] = Usages[l] * normGenS
                    usageW[l] = Usages[l] * normGenW
                np.save(outPath + mode + '_' + direction + '_' + 'usageS.npy', usageS)
                np.save(outPath + mode + '_' + direction + '_' + 'usageW.npy', usageW)
                if mode == 'square':
                    usageB = np.zeros((links, nodes, lapse))
                    for l in range(links):
                        usageB[l] = Usages[l] * normGenB
                    np.save(outPath + mode + '_' + direction + '_' + 'usageB.npy', usageB)
                Usages = None

                print('Solar')
                usageCalc(F, quantiles, usageS, nodes, links, 'solar')
                usageCalcDaily(F, quantiles, usageS, nodes, links, 'solar')
                print('Wind')
                usageCalc(F, quantiles, usageW, nodes, links, 'wind')
                usageCalcDaily(F, quantiles, usageW, nodes, links, 'wind')
                if mode == 'square':
                    print('Backup')
                    usageCalc(F, quantiles, usageB, nodes, links, 'backup')
                    usageCalcDaily(F, quantiles, usageB, nodes, links, 'backup')

            else:
                print('Solar')
                usage = np.load(outPath + mode + '_' + direction + '_' + 'usageS.npy')
                links, nodes, lapse = usage.shape
                usageCalc(F, quantiles, usage, nodes, links, 'solar')
                usageCalcDaily(F, quantiles, usage, nodes, links, 'solar')
                print('Wind')
                usage = np.load(outPath + mode + '_' + direction + '_' + 'usageW.npy')
                links, nodes, lapse = usage.shape
                usageCalc(F, quantiles, usage, nodes, links, 'wind')
                usageCalcDaily(F, quantiles, usage, nodes, links, 'wind')
                if mode == 'square':
                    print('Backup')
                    usage = np.load(outPath + mode + '_' + direction + '_' + 'usageB.npy')
                    links, nodes, lapse = usage.shape
                    usageCalc(F, quantiles, usage, nodes, links, 'backup')
                    usageCalcDaily(F, quantiles, usage, nodes, links, 'backup')


if 'plot' in task:
    if 'network' in task:
        N = EU_Nodes_usage()
        colors = ['solar', 'wind']
        for mode in modes:
            print('Mode: ' + mode)
            if mode == 'square':
                colors.append('backup')
            for direction in directions:
                print('Direction: ' + direction)
                for color in colors:
                    if 'total' in task:
                        print('Plotting total network figures')
                        drawnet_total(N, mode, direction, color)
                    elif 'day' in task:
                        print('Plotting day/night network figures')
                        drawnet_day(N, mode, direction, color)
                    else:
                        print('Plotting network figures')
                        drawnet_usage(N, mode, direction, color)

    if 'usage' in task:
        print('Plotting usage figures')
        p = Pool(len(directions))
        p.map(usagePlotter, directions)

    if 'levels' in task:
        levels = 5
        N = EU_Nodes_usage()
        colors = ['solar', 'wind']
        c = ['S', 'W']
        lnames = np.array(link_namer(N))
        nnames = np.array(node_namer(N))
        schemeNames = ['localised', 'synchronised']
        if 'norm' not in task:
            print('Plotting link levels')
            for i, scheme in enumerate(modes):
                if scheme == 'square':
                    colors.append('backup')
                name = schemeNames[i]
                quantiles = np.load('./results/quantiles_' + str(scheme) + '_70128.npy')
                for direction in directions:
                    for color in colors:
                            N_usages = np.load(outPath + 'Node_contrib_' + scheme + '_' + direction + '_' + color + '.npy')
                            link_level_bars(levels, N_usages, quantiles, name, direction, color, nnames, lnames)
        else:
            print('Plotting normed link levels')
            for i, scheme in enumerate(modes):
                if scheme == 'square':
                    colors.append('backup')
                    c.append('B')
                name = schemeNames[i]
                quantiles = np.load('./results/quantiles_' + str(scheme) + '_70128.npy')
                for direction in directions:
                    for j, color in enumerate(colors):
                        Usages = np.load(outPath + scheme + '_' + direction + '_' + 'usage' + c[j] + '.npy')
                        link_level_norm(levels, Usages, quantiles, name, direction, color, nnames, lnames)

    if 'hour' in task:
        print('Plotting hourly link levels')
        levels = 5
        N = EU_Nodes_usage()
        colors = ['solar', 'wind']
        c = ['S', 'W']
        lnames = np.array(link_namer(N))
        nnames = np.array(node_namer(N))
        schemeNames = ['localised', 'synchronised']
        for i, scheme in enumerate(modes):
            if scheme == 'square':
                colors.append('backup')
                c.append('B')
            name = schemeNames[i]
            quantiles = np.load('./results/quantiles_' + str(scheme) + '_70128.npy')
            for direction in directions:
                for j, color in enumerate(colors):
                    Usages = np.load(outPath + scheme + '_' + direction + '_' + 'usage' + c[j] + '.npy')
                    link_level_hour(levels, Usages, quantiles, name, direction, color, nnames, lnames)

if 'sanity' in task:
    for mode in modes:
        quantiles = np.load('./results/quantiles_' + str(mode) + '_70128.npy')
        for direction in directions:
            S = np.load('./results/vector/Node_contrib_' + mode + '_' + direction + '_solar.npy')
            S += np.load('./results/vector/Node_contrib_' + mode + '_' + direction + '_wind.npy')
            if mode == 'square':
                S += np.load('./results/vector/Node_contrib_' + mode + '_' + direction + '_backup.npy')
            N = np.load('./results/Node_contrib_' + mode + '_' + direction + '_70128.npy')

            error = abs(S - N) / N * 100
            weightedError = abs(S - N) / N * quantiles / np.mean(quantiles) * 100
            means = np.mean(error, axis=1)
            weightedMeans = np.mean(weightedError, axis=1)
            stds = np.std(error, axis=1)
            nodeMean = np.mean(means)
            weightedNodeMean = np.mean(weightedMeans)

            x = np.linspace(.5, 29.5, 30)
            if mode == 'linear': title = 'localised'
            if mode == 'square': title = 'synchronised'
            plt.figure()
            ax = plt.subplot()
            plt.errorbar(x, means[loadOrder], yerr=stds * 0, marker='s', lw=0, elinewidth=1)
            plt.plot([0, 30], [nodeMean, nodeMean], '--k', lw=2)
            plt.title(title + ' ' + direction + ', sum of colors vs. total network usage')
            plt.ylabel('Mean link deviation in %')
            ax.set_xticks(np.linspace(1, 30, 30))
            ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=9)
            plt.axis([0, 30, min(means) - (.1 * min(means)), max(means) + (.1 * max(means))])
            plt.legend(('individual country', 'mean of countries'), loc=2, ncol=2)
            plt.savefig(figPath + 'error/' + title + '_' + direction + '_.png', bbox_inches='tight')

            plt.figure()
            ax = plt.subplot()
            plt.errorbar(x, weightedMeans[loadOrder], yerr=stds * 0, marker='s', lw=0, elinewidth=1)
            plt.plot([0, 30], [weightedNodeMean, weightedNodeMean], '--k', lw=2)
            plt.title(title + ' ' + direction + ', sum of colors vs. total network usage')
            plt.ylabel(r'Weighed mean link deviation in % normalised to $\left\langle \mathcal{K}^T \right\rangle$')
            ax.set_xticks(np.linspace(1, 30, 30))
            ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=9)
            plt.axis([0, 30, min(weightedMeans) - (.1 * min(weightedMeans)), max(weightedMeans) + (.1 * max(weightedMeans))])
            plt.legend(('individual country', 'mean of countries'), loc=2, ncol=2)
            plt.savefig(figPath + 'error/' + 'weighted_' + title + '_' + direction + '_.png', bbox_inches='tight')
            plt.close()
