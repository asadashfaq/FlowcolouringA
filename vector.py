import sys
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
from pylab import plt
import matplotlib as mpl
from aurespf.tools import *
from EUgrid import EU_Nodes_usage, EU_Nodes_regions, EU_Nodes_superRegions
from link_namer import node_namer, link_dict
from functions import *

"""
Vector flow tracing on top of the up/down stream flow tracing algorithm.
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


def usageCalc(F, quantiles, Usages, nodes, links, name):
    """
    Calculate usages and save to file
    """
    Node_contributions = np.zeros((nodes, links))  # empty array for calculated usages
    for node in range(nodes):
        print node
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


def drawnet_usage(N=None, scheme='linear', direction='combined', color='solar'):
    """
    Make network figures of a node's usage of links for both import, export and
    combined. Adapted from drawnet() in aurespf.plotting
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
    blueDict = {'red': ((0.0, 1.0, 1.0), (0.2, 0.0, 0.0), (0.9, 0.0, 0.0), (1.0, 0.0, 0.0)),
                'green': ((0.0, 1.0, 1.0), (0.2, 0.0, 0.0), (0.9, 0.0, 0.0), (1.0, 0.0, 0.0)),
                'blue': ((0.0, 1.0, 1.0), (0.2, 1.0, 1.0), (0.9, .3, .3), (1.0, 0.3, 0.3))}

    orangeDict = {'red': ((0.0, 1.0, 1.0), (0.15, 1.0, 1.0), (0.9, 0.3, 0.3), (1.0, 0.3, 0.3)),
                  'green': ((0.0, 1.0, 1.0), (0.15, 0.65, .65), (0.9, .1, .1), (1.0, .1, .1)),
                  'blue': ((0.0, 1.0, 1.0), (0.15, 0.0, 0.0), (0.9, 0.0, 0.0), (1.0, 0.0, 0.0))}

    brownDict = {'red': ((0.0, 1.0, 1.0), (0.1, 0.57, 0.57), (0.9, 0.2, 0.2), (1.0, 0.2, 0.2)),
                 'green': ((0.0, 1.0, 1.0), (0.1, 0.36, .36), (0.9, .12, .12), (1.0, .12, .12)),
                 'blue': ((0.0, 1.0, 1.0), (0.1, 0.15, 0.15), (0.9, 0.01, 0.01), (1.0, 0.01, 0.01))}

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
        ax1.set_xlabel(r'$\mathcal{K}^T_{ln}/\mathcal{K}^T_l$')
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
        plt.savefig(figPath + "network/" + scheme + "/" + str(n.id) + '_' + str(direction) + '_' + color + ".png")
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
            if direction == 'combined':
                Usages = np.load('./linkcolouring/old_' + mode + '_copper_link_mix_import_all_alpha=same.npy')
                Usages2 = np.load('./linkcolouring/old_' + mode + '_copper_link_mix_export_all_alpha=same.npy')
                Usages += Usages2
                Usages2 = None
                Usages /= 2
            else:
                Usages = np.load('./linkcolouring/old_' + mode + '_copper_link_mix_' + direction + '_all_alpha=same.npy')

            links, nodes, lapse = Usages.shape
            usageS = np.zeros((links, nodes, lapse))
            usageW = np.zeros((links, nodes, lapse))
            for l in xrange(links):
                usageS[l] = Usages[l] * normGenS
                usageW[l] = Usages[l] * normGenW
            if mode == 'square':
                usageB = np.zeros((links, nodes, lapse))
                for l in range(links):
                    usageB[l] = Usages[l] * normGenB
            Usages = None

            print('Solar')
            usageCalc(F, quantiles, usageS, nodes, links, 'solar')
            print('Wind')
            usageCalc(F, quantiles, usageW, nodes, links, 'wind')
            if mode == 'square':
                print('Backup')
                usageCalc(F, quantiles, usageB, nodes, links, 'backup')


if 'plot' in task:
    print('Plotting network figures')
    N = EU_Nodes_usage()
    colors = ['solar', 'wind']
    for mode in modes:
        print('Mode: ' + mode)
        if mode == 'square':
            colors.append('backup')
        for direction in directions:
            print('Direction: ' + direction)
            for color in colors:
                drawnet_usage(N, mode, direction, color)
