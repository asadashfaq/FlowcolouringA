#! /usr/bin/env/python
from __future__ import division
import os
import sys
from multiprocessing import Pool
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import shapefile
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import aurespf.solvers as au
from aurespf.tools import get_q, AtoKh_old
from EUgrid import *
from figutils import *
from link_namer import *
from link_colour_less import track_flows, get_link_direction
from new_linkcolouralgorithm_less import track_link_usage_total
from functions import binMaker, bin_prob, bin_CDF, node_contrib

"""
Flow tracing on heterogeneous networks

The script can be called using the following command line arguments:
- solve:    Solve power flows
- trace:    Run scalar vector flow tracing
- usage:    Calculate link usages
- vector:   Run vector flow tracing and calculate link usages
- plot:     plot various figures. See below.

Plotting:
- plot ag:              plot alphas and gammas for various betas
- plot map:             plot maps of beta layouts
- plot network:         plot network usage for every color, scheme, node, beta
- plot network total:   plot total network usage for every scheme, node, beta
- plot network day:     plot daily network usage for every color, scheme, node, beta
- plot levels:          plot usage of links at different levels
- plot hour:            plot hourly usage of links at different levels
"""

if len(sys.argv) < 2:
    raise Exception('Not enough inputs')
else:
    task = str(sys.argv[1:])

schemes = ['linear', 'square']
directions = ['import', 'export', 'combined']
lapse = 70128
N_bins = 90
nNodes = 30
B = [0.5, 1.5, 2.5, 3.5, 4.5]  # range(11)
# B = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
meanEU = 345327.47685659607
inPath = './results/heterogen/input/'
resPath = './results/heterogen/'
figPath = './figures/heterogen/'
colorPath = './linkcolouring/heterogen/'

# Capacity factors
cfW = {'BE': 2.7778, 'FR': 3.125, 'BG': 6.6667, 'DK': 2.4390, 'HR': 5.5556, 'DE': 2.7778,
       'HU': 5.2632, 'FI': 3.7037, 'BA': 6.25, 'NL': 2.5, 'PT': 4.5455, 'NO': 3.125,
       'LV': 3.4483, 'LT': 3.333, 'LU': 4.0, 'RO': 6.6667, 'PL': 3.125, 'CH': 7.1429,
       'GR': 7.1429, 'EE': 3.5714, 'IT': 5.0, 'CZ': 4.7619, 'GB': 2.381, 'IE': 2.1739,
       'ES': 5.2632, 'RS': 6.6667, 'SK': 5.8824, 'SI': 6.6667, 'SE': 3.0303, 'AT': 5.2632}
cfS = {'AT': 5.8824, 'BA': 4.7619, 'BE': 7.1429, 'BG': 4.5455, 'CH': 5.8824, 'CZ': 6.25,
       'DE': 6.6667, 'DK': 6.6667, 'EE': 7.6923, 'ES': 4.0, 'FI': 7.6923, 'FR': 5.0,
       'GB': 6.6667, 'GR': 4.1667, 'HR': 5.0, 'HU': 5.5556, 'IE': 7.6923, 'IT': 4.3478,
       'LT': 7.1429, 'LU': 7.1429, 'LV': 7.6923, 'NL': 6.6667, 'NO': 7.6923, 'PL': 6.6667,
       'PT': 4.3478, 'RO': 5.0, 'RS': 5.2632, 'SE': 7.1429, 'SI': 5.5556, 'SK': 5.8824}


def tablePrint():
    """
    Print capacity facors for LaTeX tablePrint
    """
    def printer(id):
        return fullNames[id] + ' & ' \
            + str(np.round(1 / cfW[loadNames[id]], 2)) + ' & ' \
            + str(np.round(1 / cfS[loadNames[id]], 2))
    for i in range(10):
        print printer(i) + ' & ' + printer(i + 10) + ' & ' + printer(i + 20) + ' \\\\'


def pathCheck(path):
    """
    check if path exists. If not create path
    """
    if not os.path.exists(path):
        os.makedirs(path)


def calcAG(N, b, alpha=0.7):
    """
    Calculate individual alpha and gamma from beta factor
    """
    for n in N:
        gs = (1 / cfS[str(n.label)]) ** b * meanEU / sum([(1 / cfS[str(m.label)]) ** b * m.mean for j, m in enumerate(N)])
        gw = (1 / cfW[str(n.label)]) ** b * meanEU / sum([(1 / cfW[str(m.label)]) ** b * m.mean for j, m in enumerate(N)])
        n.set_gamma(alpha * gw + (1 - alpha) * gs)
        n.set_alpha(alpha * gw / (alpha * gw + (1 - alpha) * gs))
    return N


def solveFlows(d):
    scheme = d[0]
    b = d[1]
    N = EU_Nodes_usage()
    N = calcAG(N, b)
    N, F = au.solve(N, mode=scheme + ' copper', lapse=lapse)
    N.save_nodes(scheme, path=resPath + 'b_' + str(b) + '_')
    np.save(resPath + scheme + '_b_' + str(b) + '_flows', F)


def traceFlow(d):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """
    scheme = d[0]
    b = d[1]
    N = EU_Nodes_usage('../' + resPath + 'b_' + str(b) + '_' + scheme + '.npz')
    F = np.load(resPath + scheme + '_b_' + str(b) + '_flows.npy')

    """
    N2 is a new nodes object containing individual powermixes for import and
    export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
    """
    N2, power_mixes_total = track_flows(N, F, lapse=lapse)

    """
    track_link_usage_total tracks each nodes usage of all links. The results
    are saved to files '..._links_ex_...' and '..._links_im_...'.
    """
    boxplot, boxplotlabel = track_link_usage_total(N2, F, mode=scheme, alph='same', lapse=lapse, heterogen=b)


def calcCont(d):
    """
    Calculate nodes' contribution for a single mode and save results to file.
    """
    scheme = d[0]
    b = d[1]
    F = abs(np.load(resPath + scheme + '_b_' + str(b) + '_flows.npy'))
    quantiles = [get_q(abs(F[link, :lapse]), .99) for link in range(len(F))]

    for direction in directions:
        if direction == 'combined':
            Usages = np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_import.npy')
            Usages += np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_export.npy')
            Usages /= 2
        else:
            Usages = np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_' + direction + '.npy')

        Node_contributions = np.zeros((nNodes, len(F)))  # empty array for calculated usages
        for node in range(nNodes):
            for link in range(len(F)):
                # Stacking and sorting data
                F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
                exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
                F_matrix = np.hstack([F_vert, exp_vert])  # [flow, usage]
                F_matrix[F_matrix[:, 0].argsort()]

                H, bin_edges = binMaker(F_matrix, quantiles[link], lapse, N_bins)
                Node_contributions[node, link] = node_contrib(H, bin_edges, linkID=link)

        np.save(resPath + 'N_cont_' + scheme + '_' + direction + '_b_' + str(b) + '.npy', Node_contributions)
        if direction == 'import':
            np.save(resPath + 'quant_' + scheme + '_b_' + str(b) + '.npy', quantiles)


def vCalcCont(F, quantiles, Usages, nodes, links, name, b, scheme, direction):
    """
    Calculate usages and save to file
    """
    Node_contributions = np.zeros((nodes, links))
    for node in range(nodes):
        for link in range(links):
            F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
            exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
            F_matrix = np.hstack([F_vert, exp_vert])
            F_matrix[F_matrix[:, 0].argsort()]
            H, bin_edges = binMaker(F_matrix, quantiles[link], lapse)
            Node_contributions[node, link] = node_contrib(H, bin_edges, linkID=link)
    np.save(resPath + 'N_cont_' + scheme + '_' + direction + '_b_' + str(b) + '_' + str(name) + '.npy', Node_contributions)
    return


def vCalcContDaily(iF, quantiles, iUsages, nodes, links, name, b, scheme, direction):
    """
    Calculate usages for daytime and nighttime and save to file
    """
    days = int(lapse / 24)
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
        np.save(resPath + 'N_cont_' + scheme + '_' + direction + '_b_' + str(b) + '_' + d + '_' + str(name) + '.npy', Node_contributions)
    return


def vectorTrace(d):
    """
    Vector flow tracing
    """
    scheme = d[0]
    b = d[1]
    N = np.load(resPath + 'b_' + str(b) + '_' + scheme + '.npz', mmap_mode='r')
    F = abs(np.load(resPath + scheme + '_b_' + str(b) + '_flows.npy'))
    quantiles = [get_q(abs(F[link]), .99) for link in range(len(F))]
    nodes = 30
    names = ['solar', 'wind']
    meanLoads = np.reshape(N['mean'], (nodes, 1))
    genS = N['normsolar']
    genW = N['normwind']
    genSum = genS + genW
    genSum[np.where(genSum == 0)] = 1
    if scheme == 'square':
        genB = np.divide(N['balancing'], meanLoads)
        genSum += genB
        normGenB = genB / genSum
        names.append('backup')
    normGenS = genS / genSum
    normGenW = genW / genSum

    for direction in directions:
        print(str(direction))
        if not os.path.exists(resPath + scheme + '_' + direction + '_b_' + str(b) + '_' + 'usageS.npy'):
            if direction == 'combined':
                Usages = np.load(colorPath + scheme + '-b-' + str(b) + '_link_mix_import.npy')
                Usages += np.load(colorPath + scheme + '-b-' + str(b) + '_link_mix_export.npy')
                Usages /= 2
            else:
                Usages = np.load(colorPath + scheme + '-b-' + str(b) + '_link_mix_' + direction + '.npy')

            links, nodes, lapse = Usages.shape
            usageS = np.zeros((links, nodes, lapse))
            usageW = np.zeros((links, nodes, lapse))
            for l in xrange(links):
                usageS[l] = Usages[l] * normGenS
                usageW[l] = Usages[l] * normGenW
            np.save(colorPath + scheme + '_' + direction + '_b_' + str(b) + '_' + 'usageS.npy', usageS)
            np.save(colorPath + scheme + '_' + direction + '_b_' + str(b) + '_' + 'usageW.npy', usageW)
            if scheme == 'square':
                usageB = np.zeros((links, nodes, lapse))
                for l in range(links):
                    usageB[l] = Usages[l] * normGenB
                np.save(colorPath + scheme + '_' + direction + '_b_' + str(b) + '_' + 'usageB.npy', usageB)
            Usages = None

            print('Solar')
            vCalcCont(F, quantiles, usageS, nodes, links, 'solar', b, scheme, direction)
            vCalcContDaily(F, quantiles, usageS, nodes, links, 'solar', b, scheme, direction)
            print('Wind')
            vCalcCont(F, quantiles, usageW, nodes, links, 'wind', b, scheme, direction)
            vCalcContDaily(F, quantiles, usageW, nodes, links, 'wind', b, scheme, direction)
            if scheme == 'square':
                print('Backup')
                vCalcCont(F, quantiles, usageB, nodes, links, 'backup', b, scheme, direction)
                vCalcContDaily(F, quantiles, usageB, nodes, links, 'backup', b, scheme, direction)

        else:
            print('Solar')
            usage = np.load(colorPath + scheme + '_' + direction + '_b_' + str(b) + '_' + 'usageS.npy')
            links, nodes, lapse = usage.shape
            vCalcCont(F, quantiles, usage, nodes, links, 'solar', b, scheme, direction)
            vCalcContDaily(F, quantiles, usage, nodes, links, 'solar', b, scheme, direction)
            print('Wind')
            usage = np.load(colorPath + scheme + '_' + direction + '_b_' + str(b) + '_' + 'usageW.npy')
            links, nodes, lapse = usage.shape
            vCalcCont(F, quantiles, usage, nodes, links, 'wind', b, scheme, direction)
            vCalcContDaily(F, quantiles, usage, nodes, links, 'wind', b, scheme, direction)
            if scheme == 'square':
                print('Backup')
                usage = np.load(colorPath + scheme + '_' + direction + '_b_' + str(b) + '_' + 'usageB.npy')
                links, nodes, lapse = usage.shape
                vCalcCont(F, quantiles, usage, nodes, links, 'backup', b, scheme, direction)
                vCalcContDaily(F, quantiles, usage, nodes, links, 'backup', b, scheme, direction)


def plot_europe_map(country_weights, b=None, ax=None):
    """
    Plot a map from shapefiles with coutnries colored by gamma. Inspired by Magnus.
    """
    if ax == None:
        ax = plt.subplot(111)
    m = Basemap(llcrnrlon=-10., llcrnrlat=30., urcrnrlon=50., urcrnrlat=72.,
                projection='lcc', lat_1=40., lat_2=60., lon_0=20.,
                resolution='l', area_thresh=1000.,
                rsphere=(6378137.00, 6356752.3142))
    m.drawcoastlines(linewidth=0)
    r = shapefile.Reader(r'settings/ne_10m_admin_0_countries/ne_10m_admin_0_countries')
    all_shapes = r.shapes()
    all_records = r.records()
    shapes = []
    records = []
    for country in all_countries:
        shapes.append(all_shapes[shapefile_index[country]])
        records.append(all_records[shapefile_index[country]])

    country_count = 0
    for record, shape in zip(records, shapes):
        lons, lats = zip(*shape.points)
        data = np.array(m(lons, lats)).T
        if len(shape.parts) == 1:
            segs = [data, ]
        else:
            segs = []
            for i in range(1, len(shape.parts)):
                index = shape.parts[i - 1]
                index2 = shape.parts[i]
                segs.append(data[index:index2])
            segs.append(data[index2:])
        lines = LineCollection(segs, antialiaseds=(1,))
        lines.set_facecolor(cmap(country_weights[country_count]))
        lines.set_edgecolors('k')
        lines.set_linewidth(0.3)
        ax.add_collection(lines)
        country_count += 1
    if b: plt.text(1.5e5, 4e6, r'$\beta = ' + str(b) + r'$', fontsize=12)


def drawnet_usage(N=None, scheme='linear', direction='combined', color='solar', b=1):
    """
    Make network figures of a node's usage of links for both import, export and
    combined. A figure for each color is created.
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

    if color == 'wind':
        cmap = LinearSegmentedColormap('blue', blueDict, 1000)
    elif color == 'solar':
        cmap = LinearSegmentedColormap('orange', orangeDict, 1000)
    elif color == 'backup':
        cmap = LinearSegmentedColormap('brown', brownDict, 1000)
    elif color == '':
        cmap = LinearSegmentedColormap('blue', blueDict2, 1000)

    # Load usages for given scheme and direction
    if color != '':
        colStr = '_' + color
    else:
        colStr = color
    N_usages = np.load(resPath + '/N_cont_' + scheme + '_' + direction + '_b_' + str(b) + colStr + '.npy')
    quantiles = np.load(resPath + 'quant_' + str(scheme) + '_b_' + str(b) + '.npy')

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
            xlabel = 'Localized'
        elif scheme == 'square':
            xlabel = 'Synchronized'

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
        thisPath = figPath + 'B_' + str(b) + '/network/' + scheme + '/'
        pathCheck(thisPath)
        plt.savefig(thisPath + str(n.id) + colStr + '_' + direction + '.png', bbox_inches='tight')
        plt.close()


def drawnet_total(N=None, scheme='linear', direction='combined', color='solar', b=1):
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

    if color == 'wind':
        cmap = LinearSegmentedColormap('blue', blueDict, 1000)
    elif color == 'solar':
        cmap = LinearSegmentedColormap('orange', orangeDict, 1000)
    elif color == 'backup':
        cmap = LinearSegmentedColormap('brown', brownDict, 1000)
    elif color == '':
        cmap = LinearSegmentedColormap('blue', blueDict2, 1000)

    # Load usages for given scheme and direction
    if color != '':
        colStr = '_' + color
    else:
        colStr = color
    N_usages = np.load(resPath + '/N_cont_' + scheme + '_' + direction + '_b_' + str(b) + colStr + '.npy')
    quantiles = np.load(resPath + 'quant_' + str(scheme) + '_b_' + str(b) + '.npy')

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
        xlabel = 'Localized'
    elif scheme == 'square':
        xlabel = 'Synchronized'
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
    thisPath = figPath + 'B_' + str(b) + '/network/' + scheme + '/'
    pathCheck(thisPath)
    plt.savefig(thisPath + 'total' + colStr + '_' + direction + '.png', bbox_inches='tight')
    plt.close()


def drawnet_day(N=None, scheme='linear', direction='combined', color='solar', b=1):
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

    if color == 'wind':
        cmap = LinearSegmentedColormap('blue', blueDict, 1000)
    elif color == 'solar':
        cmap = LinearSegmentedColormap('orange', orangeDict, 1000)
    else:
        cmap = LinearSegmentedColormap('brown', brownDict, 1000)

    quantiles = np.load(resPath + 'quant_' + str(scheme) + '_b_' + str(b) + '.npy')
    for time in ['day', 'night']:
        # Load usages for given scheme and direction
        N_usages = np.load(resPath + 'N_cont_' + scheme + '_' + direction + '_b_' + str(b) + '_' + time + '_' + color + '.npy')

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
            xlabel = 'Localized'
        elif scheme == 'square':
            xlabel = 'Synchronized'
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
        thisPath = figPath + 'B_' + str(b) + '/day/' + scheme + '/'
        pathCheck(thisPath)
        plt.savefig(thisPath + color + '_' + direction + '_' + time + '.png')
        plt.close()


def link_level_bars(levels, usages, quantiles, scheme, direction, color, nnames, lnames, admat=None, b=1):
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
    else:
        cmap = 'OrRd'
    if color != '':
        colStr = '_' + color
    else:
        colStr = color

    thisPath = figPath + 'B_' + str(b) + '/levels/' + scheme + '/'
    pathCheck(thisPath)

    nodes, links = usages.shape
    usageLevels = np.zeros((nodes, levels))
    for node in range(nodes):
        nl = neighbor_levels(node, levels, admat)
        for lvl in range(levels):
            ll = link_level(nl, lvl, nnames, lnames)
            ll = np.array(ll, dtype='int')
            usageSum = sum(usages[node, ll])
            linkSum = sum(quantiles[ll])
            usageLevels[node, lvl] = usageSum / linkSum
    usages = usageLevels.transpose()
    plt.figure(figsize=(11, 3))
    ax = plt.subplot()
    plt.pcolormesh(usages[:, loadOrder], cmap=cmap)
    plt.colorbar().set_label(label=r'$U_n^{(l)}$', size=11)
    ax.set_yticks(np.linspace(.5, levels - .5, levels))
    ax.set_yticklabels(range(1, levels + 1))
    ax.yaxis.set_tick_params(width=0)
    ax.xaxis.set_tick_params(width=0)
    ax.set_xticks(np.linspace(1, nodes, nodes))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10)
    plt.ylabel('Link level')
    plt.savefig(thisPath + 'total' + colStr + '_' + direction + '.png', bbox_inches='tight')
    plt.close()


def plotTotal(scheme, color):
    """
    Plot total network usage for every color, scheme and beta
    """
    totUsageIm = np.zeros((nNodes, len(B)))
    totUsageEx = np.zeros((nNodes, len(B)))
    if color == 'solar':
        cmap = 'Oranges'
    elif color == 'wind':
        cmap = 'Blues'
    elif color == 'backup':
        cmap = 'Greys'
    if color != '':
        colStr = '_' + color
    else:
        colStr = color
        cmap = 'OrRd'
    for i, b in enumerate(B):
        usagesIm = np.load(resPath + '/N_cont_' + scheme + '_import_b_' + str(b) + colStr + '.npy')
        totUsageIm[:, i] = np.sum(usagesIm, 1) / node_mean_load
        usagesEx = np.load(resPath + '/N_cont_' + scheme + '_export_b_' + str(b) + colStr + '.npy')
        totUsageEx[:, i] = np.sum(usagesEx, 1) / node_mean_load

    plt.figure(figsize=(12, 5))
    ax1 = plt.subplot(121)
    ax1.set_xticks(np.linspace(0.5, 10.5, 11))
    ax1.set_xticklabels(range(11), fontsize=8)
    ax1.set_yticks(np.linspace(.5, 29.5, 30))
    ax1.set_yticklabels(loadNames, ha="right", va="center", fontsize=8)
    ax1.xaxis.set_tick_params(width=0)
    ax1.yaxis.set_tick_params(width=0)
    plt.pcolormesh(totUsageIm[loadOrder], cmap=cmap)
    cb1 = plt.colorbar()
    cb1.solids.set_edgecolor('face')
    cb1.set_label(label='Network usage', size=11)
    plt.xlabel(r'$\beta$')

    ax2 = plt.subplot(122)
    ax2.set_xticks(np.linspace(0.5, 10.5, 11))
    ax2.set_xticklabels(range(11), fontsize=8)
    ax2.set_yticks(np.linspace(.5, 29.5, 30))
    ax2.set_yticklabels(loadNames, ha="right", va="center", fontsize=8)
    ax2.xaxis.set_tick_params(width=0)
    ax2.yaxis.set_tick_params(width=0)
    plt.pcolormesh(totUsageEx[loadOrder], cmap=cmap)
    cb2 = plt.colorbar()
    cb2.solids.set_edgecolor('face')
    cb2.set_label(label='Network usage', size=11)
    plt.xlabel(r'$\beta$')
    plt.savefig(figPath + '/total/' + 'netUsage' + '_' + scheme + colStr + '.pdf', bbox_inches='tight')
    plt.close()


def link_level_hour(levels, usages, quantiles, scheme, direction, color, nnames, lnames, admat=None, b=1):
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
    if color != '':
        colStr = '_' + color
    else:
        colStr = color
        cmap = 'OrRd'
    links, nodes, lapse = usages.shape
    usages = np.reshape(usages, (links, nodes, lapse / 24, 24))
    totalHour = np.zeros((levels, 24))
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
        plt.colorbar().set_label(label=r'$U_n^{(l)}$', size=11)
        ax.set_yticks(np.linspace(.5, levels - .5, levels))
        ax.set_yticklabels(range(1, levels + 1))
        ax.yaxis.set_tick_params(width=0)
        ax.xaxis.set_tick_params(width=0)
        ax.set_xticks(np.linspace(.5, 23.5, 24))
        ax.set_xticklabels(np.array(np.linspace(1, 24, 24), dtype='int'), ha="center", va="top", fontsize=10)
        plt.ylabel('Link level')
        plt.axis([0, 24, 0, levels])
        plt.title(nnames[node] + ' ' + direction + ' ' + color)

        thisPath = figPath + 'B_' + str(b) + '/hourly/' + scheme + '/'
        pathCheck(thisPath)
        plt.savefig(thisPath + str(node) + colStr + '_' + direction + '.png', bbox_inches='tight')
        plt.close()

    # Plot average hourly usage
    totalHour /= nodes
    plt.figure(figsize=(9, 3))
    ax = plt.subplot()
    plt.pcolormesh(totalHour, cmap=cmap)
    plt.colorbar().set_label(label=r'$U_n^{(l)}$', size=11)
    ax.set_yticks(np.linspace(.5, levels - .5, levels))
    ax.set_yticklabels(range(1, levels + 1))
    ax.yaxis.set_tick_params(width=0)
    ax.xaxis.set_tick_params(width=0)
    ax.set_xticks(np.linspace(.5, 23.5, 24))
    ax.set_xticklabels(np.array(np.linspace(1, 24, 24), dtype='int'), ha="center", va="top", fontsize=10)
    plt.ylabel('Link level')
    plt.axis([0, 24, 0, levels])
    plt.savefig(thisPath + 'total' + colStr + '_' + direction + '.png', bbox_inches='tight')
    plt.close()


def agPlot():
    """
    Plot evolutino of alphas and gammas for increasing beta
    """
    B = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    alphas = np.zeros((nNodes, len(B)))
    gammas = np.zeros((nNodes, len(B)))
    for i, b in enumerate(B):
        N = np.load('./results/heterogen/b_' + str(b) + '_linear.npz', mmap_mode='r')
        alphas[:, i] = N['alpha']
        gammas[:, i] = N['gamma']

    plt.figure(figsize=(12, 5))
    ax1 = plt.subplot(121)
    ax1.set_xticks(np.linspace(0.5, 10.5, 11))
    ax1.set_xticklabels(B, fontsize=8)
    ax1.set_yticks(np.linspace(.5, 29.5, 30))
    ax1.set_yticklabels(loadNames, ha="right", va="center", fontsize=8)
    ax1.xaxis.set_tick_params(width=0)
    ax1.yaxis.set_tick_params(width=0)
    plt.pcolormesh(alphas[loadOrder], cmap='coolwarm_r')
    cb1 = plt.colorbar()
    cb1.solids.set_edgecolor('face')
    cb1.set_label(label=r'$\alpha$', size=13)
    plt.xlabel(r'$\beta$')

    maxGamma = np.max(gammas)
    point = 1 / maxGamma
    redgreendict = {'red': [(0.0, 1, 1), (point, 1.0, 1.0), (1.0, 0.0, 0.0)],
                    'green': [(0.0, 0.4, 0.4), (point, 1.0, 1.0), (1.0, 0.8, 0.8)],
                    'blue': [(0.0, 0.4, 0.4), (point, 1.0, 1.0), (1.0, 0.0, 0.0)]}
    cmap = LinearSegmentedColormap('redgreen', redgreendict, 1000)

    ax2 = plt.subplot(122)
    ax2.set_xticks(np.linspace(0.5, 10.5, 11))
    ax2.set_xticklabels(B, fontsize=8)
    ax2.set_yticks(np.linspace(.5, 29.5, 30))
    ax2.set_yticklabels(loadNames, ha="right", va="center", fontsize=8)
    ax2.xaxis.set_tick_params(width=0)
    ax2.yaxis.set_tick_params(width=0)
    plt.pcolormesh(gammas[loadOrder], cmap=cmap)
    cb2 = plt.colorbar()
    cb2.solids.set_edgecolor('face')
    cb2.set_label(label=r'$\gamma$', size=13)
    plt.xlabel(r'$\beta$')
    plt.savefig(figPath + 'alpha_gamma' + '.pdf', bbox_inches='tight')
    plt.close()

# Parameters for parallel calling of functions
d = []
for i in range(2 * len(B)):
    d.append([schemes[i // len(B)], B[i % len(B)]])

if 'solve' in task:
    print 'Solving power flows'
    p = Pool(4)
    p.map(solveFlows, d)

if 'trace' in task:
    print 'Flow tracing'
    p = Pool(3)
    p.map(traceFlow, d)

if 'usage' in task:
    print 'Calculating link usages'
    p = Pool(4)
    p.map(calcCont, d)

if 'vector' in task:
    print 'Vector flow tracing'
    p = Pool(4)
    p.map(vectorTrace, d)

if 'plot' in task:
    if 'ag' in task:
        print 'Plotting alphas and gammas'
        agPlot()

    if 'map' in task:
        print 'Plotting maps of beta layouts'
        scheme = 'square'
        B = [0.5, 1, 1.5, 2]
        myfig = plt.figure(figsize=(15, 6))
        maxG = np.zeros(len(B))
        gammas = np.zeros((len(B), 30))
        for i, b in enumerate(B):
            N = EU_Nodes_usage('../' + resPath + 'b_' + str(b) + '_' + scheme + '.npz')
            gammas[i] = [n.gamma for n in N]
            maxG[i] = max(gammas[i])
        maxG = max(maxG)
        point = 1 / maxG
        gammas /= maxG

        redgreendict = {'red': [(0.0, 0.8, 0.8), (point, 1.0, 1.0), (1.0, 0.0, 0.0)],
                        'green': [(0.0, 0.0, 0.0), (point, 1.0, 1.0), (1.0, 0.7, 0.7)],
                        'blue': [(0.0, 0.0, 0.0), (point, 1.0, 1.0), (1.0, 0.0, 0.0)]}
        cmap = LinearSegmentedColormap('redgreen', redgreendict, 1000)

        for j, b in enumerate(B):
                plot_europe_map(gammas[j], b, ax=plt.subplot(1, len(B), j + 1))
        cbar_ax = myfig.add_axes([0.2, 0.05, 0.6, 0.1])
        cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap, orientation='horizontal')
        cb1.solids.set_edgecolor('face')
        cb1.set_ticks([0, point, 1])
        cb1.set_ticklabels(['0', '1', str(round(maxG, 1))])
        cbar_ax.set_xlabel(r'Renewable penetration [$\gamma_n$]', fontsize=18)
        cbar_ax.xaxis.set_label_position('top')
        cbar_ax.set_xticks('none')
        plt.savefig(figPath + 'gamma_map.pdf', bbox_inches='tight')

    if 'network' in task:
        print 'Plotting network figures'
        N = EU_Nodes_usage()
        if 'day' in task:
            colors = ['solar', 'wind']
        else:
            colors = ['', 'solar', 'wind']
        for scheme in schemes:
            if scheme == 'square':
                colors.append('backup')
            for b in B:
                for direction in directions:
                    print('Direction: ' + direction)
                    for color in colors:
                        if 'total' in task:
                            print('Plotting total network figures')
                            drawnet_total(N, scheme, direction, color, b)
                        elif 'day' in task:
                            print('Plotting day/night network figures')
                            drawnet_day(N, scheme, direction, color, b)
                        else:
                            print('Plotting network figures')
                            drawnet_usage(N, scheme, direction, color, b)

    if 'levels' in task:
        levels = 4
        N = EU_Nodes_usage()
        lnames = np.array(link_namer(N))
        nnames = np.array(node_namer(N))
        schemeNames = ['localised', 'synchronised']
        print('Plotting link levels')
        for b in B:
            colors = ['', 'solar', 'wind']
            for i, scheme in enumerate(schemes):
                if scheme == 'square':
                    colors.append('backup')
                name = schemeNames[i]
                quantiles = np.load(resPath + 'quant_' + str(scheme) + '_b_' + str(b) + '.npy')
                for direction in directions:
                    for color in colors:
                        if color != '':
                            colStr = '_' + color
                        else:
                            colStr = color
                        N_usages = np.load(resPath + '/N_cont_' + scheme + '_' + direction + '_b_' + str(b) + colStr + '.npy')
                        link_level_bars(levels, N_usages, quantiles, name, direction, color, nnames, lnames, b=b)

    if 'total' in task:
        print('Plotting total network usage')
        N = np.load('./results/heterogen/b_0_linear.npz', mmap_mode='r')
        node_mean_load = N['mean']
        colors = ['', 'solar', 'wind']
        for scheme in schemes:
            if scheme == 'square':
                colors.append('backup')
            for color in colors:
                plotTotal(scheme, color)

    if 'hour' in task:
        print('Plotting hourly link level usages')
        levels = 4
        N = EU_Nodes_usage()
        lnames = np.array(link_namer(N))
        nnames = np.array(node_namer(N))
        schemeNames = ['localised', 'synchronised']
        for b in B:
            colors = ['', 'solar', 'wind']
            c = ['S', 'W']
            for i, scheme in enumerate(schemes):
                if scheme == 'square':
                    colors.append('backup')
                    c.append('B')
                name = schemeNames[i]
                quantiles = np.load(resPath + 'quant_' + str(scheme) + '_b_' + str(b) + '.npy')
                for direction in directions:
                    for j, color in enumerate(colors):
                        if color == '':
                            if direction == 'combined':
                                Usages = np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_import.npy')
                                Usages += np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_export.npy')
                                Usages /= 2
                            else:
                                Usages = np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_' + direction + '.npy')
                        else:
                            Usages = np.load('./linkcolouring/heterogen/' + scheme + '_' + direction + '_b_' + str(b) + '_' + 'usage' + c[j - 1] + '.npy')
                        link_level_hour(levels, Usages, quantiles, name, direction, color, nnames, lnames, b=b)
