from __future__ import division
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
from figutils import *

"""
Create animations from results.
"""
figPath = './figures/gif/'
resolution = 10  # number of interpolations between original results


def make_europe_graph(link_weights, node_weights, t, filename='injection/', savepath=figPath, title=None):
    G = nx.Graph()

    dcolwidth = (2 * 3.425 + 0.236)
    colwidth = (3.425)
    blue = '#134b7c'

    all_countries_ISO2 = [ISO3ISO2dict[c] for c in all_countries]

    nodelist = all_countries_ISO2
    all_links_ISO2 = []
    for l in all_links:
        l_ISO2 = copy.deepcopy(l)
        for c in all_countries:
            if c in l_ISO2:
                l_ISO2 = l_ISO2.replace(c, ISO3ISO2dict[c])
        all_links_ISO2.append(l_ISO2)

    linklist = [[link[0:2], link[-2:], all_links_ISO2.index(link)]
                for link in all_links_ISO2]

    # print linklist
    for n in nodelist:
        G.add_node(n)

    for l in linklist:
        G.add_edge(l[0], l[1], weight=link_weights[l[2]])

    cmap = LinearSegmentedColormap('redgreen', redgreendict, 1000)

    fig = plt.figure(dpi=400, figsize=(0.85 * dcolwidth, 0.85 * dcolwidth * 0.8))
    ax1 = fig.add_axes([0.05, 0.08, 0.9, 0.08])
    ax2 = fig.add_axes([-0.05, 0.15, 1.1, 0.95])

    # span = 2*np.max(np.abs(node_weights))
    # norm_node_weights = [w/span+0.5 for w in node_weights]
    # now the scale shows the actual phi normalized to unit length
    # the offset and factor of 2 is because the colormap is defined
    # from 0 to 1
    phi_length = np.sqrt(np.sum(node_weights ** 2))
    norm_node_weights = [w / (2 * phi_length) + 0.5 for w in node_weights]

    # node_colors = [cmap(w) for w in norm_node_weights]
    node_colors = [cmap(w) for w in node_weights]
    nx.draw_networkx_nodes(G, pos, node_size=400, nodelist=nodelist,
                           node_color=node_colors)
    nx.draw_networkx_labels(G, pos, font_size=10.8, font_color='k',
                            font_family='sans-serif')

    maxflow = np.max(np.abs(link_weights))
    no_arrow_limit = 0.05 * maxflow
    arrow_length = 0.05
    for l in linklist:
        if -no_arrow_limit < link_weights[l[2]] < no_arrow_limit:
            continue
        if link_weights[l[2]] >= no_arrow_limit:
            x0 = pos[l[0]][0]
            y0 = pos[l[0]][1]
            x1 = pos[l[1]][0]
            y1 = pos[l[1]][1]
        if link_weights[l[2]] <= -no_arrow_limit:
            x1 = pos[l[0]][0]
            y1 = pos[l[0]][1]
            x0 = pos[l[1]][0]
            y0 = pos[l[1]][1]
        dist = np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
        plt.arrow(x0, y0, 0.4 * (x1 - x0), 0.4 * (y1 - y0), fc=blue, ec=blue,
                  head_width=1.2 * arrow_length * (abs(link_weights[l[2]])) / maxflow,
                  head_length=arrow_length)

    nx.draw_networkx_edges(G, pos, edgelist=G.edges(), edge_color=blue)

    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap, orientation='horizontal')
    cb1.set_ticks([0, 0.5, 1])
    cb1.set_ticklabels(['-1', '0', '1'])
    # ax1.set_xlabel(r'$\Phi_n$' + ' [normalized]')
    ax1.set_xlabel(r'$P_n / \left\langleL_n\right\rangle$')
    ax1.xaxis.set_label_position('top')
    ax1.set_xticks('none')
    ax2.axis('off')
    if title != None:
        fig.suptitle(title)
    fig.savefig(savepath + filename + 'injection-' + str(t) + '.png')
    plt.close()


N = np.load('./results/square.npz', mmap_mode='r')
F = np.load('./results/square-flows.npy')
ranges = [range(24)]
for r in ranges:
    for t in r:
        node_weights = N['mismatch'][:, t] + N['balancing'][:, t]
        means = N['mean']
        node_weights = node_weights / means
        if t == 30:
            print node_weights
        # Scaling to the colorbar so [-1, 1] becomes [0, 1]
        node_weights += 1
        node_weights /= 2
        make_europe_graph(F[:, t], node_weights, t)
