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
figPath = './figures/animation/'


def interpolator(data, steps, loop=False):
    """
    Load n by t data. Interpolate linearly _steps_ matrices between the t input
    steps. If _loop_ is set, interpolate from the last t to the first t, so a
    looping .gif can be made.
    """
    n, t = data.shape
    if not loop:
        idata = np.zeros((n, (t-1) * steps + 1))
    else:
        idata = np.zeros((n, t * steps))
    index = range(0, t * steps, steps)
    for k, i in enumerate(index):
        idata[:, i] = data[:, k]
        if k > 0:
            delta = (data[:, k] - data[:, k - 1]) / steps
            for step in xrange(i - steps + 1, i, 1):
                idata[:, step] = idata[:, step - 1] + delta
    if loop:
        delta = (data[:, 0] - data[:, -1]) / steps
        for step in xrange(1, steps):
            idata[:, i + step] = idata[:, i + step - 1] + delta
    return idata


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
    ax1.set_xlabel(r'$P_n / \left\langleL_n\right\rangle$')
    ax1.xaxis.set_label_position('top')
    ax1.set_xticks('none')
    ax2.axis('off')
    if title != None:
        fig.suptitle(title)
    if t < 10:
        tstr = '00' + str(t)
    elif 10 <= t < 100:
        tstr = '0' + str(t)
    else:
        tstr = str(t)
    fig.savefig(savepath + filename + 'injection-' + tstr + '.png')
    plt.close()

N = np.load('./results/square.npz', mmap_mode='r')
F = np.load('./results/square-flows.npy')
ranges = [range(24)]
steps = 5
for r in ranges:
    node_weights = N['mismatch'][:, r] + N['balancing'][:, r] - N['curtailment'][:, r]
    means = N['mean']
    for i in range(30):
        node_weights[i, :] = node_weights[i, :] / means[i]
    node_weights = interpolator(node_weights, steps, loop=True)
    # Scaling to the colorbar so [-1, 1] becomes [0, 1]
    node_weights += 1
    node_weights /= 2
    iFlow = interpolator(F[:, r], steps, loop=True)
    for t in xrange(node_weights.shape[1]):
        make_europe_graph(iFlow[:, t], node_weights[:, t], t)
