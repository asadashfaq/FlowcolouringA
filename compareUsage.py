import sys
import numpy as np
from pylab import plt
from multiprocessing import Pool
from aurespf.tools import *
from EUgrid import EU_Nodes_usage
from link_namer import *
from functions import *

if len(sys.argv) > 1:
    task = str(sys.argv[1:])
else:
    raise Exception('Not enough inputs!')

schemes = ['linear', 'square']  # 'RND'
directions = ['import', 'export', 'combined']
figPath = './figures/compareUsage/'

# Node indices and names sorted after descending mean load
loadOrder = [18, 4, 7, 22, 24, 20, 8, 5, 2, 6, 1, 15, 0, 10, 14,
             9, 11, 12, 16, 21, 17, 19, 3, 26, 13, 29, 27, 23, 28, 25]

loadNames = np.array(['DE', 'FR', 'GB', 'IT', 'ES', 'SE', 'PL', 'NO', 'NL',
                      'BE', 'FI', 'CZ', 'AT', 'GR', 'RO', 'BG', 'PT', 'CH',
                      'HU', 'DK', 'RS', 'IE', 'BA', 'SK', 'HR', 'LT', 'EE',
                      'SI', 'LV', 'LU'], dtype='|S4')


def bars(scheme, verbose=None, norm='load'):
    """
    Figure to compare link proportional and usage proportional for a single
    scheme and put them in ./sensitivity/figures/scheme/
    """
    # Load data and results
    F = abs(np.load('./results/' + scheme + '-flows.npy'))
    quantiles = np.load('./results/quantiles_' + scheme + '_' + str(lapse) + '.npy')
    nNodes = 30

    names = node_namer(N)  # array of node labels
    links = range(len(F))
    nodes = np.linspace(0.5, 2 * nNodes - 1.5, nNodes)
    nodes_shift = nodes + .5

    for direction in directions:
        N_usages = np.load('./results/Node_contrib_' + scheme + '_' + direction + '_' + str(lapse) + '.npy')

        # Compare node transmission to mean load
        if verbose:
            print('Plotting node comparison - ' + scheme + ' - ' + direction)
        # sort node names for x-axis
        Total_usage = np.sum(N_usages, 1)
        node_ids = np.array(range(len(N))).reshape((len(N), 1))
        node_mean_load = [n.mean for n in N]

        # Vector for normalisation
        if norm == 'cap':
            normVec = np.ones(nNodes) * sum(quantiles)
        else:
            normVec = node_mean_load

        # Calculate node proportional
        EU_load = np.sum(node_mean_load)
        Total_caps = sum(quantiles)
        Node_proportional = node_mean_load / EU_load * Total_caps / normVec
        Node_proportional = np.reshape(Node_proportional, (len(Node_proportional), 1))

        # Calculate link proportional
        link_proportional = linkProportional(N, link_dic, quantiles)
        link_proportional = [link_proportional[i] / normVec[i] for i in range(nNodes)]

        # Calculate old usage proportional
        if direction == 'combined':
            old_usages = np.load('./linkcolouring/old_' + scheme + '_copper_link_mix_import_all_alpha=same.npy')
            old_usages += np.load('./linkcolouring/old_' + scheme + '_copper_link_mix_export_all_alpha=same.npy')
        else:
            old_usages = np.load('./linkcolouring/old_' + scheme + '_copper_link_mix_' + direction + '_all_alpha=same.npy')
        avg_node_usage = np.sum(np.sum(old_usages, axis=2), axis=0) / 70128.
        avg_EU_usage = np.sum(np.sum(np.sum(old_usages, axis=2), axis=0)) / 70128.
        avg_node_usage /= avg_EU_usage
        avg_node_usage /= normVec
        avg_node_usage *= 500000

        # Calculate usage and sort countries by mean load
        normed_usage = Total_usage / normVec
        normed_usage = np.reshape(normed_usage, (len(normed_usage), 1))
        node_mean_load = np.reshape(node_mean_load, (len(node_mean_load), 1))
        data = np.hstack([normed_usage, node_ids, node_mean_load, link_proportional, Node_proportional])
        data_sort = data[data[:, 2].argsort()]
        names_sort = [names[int(i)] for i in data_sort[:, 1]]
        # flip order so largest is first
        names_sort = names_sort[::-1]
        link_proportional = data_sort[:, 3][::-1]
        Node_proportional = data_sort[:, 4][::-1]
        data_sort = data_sort[:, 0][::-1]

        plt.figure(figsize=(10, 4), facecolor='w', edgecolor='k')
        ax = plt.subplot(111)
        green = '#009900'
        blue = '#000099'

        # Plot node proportional
        plt.rc('lines', lw=2)
        plt.rc('lines', dash_capstyle='round')
        plt.plot(np.linspace(0, len(N) * 2 + 2, len(N)), Node_proportional, '--k')
        # Plot link proportional
        #plt.bar(nodes, link_proportional, width=1, color=green, edgecolor='none')
        # Plot old usage proportional
        plt.bar(nodes, avg_node_usage[loadOrder], width=1, color=green, edgecolor='none')
        # Plot usage proportional
        plt.bar(nodes_shift, data_sort, width=1, color=blue, edgecolor='none')

        # Magic with ticks and labels
        ax.set_xticks(np.linspace(2, len(N) * 2 + 2, len(N) + 1))
        ax.set_xticklabels(names_sort, rotation=60, ha="right", va="top", fontsize=10.5)

        ax.xaxis.grid(False)
        ax.xaxis.set_tick_params(width=0)
        if norm == 'cap':
            ax.set_ylabel(r'$M_n/ \mathcal{K}^T$')
        else:
            # ax.set_ylabel(r'Network usage [MW$_T$/MW$_L$]')
            ax.set_ylabel(r'$M_n/\left\langle L_n \right\rangle$')
        maxes = [max(avg_node_usage), max(data_sort)]
        plt.axis([0, nNodes * 2 + .5, 0, 1.15 * max(maxes)])

        # Legend
        artists = [plt.Line2D([0, 0], [0, 0], ls='dashed', lw=2.0, c='k'), plt.Rectangle((0, 0), 0, 0, ec=green, fc=green), plt.Rectangle((0, 0), 0, 0, ec=blue, fc=blue)]
        LABS = ['$M^1$', '$M^{3}_{old}$', '$M^{3}_{new}$']
        leg = plt.legend(artists, LABS, loc='upper left', ncol=len(artists), columnspacing=0.6, borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight=1.2)
        leg.get_frame().set_alpha(0)
        leg.get_frame().set_edgecolor('white')
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=12)    # the legend text fontsize

        plt.savefig(figPath + scheme + '/network-usage-' + direction + '-' + norm + '.png', bbox_inches='tight')
        if verbose:
            print('Saved figures to ./figures/compareUsage/' + scheme + '/network-usage-' + direction + '-' + norm + '.png')


print('Plotting total network usage')
lapse = 70128
N = EU_Nodes_usage()
print('Building link dictionary')
link_dic = link_dict(N)  # dictionary of links directly connected to each node
if 'cap' in task:
    print('plotting normalised to network capacity')
    for scheme in schemes:
        bars(scheme, norm='cap')
else:
    print('Plotting')
    p = Pool(len(schemes))
    p.map(bars, schemes)
