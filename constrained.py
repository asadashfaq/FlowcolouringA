#! /usr/bin/env/python
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from aurespf.tools import get_q, get_quant_caps
from link_colour_less import track_flows, get_link_direction
from new_linkcolouralgorithm_less import track_link_usage_total
from link_namer import node_namer, link_namer, link_dict
from functions import binMaker, bin_prob, bin_CDF, node_contrib
from europe_plusgrid import *
from EUgrid import EU_Nodes

"""
Script to investigate the transition of a node's usage of the links as
the transmission is constrained

All ouput files are saved to ./constrained/
Figures are saved to ./figures/constrained/

Call the script using only one of the following command line arguments:
- trace:        Run flow tracing and save results.
- usage:        Calculate node's usage of links.
- plot:         Create figures of usage as function of transmission constraint.
                Figures for each link and each node.
- plot total:   Plotting one figure for both export schemes comparing all nodes'
                total usages of the network as a color mesh.
                Also make two line plots of the same data.
- plot area:    Stacked area plots for each link and each node
"""

if len(sys.argv) < 2:
    raise Exception('Not enough inputs')
else:
    task = str(sys.argv[1:])


def calc_usage(N, F, lapse, scheme, b):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """

    """
    N2 is a new nodes object containing individual powermixes for import and
    export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
    """
    N2, power_mixes_total = track_flows(N, F, admat=N.pathadmat, lapse=lapse)

    """
    track_link_usage_total tracks each nodes usage of all links. The results
    are saved to files '..._links_ex_...' and '..._links_im_...'.
    """
    boxplot, boxplotlabel = track_link_usage_total(N2, F, mode=scheme, lapse=lapse, constrained=b)
    return


def calc_contribution(b, verbose=False):
    """
    Calculate nodes' contribution for a single scheme and save results to file.
    Set verbose to True to follow the progress in the command line. Not recommended for multi processing.
    """
    for scheme in schemes:
        if verbose:
            print('scheme: ' + scheme)
        F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
        # Get 99% quantile of link flow
        quantiles = [get_q(abs(F[link, :lapse]), .99) for link in range(len(F))]

        # Do everything below for both import and export usages unless we want the combined case
        for direction in directions:
            if verbose:
                print('Direction: ' + direction)
            if direction == 'combined':
                Usages = np.load('./constrained/linkcolouring/' + scheme + '_link_mix_import_b_' + str(b) + '.npy')
                Usages += np.load('./constrained/linkcolouring/' + scheme + '_link_mix_export_b_' + str(b) + '.npy')
                Usages /= 2
            else:
                Usages = np.load('./constrained/linkcolouring/' + scheme + '_link_mix_' + str(direction) + '_b_' + str(b) + '.npy')
            print('Loaded ' + scheme + ' usages')

            # Calculate usages and save to file
            Node_contributions = np.zeros((nNodes, len(F)))  # empty array for calculated usages
            for node in range(nNodes):
                if verbose:
                    print node + 1, '/ ' + str(nNodes)
                for link in range(len(F)):
                    # Stacking and sorting data
                    F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
                    exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
                    F_matrix = np.hstack([F_vert, exp_vert])  # [flow, usage]
                    F_matrix[F_matrix[:, 0].argsort()]

                    H, bin_edges = binMaker(F_matrix, quantiles[link], lapse)
                    Node_contributions[node, link] = node_contrib(H, bin_edges)

            # save results to file
            np.save('./constrained/Node_contrib_' + scheme + '_' + direction + '_b_' + str(b) + '.npy', Node_contributions)
            print('Saved Node_contributions to ./constrained/Node_contrib_' + scheme + '_' + direction + '_b_' + str(b) + '.npy')
            if direction == 'import':
                np.save('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy', quantiles)
                print('Saved 99% quantiles to ./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')


def caller(b):
    for scheme in schemes:
        N = europe_plus_Nodes(load_filename='../ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '.npz')
        F = np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy')
        calc_usage(N, F, lapse, scheme, b)


def plotter():
    """
    Make a lot of figures for a single mode and put them in
    ./constrained/figures/scheme/
    This function is very ugly because it was written very fast without thinking much
    """
    # Country's usage of all links as function of b
    for node in nodes:
        plt.figure(figsize=(12, 6))
        usages = []
        totalUsagesLin = []
        scheme = 'lin'
        for b in np.linspace(0.05, 1.5, 30):
            # Load data and results
            F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
            quantiles = np.load('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')
            N_usages = np.load('./constrained/Node_contrib_' + scheme + '_combined_b_' + str(b) + '.npy')
            usages.append(.5 * N_usages[node] / quantiles)
            totalUsagesLin.append(np.sum(.5 * N_usages, 1) / node_mean_load)
        usages = np.array(usages).transpose()
        ax1 = plt.subplot(121)
        plt.pcolormesh(usages, cmap='Blues')
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        ax1.set_xticks(np.linspace(0, 30, 16))
        ax1.set_xticklabels(np.linspace(0, 1.5, 16))
        ax1.set_yticks(np.linspace(.5, 49.5, 50))
        ax1.set_yticklabels(link_names, ha="right", va="center", fontsize=6)

        usages = []
        totalUsagesSqr = []
        scheme = 'sqr'
        for b in np.linspace(0.05, 1.5, 30):
            # Load data and results
            F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
            quantiles = np.load('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')
            N_usages = np.load('./constrained/Node_contrib_' + scheme + '_combined_b_' + str(b) + '.npy')
            usages.append(.5 * N_usages[node] / quantiles)
            totalUsagesSqr.append(np.sum(.5 * N_usages, 1) / node_mean_load)
        usages = np.array(usages).transpose()
        ax2 = plt.subplot(122)
        plt.pcolormesh(usages, cmap='Blues')
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        ax2.set_xticks(np.linspace(0, 30, 16))
        ax2.set_xticklabels(np.linspace(0, 1.5, 16))
        ax2.set_yticks(np.linspace(.5, 49.5, 50))
        ax2.set_yticklabels(link_names, ha="right", va="center", fontsize=6)
        plt.savefig('./figures/constrained/node-usage-' + str(node) + '.pdf', bbox_inches='tight')
        plt.close()

    # Figure comparing nodes' total network usage as function of b
    totalUsagesLin = np.array(totalUsagesLin).transpose()
    totalUsagesLin = totalUsagesLin[node_mean_load.argsort()]
    totalUsagesSqr = np.array(totalUsagesSqr).transpose()
    totalUsagesSqr = totalUsagesSqr[node_mean_load.argsort()]
    names_sort = names[node_mean_load.argsort()]
    np.savez('./results/constrained_results', lin=totalUsagesLin, sqr=totalUsagesSqr, names=names_sort)
    plt.figure(figsize=(12, 6))
    ax3 = plt.subplot(121)
    ax3.set_xticks(np.linspace(0, 30, 16))
    ax3.set_xticklabels(np.linspace(0, 1.5, 16), fontsize=8)
    ax3.set_yticks(np.linspace(.5, 29.5, 30))
    ax3.set_yticklabels(names_sort, ha="right", va="center", fontsize=9)
    plt.pcolormesh(totalUsagesLin, cmap='Blues')
    plt.colorbar()
    plt.xlabel(r'$\beta$')
    ax4 = plt.subplot(122)
    ax4.set_xticks(np.linspace(0, 30, 16))
    ax4.set_xticklabels(np.linspace(0, 1.5, 16), fontsize=8)
    ax4.set_yticks(np.linspace(.5, 29.5, 30))
    ax4.set_yticklabels(names_sort, ha="right", va="center", fontsize=9)
    plt.pcolormesh(totalUsagesSqr, cmap='Blues')
    plt.colorbar()
    plt.xlabel(r'$\beta$')
    plt.savefig('./figures/constrained/total-network-usage.pdf', bbox_inches='tight')
    plt.close()

    # Figure plotting the same as just above, but viewed from a different angle
    betas = np.linspace(0.05, 1.5, 30)
    plt.figure(figsize=(12, 6))
    ax5 = plt.subplot(121)
    for n in range(len(names_sort)):
        if n in [0, 1, 2, 8, 13, 19, 21]:
            col = '#000099'
            alph = .7
        elif n in [27, 28, 29]:
            col = '#990000'
            alph = .7
        else:
            col = '#000000'
            alph = .3
        plt.plot(betas, totalUsagesLin[n], '-', color=col, lw=1.5, alpha=alph)
    ax5.set_xlabel(r'$\beta$')
    ax5.set_ylabel(r'Network usage [MW$_T$/MW$_L$]')
    plt.axis([0.05, 1.5, 0, np.ceil(np.max(totalUsagesLin))])
    plt.title('Total network usage, localised')
    ax6 = plt.subplot(122)
    for n in range(len(names_sort)):
        if n in [0, 1, 2, 8, 13, 19, 21]:
            col = '#000099'
            alph = .7
        elif n in [27, 28, 29]:
            col = '#990000'
            alph = .7
        else:
            col = '#000000'
            alph = .3
        plt.plot(betas, totalUsagesSqr[n], '-', color=col, lw=1.5, alpha=alph)
    ax6.set_xlabel(r'$\beta$')
    plt.axis([0.05, 1.5, 0, np.ceil(np.max(totalUsagesSqr))])
    plt.title('Total network usage, synchronised')
    plt.savefig('./figures/constrained/total-network-usage-lines.pdf', bbox_inches='tight')
    plt.close()

    # All countries usage of a single link as function of b
    for link in links:
        plt.figure(figsize=(12, 6))
        usages = []
        scheme = 'lin'
        for b in np.linspace(0.05, 1.5, 30):
            # Load data and results
            F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
            quantiles = np.load('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')
            N_usages = np.load('./constrained/Node_contrib_' + scheme + '_combined_b_' + str(b) + '.npy')
            usages.append(.5 * N_usages[:, link] / quantiles[link])
        usages = np.array(usages).transpose()
        ax7 = plt.subplot(121)
        ax7.set_xticks(np.linspace(0, 30, 16))
        ax7.set_xticklabels(np.linspace(0, 1.5, 16))
        ax7.set_yticks(np.linspace(.5, 29.5, 30))
        ax7.set_yticklabels(names, ha="right", va="center", fontsize=9)
        plt.pcolormesh(usages, cmap='Blues')
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        plt.title('Usage of link ' + link_names[link] + ', localised')

        usages = []
        scheme = 'sqr'
        for b in np.linspace(0.05, 1.5, 30):
            # Load data and results
            F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
            quantiles = np.load('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')
            N_usages = np.load('./constrained/Node_contrib_' + scheme + '_combined_b_' + str(b) + '.npy')
            usages.append(.5 * N_usages[:, link] / quantiles[link])
        usages = np.array(usages).transpose()
        ax8 = plt.subplot(122)
        ax8.set_xticks(np.linspace(0, 30, 16))
        ax8.set_xticklabels(np.linspace(0, 1.5, 16))
        ax8.set_yticks(np.linspace(.5, 29.5, 30))
        ax8.set_yticklabels(names, ha="right", va="center", fontsize=9)
        plt.pcolormesh(usages, cmap='Blues')
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        plt.title('Usage of link ' + link_names[link] + ', synchronised')
        plt.savefig('./figures/constrained/link-usage-' + str(link) + '.pdf', bbox_inches='tight')
        plt.close()

lapse = 8766  # 280512
schemes = ['lin', 'sqr']
directions = ['import', 'export', 'combined']

"""
Calculate powermixes and nodes' usages of links and save results to file.
"""
if 'trace' in task:
    print('Mode selected: flow tracing')
    b = np.linspace(0.05, 1.5, 30)
    p = Pool(8)
    p.map(caller, b)

"""
Calculate nodes' contributions and save results to file.
"""
if 'usage' in task:
    print('Mode selected: usage calculation')
    nNodes = 30
    b = np.linspace(0.05, 1.5, 30)
    p = Pool(8)
    p.map(calc_contribution, b)

"""
Create various plots of usage and save figures to ./figures/
"""
if (('plot' in task) and ('total' not in task) and ('area' not in task)):
    print('Plotting node- and link usages')
    N = EU_Nodes()
    nodes = range(len(N))
    links = range(50)
    link_dic = link_dict(N, links)
    names = np.array([str(N[i].label) for i in range(len(N))])
    link_names = link_namer(N, links)
    N = np.load('./ConstrainedFlowData/Europe_aHE_copper_DC_lin.npz', mmap_mode='r+')
    node_mean_load = N['mean']
    N = None
    plotter()

"""
ekstra plotting of total usages
"""
if (('plot' in task) and ('total' in task) and ('area' not in task)):
    print('Plotting total node- and link usages')
    results = np.load('./results/constrained_results.npz')
    totalUsagesLin = results['lin']
    totalUsagesSqr = results['sqr']
    names = results['names']

    # Figure comparing nodes' total network usage as function of b
    plt.figure(figsize=(12, 6))
    ax3 = plt.subplot(121)
    ax3.set_xticks(np.linspace(0, 30, 16))
    ax3.set_xticklabels(np.linspace(0, 1.5, 16), fontsize=8)
    ax3.set_yticks(np.linspace(.5, 29.5, 30))
    ax3.set_yticklabels(names, ha="right", va="center", fontsize=9)
    plt.pcolormesh(totalUsagesLin, cmap='Blues')
    plt.colorbar()
    plt.xlabel(r'$\beta$')
    ax4 = plt.subplot(122)
    ax4.set_xticks(np.linspace(0, 30, 16))
    ax4.set_xticklabels(np.linspace(0, 1.5, 16), fontsize=8)
    ax4.set_yticks(np.linspace(.5, 29.5, 30))
    ax4.set_yticklabels(names, ha="right", va="center", fontsize=9)
    plt.pcolormesh(totalUsagesSqr, cmap='Blues')
    plt.colorbar()
    plt.xlabel(r'$\beta$')
    plt.savefig('./figures/constrained/total-network-usage.pdf', bbox_inches='tight')

    # Figure plotting the same as just above, but viewed from a different angle
    betas = np.linspace(0.05, 1.5, 30)
    plt.figure(figsize=(12, 6))
    ax5 = plt.subplot(121)
    for n in range(len(names)):
        if n in [4, 7, 9, 14, 15]:
            col = '#000099'
            alph = .7
        elif n in [21, 27, 28, 29]:
            col = '#990000'
            alph = .7
        else:
            col = '#000000'
            alph = .3
        plt.plot(betas, totalUsagesLin[n], '-', color=col, lw=1.5, alpha=alph)
    ax5.set_xlabel(r'$\beta$')
    ax5.set_ylabel(r'Network usage [MW$_T$/MW$_L$]')
    plt.axis([0.05, 1.5, 0, np.ceil(np.max(totalUsagesLin))])
    plt.title('Total network usage, localised')
    ax6 = plt.subplot(122)
    for n in range(len(names)):
        if n in [4, 7, 9, 14, 15]:
            col = '#000099'
            alph = .7
        elif n in [21, 27, 28, 29]:
            col = '#990000'
            alph = .7
        else:
            col = '#000000'
            alph = .3
        plt.plot(betas, totalUsagesSqr[n], '-', color=col, lw=1.5, alpha=alph)
    ax6.set_xlabel(r'$\beta$')
    plt.axis([0.05, 1.5, 0, np.ceil(np.max(totalUsagesSqr))])
    plt.title('Total network usage, synchronised')
    plt.savefig('./figures/constrained/total-network-usage-lines.pdf', bbox_inches='tight')

    # normalisation across betas
    normedUsagesLin = totalUsagesLin / np.sum(totalUsagesLin, 0)
    normedUsagesSqr = totalUsagesSqr / np.sum(totalUsagesSqr, 0)
    plt.figure(figsize=(12, 6))
    ax5 = plt.subplot(121)
    for n in range(len(names)):
        if n in [4, 7, 9, 14, 15]:
            col = '#000099'
            alph = .7
        elif n in [21, 27, 28, 29]:
            col = '#990000'
            alph = .7
        else:
            col = '#000000'
            alph = .3
        plt.plot(betas, normedUsagesLin[n], '-', color=col, lw=1.5, alpha=alph)
    ax5.set_xlabel(r'$\beta$')
    ax5.set_ylabel(r'Network usage [MW$_T$/MW$_L$], normalised across $\beta$')
    plt.axis([0.05, 1.5, 0, 1.1 * np.max(normedUsagesLin)])
    plt.title('Total network usage, localised')
    ax6 = plt.subplot(122)
    for n in range(len(names)):
        if n in [4, 7, 9, 14, 15]:
            col = '#000099'
            alph = .7
        elif n in [21, 27, 28, 29]:
            col = '#990000'
            alph = .7
        else:
            col = '#000000'
            alph = .3
        plt.plot(betas, normedUsagesSqr[n], '-', color=col, lw=1.5, alpha=alph)
    ax6.set_xlabel(r'$\beta$')
    plt.axis([0.05, 1.5, 0, 1.1 * np.max(normedUsagesSqr)])
    plt.title('Total network usage, synchronised')
    plt.savefig('./figures/constrained/total-network-usage-beta-normed.pdf', bbox_inches='tight')


"""
Create stacked area plots for each link and each node
"""
if (('plot' in task) and ('area' in task)):
    print('Plotting area figures')
    N = EU_Nodes()
    F = abs(np.load('./ConstrainedFlowData/Europe_aHE_0.05q99_DC_lin_flows.npy'))
    link_dic = link_dict(N, F)
    nodes = range(len(N))
    links = range(len(F))
    names = np.array([str(N[i].label) for i in range(len(N))])
    link_names = np.array(link_namer(N, F))
    N = np.load('./ConstrainedFlowData/Europe_aHE_copper_DC_lin.npz', mmap_mode='r + ')
    node_mean_load = N['mean']
    N = None
    betas = np.linspace(0.05, 1.5, 30)

    for n in ['', 'norm']:
        if n == 'norm':
            saveString = 'normed/'
        else:
            saveString = ''

        # Country's usage of all links as function of b
        for node in nodes:
            plt.figure(figsize=(12, 5))
            usages = []
            scheme = 'lin'
            for b in np.linspace(0.05, 1.5, 30):
                # Load data and results
                F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
                quantiles = np.load('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')
                N_usages = np.load('./constrained/Node_contrib_' + scheme + '_combined_b_' + str(b) + '.npy')
                if n == 'norm':
                    usages.append(.5 * N_usages[node] / quantiles)
                else:
                    usages.append(.5 * N_usages[node])
            usages = np.array(usages).transpose()
            ax1 = plt.subplot(121)

            max_usages = np.max(usages, 1)
            sort_usages = usages[max_usages.argsort()]
            if n == 'norm':
                norm_usages = sort_usages / np.sum(sort_usages, 0)
            else:
                norm_usages = sort_usages
            sort_link_names = link_names[max_usages.argsort()]
            top_col = ["#cc66cc", "#cc6666", "#cc9966", "#66cc66", "#66cccc", "#6666cc"]

            for a in links:
                if a < len(links) - 6:
                    col = "#6e5160"
                else:
                    col = top_col[len(links) - a - 1]
                if a == 0:
                    ax1.fill_between(betas, 0, norm_usages[a, :], facecolor=col, alpha=.7)
                    baseline = norm_usages[a, :]
                else:
                    plot_data = [norm_usages[a, i] + baseline[i] for i in range(len(baseline))]
                    ax1.fill_between(betas, baseline, plot_data, facecolor=col, alpha=.7)
                    if a >= len(links) - 10:
                        plt.text(1.51, 0 + baseline[-1] + (plot_data[-1] - baseline[-1]) / 2,
                                 sort_link_names[a], color=col, va='center', fontsize=9)
                    baseline = [baseline[i] + norm_usages[a, i] for i in range(len(baseline))]

            ax1.set_xticks(np.linspace(0, 1.5, 16))
            ax1.set_xticklabels(np.linspace(0, 1.5, 16))
            plt.axis([0.05, 1.5, 0, max(np.sum(norm_usages, 0))])
            plt.xlabel(r'$\beta$')
            if n == 'norm':
                plt.ylabel(r'$\mathcal{K}^T_{ln}/\mathcal{K}^T_l$')
            else:
                plt.ylabel(r'$\mathcal{K}^T_{ln}$ [W]')

            usages = []
            scheme = 'sqr'
            for b in np.linspace(0.05, 1.5, 30):
                # Load data and results
                F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
                quantiles = np.load('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')
                N_usages = np.load('./constrained/Node_contrib_' + scheme + '_combined_b_' + str(b) + '.npy')
                if n == 'norm':
                    usages.append(.5 * N_usages[node] / quantiles)
                else:
                    usages.append(.5 * N_usages[node])
            usages = np.array(usages).transpose()
            ax2 = plt.subplot(122)

            max_usages = np.max(usages, 1)
            sort_usages = usages[max_usages.argsort()]
            if n == 'norm':
                norm_usages = sort_usages / np.sum(sort_usages, 0)
            else:
                norm_usages = sort_usages
            sort_link_names = link_names[max_usages.argsort()]
            top_col = ["#cc66cc", "#cc6666", "#cc9966", "#66cc66", "#66cccc", "#6666cc"]

            for a in links:
                if a < len(links) - 6:
                    col = "#6e5160"
                else:
                    col = top_col[len(links) - a - 1]
                if a == 0:
                    ax2.fill_between(betas, 0, norm_usages[a, :], facecolor=col, alpha=.7)
                    baseline = norm_usages[a, :]
                else:
                    plot_data = [norm_usages[a, i] + baseline[i] for i in range(len(baseline))]
                    ax2.fill_between(betas, baseline, plot_data, facecolor=col, alpha=.7)
                    if a >= len(links) - 10:
                        plt.text(1.51, 0 + baseline[-1] + (plot_data[-1] - baseline[-1]) / 2,
                                 sort_link_names[a], color=col, va='center', fontsize=9)
                    baseline = [baseline[i] + norm_usages[a, i] for i in range(len(baseline))]

            ax2.set_xticks(np.linspace(0, 1.5, 16))
            ax2.set_xticklabels(np.linspace(0, 1.5, 16))
            plt.axis([0.05, 1.5, 0, max(np.sum(norm_usages, 0))])
            plt.xlabel(r'$\beta$')
            plt.savefig('./figures/constrained/' + saveString + 'node-usage-area-' + str(node) + '.pdf', bbox_inches='tight')
            plt.close()

        # All countries usage of a single link as function of b
        for link in links:
            plt.figure(figsize=(12, 5))
            usages = []
            scheme = 'lin'
            for b in betas:
                # Load data and results
                F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
                quantiles = np.load('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')
                N_usages = np.load('./constrained/Node_contrib_' + scheme + '_combined_b_' + str(b) + '.npy')
                if n == 'norm':
                    usages.append(.5 * N_usages[:, link] / quantiles[link])
                else:
                    usages.append(.5 * N_usages[:, link])
            usages = np.array(usages).transpose()
            ax3 = plt.subplot(121)

            max_usages = np.max(usages, 1)
            sort_usages = usages[max_usages.argsort()]
            if n == 'norm':
                norm_usages = sort_usages / np.sum(sort_usages, 0)
            else:
                norm_usages = sort_usages
            sort_names = names[max_usages.argsort()]
            top_col = ["#cc66cc", "#cc6666", "#cc9966", "#66cc66", "#66cccc", "#6666cc"]

            for a in nodes:
                if a < len(nodes) - 6:
                    col = "#6e5160"
                else:
                    col = top_col[len(nodes) - a - 1]
                if a == 0:
                    ax3.fill_between(betas, 0, norm_usages[a, :], facecolor=col, alpha=.7)
                    baseline = norm_usages[a, :]
                else:
                    plot_data = [norm_usages[a, i] + baseline[i] for i in range(len(baseline))]
                    ax3.fill_between(betas, baseline, plot_data, facecolor=col, alpha=.7)
                    if a >= len(nodes) - 10:
                        plt.text(1.51, 0 + baseline[-1] + (plot_data[-1] - baseline[-1]) / 2, sort_names[a], color=col, va='center', fontsize=9)
                    baseline = [baseline[i] + norm_usages[a, i] for i in range(len(baseline))]

            ax3.set_xticks(np.linspace(0, 1.5, 16))
            ax3.set_xticklabels(np.linspace(0, 1.5, 16))
            plt.axis([0.05, 1.5, 0, max(np.sum(norm_usages, 0))])
            plt.xlabel(r'$\beta$')
            if n == 'norm':
                plt.ylabel(r'$\mathcal{K}^T_{ln}/\mathcal{K}^T_l$')
            else:
                plt.ylabel(r'$\mathcal{K}^T_{ln}$ [W]')

            usages = []
            scheme = 'sqr'
            for b in np.linspace(0.05, 1.5, 30):
                # Load data and results
                F = abs(np.load('./ConstrainedFlowData/Europe_aHE_' + str(b) + 'q99_DC_' + scheme + '_flows.npy'))
                quantiles = np.load('./constrained/quantiles_' + scheme + '_b_' + str(b) + '.npy')
                N_usages = np.load('./constrained/Node_contrib_' + scheme + '_combined_b_' + str(b) + '.npy')
                if n == 'norm':
                    usages.append(.5 * N_usages[:, link] / quantiles[link])
                else:
                    usages.append(.5 * N_usages[:, link])
            usages = np.array(usages).transpose()
            ax4 = plt.subplot(122)

            max_usages = np.max(usages, 1)
            sort_usages = usages[max_usages.argsort()]
            if n == 'norm':
                norm_usages = sort_usages / np.sum(sort_usages, 0)
            else:
                norm_usages = sort_usages
            sort_names = names[max_usages.argsort()]
            top_col = ["#cc66cc", "#cc6666", "#cc9966", "#66cc66", "#66cccc", "#6666cc"]

            for a in nodes:
                if a < len(nodes) - 6:
                    col = "#6e5160"
                else:
                    col = top_col[len(nodes) - a - 1]
                if a == 0:
                    ax4.fill_between(betas, 0, norm_usages[a, :], facecolor=col, alpha=.7)
                    baseline = norm_usages[a, :]
                else:
                    plot_data = [norm_usages[a, i] + baseline[i] for i in range(len(baseline))]
                    ax4.fill_between(betas, baseline, plot_data, facecolor=col, alpha=.7)
                    if a >= len(nodes) - 10:
                        plt.text(1.51, 0 + baseline[-1] + (plot_data[-1] - baseline[-1]) / 2,
                                 sort_names[a], color=col, va='center', fontsize=9)
                    baseline = [baseline[i] + norm_usages[a, i] for i in range(len(baseline))]

            ax4.set_xticks(np.linspace(0, 1.5, 16))
            ax4.set_xticklabels(np.linspace(0, 1.5, 16))
            plt.axis([0.05, 1.5, 0, max(np.sum(norm_usages, 0))])
            plt.xlabel(r'$\beta$')
            plt.savefig('./figures/constrained/' + saveString + 'link-usage-area-' + str(link) + '.pdf', bbox_inches='tight')
            plt.close()
