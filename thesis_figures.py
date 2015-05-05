from __future__ import division
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from aurespf.tools import get_q
from EUgrid import EU_Nodes_usage
from functions import *

"""
Special figures for my thesis.

Call modes:
- flowhist:         flow histogram to explain quantiles
- balancing:        histogram of a country's balancing
- scatter:          scatter plot to motivate calculation of usage
- scatter diagonal: includes diagonal on un-normalized scatter plots
- scatter bin:      includes a sample bin with conditional usage
- scatter usage:    includes line showing conditional usage

The optional parameters {diagonal, bin, usage} can be used separately or combined
"""

if len(sys.argv) < 2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

figPath = './figures/thesis/'


def scatter_plotter(N, F, Fmax, usage, direction, mode):
    """
    Scatter plots of nodes' import/export usages of links saved to ./figures/.
    """
    nodes = [0, 3, 18, 21, 24]
    links = usage.shape[0]
    for l in range(links):
        diag = []
        diagflow = []
        for n in nodes:
            if not os.path.exists(figPath + 'scatter/' + str(N[n].label)):
                os.makedirs(figPath + 'scatter/' + str(N[n].label))
            names = ['diagonal', r'$99\%$ quantile', 'avg. usage', 'usage']
            plt.figure()
            ax = plt.subplot(111)
            linkflow = abs(F[l, :])
            qq = get_q(abs(F[l]), .99)
            if mode == 'old':
                usages = usage[l, n, :] / Fmax[l]
                if 'diagonal' in task:
                    plt.plot([0, Fmax[l] / qq], [0, 1], '-k', lw=2, alpha=.3)
                ax.set_ylabel(r'$C_{ln}$')
                txtCoor = .04
            if mode == 'new':
                usages = usage[l, n, :] / linkflow
                names = names[1:]
                ax.set_ylabel(r'$C_{ln}$')
                txtCoor = .9 * max(linkflow / qq)

            # 99% quantile
            plt.plot([1, 1], [0, 1], ':k', alpha=.8)

            # scatter
            plt.scatter(linkflow / qq, usages, c='#000099', edgecolor='none', alpha=.2)

            if 'usage' in task:
                # Plot bin avg usage
                F_vert = np.reshape(linkflow, (len(linkflow), 1))
                exp_vert = np.reshape(usages, (len(usages), 1))
                F_matrix = np.hstack([F_vert, exp_vert])
                F_matrix[F_matrix[:, 0].argsort()]
                H, bin_edges = binMaker(F_matrix, qq, lapse=70128)
                plt.plot(bin_edges / qq, H[:, 1], '-', c='#aa0000', lw=2)

            if 'bin' in task:
                # plot example bin
                sampleBin = 30
                bLeft = sampleBin / 90
                bRight = (sampleBin + 1) / 90
                avg = np.mean(usages[np.where(np.logical_and((linkflow / qq) > bLeft, (linkflow / qq) < bRight))])
                plt.plot([bLeft, bLeft], [0, 1], '-', color='#aa0000')
                plt.plot([bRight, bRight], [0, 1], '-', color='#aa0000')
                plt.plot([1.01 * bLeft, .99 * bRight], [avg, avg], '-', lw=2, color='#aa0000')

            # Plot link name and direction in figure
            label = link_label(l, N)
            plt.text(txtCoor, .95, label)
            plt.text(txtCoor, .9, direction)

            ax.set_xlabel(r'$|F_l(t)|/\mathcal{K}_l^T$')
            plt.axis([0, Fmax[l] / qq, 0, 1])
            if mode == 'old':
                plt.savefig(figPath + 'scatter/' + str(N[n].label) + '/' + str(l) + '-' + str(direction) + '.png', bbox_inches='tight')
            if mode == 'new':
                plt.savefig(figPath + 'scatter/' + str(N[n].label) + '/' + str(l) + '-' + str(direction) + '-' + mode + '.png', bbox_inches='tight')
            plt.close()
    return

if 'scatter' in task:
    N = EU_Nodes_usage('square.npz')
    F = np.load('./results/square-flows.npy')
    Fmax = np.max(np.abs(F), 1)

    export_usage = np.load('./linkcolouring/old_square_copper_link_mix_export_all_alpha=same.npy')
    scatter_plotter(N, F, Fmax, export_usage, 'export', mode='old')
    scatter_plotter(N, F, Fmax, export_usage, 'export', mode='new')
    export_usage = []  # frees roughly 800 MB RAM before loading import usage.

    import_usage = np.load('./linkcolouring/old_square_copper_link_mix_import_all_alpha=same.npy')
    scatter_plotter(N, F, Fmax, import_usage, 'import', mode='old')
    scatter_plotter(N, F, Fmax, import_usage, 'import', mode='new')
    import_usage = []

if 'flowhist' in task:
    # Flow histogram with quantiles
    F = np.load('./results/linear-flows.npy')
    f, bins = np.histogram(F[0], bins=400, normed=True)

    qs = [get_q(F[0], q) for q in [0.005, 0.01, 0.99, 0.995]]
    topBot = [0, 0.00010]

    textHeight = 0.000103

    plt.figure(figsize=(11, 8))
    plt.subplot(211)
    plt.fill_between(bins[:-1], f, color='#0000aa', edgecolor='#0000aa')
    # plt.plot([m, m], [0, 1], '-k', lw=2)
    for q in qs:
        plt.plot([q, q], topBot, '-', lw=2, color='#aa0000')
    plt.text(qs[0] - 1080, textHeight, '0.5%')
    plt.text(qs[1] + 0, textHeight, '1%')
    plt.text(qs[2] - 950, textHeight, '99%')
    plt.text(qs[3] + 0, textHeight, '99.5%')
    plt.text(bins[0] + 400, 0.9 * 0.0002, 'a)', fontsize=16)
    plt.axis([bins[0], bins[-2], 0, 0.0002])
    plt.xlabel(r'$F_l$ [MW]')
    plt.ylabel(r'$p(F_l)$')

    f, bins = np.histogram(abs(F[0]), bins=400, normed=True)
    qs = [get_q(abs(F[0]), q) for q in [0.99, 0.995]]
    topBot = [0, 0.00010]

    plt.subplot(212)
    plt.fill_between(bins[:-1], f, color='#0000aa', edgecolor='#0000aa')
    # plt.plot([m, m], [0, 1], '-k', lw=2)
    for q in qs:
        plt.plot([q, q], topBot, '-', lw=2, color='#aa0000')
    plt.plot([max(abs(F[0])), max(abs(F[0]))], topBot, '-', lw=2, color='#aa0000')
    plt.text(qs[0] - 550, textHeight, '99%')
    plt.text(qs[1] + 0, textHeight, '99.5%')
    plt.text(max(abs(F[0])) - 700, textHeight, '100%')
    plt.text(bins[0] + 200, 0.9 * 0.0002, 'b)', color='#f6f6f6', fontsize=16)
    plt.axis([bins[0], 14000, 0, 0.0002])
    plt.xlabel(r'$F_l$ [MW]')
    plt.ylabel(r'$p(F_l)$')

    plt.savefig(figPath + 'flowHist.pdf', bbox_inches='tight')

if 'balancing' in task:
    N = EU_Nodes_usage('square.npz')
    for n in range(30):  # [21]:
        if n in [13, 23, 27, 29]:
            top = 0.0008
        elif n in [25, 28]:
            top = 0.0016
        else:
            top = 0.0003
        b, bins = np.histogram(N[n].balancing, bins=100, normed=True)
        q = get_q(N[n].balancing, .99)
        plt.figure(figsize=(7, 4))
        plt.fill_between(bins[:-1], b, color='#0000aa', edgecolor='#0000aa')
        plt.plot([q, q], [0, .4 * top], '-', lw=2, color='#aa0000')
        plt.text(q, (.4 * top) + 0.000005, '99%')
        plt.xlabel(r'$B_n$ [MW]')
        plt.ylabel(r'$p(B_n)$')
        plt.axis([0, bins[-1], 0, top])
        plt.savefig(figPath + '/balancing/' + str(n) + '-' + str(N[n].label) + '.pdf', bbox_inches='tight')
        plt.close()
