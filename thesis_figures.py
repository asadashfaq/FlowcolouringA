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
- compare:          figures to compare export schemes

The optional parameters {diagonal, bin, usage} can be used separately or combined
"""

if len(sys.argv) < 2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

figPath = './figures/thesis/'

loadOrder = [18, 4, 7, 22, 24, 20, 8, 5, 2, 6, 1, 15, 0, 10, 14,
             9, 11, 12, 16, 21, 17, 19, 3, 26, 13, 29, 27, 23, 28, 25]

loadNames = np.array(['DE', 'FR', 'GB', 'IT', 'ES', 'SE', 'PL', 'NO', 'NL',
                      'BE', 'FI', 'CZ', 'AT', 'GR', 'RO', 'BG', 'PT', 'CH',
                      'HU', 'DK', 'RS', 'IE', 'BA', 'SK', 'HR', 'LT', 'EE',
                      'SI', 'LV', 'LU'], dtype='|S4')


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

            ax.set_xlabel(r'$|f_l(t)|/\mathcal{K}_l^T$')
            plt.axis([0, Fmax[l] / qq, 0, 1])
            if mode == 'old':
                plt.savefig(figPath + 'scatter/' + str(N[n].label) + '/' + str(l) + '-' + str(direction) + '.pdf', bbox_inches='tight')
            if mode == 'new':
                plt.savefig(figPath + 'scatter/' + str(N[n].label) + '/' + str(l) + '-' + str(direction) + '-' + mode + '.pdf', bbox_inches='tight')
            plt.close()
    return

if 'scatter' in task:
    N = EU_Nodes_usage('square.npz')
    F = np.load('./results/square-flows.npy')
    Fmax = np.max(np.abs(F), 1)

    #export_usage = np.load('./linkcolouring/old_square_copper_link_mix_export_all_alpha=same.npy')
    #scatter_plotter(N, F, Fmax, export_usage, 'export', mode='old')
    #scatter_plotter(N, F, Fmax, export_usage, 'export', mode='new')
    #export_usage = []  # frees roughly 800 MB RAM before loading import usage.

    import_usage = np.load('./linkcolouring/old_square_copper_link_mix_import_all_alpha=same.npy')
    #scatter_plotter(N, F, Fmax, import_usage, 'import', mode='old')
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
        plt.xlabel(r'$G^B_n$ [MW]')
        plt.ylabel(r'$p(B_n)$')
        plt.axis([0, bins[-1], 0, top])
        plt.savefig(figPath + '/balancing/' + str(n) + '-' + str(N[n].label) + '.pdf', bbox_inches='tight')
        plt.close()

if 'compare' in task:
    schemes = ['linear', 'square']
    bars = np.zeros((2, 30))
    nodes = np.linspace(0.5, 2 * 30 - 1.5, 30)
    nodes_shift = nodes + .5
    scatter_points = nodes + 1.35
    green = '#009900'
    blue = '#000099'

    # Imports
    for i, scheme in enumerate(schemes):
        N = EU_Nodes_usage(scheme + '.npz')
        for n in N:
            P = n.mismatch + n.balancing - n.curtailment
            P[np.where(P < 0)] = 0
            bars[i, n.id] = np.mean(abs(P)) / n.mean

    print 'Imports:', np.mean(bars[0]), np.mean(bars[1]), np.mean((bars[1] - bars[0]) * 100 / bars[0])

    plt.figure(figsize=(10, 4))
    ax = plt.subplot(111)
    plt.bar(nodes, bars[0][loadOrder], width=1, color=green, edgecolor='none')
    plt.bar(nodes_shift, bars[1][loadOrder], width=1, color=blue, edgecolor='none')
    ax.set_xticks(np.linspace(2, len(N) * 2 + 2, len(N) + 1))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10.5)
    ax.xaxis.grid(False)
    ax.xaxis.set_tick_params(width=0)
    ax.set_ylabel(r'$\left\langle \max(-P_n,0) \right\rangle / \left\langle L_n \right\rangle$')
    maxes = [max(bars[0]), max(bars[1])]
    plt.axis([0, 30 * 2 + .5, 0, 1.10 * max(maxes)])

    # Legend
    artists = [plt.Rectangle((0, 0), 0, 0, ec=green, fc=green), plt.Rectangle((0, 0), 0, 0, ec=blue, fc=blue)]
    LABS = ['Localized', 'Synchronized']
    leg = plt.legend(artists, LABS, loc='upper left', ncol=len(artists), columnspacing=0.6, borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight=1.2)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=12)
    plt.savefig(figPath + 'imports.pdf', bbox_inches='tight')
    plt.close()

    # Net exports
    for i, scheme in enumerate(schemes):
        N = EU_Nodes_usage(scheme + '.npz')
        for n in N:
            bars[i, n.id] = np.mean(n.mismatch + n.balancing - n.curtailment) / n.mean

    print 'Exports:', np.mean(bars[0]), np.mean(bars[1]), np.mean((bars[1] - bars[0]) * 100 / bars[0])

    plt.figure(figsize=(10, 4))
    ax = plt.subplot(111)
    plt.plot([0, 61], [0, 0], '-k')
    plt.bar(nodes, bars[0][loadOrder], width=1, color=green, edgecolor='none')
    plt.bar(nodes_shift, bars[1][loadOrder], width=1, color=blue, edgecolor='none')
    plt.scatter(scatter_points, np.zeros(30), s=30, marker=u's', c=blue, edgecolor=blue)
    ax.set_xticks(np.linspace(2, len(N) * 2 + 2, len(N) + 1))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10.5)
    ax.xaxis.grid(False)
    ax.xaxis.set_tick_params(width=0)
    ax.set_ylabel(r'$\left\langle P_n \right\rangle / \left\langle L_n \right\rangle$')
    plt.axis([0, 30 * 2 + .5, 1.05 * min(bars[0]), 1.05 * max(bars[0])])

    # Legend
    artists = [plt.Rectangle((0, 0), 0, 0, ec=green, fc=green), plt.Rectangle((0, 0), 0, 0, ec=blue, fc=blue)]
    LABS = ['Localized', 'Synchronized']
    leg = plt.legend(artists, LABS, loc='upper left', ncol=len(artists), columnspacing=0.6, borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight=1.2)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=12)
    plt.savefig(figPath + 'exports.pdf', bbox_inches='tight')
    plt.close()

    # Balancing capacity
    for i, scheme in enumerate(schemes):
        N = EU_Nodes_usage(scheme + '.npz')
        for n in N:
            bars[i, n.id] = get_q(n.balancing, .99) / n.mean

    print 'Balancing capacity:', np.mean(bars[0]), np.mean(bars[1]), np.mean((bars[1] - bars[0]) * 100 / bars[0])

    plt.figure(figsize=(10, 4))
    ax = plt.subplot(111)
    plt.bar(nodes, bars[0][loadOrder], width=1, color=green, edgecolor='none')
    plt.bar(nodes_shift, bars[1][loadOrder], width=1, color=blue, edgecolor='none')
    ax.set_xticks(np.linspace(2, len(N) * 2 + 2, len(N) + 1))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10.5)
    ax.xaxis.grid(False)
    ax.xaxis.set_tick_params(width=0)
    ax.set_ylabel(r'$\mathcal{K}_n^B / \left\langle L_n \right\rangle$')
    maxes = [max(bars[0]), max(bars[1])]
    plt.axis([0, 30 * 2 + .5, 0, 1.10 * max(maxes)])

    # Legend
    artists = [plt.Rectangle((0, 0), 0, 0, ec=green, fc=green), plt.Rectangle((0, 0), 0, 0, ec=blue, fc=blue)]
    LABS = ['Localized', 'Synchronized']
    leg = plt.legend(artists, LABS, loc='upper left', ncol=len(artists), columnspacing=0.6, borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight=1.2)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=12)
    plt.savefig(figPath + 'backup_capacity.pdf', bbox_inches='tight')
    plt.close()

    # Balancing energy
    for i, scheme in enumerate(schemes):
        N = EU_Nodes_usage(scheme + '.npz')
        for n in N:
            bars[i, n.id] = np.mean(n.balancing) / n.mean

    print 'Balancing energy:', np.mean(bars[0]), np.mean(bars[1]), np.mean((bars[1] - bars[0]) * 100 / bars[0])

    plt.figure(figsize=(10, 4))
    ax = plt.subplot(111)
    plt.bar(nodes, bars[0][loadOrder], width=1, color=green, edgecolor='none')
    plt.bar(nodes_shift, bars[1][loadOrder], width=1, color=blue, edgecolor='none')
    ax.set_xticks(np.linspace(2, len(N) * 2 + 2, len(N) + 1))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10.5)
    ax.xaxis.grid(False)
    ax.xaxis.set_tick_params(width=0)
    ax.set_ylabel(r'$ E_n^B / \left\langle L_n \right\rangle$')
    maxes = [max(bars[0]), max(bars[1])]
    plt.axis([0, 30 * 2 + .5, 0, 1.10 * max(maxes)])

    # Legend
    artists = [plt.Rectangle((0, 0), 0, 0, ec=green, fc=green), plt.Rectangle((0, 0), 0, 0, ec=blue, fc=blue)]
    LABS = ['Localized', 'Synchronized']
    leg = plt.legend(artists, LABS, loc='upper left', ncol=len(artists), columnspacing=0.6, borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight=1.2)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=12)
    plt.savefig(figPath + 'balancing_energy.pdf', bbox_inches='tight')
    plt.close()
