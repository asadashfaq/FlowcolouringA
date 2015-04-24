import sys
import numpy as np
from multiprocessing import Pool
import aurespf.solvers as au
import aurespf.plotting as auplot
from EUgrid import *
from aurespf.tools import *
from functions import *
from link_colour_less import track_flows
from new_linkcolouralgorithm_less import track_link_usage_total

'''
Flow tracing for logistic alphas and gammas and for gammas larger than 1.

Ways to call the program:
- solve:            solve power flows
- color:            do flow tracing
- usage:            calculate usages
- solve gamma:      solve power flows for gammas larger than 1
- color gamma:      do flow tracing for gammas larger than 1
- usage gamma:      calculate usages for gammas larger than 1
- plot:             plot gammas and network usage
- plot cap:         plot nodes contribution to link capacities
- plot cap gamma:   same as above but for gammas larger than one

Example call:
python logistic.py solve
'''

if len(sys.argv) < 2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

setPath = './settings/logistic/'
resPath = './results/logistic/'
figPath = './figures/logistic/'
modes = ['linear', 'square']
directions = ['import', 'export', 'combined']
years = range(2015, 2051, 1)
lapse = 70128
alphas = np.load('./settings/logistic/alphas/alpha2050.npy')
gammas = [1.25, 1.5, 1.75, 2.]

western = [2, 4, 6, 7]
northern = [1, 5, 20, 21, 27, 28, 29]
central = [0, 18, 12, 25]
balcans = [3, 9, 10, 13, 17, 23]
southern = [11, 24, 22]
eastern = [8, 14, 15, 16]


def logSolver(year):
    '''
    Solve power flows for logistic alphas and gammas for two export schemes for
    a specific year.
    '''
    Nodes = EU_Nodes_log(year=year)
    for mode in modes:
        N, F = au.solve(Nodes, mode=mode + ' copper', lapse=lapse)
        N.save_nodes(mode + '-' + str(year), path=resPath)
        np.save(resPath + mode + '-' + str(year) + '-flows', F)


def gammaSolver(gamma):
    '''
    Solve power flows for a given gamma larger than 1.
    '''
    Nodes = EU_Nodes_custom(alphas=alphas, gammas=np.ones(30) * gamma)
    for mode in modes:
        N, F = au.solve(Nodes, mode=mode + ' copper', lapse=lapse)
        N.save_nodes(mode + '-g-' + str(gamma), path=resPath)
        np.save(resPath + mode + '-g-' + str(gamma) + '-flows', F)


def calcPowerMix(year):
    '''
    Calculate powermixes and nodes' usages of links and save results to file.
    '''
    for mode in modes:
        N = EU_Nodes_log('logistic/' + mode + '-' + str(year) + '.npz', year=year)
        F = np.load(resPath + mode + '-' + str(year) + '-flows.npy')

        '''
        N2 is a new nodes object containing individual powermixes for import and
        export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
        '''
        N2, power_mixes_total = track_flows(N, F, lapse=lapse)
        N2.save_nodes(mode + '-' + str(year) + '_pm', path=resPath)  # save node object including powermix
        N = None

        '''
        track_link_usage_total tracks each nodes usage of all links. The results
        are saved to files '..._links_ex_...' and '..._links_im_...'.
        '''
        boxplot, boxplotlabel = track_link_usage_total(N2, F, mode=mode, alph='same', lapse=lapse, logistic=year)


def calcPowerMixGamma(gamma):
    '''
    Calculate powermixes and nodes' usages of links and save results to file.
    '''
    for mode in modes:
        N = EU_Nodes_custom(alphas=alphas, gammas=np.ones(30) * gamma, load_filename='logistic/' + mode + '-g-' + str(gamma) + '.npz')
        F = np.load(resPath + mode + '-g-' + str(gamma) + '-flows.npy')

        '''
        N2 is a new nodes object containing individual powermixes for import and
        export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
        '''
        N2, power_mixes_total = track_flows(N, F, lapse=lapse)
        N2.save_nodes(mode + '-g-' + str(gamma) + '_pm', path=resPath)  # save node object including powermix
        N = None

        '''
        track_link_usage_total tracks each nodes usage of all links. The results
        are saved to files '..._links_ex_...' and '..._links_im_...'.
        '''
        boxplot, boxplotlabel = track_link_usage_total(N2, F, mode=mode, alph='same', lapse=lapse, gamma=gamma)


def calcUsage(year):
    '''
    Calculate usages and save to file
    '''
    for mode in modes:
        F = abs(np.load(resPath + mode + '-' + str(year) + '-flows.npy'))
        quantiles = [get_q(abs(F[link]), .99) for link in range(len(F))]
        for direction in directions:
            if direction == 'combined':
                Usages = np.load(resPath + 'linkcolouring/' + mode + '-' + str(year) + '_link_mix_import.npy')
                Usages += np.load(resPath + 'linkcolouring/' + mode + '-' + str(year) + '_link_mix_export.npy')
                Usages /= 2
            else:
                Usages = np.load(resPath + 'linkcolouring/' + mode + '-' + str(year) + '_link_mix_' + direction + '.npy')

            links, nodes, lapse = Usages.shape
            Node_contribs = np.zeros((nodes, links))  # empty array for calculated usages
            for node in range(nodes):
                for link in range(links):
                    # Stacking and sorting data
                    F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
                    exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
                    F_matrix = np.hstack([F_vert, exp_vert])
                    F_matrix[F_matrix[:, 0].argsort()]

                    H, bin_edges = binMaker(F_matrix, quantiles[link], lapse)
                    Node_contribs[node, link] = node_contrib(H, bin_edges, linkID=link)

            np.save(resPath + 'usages/Ncontrib-' + mode + '-' + str(year) + '_' + direction + '.npy', Node_contribs)


def calcUsageGamma(gamma):
    '''
    Calculate usages and save to file
    '''
    for mode in modes:
        F = abs(np.load(resPath + mode + '-g-' + str(gamma) + '-flows.npy'))
        quantiles = [get_q(abs(F[link]), .99) for link in range(len(F))]
        for direction in directions:
            if direction == 'combined':
                Usages = np.load(resPath + 'linkcolouring/' + mode + '-g-' + str(gamma) + '_link_mix_import.npy')
                Usages += np.load(resPath + 'linkcolouring/' + mode + '-g-' + str(gamma) + '_link_mix_export.npy')
                Usages /= 2
            else:
                Usages = np.load(resPath + 'linkcolouring/' + mode + '-g-' + str(gamma) + '_link_mix_' + direction + '.npy')

            links, nodes, lapse = Usages.shape
            Node_contribs = np.zeros((nodes, links))  # empty array for calculated usages
            for node in range(nodes):
                for link in range(links):
                    # Stacking and sorting data
                    F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
                    exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
                    F_matrix = np.hstack([F_vert, exp_vert])
                    F_matrix[F_matrix[:, 0].argsort()]

                    H, bin_edges = binMaker(F_matrix, quantiles[link], lapse)
                    Node_contribs[node, link] = node_contrib(H, bin_edges, linkID=link)

            np.save(resPath + 'usages/Ncontrib-' + mode + '-g-' + str(gamma) + '_' + direction + '.npy', Node_contribs)


def transcaps(N, mode):
    for n in N:
        n.t_share = []
        n.g_share = []
    for y in years:
        T_C = sum(au.link_q(f) for f in np.load(resPath + mode + '-' + str(y) + '-flows.npy'))
        Uim = np.load(resPath + 'usages/Ncontrib-' + mode + '-' + str(y) + '_import' + '.npy')
        Uex = np.load(resPath + 'usages/Ncontrib-' + mode + '-' + str(y) + '_export' + '.npy')
        Gs = np.load(setPath + '/gammas/gamma' + str(y - 2014) + '.npy')
        totU = [Uex[m.id].sum() + Uim[m.id].sum() for m in N]
        for n in N:
            n.t_share.append(T_C * totU[n.id] / (np.sum(totU) * n.mean))
            n.g_share.append(Gs[n.id])
    return N


def multi_TC(N, mode):  # N must be N = transcaps(EU_Nodes())
    EU_T = sum(np.array(n.t_share) * n.mean for n in N) / np.sum(n.mean for n in N)
    western = [2, 4, 6, 7]
    northern = [1, 5, 20, 21, 27, 28, 29]
    central = [0, 18, 12, 25]
    balcans = [3, 9, 10, 13, 17, 23]
    southern = [11, 24, 22]
    eastern = [8, 14, 15, 16]
    ax = plt.subplot(1, 1, 1)
    BA_T = sum(np.array(n.t_share) * n.mean * (n.id in balcans) for n in N) / np.sum(n.mean * (n.id in balcans) for n in N)
    WE_T = sum(np.array(n.t_share) * n.mean * (n.id in western) for n in N) / np.sum(n.mean * (n.id in western) for n in N)
    ax = plt.subplot(1, 1, 1)
    for n in N:
        plt.plot(years, n.t_share, alpha=0.15)
    plt.plot(years, EU_T, color=auplot.blue, lw=3, label='Europe')
    plt.plot(years, BA_T, color=auplot.orange, lw=3, label='Balkans')
    plt.plot(years, WE_T, color=auplot.green, lw=3, label='Western Europe')
    plt.plot(years, N[18].t_share, color=auplot.yellow, lw=3, label='Germany')
    plt.plot(years, N[21].t_share, color=auplot.red, lw=3, label='Denmark')
    plt.plot(years, N[10].t_share, color=auplot.lightblue, lw=3, label='Greece')
    plt.ylabel('Network usage [MW$_T$/MW$_L$]')
    plt.xlabel('year')
    if mode == 'linear':
        plt.yticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
        plt.axis([2015, 2050, 0, 3.5])
    elif mode == 'square':
        plt.yticks([1, 2, 3, 4, 5])
        plt.axis([2015, 2050, 0, 5])
    plt.xticks([2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050])
    handles, labels = ax.get_legend_handles_labels()
    leg = plt.legend(handles, labels, loc=0, prop={'size': 8}, handletextpad=0.15, columnspacing=1.0)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    plt.setp(leg.get_texts(), fontsize='small')
    plt.gcf().set_size_inches([1.5 * auplot.colwidth, 1.0 * auplot.colwidth])
    plt.gcf().set_dpi(400)
    plt.tight_layout()
    plt.savefig(figPath + 'network_usage_' + mode + '_ex.png')
    plt.close()


def multi_GS(N):  # N must be N = transcaps(EU_Nodes())
    EU_T = sum(np.array(n.g_share) * n.mean for n in N) / np.sum(n.mean for n in N)
    ax = plt.subplot(1, 1, 1)
    BA_T = sum(np.array(n.g_share) * n.mean * (n.id in balcans) for n in N) / np.sum(n.mean * (n.id in balcans) for n in N)
    WE_T = sum(np.array(n.g_share) * n.mean * (n.id in western) for n in N) / np.sum(n.mean * (n.id in western) for n in N)
    ax = plt.subplot(1, 1, 1)
    for n in N:
        plt.plot(years, n.g_share, alpha=0.15)
    plt.plot(years, EU_T, color=auplot.blue, lw=3, label='Europe')
    plt.plot(years, BA_T, color=auplot.orange, lw=3, label='Balkans')
    plt.plot(years, WE_T, color=auplot.green, lw=3, label='Western Europe')
    plt.plot(years, N[18].g_share, color=auplot.yellow, lw=3, label='Germany')
    plt.plot(years, N[21].g_share, color=auplot.red, lw=3, label='Denmark')
    plt.plot(years, N[10].g_share, color=auplot.lightblue, lw=3, label='Greece')
    plt.ylabel(r'VRES penetration $\gamma$')
    plt.xlabel('year')
    plt.axis([2015, 2050, 0, 1])
    plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    plt.xticks([2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050])
    handles, labels = ax.get_legend_handles_labels()
    leg = plt.legend(handles, labels, loc=0, prop={'size': 8}, handletextpad=0.15, columnspacing=1.0)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    plt.setp(leg.get_texts(), fontsize='small')
    plt.gcf().set_size_inches([1.5 * auplot.colwidth, 1.0 * auplot.colwidth])
    plt.gcf().set_dpi(400)
    plt.tight_layout()
    plt.savefig(figPath + 'logistic_gammas.png')
    plt.close()


def capPlotter(mode, direction):
    nodes = 30
    nodeCaps = np.zeros((nodes, len(years)))
    for i, year in enumerate(years):
        usages = np.load(resPath + 'usages/Ncontrib-' + mode + '-' + str(year) + '_' + direction + '.npy')
        F = np.load(resPath + mode + '-' + str(year) + '-flows.npy')
        quantiles = [get_q(abs(F[link]), .99) for link in range(len(F))]
        NTC = sum(quantiles)
        nodeCaps[:, i] = np.sum(usages, 1) / NTC

    EU = np.sum(nodeCaps, 0) / nodes
    BA = np.sum(nodeCaps[balcans], 0)
    WE = np.sum(nodeCaps[western], 0)
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    for node in range(nodes):
        plt.plot(years, nodeCaps[node], alpha=.15)
    plt.plot(years, EU, color=auplot.blue, lw=3, label='avg. Europe')
    plt.plot(years, BA, color=auplot.orange, lw=3, label='Balkans')
    plt.plot(years, WE, color=auplot.green, lw=3, label='Western Europe')
    plt.plot(years, nodeCaps[18], color=auplot.yellow, lw=3, label='Germany')
    plt.plot(years, nodeCaps[21], color=auplot.red, lw=3, label='Denmark')
    plt.plot(years, nodeCaps[10], color=auplot.lightblue, lw=3, label='Greece')
    plt.ylabel(r'$\mathcal{K}^T_n(t) / \mathcal{K}^T(t)$')
    plt.xlabel('year')
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5])
    handles, labels = ax.get_legend_handles_labels()
    leg = plt.legend(handles, labels, loc=0, prop={'size': 8}, handletextpad=0.15, columnspacing=1.0)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    plt.setp(leg.get_texts(), fontsize='small')
    plt.gcf().set_size_inches([1.5 * auplot.colwidth, 1.0 * auplot.colwidth])
    plt.gcf().set_dpi(400)
    plt.tight_layout()
    plt.savefig(figPath + 'caps_' + mode + '_' + direction + '.png')
    plt.close()


def capPlotterGamma(mode, direction, gammas):
    nodes = 30
    nodeCaps = np.zeros((nodes, len(gammas)))
    for i, gamma in enumerate(gammas):
        usages = np.load(resPath + 'usages/Ncontrib-' + mode + '-g-' + str(gamma) + '_' + direction + '.npy')
        F = np.load(resPath + mode + '-g-' + str(gamma) + '-flows.npy')
        quantiles = [get_q(abs(F[link]), .99) for link in range(len(F))]
        NTC = sum(quantiles)
        nodeCaps[:, i] = np.sum(usages, 1) / NTC

    EU = np.sum(nodeCaps, 0) / nodes
    BA = np.sum(nodeCaps[balcans], 0)
    WE = np.sum(nodeCaps[western], 0)
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    for node in range(nodes):
        plt.plot(gammas, nodeCaps[node], alpha=.15)
    plt.plot(gammas, EU, color=auplot.blue, lw=3, label='avg. Europe')
    plt.plot(gammas, BA, color=auplot.orange, lw=3, label='Balkans')
    plt.plot(gammas, WE, color=auplot.green, lw=3, label='Western Europe')
    plt.plot(gammas, nodeCaps[18], color=auplot.yellow, lw=3, label='Germany')
    plt.plot(gammas, nodeCaps[21], color=auplot.red, lw=3, label='Denmark')
    plt.plot(gammas, nodeCaps[10], color=auplot.lightblue, lw=3, label='Greece')
    plt.ylabel(r'$\mathcal{K}^T_n(\gamma) / \mathcal{K}^T(\gamma)$')
    plt.xlabel('gamma')
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5])
    plt.xticks([1.25, 1.5, 1.75, 2.])
    plt.xlim(xmin=1.25)
    handles, labels = ax.get_legend_handles_labels()
    leg = plt.legend(handles, labels, loc=0, prop={'size': 8}, handletextpad=0.15, columnspacing=1.0)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    plt.setp(leg.get_texts(), fontsize='small')
    plt.gcf().set_size_inches([1.5 * auplot.colwidth, 1.0 * auplot.colwidth])
    plt.gcf().set_dpi(400)
    plt.tight_layout()
    plt.savefig(figPath + 'caps_g_' + mode + '_' + direction + '.png')
    plt.close()


if 'solve' in task and not 'gamma' in task:
    p = Pool(4)
    p.map(logSolver, years)

if 'color' in task and not 'gamma' in task:
    p = Pool(3)
    p.map(calcPowerMix, years)

if 'usage' in task and not 'gamma' in task:
    p = Pool(4)
    p.map(calcUsage, years)

if 'solve' in task and 'gamma' in task:
    p = Pool(4)
    p.map(gammaSolver, gammas)

if 'color' in task and 'gamma' in task:
    p = Pool(2)
    p.map(calcPowerMixGamma, gammas)

if 'usage' in task and 'gamma' in task:
    p = Pool(4)
    p.map(calcUsageGamma, gammas)


if 'plot' in task and 'cap' not in task:
    N = EU_Nodes_log(year=2015)
    for i, mode in enumerate(modes):
        N = transcaps(N, mode)
        if i == 0: multi_GS(N)
        multi_TC(N, mode)

if 'plot' in task and 'cap' in task and 'gamma' not in task:
    for mode in modes:
        for direction in directions:
            capPlotter(mode, direction)

if 'plot' in task and 'cap' in task and 'gamma' in task:
    for mode in modes:
        for direction in directions:
            capPlotterGamma(mode, direction, gammas)
