import sys
import numpy as np
from multiprocessing import Pool
import aurespf.solvers as au
from EUgrid import EU_Nodes_log, EU_Nodes_custom
from aurespf.tools import *
from functions import *
from link_colour_less import track_flows
from new_linkcolouralgorithm_less import track_link_usage_total

"""
Flow tracing for logistic alphas and gammas and for gammas larger than 1.

Ways to call the program:
- solve:        solve power flows
- color:        do flow tracing
- usage:        calculate usages
- solve gamma:  solve power flows for gammas larger than 1
- color gamma:  do flow tracing for gammas larger than 1
- usage gamma:  calculate usages for gammas larger than 1

Example call:
python logistic.py solve
"""

if len(sys.argv) < 2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

setPath = './settings/logistic/'
resPath = './results/logistic/'
modes = ['linear', 'square']
directions = ['import', 'export', 'combined']
years = range(2015, 2055, 5)
lapse = 70128
alphas = np.ones(30) * .7
gammas = [1.25, 1.5, 1.75, 2.]


def logSolver(year):
    """
    Solve power flows for logistic alphas and gammas for two export schemes for
    a specific year.
    """
    Nodes = EU_Nodes_log(year=year)
    for mode in modes:
        N, F = au.solve(Nodes, mode=mode + ' copper', lapse=lapse)
        N.save_nodes(mode + '-' + str(year), path=resPath)
        np.save(resPath + mode + '-' + str(year) + '-flows', F)


def gammaSolver(gamma):
    """
    Solve power flows for a given gamma larger than 1.
    """
    Nodes = EU_Nodes_custom(alphas=alphas, gammas=np.ones(30) * gamma)
    for mode in modes:
        N, F = au.solve(Nodes, mode=mode + ' copper', lapse=lapse)
        N.save_nodes(mode + '-g-' + str(gamma), path=resPath)
        np.save(resPath + mode + '-g-' + str(gamma) + '-flows', F)


def calcPowerMix(year):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """
    for mode in modes:
        N = EU_Nodes_log('logistic/' + mode + '-' + str(year) + '.npz', year=year)
        F = np.load(resPath + mode + '-' + str(year) + '-flows.npy')

        """
        N2 is a new nodes object containing individual powermixes for import and
        export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
        """
        N2, power_mixes_total = track_flows(N, F, lapse=lapse)
        N2.save_nodes(mode + '-' + str(year) + '_pm', path=resPath)  # save node object including powermix
        N = None

        """
        track_link_usage_total tracks each nodes usage of all links. The results
        are saved to files '..._links_ex_...' and '..._links_im_...'.
        """
        boxplot, boxplotlabel = track_link_usage_total(N2, F, mode=mode, alph='same', lapse=lapse, logistic=year)


def calcPowerMixGamma(gamma):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """
    for mode in modes:
        N = EU_Nodes_custom(alphas=alphas, gammas=np.ones(30) * gamma, load_filename='logistic/' + mode + '-g-' + str(gamma) + '.npz')
        F = np.load(resPath + mode + '-g-' + str(gamma) + '-flows.npy')

        """
        N2 is a new nodes object containing individual powermixes for import and
        export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
        """
        N2, power_mixes_total = track_flows(N, F, lapse=lapse)
        N2.save_nodes(mode + '-g-' + str(gamma) + '_pm', path=resPath)  # save node object including powermix
        N = None

        """
        track_link_usage_total tracks each nodes usage of all links. The results
        are saved to files '..._links_ex_...' and '..._links_im_...'.
        """
        boxplot, boxplotlabel = track_link_usage_total(N2, F, mode=mode, alph='same', lapse=lapse, gamma=gamma)


def calcUsage(year):
    """
    Calculate usages and save to file
    """
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

            np.save(resPath + 'usages/Ncontrib-' + mode + '-' + str(year) + '_' + direction + '_' + '.npy', Node_contribs)


def calcUsageGamma(gamma):
    """
    Calculate usages and save to file
    """
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

            np.save(resPath + 'usages/Ncontrib-' + mode + '-g-' + str(gamma) + '_' + direction + '_' + '.npy', Node_contribs)


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
