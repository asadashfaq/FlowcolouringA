#! /usr/bin/env/python
from __future__ import division
import sys
from multiprocessing import Pool
import aurespf.solvers as au
from aurespf.tools import get_q
from EUgrid import *
from link_colour_less import track_flows, get_link_direction
from new_linkcolouralgorithm_less import track_link_usage_total
from functions import binMaker, bin_prob, bin_CDF, node_contrib

"""
Flow tracing on heterogeneous networks
"""

if len(sys.argv) < 2:
    raise Exception('Not enough inputs')
else:
    task = str(sys.argv[1:])

schemes = ['linear', 'square']
directions = ['import', 'export']
lapse = 70128
N_bins = 90
B = range(11)
meanEU = 345327.47685659607

inPath = './results/heterogen/input/'
resPath = './results/heterogen/'


def calcGamma(N, b, cfS=None, cfW=None, alpha=0.7):
    """
    Calculate gamma from beta factor
    """
    if not cfS:
        cfS = np.load('./settings/cf_solar.npy')
    if not cfW:
        cfW = np.load('./settings/cf_wind.npy')
    gammas = np.zeros(len(N))
    i = 0
    for n in N:
        gs = cfS[i] ** b * meanEU / sum([cfS[j] ** b * m.mean for j, m in enumerate(N)])
        gw = cfW[i] ** b * meanEU / sum([cfW[j] ** b * m.mean for j, m in enumerate(N)])
        n.gamma = alpha * gw + (1 - alpha) * gs
        i += 1
    return N


def solveFlows(d):
    scheme = d[0]
    b = d[1]
    N = EU_Nodes_usage()
    N = calcGamma(N, b)
    N, F = au.solve(N, mode=scheme + ' copper', lapse=lapse)
    N.save_nodes(scheme, path=resPath + 'b_' + str(b) + '_')
    np.save(resPath + scheme + '_b_' + str(b) + '_flows', F)


def traceFlow(scheme):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """
    for b in [0]:  # B:
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
        boxplot, boxplotlabel = track_link_usage_total(N2, F, mode=scheme, lapse=lapse, heterogen=b)


def calcCont(d):
    """
    Calculate nodes' contribution for a single mode and save results to file.
    """
    scheme = d[0]
    b = d[1]
    # pick one of three transmission paradigms
    N = EU_Nodes_usage('../' + resPath + 'b_' + str(b) + '_' + scheme + '.npz')
    F = abs(np.load(resPath + scheme + '_b_' + str(b) + '_flows.npy'))
    # Get 99% quantile of link flow
    quantiles = [get_q(abs(F[link, :lapse]), .99) for link in range(len(F))]

    # Do everything below for both import and export usages unless we want the combined case
    for direction in directions:
        if direction == 'combined':
            Usages = np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_import.npy')
            Usages += np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_export.npy')
            Usages /= 2
        else:
            Usages = np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_' + direction + '.npy')

        # Calculate usages and save to file
        Node_contributions = np.zeros((len(N), len(F)))  # empty array for calculated usages
        for node in range(len(N)):
            for link in range(len(F)):
                # Stacking and sorting data
                F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
                exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
                F_matrix = np.hstack([F_vert, exp_vert])  # [flow, usage]
                F_matrix[F_matrix[:, 0].argsort()]

                H, bin_edges = binMaker(F_matrix, quantiles[link], lapse, N_bins)
                Node_contributions[node, link] = node_contrib(H, bin_edges, linkID=link)

        # save results to file
        np.save(resPath + 'N_cont_' + scheme + '_' + direction + '_b_' + str(b) + '.npy', Node_contributions)
        if direction == 'import':
            np.save(resPath + 'quant_' + scheme + '_b_' + str(b) + '.npy', quantiles)


if 'solve' in task:
    d = []
    for i in range(2 * len(B)):
        d.append([schemes[i // len(B)], B[i % len(B)]])
    p = Pool(4)
    p.map(solveFlows, d)

if 'trace' in task:
    p = Pool(3)
    p.map(traceFlow, schemes)

if 'cont' in task:
    d = []
    for i in range(2 * len(B)):
        d.append([schemes[i // len(B)], B[i % len(B)]])
    p = Pool(4)
    p.map(calcCont, d)
