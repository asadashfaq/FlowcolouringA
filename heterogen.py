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


def vCalcCont(F, quantiles, Usages, nodes, links, name, b):
    """
    Calculate usages and save to file
    """
    Node_contributions = np.zeros((nodes, links))
    for node in range(nodes):
        for link in range(links):
            # Stacking and sorting data
            F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
            exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
            F_matrix = np.hstack([F_vert, exp_vert])
            F_matrix[F_matrix[:, 0].argsort()]

            H, bin_edges = binMaker(F_matrix, quantiles[link], lapse)
            Node_contributions[node, link] = node_contrib(H, bin_edges, linkID=link)

    np.save(resPath + 'N_cont_' + scheme + '_' + direction + '_b_' + str(b) + '_' + str(name) + '.npy', Node_contributions)
    return


def vectorTrace(d):
    """
    Vector flow tracing
    """
    scheme = d[0]
    b = d[1]
    N = np.load(resPath + 'b_' + str(b) + '_' + scheme, mmap_mode='r')
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
        if not os.path.exists(resPath + scheme + '_' + direction + '_' + 'usageS.npy'):
            if direction == 'combined':
                Usages = np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_import.npy')
                Usages += np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_export.npy')
                Usages /= 2
            else:
                Usages = np.load('./linkcolouring/heterogen/' + scheme + '-b-' + str(b) + '_link_mix_' + direction + '.npy')

            links, nodes, lapse = Usages.shape
            usageS = np.zeros((links, nodes, lapse))
            usageW = np.zeros((links, nodes, lapse))
            for l in xrange(links):
                usageS[l] = Usages[l] * normGenS
                usageW[l] = Usages[l] * normGenW
            np.save(resPath + scheme + '_' + direction + '_' + 'usageS.npy', usageS)
            np.save(resPath + scheme + '_' + direction + '_' + 'usageW.npy', usageW)
            if scheme == 'square':
                usageB = np.zeros((links, nodes, lapse))
                for l in range(links):
                    usageB[l] = Usages[l] * normGenB
                np.save(resPath + scheme + '_' + direction + '_' + 'usageB.npy', usageB)
            Usages = None

            print('Solar')
            vCalcCont(F, quantiles, usageS, nodes, links, 'solar', b)
            print('Wind')
            vCalcCont(F, quantiles, usageW, nodes, links, 'wind', b)
            if scheme == 'square':
                print('Backup')
                vCalcCont(F, quantiles, usageB, nodes, links, 'backup', b)

        else:
            print('Solar')
            usage = np.load(resPath + scheme + '_' + direction + '_' + 'usageS.npy')
            links, nodes, lapse = usage.shape
            vCalcCont(F, quantiles, usage, nodes, links, 'solar', b)
            print('Wind')
            usage = np.load(resPath + scheme + '_' + direction + '_' + 'usageW.npy')
            links, nodes, lapse = usage.shape
            vCalcCont(F, quantiles, usage, nodes, links, 'wind', b)
            if scheme == 'square':
                print('Backup')
                usage = np.load(resPath + scheme + '_' + direction + '_' + 'usageB.npy')
                links, nodes, lapse = usage.shape
                vCalcCont(F, quantiles, usage, nodes, links, 'backup', b)


d = []
for i in range(2 * len(B)):
    d.append([schemes[i // len(B)], B[i % len(B)]])

if 'solve' in task:
    p = Pool(4)
    p.map(solveFlows, d)

if 'trace' in task:
    p = Pool(4)
    p.map(traceFlow, d)

if 'cont' in task:
    p = Pool(4)
    p.map(calcCont, d)

if 'vector' in task:
    p = Pool(4)
    p.map(vectorTrace, d)
