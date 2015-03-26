import sys
import numpy as np
from pylab import plt
from aurespf.tools import *
from EUgrid import EU_Nodes_usage, EU_Nodes_regions, EU_Nodes_superRegions
from link_namer import node_namer, link_dict
from functions import *

"""
Vector flow tracing on top of the up/down stream flow tracing algorithm.
"""

"""
Initialisation
"""
if len(sys.argv) < 2:
    raise Exception('Not enough inputs!')
else:
    task = str(sys.argv[1:])

modes = ['square', 'square']
directions = ['import', 'export', 'combined']
outPath = './results/vector/'
figPath = './figures/vector/'


def usageCalc(F, quantiles, Usages, nodes, links, name):
    """
    Calculate usages and save to file
    """
    Node_contributions = np.zeros((nodes, links))  # empty array for calculated usages
    for node in range(nodes):
        print node
        for link in range(links):
            # Stacking and sorting data
            F_vert = np.reshape(F[link, :lapse], (len(F[link, :lapse]), 1))
            exp_vert = np.reshape(Usages[link, node, :lapse], (len(Usages[link, node, :lapse]), 1))
            F_matrix = np.hstack([F_vert, exp_vert])
            F_matrix[F_matrix[:, 0].argsort()]

            H, bin_edges = binMaker(F_matrix, quantiles[link], lapse)
            Node_contributions[node, link] = node_contrib(H, bin_edges, linkID=link)

    np.save(outPath + 'Node_contrib_' + mode + '_' + direction + '_' + str(name) + '.npy', Node_contributions)
    return


if 'trace' in task:
    print('tracing')
    for mode in modes:
        print(str(mode))
        N = np.load('./results/' + mode + '_pm.npz', mmap_mode='r')
        F = abs(np.load('./results/' + mode + '-flows.npy'))
        quantiles = [get_q(abs(F[link]), .99) for link in range(len(F))]
        nodes = 30
        names = ['solar', 'wind']
        meanLoads = np.reshape(N['mean'], (nodes, 1))
        genS = N['normsolar']
        genW = N['normwind']
        genSum = genS + genW
        genSum[np.where(genSum == 0)] = 1
        if mode == 'square':
            genB = np.divide(N['balancing'], meanLoads)
            genSum += genB
            normGenB = genB / genSum
            names.append('backup')
        normGenS = genS / genSum
        normGenW = genW / genSum

        for direction in directions:
            print(str(direction))
            if direction == 'combined':
                Usages = np.load(
                    './linkcolouring/old_' + mode + '_copper_link_mix_import_all_alpha=same.npy')
                Usages2 = np.load(
                    './linkcolouring/old_' + mode + '_copper_link_mix_export_all_alpha=same.npy')
                Usages += Usages2
                Usages2 = None
                Usages /= 2
            else:
                Usages = np.load('./linkcolouring/old_' + mode + '_copper_link_mix_' + direction + '_all_alpha=same.npy')

            links, nodes, lapse = Usages.shape
            usageS = np.zeros((links, nodes, lapse))
            usageW = np.zeros((links, nodes, lapse))
            for l in xrange(links):
                usageS[l] = Usages[l] * normGenS
                usageW[l] = Usages[l] * normGenW
            if mode == 'square':
                usageB = np.zeros((links, nodes, lapse))
                for l in range(links):
                    usageB[l] = Usages[l] * normGenB
            Usages = None

            print('Solar')
            usageCalc(F, quantiles, usageS, nodes, links, 'solar')
            print('Wind')
            usageCalc(F, quantiles, usageW, nodes, links, 'wind')
            if mode == 'square':
                print('Backup')
                usageCalc(F, quantiles, usageB, nodes, links, 'backup')
