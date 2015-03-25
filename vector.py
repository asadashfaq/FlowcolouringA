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

modes = ['linear', 'square']
directions = ['import', 'export', 'combined']
outPath = './results/vector/'
figPath = './figures/vector/'

if 'trace' in task:
    for mode in modes:
        N = np.load('./results/' + mode + '_pm.npz', mmap_mode='r')
        loads = N['load']
        nodes = loads.shape[0]
        meanLoads = np.reshape(N['mean'], (nodes, 1))
        genS = N['normsolar']
        genW = N['normwind']
        genSum = genS + genW
        genSum[np.where(genSum == 0)] = 1
        if mode == 'square':
            genB = np.divide(N['balancing'], meanLoads)
            genSum += genB
            normGenB = genB / genSum
        normGenS = genS / genSum
        normGenW = genW / genSum

        for direction in directions:
            if direction == 'combined':
                Usages1 = np.load(
                    './linkcolouring/old_' + mode + '_copper_link_mix_import_all_alpha=same.npy')
                Usages2 = np.load(
                    './linkcolouring/old_' + mode + '_copper_link_mix_export_all_alpha=same.npy')
                Usages = (Usages1 + Usages2) * .5
                Usages1, Usages2 = None, None
            else:
                Usages = np.load('./linkcolouring/old_' + mode + '_copper_link_mix_' + direction + '_all_alpha=same.npy')        links, nodes, lapse = Usages.shape

            usageS = np.zeros((links, nodes, lapse))
            usageW = np.zeros((links, nodes, lapse))
            for l in xrange(links):
                usageS[l] = Usages[l] * normGenS
                usageW[l] = Usages[l] * normGenW
            if mode == 'square':
                usageB = np.zeros((links, nodes, lapse))
                for l in links:
                    usageB[l] = Usages[l] * normGenB
