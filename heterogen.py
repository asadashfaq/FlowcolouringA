#! /usr/bin/env/python
from __future__ import division
import sys
from EUgrid import *
from link_colour_less import track_flows, get_link_direction
from new_linkcolouralgorithm_less import track_link_usage_total

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
B = [0, 0.5, 1, 1.5, 2, 4, 6, 8, 10, 12]

inPath = './results/heterogen/input/'
outPath = './results/heterogen/'


def calc_usage(scheme):
    """
    Calculate powermixes and nodes' usages of links and save results to file.
    """
    for b in [1, 2, 4, 8, 16]:
        N = EU_Nodes_usage('../' + inPath + 'N_cor_' + scheme + '_b=' + str(b) + '_g=1_a=0.7.npz')
        F = np.load(inPath + 'F_cor_' + scheme + '_b=' + str(b) + '_g=1_a=0.7.npy')

        """
        N2 is a new nodes object containing individual powermixes for import and
        export in the variables N2[n].power_mix and N2[n].power_mix_ex respectively.
        """
        N2, power_mixes_total = track_flows(N, F, lapse=lapse)
        N2.save_nodes(scheme + '_pm')  # save node object including powermix

        """
        track_link_usage_total tracks each nodes usage of all links. The results
        are saved to files '..._links_ex_...' and '..._links_im_...'.
        """
        boxplot, boxplotlabel = track_link_usage_total(N2, F, mode=scheme, lapse=lapse, heterogen=b)
        return


if 'trace' in task:
    calc_usage('lin')
