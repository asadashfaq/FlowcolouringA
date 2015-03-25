#! /usr/bin/env/python
from __future__ import division
from link_colour_less import get_link_direction
import numpy as np

"""
A collection of commonly used functions.
"""


def link_label(num, N):
    """
    Translate link number into a string like: 'Country1-Country2'.
    """
    label = get_link_direction(num, N)
    return str(label[0].label) + '-' + str(label[1].label)


def getLengths(lengths):
    if lengths == 'countries':
        lengths = np.load('./settings/country_link_length.npy')
    elif lengths == 'regions':
        lengths = np.load('./settings/region_link_length.npy')
    elif lengths == 'superRegions':
        lengths = np.load('./settings/superRegion_link_length.npy')
    return lengths


def linkProportional(N, link_dic, quantiles, lengths=None):
    """
    Calculate a node's link proportional.
    For each node: find directly connected links and sum half their capacity
    """
    if lengths:  # includes modelling of link lengths
        if ((type(lengths) != list) and (not hasattr(lengths, 'all'))):  # load lengths if none are given
            lengths = getLengths(lengths)
        quantiles = quantiles * lengths
    link_proportional = np.zeros(len(N))
    for n in N:
        node_id = n.id
        node_label = str(n.label)
        for link in link_dic[node_label]:
            link_proportional[node_id] += quantiles[link] * .5
    link_proportional = np.reshape(
        link_proportional, (len(link_proportional), 1))
    return link_proportional


def binMaker(F_matrix, q, lapse, nbins=None, verbose=False):
    """
    Histogram to calculate each nodes average usage of a link. The last bin
    should be placed such that its right side is aligned with the 99% quantile
    of the flow on the link. This value is set as the capacity of the link.
    """
    if (verbose and (nbins in [10, 50, 100, 150, 200])):
        start = time.time()
    bin_max = np.ceil(q)  # last bin ends at 99% quantile
    if not nbins:
        nbins = 90  # this number is not at all arbitrary
    bin_size = bin_max / nbins
    bin_edges = np.linspace(bin_size, bin_max, nbins)
                            # value at the right side of each bin, last bin is
                            # 99% quantile

    flows = F_matrix[:, 0]
    usages = F_matrix[:, 1]
    H = np.zeros((nbins, 2))  # [number of events, mean usage]
    for b in range(nbins):
        if b == nbins - 2:
            usage = usages[
                np.where(np.logical_and(flows >= b * bin_size, flows < (b + 1) * bin_size))]
            # Take care of events beyond the given quantile. Scale them and
            # place them in the last bin
            indices = np.where(flows >= (b + 1) * bin_size)
            np.append(
                usage, usages[indices] / flows[indices] * (b + 1) * bin_size)
        else:
            usage = usages[
                np.where(np.logical_and(flows >= b * bin_size, flows < (b + 1) * bin_size))]
        events = len(usage)
        usage = sum(usage)
        if events > 0:
            H[b, 0] = events
            H[b, 1] = usage / events
        else:
            H[b, 0] = 0
            H[b, 1] = 0
    if (verbose and (nbins in [10, 50, 100, 150, 200])):
        end = time.time()
        print 'BinMaker with ' + str(nbins) + ' bins took ' + str(end - start) + ' seconds'
    return H, bin_edges


def bin_prob(bin_id, H):
    """
    Probability of occurence of given flow
    """
    return H[bin_id, 0] / sum(H[:, 0])


def bin_CDF(bin_id, H):
    """
    Cumulative distribution function to be used on the bins created in the
    binMaker.
    """
    i = 0
    P = 0
    while i <= bin_id:
        P += H[i, 0]
        i += 1
    P = P / sum(H[:, 0])
    return P


def node_contrib(H, bin_edges, linkID=None, lengths=None):
    """
    Calculate a node's contribution to a specific links capacity
    """
    flows = np.append([0], bin_edges)
    bin_size = flows[2] - flows[1]
    nbins = len(bin_edges)
    C = 0  # total contribution
    for i in range(nbins - 1):
        c1, c2 = 0, 0  # partial contributions
        if i == 0:
            c1 += (flows[i + 1] - flows[i])
        else:
            c1 += (flows[i + 1] - flows[i]) / (1 - bin_CDF(i, H))
        l = i + 1
        while l < nbins:
            c2 += bin_prob(l, H) * H[l, 1] / flows[l]
            l += 1
        C += c1 * c2
    if not lengths:
        return C
    else:  # includes modelling of link lengths
        if ((type(lengths) != list) and (not hasattr(lengths, 'all'))):  # load lengths if none are given
            lengths = getLengths(lengths)
        return C * lengths[linkID]
"""
The following three functions are for manipulating adjacency matrices.
They are especially used in newRegions.py.
"""


def idFinder(name, regions):
    i = 0
    for j in regions:
        if j == name + '.npz':
            return i
        i += 1
    raise Exception('Region not found')


def removeLink(name1, name2, admat, regions):
    id1 = idFinder(name1, regions)
    id2 = idFinder(name2, regions)
    if (admat[id1, id2] == 0) and (admat[id2, id1] == 0):
        raise Exception('Link not found')
    elif (admat[id1, id2] == 0) or (admat[id2, id1] == 0):
        raise Exception('Adjacency matrix is broken')
    else:
        admat[id1, id2] = 0
        admat[id2, id1] = 0


def addLink(name1, name2, admat, regions):
    id1 = idFinder(name1, regions)
    id2 = idFinder(name2, regions)
    admat[id1, id2] = 1
    admat[id2, id1] = 1
    return admat

"""
Functions for calculating distance correlation
See: http://jpktd.blogspot.com/2012/06/non-linear-dependence-measures-distance.html
"""


def dist(x, y):
    # 1d only
    return np.abs(x[:, None] - y)


def d_n(x):
    d = dist(x, x)
    dn = d - d.mean(0) - d.mean(1)[:, None] + d.mean()
    return dn


def dcov_all(x, y):
    dnx = d_n(x)
    dny = d_n(y)

    denom = np.product(dnx.shape)
    dc = (dnx * dny).sum() / denom
    dvx = (dnx ** 2).sum() / denom
    dvy = (dny ** 2).sum() / denom
    dr = dc / (np.sqrt(dvx) * np.sqrt(dvy))
    return dc, dr, dvx, dvy
