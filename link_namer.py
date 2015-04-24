#! /usr/bin/env/python
from __future__ import division
import re
import os
import numpy as np
import aurespf.solvers as au
from EUgrid import *

"""
Various functions regarding names of nodes and links.
"""


def link_direction(num, N):
    """
    Adapted from Anders' version. Compatible with names of variable length
    """
    a, b, c, d, e = au.AtoKh(N)
    link_label = re.split('\W+', e[num][0])
    for k in N:
        if str(k.label) == link_label[0]:
            start_node = k
        if str(k.label) == link_label[2]:
            end_node = k
    return [start_node, end_node]


def link_label(num, N):
    """
    Translate link number into a string like: 'Country1-Country2'.
    """
    label = link_direction(num, N)
    return str(label[0].label) + '-' + str(label[1].label)


def link_namer(N=None, F=None):
    """
    Return a list of link labels: 'from-to'
    """
    if os.path.exists('lnames.npy'):
        return np.load('lnames.npy').tolist()
    else:
        if not N:
            N = EU_Nodes_usage('linear.npz')
        if not F:
            F = np.load('./results/linear-flows.npy')
        lnames = []
        for link in range(len(F)):
            lnames.append(link_label(link, N))
        return lnames


def node_namer(N=None):
    """
    Return a list of node labels eg. 'DK'
    """
    if not N:
        N = EU_Nodes_usage('linear.npz')
    nnames = []
    for node in N:
        nnames.append(node.label)
    return nnames


def link_id(n1, n2, ln, silent=False):
    """
    Find link id from attached nodes' labels
    """
    ln = np.array(ln)
    name1 = str(n1) + '-' + str(n2)
    name2 = str(n2) + '-' + str(n1)
    try:
        ID = np.where(ln == name1)[0][0]
        return ID
    except:
        try:
            ID = np.where(ln == name2)[0][0]
            return ID
        except:
            if not silent:
                print('link_id error: ' + str(n1) + ' and ' + str(n2) + ' not connected')
            return


def link_ids(n1, n2, ln, silent=False):
    """
    Similar to link_id but takes a list of labels as n2
    """
    ids = [link_id(n1, i, ln, silent) for i in n2]
    ids = np.array(ids)
    return ids[np.where(ids)]


def link_dict(N=None, F=None):
    """
    Return a dictionary of node labels and IDs of their directly connected links
    eg. 'AT':[0,1,2,3,4,5]
    """
    if not N:
        N = EU_Nodes_usage('linear.npz')
    if not F:
        F = np.load('./results/linear-flows.npy')
    ldict = {}
    for link in range(len(F)):
        label = link_direction(link, N)
        for l in label:
            if str(l.label) in ldict:
                ldict[str(l.label)].append(link)
            else:
                ldict.update({str(l.label): [link]})
    return ldict


def get_neighbors(node, admat):
    """
    Function to find a node's directly connected neighbors. Input shold be an
    integer and an adjacency matrix as an array.
    """
    temp = np.where(admat[node] != 0)[0]
    return temp[np.where(temp != node)]


def get_multi_neighbors(nodes, admat):
    """
    Similar to get_neighbors but supports multiple input nodes and returns
    unique neighbors to these.
    """
    temp = [get_neighbors(i, admat) for i in nodes]
    temp = np.concatenate(temp)
    return np.unique(temp)


def neighbor_levels(node, levels, admat):
    """
    Get a list of neighbors at different levels. Level 0 is the node it self,
    level 1 is the first neighbors and so on.
    """
    neighbors = [[node]]
    for lvl in range(levels):
        temp = get_multi_neighbors(neighbors[lvl], admat)
        ignoreList = np.concatenate(neighbors)
        ignore = [np.where(temp == i)[0].tolist() for i in ignoreList]
        temp2 = np.delete(temp, np.concatenate(ignore)).tolist()
        if len(temp2) == 0:
            return neighbors
        else:
            neighbors.append(temp2)
    return neighbors


def link_level(nl, lvl, nnames, lnames):
    """
    Return a list of ids for the links that connect a node's neighbors at level
    lvl - 1 with those at level lvl. Takes neighbor_levels as first input.
    Examples:
    lvl=0 gives the links from a node to it's first neighbors.
    lvl=1 gives the between a node's first and secondary neighbors.
    """
    ll = [link_ids(nnames[i], nnames[nl[lvl + 1]], lnames, silent=True) for i in nl[lvl]]
    return np.concatenate(ll)


def get_diameter(admat=None):
    """
    Find diameter of network
    """
    if not admat:
        admat = np.genfromtxt('./settings/eadmat.txt')
    nodes = admat.shape[0]
    d = []
    for n in range(nodes):
        d.append(len(neighbor_levels(n, 20, admat)))
    return min(d) - 1, max(d) - 1
