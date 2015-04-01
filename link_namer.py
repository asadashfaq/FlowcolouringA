#! /usr/bin/env/python
from __future__ import division
import re
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


def link_id(n1, n2, ln):
    """
    Find link id from attached nodes' labels
    """
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
            print('link_id error: ' + str(n1) + ' and ' + str(n2) + ' not connected')
            return


def link_ids(n1, n2, ln):
    """
    Similar to link_id but takes a list of labels as n2
    """
    ids = [link_id(n1, i, ln) for i in n2]
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
