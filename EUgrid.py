#! /usr/bin/env python

import aurespf.solvers as au
import numpy as np
au = reload(au)
alphas = [0.674268, 0.795557, 0.716455, 0.68291, 0.755273, 0.849023, 0.701074, 0.787354, 0.739453, 0.659473, 0.641748, 0.660791, 0.690674, 0.639551, 0.68862, 0.713086, 0.662549, 0.695068, 0.716016, 0.754102, 0.817676, 0.731543, 0.646582, 0.650977, 0.696533, 0.70708, 0.706787, 0.812915, 0.810791, 0.789551]

homo_alphas=[]
for i in range(30):
	for j in range(11):
		homo_alphas=np.append(homo_alphas,j/float(10))

files = ['AT.npz', 'FI.npz', 'NL.npz', 'BA.npz', 'FR.npz', 'NO.npz', 'BE.npz', 'GB.npz', 'PL.npz', 'BG.npz', 'GR.npz', 'PT.npz', 'CH.npz', 'HR.npz', 'RO.npz', 'CZ.npz', 'HU.npz', 'RS.npz', 'DE.npz', 'IE.npz', 'SE.npz', 'DK.npz', 'IT.npz', 'SI.npz', 'ES.npz', 'LU.npz', 'SK.npz', 'EE.npz', 'LV.npz', 'LT.npz']

def EU_Nodes(load_filename=None, full_load=False):
    return au.Nodes(admat='./settings/eadmat.txt', path='./data/', prefix = "ISET_country_", files=files, load_filename=load_filename, full_load=full_load, alphas=alphas, gammas=np.ones(30))

def EU_Nodes_homo(load_filename=None, full_load=False,alph=0):# for homogenous alphas, use alph from 0 to 1 included
    return au.Nodes(admat='./settings/eadmat.txt', path='./data/', prefix = "ISET_country_", files=files, load_filename=load_filename, full_load=full_load, alphas=homo_alphas[alph*10::11], gammas=np.ones(30))

def EU_Nodes_gamma2(load_filename=None, full_load=False):
    return au.Nodes(admat='./settings/eadmat.txt', path='./data/', prefix = "ISET_country_", files=files, load_filename=load_filename, full_load=full_load, alphas=alphas, gammas=np.ones(30)*2)
def EU_Nodes_log(load_filename=None, full_load=False,year=0): # used to create a nodes object with the alphas and gammas calculated from Sarahs logistic fit function
    return au.Nodes(admat='./settings/eadmat.txt', path='./data/', prefix = "ISET_country_", files=files, load_filename=load_filename, full_load=full_load, alphas=np.load("./alphas/alpha"+str(year)+".npy"), gammas=np.load("./gammas/gamma"+str(year-2014)+".npy"))
def EU_Nodes_Ben(load_filename=None, full_load=False,alphas=None,gammas=None):
    return au.Nodes(admat='./settings/eadmat.txt', path='./data/', prefix = "ISET_country_", files=files, load_filename=load_filename, full_load=full_load, alphas=alphas, gammas=np.ones(30))
def linfo():
    link_info = open("line_info")
    LI = []
    for l in link_info:
        LI.append([l[0:8], l[12:14], int(l[15:-1])])
    return LI
