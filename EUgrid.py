#! /usr/bin/env python

import aurespf.solvers as au
import numpy as np
au = reload(au)
alphas = [0.674268, 0.795557, 0.716455, 0.68291, 0.755273, 0.849023, 0.701074, 0.787354, 0.739453, 0.659473, 0.641748, 0.660791, 0.690674, 0.639551, 0.68862, 0.713086, 0.662549, 0.695068, 0.716016, 0.754102, 0.817676, 0.731543, 0.646582, 0.650977, 0.696533, 0.70708, 0.706787, 0.812915, 0.810791, 0.789551]

homo_alphas=[]
for i in range(30):
	for j in range(11):
		homo_alphas=np.append(homo_alphas,j/float(10))

# Defining names of data files
# List of countries
files = ['AT.npz', 'FI.npz', 'NL.npz', 'BA.npz', 'FR.npz', 'NO.npz', 'BE.npz', 'GB.npz', 'PL.npz', 'BG.npz', 'GR.npz', 'PT.npz', 'CH.npz', 'HR.npz', 'RO.npz', 'CZ.npz', 'HU.npz', 'RS.npz', 'DE.npz', 'IE.npz', 'SE.npz', 'DK.npz', 'IT.npz', 'SI.npz', 'ES.npz', 'LU.npz', 'SK.npz', 'EE.npz', 'LV.npz', 'LT.npz']

# List of regions
regions = ['AT.npz', 'BA.npz', 'BA_off.npz', 'BE.npz', 'BE_off.npz', 'BG.npz', 'BG_off.npz', 'CH.npz', 'CZ.npz', 'DE_M.npz', 'DE_NE.npz', 'DE_NE_off.npz', 'DE_NW.npz', 'DE_NW_off.npz', 'DE_SE.npz', 'DE_SW.npz', 'DE_W.npz', 'DK_E.npz', 'DK_E_off.npz', 'DK_W.npz', 'DK_W_off.npz', 'ES_NW.npz', 'ES_NW_off.npz', 'ES_SE.npz', 'ES_SE_off.npz', 'ES_SW.npz', 'ES_SW_off.npz', 'FI_N.npz', 'FI_N_off.npz', 'FI_S.npz', 'FI_S_off.npz', 'FR_NE.npz', 'FR_NW.npz', 'FR_NW_off.npz', 'FR_SE.npz', 'FR_SE_off.npz', 'FR_SW.npz', 'FR_SW_off.npz', 'GB_E_off.npz', 'GB_N.npz', 'GB_N_off.npz', 'GB_S.npz', 'GB_W_off.npz', 'GR.npz', 'GR_off.npz', 'HR.npz', 'HR_off.npz', 'HU.npz', 'IE.npz', 'IE_N.npz', 'IE_off.npz', 'IT_N.npz', 'IT_N_off.npz', 'IT_S.npz', 'IT_S_off.npz', 'IT_Sar.npz', 'IT_Sic.npz', 'LU.npz', 'NL.npz', 'NL_off.npz', 'NO_M.npz', 'NO_N.npz', 'NO_S.npz', 'NO_S_off.npz', 'PL.npz', 'PL_off.npz', 'PT_M.npz', 'PT_M_off.npz', 'PT_N.npz', 'PT_N_off.npz', 'PT_S.npz', 'PT_S_off.npz', 'RO.npz', 'RO_off.npz', 'RS.npz', 'SE_M.npz', 'SE_M_off.npz', 'SE_N.npz', 'SE_N_off.npz', 'SE_S.npz', 'SE_S_off.npz', 'SI.npz', 'SK.npz', 'EE.npz', 'LV.npz', 'LT.npz']

# List of regions where offshore windfarms have been merged with their neighbour
mergedRegions = ['AT.npz', 'BA_merge.npz', 'BE_merge.npz', 'BG_merge.npz', 'CH.npz', 'CZ.npz', 'DE_M.npz', 'DE_NE_merge.npz', 'DE_NW_merge.npz', 'DE_SE.npz', 'DE_SW.npz', 'DE_W.npz', 'DK_E_merge.npz', 'DK_W_merge.npz', 'ES_NW_merge.npz', 'ES_SE_merge.npz', 'ES_SW_merge.npz', 'FI_N_merge.npz', 'FI_S_merge.npz', 'FR_NE.npz', 'FR_NW_merge.npz', 'FR_SE_merge.npz', 'FR_SW_merge.npz', 'GB_N_merge.npz', 'GB_S_merge.npz', 'GR_merge.npz', 'HR_merge.npz', 'HU.npz', 'IE_N.npz', 'IE_merge.npz', 'IT_N_merge.npz', 'IT_S_merge.npz', 'IT_Sar.npz', 'IT_Sic.npz', 'LU.npz', 'NL_merge.npz', 'NO_M.npz', 'NO_N.npz', 'NO_S_merge.npz', 'PL_merge.npz', 'PT_M_merge.npz', 'PT_N_merge.npz', 'PT_S_merge.npz', 'RO_merge.npz', 'RS.npz', 'SE_M_merge.npz', 'SE_N_merge.npz', 'SE_S_merge.npz', 'SI.npz', 'SK.npz', 'EE.npz', 'LV.npz', 'LT.npz']

# List of superregions made from merging regions
superRegions = ['DE.npz', 'EE.npz', 'FR.npz', 'HI.npz', 'IT.npz', 'SC.npz', 'SG.npz', 'UK.npz']

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

# Nodes object for Bo's usage calculation, alphas=.7, gammas=1
def EU_Nodes_usage(load_filename=None, full_load=False):
    return au.Nodes(admat='./settings/eadmat.txt', path='./data/', prefix = "ISET_country_", files=files, load_filename=load_filename, full_load=full_load, alphas=np.ones(30)*.7, gammas=np.ones(30))

# Nodes object containing 53 regions, alphas=.7, gammas=1
def EU_Nodes_regions(load_filename=None, full_load=False):
    return au.Nodes(admat='./settings/mergedadmat.txt', path='./data/regions/', prefix = "ISET_region_", files=mergedRegions, load_filename=load_filename, full_load=full_load, alphas=np.ones(len(mergedRegions))*.7, gammas=np.ones(len(mergedRegions)))

# Nodes object containing 8 super regions, alphas=.7, gammas=1
def EU_Nodes_superRegions(load_filename=None, full_load=False):
    return au.Nodes(admat='./settings/superregionadmat.txt', path='./data/superregions/', prefix = "", files=superRegions, load_filename=load_filename, full_load=full_load, alphas=np.ones(len(superRegions))*.7, gammas=np.ones(len(superRegions)))
	
def linfo():
    link_info = open("line_info")
    LI = []
    for l in link_info:
        LI.append([l[0:8], l[12:14], int(l[15:-1])])
    return LI
