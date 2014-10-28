#! /usr/bin/env python
import numpy as np
from functions import idFinder, removeLink, addLink

"""
This script expands the old 50+33 node region network to a 53+33 node network.
The three Baltic countries are added and some links from Luxembourg are
destroyed. A new adjacency matrix is created.
This script should only be run once. If it has already been run it will tell you
so and stop, so don't worry about breaking anything.
After this script is run you should run _mergeRegions.py_.
"""

# A list of all regions including the new ones at the end.
regions = ['AT.npz', 'BA.npz', 'BA_off.npz', 'BE.npz', 'BE_off.npz', 'BG.npz',
'BG_off.npz', 'CH.npz', 'CZ.npz', 'DE_M.npz', 'DE_NE.npz', 'DE_NE_off.npz',
'DE_NW.npz', 'DE_NW_off.npz', 'DE_SE.npz', 'DE_SW.npz', 'DE_W.npz', 'DK_E.npz',
'DK_E_off.npz', 'DK_W.npz', 'DK_W_off.npz', 'ES_NW.npz', 'ES_NW_off.npz',
'ES_SE.npz', 'ES_SE_off.npz', 'ES_SW.npz', 'ES_SW_off.npz', 'FI_N.npz',
'FI_N_off.npz', 'FI_S.npz', 'FI_S_off.npz', 'FR_NE.npz', 'FR_NW.npz',
'FR_NW_off.npz', 'FR_SE.npz', 'FR_SE_off.npz', 'FR_SW.npz', 'FR_SW_off.npz',
'GB_E_off.npz', 'GB_N.npz', 'GB_N_off.npz', 'GB_S.npz', 'GB_W_off.npz',
'GR.npz', 'GR_off.npz', 'HR.npz', 'HR_off.npz', 'HU.npz', 'IE.npz', 'IE_N.npz',
'IE_off.npz', 'IT_N.npz', 'IT_N_off.npz', 'IT_S.npz', 'IT_S_off.npz',
'IT_Sar.npz', 'IT_Sic.npz', 'LU.npz', 'NL.npz', 'NL_off.npz', 'NO_M.npz',
'NO_N.npz', 'NO_S.npz', 'NO_S_off.npz', 'PL.npz', 'PL_off.npz', 'PT_M.npz',
'PT_M_off.npz', 'PT_N.npz', 'PT_N_off.npz', 'PT_S.npz', 'PT_S_off.npz',
'RO.npz', 'RO_off.npz', 'RS.npz', 'SE_M.npz', 'SE_M_off.npz', 'SE_N.npz',
'SE_N_off.npz', 'SE_S.npz', 'SE_S_off.npz', 'SI.npz', 'SK.npz', 'EE.npz',
'LV.npz', 'LT.npz']

print('You are about to create a new and improved region and super region network.')
print('This script should only be run once - ever.')

admat = np.genfromtxt('./settings/regionadmat.txt')

if admat.shape[0] != admat.shape[1]:
    raise Exception('Adjacency matrix dimensions are broken')
if ((admat.shape[0] == 86) and (len(np.nonzero(admat)[0])/2 == 128)):
    raise Exception('The network is already up to date')

oldDim = admat.shape[0]
newRegions = ['EE', 'LV', 'LT']
newLinks = [['NL', 'NO_S'], ['EE', 'FI_S'], ['EE', 'LV'], ['LV', 'LT']]
oldLinks = [['BE', 'LU'], ['FR_NE', 'LU']]

newAdmat = np.zeros((oldDim+len(newRegions), oldDim+len(newRegions)))
newAdmat[:oldDim, :oldDim] = admat

for link in oldLinks:
    removeLink(link[0], link[1], newAdmat, regions)

for newLink in newLinks:
    newAdmat = addLink(newLink[0], newLink[1], newAdmat)

np.savetxt('./settings/regionadmat.txt',newAdmat)
print('The region network has now been updated. Run _mergeRegions.py_.')
