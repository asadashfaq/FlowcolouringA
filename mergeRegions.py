import numpy as np

"""
Script to merge wind time series for offshore regions (wind farms) to their
neighbouring onshore region
"""

# List of all regions
regions = ['AT', 'BA', 'BA_off', 'BE', 'BE_off', 'BG', 'BG_off', 'CH', 'CZ', 'DE_M', 'DE_NE', 'DE_NE_off', 'DE_NW', 'DE_NW_off', 'DE_SE', 'DE_SW', 'DE_W', 'DK_E', 'DK_E_off', 'DK_W', 'DK_W_off', 'ES_NW', 'ES_NW_off', 'ES_SE', 'ES_SE_off', 'ES_SW', 'ES_SW_off', 'FI_N', 'FI_N_off', 'FI_S', 'FI_S_off', 'FR_NE', 'FR_NW', 'FR_NW_off', 'FR_SE', 'FR_SE_off', 'FR_SW', 'FR_SW_off', 'GB_E_off', 'GB_N', 'GB_N_off', 'GB_S', 'GB_W_off', 'GR', 'GR_off', 'HR', 'HR_off', 'HU', 'IE', 'IE_N', 'IE_off', 'IT_N', 'IT_N_off', 'IT_S', 'IT_S_off', 'IT_Sar', 'IT_Sic', 'LU', 'NL', 'NL_off', 'NO_M', 'NO_N', 'NO_S', 'NO_S_off', 'PL', 'PL_off', 'PT_M', 'PT_M_off', 'PT_N', 'PT_N_off', 'PT_S', 'PT_S_off', 'RO', 'RO_off', 'RS', 'SE_M', 'SE_M_off', 'SE_N', 'SE_N_off', 'SE_S', 'SE_S_off', 'SI', 'SK']

# List of regions that have an offshore region attached
offRegions = ['BA', 'BE', 'BG', 'DE_NE', 'DE_NW', 'DK_E', 'DK_W', 'ES_NW', 'ES_SE', 'ES_SW', 'FI_N', 'FI_S', 'FR_NW', 'FR_SE', 'FR_SW', 'GB_E', 'GB_N', 'GB_W', 'GR', 'HR', 'IE', 'IT_N', 'IT_S', 'NL', 'NO_S', 'PL', 'PT_M', 'PT_N', 'PT_S', 'RO', 'SE_M', 'SE_N', 'SE_S']

prefix = 'ISET_region_'     # prefix of data file
offEnding = '_off'          # file ending of offshore regions
saveFileEnding = '_merge'   # file ending of saved merged data


# Merge winddata from offshore regions and save to new data file
for region in offRegions:
    if region == 'GB_W': continue
    # Load data
    if region != 'GB_E':
        onshore = np.load('./data/regions/'+prefix+region+'.npz')
        offshore = np.load('./data/regions/'+prefix+region+offEnding+'.npz')
        offshoreGW = offshore['Gw']
    else:
        onshore = np.load('./data/regions/'+prefix+'GB_S.npz')
        offshoreE = np.load('./data/regions/'+prefix+region+offEnding+'.npz')
        offshoreW = np.load('./data/regions/'+prefix+'GB_W'+offEnding+'.npz')
        offshoreGW = offshoreE['Gw']+offshoreW['Gw']
        region = 'GB_S'

    # Merge data
    Gw = onshore['Gw']+offshoreGW
    Gs = onshore['Gs']
    L = onshore['L']
    t = onshore['t']
    datalabel = onshore['datalabel']
    datetime_offset = onshore['datetime_offset']

    # Save merged data to new file
    np.savez('./data/regions/'+prefix+region+saveFileEnding+'.npz',
            Gw=Gw, Gs=Gs, L=L, t=t, datalabel=datalabel, datetime_offset=datetime_offset)


# Create a new adjacency matrix for the merged regions
admat = np.genfromtxt('./settings/regionadmat.txt')

index = 0
offIndex = 0
indices = []
for region in regions:
    if offIndex >= len(offRegions): continue
    offRegion = str(offRegions[offIndex])+'_off'
    if region == offRegion:
        offIndex += 1
        indices.append(index)
    index += 1

revIndices = indices[::-1]
for i in revIndices:
    admat = np.delete(admat,i,0)
    admat = np.delete(admat,i,1)

np.savetxt('./settings/mergedadmat.txt',admat)
