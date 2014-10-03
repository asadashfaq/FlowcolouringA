import numpy as np

"""
Script to merge wind time series for offshore regions (wind farms) to their
neighbouring onshore region
"""

# List of regions that have an offshore region attached
offRegions = ['BA', 'BE', 'BG', 'DE_NE', 'DE_NW', 'DK_E', 'DK_W', 'ES_NW', 'ES_SE', 'ES_SW', 'FI_N', 'FI_S', 'FR_NW', 'FR_SE', 'FR_SW', 'GB_E', 'GB_N', 'GB_W', 'GR', 'HR', 'IE', 'IT_N', 'IT_S', 'NL', 'NO_S', 'PL', 'PT_M', 'PT_N', 'PT_S', 'RO', 'SE_M', 'SE_N', 'SE_S']
prefix = 'ISET_region_'     # prefix of data file
offEnding = '_off'          # file ending of offshore regions
saveFileEnding = '_merge'   # file ending of saved merged data


# For each of the regions:
    # load onshore/offshore regions
    # merge wind data
    # save to new merged file

# for region in offRegions:
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
