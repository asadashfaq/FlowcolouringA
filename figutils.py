import numpy as np

"""
Useful stuff for plotting functions
"""

# Node indices and names sorted after descending mean load
loadOrder = [18, 4, 7, 22, 24, 20, 8, 5, 2, 6, 1, 15, 0, 10, 14,
             9, 11, 12, 16, 21, 17, 19, 3, 26, 13, 29, 27, 23, 28, 25]

loadNames = np.array(['DE', 'FR', 'GB', 'IT', 'ES', 'SE', 'PL', 'NO', 'NL',
                      'BE', 'FI', 'CZ', 'AT', 'GR', 'RO', 'BG', 'PT', 'CH',
                      'HU', 'DK', 'RS', 'IE', 'BA', 'SK', 'HR', 'LT', 'EE',
                      'SI', 'LV', 'LU'], dtype='|S2')

fullNames = ['Germany', 'France', 'Great Britain', 'Italy', 'Spain', 'Sweden', 'Poland',
             'Norway', 'Netherlands', 'Belgium', 'Finland', 'Czech Republic', 'Austria',
             'Greece', 'Romania', 'Bulgaria', 'Portugal', 'Switzerland', 'Hungary',
             'Denmark', 'Serbia', 'Ireland', 'Bosnia & Herz.', 'Slovakia', 'Croatia',
             'Lithuania', 'Estonia', 'Slovenia', 'Latvia', 'Luxembourg']

all_countries = ['AUT', 'FIN', 'NLD', 'BIH', 'FRA', 'NOR', 'BEL', 'GBR', 'POL', 'BGR',
                 'GRC', 'PRT', 'CHE', 'HRV', 'ROU', 'CZE', 'HUN', 'SRB', 'DEU', 'IRL',
                 'SWE', 'DNK', 'ITA', 'SVN', 'ESP', 'LUX', 'SVK', 'EST', 'LVA', 'LTU']

# Dictionary with index of the countries in the shapefiles
shapefile_index = {'AUT': 16, 'BEL': 19, 'BGR': 23, 'BIH': 26, 'CHE': 40, 'CZE': 60,
                   'DEU': 61, 'DNK': 64, 'ESP': 71, 'EST': 72, 'FIN': 74, 'FRA': 77,
                   'GBR': 81, 'GRC': 90, 'HRV': 99, 'HUN': 101, 'IRL': 107, 'ITA': 112,
                   'LTU': 136, 'LUX': 137, 'LVA': 138, 'NLD': 168, 'NOR': 169, 'POL': 182,
                   'PRT': 185, 'ROU': 190, 'SRB': 210, 'SVK': 213, 'SVN': 214, 'SWE': 215}

# Define position of nodes for network figures
pos = {}
pos['AT'] = [0.55, 0.45]
pos['FI'] = [.95, 1.1]
pos['NL'] = [0.40, 0.85]
pos['BA'] = [0.65, 0.15]
pos['FR'] = [0.15, 0.60]
pos['NO'] = [0.5, 1.1]
pos['BE'] = [0.275, 0.775]
pos['GB'] = [0.10, 1.05]
pos['PL'] = [0.75, 0.8]
pos['BG'] = [0.9, 0.0]
pos['GR'] = [0.7, 0.0]
pos['PT'] = [0.0, 0.15]
pos['CH'] = [0.4, 0.45]
pos['HR'] = [0.75, 0.3]
pos['RO'] = [1.0, 0.15]
pos['CZ'] = [0.75, 0.60]
pos['HU'] = [1.0, 0.45]
pos['RS'] = [0.85, 0.15]
pos['DE'] = [0.45, 0.7]
pos['IE'] = [0.0, 0.95]
pos['SE'] = [0.75, 1.0]
pos['DK'] = [0.5, 0.875]
pos['IT'] = [0.4, 0.2]
pos['SI'] = [0.55, 0.3]
pos['ES'] = [0.15, 0.35]
pos['LU'] = [0.325, 0.575]
pos['SK'] = [0.90, 0.55]
pos['EE'] = [1.0, 0.985]
pos['LV'] = [0.975, 0.87]
pos['LT'] = [0.925, 0.77]


# Dictionaries for color maps
blueDict = {'red': ((0.0, 1.0, 1.0), (0.15, 0.0, 0.0), (1.0, 0.0, 0.0)),
            'green': ((0.0, 1.0, 1.0), (0.15, 0.0, 0.0), (1.0, 0.0, 0.0)),
            'blue': ((0.0, 1.0, 1.0), (0.15, 1.0, 1.0), (1.0, 0.3, 0.3))}

blueDict2 = {'red': ((0.0, 1.0, 1.0), (.15, 0, 0), (.4, 0.5, 0.5), (1.0, 1, 1)),
             'green': ((0.0, 1.0, 1.0), (.15, 0, 0), (.4, 0, 0), (1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0), (.15, 1, 1), (.4, 0.5, 0.5), (1.0, 0, 0))}

orangeDict = {'red': ((0.0, 1.0, 1.0), (0.1, 1.0, 1.0), (1.0, 0.3, 0.3)),
              'green': ((0.0, 1.0, 1.0), (0.1, 0.65, .65), (1.0, .1, .1)),
              'blue': ((0.0, 1.0, 1.0), (0.1, 0.0, 0.0), (1.0, 0.0, 0.0))}

brownDict = {'red': ((0.0, 1.0, 1.0), (0.05, 0.57, 0.57), (1.0, 0.2, 0.2)),
             'green': ((0.0, 1.0, 1.0), (0.05, 0.36, .36), (1.0, .12, .12)),
             'blue': ((0.0, 1.0, 1.0), (0.05, 0.15, 0.15), (1.0, 0.01, 0.01))}
