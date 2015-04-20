import numpy as np
from multiprocessing import Pool
import aurespf.solvers as au
from EUgrid import EU_Nodes_log
from aurespf.tools import *
from functions import *

"""
Flow tracing for logistic alphas and gammas and for gammas larger than 1.
"""

setPath = './settings/logistic/'
resPath = './results/logistic/'

modes = ['linear', 'square']
years = range(2015, 2055, 5)


def logSolver(year):
    Nodes = EU_Nodes_log(year=year)
    for mode in modes:
        N, F = au.solve(Nodes, mode=mode + ' copper verbose')
        N.save_nodes(mode + '-' + str(year), path=resPath)
        np.save(resPath + mode + '-' + str(year) + '-flows', F)


p = Pool(4)
p.map(logSolver, years)
