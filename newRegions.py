#! /usr/bin/env python
import numpy as np
from EUgrid import regions

def idFinder(name, regions):
    i = 0
    for j in regions:
        if j == name:
            return j
        i += 1
    raise Exception('Region not found')


admat = np.genfromtxt('./settings/regionadmat.txt')


