#!/usr/bin/env python3

import sys
import numpy as np
import ebcm
import ioput 


### PARAMETERS

# Reading parameters from input

betaG = float(sys.argv[1])
lambdaH = float(sys.argv[2])
lambdaW = float(sys.argv[3])
gamma = float(sys.argv[4])

eps = float(sys.argv[5])
tstop = float(sys.argv[6])


# Household and workplace size distributions

pH = np.array([0.37, 0.32, 0.13, 0.12, 0.04, 0.02])

pW = np.append([0.2192, 0.12, 0.08, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03], np.repeat(0.02, 3))
pW  = np.append(pW, np.repeat(0.015, 4))
pW = np.append(pW, np.repeat(0.001,25))
pW = np.append(pW, np.repeat(0.0001, 8))
pW = np.append(pW, [0.175])

# Where to save the output

pwd = '../article-data/ebcm/'


### MAIN

time, susceptibles, infected = ebcm.ebcm(betaG, lambdaH, lambdaW, gamma, eps, pH, pW, tstop)


filename = 'ebcm-bG%0.5f-bH%0.5f-bW%0.5f-nu%0.5f-eps%0.5f-Tmax%i' %(betaG, lambdaH, lambdaW, gamma, eps, tstop)
filename = filename.replace('.','p')+'.txt'
filename = pwd + filename

ioput.save_reduction(time, susceptibles, infected, filename, pH, pW, betaG, lambdaH, lambdaW, gamma, eps, tstop)