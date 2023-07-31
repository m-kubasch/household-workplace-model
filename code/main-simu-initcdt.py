#!/usr/bin/env python3

import sys
import numpy as np
import agentbasedmodel as abm
import ioput 

rng = np.random.default_rng()

### PARAMETERS

# Reading parameters from input

K = int(sys.argv[1])

betaG = float(sys.argv[2])
lambdaH = float(sys.argv[3])
lambdaW = float(sys.argv[4])
gamma = float(sys.argv[5])

tstop = float(sys.argv[6])
threshold = float(sys.argv[7])

job_id = int(sys.argv[8])
array_nb = int(sys.argv[9])


# Household and workplace size distributions

pH = np.array([0.37, 0.32, 0.13, 0.12, 0.04, 0.02])

pW = np.append([0.2192, 0.12, 0.08, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03], np.repeat(0.02, 3))
pW  = np.append(pW, np.repeat(0.015, 4))
pW = np.append(pW, np.repeat(0.001,25))
pW = np.append(pW, np.repeat(0.0001, 8))
pW = np.append(pW, [0.175])


# Where to save the output
hpwd = '../article-data/simulations/initial-condition/household-states/'
wpwd = '../article-data/simulations/initial-condition/workplace-states/'

### MAIN

GH, GW = abm.generate_population_hw(pH, pW, K, rng)
data = abm.epidemic_ghw_distribution(K, GH, GW, betaG, lambdaH, lambdaW, gamma, tstop, threshold, rng)


title = 'household-stateDist-1I0-pXinsee-bG%0.5f-bH%0.5f-bW%0.5f-nu%0.5f-N%i-Tmax%i-threshold%0.4f-job%i-rep%i' %(betaG, lambdaH, lambdaW, gamma, K, tstop, threshold, job_id, array_nb)
hfilename = hpwd + title.replace('.','p')+'.txt'

title = 'workplace-stateDist-1I0-pXinsee-bG%0.5f-bH%0.5f-bW%0.5f-nu%0.5f-N%i-Tmax%i-threshold%0.4f-job%i-rep%i' %(betaG, lambdaH, lambdaW, gamma, K, tstop, threshold, job_id, array_nb)
wfilename = wpwd + title.replace('.','p')+'.txt'

ioput.save_simulated_initcdt(data, threshold, hfilename, wfilename, K, pH, pW, betaG, lambdaH, lambdaW, gamma, tstop)