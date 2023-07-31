#!/usr/bin/env python3

import sys
import numpy as np
import agentbasedmodel as abm
import largepoplimit as limit
import ebcm
import comptime
import ioput

RNG = np.random.default_rng()



def simulation(K, pH, pW, betaG, lambdaH, lambdaW, gamma, tstop, eps, rng) :
	"""Performing one simulation of the ABM. """
	GH, GW = abm.generate_population_hw(pH, pW, K, rng)
	data = abm.epidemic_ghw_timeshots(K, GH, GW, betaG, lambdaH, lambdaW, gamma, tstop, eps, rng)
	return()

### PARAMETERS

mode = sys.argv[1]

# Three possible modes : 
#	simulation : one run of the agent-based model (generating population structure and simulating the epidemic)
#	largepop : solving the large population limit
# 	ebcm : solving the edge-based model


betaG = float(sys.argv[2])
lambdaH = float(sys.argv[3])
lambdaW = float(sys.argv[4])
gamma = float(sys.argv[5])

eps = float(sys.argv[6])
tstop = float(sys.argv[7])

if mode == 'simulation' :
	K = int(sys.argv[8])
	job_id = int(sys.argv[9])
	array_nb = int(sys.argv[10])
	
else : 
	job_id = int(sys.argv[8])
	array_nb = int(sys.argv[9])

# Household and workplace size distributions
	
pH = np.array([0.37, 0.32, 0.13, 0.12, 0.04, 0.02])

pW = np.append([0.2192, 0.12, 0.08, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03], np.repeat(0.02, 3))
pW  = np.append(pW, np.repeat(0.015, 4))
pW = np.append(pW, np.repeat(0.001,25))
pW = np.append(pW, np.repeat(0.0001, 8))
pW = np.append(pW, [0.175])



### MAIN
	
if mode == 'simulation' :
	func = simulation
	param = [K, pH, pW, betaG, lambdaH, lambdaW, gamma, tstop, eps, RNG]
	filename = 'simulation-bG%0.5f-bH%0.5f-bW%0.5f-nu%0.5f-eps%0.5f-N%i-Tmax%i-job%i-rep%i' %(betaG, lambdaH, lambdaW, gamma, eps, K, tstop, job_id, array_nb)
	
elif mode == 'largepop' :
	func = limit.largepop
	param = [betaG, lambdaH, lambdaW, gamma, eps, pH, pW, tstop]
	filename = 'largepop-bG%0.5f-bH%0.5f-bW%0.5f-nu%0.5f-eps%0.5f-Tmax%i-job%i-rep%i' %(betaG, lambdaH, lambdaW, gamma, eps, tstop, job_id, array_nb)
	
elif mode == 'ebcm' :
	func = ebcm.ebcm
	param = [betaG, lambdaH, lambdaW, gamma, eps, pH, pW, tstop]
	filename = 'ebcm-bG%0.5f-bH%0.5f-bW%0.5f-nu%0.5f-eps%0.5f-Tmax%i-job%i-rep%i' %(betaG, lambdaH, lambdaW, gamma, eps, tstop, job_id, array_nb)

filename = 'ctime-' + filename.replace('.','p') + '.txt'


reftime, ftime = comptime.timefct(func, param)

ioput.save_ctime(reftime, ftime, filename, pH, pW, betaG, lambdaH, lambdaW, gamma, eps, tstop)


	