#!/usr/bin/env python3

import numpy as np
import largepoplimit as limit
import ioput 



def reduction(parameters, pwd = '') :
	""" Computes the reduced model (large population limit) for the desired parameters, and saves the output.
	
	INPUT
		parameters = (betaG, lambdaH, lambdaW, gamma, epsilon, Tmax)
			betaG : one-to-all infectious contact rate in the general population
			lambdaH : one-to-one infectious contact rate within households
			lambdaW : one-to-one infectious contact rate within workplaces
			gamma : recovery rate
			epsilon : proportion of infected at time 0
			Tmax : time up to which the solution is computed
		pwd = folder where to save the output; default = parent folder
			
	OUTPUT
		None
	"""
	
	betaG, lambdaH, lambdaW, gamma, eps, tstop = parameters
	
	# Household and workplace size distributions
	
	pH = np.array([0.37, 0.32, 0.13, 0.12, 0.04, 0.02])
	
	pW = np.append([0.2192, 0.12, 0.08, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03], np.repeat(0.02, 3))
	pW  = np.append(pW, np.repeat(0.015, 4))
	pW = np.append(pW, np.repeat(0.001,25))
	pW = np.append(pW, np.repeat(0.0001, 8))
	pW = np.append(pW, [0.175])
	
	
	### MAIN
	
	time, susceptibles, infected = limit.largepop(betaG, lambdaH, lambdaW, gamma, eps, pH, pW, tstop)
	
	
	filename = 'largepop-bG%0.5f-bH%0.5f-bW%0.5f-nu%0.5f-eps%0.5f-Tmax%i' %(betaG, lambdaH, lambdaW, gamma, eps, tstop)
	filename = filename.replace('.','p')+'.txt'
	filename = pwd + filename
	
	ioput.save_reduction(time, susceptibles, infected, filename, pH, pW, betaG, lambdaH, lambdaW, gamma, eps, tstop)
	
	return()



### MAIN 

parameters = [(0.125, 1.5, 0.00115, 0.125, 0.001, 60), 
	(0.125, 1.5, 0.00115, 0.125, 0.01, 60), 
	(0.125, 1.5, 0.00115, 0.125, 0.05, 60),
	(0.03, 0.05, 0.0015, 0.125, 0.001, 400),
	(0.03, 0.05, 0.0015, 0.125, 0.01, 250),
	(0.085, 0.1, 0.001, 0.125, 0.01, 120),
	(0.03, 0.05, 0.0015, 0.125, 0.005, 150),
	(0.035, 0.07, 0.0016, 0.125, 0.005, 150),
	(0.045, 0.09, 0.0018, 0.125, 0.005, 150),
	(0.05, 0.15, 0.002, 0.125, 0.005, 150),
	(0.06, 0.2, 0.0022, 0.125, 0.005, 150),
	(0.06, 0.06, 0.00075, 0.125, 0.005, 150),
	(0.07, 0.07, 0.0008, 0.125, 0.005, 150),
	(0.085, 0.1, 0.001, 0.125, 0.005, 150),
	(0.1, 0.15, 0.0011, 0.125, 0.005, 150),
	(0.125, 1.5, 0.00115, 0.125, 0.005, 150)]

for param in parameters :
	reduction(param, '../article-data/limit/')
	