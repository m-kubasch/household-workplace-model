#!/usr/bin/env python3
#!/usr/bin/env python3

import numpy as np
import os 

import largepoplimit as limit
import ioput


### PARAMETERS 

nbrep = 5000 # Maximal number of simulations for which an initial condition is available

# hbase is such that all files with household data are of the form (hbase)-rep(int).txt 
# wbase : idem for workplace data

hbase = '../article-data/simulations/initial-condition/household-states/household-stateDist-1I0-pXinsee-bG0p08500-bH0p10000-bW0p00100-nu0p12500-N10000-Tmax150-threshold0p0100-job313'
wbase = '../article-data/simulations/initial-condition/workplace-states/workplace-stateDist-1I0-pXinsee-bG0p08500-bH0p10000-bW0p00100-nu0p12500-N10000-Tmax150-threshold0p0100-job313'


# where to save the result

pwd = '../article-data/limit/'


# Model parameters

betaG = 0.085
lambdaH = 0.1
lambdaW = 0.001
gamma = 0.125
tstop = 150

# Structure size distributions

pH = np.array([0.37, 0.32, 0.13, 0.12, 0.04, 0.02])
pW = np.append([0.2192, 0.12, 0.08, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03], np.repeat(0.02, 3))
pW  = np.append(pW, np.repeat(0.015, 4))
pW = np.append(pW, np.repeat(0.001,25))
pW = np.append(pW, np.repeat(0.0001, 8))
pW = np.append(pW, [0.175])

### MAIN

nHmax = pH.size 
nWmax = pW.size 

# Number of variables in the large population limit dynamical system
vecsize = nHmax*(nHmax + 1)/2 + nWmax*(nWmax + 1)/2
vecsize = round(vecsize)


# Checking which simulations have attaigned the necessary threshold for computing the initial condition (= runs for which the output has been saved)
repnbs = []
for j in range(1,nbrep+1) :
	currentname = hbase + '-rep%i.txt' % j
	if os.path.exists(currentname) :
		repnbs.append(j)
repsize = len(repnbs)

# y will contain all the observed initial conditions
y = np.zeros((repsize, vecsize))
k = 0
# Computing the initial condition for each run and storing it in y
for j in repnbs :
	hfile = hbase + '-rep%i.txt' % j
	wfile = wbase + '-rep%i.txt' % j
	yH = ioput.read_simulated_initcdt(hfile, vecsize, limit.coeffH, args = [nHmax], init = True)
	yW = ioput.read_simulated_initcdt(wfile, vecsize, limit.coeffW, args = [nHmax, nWmax])
	y[k,:] = yH + yW
	k += 1

# average initial condition
y0 = np.mean(y, axis=0)

# computing the reduced model
time, susceptibles, infected = limit.largepopY0(betaG, lambdaH, lambdaW, gamma, pH, pW, tstop, y0)

# saving the output
filename = 'largepop-Y0-bG%0.5f-bH%0.5f-bW%0.5f-nu%0.5f-Tmax%i' %(betaG, lambdaH, lambdaW, gamma, tstop)
filename = filename.replace('.','p')+'.txt'
filename = pwd + filename

ioput.save_reduction_y0(time, susceptibles, infected, filename, pH, pW, betaG, lambdaH, lambdaW, gamma, tstop, y0)
