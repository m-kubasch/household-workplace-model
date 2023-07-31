#!/usr/bin/env python3

import ioput

import numpy as np
import time

def reference(N) :
	"""Short reference function for normalizing computation time fluctuations. """
	k = 0
	for j in range(N) :
		k += 1
	return(k)


def timefct(func, param, n=1000000000) :
	"""Computation time for one execution of the reference function reference(n), and of the function of interest func(*param). """
	
	start = time.time()
	reference(n)
	end = time.time()
	reftime = end-start
	
	start = time.time()
	func(*param)
	end = time.time()
	ftime = end-start
	
	return(reftime, ftime)
	

def syntheticInfo(filelist, nbrep) :
	"""Computes the ratio of the average compution time for the reference function and function of interest, for each set of data in filelist, each dataset containing nbrep runs. 
	
	INPUT
		filelist, nbrep : considered files are of the form fbase-rep1.txt, ..., fbase-rep{nbrep}.txt for fbase in filelist.
		
	OUTPUT 
		averageRatio : list of ratios of the average compution times
	"""
	averageRatio = []
	for fbase in filelist : 
		rfc, ftime = ioput.read_timingdata(fbase, nbrep)
		averageRatio.append(np.mean(ftime/rfc))
	averageRatio = np.array(averageRatio)
	return(averageRatio)

def timeratio(simlist, redlist, nbrep) :
	"""Returns the ratio of average normalized computation times for simulations and reductions.
	
	INPUT
		similist, redlist : datasets of computation time for simulations and reduced model, respectively
		nbrep : number of repeats in each dataset
		-- considered files are of the form fbase-rep1.txt, ..., fbase-rep{nbrep}.txt for fbase in simlist and redlist, respectively
	
	OUTPUT
		totalRatio : list of normalized average computation times (reduced model / simulation)
	"""
	
	simRatios = syntheticInfo(simlist, nbrep)
	redRatios = syntheticInfo(redlist, nbrep)
	totalRatio = redRatios/simRatios
	return(totalRatio)

