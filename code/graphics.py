#!/usr/bin/env python3

import ioput
import numpy as np
import matplotlib.pyplot as plt

def plotsimu(ax, t, S, I, threshold = 0, scol = 'yellowgreen', icol = 'lightcoral', alpha = 0.3, lst = 'solid', negtime = False) :
	"""Plot all epidemic trajectories of the simulated model stored in a matrix, each row containing the proportion of susceptibles (S) (resp. infected (I)) in the population for each timestep of t. 
	Time is shifted so that at time 0, the proportion of infected exceeds the indicated threshold for the first time.
	If negtime is True, then the trajectories up to time 0 will also be plotted."""
	
	nbsteps = t.size
	nrep = S.shape[0]
	for j in range(nrep) :
		i = 0
		trying = True
		while I[j,i] < threshold and trying :
			i += 1
			if i == nbsteps - 1 :
				trying = False

		if trying : 
			if negtime :
				ax.plot(t - t[i], S[j,:], scol, linestyle = lst, alpha=alpha)
				ax.plot(t - t[i], I[j,:], icol, linestyle = lst, alpha=alpha)
			else :
				ax.plot(t[i:] - t[i], S[j,i:], scol, linestyle = lst, alpha=alpha)
				ax.plot(t[i:] - t[i], I[j,i:], icol, linestyle = lst, alpha=alpha)
	return()


def plotreduction(ax, t, S, I, threshold = 0, scol = 'green', icol = 'red', lst = 'solid', negtime = False) :
	"""Plot one trajectory of a deterministic model, with timeshift so that at time 0, the proportion of infected exceeds the indicated threshold for the first time. 
	If negtime is True, then the trajectories up to time 0 will also be plotted."""
	
	i = 0
	while I[i] < threshold :
		i += 1
	
	if negtime :
		ax.plot(t - t[i], S, scol, linestyle = lst)
		ax.plot(t - t[i], I, icol, linestyle = lst)
	else :
		ax.plot(t[i:] - t[i], S[i:], scol, linestyle = lst)
		ax.plot(t[i:] - t[i], I[i:], icol, linestyle = lst)
	
	return()


def plotscenarios(ax, filenames, r0, colorS, colorI) : 
	""" Plotting several trajectories, as in Figure A1. 
	Compute and return the times at which those trajectories fall below the threshold of 1% of infetcted, after the epidemic peak. """
	
	stoptimes = []
	stopvalues = []
	
	for fname, R0 in zip(filenames, r0) :
		_,t,s,i = ioput.read_reduction(fname)
		stopindex = np.max(np.where(i >= 0.01))
		stoptimes.append(t[int(stopindex)])
		stopvalues.append(i[int(stopindex)])
		
		ax.plot(t, s, label = 'S %0.2f' % R0, c = colorS[R0])
		ax.plot(t, i, label = 'I %0.2f' % R0, c = colorI[R0])
		
	ax.scatter(stoptimes, stopvalues, marker = 'P', color = 'black', zorder=2)
		
	return(stoptimes)
	