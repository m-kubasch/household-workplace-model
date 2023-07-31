#!/usr/bin/env python3

import numpy as np
from scipy.integrate import odeint
import math

def binom(n, k):
	""" Computes n choose k """
	return math.factorial(n) // math.factorial(k) // math.factorial(n - k)


def coeffH(s, i, nHmax) :
	""" Computes the position describing a household in state (s,i) given the maximal household size nHmax """
	m = s + i
	c = 1 + nHmax * (nHmax + 1)/2 - m*(m+1)/2 + i + 1
	c = round(c)
	return(c)


def coeffW(s, i, nHmax, nWmax) :
	""" Computes the position describing a workplace in state (s,i) given the maximal household and workplace sizes nHmax and nWmax """
	m = s + i 
	c = nHmax * (nHmax + 1)/2 + nWmax * (nWmax + 1)/2 - m*(m+1)/2 + i + 1
	c = round(c) 
	return(c)


def limit(y, t, betaG, lambdaH, lambdaW, gamma, nHmax, nWmax, NH, NW, vecsize, hstart, hstop, wstart, wstop, prodSIh, prodSIw) :
	""" Computes dy/dt as defined in the large population limit dynamical system.
	
	INPUT
		y, t : current state at time t

		betaG : one-to-all infectious contact rate in the general population
		lambdaH : one-to-one infectious contact rate within households
		lambdaW : one-to-one infectious contact rate within workplaces
		gamma : recovery rate
		
		nHmax : maximum household size
		nWmax : maximum workplace size

		vecsize : number of variables in the dynamical system
		hstart, hstop : indexes delimiting positions corresponding to households
		wstart, wstop : indexes delimiting positions corresponding to workplaces
		prodSIh : prodcut S*I for each epidemic state (S,I) considered for households
		prodSIw : prodcut S*I for each epidemic state (S,I) considered for workplaces
		
	OUTPUT
		dy : dy/dt at time t as given by the dynmical system
	
	"""
	
	
	dy = np.zeros(vecsize)
	S = y[0]
	I = y[1]
	tauH = lambdaH/NH * np.dot(prodSIh, y[hstart:hstop])
	tauW = lambdaW/NW * np.dot(prodSIw, y[wstart:wstop])
	tauG = betaG * I
	
	# Derivatives of the proportion of susceptibles and infected
	dy[0] = -(tauH + tauW + tauG*S)
	dy[1] = -dy[0] - gamma*I
	
	# Derivatives for households - parsing all epidemic states
	for m in range(2, nHmax + 1) :
		for i in range(m) :
			s = m-i
			c = coeffH(s, i, nHmax)
			nh = y[c]
			dnh = - (lambdaH * s * i + tauW/S * s + tauG * s + gamma * i) * nh
			if i >= 1 : 
				sI = s+1
				iI = i-1
				cI = coeffH(sI, iI, nHmax)
				nhI = y[cI]
				dnh += (lambdaH * sI * iI + tauW/S * sI + tauG * sI) * nhI
			if i + s < nHmax :
				iG = i + 1
				cG = coeffH(s, iG, nHmax)
				nhG = y[cG]
				dnh += gamma * iG * nhG 
			dy[c] = dnh
			
	# Derivatives for workplaces - parsing all epidemic states
	for m in range(2, nWmax + 1) :
		for i in range(m) :
			s = m-i
			c = coeffW(s, i, nHmax, nWmax)
			nw = y[c]
			dnw = -(lambdaW * s * i + tauH/S * s + tauG * s + gamma * i) * nw
			if i >= 1 :
				sI = s+1
				iI = i-1
				cI = coeffW(sI, iI, nHmax, nWmax)
				nwI = y[cI]
				dnw += (lambdaW * sI * iI + tauH/S * sI + tauG * sI) * nwI
			if i + s < nWmax :
				iG = i + 1
				cG = coeffW(s, iG, nHmax, nWmax)
				nwG = y[cG]
				dnw += gamma * iG * nwG 
			dy[c] = dnw
			
	return(dy)


def largepop(betaG, lambdaH, lambdaW, gamma, epsilon, pH, pW, Tmax) :
	""" Computes the solution to the large population limit dynamical system, starting from a given proportion of infected chosen uniformly at random in an otherwise susceptible population at time 0. 

	INPUT 
		betaG : one-to-all infectious contact rate in the general population
		lambdaH : one-to-one infectious contact rate within households
		lambdaW : one-to-one infectious contact rate within workplaces
		gamma : recovery rate
		
		epsilon : proportion of infected at time 0
		
		pH : household size distribution
		pW : workplace size distribution 
		
		Tmax : time up to which the solution is computed

	OUTPUT 
		t, s, i : s[k], i[k] proportions of susceptibles and infected at time t[k].
	"""
	
	
	# Times at which the solution to the dynamical system is evaluated 
	spacing = 1000
	t  = np.linspace(0, Tmax, spacing)
	
	# Maximal household and workplace sizes under consideration 
	nHmax = pH.size 
	nWmax = pW.size 
	
	# Average household and workplace size
	NH = np.dot(np.arange(1,nHmax+1), pH)
	NW = np.dot(np.arange(1,nWmax+1), pW)
	
	# Number of variables in the dynamic system
	vecsize = nHmax*(nHmax + 1)/2 + nWmax*(nWmax + 1)/2
	vecsize = round(vecsize)
	
	# Households = y[hstart:hstop]
	hstart = 2
	hstop = round(nHmax*(nHmax + 1)/2) + 1
	
	# Workplaces = y[wstart:wstop]
	wstart = hstop
	wstop = wstart + round(nWmax*(nWmax + 1)/2)
	
	# Vector of products (S * I) for each epidemic state and Initial condition y0
	prodSI = np.zeros(vecsize)
	
	y0 = np.zeros(vecsize)
	y0[0] = 1 - epsilon # Proportion of susceptibles
	y0[1] = epsilon # Proportion of infected

	# Households - parsing all epidemic states
	for m in range(2, nHmax + 1) :
		for i in range(m) :
			s = m-i
			c = coeffH(s, i, nHmax)
			y0[c] = pH[m-1] * binom(m, s) * (1 - epsilon) ** s * epsilon ** i
			prodSI[c] = s*i
	
	# Workplaces - parsing all epidemic states
	for m in range(2, nWmax + 1) :
		for i in range(m) :
			s = m-i
			c = coeffW(s, i, nHmax, nWmax)
			y0[c] = pW[m-1] * binom(m, s) * (1 - epsilon) ** s * epsilon ** i
			prodSI[c] = s*i
			
	
	# Seperating households from workplaces
	prodSIh = prodSI[hstart:hstop]
	prodSIw = prodSI[wstart:wstop]
	
	# Solving the dynamical system
	solv = odeint(limit, y0, t, args = (betaG, lambdaH, lambdaW, gamma, nHmax, nWmax, NH, NW, vecsize, hstart, hstop, wstart, wstop, prodSIh, prodSIw))
	
	# Extract the proportion of susceptibles and infected
	s = solv[:, 0]
	i = solv[:, 1]
	
	return(t, s, i)


def largepopY0(betaG, lambdaH, lambdaW, gamma, pH, pW, Tmax, y0) :
	""" Computes the solution to the large population limit dynamical system, starting from a given initial condition at time 0. 

	INPUT 
		betaG : one-to-all infectious contact rate in the general population
		lambdaH : one-to-one infectious contact rate within households
		lambdaW : one-to-one infectious contact rate within workplaces
		gamma : recovery rate
		
		pH : household size distribution
		pW : workplace size distribution 
		
		Tmax : time up to which the solution is computed

		y0 : initial condition

	OUTPUT 
		t, s, i : s[k], i[k] proportions of susceptibles and infected at time t[k].
	"""
	
	
	# Times at which the solution to the dynamical system is evaluated 
	spacing = 1000
	t  = np.linspace(0, Tmax, spacing)
	
	# Maximal household and workplace sizes under consideration 
	nHmax = pH.size 
	nWmax = pW.size 
	
	# Average household and workplace size
	NH = np.dot(np.arange(1,nHmax+1), pH)
	NW = np.dot(np.arange(1,nWmax+1), pW)
	
	# Number of variables in the dynamic system
	vecsize = nHmax*(nHmax + 1)/2 + nWmax*(nWmax + 1)/2
	vecsize = round(vecsize)
	
	# Households = y[hstart:hstop]
	hstart = 2
	hstop = round(nHmax*(nHmax + 1)/2) + 1
	
	# Workplaces = y[wstart:wstop]
	wstart = hstop
	wstop = wstart + round(nWmax*(nWmax + 1)/2)
	
	# Vector of products (S * I) for each epidemic state and Initial condition y0
	prodSI = np.zeros(vecsize)
	
	# Households - parsing all epidemic states
	for m in range(2, nHmax + 1) :
		for i in range(m) :
			s = m-i
			c = coeffH(s, i, nHmax)
			prodSI[c] = s*i
			
	# Workplaces - parsing all epidemic states
	for m in range(2, nWmax + 1) :
		for i in range(m) :
			s = m-i
			c = coeffW(s, i, nHmax, nWmax)
			prodSI[c] = s*i
			
			
	# Seperating households from workplaces
	prodSIh = prodSI[hstart:hstop]
	prodSIw = prodSI[wstart:wstop]
	
	# Solving the dynamical system
	solv = odeint(limit, y0, t, args = (betaG, lambdaH, lambdaW, gamma, nHmax, nWmax, NH, NW, vecsize, hstart, hstop, wstart, wstop, prodSIh, prodSIw))
	
	# Extract the proportion of susceptibles and infected
	s = solv[:, 0]
	i = solv[:, 1]
	
	return(t, s, i)