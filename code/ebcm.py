#!/usr/bin/env python3

import numpy as np
from scipy.integrate import odeint


def coeffH(s, i, h, hmax, wmax) :
	""" Computes the position describing a household of size h in state (s,i) given the maximal household and workplace sizes (hmax and wmax, resp.) """
	c = hmax + wmax + (h-1)*h*( (2*h-1)/12 + 1/4 ) - (h-1) + h*(h+1)/2 - (s+i)*(s+i+1)/2 + i
	c = round(c)
	return(c)


def coeffW(s, i, w, hmax, wmax) :
	""" Computes the position describing a workplace of size w in state (s,i) given the maximal household and workplace sizes (hmax and wmax, resp.) """
	c = wmax + hmax*(hmax+1)*( (2*hmax+1)/12 + 1/4 ) + (w-1)*w*( (2*w-1)/12 + 1/4 ) - (w-1) + w*(w+1)/2 - (s+i)*(s+i+1)/2 + i
	c = round(c) 
	return(c)


def statespace(n) :
	""" All the states (s,i) to consider for a structure of size n. """
	states = []
	for s in range(1, n) :
		for i in range(1, n) :
			if s + i <= n :
				states.append((s,i))
	for s in range(2,n+1) :
		states.append((s, 0))
	return(states)



def dy_ebcm(y, t, betaG, lambdaH, lambdaW, gamma, pHb, pWb, hmax, wmax, athetaH1, athetaW1) :
	""" Computes dy/dt as defined in the edge-based compartimental model.
	
	INPUT
		y, t : current state at time t

		betaG : one-to-all infectious contact rate in the general population
		lambdaH : one-to-one infectious contact rate within households
		lambdaW : one-to-one infectious contact rate within workplaces
		gamma : recovery rate
		
		pHb, pWb : size-biased household and workplace size distributions
		
		hmax, wmax : maximum household and workplace sizes

		athetaH1, athetaW1 : values of theta^H_1 and theta^W_1 which are constant over time
		
	OUTPUT
		dy : dy/dt at time t as given by the dynmical system
	
	"""
	
	
	dy = np.zeros(y.size)
	I = y[0]
	thetaG = y[1]
	thetaH = np.append(athetaH1, y[2:hmax+1]) # inserting thetaH1
	thetaW = np.append(athetaW1, y[hmax+1:hmax+wmax]) # inserting thetaW1
	
	sumptW = np.array([pHb[h] * np.dot(pWb, thetaW) for h in range(hmax)])
	sumptH = np.array([pWb[w] * np.dot(pHb,thetaH) for w in range(wmax)])
	
	nH = thetaG * np.array([thetaH[h] * sumptW[h] for h in range(hmax)])
	nW = thetaG * np.array([thetaW[w] * sumptH[w] for w in range(wmax)])
	
	SHW = np.dot(pHb, thetaH)*np.dot(pWb, thetaW)
	S = thetaG * SHW
	
	TG = betaG * S * I
	
	TH = lambdaH * np.array([np.sum(np.array([s*i*y[coeffH(s, i, h, hmax, wmax)] for s,i in statespace(h)])) for h in range(2, hmax+1)])
	TW = lambdaW * np.array([np.sum(np.array([s*i*y[coeffW(s, i, w, hmax, wmax)] for s,i in statespace(w)])) for w in range(2, wmax+1)])
	
	
	# within-structure infections for households/workplaces of size 1 are impossible
	TH = np.append(np.zeros(1, dtype = 'int'), TH)
	TW = np.append(np.zeros(1, dtype = 'int'), TW)
	
	tauH = np.array([
		TG * thetaH[h] * sumptW[h]/SHW
		+ pHb[h] * np.sum(np.array([TW[w] * pWb[w] * thetaH[h]/sumptH[w] for w in range(wmax)])) 
		for h in range(hmax)
	])
	
	tauW = np.array([
		TG * thetaW[w] * sumptH[w]/SHW
		+ pWb[w] * np.sum(np.array([TH[h] * pHb[h] * thetaW[w]/sumptW[h] for h in range(hmax)])) 
		for w in range(wmax)
	])
	
	
	# derivatives of the probability of escaping infection
	
	# dt_thetaH for households of size >= 2
	for h in range(hmax - 1) :
		dy[2+h] = - thetaH[h + 1] * TH[h + 1] / nH[h + 1]
		
	# dt_thetaW for workplaces of size >= 2
	for w in range(wmax - 1) :
		dy[hmax+1+w] = -thetaW[w + 1] * TW[w + 1] / nW[w + 1]
		
	# dt_thetaG
	dy[1] = - thetaG * betaG * I
	
	# derivatives for households - parsing all epidemic states
	for h in range(2, hmax + 1) :
		tauHh = tauH[h-1]
		nHh = nH[h-1]
		
		for s, i in statespace(h) :
			c = coeffH(s, i, h, hmax, wmax)
			dy[c] = -(lambdaH * s + gamma) * i * y[c] - tauHh / nHh * s * y[c] 
			if i >= 1 :
				p = coeffH(s+1, i-1, h, hmax, wmax)
				dy[c] += lambdaH * (s+1) * (i-1) * y[p] + tauHh / nHh * (s+1) * y[p] 
			if h - s - i >= 1 :
				g = coeffH(s, i+1, h, hmax, wmax)
				dy[c] += gamma * (i+1) * y[g]
				
	# derivatives for workplaces - parsing all epidemic states
	for w in range(2, wmax + 1) :
		tauWw = tauW[w-1]
		nWw = nW[w-1]
		
		for s, i in statespace(w) :
			c = coeffW(s, i, w, hmax, wmax)
			dy[c] = -(lambdaW * s + gamma) * i * y[c] - tauWw / nWw * s * y[c] 
			if i >= 1 : 
				p = coeffW(s+1, i-1, w, hmax, wmax)
				dy[c] += lambdaW * (s+1) * (i-1) * y[p] + tauWw / nWw * (s+1) * y[p] 
			if w - s - i >= 1 :
				g = coeffW(s, i+1, w, hmax, wmax)
				dy[c] += gamma * (i+1) * y[g]
				
				
	# thetaH1 and thetaW1 are constant
	dt_thetaH = np.append(np.zeros(1), dy[2:hmax+1])
	dt_thetaW = np.append(np.zeros(1), dy[hmax+1:hmax+wmax])
	dt_thetaG = dy[1]
	
	# derivative of the proportion of infected
	dy[0] = - gamma * I - dt_thetaG * SHW - thetaG * (np.dot(pHb,dt_thetaH) * np.dot(pWb, thetaW) + np.dot(pHb, thetaH) * np.dot(pWb, dt_thetaW))
	
	return(dy)


def ebcm(betaG, lambdaH, lambdaW, gamma, epsilon, pH, pW, Tmax) :
	""" Computes the solution to the edge-based compartmental model. 

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
	
	# Size-biased distributions
	pHb = np.array([(k+1)*p for k, p in enumerate(pH)])
	pWb = np.array([(k+1)*p for k, p in enumerate(pW)])
	pHb = pHb/np.sum(pHb)
	pWb = pWb/np.sum(pWb)
	
	# Maximal household and workplace sizes under consideration
	hmax = pHb.size
	wmax = pWb.size
	
	# Some constants for computations 
	epsbar = 1 - epsilon
	thetaH1 = epsbar
	thetaW1 = epsbar
	athetaH1 = np.array(thetaH1)
	athetaW1 = np.array(thetaW1)
	
	# Number of variables in the dynamic system
	vecsize = hmax + wmax + hmax*(hmax+1)*((2*hmax+1)/12 + 1/4) - hmax + wmax*(wmax+1)*((2*wmax+1)/12 + 1/4) - wmax
	vecsize = round(vecsize)
	
	# Correspondance between structure types and their position in the vector y
	dstatespace = {n : statespace(n) for n in range(2, max(hmax, wmax)+1)}
	dcoeffH = {(s,i,h) : coeffH(s, i, h, hmax, wmax) for h in range(2, hmax+1) for s,i in dstatespace[h]}
	dcoeffW = {(s,i,w) : coeffW(s, i, w, hmax, wmax) for w in range(2, wmax+1) for s,i in dstatespace[w]}
	
	
	# Computing the initial condition y(0)
	y0 = np.zeros(vecsize)
	y0[0] = epsilon # I
	y0[1:hmax+wmax] = epsbar # thetaG, thetaHh, thetaWw
	
	# Households - parsing all possible types
	for h in range(2, hmax + 1) :
		for s, i in dstatespace[h] :
			if s+i == h:
				c = dcoeffH[(s, i, h)]
				norm = 1
				if h == s :
					norm = h
				y0[c] = pHb[h-1]/norm * (epsilon**i) * (epsbar**s)

	# Workplaces - parsing all possible types 
	for w in range(2, wmax + 1) :
		for s, i in dstatespace[w] :
			if s+i == w:
				c = dcoeffW[(s, i, w)]
				norm = 1
				if w == s :
					norm = w
				y0[c] = pWb[w-1]/norm * (epsilon**i) * (epsbar**s)
	
	# Solving the dynamical system 
	solv = odeint(dy_ebcm, y0, t, args = (betaG, lambdaH, lambdaW, gamma, pHb, pWb, hmax, wmax, athetaH1, athetaW1))
	
	I = solv[:, 0]
	thetaG = solv[:, 1]
	thetaH = solv[:, 2:hmax+1]
	thetaW = solv[:, hmax+1:hmax+wmax]
	
	# Computing the proportion of susceptibles over time
	vh = pHb[0]*thetaH1
	vw = pWb[0]*thetaW1
	S = np.array([thetaG[j]*(vh + np.dot(thetaH[j,:], pHb[1:]))*(vw + np.dot(thetaW[j,:], pWb[1:])) for j in range(spacing)])
	
	return(t, S, I)