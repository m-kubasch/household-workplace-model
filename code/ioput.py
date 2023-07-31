#!/usr/bin/env python3

import numpy as np
import os


### OUTPUT

def save_simulation(data, filename, K, pH, pW, betaG, lbdaH, lbdaW, gamma, eps, Tmax, explosion = 0.005) :
	""" Save simulated data in a .txt file. 

	INPUT 
		data : output from simulation
		filename : desired name of the output file
		K, pH, pW : population size, housheold and workplace size distributions used for the simulation
		betaG, lbdaH, lbdaW, gamma : epidemic parameters used for the simulation
		eps : initial proportion of infected 
		Tmax : final time
		explosion : data is only solved if the maximal proportion of infected exceeds this threshold value

	OUTPUT FORMAT 
		Line 1: parameters
		Line 2: household size distribution
		Line 3: workplace size distribution
		Line 4: time 
		Line 5: proportion of susceptibles
		Line 6: proportion of infected 
		Line 7: proporiton of recovered
	"""
	
	tps, pops = data
	nbS = np.array([pop[0] for pop in pops])/K
	nbI = np.array([pop[1] for pop in pops])/K
	nbR = np.array([pop[2] for pop in pops])/K

	if np.max(nbI) > explosion:
		outfile = open(filename, "w")
		outfile.write( "pop size : %i, betaG : %0.5f, lambdaH : %0.5f, lambdaW : %0.5f, gamma : %0.5f, epsilon : %0.5f, Tmax : %0.2f\n" %(K, betaG, lbdaH, lbdaW, gamma, eps, Tmax) )
		outfile.writelines(["%s " %p for p in pH])
		outfile.write("\n")
		outfile.writelines(["%s " %p for p in pW])
		outfile.write("\n")
		outfile.writelines(["%s " %t for t in tps])
		outfile.write("\n")
		outfile.writelines(["%s " %s for s in nbS])
		outfile.write("\n")
		outfile.writelines(["%s " %s for s in nbI])
		outfile.write("\n")
		outfile.writelines(["%s " %s for s in nbR])
		outfile.write("\n")
		outfile.close()

	return()


def save_simulated_initcdt(data, threshold, hfilename, wfilename, K, pH, pW, betaG, lbdaH, lbdaW, gamma, Tmax) :
	""" Save simulated initial condition data in a .txt file. 

	INPUT 
		data : output from simulation
		threshold : proportion of infected that needs to be attained for the initial condition
		hfilename, wfilename : desired name of the output files for household and workplca type dsitributions, respectively
		K, pH, pW : population size, housheold and workplace size distributions used for the simulation
		betaG, lbdaH, lbdaW, gamma : epidemic parameters used for the simulation
		Tmax : final time

	OUTPUT FORMAT 
		Line 1: parameters
		Line 2: household size distribution
		Line 3: workplace size distribution
		Line 4: "Proportion of infected at cutoff : ", proportion of infected 
		Line 5: "Proportion of susceptibles at cutoff : ", proportion of susceptibles
		Line 6+: "(s,i)" : proportion of structures (housheolds or wrokplaces, resp.) in state (s,i)
	"""
	
	prop_inf, prop_susc, hstate_dist, wstate_dist = data
	
	if prop_inf >= threshold:
		outfile = open(hfilename, "w")
		outfile.write( "pop size : %i, betaG : %0.5f, lambdaH : %0.5f, lambdaW : %0.5f, gamma : %0.5f, Tmax : %0.2f, threshold of infected : %0.3f\n" %(K, betaG, lbdaH, lbdaW, gamma, Tmax, threshold) )
		outfile.writelines(["%s " %p for p in pH])
		outfile.write("\n")
		outfile.writelines(["%s " %p for p in pW])
		outfile.write("\n")
		outfile.write("Proportion of infected at cutoff : %0.4f\n" %prop_inf)
		outfile.write("Proportion of susceptibles at cutoff : %0.4f\n" %prop_susc)
		for key, value in hstate_dist.items() :
			outfile.write("%s : %i\n" %(key, value))
		outfile.close()
		
		outfile = open(wfilename, "w")
		outfile.write( "pop size : %i, betaG : %0.5f, lambdaH : %0.5f, lambdaW : %0.5f, gamma : %0.5f, Tmax : %0.2f, threshold of infected : %0.3f\n" %(K, betaG, lbdaH, lbdaW, gamma, Tmax, threshold) )
		outfile.writelines(["%s " %p for p in pH])
		outfile.write("\n")
		outfile.writelines(["%s " %p for p in pW])
		outfile.write("\n")
		outfile.write("Proportion of infected at cutoff : %0.4f\n" %prop_inf)
		outfile.write("Proportion of susceptibles at cutoff : %0.4f\n" %prop_susc)
		for key, value in wstate_dist.items() :
			outfile.write("%s : %i\n" %(key, value))
		outfile.close()
		
	return()


def save_reduction(time, s, i, filename, pH, pW, betaG, lbdaH, lbdaW, gamma, eps, Tmax) :
	""" Save solution of dynamical system in a .txt file. 
	
	INPUT 
		time : time steps at which the reduced model is evaluated
		s, i : output from reduced model - proporiton of usceptibles and infected, respectively
		filename : desired name of the output file
		pH, pW : housheold and workplace size distributions used for the simulation
		betaG, lbdaH, lbdaW, gamma : epidemic parameters used for the simulation
		eps : initial proportion of infected 
		Tmax : final time
	
	OUTPUT FORMAT 
		Line 1: parameters
		Line 2: household size distribution
		Line 3: workplace size distribution
		Line 4: time 
		Line 5: proportion of susceptibles
		Line 6: proportion of infected 
	"""
	
	outfile = open(filename, "w")
	outfile.write( "betaG : %0.5f, lambdaH : %0.5f, lambdaW : %0.5f, gamma : %0.5f, epsilon : %0.5f, Tmax : %0.2f\n" %(betaG, lbdaH, lbdaW, gamma, eps, Tmax) )
	outfile.writelines(["%s " %p for p in pH])
	outfile.write("\n")
	outfile.writelines(["%s " %p for p in pW])
	outfile.write("\n")
	outfile.writelines(["%s " %x for x in time])
	outfile.write("\n")
	outfile.writelines(["%s " %x for x in s])
	outfile.write("\n")
	outfile.writelines(["%s " %x for x in i])
	outfile.write("\n")
	outfile.close()
	
	return() 

def save_reduction_y0(time, s, i, filename, pH, pW, betaG, lbdaH, lbdaW, gamma, Tmax, y0) :
	""" Save solution of dynamical system and its initial condition in a .txt file. 
	
	INPUT 
		time : time steps at which the reduced model is evaluated
		s, i : output from reduced model - proporiton of usceptibles and infected, respectively
		filename : desired name of the output file
		pH, pW : housheold and workplace size distributions used for the simulation
		betaG, lbdaH, lbdaW, gamma : epidemic parameters used for the simulation
		Tmax : final time
		y0 : initial condition
	
	OUTPUT FORMAT 
		Line 1: parameters
		Line 2: household size distribution
		Line 3: workplace size distribution
		Line 4: initial condition
		Line 5: time 
		Line 6: proportion of susceptibles
		Line 7: proportion of infected 
	"""
	
	outfile = open(filename, "w")
	outfile.write( "betaG : %0.5f, lambdaH : %0.5f, lambdaW : %0.5f, gamma : %0.5f, Tmax : %0.2f\n" %(betaG, lbdaH, lbdaW, gamma, Tmax) )
	outfile.writelines(["%s " %p for p in pH])
	outfile.write("\n")
	outfile.writelines(["%s " %p for p in pW])
	outfile.write("\n")
	outfile.writelines(["%s " %y for y in y0])
	outfile.write("\n")
	outfile.writelines(["%s " %x for x in time])
	outfile.write("\n")
	outfile.writelines(["%s " %x for x in s])
	outfile.write("\n")
	outfile.writelines(["%s " %x for x in i])
	outfile.write("\n")
	outfile.close()
	
	return() 

def save_ctime(reftime, ftime, filename, pH, pW, betaG, lbdaH, lbdaW, gamma, eps, Tmax) :
	""" Save computation time data in a .txt file. 
	
	INPUT 
		reftime, ftime : computation time for executing the reference function / function of interest
		filename : desired name of the output file
		pH, pW : housheold and workplace size distributions used for the simulation
		betaG, lbdaH, lbdaW, gamma : epidemic parameters used for the simulation
		eps : initial proportion of infected 
		Tmax : final time
	
	OUTPUT FORMAT 
		Line 1: parameters
		Line 2: household size distribution
		Line 3: workplace size distribution
		Line 4: 'reference time : ' reftime 
		Line 5: 'main time : ' ftime
	"""
	
	outfile = open(filename, "w")
	outfile.write( "betaG : %0.5f, lambdaH : %0.5f, lambdaW : %0.5f, gamma : %0.5f, epsilon : %0.5f, Tmax : %0.2f\n" %(betaG, lbdaH, lbdaW, gamma, eps, Tmax) )
	outfile.writelines(["%s " %p for p in pH])
	outfile.write("\n")
	outfile.writelines(["%s " %p for p in pW])
	outfile.write("\n")
	outfile.write("reference time : %0.5f\n" % reftime)
	outfile.write("main time : %0.5f\n" % ftime)
	outfile.close()
	
	return() 


### INPUT

def read_simarray(root, nbrep) :
	""" Extract information from simulation files.
		
		INPUT : 
			root, nbrep : root = 'pwd/base' such that all simulation files are fo the form 'pwd/base-repX.txt' for 1 <= X <= nbrep
		
		OUTPUT :
			parameters : dictionnary containing parameters used for simulation
			simtemps : list of times at which the simulation is evaluated
			S, I : proportions of susceptibels and infected at each time, for each simulation
	
	"""
	repnbs = []
	# Checking which simulations actually have succeeded 
	for j in range(nbrep) :
		currentname = root + '-rep%i.txt' % j
		if os.path.exists(currentname) :
			repnbs.append(j)
	repsize = len(repnbs)
	
	# Extracting information from those that have succeeded
	currentname = root + '-rep%i.txt' % repnbs[0]
	file = open(currentname,"r")
	data = file.readlines()
	temps = []
	
	# Parameters - common for all simulations
	parameters = {}
	raw_param = data[0].split(',')
	for string in raw_param :
		name, val = string.split(' : ')
		name = name.replace(' ', '')
		val = float(val)
		parameters[name] = val
	pH = np.array([float(p) for p in data[1].split(" ")[:-1]])
	pW = np.array([float(p) for p in data[2].split(" ")[:-1]])
	parameters['pH'] = pH
	parameters['pW'] = pW
	
	# Extracting the times at which simulations are evaluated - common for all simulations
	simtemps = np.array([float(t) for t in data[3].split(" ")[:-1]])
	
	# Proportions of susceptibles and infected 
	S = np.zeros((repsize, simtemps.size))
	I = np.zeros((repsize, simtemps.size))
	S[0,:] = [float(s) for s in data[4].split(" ")[:-1]]
	I[0,:] = [float(i) for i in data[5].split(" ")[:-1]]
	file.close()
	
	for k in range(1, repsize) :
		currentname = root + '-rep%i.txt' % repnbs[k]
		file = open(currentname,"r")
		data = file.readlines()
		S[k-1,:] = [float(s) for s in data[4].split(" ")[:-1]]
		I[k-1,:] = [float(i) for i in data[5].split(" ")[:-1]]
		file.close()
		
	return(parameters, simtemps, S, I)



def read_simulated_initcdt(fname, vsize, cfct, args, init=False, params = False):
	"""Extract information form output files generated using save_simulated_initcdt.

	INPUT
		fname : input filename
		vsize : number of variables of the large population limit dynamical system for which the initial value is computed
		cfct : function associating the correct vector index to each structure type (s,i)
		args : other arguments for cfunc (except for s and i)
		init : if True, the proportion of suscpetibles and infected is stored in the initial condition (needs only to be done once for either households or workplaces)
		params : if True, extract model parameters and size distributions

	OUTPUT 
		y0 : initial condtion obtained in this simulation for the structure under consideration (either households or workplaces)
	"""
	
	rawfile = open(fname,"r")
	data = rawfile.readlines()
	propI = float(data[3].split(' : ')[1][:-1])
	propS = float(data[4].split(' : ')[1][:-1])
	states = {}
	for line in data[5:]:
		state, count = line.split(':')
		s, i = state.split(',')
		s = int(s[1:]) # on vire la parenthese
		i = int(i[:-2]) # on vire la parenthese et l'espace
		count = int(count[:-1]) 
		key = (s,i) 
		if key in states : 
			states[key] += count
		else :
			states[key] = count
	rawfile.close()
	nbcliques = np.sum(list(states.values()))
	
	y0 = np.zeros(vsize)
	if init :
		y0[1] = propI 
		y0[0] = propS 
	for key, value in states.items() :
		s, i = key
		if s >= 2 or s*i >= 1 :
			c = cfct(s,i,*args)
			y0[c] = value/nbcliques
	return(y0)


def read_reduction(filename):
	""" Extract data from reduced model file.
	
	INPUT :
		filename : name of the file containing the data

	OUTPUT : 
		parameters : dictionnary containing model parameters 
		t : list of times at which the reduced model is evaluated
		s, i : proportions of suceptibels and infected at each time
	
	"""
	file = open(filename,"r")
	data = file.readlines()
	
	# Parameters
	raw_param = data[0].split(',')
	parameters = {}
	for string in raw_param :
		name, val = string.split(' : ')
		name = name.replace(' ', '')
		val = float(val)
		parameters[name] = val
	pH = np.array([float(p) for p in data[1].split(" ")[:-1]])
	pW = np.array([float(p) for p in data[2].split(" ")[:-1]])
	parameters['pH'] = pH
	parameters['pW'] = pW
	
	# Reduction
	t = np.array([float(t) for t in data[3].split(" ")[:-1]])
	s = np.array([float(s) for s in data[4].split(" ")[:-1]])
	i = np.array([float(i) for i in data[5].split(" ")[:-1]])
	
	file.close()
	return(parameters, t, s, i)

def read_reduction_y0(filename):
	""" Extract data from reduced model file with saved initial condition.
	
	INPUT :
		filename : name of the file containing the data

	OUTPUT : 
		parameters : dictionnary containing model parameters
		t : list of times at which the reduced model is evaluated
		s, i : proportions of suceptibels and infected at each time
	
	"""
	file = open(filename,"r")
	data = file.readlines()
	
	# Parameters
	raw_param = data[0].split(',')
	parameters = {}
	for string in raw_param :
		name, val = string.split(' : ')
		val = float(val)
		parameters[name] = val
	pH = np.array([float(p) for p in data[1].split(" ")[:-1]])
	pW = np.array([float(p) for p in data[2].split(" ")[:-1]])
	y0 = np.array([float(y) for y in data[3].split(" ")[:-1]])
	parameters['pH'] = pH
	parameters['pW'] = pW
	parameters['y0'] = y0
	
	# Reduction
	t = np.array([float(t) for t in data[4].split(" ")[:-1]])
	s = np.array([float(s) for s in data[5].split(" ")[:-1]])
	i = np.array([float(i) for i in data[6].split(" ")[:-1]])
	
	file.close()
	return(parameters, t, s, i)


def read_timingdata(fbase, nbrep) :
	""" Extract data from computation time file.
	
	INPUT :
		fbase : root of the file names containing the data
		nbrep : number of repetitions of the script 
		-- Considered files are of the form fbase-rep1.txt, ..., fbase-rep{nbrep}.txt

	OUTPUT : 
		rfc, ftime : lists of computing times for each run, for the reference function and function of interest.
	
	"""
	rfc = [] # temps mis par la fonction de reference
	ftime = []  # temps mis pour faire la partie d'interet
	
	for j in range(1, nbrep+1) :
		fname = fbase + '-rep%i.txt' %j
		file = open(fname, "r")
		rawdata = file.readlines()
		file.close()
		
		t = float(rawdata[3].split(' ')[-1])
		rfc.append(t)

		t = float(rawdata[4].split(' ')[-1])
		ftime.append(t)
	
	ftime = np.array(ftime)
	rfc = np.array(rfc)
	
	return(rfc, ftime)