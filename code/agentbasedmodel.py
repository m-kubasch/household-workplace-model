#!/usr/bin/env python3

import numpy as np
import collections

def generate_population_hw(pH, pW, N, seed) :
	""" Generates two adjacency matrixes representing social structures (households and workplaces).
	
	INPUT
		pH : household size distribution (np.array whose k-th coefficient is the poporiton of households of size k)
		pW : workplace size distribution
		N : population size
		seed : seed for RNG (or RNG itself)
	
	OUTPUT
		H, W : household (resp. workplace) adjacency matrix 
	"""
	rng = np.random.default_rng(seed)
	hmax = pH.size
	wmax = pW.size
	struc = np.zeros((N, 2))
	
	# households
	remaining = N
	hnumber = 1
	index = 0
	while remaining > 0 : 
		size = 1 + rng.choice(hmax, p = pH)
		if size > remaining :
			size = remaining 
		struc[index : index+size,0] = np.repeat(hnumber, size)
		index += size
		hnumber += 1
		remaining += -size
		
	# workplaces
	wnumber = 1
	index = 0
	remaining = N
	while remaining > 0 : 
		size = 1 + rng.choice(wmax, p = pW)
		if size > remaining :
			size = remaining 
		struc[index : index+size,1] = np.repeat(wnumber, size)
		index += size
		wnumber += 1
		remaining += -size
	
	# random pairings between households and workplaces
	rng.shuffle(struc[:,0])
	rng.shuffle(struc[:,1])  
	households = struc[:,0]
	workplaces = struc[:,1]
	
	# adjacency matrixes for households and workplaces
	H = np.zeros((N, N))
	hmax = np.max(households)
	for h in np.arange(1, hmax+1) :
		members = np.where(households == h)[0]
		for a in members :
			for b in members :
				if a != b :
					H[a,b] = 1
	W = np.zeros((N, N))
	wmax = np.max(workplaces)
	for w in np.arange(1, wmax+1) :
		members = np.where(workplaces == w)[0]
		for a in members :
			for b in members :
				if a != b :
					W[a,b] = 1
	return(H,W)


def neighbours(A, v) :
	""" Returns the neighbours of a node v in a graph of adjacency matrix A """
	adj = A[v,:] 
	nbrs = np.argwhere(adj > 0).flatten()
	return(nbrs)


def composition(state) :
	""" Returns the number of susceptibles, infected and recovered observed in state """
	nbS = np.sum(state == "S")
	nbI = np.sum(state == "I")
	nbR = np.sum(state == "R")
	return(nbS, nbI, nbR)


def typedist(pop, G):
	""" Returns the number of cliques observed in each epidemic state (S,I).

	INPUT
		pop : np.array of length N, indicating the epidemic state (S/I/R) of each node.
		G : graph formed by social structures (households/workplaces) -- adjacency matrix

	OUTPUT
		epidemic_states : dictionnary whose keys represent epidemic states (S,I), and the associated values are the number of occurrences. Only observed epidemic states are taken into account.
	"""
	
	cliques = []
	epidemic_states = {}
	N = len(pop)
	
	# Extracting the cliques from the adjacency matrix
	unknown_group = np.arange(N)
	trying = True
	while np.size(unknown_group) > 0:
		ind = unknown_group[0]
		neighbrs = np.where(G[ind,:] == 1)[0]
		cliques.append(np.append(neighbrs, ind))
		unknown_group = np.setdiff1d(unknown_group, cliques[-1])
		
	# For each clique, compute its current epidemic state and update dictionnary accordingly
	for clque in cliques:
		state = collections.Counter([pop[i] for i in clque])
		key = "(%i,%i)" %(state["S"],state["I"])
		if key not in epidemic_states :
			epidemic_states[key] = 1
		else :
			epidemic_states[key] += 1
			
	return(epidemic_states)



def epidemic_ghw_timeshots(N, H, W, betaG, lambdaH, lambdaW, nu, Tmax, initial_prop, seed) : 
	""" SSA algorithm (Gillespie's algorithm) for the SIR household-workplace model. 
		Initial condition : a fraction of uniformly chosen individuals is infected at random
		Aim : keep track of the proportion of susceptibles, infected and recovered in the population over time.
	
	INPUT
		N : population size
		H : adjacency matrix representing households
		W : adjacency matrix representing workplaces
		betaG : one-to-all infectious contact rate in the general population
		lambdaH : one-to-one infectious contact rate within households
		lambdaW : one-to-one infectious contact rate within workplaces
		nu : removal rate
		Tmax : final time at which the simulation is stopped
		initial_prop : initial proportion of infected
		seed : seed for RNG (or RNG itself)

	OUTPUT
		times, snapshots : np.arrays such that snapshots[j] yields the number of S/I/R in the population at time times[j].
	"""
	
	
	rng = np.random.default_rng(seed)
	
	# Initial condition
	t = 0
	pop = np.repeat("S", N)
	nb_infected = int(initial_prop*N) # nombre d'infetces
	infected = rng.permutation(N)[:nb_infected]
	nb_susceptibles = N - nb_infected # reste de la pop : susceptibles
	susceptibles = np.array([j for j in range(N) if j not in infected])
	pop[infected] = "I"
	
	# Directed edges connecting infected and susceptible individuals i.e. along which infections may occur
	edgeH = np.array([(x, y) for x in infected for y in neighbours(H, x) if pop[y] == "S"], dtype="int,int")
	edgeW = np.array([(x, y) for x in infected for y in neighbours(W, x) if pop[y] == "S"], dtype="int,int")
	nb_edgeH = edgeH.size
	nb_edgeW = edgeW.size
	
	# times : list of points in time at which the current state of the population is saved ('timeshot')
	maxsteps = 1000
	times = np.linspace(0,Tmax, maxsteps)
	snapshot = [composition(pop)]
	nbsteps = 1
	
	# Starting simulation - it is stopped when maximum time is reached, or when there are no infected individuals
	while t < Tmax and nb_infected != 0:
		
		# Jump rates
		tauxG = betaG*nb_susceptibles*nb_infected/N
		tauxH = lambdaH*nb_edgeH
		tauxW = lambdaW*nb_edgeW
		tauxR = nu*nb_infected
		total = tauxG + tauxH + tauxW + tauxR
		
		# Time of the next jump
		u = rng.uniform()
		dt = rng.exponential(1/total)
		t += dt

		# Save population state if the time after the jump exceeds the next timeshot
		if t >= times[nbsteps] :
			trying = True 
			while trying :
				nbsteps += 1
				snapshot.append(composition(pop))
				if nbsteps == maxsteps :
					trying = False
				elif t < times[nbsteps] :
					trying = False
		
		# Updating the population state
		# Option 1 : infection in the general population
		if u < tauxG/total : 
			
			# A susceptible is infected 
			new_infected = rng.choice(susceptibles)
			infected = np.append(infected, new_infected)
			pop[new_infected] = "I"
			susceptibles = susceptibles[susceptibles != new_infected]
			
			# Updating the set of directed edges along which infections may happen 
			edgeH = np.array([arc for arc in edgeH if arc[1] != new_infected], dtype="int,int")
			edgeW = np.array([arc for arc in edgeW if arc[1] != new_infected], dtype="int,int")
			newH = np.array([(new_infected, y) for y in neighbours(H, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeH = np.concatenate((edgeH, newH))
			newW = np.array([(new_infected, y) for y in neighbours(W, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeW = np.concatenate((edgeW, newW))
		
		# Option 2 : within-household infection
		elif u < (tauxG+tauxH)/total : 

			# Directed edge responsible for the infection 
			infectious_edge_id = rng.integers(nb_edgeH)
			new_infected = edgeH[infectious_edge_id][1]
			
			# Contamination of the associated susceptible 
			susceptibles = susceptibles[susceptibles != new_infected]
			infected = np.append(infected, new_infected)
			pop[new_infected] = "I"
			
			# Updating the set of directed edges along which infections may happen
			edgeH = np.array([arc for arc in edgeH if arc[1] != new_infected], dtype="int,int")
			newH = np.array([(new_infected, y) for y in neighbours(H, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeH = np.concatenate((edgeH, newH))
			edgeW = np.array([arc for arc in edgeW if arc[1] != new_infected], dtype="int,int")
			newW = np.array([(new_infected, y) for y in neighbours(W, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeW = np.concatenate((edgeW, newW))
		
		# Option 3 : within-workplace infection
		elif u < (tauxG + tauxH + tauxW)/total : 
			
			# Directed edge responsible for the infection 
			infectious_edge_id = rng.integers(nb_edgeW)
			new_infected = edgeW[infectious_edge_id][1]
			
			# Contamination of the associated susceptible
			susceptibles = susceptibles[susceptibles != new_infected]
			infected = np.append(infected, new_infected)
			pop[new_infected] = "I"
			
			# Updating the set of directed edges along which infections may happen
			edgeH = np.array([arc for arc in edgeH if arc[1] != new_infected], dtype="int,int")
			newH = np.array([(new_infected, y) for y in neighbours(H, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeH = np.concatenate((edgeH, newH)) 
			edgeW = np.array([arc for arc in edgeW if arc[1] != new_infected], dtype="int,int")
			newW = np.array([(new_infected, y) for y in neighbours(W, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeW = np.concatenate((edgeW, newW))
		
		# Option 4 : removal 
		else : # Guerison
		
			# An infected recovers
			removed_id = rng.integers(nb_infected)
			removed = infected[removed_id]
			infected = np.delete(infected, removed_id)
			pop[removed] = "R"
		
			# Updating the set of directed edges along which infections may happen
			edgeH = np.array([arc for arc in edgeH if arc[0] != removed], dtype="int,int")
			edgeW = np.array([arc for arc in edgeW if arc[0] != removed], dtype="int,int")
		
		# Updating descriptive quantities
		nb_susceptibles = susceptibles.size
		nb_infected = infected.size
		nb_edgeH = edgeH.size
		nb_edgeW = edgeW.size

	# Final timeshots if the disease dies out
	if t < Tmax and infected.size == 0 :
		final = composition(pop)
		for j in range(nbsteps, maxsteps) :
			snapshot.append(final)
	
	return(times, snapshot)	




def epidemic_ghw_distribution(N, H, W, betaG, lambdaH, lambdaW, nu, Tmax, final_prop, seed) : 
	""" SSA algorithm for the SIR household-workplace model. 
		Initial condition : 1 infected
		Aim : compute epidemic state distributions for households and workplaces when a given threshold of infected is reached.

		* The code being almost identical to epidemic_ghw_timeshots, only changes are commented. *
	
	INPUT
		N : population size
		H : adjacency matrix representing households
		W : adjacency matrix representing workplaces
		betaG : one-to-all infectious contact rate in the general population
		lambdaH : one-to-one infectious contact rate within households
		lambdaW : one-to-one infectious contact rate within workplaces
		nu : removal rate
		Tmax : final time at which the simulation is stopped
		final_prop : threshold proportion of infected
		seed : seed for RNG (or RNG itself)

	OUTPUT
		prop_inf, prop_susc : proportions of infected and susceptibles when the threshold is reached
		hstate_dist, wstate_dist : dictionnaries describing the epidemic state distributions for households and workplaces, respectively (generated by typedist)

		NB : if the epidemic dies out before reaching the threshold, prop_susc = -1 and the epidemic state distributions are empty dictionnaries.
	"""
	
	rng = np.random.default_rng(seed)
	
	# Initial condition - starting from a single infected
	t = 0
	pop = np.repeat("S", N)
	nb_infected = 1
	infected = rng.permutation(N)[:nb_infected]
	prop_inf = nb_infected/N
	nb_susceptibles = N - nb_infected 
	susceptibles = np.array([j for j in range(N) if j not in infected])
	pop[infected] = "I"
	
	edgeH = np.array([(x, y) for x in infected for y in neighbours(H, x) if pop[y] == "S"], dtype="int,int")
	edgeW = np.array([(x, y) for x in infected for y in neighbours(W, x) if pop[y] == "S"], dtype="int,int")
	nb_edgeH = edgeH.size
	nb_edgeW = edgeW.size
	
	
	while t < Tmax and nb_infected != 0 and prop_inf < final_prop:
		
		tauxG = betaG*nb_susceptibles*nb_infected/N
		tauxH = lambdaH*nb_edgeH
		tauxW = lambdaW*nb_edgeW
		tauxR = nu*nb_infected
		total = tauxG + tauxH + tauxW + tauxR
		u = rng.uniform()
		dt = rng.exponential(1/total)
		t += dt
		
		if u < tauxG/total :
			
			new_infected = rng.choice(susceptibles)
			infected = np.append(infected, new_infected)
			pop[new_infected] = "I"
			susceptibles = susceptibles[susceptibles != new_infected]

			edgeH = np.array([arc for arc in edgeH if arc[1] != new_infected], dtype="int,int")
			edgeW = np.array([arc for arc in edgeW if arc[1] != new_infected], dtype="int,int")
			newH = np.array([(new_infected, y) for y in neighbours(H, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeH = np.concatenate((edgeH, newH))
			newW = np.array([(new_infected, y) for y in neighbours(W, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeW = np.concatenate((edgeW, newW))
		
		elif u < (tauxG+tauxH)/total : 

			infectious_edge_id =rng.integers(nb_edgeH)
			new_infected = edgeH[infectious_edge_id][1]
			susceptibles = susceptibles[susceptibles != new_infected]
			infected = np.append(infected, new_infected)
			pop[new_infected] = "I"
			
			edgeH = np.array([arc for arc in edgeH if arc[1] != new_infected], dtype="int,int")
			newH = np.array([(new_infected, y) for y in neighbours(H, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeH = np.concatenate((edgeH, newH))
			edgeW = np.array([arc for arc in edgeW if arc[1] != new_infected], dtype="int,int")
			newW = np.array([(new_infected, y) for y in neighbours(W, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeW = np.concatenate((edgeW, newW))
		
		elif u < (tauxG + tauxH + tauxW)/total : 
		
			infectious_edge_id = rng.integers(nb_edgeW)
			new_infected = edgeW[infectious_edge_id][1]
			susceptibles = susceptibles[susceptibles != new_infected]
			infected = np.append(infected, new_infected)
			pop[new_infected] = "I"
			
			edgeH = np.array([arc for arc in edgeH if arc[1] != new_infected], dtype="int,int")
			newH = np.array([(new_infected, y) for y in neighbours(H, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeH = np.concatenate((edgeH, newH))
			edgeW = np.array([arc for arc in edgeW if arc[1] != new_infected], dtype="int,int")
			newW = np.array([(new_infected, y) for y in neighbours(W, new_infected) if pop[y] == "S"], dtype="int,int")
			edgeW = np.concatenate((edgeW, newW))
		
		else : 
		
			removed_id = rng.integers(nb_infected)
			removed = infected[removed_id]
			infected = np.delete(infected, removed_id)
			pop[removed] = "R"
		
			edgeH = np.array([arc for arc in edgeH if arc[0] != removed], dtype="int,int")
			edgeW = np.array([arc for arc in edgeW if arc[0] != removed], dtype="int,int")
		
		
		nb_susceptibles = susceptibles.size
		nb_infected = infected.size
		prop_inf = nb_infected/N
		nb_edgeH = edgeH.size
		nb_edgeW = edgeW.size

	
	hstate_dist = {}
	wstate_dist = {}
	prop_susc = -1
	
	# If the desired threshold of infected is reached, compute the epidemic state distributions for households and workplaces
	if prop_inf >= final_prop:
		hstate_dist = typedist(pop, H)
		wstate_dist = typedist(pop, W)
		prop_susc = np.sum(pop == "S")/N
		
	return(prop_inf, prop_susc, hstate_dist, wstate_dist)	




