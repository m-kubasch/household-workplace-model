#!/usr/bin/env python3

import ioput 
import graphics
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

### FIGURE 2


# Input files

# R0 2p5, eps 0p001 - Pannel (a)
fnameA = '../article-data/simulations/comparison-limit/r02p5-eps0p001/simulation-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00100-N10000-Tmax60-job692134'
rnameA = '../article-data/limit/largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00100-Tmax60.txt'


# R0 1p2, eps 0p001 - Pannel (b)
fnameB = '../article-data/simulations/comparison-limit/r01p2-eps0p001/simulation-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00100-N10000-Tmax400-job692186'
rnameB = '../article-data/limit/largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00100-Tmax400.txt'

#R0 2p5, eps 0p01 - Pannel (c)
fnameC = '../article-data/simulations/comparison-limit/r02p5-eps0p01/simulation-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p01000-N10000-Tmax60-job692160'
rnameC = '../article-data/limit/largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p01000-Tmax60.txt'

#R0 1p2, eps 0p01 - Pannel (d)
fnameD = '../article-data/simulations/comparison-limit/r01p2-eps0p01/simulation-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p01000-N10000-Tmax400-job692185'
rnameD = '../article-data/limit/largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p01000-Tmax250.txt'

# Parameters

nbsim = 50 # Number of repetitions per simulation
simscol = 'yellowgreen'
simicol = 'lightcoral'
redscol = 'green'
redicol = 'red'
alpha = 0.3


# Output 
pwd = '../article-data/figures/'
ftitle = pwd + 'fig-comparison-limit.png'

# Disposition of the subplots

disposition = [(0,0,fnameA,rnameA),(0,1,fnameB,rnameB),(1,0,fnameC,rnameC),(1,1,fnameD,rnameD)]

# pannel names
dic = {(0,0) : '(a)', (0,1) : '(b)', (1,0) : '(c)', (1,1) : '(d)'}


fig, axs = plt.subplots(nrows=2, ncols=2, figsize = (10,10))

# SUBPLOTS

for x, y, fname, rname in disposition :
	
	# Plotting the data
	
	p, simtime, simS, simI = ioput.read_simarray(fname, nbsim)
	p, t, S, I = ioput.read_reduction(rname)
	
	graphics.plotsimu(axs[x,y], simtime, simS, simI, scol = simscol, icol = simicol, threshold = 0.005)
	graphics.plotreduction(axs[x,y], t, S, I, scol = redscol, icol = redicol, threshold = 0.005)

	
	if y==0 :
		Tmax=50
	if y==1 :
		Tmax=200
	
	axs[x,y].set_xlim(0, Tmax)
	if x==0 :
		axs[x,y].set_xticks([], []) 
	if x==1 :
		axs[x,y].set_xlabel('Time')
	if y==0 :
		axs[x,y].set_ylabel('Proportion of individuals')
	if y==1 :
		axs[x,y].set_yticks([], [])
		
	if x==1 and y==0 :
		legend_elements = [Rectangle((0,0), 0, 0, color='none', label = 'Stochastic simulations'),
			Line2D([0], [0], color=simscol, lw=2, alpha = alpha, label='Proportion of susceptibles'),
			Line2D([0], [0], color=simicol, lw=2, alpha = alpha, label='Proportion of infected'),
			Rectangle((0,0), 0, 0, color='none', label = 'Large population approximation'),
			Line2D([0], [0], color=redscol,lw=2, label='Proportion of susceptibles'),
			Line2D([0], [0], color=redicol,lw=2, label='Proportion of infected')]

		leg = axs[x,y].legend(handles=legend_elements, bbox_to_anchor=(0.35, -0.12), loc='upper left', ncol = 2)
		leg.get_texts()[0].set_weight('bold')
		leg.get_texts()[3].set_weight('bold')
		
	if x == 0 :
		eps = 0.001
	if x == 1 :
		eps = 0.01
	if y == 0 :
		r0 = 2.5
	if y == 1 :
		r0 = 1.2
		
	axs[x,y].set_title(r"%s $R_0=$%0.2f, $\varepsilon=$%0.3f" %(dic[(x,y)], r0,eps))
	
plt.subplots_adjust(wspace=0.07,hspace=0.08)
plt.savefig(ftitle, dpi=400, bbox_extra_artists=(leg,), bbox_inches='tight')

