#!/usr/bin/env python3

import ioput 
import graphics
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

### FIGURE 1 SUPPLEMENTARY MATERIAL

# Input files

# epsilon = 0.001 - Pannel (a)

fsimA = '../article-data/simulations/comparison-limit/r02p5-eps0p001/simulation-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00100-N10000-Tmax60-job692134'
flimA = '../article-data/limit/largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00100-Tmax60.txt'
febcmA = '../article-data/ebcm/ebcm-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00100-Tmax60.txt'

# epsilon = 0.01 - Pannel (b)

fsimB = '../article-data/simulations/comparison-limit/r02p5-eps0p01/simulation-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p01000-N10000-Tmax60-job692160'
flimB = '../article-data/limit/largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p01000-Tmax60.txt'
febcmB = '../article-data/ebcm/ebcm-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p01000-Tmax60.txt'

# epsilon = 0.05 - Pannel (c)

fsimC = '../article-data/simulations/comparison-ebcm/simulation-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p05000-N10000-Tmax60-job692184'
flimC = '../article-data/limit/largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p05000-Tmax60.txt'
febcmC = '../article-data/ebcm/ebcm-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p05000-Tmax60.txt'



# Parameters

nbsim = 50 # Number of repetitions per simulation
simscol = 'yellowgreen'
simicol = 'lightcoral'
limscol = 'green'
limicol = 'red'
ebcmscol = 'blue'
ebcmicol =  'black'
alpha = 0.3 
Tmax = 50 # Time up to which the plot is shown

# Output 
pwd = '../article-data/figures/'
ftitle = pwd + 'fig-ebcm.png'


# Disposition of the subplots
disposition = [(0,0,flimA,febcmA,fsimA), (0,1,flimB,febcmB,fsimB), (1,0,flimC,febcmC,fsimC)]

# Pannel names
dic = {(0,0) : '(a)', (0,1) : '(b)', (1,0) : '(c)', (1,1) : '(d)'}


fig, axs = plt.subplots(nrows=2, ncols=2, figsize = (10,10))

# SUBPLOTS

for x, y, flim, febcm, fsim in disposition :
	
	# Plotting the data
	
	p, simtime, simS, simI = ioput.read_simarray(fsim, nbsim)
	p, t, S, I = ioput.read_reduction(flim)
	pebcm, tebcm, Sebcm, Iebcm = ioput.read_reduction(febcm)
	
	graphics.plotsimu(axs[x,y], simtime, simS, simI, scol = simscol, icol = simicol, alpha = alpha, threshold = 0.005)
	graphics.plotreduction(axs[x,y], t, S, I, scol = limscol, icol = limicol, threshold = 0.005)
	graphics.plotreduction(axs[x,y], tebcm, Sebcm, Iebcm, scol = ebcmscol, icol = ebcmicol, lst = 'dashed', threshold = 0.005)
	
	
	axs[x,y].set_xlim(0, Tmax)
	if x==0 and y == 0 :
		axs[x,y].set_xticks([], []) 
	if (x==1 and y == 0) or (x==0 and y == 1):
		axs[x,y].set_xlabel('Time')
	if y==0 :
		axs[x,y].set_ylabel('Proportion of individuals')
	if y==1 :
		axs[x,y].set_yticks([], [])
		
	if x == 0 :
		if y == 0 :
			eps = 0.001
		if y == 1 :
			eps = 0.01
	if x == 1 :
		eps = 0.05
	if not (x == 1 and y == 1):
		axs[x,y].set_title(r"%s $\varepsilon=$%0.3f" %(dic[(x,y)],eps))
		
		
legend_elements = [	Rectangle((0,0), 0, 0, color='none', label = 'Stochastic simulations'),
	Line2D([0], [0], color=simscol, lw=2, alpha = alpha, label='Proportion of susceptibles'),
	Line2D([0], [0], color=simicol, lw=2, alpha = alpha, label='Proportion of infected'),
	Rectangle((0,0), 0, 0, color='none', label = ''),
	Rectangle((0,0), 0, 0, color='none', label = 'Large population approximation'),
	Line2D([0], [0], color = limscol, lw=2, label='Proportion of susceptibles'),
	Line2D([0], [0], color = limicol, lw=2, label='Proportion of infected'),
	Rectangle((0,0), 0, 0, color='none', label = ''),
	Rectangle((0,0), 0, 0, color='none', label = 'EBCM'),
	Line2D([0], [0], color = ebcmscol, linestyle = 'dashed', lw=2, label='Proportion of susceptibles'),
	Line2D([0], [0], color = ebcmicol, linestyle = 'dashed', lw=2, label='Proportion of infected')]

leg = axs[1,1].legend(handles=legend_elements, mode = 'expand', loc='center left')
leg.get_texts()[0].set_weight('bold')
leg.get_texts()[4].set_weight('bold')
leg.get_texts()[8].set_weight('bold')
axs[1,1].xaxis.set_visible(False)
axs[1,1].yaxis.set_visible(False)
axs[1,1].patch.set_facecolor('none')
for spine in ['top', 'right', 'left', 'bottom']:
	axs[1,1].spines[spine].set_visible(False)
	
	
plt.subplots_adjust(wspace=0.07,hspace=0.08)
plt.savefig(ftitle, dpi=400, bbox_extra_artists=(leg,), bbox_inches='tight')

	
	
	