#!/usr/bin/env python3

import ioput 
import graphics

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

### FIGURE 3

### Input files

fname = '../article-data/simulations/initial-condition/trajectories/simulation-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p00010-N10000-Tmax150-job692187'
rname = '../article-data/limit/largepop-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p01000-Tmax120.txt'
rnameY0 = '../article-data/limit/largepop-Y0-bG0p08500-bH0p10000-bW0p00100-nu0p12500-Tmax150.txt'

### Output file
outfile = '../article-data/figures/fig-y0.pdf'

### Parameters

# Number of repetitions for simulations
nbsim = 100 

# Colors and transparency for simulations
simscol = 'yellowgreen'
simicol = 'lightcoral'
alpha = 0.3

# Colors for reduced model woth infered initial condition
scol1 = 'darkgreen'
icol1 = 'firebrick'

# Colors for reduced model with uniform initial condition
scol2 = 'green'
icol2 = 'red'

# Threshold for aligning simulations in time
threshold = 0.03

### Reading the data

p, simtime, simS, simI = ioput.read_simarray(fname, nbsim)
N = p['popsize']
I0 = N * p['epsilon']
p, t, S, I = ioput.read_reduction(rname)
py0, ty0, Sy0, Iy0 = ioput.read_reduction_y0(rnameY0)

### Plotting

fig = plt.figure(figsize = (7,6))
ax = fig.add_subplot(111)

plt.axvline(color='black', alpha = 1, linewidth = 0.1)
plt.axhline(y=threshold, color='black', alpha = 1, linewidth = 0.2)
plt.yticks(list(plt.yticks()[0]) + [threshold])


legend_elements = [Rectangle((0,0), 0, 0, color='none', label = "Stochastic simulations (SSA)"),
			Rectangle((0,0), 0, 0, color='none', label = r"$K = %i, I(0)= %i$" %(N, I0)),
			Line2D([0], [0], color='yellowgreen', lw=2, alpha = 0.3, label='Proportion of susceptibles'),
			Line2D([0], [0], color='lightcoral', lw=2, alpha = 0.3, label='Proportion of infected'),
			Rectangle((0,0), 0, 0, color='none', label = ""),
			Rectangle((0,0), 0, 0, color='none', label = ""),
			Rectangle((0,0), 0, 0, color='none', label = r"Large population approximation"),
			#Rectangle((0,0), 0, 0, color='none', label = ""),
			Rectangle((0,0), 0, 0, color='none', label = r"Initial condition (7), $\varepsilon=0.01$"),
			Line2D([0], [0], color=scol2,lw=2, linestyle='dashed', label='Proportion of susceptibles'),
			Line2D([0], [0], color=icol2,lw=2, linestyle='dashed', label='Proportion of infected'),
			Rectangle((0,0), 0, 0, color='none', label = ""),
			Rectangle((0,0), 0, 0, color='none', label = r"Inferred initial condition"),
			Line2D([0], [0], color=scol1,lw=2, label='Proportion of susceptibles'),
			Line2D([0], [0], color=icol1,lw=2, label='Proportion of infected')
]
leg = ax.legend(handles=legend_elements, bbox_to_anchor=(1.01, 1), loc='upper left', ncol = 1)
leg.get_texts()[0].set_weight('bold')
leg.get_texts()[6].set_weight('bold')
leg.get_texts()[7].set_style('italic')
leg.get_texts()[11].set_style('italic')

graphics.plotsimu(ax, simtime, simS, simI, scol = simscol, icol = simicol, threshold = threshold, negtime = True)
graphics.plotreduction(ax, t, S, I, scol = scol2, icol = icol2, threshold = threshold, lst = 'dashed', negtime = True)
graphics.plotreduction(ax, ty0, Sy0, Iy0, scol = scol1, icol = icol1, threshold=threshold, negtime = True)

ax.set_xlim(-50, 100)
ax.set_ylabel('Proportion of individuals')
ax.set_xlabel('Time')

plt.savefig(outfile, dpi=400, bbox_extra_artists=(leg,), bbox_inches='tight')