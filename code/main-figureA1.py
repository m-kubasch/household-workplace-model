#!/usr/bin/env python3

import ioput 
import graphics

import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple


### FIGURE A1

# Input files

filePGfR01p2 = '../article-data/limit/largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-Tmax150.txt'
filePGfR01p4 = '../article-data/limit/largepop-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-Tmax150.txt'
filePGfR01p7 = '../article-data/limit/largepop-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-Tmax150.txt'
filePGfR02p0 = '../article-data/limit/largepop-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-Tmax150.txt'
filePGfR02p5 = '../article-data/limit/largepop-bG0p06000-bH0p20000-bW0p00220-nu0p12500-eps0p00500-Tmax150.txt'

filePGmR01p2 = '../article-data/limit/largepop-bG0p06000-bH0p06000-bW0p00075-nu0p12500-eps0p00500-Tmax150.txt'
filePGmR01p4 = '../article-data/limit/largepop-bG0p07000-bH0p07000-bW0p00080-nu0p12500-eps0p00500-Tmax150.txt'
filePGmR01p7 = '../article-data/limit/largepop-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p00500-Tmax150.txt'
filePGmR02p0 = '../article-data/limit/largepop-bG0p10000-bH0p15000-bW0p00110-nu0p12500-eps0p00500-Tmax150.txt'
filePGmR02p5 = '../article-data/limit/largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-Tmax150.txt'


# Output file

pwd = "../article-data/figures/"
fname = pwd + "fig-runtime-scenarios.pdf"

filePGf = [filePGfR01p2, filePGfR01p4, filePGfR01p7, filePGfR02p0, filePGfR02p5]
filePGm = [filePGmR01p2, filePGmR01p4, filePGmR01p7, filePGmR02p0, filePGmR02p5]

r0 = [1.2, 1.4, 1.7, 2.0, 2.5]

stoptimes = np.zeros((2,5))

# General layout 

fig, axs = plt.subplots(nrows=2, ncols=1, figsize = (9,7))

colorS = dict(zip(r0, ["#5CE684", "#2BDF5F", "#1BB648", "#138534", "#0C5421"]))
colorI = dict(zip(r0, ["#EE6A6A", "#E83737", "#CF1717", "#9C1111", "#690B0B"]))

# Plotting the trajectories for each scenario
stoptimes[0,:] = graphics.plotscenarios(axs[0], filePGf, r0, colorS, colorI)
stoptimes[1,:] = graphics.plotscenarios(axs[1], filePGm, r0, colorS, colorI)

# Decorating the axis
axs[0].set_xticks([],[])
axs[0].set_ylabel('Proportion of individuals', fontsize = 12)
axs[0].set_xlim(0,150)
axs[0].set_title(r"$(p_G,p_H,p_W) = (0.2, 0.4, 0.4)$")

axs[1].set_xlim(0,150)
axs[1].set_xlabel('Time', fontsize = 12)
axs[1].set_ylabel('Proportion of individuals', fontsize = 12)
axs[1].set_title(r"$(p_G,p_H,p_W) = (0.4, 0.4, 0.2)$")

# Adding the legend
legend_elements = [Line2D([0], [0], color=colorS[R0], lw=2, label=r"$R_0: %0.1f$ - susceptibles" % R0) for R0 in r0]
legend_elements += [Line2D([0], [0], color=colorI[R0], lw=2, label=r"$R_0: %0.1f$ - infected" % R0) for R0 in r0]
handles = [(Rectangle((0,0), 0, 0, color=colorS[R0]),Rectangle((0,0), 0, 0, color=colorI[R0])) for R0 in r0]
labels = [r"%0.1f" % R0 for R0 in r0]
leg = axs[0].legend(handles, labels, bbox_to_anchor=(1.01, 1), loc='upper left',handler_map={tuple: HandlerTuple(ndivide=None)}, title = r"S  I    $\mathbf{R_0}$", fontsize = 12, title_fontproperties = {'size' : 12, 'weight' : 'semibold'})
plt.subplots_adjust(hspace=0.12)

# Saving the figure
plt.savefig(fname, bbox_extra_artists=(leg,), bbox_inches='tight')

## Uncomment the following line to see printed the different times at which the trajectories fall below the threshold of 1%, for each scenario.
#print(stoptimes) 

