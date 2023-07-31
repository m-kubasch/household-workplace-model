#!/usr/bin/env python3

import comptime

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np

### FIGURE A2

# List of considered reproduction numbers
r0 = [1.2, 1.4, 1.7, 2.0, 2.5]

# Initializing the data structures
data = np.zeros((6, 5))
data75 = np.zeros((3, 4))

# Number of repeats per parameter set
nbrep = 99

# Where to save the figure
pwd = '../article-data/figures/'
figname = pwd + 'fig-runtime-ratios.pdf'


# All the input files

### (pG, pH, pW) = (0.2, 0.4, 0.4) -- REPEAT 1

sim1p2 = '../article-data/ctime/pG0p2/R01p2/rep1/ctime-simulation-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-N10000-Tmax130-job695025'
sim1p4 = '../article-data/ctime/pG0p2/R01p4/rep1/ctime-simulation-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-N10000-Tmax130-job695026'
sim1p7= '../article-data/ctime/pG0p2/R01p7/rep1/ctime-simulation-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-N10000-Tmax105-job695027'
sim2p0 = '../article-data/ctime/pG0p2/R02p0/rep1/ctime-simulation-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-N10000-Tmax85-job695028'
sim2p5 = '../article-data/ctime/pG0p2/R02p5/rep1/ctime-simulation-bG0p06000-bH0p20000-bW0p00220-nu0p12500-eps0p00500-N10000-Tmax75-job695029'

red1p2 = '../article-data/ctime/pG0p2/R01p2/rep1/ctime-largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-Tmax130-job694017'
red1p4 = '../article-data/ctime/pG0p2/R01p4/rep1/ctime-largepop-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-Tmax130-job694019'
red1p7 = '../article-data/ctime/pG0p2/R01p7/rep1/ctime-largepop-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-Tmax105-job694021'
red2p0 = '../article-data/ctime/pG0p2/R02p0/rep1/ctime-largepop-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-Tmax85-job694023'
red2p5 = '../article-data/ctime/pG0p2/R02p5/rep1/ctime-largepop-bG0p06000-bH0p20000-bW0p00220-nu0p12500-eps0p00500-Tmax75-job694025'


simlist = [sim1p2, sim1p4, sim1p7, sim2p0, sim2p5]
redlist = [red1p2, red1p4, red1p7, red2p0, red2p5]
r = comptime.timeratio(simlist, redlist, nbrep)
data[0,:] = r

### (pG, pH, pW) = (0.2, 0.4, 0.4) -- REPEAT 2

sim1p2 = '../article-data/ctime/pG0p2/R01p2/rep2/ctime-simulation-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-N10000-Tmax130-job697824'
sim1p4 = '../article-data/ctime/pG0p2/R01p4/rep2/ctime-simulation-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-N10000-Tmax130-job697825'
sim1p7= '../article-data/ctime/pG0p2/R01p7/rep2/ctime-simulation-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-N10000-Tmax105-job697826'
sim2p0 = '../article-data/ctime/pG0p2/R02p0/rep2/ctime-simulation-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-N10000-Tmax85-job697827'
sim2p5 = '../article-data/ctime/pG0p2/R02p5/rep2/ctime-simulation-bG0p06000-bH0p20000-bW0p00220-nu0p12500-eps0p00500-N10000-Tmax75-job697828'

red1p2 = '../article-data/ctime/pG0p2/R01p2/rep2/ctime-largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-Tmax130-job696422'
red1p4 = '../article-data/ctime/pG0p2/R01p4/rep2/ctime-largepop-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-Tmax130-job696424'
red1p7 = '../article-data/ctime/pG0p2/R01p7/rep2/ctime-largepop-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-Tmax105-job696426'
red2p0 = '../article-data/ctime/pG0p2/R02p0/rep2/ctime-largepop-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-Tmax85-job696428'
red2p5 = '../article-data/ctime/pG0p2/R02p5/rep2/ctime-largepop-bG0p06000-bH0p20000-bW0p00220-nu0p12500-eps0p00500-Tmax75-job696430'


simlist = [sim1p2, sim1p4, sim1p7, sim2p0, sim2p5]
redlist = [red1p2, red1p4, red1p7, red2p0, red2p5]
r = comptime.timeratio(simlist, redlist, nbrep)
data[1,:] = r


### (pG, pH, pW) = (0.2, 0.4, 0.4) -- REPEAT 3

sim1p2 = '../article-data/ctime/pG0p2/R01p2/rep3/ctime-simulation-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-N10000-Tmax130-job698827'
sim1p4 = '../article-data/ctime/pG0p2/R01p4/rep3/ctime-simulation-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-N10000-Tmax130-job698837'
sim1p7= '../article-data/ctime/pG0p2/R01p7/rep3/ctime-simulation-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-N10000-Tmax105-job698845'
sim2p0 = '../article-data/ctime/pG0p2/R02p0/rep3/ctime-simulation-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-N10000-Tmax85-job698853'
sim2p5 = '../article-data/ctime/pG0p2/R02p5/rep3/ctime-simulation-bG0p06000-bH0p20000-bW0p00220-nu0p12500-eps0p00500-N10000-Tmax75-job698857'

red1p2 = '../article-data/ctime/pG0p2/R01p2/rep3/ctime-largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-Tmax130-job698828'
red1p4 = '../article-data/ctime/pG0p2/R01p4/rep3/ctime-largepop-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-Tmax130-job698830'
red1p7 = '../article-data/ctime/pG0p2/R01p7/rep3/ctime-largepop-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-Tmax105-job698834'
red2p0 = '../article-data/ctime/pG0p2/R02p0/rep3/ctime-largepop-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-Tmax85-job698838'
red2p5 = '../article-data/ctime/pG0p2/R02p5/rep3/ctime-largepop-bG0p06000-bH0p20000-bW0p00220-nu0p12500-eps0p00500-Tmax75-job698842'


simlist = [sim1p2, sim1p4, sim1p7, sim2p0, sim2p5]
redlist = [red1p2, red1p4, red1p7, red2p0, red2p5]
r = comptime.timeratio(simlist, redlist, nbrep)
data[2,:] = r



### (pG, pH, pW) = (0.2, 0.4, 0.4), Tmax = 75 -- REPEAT 1
sim1p2 = '../article-data/ctime/pG0p2/R01p2/rep1/ctime-simulation-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-N10000-Tmax75-job698833'
sim1p4 = '../article-data/ctime/pG0p2/R01p4/rep1/ctime-simulation-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-N10000-Tmax75-job698841'
sim1p7 = '../article-data/ctime/pG0p2/R01p7/rep1/ctime-simulation-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-N10000-Tmax105-job695027'
sim2p0 = '../article-data/ctime/pG0p2/R02p0/rep1/ctime-simulation-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-N10000-Tmax75-job698855'

red1p2 = '../article-data/ctime/pG0p2/R01p2/rep1/ctime-largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-Tmax75-job694018'
red1p4 = '../article-data/ctime/pG0p2/R01p4/rep1/ctime-largepop-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-Tmax75-job694020'
red1p7 = '../article-data/ctime/pG0p2/R01p7/rep1/ctime-largepop-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-Tmax75-job694022'
red2p0 = '../article-data/ctime/pG0p2/R02p0/rep1/ctime-largepop-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-Tmax75-job694024'


simlist = [sim1p2, sim1p4, sim1p7, sim2p0]
redlist = [red1p2, red1p4, red1p7, red2p0]
data75[0,:] = comptime.timeratio(simlist, redlist, nbrep)


### (pG, pH, pW) = (0.2, 0.4, 0.4), Tmax = 75 -- REPEAT 2
sim1p2 = '../article-data/ctime/pG0p2/R01p2/rep2/ctime-simulation-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-N10000-Tmax75-job698835'
sim1p4 = '../article-data/ctime/pG0p2/R01p4/rep2/ctime-simulation-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-N10000-Tmax75-job698843'
sim1p7 = '../article-data/ctime/pG0p2/R01p7/rep2/ctime-simulation-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-N10000-Tmax75-job698851'
sim2p0 = '../article-data/ctime/pG0p2/R02p0/rep2/ctime-simulation-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-N10000-Tmax75-job698856'

red1p2 = '../article-data/ctime/pG0p2/R01p2/rep2/ctime-largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-Tmax75-job696423'
red1p4 = '../article-data/ctime/pG0p2/R01p4/rep2/ctime-largepop-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-Tmax75-job696425'
red1p7 = '../article-data/ctime/pG0p2/R01p7/rep2/ctime-largepop-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-Tmax75-job696427'
red2p0 = '../article-data/ctime/pG0p2/R02p0/rep2/ctime-largepop-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-Tmax75-job696429'

simlist = [sim1p2, sim1p4, sim1p7, sim2p0]
redlist = [red1p2, red1p4, red1p7, red2p0]
data75[1,:] = comptime.timeratio(simlist, redlist, nbrep)

### (pG, pH, pW) = (0.2, 0.4, 0.4), Tmax = 75 -- REPEAT 3
sim1p2 = '../article-data/ctime/pG0p2/R01p2/rep3/ctime-simulation-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-N10000-Tmax75-job698831'
sim1p4 = '../article-data/ctime/pG0p2/R01p4/rep3/ctime-simulation-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-N10000-Tmax75-job698839'
sim1p7 = '../article-data/ctime/pG0p2/R01p7/rep3/ctime-simulation-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-N10000-Tmax75-job698847'
sim2p0 = '../article-data/ctime/pG0p2/R02p0/rep3/ctime-simulation-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-N10000-Tmax75-job698854'

red1p2 = '../article-data/ctime/pG0p2/R01p2/rep3/ctime-largepop-bG0p03000-bH0p05000-bW0p00150-nu0p12500-eps0p00500-Tmax75-job698829'
red1p4 = '../article-data/ctime/pG0p2/R01p4/rep3/ctime-largepop-bG0p03500-bH0p07000-bW0p00160-nu0p12500-eps0p00500-Tmax75-job698832'
red1p7 = '../article-data/ctime/pG0p2/R01p7/rep3/ctime-largepop-bG0p04500-bH0p09000-bW0p00180-nu0p12500-eps0p00500-Tmax75-job698836'
red2p0 = '../article-data/ctime/pG0p2/R02p0/rep3/ctime-largepop-bG0p05000-bH0p15000-bW0p00200-nu0p12500-eps0p00500-Tmax75-job698840'


simlist = [sim1p2, sim1p4, sim1p7, sim2p0]
redlist = [red1p2, red1p4, red1p7, red2p0]
data75[2,:] = comptime.timeratio(simlist, redlist, nbrep)



### (pG, pH, pW) = (0.4, 0.4, 0.2) -- REPEAT 1

sim1p2 = '../article-data/ctime/pG0p4/R01p2/rep1/ctime-simulation-bG0p06000-bH0p06000-bW0p00075-nu0p12500-eps0p00500-N10000-Tmax145-job695030'
sim1p4 = '../article-data/ctime/pG0p4/R01p4/rep1/ctime-simulation-bG0p07000-bH0p07000-bW0p00080-nu0p12500-eps0p00500-N10000-Tmax130-job695031'
sim1p7= '../article-data/ctime/pG0p4/R01p7/rep1/ctime-simulation-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p00500-N10000-Tmax95-job695032'
sim2p0 = '../article-data/ctime/pG0p4/R02p0/rep1/ctime-simulation-bG0p10000-bH0p15000-bW0p00110-nu0p12500-eps0p00500-N10000-Tmax80-job695033'
sim2p5 = '../article-data/ctime/pG0p4/R02p5/rep1/ctime-simulation-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-N10000-Tmax55-job695034'

red1p2 = '../article-data/ctime/pG0p4/R01p2/rep1/ctime-largepop-bG0p06000-bH0p06000-bW0p00075-nu0p12500-eps0p00500-Tmax145-job694026'
red1p4 = '../article-data/ctime/pG0p4/R01p4/rep1/ctime-largepop-bG0p07000-bH0p07000-bW0p00080-nu0p12500-eps0p00500-Tmax130-job694027'
red1p7 = '../article-data/ctime/pG0p4/R01p7/rep1/ctime-largepop-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p00500-Tmax95-job694028'
red2p0 = '../article-data/ctime/pG0p4/R02p0/rep1/ctime-largepop-bG0p10000-bH0p15000-bW0p00110-nu0p12500-eps0p00500-Tmax80-job694029'
red2p5 = '../article-data/ctime/pG0p4/R02p5/rep1/ctime-largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-Tmax55-job694030'


simlist = [sim1p2, sim1p4, sim1p7, sim2p0, sim2p5]
redlist = [red1p2, red1p4, red1p7, red2p0, red2p5]
r = comptime.timeratio(simlist, redlist, nbrep)
data[3,:] = r

### (pG, pH, pW) = (0.4, 0.4, 0.2) -- REPEAT 2

sim1p2 = '../article-data/ctime/pG0p4/R01p2/rep2/ctime-simulation-bG0p06000-bH0p06000-bW0p00075-nu0p12500-eps0p00500-N10000-Tmax145-job697829'
sim1p4 = '../article-data/ctime/pG0p4/R01p4/rep2/ctime-simulation-bG0p07000-bH0p07000-bW0p00080-nu0p12500-eps0p00500-N10000-Tmax130-job697830'
sim1p7= '../article-data/ctime/pG0p4/R01p7/rep2/ctime-simulation-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p00500-N10000-Tmax95-job697831'
sim2p0 = '../article-data/ctime/pG0p4/R02p0/rep2/ctime-simulation-bG0p10000-bH0p15000-bW0p00110-nu0p12500-eps0p00500-N10000-Tmax80-job697832'
sim2p5 = '../article-data/ctime/pG0p4/R02p5/rep2/ctime-simulation-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-N10000-Tmax55-job697833'

red1p2 = '../article-data/ctime/pG0p4/R01p2/rep2/ctime-largepop-bG0p06000-bH0p06000-bW0p00075-nu0p12500-eps0p00500-Tmax145-job696431'
red1p4 = '../article-data/ctime/pG0p4/R01p4/rep2/ctime-largepop-bG0p07000-bH0p07000-bW0p00080-nu0p12500-eps0p00500-Tmax130-job696432'
red1p7 = '../article-data/ctime/pG0p4/R01p7/rep2/ctime-largepop-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p00500-Tmax95-job696433'
red2p0 = '../article-data/ctime/pG0p4/R02p0/rep2/ctime-largepop-bG0p10000-bH0p15000-bW0p00110-nu0p12500-eps0p00500-Tmax80-job696434'
red2p5 = '../article-data/ctime/pG0p4/R02p5/rep2/ctime-largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-Tmax55-job696435'


simlist = [sim1p2, sim1p4, sim1p7, sim2p0, sim2p5]
redlist = [red1p2, red1p4, red1p7, red2p0, red2p5]
r = comptime.timeratio(simlist, redlist, nbrep)
data[4,:] = r


### (pG, pH, pW) = (0.4, 0.4, 0.2) -- REPEAT 3

sim1p2 = '../article-data/ctime/pG0p4/R01p2/rep3/ctime-simulation-bG0p06000-bH0p06000-bW0p00075-nu0p12500-eps0p00500-N10000-Tmax145-job698858'
sim1p4 = '../article-data/ctime/pG0p4/R01p4/rep3/ctime-simulation-bG0p07000-bH0p07000-bW0p00080-nu0p12500-eps0p00500-N10000-Tmax130-job698859'
sim1p7= '../article-data/ctime/pG0p4/R01p7/rep3/ctime-simulation-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p00500-N10000-Tmax95-job698860'
sim2p0 = '../article-data/ctime/pG0p4/R02p0/rep3/ctime-simulation-bG0p10000-bH0p15000-bW0p00110-nu0p12500-eps0p00500-N10000-Tmax80-job698861'
sim2p5 = '../article-data/ctime/pG0p4/R02p5/rep3/ctime-simulation-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-N10000-Tmax55-job698862'

red1p2 = '../article-data/ctime/pG0p4/R01p2/rep3/ctime-largepop-bG0p06000-bH0p06000-bW0p00075-nu0p12500-eps0p00500-Tmax145-job698844'
red1p4 = '../article-data/ctime/pG0p4/R01p4/rep3/ctime-largepop-bG0p07000-bH0p07000-bW0p00080-nu0p12500-eps0p00500-Tmax130-job698846'
red1p7 = '../article-data/ctime/pG0p4/R01p7/rep3/ctime-largepop-bG0p08500-bH0p10000-bW0p00100-nu0p12500-eps0p00500-Tmax95-job698848'
red2p0 = '../article-data/ctime/pG0p4/R02p0/rep3/ctime-largepop-bG0p10000-bH0p15000-bW0p00110-nu0p12500-eps0p00500-Tmax80-job698850'
red2p5 = '../article-data/ctime/pG0p4/R02p5/rep3/ctime-largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-Tmax55-job698852'

## WAITING FOR DATA
simlist = [sim1p2, sim1p4, sim1p7, sim2p0, sim2p5]
redlist = [red1p2, red1p4, red1p7, red2p0, red2p5]
r = comptime.timeratio(simlist, redlist, nbrep)
data[5,:] = r




##### PLOTTING THE FIGURE

fig, ax = plt.subplots(figsize = (7,5))

plt.scatter(np.array([r0,r0,r0]), data[3:,:], s = 40, color = 'blue', alpha = 0.5, marker = 'P', label = r"$(p_G,p_H,p_W) = (0.4, 0.4, 0.2)$")
plt.scatter(np.array([r0,r0,r0]), data[:3,:], s = 40, color = 'red', alpha = 0.5, label = r"$(p_G,p_H,p_W) = (0.2, 0.4, 0.4)$")
plt.scatter(np.array([r0[:-1],r0[:-1],r0[:-1]]), data75, s = 40, color = 'darkred', alpha = 0.5, marker = 'd', label = r"$(p_G,p_H,p_W) = (0.2, 0.4, 0.4)$, $T=75$")
plt.hlines(1, 1, 2.7, 'black', linewidth = 0.5, linestyle = 'dashed') 
plt.legend()

plt.yscale('log')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
plt.xlabel(r"$R_0$", fontsize = 12) 
plt.xticks(r0)
plt.ylabel("Ratio of average normalised runtimes\nReduced model over SSA (log scale)", fontsize = 12)
plt.savefig(figname, bbox_inches='tight')