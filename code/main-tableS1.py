#!/usr/bin/env python3

import ioput

import numpy as np
import matplotlib.pyplot as plt


### TABLE S1 SUPPLEMENTARY MATERIAL

# input files

febcm1 = '../article-data/ctime/ebcm/ctime-ebcm-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-Tmax30-job702612'
febcm2 = '../article-data/ctime/ebcm/ctime-ebcm-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-Tmax30-job702755'
flim = '../article-data/ctime/ebcm/ctime-largepop-bG0p12500-bH1p50000-bW0p00115-nu0p12500-eps0p00500-Tmax30-job702760'

nebcm = 5
nsim = 50


# raw data 

ref1, ebcm1 = ioput.read_timingdata(febcm1, nebcm)
ref2, ebcm2 = ioput.read_timingdata(febcm2, nebcm)
ref = np.append(ref1, ref2)
ebcm = np.append(ebcm1, ebcm2)
ebcm_ratio = ebcm/ref

reflim, lim = ioput.read_timingdata(flim, nsim)
lim_ratio = lim/reflim

# Values of interest 

print("EBCM : average %0.3f, min : %0.2f, max : %0.2f" % (np.mean(ebcm_ratio), np.min(ebcm_ratio), np.max(ebcm_ratio)))
print("Large population limit : average %0.3f, min : %0.2f, max : %0.2f" % (np.mean(lim_ratio), np.min(lim_ratio), np.max(lim_ratio)))





