# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:25:35 2022

@author: anumi
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

expt = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '21-12-20_M2_P1_C3_disc20um'

file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder+'_PY.csv'
data = pd.read_csv(file, sep=';')
bead_dia = 4.503 #in um

xyz_dist = data['D3'] - bead_dia
xy_dist = data['D2'] - bead_dia
dz = data['dz']
#t = (data['T']*1000)/60
t = np.linspace(0,len(data['T']),len(data['T']))
Nan_thresh = 3


outlier = np.where(xyz_dist > Nan_thresh)[0]
xyz_dist[outlier] = np.nan
xy_dist[outlier] = np.nan
dz[outlier] = np.nan

plt.figure(figsize=(20,20))
plt.rcParams.update({'font.size': 22})

ax1 = plt.subplot(311)
plt.plot(t, xyz_dist)
plt.title('3D Distance vs. Time')
#plt.xlim(0,25)

plt.rcParams.update({'font.size': 22})

# share x only
ax2 = plt.subplot(312)
plt.plot(t, xy_dist)
plt.title('2D Distance (XY) vs. Time (mins)')
# make these tick labels invisible

# share x and y
ax3 = plt.subplot(313)
plt.title('Dz vs. Time (mins)')
plt.plot(t, dz)
plt.show()