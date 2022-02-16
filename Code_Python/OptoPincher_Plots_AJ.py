# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:25:35 2022

@author: anumi
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#%% Plotting alll three graphs (3D, 2D and Dz)

expt = '20220203_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '22-02-03_M5_P1_C3_disc20um'
date = '22.02.03'

file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder+'_PY.csv'
data = pd.read_csv(file, sep=';')
bead_dia = 4.503 #in um

xyz_dist = data['D3'] - bead_dia
xy_dist = data['D2'] - bead_dia
dz = data['dz']
t = (data['T']*1000)/60
#t = np.linspace(0,len(data['T']),len(data['T']))
Nan_thresh = 3

outlier = np.where(xyz_dist > Nan_thresh)[0]
xyz_dist[outlier] = np.nan
xy_dist[outlier] = np.nan
dz[outlier] = np.nan

fig = plt.figure(figsize=(20,20))
fig.suptitle(folder, fontsize=16)
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

plt.savefig('D:/Anumita/PincherPlots/'+folder+'_DistancevTime.jpg')


#%% Plotting just 3D graphs

expt = '20220203_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '22-02-03_M5_P1_C3_disc20um'
date = '22.02.03'

file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder+'_PY.csv'
data = pd.read_csv(file, sep=';')
bead_dia = 4.503 #in um

xyz_dist = data['D3'] - bead_dia
xy_dist = data['D2'] - bead_dia
dz = data['dz']
t = (data['T']*1000)/60
#t = np.linspace(0,len(data['T']),len(data['T']))
Nan_thresh = 3

outlier = np.where(xyz_dist > Nan_thresh)[0]
xyz_dist[outlier] = np.nan
xy_dist[outlier] = np.nan
dz[outlier] = np.nan

fig = plt.figure(figsize=(20,20))
fig.suptitle(folder, fontsize=16)
plt.rcParams.update({'font.size': 22})
plt.ylabel('Thickness (nm)')
plt.xlabel('Time (mins)')

plt.plot(t, xyz_dist)
plt.title('3D Distance vs. Time')

plt.savefig('D:/Anumita/PincherPlots/'+folder+'_3DistancevTime.jpg')


# %%
expt1 = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
folder1 = '21-12-20_M1_P1_C3_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder1+'_PY.csv'
data1 = pd.read_csv(file, sep=';')
t1 = np.linspace(0,len(data1['T']),len(data1['T']))
xy_dist1 = data1['D2'] - bead_dia

expt2 = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
folder2 = '21-12-20_M2_P1_C3_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder2+'_PY.csv'
data2 =  pd.read_csv(file, sep=';')
t2 = np.linspace(0,len(data2['T']),len(data2['T']))
xy_dist2 = data2['D2'] - bead_dia

Nan_thresh = 3

outlier = np.where(xy_dist2 > Nan_thresh)[0]
xy_dist2[outlier] = np.nan

# expt3 = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
# folder3 = '21-12-20_M3_P1_C2_disc20um'
# file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder3+'_PY.csv'
# data3 = pd.read_csv(file, sep=';')
# t3 = np.linspace(0,len(data3['T']),len(data3['T']))
# xy_dist3 = data3['D2'] - bead_dia

plt.figure(figsize=(20,20))
plt.rcParams.update({'font.size': 35})
plt.title('M1-2_C3_DistancevsTime')

plt.plot(t1, xy_dist1, label="30s")
plt.plot(t2, xy_dist2, label='60s')
# plt.plot(t3, xy_dist3, label='90s')
plt.legend()
plt.show()
plt.savefig('D:/Anumita/PincherPlots/C3_DistancevTime.jpg')


