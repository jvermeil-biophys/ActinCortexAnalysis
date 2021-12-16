# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 15:11:54 2021

@author: anumi
"""

import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt

file = 'D:/Anumita/20211210_100xoil_3t3optorhoa_4.5beads/Chamber1_Cell1_15mT_50msExp_15sActivation_NoScalefactor/test.LOG'
data = pd.read_csv(file, sep=',', skiprows=[0,3])

col_planeNo =  np.asarray(data[data.columns[0]])
col_time = np.asarray(data[data.columns[1]])
col_plane =  np.asarray(data[data.columns[2]])

col_ind = np.where(col_time == col_time[0])
noOfPlanes = 3

times = []
planes = []
planeNos = []

# %% Making a text file for the main tracker analysis

for i in col_ind[0]:
    for j in range(1,noOfPlanes+1):
        time = col_time[i+j]
        plane = col_plane[i+j]
        planeNo =  col_planeNo[i+j]
        np.asarray(times.append(time[2:-1]))
        np.asarray(planes.append(plane))
        np.asarray(planeNos.append(planeNo))
            
#Creating a fake magnetic field column
field = [15.0]*len(times)
field = np.asarray(field)

#writing the data in a new txt file
all_data = np.asarray([field, times, field, planes])
np.savetxt('D:/Anumita/20211210_100xoil_3t3optorhoa_4.5beads/Chamber1_Cell1_15mT_50msExp_15sActivation_NoScalefactor/PTLResults.LOG', all_data.T, fmt='%s, %s, %s, %s', delimiter=' ')

# %% Categorising all the 1st planes to get an idea of the time between two aquisitions 

p1_times = []
frames = []
delta_sec = []
for i in range(len(planeNos)):
    if planeNos[i] == '1':
        p1_times.append(datetime.datetime.strptime(times[i], '%H:%M:%S.%f'))
        frames.append(i)
p1_times = np.asarray(p1_times)
frames = np.asarray(frames)

delta = np.roll(p1_times,-1) - p1_times

for i in delta:
    delta = i.total_seconds()*1000
    delta_sec.append(delta)

plt.figure(figsize=(20,20))
plt.rcParams.update({'font.size': 22})
plt.ylabel('Time interval (ms)')
plt.xlabel('Frame No.')
frame = np.linspace(0, len(delta_sec[0:-2]), len(delta_sec[0:-2]))
plt.plot(frames[0:-2], delta_sec[0:-2])

