# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 15:11:54 2021

@author: anumi
"""

import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt


expt = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '21-12-20_M3_P1_C3_disc20um'

out_path = 'D:/Anumita/optoPincher Experiments/'+expt+'/'+folder+'/PTLResults.txt'

file = 'D:/Anumita/optoPincher Experiments/'+expt+'/'+folder+'/test.LOG'
data = pd.read_csv(file, sep=',', skiprows=[0,3,7,8])

col_planeNo =  np.asarray(data[data.columns[0]])
col_time = np.asarray(data[data.columns[1]])
col_plane =  np.asarray(data[data.columns[2]])

col_ind = np.where(col_time == col_time[0])[0]
noOfPlanes = 3
noOfCh = 2
totalFrames = 6000

times = []
planes = []
planeNos = []


# %% Making a text file for the main tracker analysis

for i in col_ind[0:-1]:
    for j in range(1,noOfPlanes+1):
        ind = j
        time = col_time[i+ind]
        plane = col_plane[i+ind][1:]
        planeNo =  col_planeNo[i+ind]
        split_time = time[2:-1].split(':')
        time_sec = 3600*int(split_time[0]) + 60*int(split_time[1]) + float(split_time[2])
        np.asarray(times.append(time_sec))
        np.asarray(planes.append(plane))
        np.asarray(planeNos.append(planeNo))

    if (col_planeNo[i+noOfPlanes+1] == '1') == True:
        for j in range(1,noOfPlanes+1):
            print('option1')
            ind = noOfPlanes+j
            time = col_time[i+ind]
            plane = col_plane[i+ind][1:]
            planeNo =  col_planeNo[i+ind]
            split_time = time[2:-1].split(':')
            time_sec = 3600*int(split_time[0]) + 60*int(split_time[1]) + float(split_time[2])
            np.asarray(times.append(time_sec))
            np.asarray(planes.append(plane))
            np.asarray(planeNos.append(planeNo))
        
    if (col_planeNo[i+noOfCh*noOfPlanes+1] == '1') == True:
        for j in range(1,noOfPlanes+1):
            print('option2')
            ind = noOfCh*noOfPlanes+j
            time = col_time[i+ind]
            plane = col_plane[i+ind][1:]
            planeNo =  col_planeNo[i+ind]
            split_time = time[2:-1].split(':')
            time_sec = 3600*int(split_time[0]) + 60*int(split_time[1]) + float(split_time[2])
            np.asarray(times.append(time_sec))
            np.asarray(planes.append(plane))
            np.asarray(planeNos.append(planeNo))

if len(times) < totalFrames:
    for i in range(totalFrames - len(times)):
        np.asarray(times.append('nan'))
        np.asarray(planes.append('nan'))
        np.asarray(planeNos.append('nan'))
    
    
#Creating a fake magnetic field column
field = [15.0]*len(times)
field = np.asarray(field)

#writing the data in a new txt file
all_data = np.asarray([field, times, field, planes])
np.savetxt(out_path, all_data.T, fmt='%s, %s, %s, %s')

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

