# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 15:11:54 2021

@author: anumi
"""

import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import os


date = '22.05.09'
expt = '20220509_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '22-05-09_M1_P1_C1_disc20um'

out_path = 'D:/Anumita/MagneticPincherData/Raw/'+date+'/'+folder+'_Field.txt'
rawDir = 'F:/Cortex Experiments/OptoPincher Experiments/'+expt+'/'+date

allFiles = os.listdir(rawDir)

for f in allFiles:
    file = 'F:/Cortex Experiments/OptoPincher Experiments/'+expt+'/'+date+'/'+f+'/test.LOG'
    out_path = 'D:/Anumita/MagneticPincherData/Raw/'+date+'/'+f+'_Field.txt'
    data = pd.read_csv(file, sep=',', skiprows=[0,3,7,8])
    
    col_planeNo =  np.asarray(data[data.columns[0]])
    col_time = np.asarray(data[data.columns[1]])
    col_plane =  np.asarray(data[data.columns[2]])
    
    col_ind = np.where(col_time == col_time[0])[0]
    noOfPlanes = 3
    noOfCh = 2
    totalFrames = 4500
    
    times = []
    planes = []
    planeNos = []
    
    for i in range(len(col_planeNo)):
        if (col_planeNo[i] == '1') == True:
            for j in range(0,noOfPlanes):
                ind = i+j
                time = col_time[ind]
                plane = col_plane[ind][1:]
                planeNo =  col_planeNo[ind]
                split_time = time[2:-1].split(':')
                time_sec = 3600*int(split_time[0]) + 60*int(split_time[1]) + float(split_time[2])
                np.asarray(times.append(time_sec))
                np.asarray(planes.append(plane))
                np.asarray(planeNos.append(planeNo))
        else:
            continue

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
    np.savetxt(out_path, all_data.T, fmt='%s\t%s\t%s\t%s')

# %% Making a text file for the main tracker analysis



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

# %%
import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt


expt = '20220301_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '22-03-01_M2_P1_C7_disc20um'
date = '22.03.01'

out_path = 'D:/Anumita/MagneticPincherData/Raw/'+date+'/'+folder
extDataDir = 'E:/Cortex Experiments/optoPincher Experiments/'+expt+'/'+date+'/'+folder

def CreateFieldFile(self, extDataDir, out_path):
    noOfPlanes = self.Nuplet
    fieldValue = self.MagField
    out_path = out_path+'_Field.txt'
    file = extDataDir+'/test.LOG'
    
    #Reads .LOG file, skips the first couple of rows which are not useful
    data = pd.read_csv(file, sep=',', skiprows=[0,3,7,8])
    
    col_planeNo =  np.asarray(data[data.columns[0]])
    col_time = np.asarray(data[data.columns[1]])
    col_plane =  np.asarray(data[data.columns[2]])
    #totalFrames = len(out_path+'/SplitTriplets')/noOfPlanes
    times = []
    planes = []
    planeNos = []
    
    for i in range(len(col_planeNo)):
        #Finds the row for the first plane of each acquisition and appends
        #the next couple of rows to an empty array which are all written to a .txt file 
        #at the end in the right format
        if (col_planeNo[i] == '1') == True: 
            for j in range(0,noOfPlanes):
                ind = i+j
                time = col_time[ind]
                plane = col_plane[ind][1:]
                planeNo =  col_planeNo[ind]
                split_time = time[2:-1].split(':')
                time_sec = 3600*int(split_time[0]) + 60*int(split_time[1]) + float(split_time[2])
                np.asarray(times.append(time_sec))
                np.asarray(planes.append(plane))
                np.asarray(planeNos.append(planeNo))
        else:
            continue
      
    #Creating a fake magnetic field column
    field = [fieldValue]*len(times)
    field = np.asarray(field)
    
    #writing the data in a new txt file
    all_data = np.asarray([field, times, field, planes])
    np.savetxt(out_path, all_data.T, fmt='%s,%s,%s,%s')
    
    
# %%
date = '22.05.09'
expt = '20220509_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '22-05-09_M1_P1_C1_disc20um'

out_path = 'D:/Anumita/MagneticPincherData/Raw/'+date+'/'+folder
extDataDir = 'E:/Cortex Experiments/optoPincher Experiments/'+expt+'/'+date+'/'+folder

CreateFieldFile(extDataDir, out_path)