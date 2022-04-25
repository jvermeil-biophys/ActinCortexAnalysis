# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:25:35 2022

@author: anumi
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os


#%% Global constants

bead_dia = 4.503


#%% Utility functions

def getExperimentalConditions(experimentalDataDir, save = False, sep = ';'):
    """"
    Import the table with all the conditions in a clean way.
    It is a tedious function to read because it's doing a boring job:
    Converting strings into numbers when possible
    Converting commas into dots to correct for the French decimal notation
    Converting semicolon separated values into lists when needed
    Etc
    """
    #### 0. Import the table
    experimentalDataFile = 'ExperimentalConditions.csv'
    experimentalDataFilePath = os.path.join(experimentalDataDir, experimentalDataFile)
    expConditionsDF = pd.read_csv(experimentalDataFilePath, sep=sep, header=0)
    print(BLUE + 'Importing Experimental Conditions' + NORMAL)
    print(BLUE + 'Extracted a table with ' + str(expConditionsDF.shape[0]) + ' lines and ' + str(expConditionsDF.shape[1]) + ' columns' + NORMAL)
    #### 1. Clean the table
    
    #### 1.1 Remove useless columns
    for c in expConditionsDF.columns:
        if 'Unnamed' in c:
            expConditionsDF = expConditionsDF.drop([c], axis=1)
        if '.1' in c:
            expConditionsDF = expConditionsDF.drop([c], axis=1)
    expConditionsDF = expConditionsDF.convert_dtypes()

    #### 1.2 Convert commas into dots
    listTextColumns = []
    for col in expConditionsDF.columns:
        try:
            if expConditionsDF[col].dtype == 'string':
                listTextColumns.append(col)
        except:
            pass
    expConditionsDF[listTextColumns] = expConditionsDF[listTextColumns].apply(lambda x: x.str.replace(',','.'))

    #### 1.3 Format 'scale'
    expConditionsDF['scale pixel per um'] = expConditionsDF['scale pixel per um'].astype(float)
    
    #### 1.4 Format 'optical index correction'
    try: # In case the format is 'n1/n2'
        expConditionsDF['optical index correction'] = \
                  expConditionsDF['optical index correction'].apply(lambda x: x.split('/')[0]).astype(float) \
                / expConditionsDF['optical index correction'].apply(lambda x: x.split('/')[1]).astype(float)
        print(ORANGE + 'optical index correction : format changed' + NORMAL)
    except:
        pass
    
    #### 1.5 Format 'magnetic field correction'
    expConditionsDF['magnetic field correction'] = expConditionsDF['magnetic field correction'].astype(float)
    
    #### 1.6 Format 'with fluo images'
    expConditionsDF['with fluo images'] = expConditionsDF['with fluo images'].astype(bool)

    # #### 1.7 Format 'ramp field'
    # try:
    #     print(ORANGE + 'ramp field : converted to list successfully' + NORMAL)
    #     expConditionsDF['ramp field'] = \
    #     expConditionsDF['ramp field'].apply(lambda x: [x.split(';')[0], x.split(';')[1]] if not pd.isnull(x) else [])
    # except:
    #     pass

    #### 1.8 Format 'date'
    dateExemple = expConditionsDF.loc[expConditionsDF.index[1],'date']
    if re.match(dateFormatExcel, dateExemple):
        print(ORANGE + 'dates : format corrected' + NORMAL)
        expConditionsDF.loc[:,'date'] = expConditionsDF.loc[:,'date'].apply(lambda x: x.split('/')[0] + '-' + x.split('/')[1] + '-' + x.split('/')[2][2:])        
    elif re.match(dateFormatExcel2, dateExemple):
        print(ORANGE + 'dates : format corrected' + NORMAL)
        expConditionsDF.loc[:,'date'] = expConditionsDF.loc[:,'date'].apply(lambda x: x.split('-')[0] + '-' + x.split('-')[1] + '-' + x.split('-')[2][2:])  
        
    #### 1.9 Format activation fields
    expConditionsDF['first activation'] = expConditionsDF['first activation'].astype(np.float)
    expConditionsDF['activation frequency'] = expConditionsDF['activation frequency'].astype(np.float)


    #### 2. Save the table, if required
    if save:
        saveName = 'ExperimentalConditions.csv'
        savePath = os.path.join(experimentalDataDir, saveName)
        expConditionsDF.to_csv(savePath, sep=';')

    #### 3. Generate additionnal field that won't be saved
    
    def str2int(s):
        try:
            x = int(s)
        except:
            x = np.nan
        return(x)
    
    def str2float(s):
        try:
            x = float(s)
        except:
            x = np.nan
        return(x)
    
    #### 3.1 Make 'manipID'
    expConditionsDF['manipID'] = expConditionsDF['date'] + '_' + expConditionsDF['manip']

def findInfosInFileName(f, infoType):
    """
    Return a given type of info from a file name.
    Inputs : f (str), the file name.
             infoType (str), the type of info wanted.
             infoType can be equal to : 
             * 'M', 'P', 'C' -> will return the number of manip (M), well (P), or cell (C) in a cellID.
             ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'C', the function will return 8.
             * 'manipID'     -> will return the full manip ID.
             ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'manipID', the function will return '21-01-18_M2'.
             * 'cellID'     -> will return the full cell ID.
             ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'cellID', the function will return '21-01-18_M2_P1_C8'.
    """
    if infoType in ['M', 'P', 'C']:
        acceptedChar = [str(i) for i in range(10)] + ['.', '-']
        string = '_' + infoType
        iStart = re.search(string, f).end()
        i = iStart
        infoString = '' + f[i]
        while f[i+1] in acceptedChar and i < len(f)-1:
            i += 1
            infoString += f[i]
            
    elif infoType == 'date':
        datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
        date = f[datePos.start():datePos.end()]
        infoString = date
    
    elif infoType == 'manipID':
        datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
        date = f[datePos.start():datePos.end()]
        manip = 'M' + findInfosInFileName(f, 'M')
        infoString = date + '_' + manip
        
    elif infoType == 'cellID':
        datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
        date = f[datePos.start():datePos.end()]
        infoString = date + '_' + 'M' + findInfosInFileName(f, 'M') + \
                            '_' + 'P' + findInfosInFileName(f, 'P') + \
                            '_' + 'C' + findInfosInFileName(f, 'C')
                            
    elif infoType == 'substrate':
        try:
            pos = re.search(r"disc[\d]*um", f)
            infoString = f[pos.start():pos.end()]
        except:
            infoString = ''
                 
                             
    return(infoString)

def isFileOfInterest(f, manips, wells, cells):
    """
    Determine if a file f correspond to the given criteria.
    More precisely, return a boolean saying if the manip, well and cell number are in the given range.
    f is a file name. Each of the fields 'manips', 'wells', 'cells' can be either a number, a list of numbers, or 'all'.
    Example : if f = '21-01-18_M2_P1_C8.tif'
    * manips = 'all', wells = 'all', cells = 'all' -> the function return True.
    * manips = 1, wells = 'all', cells = 'all' -> the function return False.
    * manips = [1, 2], wells = 'all', cells = 'all' -> the function return True.
    * manips = [1, 2], wells = 2, cells = 'all' -> the function return False.
    * manips = [1, 2], wells = 1, cells = [5, 6, 7, 8] -> the function return True.
    Note : if manips = 'all', the code will consider that wells = 'all', cells = 'all'.
           if wells = 'all', the code will consider that cells = 'all'.
           This means you can add filters only in this order : manips > wells > cells.
    """
    test = False
    if f.endswith(".tif"):
        if manips == 'all':
            test = True
        else:
            try:
                manips_str = [str(i) for i in manips]
            except:
                manips_str = [str(manips)]
            infoM = findInfosInFileName(f, 'M')
            if infoM in manips_str:
                if wells == 'all':
                    test = True
                else:
                    try:
                        wells_str = [str(i) for i in wells]
                    except:
                        wells_str = [str(wells)]
                    infoP = findInfosInFileName(f, 'P')
                    if infoP in wells_str:
                        if cells == 'all':
                            test = True
                        else:
                            try:
                                cells_str = [str(i) for i in cells]
                            except:
                                cells_str = [str(cells)]
                            infoC = findInfosInFileName(f, 'C')
                            if infoC in cells_str:
                                test = True
    return(test)


def getFieldFile(f):
    return
    
####Plot median thickness before and after activation (box plots)
def cfComparison(f):
    return

####Plot only 3D distance vs Time


####Plot 3D distance, 2D distance and Dz vs. Time


#%% Tests

path = 'D:/Anumita/MagneticPincherData/Raw/22.03.31/'
expt = '22-03-31_M6_P1_C2_disc20um_L40'
file = path+expt+'_Field.txt'

def findFirstActivation(f):
    allData = pd.read_csv(file, sep='\t') #Read the _Field.txt file
    maxZidx = (allData[allData.columns[3]].argmax()) #Finding the index of the max Z
    maxZ = allData[allData.columns[3]][maxZidx] #To check if the value is correct
    return maxZidx, maxZ


maxZ = findFirstActivation(file)

#%% Plotting all three graphs (3D, 2D and Dz)

expt = '20220331_100xoil_3t3optorhoa_4.5beads_15mT_Mechanics'
folder = '22-03-31_M7_P1_C2_disc20um_L40'
date = '22.03.31'

file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder+'_PY.csv'
data = pd.read_csv(file, sep=';')
 #in um

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

plt.style.use('dark_background')

fig = plt.figure(figsize=(20,20))
fig.suptitle(folder, fontsize=16)
plt.rcParams.update({'font.size': 25})

ax1 = plt.subplot(311)
plt.plot(t, xyz_dist)
plt.axvline(x = 5, color = 'r', label = 'Activation begins')
plt.axvline(x = 8, color = 'r', label = 'Activation begins')

plt.title('3D Distance vs. Time')
#plt.xlim(0,25)


# share x only
ax2 = plt.subplot(312)
plt.plot(t, xy_dist)
plt.axvline(x = 5, color = 'r', label = 'Activation begins')
plt.axvline(x = 8, color = 'r', label = 'Activation begins')

plt.title('2D Distance (XY) vs. Time (mins)')
# make these tick labels invisible

# share x and y
ax3 = plt.subplot(313)
plt.axvline(x = 5, color = 'r', label = 'Activation begins')
plt.axvline(x = 8, color = 'r', label = 'Activation begins')

plt.title('Dz vs. Time (mins)')
plt.plot(t, dz)

plt.show()

plt.savefig('D:/Anumita/PincherPlots/'+folder+'_DistancevTime.jpg')


#%% Plotting just 3D graphs

expt = '20220322_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '22-03-22_M2_P1_C5_disc20um'
date = '22.03.22'

file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder+'_PY.csv'
data = pd.read_csv(file, sep=';')
bead_dia = 4.503 #in um

xyz_dist = data['D3'] - bead_dia
xy_dist = data['D2'] - bead_dia
dz = data['dz']
t = (data['T']*1000)/60
#t = np.linspace(0,len(data['T']),len(data['T']))
# Nan_thresh = 3

# outlier = np.where(xyz_dist > Nan_thresh)[0]
# xyz_dist[outlier] = np.nan
# xy_dist[outlier] = np.nan
# dz[outlier] = np.nan

plt.style.use('dark_background')

fig = plt.figure(figsize=(30,10))
fig.suptitle(folder, fontsize=16)
plt.rcParams.update({'font.size': 40})
plt.axvline(x = 5, color = 'r', label = 'Activation begins')
plt.xlim(0.5,8)

plt.ylabel('Thickness (um)')
plt.xlabel('Time (mins)')

plt.plot(t, xyz_dist)
plt.title('3D Distance vs. Time')

plt.savefig('D:/Anumita/PincherPlots/'+folder+'_3DistancevTime.jpg')


# %%
expt1 = '20220322_100xoil_3t3optorhoa_4.5beads_15mT'
folder1 = '22-03-22_M4_P1_C5_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder1+'_PY.csv'
data1 = pd.read_csv(file, sep=';')
#t1 = np.linspace(0,len(data1['T']),len(data1['T']))
t1 = (data1['T']*1000/60)
xy_dist1 = data1['D3'] - bead_dia

expt2 = '20220322_100xoil_3t3optorhoa_4.5beads_15mT'
folder2 = '22-03-22_M3_P1_C5_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder2+'_PY.csv'
data2 =  pd.read_csv(file, sep=';')
# t2 = np.linspace(0,len(data2['T']),len(data2['T']))
t2 = (data2['T']*1000/60)
xy_dist2 = data2['D3'] - bead_dia

Nan_thresh1 = 3
Nan_thresh2 = 3

outlier = np.where(xy_dist2 > Nan_thresh2)[0]
xy_dist2[outlier] = np.nan

outlier = np.where(xy_dist1 > Nan_thresh1)[0]
xy_dist1[outlier] = np.nan

# expt3 = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
# folder3 = '21-12-20_M3_P1_C2_disc20um'
# file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder3+'_PY.csv'
# data3 = pd.read_csv(file, sep=';')
# t3 = np.linspace(0,len(data3['T']),len(data3['T']))
# xy_dist3 = data3['D2'] - bead_dia

plt.style.use('dark_background')


fig= plt.figure(figsize=(30,10))
fig.suptitle(folder, fontsize=16)

# right_side = fig.spines["right"]
# right_side.set_visible(False)
# top_side = fig.spines["top"]
# top_side.set_visible(False)


plt.rcParams.update({'font.size':35})
plt.title('3D Distance vs Time')
plt.ylabel('Thickness (um)')
plt.xlabel('Time (secs)')

plt.plot(t1, xy_dist1, label="Activation away from beads", color = 'orange')
plt.plot(t2, xy_dist2, label='Activation at beads', color = 'royalblue')
plt.axvline(x = 5, color = 'r')

# plt.plot(t3, xy_dist3, label='90s')
plt.legend()
plt.show()
plt.savefig('D:/Anumita/PincherPlots/C3_DistancevTime.jpg')

# %% All curves

expt1 = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
folder1 = '21-12-20_M1_P1_C1_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder1+'_PY.csv'
data1 = pd.read_csv(file, sep=';')
#t1 = np.linspace(0,len(data1['T']),len(data1['T']))
t1 = (data1['T']*1000)
xy_dist1 = data1['D3'] - bead_dia

expt2 = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
folder2 = '21-12-20_M1_P1_C3_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder2+'_PY.csv'
data2 =  pd.read_csv(file, sep=';')
# t2 = np.linspace(0,len(data2['T']),len(data2['T']))
t2 = (data2['T']*1000)
xy_dist2 = data2['D3'] - bead_dia

expt3 = '20220203_100xoil_3t3optorhoa_4.5beads_15mT'
folder3 = '22-02-03_M4_P1_C1_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder3+'_PY.csv'
data3 = pd.read_csv(file, sep=';')
#t1 = np.linspace(0,len(data1['T']),len(data1['T']))
t3 = (data3['T']*1000)
xy_dist3 = data3['D3'] - bead_dia

expt4 = '20220203_100xoil_3t3optorhoa_4.5beads_15mT'
folder4 = '22-02-03_M3_P1_C1_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder4+'_PY.csv'
data4 =  pd.read_csv(file, sep=';')
# t2 = np.linspace(0,len(data2['T']),len(data2['T']))
t4 = (data4['T']*1000)
xy_dist4 = data4['D3'] - bead_dia

expt5 = '20220203_100xoil_3t3optorhoa_4.5beads_15mT'
folder5 = '22-02-03_M5_P1_C3_disc20um'
file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder5+'_PY.csv'
data5 =  pd.read_csv(file, sep=';')
# t2 = np.linspace(0,len(data2['T']),len(data2['T']))
t5 = (data5['T']*1000)
xy_dist5 = data5['D3'] - bead_dia


Nan_thresh1 = 3
Nan_thresh2 = 3

outlier = np.where(xy_dist2 > Nan_thresh2)[0]
xy_dist2[outlier] = np.nan

outlier = np.where(xy_dist1 > Nan_thresh1)[0]
xy_dist1[outlier] = np.nan

# expt3 = '20211220_100xoil_3t3optorhoa_4.5beads_15mT'
# folder3 = '21-12-20_M3_P1_C2_disc20um'
# file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+folder3+'_PY.csv'
# data3 = pd.read_csv(file, sep=';')
# t3 = np.linspace(0,len(data3['T']),len(data3['T']))
# xy_dist3 = data3['D2'] - bead_dia

plt.style.use('dark_background')


fig= plt.figure(figsize=(20,20))
# fig.suptitle(fontsize=16)

# right_side = fig.spines["right"]
# right_side.set_visible(False)
# top_side = fig.spines["top"]
# top_side.set_visible(False)


plt.rcParams.update({'font.size':35})
plt.title('3D Distance vs Time')
plt.ylabel('Thickness (nm)')
plt.xlabel('Time (secs)')

plt.plot(t1, xy_dist1, label=folder1, color = 'red')
plt.plot(t2, xy_dist2, label=folder2, color = 'blue')
plt.plot(t3, xy_dist3, label=folder3, color = 'orange')
plt.plot(t4, xy_dist4, label=folder4, color = 'pink')
plt.plot(t5, xy_dist5, label=folder5, color='yellow')
# plt.plot(t3, xy_dist3, label='90s')
plt.legend()
plt.show()
plt.savefig('D:/Anumita/PincherPlots/All_DistancevTime.jpg')


# %% Plotting with fluorescence recruitment values

expt = '20220301_100xoil_3t3optorhoa_4.5beads_15mT'
folder = '22-03-01_M1_P1_C3_disc20um'
date = '22.03.01'

file = 'D:/Anumita/MagneticPincherData/Raw/'+date+'/'+folder+'_Values.csv'
data = pd.read_csv(file, sep=',')

radius = np.asarray(data['Radius_[pixels]'])

# for i in range(np.shape(data)[1]):
#     plt

# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

bead_dia = 4.503

# %% Plotting combined graphs

path = 'D:/Anumita/CombinePlotsFolder/'
actType = 'Global'
bead_dia = 4.503

if actType == 'Away':
    path = path+'AwayBeads/'
    files = os.listdir(path)  
    for file in files:
        data =  pd.read_csv(path+file, sep=';')
        xyz_dist = data['D3'] - bead_dia
        t = (data['T']*1000)/60
        plt.plot(t, xyz_dist)
    plt.show()
elif actType == 'At':
    path = path+'AtBeads/'
    files = os.listdir(path) 
    for file in files:
        data =  pd.read_csv(path+file, sep=';')
        xyz_dist = data['D3'] - bead_dia
        t = (data['T']*1000)/60
        plt.plot(t, xyz_dist)
    plt.show()
elif actType == 'Global':
    path = path+'Global/'
    files = os.listdir(path) 
    for file in files:
        data =  pd.read_csv(path+file, sep=';')
        xyz_dist = data['D3'] - bead_dia
        t = (data['T']*1000)/60
        plt.plot(t, xyz_dist, label = file)
    plt.show()
    plt.legend()
               