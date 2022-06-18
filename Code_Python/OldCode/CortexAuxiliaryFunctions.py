# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 10:04:47 2021

@author: joseph
"""

#%%

import numpy as np
import pandas as pd
import os
from cycler import cycler
from copy import copy
import itertools
from statannot import add_stat_annotation
from datetime import date
import scipy.io as sio
# from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
pd.set_option('mode.chained_assignment',None)


# In[ ]:

import re
dateFormatExcel = re.compile('\d{2}/\d{2}/\d{4}')
dateFormatOk = re.compile('\d{2}-\d{2}-\d{2}')

# In[ ]:

from bokeh.io import output_notebook, show
output_notebook()
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Range1d
from bokeh.transform import factor_cmap
from bokeh.palettes import Category10
from bokeh.layouts import gridplot

# In[ ]:

import matplotlib
get_ipython().run_line_magic('matplotlib', 'widget')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'figure.autolayout': True})

# %matplotlib inline
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# prop_cycle = plt.rcParams['axes.prop_cycle']
# colors = prop_cycle.by_key()['color']
new_color_cycle = cycler(color=['#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
plt.rcParams['axes.prop_cycle'] = new_color_cycle


# In[ ]:

# mainDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis"
mainDir = "C://Users//josep//Desktop//ActinCortexAnalysis"
experimentalDataDir = os.path.join(mainDir, "ExperimentalData")
dataDir = os.path.join(mainDir, "DataAnalysis")
figDir = os.path.join(dataDir, "Figures")
todayFigDir = os.path.join(figDir, "Historique//" + str(date.today()))
# os.mkdir(todayFigDir)
timeSeriesDataDir = os.path.join(dataDir, "TimeSeriesData")
allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
allTimeSeriesDataFiles

# In[ ]:

def get_R2(Y1, Y2):
    meanY = np.mean(Y1)
    meanYarray = meanY*np.ones(len(Y1))
    SST = np.sum((Y1-meanYarray)**2)
    SSE = np.sum((Y2-meanYarray)**2)
    R2 = SSE/SST
    return(R2)

def getDictAggMean(df):
    dictAggMean = {}
    for c in df.columns:
#         t = df[c].dtype
#         print(c, t)
        try :
            if np.array_equal(df[c], df[c].astype(bool)):
                dictAggMean[c] = 'first'
            else:
                try:
                    np.mean(df[c])
                    dictAggMean[c] = 'mean'
                except:
                    dictAggMean[c] = 'first'
        except:
                dictAggMean[c] = 'first'
    return(dictAggMean)


#%%

    #----------------------#
    # TimeSeries functions #
    #----------------------#

# In[ ]:

def getCellTimeSeriesData(cellID, timeSeriesDataDir=timeSeriesDataDir):
    allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
    fileFound = False
    nFile = len(allTimeSeriesDataFiles)
    iFile = 0
    while (not fileFound) and (iFile < nFile):
        f = allTimeSeriesDataFiles[iFile]
        if f.startswith(cellID):
            timeSeriesDataFilePath = os.path.join(timeSeriesDataDir, f)
            timeSeriesDataFrame = pd.read_csv(timeSeriesDataFilePath, sep=';')
            fileFound = True
        iFile += 1
    if not fileFound:
        timeSeriesDataFrame = pd.DataFrame([])
    return(timeSeriesDataFrame)

def plotCellTimeSeriesData(cellID):
    X = 'T'
    Y = np.array(['B', 'F', 'dx', 'dy', 'dz', 'D2', 'D3'])
    units = np.array([' (mT)', ' (pN)', ' (µm)', ' (µm)', ' (µm)', ' (µm)', ' (µm)'])
    timeSeriesDataFrame = getCellTimeSeriesData(cellID)
    if not timeSeriesDataFrame.size == 0:
#         plt.tight_layout()
#         fig.show() # figsize=(20,20)
        axes = timeSeriesDataFrame.plot(x=X, y=Y, kind='line', ax=None, subplots=True, sharex=True, sharey=False, layout=None,                        figsize=(8,10), use_index=True, title = cellID + '- Time dependant data', grid=None, legend=False, style=None, logx=False, logy=False,                        loglog=False, xticks=None, yticks=None, xlim=None, ylim=None, rot=None, fontsize=None, colormap=None,                        table=False, yerr=None, xerr=None, secondary_y=False, sort_columns=False)
        plt.gcf().tight_layout()
        for i in range(len(Y)):
            axes[i].set_ylabel(Y[i] + units[i])
        
    else:
        print('cell not found')
        
def addExcludedCell(cellID, motive):
    f = open(os.path.join(experimentalDataDir, 'ExcludedCells.txt'), 'r')
    lines = f.readlines()
    nLines = len(lines)
    excludedCellsList = []
    for iLine in range(nLines):
        line = lines[iLine]
        splitLine = line[:-1].split(',')
        excludedCellsList.append(splitLine[0])
    if cellID in excludedCellsList:
        newlines = copy(lines)
        iLineOfInterest = excludedCellsList.index(cellID)
        if motive not in newlines[iLineOfInterest][:-1].split(','):
            newlines[iLineOfInterest] = newlines[iLineOfInterest][:-1] + ',' + motive + '\n'            
    else:
        newlines = copy(lines)
        newlines.append('' + cellID + ',' + motive + '\n')
    f.close()
    f = open(os.path.join(experimentalDataDir, 'ExcludedCells.txt'), 'w')
    f.writelines(newlines)
    
def getExcludedCells():
    f = open(os.path.join(experimentalDataDir, 'ExcludedCells.txt'), 'r')
    lines = f.readlines()
    nLines = len(lines)
    excludedCellsDict = {}
    for iLine in range(nLines):
        line = lines[iLine]
        splitLine = line[:-1].split(',')
        excludedCellsDict[splitLine[0]] = splitLine[1:]
    return(excludedCellsDict)


#%%

    #------------------------#
    # GlobalTables functions #
    #------------------------#


#%%

# (1) Experimental Conditions


# In[ ]:


def getExperimentalConditions(save = False, experimentalDataDir=experimentalDataDir):
    # Getting the table
    experimentalDataFile = 'ExperimentalConditions.csv'
    experimentalDataFilePath = os.path.join(experimentalDataDir, experimentalDataFile)
    expConditionsDF = pd.read_csv(experimentalDataFilePath, sep=';',header=0)
    print('Extracted a table with ' + str(expConditionsDF.shape[0]) + ' lines and ' + str(expConditionsDF.shape[1]) + ' columns.')
    # Cleaning the table
    try:
        for c in expConditionsDF.columns:
            if 'Unnamed' in c:
                expConditionsDF = expConditionsDF.drop([c], axis=1)
        expConditionsDF = expConditionsDF.convert_dtypes()

        listTextColumns = []
        for col in expConditionsDF.columns:
            if expConditionsDF[col].dtype == 'string':
                listTextColumns.append(col)

        expConditionsDF[listTextColumns] = expConditionsDF[listTextColumns].apply(lambda x: x.str.replace(',','.'))

        expConditionsDF['scale pixel per um'] = expConditionsDF['scale pixel per um'].astype(float)
        expConditionsDF['optical index correction'] =                   expConditionsDF['optical index correction'].apply(lambda x: x.split('/')[0]).astype(float)                 / expConditionsDF['optical index correction'].apply(lambda x: x.split('/')[1]).astype(float)
        expConditionsDF['magnetic field correction'] = expConditionsDF['magnetic field correction'].astype(float)
        expConditionsDF['with fluo images'] = expConditionsDF['with fluo images'].astype(bool)

        expConditionsDF['ramp field'] =         expConditionsDF['ramp field'].apply(lambda x: [x.split(';')[0], x.split(';')[1]] if not pd.isnull(x) else [])


        dateExemple = expConditionsDF.loc[expConditionsDF.index[1],'date']

        if re.match(dateFormatExcel, dateExemple):
            print('dates corrected')
            expConditionsDF.loc[1:,'date'] = expConditionsDF.loc[1:,'date'].apply(lambda x: x.split('/')[0] + '-' + x.split('/')[1] + '-' + x.split('/')[2][2:])        
        
    except:
        print('Unexpected bug with the cleaning step')

    if save:
        saveName = 'ExperimentalConditions.csv'
        savePath = os.path.join(experimentalDataDir, saveName)
        expConditionsDF.to_csv(savePath, sep=';')

    expConditionsDF['manipID'] = expConditionsDF['date'] + '_' + expConditionsDF['manip']
    
    return(expConditionsDF)


#%%

# (2) Constant field


# In[ ]:


def analyseTimeSeries_ctField(tS_df):
    results = {}
    results['duration'] = np.max(tS_df['T'])
    results['medianRawB'] = np.median(tS_df.B)
    results['medianThickness'] = np.median(tS_df.D3)
    results['1stDThickness'] = np.percentile(tS_df.D3, 10)
    results['9thDThickness'] = np.percentile(tS_df.D3, 90)
    results['fluctuAmpli'] = results['9thDThickness'] - results['1stDThickness']
    results['validated'] = (results['1stDThickness'] > 0)
    X, Y = tS_df['T'], tS_df['D3']
    p, residuals, rank, singular_values, rcond = np.polyfit(X, Y, deg=5, full=True)
    Y2 = np.zeros(len(X))
    for i in range(len(X)):
        deg = len(p)-1
        for k in range(deg+1):
            Y2[i] += p[k]*(X[i]**(deg-k))
    results['R2_polyFit'] = get_R2(Y, Y2)
    return(results)

def createCtFieldDataDict(list_ctFieldFiles, timeSeriesDataDir=timeSeriesDataDir):
    tableDict = {}
    tableDict['date'], tableDict['cellName'], tableDict['cellID'], tableDict['manipID'] = [], [], [], []
    tableDict['duration'], tableDict['medianRawB'], tableDict['medianThickness'] = [], [], []
    tableDict['1stDThickness'], tableDict['9thDThickness'], tableDict['fluctuAmpli'] = [], [], []
    tableDict['R2_polyFit'], tableDict['validated'] = [], []
    for f in list_ctFieldFiles:
        split_f = f.split('_')
        tableDict['date'].append(split_f[0])
        tableDict['cellName'].append(split_f[1] + '_' + split_f[2] + '_' + split_f[3])
        tableDict['cellID'].append(split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3])
        tableDict['manipID'].append(split_f[0] + '_' + split_f[1])
        tS_DataFilePath = os.path.join(timeSeriesDataDir, f)
        current_tS_df = pd.read_csv(tS_DataFilePath, ';')
        current_resultDict = analyseTimeSeries_ctField(current_tS_df)
        for k in current_resultDict.keys():
            tableDict[k].append(current_resultDict[k])
    return(tableDict)

def computeGlobalTable_ctField(task = 'fromScratch', fileName = 'Global_CtFieldData', save = False, dataDir=dataDir, timeSeriesDataDir=timeSeriesDataDir):
    ctFieldTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) \
                                  if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                                      and ('thickness' in f))]
    if task == 'fromScratch':
        # create a dict containing the data
        tableDict = createCtFieldDataDict(ctFieldTimeSeriesDataFiles)
        # create the table
        CtField_DF = pd.DataFrame(tableDict)
        
    elif task == 'updateExisting':
        # get existing table
        try:
            savePath = os.path.join(dataDir, (fileName + '.csv'))
            existing_CtField_DF = pd.read_csv(savePath, sep=';')
            for c in existing_CtField_DF.columns:
                if 'Unnamed' in c:
                    existing_CtField_DF = existing_CtField_DF.drop([c], axis=1)
        except:
            print('No existing table found')
        # find which of the time series files are new
        new_ctFieldTimeSeriesDataFiles = []
        for f in ctFieldTimeSeriesDataFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            if currentCellID not in existing_CtField_DF.cellID.values:
                new_ctFieldTimeSeriesDataFiles.append(f)
        new_tableDict = createCtFieldDataDict(new_ctFieldTimeSeriesDataFiles)
        # create the table with new data
        new_CtField_DF = pd.DataFrame(new_tableDict)
        # fuse the two
        new_CtField_DF.index += existing_CtField_DF.shape[0]
        CtField_DF = pd.concat([existing_CtField_DF, new_CtField_DF])
    
    dateExemple = CtField_DF.loc[CtField_DF.index[0],'date']
    if re.match(dateFormatExcel, dateExemple):
        CtField_DF.loc[:,'date'] = CtField_DF.loc[:,'date'].apply(lambda x: x.split('/')[0] + '-' + x.split('/')[1] + '-' + x.split('/')[2][2:])
    
    if save:
        saveName = fileName + '.csv'
        savePath = os.path.join(dataDir, saveName)
        CtField_DF.to_csv(savePath, sep=';')
        
    return(CtField_DF)

def getGlobalTable_ctField(fileName = 'Global_CtFieldData', dataDir=dataDir):
    try:
        savePath = os.path.join(dataDir, (fileName + '.csv'))
        CtField_DF = pd.read_csv(savePath, sep=';')
        for c in CtField_DF.columns:
            if 'Unnamed' in c:
                CtField_DF = CtField_DF.drop([c], axis=1)
        print('Extracted a table with ' + str(CtField_DF.shape[0]) + ' lines and ' + str(CtField_DF.shape[1]) + ' columns.')
        
    except:
        print('No existing table found')
        
    dateExemple = CtField_DF.loc[CtField_DF.index[0],'date']
    if re.match(dateFormatExcel, dateExemple):
        print('dates corrected')
        CtField_DF.loc[:,'date'] = CtField_DF.loc[:,'date'].apply(lambda x: x.split('/')[0] + '-' + x.split('/')[1] + '-' + x.split('/')[2][2:])
#         mecaDF['ManipID'] = mecaDF['ExpDay'] + '_' + mecaDF['CellName'].apply(lambda x: x.split('_')[0])
    return(CtField_DF)

#%%

# (3) Mechanics


# In[ ]:


def analyseTimeSeries_meca(tS_df):
    results = {}
    ### TBC
    return(results)

def createMecaDataDict(list_ctFieldFiles):
    tableDict = {}
    tableDict['date'], tableDict['cellName'], tableDict['cellID'] = [], [], []
#     tableDict[''], tableDict[''], tableDict[''] = [], [], []
#     tableDict[''], tableDict[''], tableDict[''] = [], [], []
#     tableDict[''] = []
# TBC
    for f in list_ctFieldFiles:
        split_f = f.split('_')
        tableDict['date'].append(split_f[0])
        tableDict['cellName'].append(split_f[1] + '_' + split_f[2] + '_' + split_f[3])
        tableDict['cellID'].append(split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3])
        tS_DataFilePath = os.path.join(timeSeriesDataDir, f)
        current_tS_df = pd.read_csv(tS_DataFilePath, ';')
        current_resultDict = analyseTimeSeries_meca(current_tS_df)
        for k in current_resultDict.keys():
            tableDict[k].append(current_resultDict[k])
    return(tableDict)

def computeGlobalTable_meca(task = 'fromScratch', fileName = 'Global_CtFieldData', save = True):
    mecaTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir)                                   if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv")                                       and ('R40' in f))] # Change to allow different formats in the future
    print(ctFieldTimeSeriesDataFiles)
    if task == 'fromScratch':
        # create a dict containing the data
        tableDict = createMecaDataDict(mecaTimeSeriesDataFiles)
        # create the table
        meca_DF = pd.DataFrame(tableDict)
        
    elif task == 'updateExisting':
        # get existing table
        try:
            savePath = os.path.join(dataDir, (fileName + '.csv'))
            existing_meca_DF = pd.read_csv(savePath, sep=';')
        except:
            print('No existing table found')
        # find which of the time series files are new
        new_mecaTimeSeriesDataFiles = []
        for f in mecaTimeSeriesDataFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            if currentCellID not in existing_meca_DF.cellID:
                new_mecaTimeSeriesDataFiles.append(f)
        new_tableDict = createCtFieldDataDict(new_mecaTimeSeriesDataFiles)
        # create the table with new data
        new_meca_DF = pd.dataframe(tableDict)
        # fuse the two
        meca_DF = pd.concat([existing_meca_DF, new_meca_DF])
        
    if save:
        saveName = fileName + '.csv'
        savePath = os.path.join(dataDir, saveName)
        meca_DF.to_csv(savePath, sep=';')
        
    return(meca_DF)
              
def getGlobalTable_meca(fileName = 'Global_mecaData'):
    try:
        savePath = os.path.join(dataDir, (fileName + '.csv'))
        meca_DF = pd.read_csv(savePath, sep=';')
        print('Extracted a table with ' + str(meca_DF.shape[0]) + ' lines and ' + str(meca_DF.shape[1]) + ' columns.')
    except:
        print('No existing table found')
    for c in meca_DF.columns:
            if 'Unnamed' in c:
                meca_DF = meca_DF.drop([c], axis=1)
    if not ('manipID' in meca_DF.columns):
        meca_DF['manipID'] = meca_DF['ExpDay'] + '_' + meca_DF['CellID'].apply(lambda x: x.split('_')[0])
    dateExemple = meca_DF.loc[meca_DF.index[0],'ExpDay']
    if re.match(dateFormatExcel, dateExemple):
        print('bad date')
    return(meca_DF)


#%%

# (4) Fluorescence


# In[ ]:


def getFluoData(save = False):
    # Getting the table
    fluoDataFile = 'FluoQuantification.csv'
    fluoDataFilePath = os.path.join(dataDir, fluoDataFile)
    fluoDF = pd.read_csv(fluoDataFilePath, sep=';',header=0)
    print('Extracted a table with ' + str(fluoDF.shape[0]) + ' lines and ' + str(fluoDF.shape[1]) + ' columns.')
    # Cleaning the table
    try:
        for c in fluoDF.columns:
            if 'Unnamed' in c:
                fluoDF = fluoDF.drop([c], axis=1)
        
    except:
        print('Unexpected bug with the cleaning step')

    if save:
        saveName = 'FluoQuantification.csv'
        savePath = os.path.join(dataDir, saveName)
        fluoDF.to_csv(savePath, sep=';')
    return(fluoDF)
