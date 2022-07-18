# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 13:07:45 2022

@author: JosephVermeil
"""

# %% (0) Imports and settings

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as st
import statsmodels.api as sm

import matplotlib
import matplotlib.pyplot as plt

import re
import os
import sys
import time
import random
import itertools

from copy import copy
from cycler import cycler
from datetime import date
from scipy.optimize import curve_fit


# Local imports
COMPUTERNAME = os.environ['COMPUTERNAME']
if COMPUTERNAME == 'ORDI-JOSEPH':
    mainDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis"
    rawDir = "D://MagneticPincherData"
    ownCloudDir = "C://Users//JosephVermeil//ownCloud//ActinCortexAnalysis"
    experimentalDataDir = os.path.join(mainDir, "Data_Experimental_JV")
elif COMPUTERNAME == 'LARISA':
    mainDir = "C://Users//Joseph//Desktop//ActinCortexAnalysis"
    rawDir = "F://JosephVermeil//MagneticPincherData"    
    ownCloudDir = "C://Users//Joseph//ownCloud//ActinCortexAnalysis"
    experimentalDataDir = os.path.join(mainDir, "Data_Experimental_JV")
elif COMPUTERNAME == '':
    mainDir = "C://Users//josep//Desktop//ActinCortexAnalysis"
    ownCloudDir = "C://Users//josep//ownCloud//ActinCortexAnalysis"
    experimentalDataDir = os.path.join(mainDir, "Data_Experimental_JV")

# Add the folder to path

sys.path.append(mainDir + "//Code_Python")
import utilityFunctions_JV as jvu

# %%% Smaller settings

pd.set_option('mode.chained_assignment',None)
pd.set_option('display.max_columns', None)

dateFormatExcel = re.compile(r'\d{2}/\d{2}/\d{4}')
dateFormatOk = re.compile(r'\d{2}-\d{2}-\d{2}')

SMALLER_SIZE = 10
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


####

dictSubstrates = {}
for i in range(0,105,5):
    dictSubstrates['disc' + str(i) + 'um'] = str(i) + 'um fibronectin discs'

# %%% Color and marker lists

colorList10 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
markerList10 = ['o', 's', 'D', '>', '^', 'P', 'X', '<', 'v', 'p']

bigPalette1 = sns.color_palette("tab20b")
bigPalette1_hex = bigPalette1.as_hex()

bigPalette2 = sns.color_palette("tab20c")
bigPalette2_hex = bigPalette2.as_hex()

colorList30 = []
for ii in range(2, -1, -1):
    colorList30.append(bigPalette2_hex[4*0 + ii]) # blue
    colorList30.append(bigPalette2_hex[4*1 + ii]) # orange
    colorList30.append(bigPalette2_hex[4*2 + ii]) # green
    colorList30.append(bigPalette1_hex[4*3 + ii]) # red
    colorList30.append(bigPalette2_hex[4*3 + ii]) # purple
    colorList30.append(bigPalette1_hex[4*2 + ii]) # yellow-brown
    colorList30.append(bigPalette1_hex[4*4 + ii]) # pink
    colorList30.append(bigPalette1_hex[4*0 + ii]) # navy    
    colorList30.append(bigPalette1_hex[4*1 + ii]) # yellow-green
    colorList30.append(bigPalette2_hex[4*4 + ii]) # gray
    

# %% (1) Directories adress

# experimentalDataDir = os.path.join(mainDir, "Data_Experimental")
dataDir = os.path.join(mainDir, "Data_Analysis")
timeSeriesDataDir = os.path.join(dataDir, "TimeSeriesData")

figDir = os.path.join(dataDir, "Figures")
todayFigDir = os.path.join(figDir, "Historique//" + str(date.today()))

figDirLocal = os.path.join(rawDir, "Figures")
todayFigDirLocal = os.path.join(figDirLocal, "Historique//" + str(date.today()))

ownCloudFigDir = os.path.join(ownCloudDir, "Data_Analysis", "Figures")
ownCloudTodayFigDir = os.path.join(ownCloudFigDir, "Historique//" + str(date.today()))

# %% (2) Utility functions

def get_R2(Y1, Y2):
    meanY = np.mean(Y1)
    meanYarray = meanY*np.ones(len(Y1))
    SST = np.sum((Y1-meanYarray)**2)
    SSE = np.sum((Y2-meanYarray)**2)
    R2 = SSE/SST
    return(R2)

def get_Chi2(Ymeas, Ymodel, dof, S):
    #### To be validated soon !
    residuals = Ymeas-Ymodel
    Chi2 = np.sum((residuals/S)**2)
    Chi2_dof = Chi2/dof
    return(Chi2_dof)

def getDictAggMean(df):
    dictAggMean = {}
    for c in df.columns:
    #         t = df[c].dtype
    #         print(c, t)
            try :
                if np.array_equal(df[c], df[c].astype(bool)):
                    dictAggMean[c] = 'min'
                else:
                    try:
                        if not c.isnull().all():
                            np.mean(df[c])
                            dictAggMean[c] = 'mean'
                    except:
                        dictAggMean[c] = 'first'
            except:
                    dictAggMean[c] = 'first'
    return(dictAggMean)

def findFirst(x, A):
    idx = (A==x).view(bool).argmax()
    return(idx)

def fitLine(X, Y):
    X = sm.add_constant(X)
    model = sm.OLS(Y, X)
    results = model.fit()
    params = results.params # Y=a*X+b ; params[0] = b,  params[1] = a
#     print(dir(results))
#     R2 = results.rsquared
#     ci = results.conf_int(alpha=0.05)
#     CovM = results.cov_params()
#     p = results.pvalues

# This is how are computed conf_int:
#
#     bse = results.bse
#     dist = stats.t
#     alpha = 0.05
#     q = dist.ppf(1 - alpha / 2, results.df_resid)
#     params = results.params
#     lower = params - q * bse
#     upper = params + q * bse
#     print(lower, upper)
    
    return(results.params, results)

def archiveFig(fig, ax, name='auto', figDir = todayFigDir, figSubDir='', dpi = 100):
    
    if not os.path.exists(figDir):
        os.makedirs(figDir)
    
    saveDir = os.path.join(figDir, figSubDir)
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    
    if name != 'auto':
        fig.savefig(os.path.join(saveDir, name + '.png'), dpi=dpi)
    
    else:
        suptitle = fig._suptitle.get_text()
        if len(suptitle) > 0:
            name = suptitle
            fig.savefig(os.path.join(saveDir, name + '.png'), dpi=dpi)
        
        else:
            try:
                N = len(ax)
                ax = ax[0]
            except:
                N = 1
                ax = ax
                
            xlabel = ax.get_xlabel()
            ylabel = ax.get_ylabel()
            if len(xlabel) > 0 and len(ylabel) > 0:
                name = ylabel + ' Vs ' + xlabel
                if N > 1:
                    name = name + '___etc'
                fig.savefig(os.path.join(saveDir, name + '.png'), dpi=dpi)
            
            else:
                title = ax.get_title()
                if len(title) > 0:
                    if N > 1:
                        name = name + '___etc'
                    fig.savefig(os.path.join(saveDir, name + '.png'), dpi=dpi)
                
                else:
                    figNum = plt.gcf().number
                    name = 'figure ' + str(figNum) 
                    fig.savefig(os.path.join(saveDir, name + '.png'), dpi=dpi)
                    
# def findInfosInFileName(f, infoType):
#     """
#     Return a given type of info from a file name.
#     Inputs : f (str), the file name.
#               infoType (str), the type of info wanted.
#               infoType can be equal to : 
#               * 'M', 'P', 'C' -> will return the number of manip (M), well (P), or cell (C) in a cellID.
#               ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'C', the function will return 8.
#               * 'manipID'     -> will return the full manip ID.
#               ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'manipID', the function will return '21-01-18_M2'.
#               * 'cellID'     -> will return the full cell ID.
#               ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'cellID', the function will return '21-01-18_M2_P1_C8'.
#     """
#     if infoType in ['M', 'P', 'C']:
#         acceptedChar = [str(i) for i in range(10)] + ['.', '-']
#         string = '_' + infoType
#         iStart = re.search(string, f).end()
#         i = iStart
#         infoString = '' + f[i]
#         while f[i+1] in acceptedChar and i < len(f)-1:
#             i += 1
#             infoString += f[i]
            
#     elif infoType == 'date':
#         datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
#         date = f[datePos.start():datePos.end()]
#         infoString = date
    
#     elif infoType == 'manipID':
#         datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
#         date = f[datePos.start():datePos.end()]
#         manip = 'M' + findInfosInFileName(f, 'M')
#         infoString = date + '_' + manip
        
#     elif infoType == 'cellID':
#         datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
#         date = f[datePos.start():datePos.end()]
#         infoString = date + '_' + 'M' + findInfosInFileName(f, 'M') + \
#                             '_' + 'P' + findInfosInFileName(f, 'P') + \
#                             '_' + 'C' + findInfosInFileName(f, 'C')

    
#     return(infoString)
                    
                    
# %% (3) TimeSeries functions

def getCellTimeSeriesData(cellID, fromPython = True):
    if fromPython:
        allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) 
                              if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) 
                                  and f.endswith("PY.csv"))]
    else:
        allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) 
                              if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) 
                                  and f.endswith(".csv") and not f.endswith("PY.csv"))]
    fileFound = False
    nFile = len(allTimeSeriesDataFiles)
    iFile = 0
    while (not fileFound) and (iFile < nFile):
        f = allTimeSeriesDataFiles[iFile]
        if f.startswith(cellID + '_'):
            timeSeriesDataFilePath = os.path.join(timeSeriesDataDir, f)
            timeSeriesDataFrame = pd.read_csv(timeSeriesDataFilePath, sep=';')
            fileFound = True
        iFile += 1
    if not fileFound:
        timeSeriesDataFrame = pd.DataFrame([])
    else:
        for c in timeSeriesDataFrame.columns:
                if 'Unnamed' in c:
                    timeSeriesDataFrame = timeSeriesDataFrame.drop([c], axis=1)
    return(timeSeriesDataFrame)

def plotCellTimeSeriesData(cellID, fromPython = True):
    X = 'T'
    Y = np.array(['B', 'F', 'dx', 'dy', 'dz', 'D2', 'D3'])
    units = np.array([' (mT)', ' (pN)', ' (µm)', ' (µm)', ' (µm)', ' (µm)', ' (µm)'])
    timeSeriesDataFrame = getCellTimeSeriesData(cellID, fromPython)
    print(timeSeriesDataFrame.shape)
    # my_default_color_cycle = cycler(color=my_default_color_list)
    # plt.rcParams['axes.prop_cycle'] = my_default_color_cycle
    if not timeSeriesDataFrame.size == 0:
#         plt.tight_layout()
#         fig.show() # figsize=(20,20)
        axes = timeSeriesDataFrame.plot(x=X, y=Y, kind='line', ax=None, subplots=True, sharex=True, sharey=False, layout=None, \
                       figsize=(8,10), use_index=True, title = cellID + ' - f(t)', grid=None, legend=False, style=None, logx=False, logy=False, \
                       loglog=False, xticks=None, yticks=None, xlim=None, ylim=None, rot=None, fontsize=None, colormap=None, \
                       table=False, yerr=None, xerr=None, secondary_y=False, sort_columns=False)
        plt.gcf().tight_layout()
        for i in range(len(Y)):
            axes[i].set_ylabel(Y[i] + units[i])
        # plt.gcf().show()
        plt.show()
    else:
        print('cell not found')
    # plt.rcParams['axes.prop_cycle'] = my_default_color_cycle
        
def getCellTrajData(cellID, Ntraj = 2):
    trajDir = os.path.join(timeSeriesDataDir, 'Trajectories')
    allTrajFiles = [f for f in os.listdir(trajDir) 
                    if (os.path.isfile(os.path.join(trajDir, f)) 
                        and f.endswith(".csv"))]
    fileFound = 0
    nFile = len(allTrajFiles)
    iFile = 0
    listTraj = []
    while (fileFound < Ntraj) and (iFile < nFile):
        f = allTrajFiles[iFile]
        if f.startswith(cellID):
            fileFound += 1
            trajFilePath = os.path.join(trajDir, f)
            trajDataFrame = pd.read_csv(trajFilePath, sep='\t')
            for c in trajDataFrame.columns:
                if 'Unnamed' in c:
                    trajDataFrame = trajDataFrame.drop([c], axis=1)
            
            pos = 'na'
            if 'In' in f:
                pos = 'in'
            elif 'Out' in f:
                pos = 'out'
            
            dictTraj = {}
            dictTraj['df'] = trajDataFrame
            dictTraj['pos'] = pos
            dictTraj['path'] = trajFilePath
            
            listTraj.append(dictTraj)
            
        iFile += 1

    return(listTraj)
        
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

# %% (4) GlobalTables functions

# %%% (4.1) Exp conditions

# %%% (4.2) Constant Field experiments

listColumnsCtField = ['date','cellName','cellID','manipID',\
                      'duration','medianRawB','medianThickness',\
                      '1stDThickness','9thDThickness','fluctuAmpli',\
                      'R2_polyFit','validated']

def analyseTimeSeries_ctField(f, tsDf, expDf):
    results = {}
    
    thisManipID = jvu.findInfosInFileName(f, 'manipID')
    thisExpDf = expDf.loc[expDf['manipID'] == thisManipID]
    # Deal with the asymmetric pair case : the diameter can be for instance 4503 (float) or '4503_2691' (string)
    diameters = thisExpDf.at[thisExpDf.index.values[0], 'bead diameter'].split('_')
    if len(diameters) == 2:
        DIAMETER = (int(diameters[0]) + int(diameters[1]))/2.
    else:
        DIAMETER = int(diameters[0])
    
    results['duration'] = np.max(tsDf['T'])
    results['medianRawB'] = np.median(tsDf.B)
    results['medianThickness'] = (1000*np.median(tsDf.D3))-DIAMETER
    results['1stDThickness'] = (1000*np.percentile(tsDf.D3, 10))-DIAMETER
    results['9thDThickness'] = (1000*np.percentile(tsDf.D3, 90))-DIAMETER
    results['fluctuAmpli'] = (results['9thDThickness'] - results['1stDThickness'])
    results['validated'] = (results['1stDThickness'] > 0)

    
    # R2 polyfit to see the regularity. Order = 5 !
    X, Y = tsDf['T'], tsDf['D3']
    p, residuals, rank, singular_values, rcond = np.polyfit(X, Y, deg=5, full=True)
    Y2 = np.zeros(len(X))
    for i in range(len(X)):
        deg = len(p)-1
        for k in range(deg+1):
            Y2[i] += p[k]*(X[i]**(deg-k))
    results['R2_polyFit'] = jvu.get_R2(Y, Y2)
    return(results)



def createDataDict_ctField(list_ctFieldFiles):
    tableDict = {}
    tableDict['date'], tableDict['cellName'], tableDict['cellID'], tableDict['manipID'] = [], [], [], []
    tableDict['duration'], tableDict['medianRawB'], tableDict['medianThickness'] = [], [], []
    tableDict['1stDThickness'], tableDict['9thDThickness'], tableDict['fluctuAmpli'] = [], [], []
    tableDict['R2_polyFit'], tableDict['validated'] = [], []
    expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
    for f in list_ctFieldFiles:
        split_f = f.split('_')
        tableDict['date'].append(split_f[0])
        tableDict['cellName'].append(split_f[1] + '_' + split_f[2] + '_' + split_f[3])
        tableDict['cellID'].append(split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3])
        tableDict['manipID'].append(split_f[0] + '_' + split_f[1])
        tS_DataFilePath = os.path.join(timeSeriesDataDir, f)
        current_tsDf = pd.read_csv(tS_DataFilePath, ';')
        current_resultDict = analyseTimeSeries_ctField(f, current_tsDf, expDf)
        for k in current_resultDict.keys():
            tableDict[k].append(current_resultDict[k])
    return(tableDict)

# def computeGlobalTable_meca(task = 'fromScratch', fileName = 'Global_MecaData', save = False, PLOT = False, \
#                             source = 'Matlab', listColumnsMeca=listColumnsMeca):
#     """
#     Compute the GlobalTable_meca from the time series data files.
#     Option task='fromScratch' will analyse all the time series data files and construct a new GlobalTable from them regardless of the existing GlobalTable.
#     Option task='updateExisting' will open the existing GlobalTable and determine which of the time series data files are new ones, and will append the existing GlobalTable with the data analysed from those new fils.
#     listColumnsMeca have to contain all the fields of the table that will be constructed.
#     """
#     top = time.time()
    
# #     list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
# #                       if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
# #                       and ('R40' in f))] # Change to allow different formats in the future
    
#     suffixPython = '_PY'
#     if source == 'Matlab':
#         list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
#                       if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
#                       and ('R40' in f) and not (suffixPython in f))]
        
#     elif source == 'Python':
#         list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
#                       if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
#                       and (('R40' in f) or ('L40' in f)) and (suffixPython in f))]
#         # print(list_mecaFiles)

def computeGlobalTable_ctField(task = 'fromScratch', fileName = 'Global_CtFieldData', save = False,
                               source = 'Matlab'):
    """
    Compute the GlobalTable_ctField from the time series data files.
    > Option task='fromScratch' will analyse all the time series data files and construct 
    a new GlobalTable from them regardless of the existing GlobalTable.
    > Option task='updateExisting' will open the existing GlobalTable and determine 
    which of the time series data files are new ones, and will append the existing GlobalTable 
    with the data analysed from those new files.
    > Else, having task= a date, a cellID or a manipID will create a globalTable with this source only.
    """
    ctFieldFiles = [f for f in os.listdir(timeSeriesDataDir) \
                                  if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                                      and ('thickness' in f))]
        
#     print(ctFieldFiles)

    suffixPython = '_PY'
    if source == 'Matlab':
        ctFieldFiles = [f for f in os.listdir(timeSeriesDataDir) \
                      if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                      and ('thickness' in f) and not (suffixPython in f))]
        
    elif source == 'Python':
        ctFieldFiles = [f for f in os.listdir(timeSeriesDataDir) \
                      if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                      and ('thickness' in f) and (suffixPython in f))]
        # print(list_mecaFiles)


    if task == 'fromScratch':
        # create a dict containing the data
        tableDict = createDataDict_ctField(ctFieldFiles) # MAIN SUBFUNCTION
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
        new_ctFieldFiles = []
        for f in ctFieldFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            if currentCellID not in existing_CtField_DF.cellID.values:
                new_ctFieldFiles.append(f)
        new_tableDict = createDataDict_ctField(new_ctFieldFiles) # MAIN SUBFUNCTION
        # create the table with new data
        new_CtField_DF = pd.DataFrame(new_tableDict)
        # fuse the two
        new_CtField_DF.index += existing_CtField_DF.shape[0]
        CtField_DF = pd.concat([existing_CtField_DF, new_CtField_DF])
        
    else: # If task is neither 'fromScratch' nor 'updateExisting'
    # Then task can be a substring that can be in some timeSeries file !
    # It will create a table with only these files, WITH OR WITHOUT SAVING IT !
    # And it can plot figs from it.
        # save = False
        new_ctFieldFiles = []
        for f in new_ctFieldFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            if task in currentCellID:
                new_ctFieldFiles.append(f)
        # create the dict with new data
        new_tableDict = createDataDict_ctField(new_ctFieldFiles) # MAIN SUBFUNCTION
        # create the dataframe from it
        CtField_DF = pd.DataFrame(new_tableDict)
    
    dateExemple = CtField_DF.loc[CtField_DF.index[0],'date']
    if re.match(dateFormatExcel, dateExemple):
        CtField_DF.loc[:,'date'] = CtField_DF.loc[:,'date'].apply(lambda x: x.split('/')[0] + '-' + x.split('/')[1] + '-' + x.split('/')[2][2:])
    
    for c in CtField_DF.columns:
        if 'Unnamed' in c:
            CtField_DF = CtField_DF.drop([c], axis=1)
    
    if save:
        saveName = fileName + '.csv'
        savePath = os.path.join(dataDir, saveName)
        CtField_DF.to_csv(savePath, sep=';')
        
    return(CtField_DF)



def getGlobalTable_ctField(fileName = 'Global_CtFieldData'):
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

# %%% (4.3) Compressions experiments

#### Workflow
# * analyseTimeSeries_meca() analyse 1 file and return the dict (with the results of the analysis)
# * createMecaDataDict() call the previous function on the given list of files and concatenate the results
# * computeGlobalTable_meca() call the previous function and convert the dict to a DataFrame


listColumnsMeca = ['date','cellName','cellID','manipID',
                   'compNum','compDuration','compStartTime','compAbsStartTime','compStartTimeThisDay',
                   'initialThickness','minThickness','maxIndent','previousThickness',
                   'surroundingThickness','surroundingDx','surroundingDz',
                   'validatedThickness', 'jumpD3',
                   'minForce', 'maxForce', 'minStress', 'maxStress', 'minStrain', 'maxStrain',
                   'ctFieldThickness','ctFieldFluctuAmpli','ctFieldDX','ctFieldDZ',
                   'bestH0',
                   'H0Chadwick','EChadwick','R2Chadwick','EChadwick_CIWidth',
                   'hysteresis',
                   'critFit', 'validatedFit','comments'] # 'fitParams',

#### SETTING ! Fit Selection R2 & Chi2
dictSelectionCurve = {'R2' : 0.6, 'Chi2' : 10, 'Error' : 0.02}

#### SETTING ! Change the region fits NAMES here 

#### >>> OPTION 1 - MANY FITS, FEW PLOTS

# fit_intervals = [S for S in range(0,450,50)] + [S for S in range(500,1100,100)] + [S for S in range(1500,2500,500)]
# regionFitsNames = []
# for ii in range(len(fit_intervals)-1):
#     for jj in range(ii+1, len(fit_intervals)):
#         regionFitsNames.append(str(fit_intervals[ii]) + '<s<' + str(fit_intervals[jj]))
        
# fit_toPlot = ['50<s<200', '200<s<350', '350<s<500', '500<s<700', '700<s<1000', '1000<s<1500', '1500<s<2000']


#### >>> OPTION 2 - TO LIGHTEN THE COMPUTATION
# regionFitsNames = fit_toPlot


#### >>> OPTION 3 - OLIVIA'S IDEA
# fitMin = [S for S in range(25,1225,50)]
# fitMax = [S+150 for S in fitMin]
# fitCenters = np.array([S+75 for S in fitMin])
# regionFitsNames = [str(fitMin[ii]) + '<s<' + str(fitMax[ii]) for ii in range(len(fitMin))]
# fit_toPlot = [regionFitsNames[ii] for ii in range(0, len(regionFitsNames), 2)]

# mask_fitToPlot = np.array(list(map(lambda x : x in fit_toPlot, regionFitsNames)))

# for rFN in regionFitsNames:
#     listColumnsMeca += ['KChadwick_'+rFN, 'K_CIW_'+rFN, 'R2Chadwick_'+rFN, 'K2Chadwick_'+rFN, 
#                         'H0Chadwick_'+rFN, 'Npts_'+rFN, 'validatedFit_'+rFN] 



#### >>> OPTION 4 - Longe ranges
# intervalSize = 250
# fitMin = [S for S in range(25,975,50)] + [S for S in range(975,2125,100)]
# fitMax = [S+intervalSize for S in fitMin]
# fitCenters = np.array([S+0.5*intervalSize for S in fitMin])
# regionFitsNames = [str(fitMin[ii]) + '<s<' + str(fitMax[ii]) for ii in range(len(fitMin))]
# fit_toPlot = [regionFitsNames[ii] for ii in range(0, len(regionFitsNames), 4)]

# mask_fitToPlot = np.array(list(map(lambda x : x in fit_toPlot, regionFitsNames)))

# for rFN in regionFitsNames:
#     listColumnsMeca += ['KChadwick_'+rFN, 'K_CIW_'+rFN, 'R2Chadwick_'+rFN, 'K2Chadwick_'+rFN, 
#                         'H0Chadwick_'+rFN, 'Npts_'+rFN, 'validatedFit_'+rFN] 
    # 'H0Chadwick_'+rFN, 'EChadwick_'+rFN
    
#### >>> OPTION 5 - Effect of range width

fitC =  np.array([S for S in range(100, 1150, 50)])
fitW = [150, 200, 250, 300]

fitCenters = np.array([[int(S) for S in fitC] for w in fitW]).flatten()
fitWidth = np.array([[int(w) for S in fitC] for w in fitW]).flatten()

fitMin = np.array([[int(S-(w/2)) for S in fitC] for w in fitW]).flatten()
fitMax = np.array([[int(S+(w/2)) for S in fitC] for w in fitW]).flatten()

fitCenters, fitWidth = fitCenters[fitMin>0], fitWidth[fitMin>0]
fitMin, fitMax = fitMin[fitMin>0], fitMax[fitMin>0]

regionFitsNames = ['S='  + str(fitCenters[ii]) + '+/-' + str(int(fitWidth[ii]//2)) for ii in range(len(fitCenters))]

fit_toPlot = [regionFitsNames[ii] for ii in range(0, len(fitC), 2)]
mask_fitToPlot = np.array(list(map(lambda x : x in fit_toPlot, regionFitsNames)))

for rFN in regionFitsNames:
    listColumnsMeca += ['KChadwick_'+rFN, 'K_CIW_'+rFN, 'R2Chadwick_'+rFN, 'K2Chadwick_'+rFN, 
                        'H0Chadwick_'+rFN, 'Npts_'+rFN, 'validatedFit_'+rFN] 
    # 'H0Chadwick_'+rFN, 'EChadwick_'+rFN




def compressionFitChadwick(hCompr, fCompr, DIAMETER):
    
    error = False
    
    def chadwickModel(h, E, H0):
        R = DIAMETER/2
        f = (np.pi*E*R*((H0-h)**2))/(3*H0)
        return(f)

    def inversedChadwickModel(f, E, H0):
        R = DIAMETER/2
        h = H0 - ((3*H0*f)/(np.pi*E*R))**0.5
        return(h)

    # some initial parameter values - must be within bounds
    initH0 = max(hCompr) # H0 ~ h_max
    initE = (3*max(hCompr)*max(fCompr))/(np.pi*(DIAMETER/2)*(max(hCompr)-min(hCompr))**2) # E ~ 3*H0*F_max / pi*R*(H0-h_min)²
    # initH0, initE = initH0*(initH0>0), initE*(initE>0)
    
    initialParameters = [initE, initH0]

    # bounds on parameters - initial parameters must be within these
    lowerBounds = (0, 0)
    upperBounds = (np.Inf, np.Inf)
    parameterBounds = [lowerBounds, upperBounds]


    try:
        params, covM = curve_fit(inversedChadwickModel, fCompr, hCompr, initialParameters, bounds = parameterBounds)

        # Previously I fitted with y=F and x=H, but it didn't work so well cause H(t) isn't monotonous:
        # params, covM = curve_fit(chadwickModel, hCompr, fCompr, initialParameters, bounds = parameterBounds)
        # Fitting with the 'inverse Chadwick model', with y=H and x=F is more convenient

        E, H0 = params
        hPredict = inversedChadwickModel(fCompr, E, H0)
        err = dictSelectionCurve['Error']
        
        comprMat = np.array([hCompr, fCompr]).T
        comprMatSorted = comprMat[comprMat[:, 0].argsort()]
        hComprSorted, fComprSorted = comprMatSorted[:, 0], comprMatSorted[:, 1]
        fPredict = chadwickModel(hComprSorted, E, H0)
        
        # Stress and strain
        deltaCompr = (H0 - hCompr)/1000 # µm
        stressCompr = fCompr / (np.pi * (DIAMETER/2000) * deltaCompr)
        strainCompr = deltaCompr / (3*(H0/1000)) 
        strainPredict = stressCompr / (E*1e6) #((H0 - hPredict)/1000) / (3*(H0/1000))
        
        # residuals_h = hCompr-hPredict
        # residuals_f = fComprSorted-fPredict

        alpha = 0.975
        dof = len(fCompr)-len(params)
        q = st.t.ppf(alpha, dof) # Student coefficient
        R2 = jvu.get_R2(hCompr, hPredict)
        
        Chi2 = jvu.get_Chi2(strainCompr, strainPredict, dof, err)

        varE = covM[0,0]
        seE = (varE)**0.5
        E, seE = E*1e6, seE*1e6
        confIntE = [E-q*seE, E+q*seE]
        confIntEWidth = 2*q*seE

        varH0 = covM[1,1]
        seH0 = (varH0)**0.5
        confIntH0 = [H0-q*seH0, H0+q*seH0]
        confIntH0Width = 2*q*seH0
        
        
    except:
        error = True
        E, H0, hPredict, R2, Chi2, confIntE, confIntH0 = -1, -1, np.ones(len(hCompr))*(-1), -1, -1, [-1,-1], [-1,-1]
    
    return(E, H0, hPredict, R2, Chi2, confIntE, confIntH0, error)





def compressionFitChadwick_StressStrain(hCompr, fCompr, H0, DIAMETER):
    
    error = False
    
    # def chadwickModel(h, E, H0):
    #     R = DIAMETER/2
    #     f = (np.pi*E*R*((H0-h)**2))/(3*H0)
    #     return(f)

    # def inversedChadwickModel(f, E, H0):
    #     R = DIAMETER/2
    #     h = H0 - ((3*H0*f)/(np.pi*E*R))**0.5
    #     return(h)
    
    def computeStress(f, h, H0):
        R = DIAMETER/2000
        delta = (H0 - h)/1000
        stress = f / (np.pi * R * delta)
        return(stress)
        
    def computeStrain(h, H0):
        delta = (H0 - h)/1000
        strain = delta / (3 * (H0/1000))
        return(strain)
    
    def constitutiveRelation(strain, K, stress0):
        stress = (K * strain) + stress0
        return(stress)
    
    def inversedConstitutiveRelation(stress, K, strain0):
        strain = (stress / K) + strain0
        return(strain)

    # some initial parameter values - must be within bounds
    initK = (3*max(hCompr)*max(fCompr))/(np.pi*(DIAMETER/2)*(max(hCompr)-min(hCompr))**2) # E ~ 3*H0*F_max / pi*R*(H0-h_min)²
    init0 = 0
    
    initialParameters = [initK, init0]

    # bounds on parameters - initial parameters must be within these
    lowerBounds = (0, -np.Inf)
    upperBounds = (np.Inf, np.Inf)
    parameterBounds = [lowerBounds, upperBounds]


    try:
        strainCompr = computeStrain(hCompr, H0)
        stressCompr = computeStress(fCompr, hCompr, H0)
        
        params, covM = curve_fit(inversedConstitutiveRelation, stressCompr, strainCompr, initialParameters, bounds = parameterBounds)

        K, strain0 = params
        strainPredict = inversedConstitutiveRelation(stressCompr, K, strain0)
        err = dictSelectionCurve['Error']

        alpha = 0.975
        dof = len(stressCompr)-len(params)
        
        R2 = jvu.get_R2(strainCompr, strainPredict)
        
        Chi2 = jvu.get_Chi2(strainCompr, strainPredict, dof, err)

        varK = covM[0,0]
        seK = (varK)**0.5
        
        q = st.t.ppf(alpha, dof) # Student coefficient
        # K, seK = K*1e6, seK*1e6
        confIntK = [K-q*seK, K+q*seK]
        confIntKWidth = 2*q*seK
        
        
    except:
        error = True
        K, strainPredict, R2, Chi2, confIntK = -1, np.ones(len(strainCompr))*(-1), -1, -1, [-1,-1]
    
    return(K, strainPredict, R2, Chi2, confIntK, error)



def analyseTimeSeries_meca(f, tsDF, expDf, listColumnsMeca, task, PLOT, PLOT_SHOW):
    
    plotSmallElements = True
    
    #### (0) Import experimental infos
    split_f = f.split('_')
    tsDF.dx, tsDF.dy, tsDF.dz, tsDF.D2, tsDF.D3 = tsDF.dx*1000, tsDF.dy*1000, tsDF.dz*1000, tsDF.D2*1000, tsDF.D3*1000
    thisManipID = jvu.findInfosInFileName(f, 'manipID')
    thisExpDf = expDf.loc[expDf['manipID'] == thisManipID]

    # Deal with the asymmetric pair case : the diameter can be for instance 4503 (float) or '4503_2691' (string)
    diameters = thisExpDf.at[thisExpDf.index.values[0], 'bead diameter'].split('_')
    if len(diameters) == 2:
        DIAMETER = (int(diameters[0]) + int(diameters[1]))/2.
    else:
        DIAMETER = int(diameters[0])
    
    EXPTYPE = str(thisExpDf.at[thisExpDf.index.values[0], 'experimentType'])

    results = {}
    for c in listColumnsMeca:
        results[c] = []

    Ncomp = max(tsDF['idxAnalysis'])
    loopStruct = thisExpDf.at[thisExpDf.index.values[0], 'loop structure'].split('_')
    nUplet = thisExpDf.at[thisExpDf.index.values[0], 'normal field multi images']
    
    if 'compression' in EXPTYPE:
    
        loop_totalSize = int(loopStruct[0])
        if len(loopStruct) >= 2:
            loop_rampSize = int(loopStruct[1])
        else:
            loop_rampSize = 0
        if len(loopStruct) >= 3:
            loop_excludedSize = int(loopStruct[2])
        else:
            loop_excludedSize = 0
        loop_ctSize = int((loop_totalSize - (loop_rampSize+loop_excludedSize))/nUplet)

    Ncomp = max(tsDF['idxAnalysis'])
    normalField = thisExpDf.at[thisExpDf.index.values[0], 'normal field']
    normalField = int(normalField)
    
    if 'compression' in EXPTYPE:
        compField = thisExpDf.at[thisExpDf.index.values[0], 'ramp field'].split('_')
        minCompField = float(compField[0])
        maxCompField = float(compField[1])

    #### (1) Get global values
    # These values are computed once for the whole cell D3 time series, but since the table has 1 line per compression, 
    # that same value will be put in the table for each line corresponding to that cell
    ctFieldH = (tsDF.loc[tsDF['idxAnalysis'] == 0, 'D3'].values - DIAMETER)
    ctFieldDX = np.median(tsDF.loc[tsDF['idxAnalysis'] == 0, 'dx'].values - DIAMETER)
    ctFieldDZ = np.median(tsDF.loc[tsDF['idxAnalysis'] == 0, 'dz'].values)
    ctFieldThickness   = np.median(ctFieldH)
    ctFieldFluctuAmpli = np.percentile(ctFieldH, 90) - np.percentile(ctFieldH,10)
    
    #### PLOT [1/4]
    # First part of the plot [mainly ax1 and ax1bis]
    if PLOT:
        # 1st plot - fig1 & ax1, is the F(t) curves with the different compressions colored depending of the analysis success
        fig1, ax1 = plt.subplots(1,1,figsize=(tsDF.shape[0]*(1/100),5))
        color = 'blue'
        ax1.set_xlabel('t (s)')
        ax1.set_ylabel('h (nm)', color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        
        nColsSubplot = 5
        nRowsSubplot = ((Ncomp-1) // nColsSubplot) + 1
        
        # 2nd plot - fig2 & ax2, gather all the F(h) curves, and will be completed later in the code
        fig2, ax2 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
        
        # 3rd plot - fig3 & ax3, gather all the stress-strain curves, and will be completed later in the code
        fig3, ax3 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
        
        # 4th plot - fig4 & ax4, gather all the F(h) curves with local fits, and will be completed later in the code
        fig4, ax4 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
        alreadyLabeled4 = []
        
        # 5th plot - fig5 & ax5, gather all the stress-strain curves with local fits, and will be completed later in the code
        fig5, ax5 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
        alreadyLabeled5 = []
        
        # 6th plot - fig6 & ax6
        fig6, ax6 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
        
        if plotSmallElements:
            # 7th plot - fig7 & ax7, gather all the "small elements", and will be completed later in the code
            fig7, ax7 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))



    for i in range(1, Ncomp+1):#Ncomp+1):

        #### (1) Identifiers
        results['date'].append(split_f[0])
        results['cellName'].append(split_f[1] + '_' + split_f[2] + '_' + split_f[3])
        results['cellID'].append(split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3])
        results['manipID'].append(split_f[0] + '_' + split_f[1])
        
        #### Jump correction
        if EXPTYPE == 'compressionsLowStart' or normalField == minCompField:
            colToCorrect = ['dx', 'dy', 'dz', 'D2', 'D3']
            maskCompAndPrecomp = np.abs(tsDF['idxAnalysis']) == i
            iStart = jvu.findFirst(np.abs(tsDF['idxAnalysis']), i)
            for c in colToCorrect:
                jump = np.mean(tsDF[c].values[iStart:iStart+3]) - tsDF[c].values[iStart-1]
                tsDF.loc[maskCompAndPrecomp, c] -= jump
                if c == 'D3':
                    D3corrected = True
                    jumpD3 = jump
        else:
            D3corrected = False
            jumpD3 = 0
            
        #### (2) Segment the compression n°i
        thisCompDf = tsDF.loc[tsDF['idxAnalysis'] == i,:]
        iStart = (jvu.findFirst(tsDF['idxAnalysis'], i))
        iStop = iStart+thisCompDf.shape[0]

        # Easy-to-get parameters
        results['compNum'].append(i)
        results['compDuration'].append(thisExpDf.at[thisExpDf.index.values[0], 'compression duration'])
        results['compStartTime'].append(thisCompDf['T'].values[0])
        results['compAbsStartTime'].append(thisCompDf['Tabs'].values[0])
        results['compStartTimeThisDay'].append(thisCompDf['Tabs'].values[0])

        listB = thisCompDf.B.values
        
        # Test to check if most of the compression have not been deleted due to bad image quality 
        highBvalues = (listB > (maxCompField + minCompField)/2)
        N_highBvalues = np.sum(highBvalues)
        testHighVal = (N_highBvalues > 20)

        # Test to check if the range of B field is large enough
        minB, maxB = min(listB), max(listB)
        testRangeB = ((maxB-minB) > 0.7*(maxCompField - minCompField))
        thresholdB = (maxB-minB)/50
        thresholdDeltaB = (maxB-minB)/400

        # Is the curve ok to analyse ?
        doThisCompAnalysis = testHighVal and testRangeB # Some criteria can be added here

        if doThisCompAnalysis:
            #### (3) Inside the compression n°i, delimit the compression and relaxation phases

            # Delimit the start of the increase of B (typically the moment when the field decrease from 5 to 3)
            # and the end of its decrease (typically when it goes back from 3 to 5)
            
            try:
                # Correct for bugs in the B data
                if 'compressions' in EXPTYPE:
                    for k in range(1,len(listB)):
                        B = listB[k]
                        if B > 1.25*maxCompField:
                            listB[k] = listB[k-1]

                offsetStart, offsetStop = 0, 0
                minB, maxB = min(listB), max(listB)
                thresholdB = (maxB-minB)/50
                thresholdDeltaB = (maxB-minB)/400 # NEW CONDITION for the beginning of the compression : 
                # remove the first points where the steps in B are very very small

                k = 0
                while (listB[k] < minB+thresholdB) or (listB[k+1]-listB[k] < thresholdDeltaB):
                    offsetStart += int((listB[k] < minB+thresholdB) or ((listB[k+1]-listB[k]) < thresholdDeltaB))
                    k += 1

                k = 0
                while (listB[-1-k] < minB+thresholdB):
                    offsetStop += int(listB[-1-k] < minB+thresholdB)
                    k += 1

                jStart = offsetStart # Beginning of compression
                jMax = np.argmax(thisCompDf.B) # End of compression, beginning of relaxation
                jStop = thisCompDf.shape[0] - offsetStop # End of relaxation
            
            except:
                print(listB)
                print(testRangeB, thresholdB, thresholdDeltaB)

            # Four arrays
            hCompr = (thisCompDf.D3.values[jStart:jMax+1] - DIAMETER)
            hRelax = (thisCompDf.D3.values[jMax+1:jStop] - DIAMETER)
            fCompr = (thisCompDf.F.values[jStart:jMax+1])
            fRelax = (thisCompDf.F.values[jMax+1:jStop])


            # Refinement of the compression delimitation.
            # Remove the 1-2 points at the begining where there is just the viscous relaxation of the cortex
            # because of the initial decrease of B and the cortex thickness increases.
            offsetStart2 = 0
            k = 0
            while (k<len(hCompr)-10) and (hCompr[k] < np.max(hCompr[k+1:min(k+10, len(hCompr))])):
                offsetStart2 += 1
                k += 1
            # Better compressions arrays
            hCompr = hCompr[offsetStart2:]
            fCompr = fCompr[offsetStart2:]

            # Get the points of constant field preceding and surrounding the current compression
            # Ex : if the labview code was set so that there is 6 points of ct field before and after each compression,
            # previousPoints will contains D3[iStart-12:iStart]
            # surroundingPoints will contains D3[iStart-6:iStart] and D3[iStop:iStop+6]
            previousPoints = (tsDF.D3.values[max(0,iStart-(loop_ctSize)):iStart]) - DIAMETER
            surroundingPoints = np.concatenate([tsDF.D3.values[max(0,iStart-(loop_ctSize//2)):iStart],tsDF.D3.values[iStop:iStop+(loop_ctSize//2)]]) - DIAMETER
            surroundingPointsX = np.concatenate([tsDF.dx.values[max(0,iStart-(loop_ctSize//2)):iStart],tsDF.dx.values[iStop:iStop+(loop_ctSize//2)]]) - DIAMETER
            surroundingPointsZ = np.concatenate([tsDF.dz.values[max(0,iStart-(loop_ctSize//2)):iStart],tsDF.dz.values[iStop:iStop+(loop_ctSize//2)]])
            
            # Parameters relative to the thickness ( = D3-DIAMETER)
            results['initialThickness'].append(np.mean(hCompr[0:3]))
            results['minThickness'].append(np.min(hCompr))
            results['maxIndent'].append(results['initialThickness'][-1] - results['minThickness'][-1])
            results['previousThickness'].append(np.median(previousPoints))
            results['surroundingThickness'].append(np.median(surroundingPoints))
            results['surroundingDx'].append(np.median(surroundingPointsX))
            results['surroundingDz'].append(np.median(surroundingPointsZ))
            results['ctFieldDX'].append(ctFieldDX)
            results['ctFieldDZ'].append(ctFieldDZ)
            results['ctFieldThickness'].append(ctFieldThickness)
            results['ctFieldFluctuAmpli'].append(ctFieldFluctuAmpli)
            results['jumpD3'].append(jumpD3)

            validatedThickness = np.min([results['initialThickness'],results['minThickness'],results['previousThickness'],\
                                        results['surroundingThickness'],results['ctFieldThickness']]) > 0
            results['validatedThickness'].append(validatedThickness)
            
            results['minForce'].append(np.min(fCompr))
            results['maxForce'].append(np.max(fCompr))




            #### (4) Fit with Chadwick model of the force-thickness curve
            
            #### Classic, all curve, Chadwick fit
            E, H0, hPredict, R2, Chi2, confIntE, confIntH0, fitError = compressionFitChadwick(hCompr, fCompr, DIAMETER) # IMPORTANT SUBFUNCTION

            R2CRITERION = dictSelectionCurve['R2']
            CHI2CRITERION = dictSelectionCurve['Chi2']
            critFit = 'R2 > ' + str(R2CRITERION)
            results['critFit'].append(critFit)
            validatedFit = ((R2 > R2CRITERION) and (Chi2 < CHI2CRITERION))


            if not fitError:
                results['validatedFit'].append(validatedFit)
                if validatedFit:
                    results['comments'].append('ok')
                else:
                    results['comments'].append('R2 < ' + str(R2CRITERION))

                confIntEWidth = abs(confIntE[0] - confIntE[1])
                results['H0Chadwick'].append(H0)
                results['EChadwick'].append(E)
                results['R2Chadwick'].append(R2)
                results['EChadwick_CIWidth'].append(confIntEWidth)
                
            elif fitError:
                validatedFit = False
                results['H0Chadwick'].append(np.nan)
                results['EChadwick'].append(np.nan)
                results['R2Chadwick'].append(np.nan)
                results['EChadwick_CIWidth'].append(np.nan)
                results['validatedFit'].append(validatedFit)
                results['comments'].append('fitFailure')
                
                
            #### (4.0) Find the best H0
            findH0_NbPts = 15
            findH0_E, bestH0, findH0_hPredict, findH0_R2, findH0_Chi2, findH0_confIntE, findH0_confIntH0, findH0_fitError = \
                compressionFitChadwick(hCompr[:findH0_NbPts], fCompr[:findH0_NbPts], DIAMETER)
                
            maxH0 = max(H0, bestH0)
            if max(hCompr) > maxH0:
                validatedFit = False
                findH0_fitError = True
            
            if not findH0_fitError:
                results['bestH0'].append(bestH0)
                
                #### Stress-strain computation
                deltaCompr = (maxH0 - hCompr)/1000
                stressCompr = fCompr / (np.pi * (DIAMETER/2000) * deltaCompr)
                strainCompr = deltaCompr / (3*(maxH0/1000))
                # trueStrainCompr = np.log((bestH0)/(hCompr*(hCompr>0)))
                validDelta = (deltaCompr > 0)
                
                results['minStress'].append(np.min(stressCompr))
                results['maxStress'].append(np.max(stressCompr))
                results['minStrain'].append(np.min(strainCompr))
                results['maxStrain'].append(np.max(strainCompr))
                
            elif findH0_fitError:
                results['bestH0'].append(np.nan)
                
                results['minStress'].append(np.nan)
                results['maxStress'].append(np.nan)
                results['minStrain'].append(np.nan)
                results['maxStrain'].append(np.nan)
                

            
            


            #### (4.1) Fits on specific regions of the curve
            
            list_strainPredict_fitToPlot = [[] for kk in range(len(fit_toPlot))]
            
            dictRegionFit = {'regionFitNames' : [], 'K' : [], 'R2' : [],  'K2' : [], 'H0' : [],
                             'fitError' : [], 'validatedFit' : [], 'Npts' : [], 'K_CIW' : []} 
            # 'E' : [], 'H0' : []
            
            #### SETTING ! Setting of the region fits
           
            
            # 

            if not findH0_fitError:
                
                # fitConditions = []
                # for ii in range(len(fit_intervals)-1):
                #     for jj in range(ii+1, len(fit_intervals)):
                #         fitConditions.append((stressCompr > fit_intervals[ii]) & (stressCompr < fit_intervals[jj]))
                
                #### >>> OPTION TO LIGHTEN THE COMPUTATION
                fitConditions = []
                for kk in range(len(regionFitsNames)):
                    ftP = regionFitsNames[kk]
                    lowS, highS = int(fitMin[kk]), int(fitMax[kk])
                    fitConditions.append((stressCompr > lowS) & (stressCompr < highS))
                
                N_Fits = len(fitConditions)
                
                for ii in range(N_Fits):
                    regionFitName = regionFitsNames[ii]
                    mask_region = fitConditions[ii]
                    Npts_region = np.sum(mask_region)
                    
                    if Npts_region > 5:
                        fCompr_region = fCompr[mask_region]
                        hCompr_region = hCompr[mask_region]
                        
                        # Vérifier la concordance !
                        K2_region, H0_region, hPredict_region, R2_region, Chi2_region, confIntE_region, confIntH0_region, fitError_region = \
                                      compressionFitChadwick(hCompr_region, fCompr_region, DIAMETER)
                        
                        K_region, strainPredict_region, R2_region, Chi2_region, confIntK_region, fitError_region = \
                            compressionFitChadwick_StressStrain(hCompr_region, fCompr_region, maxH0, DIAMETER)
                            
                        confIntWidthK_region = np.abs(confIntK_region[0] - confIntK_region[1])
                            
                        
                    else:
                        fitError_region = True
                        
                    if (regionFitName in fit_toPlot) and not fitError_region:
                        kk = np.argmax(np.array(fit_toPlot) == regionFitName)
                        list_strainPredict_fitToPlot[kk] = strainPredict_region
                    
                    if not fitError_region:
                        R2CRITERION = dictSelectionCurve['R2']
                        CHI2CRITERION = dictSelectionCurve['Chi2']
                        validatedFit_region = ((R2_region > R2CRITERION) and 
                                                (Chi2_region < CHI2CRITERION))
                    else:
                        validatedFit_region, fitError_region = False, True
                        K_region, confIntK_region, R2_region  = np.nan, [np.nan, np.nan], np.nan
                        K2_region, H0_region = np.nan, np.nan
                        confIntWidthK_region = np.nan
                        # E_region  = np.nan
                        
                    # Fill dictRegionFit (reinitialized for each cell)
                    dictRegionFit['regionFitNames'].append(regionFitName)
                    dictRegionFit['Npts'].append(Npts_region)
                    dictRegionFit['K'].append(K_region)
                    dictRegionFit['K_CIW'].append(confIntWidthK_region)
                    dictRegionFit['R2'].append(R2_region)
                    dictRegionFit['K2'].append(K2_region)
                    dictRegionFit['H0'].append(H0_region)
                    # dictRegionFit['E'].append(E_region)
                    dictRegionFit['fitError'].append(fitError_region)
                    dictRegionFit['validatedFit'].append(validatedFit_region)


                    # Fill results (memorize everything)
                    rFN = regionFitName
                    results['Npts_'+rFN].append(Npts_region)
                    results['KChadwick_'+rFN].append(K_region)
                    results['K_CIW_'+rFN].append(confIntWidthK_region)
                    results['R2Chadwick_'+rFN].append(R2_region)
                    results['K2Chadwick_'+rFN].append(K2_region)
                    results['H0Chadwick_'+rFN].append(H0_region)
                    # results['EChadwick_'+rFN].append(dictRegionFit['E'][ii])
                    results['validatedFit_'+rFN].append(validatedFit_region)
                    
            elif findH0_fitError:
                for rFN in regionFitsNames:
                    dictRegionFit['regionFitNames'].append(rFN)
                    dictRegionFit['Npts'].append(np.nan)
                    dictRegionFit['K'].append(np.nan)
                    dictRegionFit['K_CIW'].append(np.nan)
                    dictRegionFit['R2'].append(np.nan)
                    dictRegionFit['K2'].append(np.nan)
                    dictRegionFit['H0'].append(np.nan)
                    # dictRegionFit['E'].append(E_region)
                    dictRegionFit['fitError'].append(True)
                    dictRegionFit['validatedFit'].append(False)
                                        
                    results['Npts_'+rFN].append(np.nan)
                    results['KChadwick_'+rFN].append(np.nan)
                    results['K_CIW_'+rFN].append(np.nan)
                    results['R2Chadwick_'+rFN].append(np.nan)
                    results['K2Chadwick_'+rFN].append(np.nan)
                    results['H0Chadwick_'+rFN].append(np.nan)
                    # results['EChadwick_'+rFN].append(np.nan)
                    results['validatedFit_'+rFN].append(False)    
            
            for k in dictRegionFit.keys():
                dictRegionFit[k] = np.array(dictRegionFit[k])
            
            
            #### PLOT [2/4]
            # Complete fig 1, 2, 3 with the results of the fit
            if PLOT:
                
                #### fig1
                if not fitError:
                    ax1.plot(thisCompDf['T'].values, thisCompDf['D3'].values-DIAMETER, color = 'chartreuse', linestyle = '-', linewidth = 1.25)
                    
                    # if validatedFit:
                    #     ax1.plot(thisCompDf['T'].values, thisCompDf['D3'].values-DIAMETER, color = 'chartreuse', linestyle = '-', linewidth = 1.25)
                    # else:
                    #     ax1.plot(thisCompDf['T'].values, thisCompDf['D3'].values-DIAMETER, color = 'gold', linestyle = '-', linewidth = 1.25)
                else:
                    ax1.plot(thisCompDf['T'].values, thisCompDf['D3'].values-DIAMETER, color = 'crimson', linestyle = '-', linewidth = 1.25)
                
                # Display jumpD3 >>> DISABLED
                # if jumpD3 != 0:
                #     x = np.mean(thisCompDf['T'].values)
                #     y = np.mean(thisCompDf['D3'].values-DIAMETER) * 1.3
                #     ax1.text(x, y, '{:.2f}'.format(jumpD3), ha = 'center')

                fig1.suptitle(results['cellID'][-1])
                

                #### fig2 & fig3
                colSp = (i-1) % nColsSubplot
                rowSp = (i-1) // nColsSubplot
                # ax2[i-1] with the 1 line plot
                if nRowsSubplot == 1:
                    thisAx2 = ax2[colSp]
                    thisAx3 = ax3[colSp]
                elif nRowsSubplot >= 1:
                    thisAx2 = ax2[rowSp,colSp]
                    thisAx3 = ax3[rowSp,colSp]

                thisAx2.plot(hCompr,fCompr,'b-', linewidth = 0.8)
                thisAx2.plot(hRelax,fRelax,'r-', linewidth = 0.8)
                titleText = results['cellID'][-1] + '__c' + str(i)
                legendText = ''
                thisAx2.set_xlabel('h (nm)')
                thisAx2.set_ylabel('f (pN)')


                if not fitError:
                    legendText += 'H0 = {:.1f}nm\nE = {:.2e}Pa\nR2 = {:.3f}\nChi2 = {:.1f}'.format(H0, E, R2, Chi2)
                    thisAx2.plot(hPredict,fCompr,'k--', linewidth = 0.8, label = legendText, zorder = 2)
                    thisAx2.legend(loc = 'upper right', prop={'size': 6})
                    # if not validatedFit:
                    #     titleText += '\nNON VALIDATED'
                        
                    if not findH0_fitError:
                        thisAx3.plot(stressCompr, strainCompr, 'go', ms = 3)
                        # thisAx3.plot(stressCompr, trueStrainCompr, 'bo', ms = 2)
                        thisAx3.plot([np.percentile(stressCompr,10), np.percentile(stressCompr,90)], 
                                     [np.percentile(stressCompr,10)/E, np.percentile(stressCompr,90)/E],
                                     'k--', linewidth = 1.2, label = legendText)
                        thisAx3.legend(loc = 'lower right', prop={'size': 6})
                    thisAx3.set_xlabel('Stress (Pa)')
                    thisAx3.set_ylabel('Strain')
                    
                else:
                    titleText += '\nFIT ERROR'
                    
                if not findH0_fitError:
                    # Computations to display the fit of the bestH0
                    min_f = np.min(fCompr)
                    low_f = np.linspace(0, min_f, 20)
                    R = DIAMETER/2
                    low_h = bestH0 - ((3*bestH0*low_f)/(np.pi*(findH0_E/1e6)*R))**0.5

                    
                    legendText2 = 'bestH0 = {:.2f}nm'.format(bestH0)
                    plot_startH = np.concatenate((low_h, findH0_hPredict[:]))
                    plot_startF = np.concatenate((low_f, fCompr[:findH0_NbPts]))
                    thisAx2.plot(plot_startH[0], plot_startF[0], ls = '', 
                                  marker = 'o', color = 'skyblue', markersize = 5, label = legendText2)
                    thisAx2.plot(plot_startH, plot_startF, ls = '--', color = 'skyblue', linewidth = 1.2, zorder = 1)
                    thisAx2.legend(loc = 'upper right', prop={'size': 6})
                
                
                multiAxes = [thisAx2, thisAx3]
                
                for ax in multiAxes:
                    ax.title.set_text(titleText)
                    for item in ([ax.title, ax.xaxis.label, \
                                  ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
                        item.set_fontsize(9)
                        
                        
                        
                        
                #### fig4 & fig5
                
                Npts_fitToPlot = dictRegionFit['Npts'][mask_fitToPlot]
                K_fitToPlot = dictRegionFit['K'][mask_fitToPlot]
                K_CIW_fitToPlot =dictRegionFit['K_CIW'][mask_fitToPlot]
                K2_fitToPlot = dictRegionFit['K2'][mask_fitToPlot]
                H0_fitToPlot = dictRegionFit['H0'][mask_fitToPlot]
                R2_fitToPlot = dictRegionFit['R2'][mask_fitToPlot]
                fitError_fitToPlot = dictRegionFit['fitError'][mask_fitToPlot]
                validatedFit_fitToPlot = dictRegionFit['validatedFit'][mask_fitToPlot]
                
                if not findH0_fitError:
                    fitConditions_fitToPlot = np.array(fitConditions)[mask_fitToPlot]
                
                
                # ax2[i-1] with the 1 line plot
                if nRowsSubplot == 1:
                    thisAx4 = ax4[colSp]
                    thisAx5 = ax5[colSp]
                elif nRowsSubplot >= 1:
                    thisAx4 = ax4[rowSp,colSp]
                    thisAx5 = ax5[rowSp,colSp]
                
                main_color = 'k' # colorList10[0]
                
                thisAx4.plot(hCompr,fCompr, color = main_color, marker = 'o', markersize = 2, ls = '', alpha = 0.8) # 'b-'
                # thisAx4.plot(hRelax,fRelax, color = main_color, ls = '-', linewidth = 0.8, alpha = 0.5) # 'r-'
                titleText = results['cellID'][-1] + '__c' + str(i)
                thisAx4.set_xlabel('h (nm)')
                thisAx4.set_ylabel('f (pN)')

                for k in range(len(fit_toPlot)):
                    if not fitError_fitToPlot[k]:

                        fit = fit_toPlot[k]
                        
                        Npts_fit = Npts_fitToPlot[k]
                        K_fit = K_fitToPlot[k]
                        K2_fit = K2_fitToPlot[k]
                        H0_fit = H0_fitToPlot[k]
                        R2_fit = R2_fitToPlot[k]
                        fitError_fit = fitError_fitToPlot[k]
                        validatedFit_fit = validatedFit_fitToPlot[k]
                        fitConditions_fit = fitConditions_fitToPlot[k]
                        
                        R = DIAMETER/2
                        fCompr_fit = fCompr[fitConditions_fit]
                        
                        hPredict_fit2 = H0_fit - ((3*H0_fit*fCompr_fit)/(np.pi*(K2_fit/1e6)*R))**0.5
                        
                        color = colorList30[k]
                        legendText4 = ''
                        
                        # hPredict_fit0 = bestH0 - ((3*bestH0*fCompr_fit)/(np.pi*(K_fit/1e6)*R))**0.5
                        # strainPredict_fit = list_strainPredict_fitToPlot[k]
                        # hPredict_fit = bestH0 * (1 - 3*strainPredict_fit)
                        # legendText4 += 'Range ' + fit + '\n'
                        # legendText4 += 'K = {:.2e}Pa'.format(K_fit)
                        # thisAx4.plot(hPredict_fit, fCompr_fit, color = color, ls = '--', linewidth = 1.8, label = legendText4)
                        
                        if fit not in alreadyLabeled4:
                            alreadyLabeled4.append(fit)
                            legendText4 += '' + fit + ''
                            thisAx4.plot(hPredict_fit2, fCompr_fit, color = color, ls = '-', linewidth = 1.8, label = legendText4)
                        elif fit in alreadyLabeled4:
                            thisAx4.plot(hPredict_fit2, fCompr_fit, color = color, ls = '-', linewidth = 1.8)
                    else:
                        pass
                
                
                
                
                thisAx5.set_xlabel('Strain')
                thisAx5.set_ylabel('Stress (Pa)')
                if not fitError and not findH0_fitError:
                    thisAx5.plot(strainCompr, stressCompr, color = main_color, marker = 'o', markersize = 2, ls = '', alpha = 0.8)
                    
                    for k in range(len(fit_toPlot)):
                        fit = fit_toPlot[k]
                        
                        Npts_fit = Npts_fitToPlot[k]
                        K_fit = K_fitToPlot[k]
                        R2_fit = R2_fitToPlot[k]
                        fitError_fit = fitError_fitToPlot[k]
                        validatedFit_fit = validatedFit_fitToPlot[k]
                        fitConditions_fit = fitConditions_fitToPlot[k]
                        
                        stressCompr_fit = stressCompr[fitConditions_fit]
                        strainPredict_fit = list_strainPredict_fitToPlot[k]
                        
                        color = colorList30[k]
                        legendText5 = ''
                        
                        if not fitError_fit:
                            if fit not in alreadyLabeled5:
                                alreadyLabeled5.append(fit)
                                legendText5 += '' + fit + ''
                                thisAx5.plot(strainPredict_fit, stressCompr_fit,  
                                             color = color, ls = '-', linewidth = 1.8, label = legendText5)
                            elif fit in alreadyLabeled5:
                                thisAx5.plot(strainPredict_fit, stressCompr_fit,  
                                             color = color, ls = '-', linewidth = 1.8)
                        else:
                            pass

                multiAxes = [thisAx4, thisAx5]
                
                for ax in multiAxes:
                    ax.title.set_text(titleText)
                    for item in ([ax.title, ax.xaxis.label, \
                                  ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
                        item.set_fontsize(9)
                        
                        
                #### fig6
                
                fitCentersPlot = fitCenters[mask_fitToPlot]
                
                if nRowsSubplot == 1:
                    thisAx6 = ax6[colSp]
                elif nRowsSubplot >= 1:
                    thisAx6 = ax6[rowSp,colSp]
                
                
                thisAx6.set_xlabel('sigma (Pa)')
                thisAx6.set_xlim([0, 1200])
                thisAx6.set_ylabel('K (kPa)')
                
                relativeError = np.zeros(len(K_fitToPlot))
                if not fitError and not findH0_fitError:
                    
                    for k in range(len(fit_toPlot)):
                        fit = fit_toPlot[k]
                        
                        K_fit = K_fitToPlot[k]
                        K_CIW_fit = K_CIW_fitToPlot[k]
                        fitError_fit = fitError_fitToPlot[k]
                        validatedFit_fit = validatedFit_fitToPlot[k]
                        fitConditions_fit = fitConditions_fitToPlot[k]
                        
                        stressCompr_fit = stressCompr[fitConditions_fit]
                        strainPredict_fit = list_strainPredict_fitToPlot[k]
    
                        color = colorList30[k]
                        
                        if not fitError_fit:
                            
                            Err = K_CIW_fit
                            relativeError[k] = (Err/K_fit)
                            mec = None
                            thisAx6.errorbar([fitCentersPlot[k]], [K_fit/1000], yerr = [(Err/2)/1000],
                                          color = color, marker = 'o', ms = 5, mec = mec)                           
                            
                        
                    multiAxes = [thisAx6]
                    
                    for ax in multiAxes:
                        ax.title.set_text(titleText)
                        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + \
                                     ax.get_xticklabels() + ax.get_yticklabels()):
                            item.set_fontsize(9)
                            
                
                        
                #### fig7
                if plotSmallElements:

                    if nRowsSubplot == 1:
                        thisAx7 = ax7[colSp]
                    elif nRowsSubplot >= 1:
                        thisAx7 = ax7[rowSp,colSp]
                        
                    thisAx7.set_xlabel('epsilon')
                    thisAx7.set_ylabel('ratios')
                    legendText7 = ''
                    
                    # def def2delta(e):
                    #     d = 3*bestH0*e
                    #     return(d)
                    
                    # def delta2def(d):
                    #     e = d/(3*bestH0)
                    #     return(e)
                    
                    # secax = thisAx7.secondary_xaxis('top', functions=(def2delta, delta2def))
                    # secax.set_xlabel('delta (nm)')
                    
                    if not fitError and not findH0_fitError:
                        
                        A = 2* ((deltaCompr*(DIAMETER/2000))**0.5)
                        largeX = A/(maxH0/1000)
                        smallX = deltaCompr/(DIAMETER/1000)
                        # legendText7 = 'bestH0 = {:.2f}nm'.format(bestH0)
                        
                        thisAx7.plot(strainCompr, largeX, color = 'red',     ls = '', marker = '+', label = 'a/H0',     markersize = 3)#, mec = 'k', mew = 0.5)
                        thisAx7.plot(strainCompr, smallX, color = 'skyblue', ls = '', marker = '+', label = 'delta/2R', markersize = 3)#, mec = 'k', mew = 0.5)
                        thisAx7.legend(loc = 'upper left', prop={'size': 5})
                        thisAx7.set_yscale('log')
                        thisAx7.set_ylim([5e-3,10])
                        minPlot, maxPlot = thisAx7.get_xlim()
                        thisAx7.set_xlim([0,maxPlot])
                        thisAx7.plot([0,maxPlot], [1, 1], color = 'k', ls = '--', lw = 0.5)
                        
                        # # thisAx7bis.tick_params(axis='y', labelcolor='b')
                        # thisAx7bis = thisAx7.twinx()
                        # color = 'firebrick'
                        # thisAx7bis.set_ylabel('F (pN)', color=color)
                        # thisAx7bis.plot(strainCompr, fCompr, color=color, lw = 0.5)
                        # thisAx7bis.set_ylim([0,1400])
                        # thisAx7bis.tick_params(axis='y', labelrotation = 50, labelsize = 10)
                        # # thisAx7bis.tick_params(axis='y', labelcolor=color)
                        # # thisAx7bis.set_yticks([0,500,1000,1500])
                        # # minh = np.min(tsDF['D3'].values-DIAMETER)
                        
                        
                        epsLim = (largeX < 1.0)
                        if len(strainCompr[epsLim]) > 0:
                            strainLimit = np.max(strainCompr[epsLim])
                            minPlot, maxPlot = thisAx5.get_ylim()
                            # thisAx5.plot([strainLimit, strainLimit], [minPlot, maxPlot], color = 'gold', ls = '--')
                            minPlot, maxPlot = thisAx7.get_ylim()
                            thisAx7.plot([strainLimit, strainLimit], [minPlot, maxPlot], color = 'gold', ls = '--')
                        
                                
                        multiAxes = [thisAx7] #, thisAx7bis]
                        
                        for ax in multiAxes:
                            ax.title.set_text(titleText + '\nbestH0 = {:.2f}nm'.format(maxH0))
                            for item in ([ax.title, ax.xaxis.label, \
                                          ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
                                item.set_fontsize(9)


            #### (5) hysteresis (its definition may change)
            try:
                results['hysteresis'].append(hCompr[0] - hRelax[-1])
            except:
                print('hysteresis computation failed, see hCompr and hRelax below:')
                print(hCompr, hRelax)
                print('Details of the delimitation of the curve')
                print(jStart, jMax, jStop, offsetStop, thisCompDf.shape[0])
                print(offsetStart2)
                print(thisCompDf.B.values)
                results['hysteresis'].append(np.nan)
        
        #### (6) Deal with the non analysed compressions
        else: # The compression curve was detected as not suitable for analysis
            generalBug = True
            validatedThickness = False
            results['initialThickness'].append(np.nan)
            results['minThickness'].append(np.nan)
            results['maxIndent'].append(np.nan)
            results['previousThickness'].append(np.nan)
            results['surroundingThickness'].append(np.nan)
            results['surroundingDx'].append(np.nan)
            results['surroundingDz'].append(np.nan)
            results['ctFieldDX'].append(np.nan)
            results['ctFieldDZ'].append(np.nan)
            results['ctFieldThickness'].append(np.nan)
            results['ctFieldFluctuAmpli'].append(np.nan)
            results['jumpD3'].append(np.nan) 
            results['validatedThickness'].append(validatedThickness)
            validatedFit = False
            results['critFit'].append('Not relevant')
            results['bestH0'].append(np.nan)
            results['H0Chadwick'].append(np.nan)
            results['EChadwick'].append(np.nan)
            results['R2Chadwick'].append(np.nan)
            results['EChadwick_CIWidth'].append(np.nan)
            results['validatedFit'].append(validatedFit)

            results['minForce'].append(np.nan)
            results['maxForce'].append(np.nan)
            results['minStress'].append(np.nan)
            results['maxStress'].append(np.nan)
            results['minStrain'].append(np.nan)
            results['maxStrain'].append(np.nan)

            results['comments'].append('Unspecified bug in the code')
            results['hysteresis'].append(np.nan)

            
            for rFN in regionFitsNames:
                results['Npts_'+rFN].append(np.nan)
                results['KChadwick_'+rFN].append(np.nan)
                results['K_CIW_'+rFN].append(np.nan)
                results['R2Chadwick_'+rFN].append(np.nan)
                results['K2Chadwick_'+rFN].append(np.nan)
                results['H0Chadwick_'+rFN].append(np.nan)
                # results['EChadwick_'+rFN].append(np.nan)
                results['validatedFit_'+rFN].append(False)

        if not doThisCompAnalysis:
            print('Curve not suitable for analysis !')
            print(results['cellID'][-1])
            print('Compression no ' + str(i))

    
    #### PLOT [3/4]
    
    if PLOT:
        #### fig1
        color = 'blue'
        ax1.plot(tsDF['T'].values, tsDF['D3'].values-DIAMETER, color = color, ls = '-', linewidth = 1, zorder = 1)
        (axm, axM) = ax1.get_ylim()
        ax1.set_ylim([min(0,axm), axM])
        if (max(tsDF['D3'].values-DIAMETER) > 200):
            ax1.set_yticks(np.arange(0, max(tsDF['D3'].values-DIAMETER), 100))
        
        twinAxis = True
        if twinAxis:
            ax1.tick_params(axis='y', labelcolor='b')
            ax1bis = ax1.twinx()
            color = 'firebrick'
            ax1bis.set_ylabel('F (pN)', color=color)
            ax1bis.plot(tsDF['T'].values, tsDF['F'].values, color=color)
            ax1bis.tick_params(axis='y', labelcolor=color)
            ax1bis.set_yticks([0,500,1000,1500])
            minh = np.min(tsDF['D3'].values-DIAMETER)
            ratio = min(1/abs(minh/axM), 5)
#             print(ratio)
            (axmbis, axMbis) = ax1bis.get_ylim()
            ax1bis.set_ylim([0, max(axMbis*ratio, 3*max(tsDF['F'].values))])
        
        
        #### fig3
        # Rescale fig3 axes
        eMin, eMax = 1, 0
        sMin, sMax = 1000, 0
        for i in range(1, Ncomp+1):
            colSp = (i-1) % nColsSubplot
            rowSp = (i-1) // nColsSubplot
            # ax2[i-1] with the 1 line plot
            if nRowsSubplot == 1:
                thisAx3 = ax3[colSp]
            elif nRowsSubplot >= 1:
                thisAx3 = ax3[rowSp,colSp]
                
            if thisAx3.get_ylim()[1] > eMax:
                eMax = thisAx3.get_ylim()[1]
            if thisAx3.get_xlim()[1] > sMax:
                sMax = thisAx3.get_xlim()[1]
                
        for i in range(1, Ncomp+1):
            colSp = (i-1) % nColsSubplot
            rowSp = (i-1) // nColsSubplot
            # ax2[i-1] with the 1 line plot
            if nRowsSubplot == 1:
                thisAx3 = ax3[colSp]
            elif nRowsSubplot >= 1:
                thisAx3 = ax3[rowSp,colSp]
                
            sMax = min(sMax, 2100)
            eMax = min(eMax, 0.5)
            
            thisAx3.set_xlim([0, sMax])
            thisAx3.set_ylim([0, eMax])
                
        #### fig4
        NL = len(alreadyLabeled4)
        titleVoid = ' '
        if NL > 16:
            titleVoid += '\n '
        fig4.suptitle(titleVoid)
        fig4.legend(loc='upper center', bbox_to_anchor=(0.5,1), ncol = min(8, NL))
                
        #### fig5
        
        NL = len(alreadyLabeled5)
        titleVoid = ' '
        if NL > 16:
            titleVoid += '\n '
        fig5.suptitle(titleVoid)
        fig5.legend(loc='upper center', bbox_to_anchor=(0.5,1), ncol = min(8, NL))
        
        # Rescale fig5 axes
        eMin, eMax = 1, 0
        sMin, sMax = 1000, 0
        for i in range(1, Ncomp+1):
            colSp = (i-1) % nColsSubplot
            rowSp = (i-1) // nColsSubplot
            # ax2[i-1] with the 1 line plot
            if nRowsSubplot == 1:
                thisAx5 = ax5[colSp]
            elif nRowsSubplot >= 1:
                thisAx5 = ax5[rowSp,colSp]
                
            if thisAx5.get_xlim()[1] > eMax:
                eMax = thisAx5.get_xlim()[1]
            if thisAx5.get_ylim()[1] > sMax:
                sMax = thisAx5.get_ylim()[1]
                
        for i in range(1, Ncomp+1):
            colSp = (i-1) % nColsSubplot
            rowSp = (i-1) // nColsSubplot
            # ax2[i-1] with the 1 line plot
            if nRowsSubplot == 1:
                thisAx5 = ax5[colSp]
            elif nRowsSubplot >= 1:
                thisAx5 = ax5[rowSp,colSp]
            
            sMax = min(sMax, 2100)
            eMax = min(eMax, 0.5)
            thisAx5.set_ylim([0, sMax])
            thisAx5.set_xlim([0, eMax])
                
        #### fig6
        # Rescale fig6 axes
        
        sMax = 0
        Kmax = 0

        for i in range(1, Ncomp+1):
            colSp = (i-1) % nColsSubplot
            rowSp = (i-1) // nColsSubplot
            if nRowsSubplot == 1:
                thisAx6 = ax6[colSp]
            elif nRowsSubplot >= 1:
                thisAx6 = ax6[rowSp,colSp]

            if thisAx6.get_xlim()[1] > sMax:
                sMax = thisAx6.get_xlim()[1]
            if thisAx6.get_ylim()[1] > Kmax:
                Kmax = thisAx6.get_xlim()[1]
        
        for i in range(1, Ncomp+1):
            colSp = (i-1) % nColsSubplot
            rowSp = (i-1) // nColsSubplot
            if nRowsSubplot == 1:
                thisAx6 = ax6[colSp]
            elif nRowsSubplot >= 1:
                thisAx6 = ax6[rowSp,colSp]
                
            Kmax = min(Kmax, 20)

            thisAx6.set_xlim([0, sMax])
            thisAx6.set_ylim([0, Kmax])
                
                
        Allfigs = [fig1,fig2,fig3,fig4,fig5,fig6,fig7]
        
        for fig in Allfigs:
            fig.tight_layout()

    
    
    #### PLOT [4/4]
    # Save the figures
    if PLOT:
        dpi1 = 150
        dpi2 = 150
        # figDir = todayFigDir # Already by default
        # figSubDir = 'MecaAnalysis_allCells'
        # archiveFig(fig1, ax1, name=results['cellID'][-1] + '_01_h(t)', dpi = dpi1)
        # archiveFig(fig2, ax2, name=results['cellID'][-1] + '_02_F(h)', dpi = dpi1)
        # archiveFig(fig3, ax3, name=results['cellID'][-1] + '_03_sig(eps)', dpi = dpi1)
        # archiveFig(fig4, ax4, name=results['cellID'][-1] + '_04_F(h)_regionFits', dpi = dpi1)
        # archiveFig(fig5, ax5, name=results['cellID'][-1] + '_05_sig(eps)_regionFits', dpi = dpi1)
        # archiveFig(fig6, ax6, name=results['cellID'][-1] + '_06_K(s)', dpi = dpi1)
        # archiveFig(fig7, ax7, name=results['cellID'][-1] + '_07_smallElements', dpi = dpi1)
        

        figDir = os.path.join(todayFigDirLocal, 'MecaAnalysis_allCells_' + task)
        jvu.archiveFig(fig1, ax1, figDir, name=results['cellID'][-1] + '_01_h(t)', dpi = dpi1)
        jvu.archiveFig(fig2, ax2, figDir, name=results['cellID'][-1] + '_02_F(h)', dpi = dpi1)
        jvu.archiveFig(fig3, ax3, figDir, name=results['cellID'][-1] + '_03_sig(eps)', dpi = dpi1)
        jvu.archiveFig(fig4, ax4, figDir, name=results['cellID'][-1] + '_04_F(h)_regionFits', dpi = dpi1)
        jvu.archiveFig(fig5, ax5, figDir, name=results['cellID'][-1] + '_05_sig(eps)_regionFits', dpi = dpi1)
        jvu.archiveFig(fig6, ax6, figDir, name=results['cellID'][-1] + '_06_K(s)', dpi = dpi1)
        jvu.archiveFig(fig7, ax7, figDir, name=results['cellID'][-1] + '_07_smallElements', dpi = dpi1)
        

        figDir = os.path.join(ownCloudTodayFigDir, 'MecaAnalysis_allCells_' + task)
        jvu.archiveFig(fig1, ax1, figDir, name=results['cellID'][-1] + '_01_h(t)', dpi = dpi2)
        jvu.archiveFig(fig2, ax2, figDir, name=results['cellID'][-1] + '_02_F(h)', dpi = dpi2)
        jvu.archiveFig(fig3, ax3, figDir, name=results['cellID'][-1] + '_03_sig(eps)', dpi = dpi2)
        jvu.archiveFig(fig4, ax4, figDir, name=results['cellID'][-1] + '_04_F(h)_regionFits', dpi = dpi2)
        jvu.archiveFig(fig5, ax5, figDir, name=results['cellID'][-1] + '_05_sig(eps)_regionFits', dpi = dpi2)
        jvu.archiveFig(fig6, ax6, figDir, name=results['cellID'][-1] + '_06_K(s)', dpi = dpi2)
        jvu.archiveFig(fig7, ax7, figDir, name=results['cellID'][-1] + '_07_smallElements', dpi = dpi2)
        
        
        if PLOT_SHOW:
            Allfigs = [fig1,fig2,fig3,fig4,fig5,fig6,fig7]
            # for fig in Allfigs:
                # fig.tight_layout()
                # fig.show()
            plt.show()
        else:
            plt.close('all')

    return(results)



def createDataDict_meca(list_mecaFiles, listColumnsMeca, task, PLOT):
    """
    Subfunction of computeGlobalTable_meca
    Create the dictionnary that will be converted in a pandas table in the end.
    """
    expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
    tableDict = {}
    Nfiles = len(list_mecaFiles)
    PLOT_SHOW = (Nfiles==1)
    PLOT_SHOW = 0
    if not PLOT_SHOW:
        plt.ioff()
    for c in listColumnsMeca:
        tableDict[c] = []
    for f in list_mecaFiles: #[:10]:
        tS_DataFilePath = os.path.join(timeSeriesDataDir, f)
        current_tsDF = pd.read_csv(tS_DataFilePath, sep = ';')
         # MAIN SUBFUNCTION
        current_resultDict = analyseTimeSeries_meca(f, current_tsDF, expDf, 
                                                    listColumnsMeca, task,
                                                    PLOT, PLOT_SHOW)
        for k in current_resultDict.keys():
            tableDict[k] += current_resultDict[k]
#     for k in tableDict.keys():
#         print(k, len(tableDict[k]))
    return(tableDict)


def update_uiDf(ui_fileName, mecaDf):
    """
    """
    listColumnsUI = ['date','cellName','cellID','manipID','compNum',
                     'UI_Valid','UI_Comments']

    try:
        print('Imported existing UI table')
        savePath = os.path.join(dataDir, (ui_fileName + '.csv'))
        uiDf = pd.read_csv(savePath, sep='\t', )
        fromScratch = False
    except:
        print('No existing UI table found with this name')
        fromScratch = True

    new_uiDf = mecaDf[listColumnsUI[:5]]
    if not fromScratch:
        existingCellId = uiDf['cellID'].values
        new_uiDf = new_uiDf.loc[new_uiDf['cellID'].apply(lambda x : x not in existingCellId), :]
    
    nrows = new_uiDf.shape[0]
    new_uiDf['UI_Valid'] = np.ones(nrows, dtype = bool)
    new_uiDf['UI_Comments'] = np.array(['' for i in range(nrows)])
    
    if not fromScratch:
        new_uiDf = pd.concat([uiDf, new_uiDf], axis = 0, ignore_index=True)
        
    savePath = os.path.join(dataDir, (ui_fileName + '.csv'))
    new_uiDf.sort_values(by=['cellID', 'compNum'], inplace = True)
    new_uiDf.to_csv(savePath, sep='\t', index = False)



def computeGlobalTable_meca(task = 'fromScratch', fileName = 'Global_MecaData', 
                            save = False, PLOT = False, \
                            source = 'Matlab', listColumnsMeca=listColumnsMeca,
                            ui_fileName = 'UserManualSelection_MecaData'):
    """
    Compute the GlobalTable_meca from the time series data files.
    Option task='fromScratch' will analyse all the time series data files and construct a new GlobalTable from them regardless of the existing GlobalTable.
    Option task='updateExisting' will open the existing GlobalTable and determine which of the time series data files are new ones, and will append the existing GlobalTable with the data analysed from those new fils.
    listColumnsMeca have to contain all the fields of the table that will be constructed.
    """
    top = time.time()
    
#     list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
#                       if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
#                       and ('R40' in f))] # Change to allow different formats in the future
    
    suffixPython = '_PY'
    if source == 'Matlab':
        list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
                      if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                      and (('R40' in f) or ('L40' in f)) and not (suffixPython in f))]
        
    elif source == 'Python':
        list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
                      if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                      and (('R40' in f) or ('L40' in f)) and (suffixPython in f))]
        # print(list_mecaFiles)
    
#     print(list_mecaFiles)
    
    if task == 'fromScratch':
        # create a dict containing the data
        tableDict = createDataDict_meca(list_mecaFiles, listColumnsMeca, task, PLOT) # MAIN SUBFUNCTION
        
        # create the dataframe from it
        mecaDf = pd.DataFrame(tableDict)
        
        update_uiDf(ui_fileName, mecaDf)
        
        # last step: now that the dataFrame is complete, one can use "compStartTimeThisDay" col to compute the start time of each compression relative to the first one done this day.
        allDates = list(mecaDf['date'].unique())
        for d in allDates:
            subDf = mecaDf.loc[mecaDf['date'] == d]
            experimentStartTime = np.min(subDf['compStartTimeThisDay'])
            mecaDf['compStartTimeThisDay'].loc[mecaDf['date'] == d] = mecaDf['compStartTimeThisDay'] - experimentStartTime
        
    elif task == 'updateExisting':
        # get existing table
        try:
            savePath = os.path.join(dataDir, (fileName + '.csv'))
            existing_mecaDf = pd.read_csv(savePath, sep=';')
        except:
            print('No existing table found')
            
        # find which of the time series files are new
        new_list_mecaFiles = []
        for f in list_mecaFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            if currentCellID not in existing_mecaDf.cellID.values:
                new_list_mecaFiles.append(f)
                
        # create the dict with new data
        new_tableDict = createDataDict_meca(new_list_mecaFiles, listColumnsMeca, task, PLOT) # MAIN SUBFUNCTION
        # create the dataframe from it
        new_mecaDf = pd.DataFrame(new_tableDict)
        # fuse the existing table with the new one
        mecaDf = pd.concat([existing_mecaDf, new_mecaDf])
        
        update_uiDf(ui_fileName, mecaDf)
        
    else: # If task is neither 'fromScratch' nor 'updateExisting'
    # Then task can be a substring that can be in some timeSeries file !
    # It will create a table with only these files !
        # save = False
        task_list = task.split(' & ')
        new_list_mecaFiles = []
        for f in list_mecaFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            for t in task_list:
                if t in currentCellID:
                    new_list_mecaFiles.append(f)
                    break
        # create the dict with new data
        new_tableDict = createDataDict_meca(new_list_mecaFiles, listColumnsMeca, task, PLOT) # MAIN SUBFUNCTION
        # create the dataframe from it
        mecaDf = pd.DataFrame(new_tableDict)
        
        update_uiDf(ui_fileName, mecaDf)
    
    for c in mecaDf.columns:
            if 'Unnamed' in c:
                mecaDf = mecaDf.drop([c], axis=1)
    
    if save:
        saveName = fileName + '.csv'
        savePath = os.path.join(dataDir, saveName)
        mecaDf.to_csv(savePath, sep=';')
    
    delta = time.time() - top
    print(delta)
    
    return(mecaDf)
            

    
def getGlobalTable_meca(fileName):
    try:
        savePath = os.path.join(dataDir, (fileName + '.csv'))
        mecaDf = pd.read_csv(savePath, sep=';')
        print('Extracted a table with ' + str(mecaDf.shape[0]) + ' lines and ' + str(mecaDf.shape[1]) + ' columns.')
    except:
        print('No existing table found')
        
    for c in mecaDf.columns:
        if 'Unnamed' in c:
            mecaDf = mecaDf.drop([c], axis=1)
        # if 'K_CIW_' in c:    
        #     mecaDf[c].apply(lambda x : x.strip('][').split(', ')).apply(lambda x : [float(x[0]), float(x[1])])
    
    if 'ExpDay' in mecaDf.columns:
        dateExemple = mecaDf.loc[mecaDf.index[0],'ExpDay']
        if not ('manipID' in mecaDf.columns):
            mecaDf['manipID'] = mecaDf['ExpDay'] + '_' + mecaDf['CellID'].apply(lambda x: x.split('_')[0])
            
    elif 'date' in mecaDf.columns:
        dateExemple = mecaDf.loc[mecaDf.index[0],'date']
        if re.match(dateFormatExcel, dateExemple):
            print('bad date')
        
    if not ('manipID' in mecaDf.columns):
        mecaDf['manipID'] = mecaDf['date'] + '_' + mecaDf['cellName'].apply(lambda x: x.split('_')[0])

        
    return(mecaDf)

# %%% (4.4) Fluorescence data

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

# %%% (4.5) Oscillations

def analyseTimeSeries_sinus(f, tsDF, expDf, listColumns, PLOT, PLOT_SHOW):
    
    plotSmallElements = True
    
    #### (0) Import experimental infos
    split_f = f.split('_')
    tsDF.dx, tsDF.dy, tsDF.dz, tsDF.D2, tsDF.D3 = tsDF.dx*1000, tsDF.dy*1000, tsDF.dz*1000, tsDF.D2*1000, tsDF.D3*1000
    thisManipID = jvu.findInfosInFileName(f, 'manipID')
    thisExpDf = expDf.loc[expDf['manipID'] == thisManipID]

    # Deal with the asymmetric pair case : the diameter can be for instance 4503 (float) or '4503_2691' (string)
    diameters = thisExpDf.at[thisExpDf.index.values[0], 'bead diameter'].split('_')
    if len(diameters) == 2:
        DIAMETER = (int(diameters[0]) + int(diameters[1]))/2.
    else:
        DIAMETER = int(diameters[0])
    
    EXPTYPE = str(thisExpDf.at[thisExpDf.index.values[0], 'experimentType'])

    results = {}
    for c in listColumnsMeca:
        results[c] = []

    nUplet = thisExpDf.at[thisExpDf.index.values[0], 'normal field multi images']
    
    if 'sinus' in EXPTYPE:
        pass


def createDataDict_sinus(listFiles, listColumns, PLOT):
    """
    Subfunction of computeGlobalTable_meca
    Create the dictionnary that will be converted in a pandas table in the end.
    """
    expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
    tableDict = {}
    Nfiles = len(listFiles)
    PLOT_SHOW = (Nfiles==1)
    PLOT_SHOW = 0
    if not PLOT_SHOW:
        plt.ioff()
    for c in listColumns:
        tableDict[c] = []
    for f in listFiles: #[:10]:
        tS_DataFilePath = os.path.join(timeSeriesDataDir, f)
        current_tsDF = pd.read_csv(tS_DataFilePath, sep = ';')
         # MAIN SUBFUNCTION
        current_resultDict = analyseTimeSeries_sinus(f, current_tsDF, expDf, 
                                                    listColumns, PLOT, PLOT_SHOW)
        for k in current_resultDict.keys():
            tableDict[k] += current_resultDict[k]
#     for k in tableDict.keys():
#         print(k, len(tableDict[k]))
    return(tableDict)

def computeGlobalTable_sinus(task = 'fromScratch', fileName = 'Global_Sinus', save = False, PLOT = False, \
                            source = 'Matlab', listColumns=listColumnsMeca):
    """
    Compute the GlobalTable_Sinus from the time series data files.
    Option task='fromScratch' will analyse all the time series data files and construct a new GlobalTable from them regardless of the existing GlobalTable.
    Option task='updateExisting' will open the existing GlobalTable and determine which of the time series data files are new ones, and will append the existing GlobalTable with the data analysed from those new fils.
    listColumns have to contain all the fields of the table that will be constructed.
    """
    top = time.time()
    
#     list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
#                       if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
#                       and ('R40' in f))] # Change to allow different formats in the future
    
    suffixPython = '_PY'
    if source == 'Matlab':
        list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
                      if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                      and ('R40' in f) and not (suffixPython in f))]
        
    elif source == 'Python':
        list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
                      if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                      and ('sin' in f) and (suffixPython in f))]
        # print(list_mecaFiles)
    
#     print(list_mecaFiles)
    
    if task == 'fromScratch':
        # create a dict containing the data
        tableDict = createDataDict_sinus(list_mecaFiles, listColumns, PLOT) # MAIN SUBFUNCTION
        # create the dataframe from it
        DF = pd.DataFrame(tableDict)
        
        # last step: now that the dataFrame is complete, one can use "compStartTimeThisDay" col to compute the start time of each compression relative to the first one done this day.
        allDates = list(DF['date'].unique())
        for d in allDates:
            subDf = DF.loc[DF['date'] == d]
            experimentStartTime = np.min(subDf['compStartTimeThisDay'])
            DF['compStartTimeThisDay'].loc[DF['date'] == d] = DF['compStartTimeThisDay'] - experimentStartTime
        
    elif task == 'updateExisting':
        # get existing table
        try:
            savePath = os.path.join(dataDir, (fileName + '.csv'))
            existing_mecaDf = pd.read_csv(savePath, sep=';')
        except:
            print('No existing table found')
            
        # find which of the time series files are new
        new_list_mecaFiles = []
        for f in list_mecaFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            if currentCellID not in existing_mecaDf.cellID.values:
                new_list_mecaFiles.append(f)
                
        # create the dict with new data
        new_tableDict = createDataDict_meca(new_list_mecaFiles, listColumnsMeca, PLOT) # MAIN SUBFUNCTION
        # create the dataframe from it
        new_mecaDf = pd.DataFrame(new_tableDict)
        # fuse the existing table with the new one
        DF = pd.concat([existing_mecaDf, new_mecaDf])
        
    else: # If task is neither 'fromScratch' nor 'updateExisting'
    # Then task can be a substring that can be in some timeSeries file !
    # It will create a table with only these files, WITHOUT SAVING IT !
    # But it can plot figs from it.
        # save = False
        task_list = task.split(' & ')
        new_list_mecaFiles = []
        for f in list_mecaFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            for t in task_list:
                if t in currentCellID:
                    new_list_mecaFiles.append(f)
                    break
        # create the dict with new data
        new_tableDict = createDataDict_meca(new_list_mecaFiles, listColumnsMeca, PLOT) # MAIN SUBFUNCTION
        # create the dataframe from it
        DF = pd.DataFrame(new_tableDict)
    
    for c in DF.columns:
            if 'Unnamed' in c:
                DF = DF.drop([c], axis=1)
    
    if save:
        saveName = fileName + '.csv'
        savePath = os.path.join(dataDir, saveName)
        DF.to_csv(savePath, sep=';')
    
    delta = time.time() - top
    print(delta)
    
    return(DF)



# %% (5) General import functions

# %%% Utility functions

def removeColumnsDuplicate(df):
    cols = df.columns.values
    for c in cols:
        if c.endswith('_x'):
            df = df.rename(columns={c: c[:-2]})
        elif c.endswith('_y'):
            df = df.drop(columns=[c])
    return(df)
    
# %%% Main function

def getGlobalTable(kind, experimentalDataDir = experimentalDataDir):
    if kind == 'ctField':
        GlobalTable = getGlobalTable_ctField()
        expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
        fluoDf = getFluoData()
        GlobalTable = pd.merge(expDf, GlobalTable, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = pd.merge(GlobalTable, fluoDf, how="left", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = removeColumnsDuplicate(GlobalTable)
        print('Merged table has ' + str(GlobalTable.shape[0]) + ' lines and ' + str(GlobalTable.shape[1]) + ' columns.')
        
        # print(GlobalTable_ctField.head())
    
    elif kind == 'ctField_py':
        GlobalTable = getGlobalTable_ctField('Global_CtFieldData_Py')
        expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
        fluoDf = getFluoData()
        GlobalTable = pd.merge(expDf, GlobalTable, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = pd.merge(GlobalTable, fluoDf, how="left", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = removeColumnsDuplicate(GlobalTable)
        print('Merged table has ' + str(GlobalTable.shape[0]) + ' lines and ' + str(GlobalTable.shape[1]) + ' columns.')
        
        # print(GlobalTable_ctField.head())
        
        # return(GlobalTable_ctField)

    elif kind == 'meca_matlab':
        GlobalTable = getGlobalTable_meca('Global_MecaData')
        expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
        fluoDf = getFluoData()
        GlobalTable = pd.merge(GlobalTable, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = pd.merge(GlobalTable, fluoDf, how="left", left_on='CellName', right_on='cellID'
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable.shape[0]) + ' lines and ' + str(GlobalTable.shape[1]) + ' columns.')
        
        # print(GlobalTable.tail())
        GlobalTable = removeColumnsDuplicate(GlobalTable)
        # return(GlobalTable_meca_Matlab)


    elif kind == 'meca_py':
        GlobalTable = getGlobalTable_meca('Global_MecaData_Py')
        expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
        fluoDf = getFluoData()
        GlobalTable = pd.merge(GlobalTable, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = pd.merge(GlobalTable, fluoDf, how="left", left_on='cellID', right_on='cellID'
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable.shape[0]) + ' lines and ' + str(GlobalTable.shape[1]) + ' columns.')
        
        # print(GlobalTable_meca_Py.tail())
        GlobalTable = removeColumnsDuplicate(GlobalTable)
        # return(GlobalTable_meca_Py)


    elif kind == 'meca_py2':
        GlobalTable = getGlobalTable_meca('Global_MecaData_Py2')
        expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
        fluoDf = getFluoData()
        GlobalTable = pd.merge(GlobalTable, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = pd.merge(GlobalTable, fluoDf, how="left", left_on='cellID', right_on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable.shape[0]) + ' lines and ' + str(GlobalTable.shape[1]) + ' columns.')

        # print(GlobalTable_meca_Py2.tail())
        GlobalTable = removeColumnsDuplicate(GlobalTable)
        # return(GlobalTable_meca_Py2)
    
    elif kind == 'meca_nonLin':
        GlobalTable = getGlobalTable_meca('Global_MecaData_NonLin_Py')
        expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
        fluoDf = getFluoData()
        GlobalTable = pd.merge(GlobalTable, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = pd.merge(GlobalTable, fluoDf, how="left", left_on='cellID', right_on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable.shape[0]) + ' lines and ' + str(GlobalTable.shape[1]) + ' columns.')

        # print(GlobalTable_meca_nonLin.tail())
        GlobalTable = removeColumnsDuplicate(GlobalTable)
        # return(GlobalTable_meca_nonLin)
    
    elif kind == 'meca_MCA':
        GlobalTable = getGlobalTable_meca('Global_MecaData_MCA')
        expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
        fluoDf = getFluoData()
        GlobalTable = pd.merge(GlobalTable, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = pd.merge(GlobalTable, fluoDf, how="left", left_on='cellID', right_on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable.shape[0]) + ' lines and ' + str(GlobalTable.shape[1]) + ' columns.')

        # print(GlobalTable_meca_nonLin.tail())
        GlobalTable = removeColumnsDuplicate(GlobalTable)
        # return(GlobalTable)
    
    else:
        GlobalTable = getGlobalTable_meca(kind)
        expDf = jvu.getExperimentalConditions(experimentalDataDir, suffix = '_JV')
        fluoDf = getFluoData()
        GlobalTable = pd.merge(GlobalTable, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable = pd.merge(GlobalTable, fluoDf, how="left", left_on='cellID', right_on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable.shape[0]) + ' lines and ' + str(GlobalTable.shape[1]) + ' columns.')
    
        # print(GlobalTable_meca_nonLin.tail())
        GlobalTable = removeColumnsDuplicate(GlobalTable)
    
    if 'substrate' in GlobalTable.columns:
        vals_substrate = GlobalTable['substrate'].values
        if 'diverse fibronectin discs' in vals_substrate:
            try:
                cellIDs = GlobalTable[GlobalTable['substrate'] == 'diverse fibronectin discs']['cellID'].values
                listFiles = [f for f in os.listdir(timeSeriesDataDir) \
                              if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
                for Id in cellIDs:
                    for f in listFiles:
                        if Id == jvu.findInfosInFileName(f, 'cellID'):
                            thisCellSubstrate = jvu.findInfosInFileName(f, 'substrate')
                            thisCellSubstrate = dictSubstrates[thisCellSubstrate]
                            if not thisCellSubstrate == '':
                                GlobalTable.loc[GlobalTable['cellID'] == Id, 'substrate'] = thisCellSubstrate
                print('Automatic determination of substrate type SUCCEDED !')
                
            except:
                print('Automatic determination of substrate type FAILED !')
    
    return(GlobalTable)
