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
import time
import random
import itertools

from copy import copy
from cycler import cycler
from datetime import date
from scipy.optimize import curve_fit

# %%% Smaller settings

pd.set_option('mode.chained_assignment',None)
pd.set_option('display.max_columns', None)

dateFormatExcel = re.compile('\d{2}/\d{2}/\d{4}')
dateFormatOk = re.compile('\d{2}-\d{2}-\d{2}')

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

my_default_color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
my_default_marker_list = ['o', 's', 'D', '>', '^', 'P', 'X', '<', 'v', 'p']

# %% (1) Directories adress

COMPUTERNAME = os.environ['COMPUTERNAME']
if COMPUTERNAME == 'ORDI-JOSEPH':
    mainDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis"
elif COMPUTERNAME == 'LARISA':
    mainDir = "C://Users//Joseph//Desktop//ActinCortexAnalysis"
elif COMPUTERNAME == '':
    mainDir = "C://Users//josep//Desktop//ActinCortexAnalysis"

# Add the folder to path
import sys
sys.path.append(mainDir + "//Code_Python")
from getExperimentalConditions import getExperimentalConditions


experimentalDataDir = os.path.join(mainDir, "Data_Experimental")
dataDir = os.path.join(mainDir, "Data_Analysis")
figDir = os.path.join(dataDir, "Figures")
todayFigDir = os.path.join(figDir, "Historique//" + str(date.today()))
timeSeriesDataDir = os.path.join(dataDir, "TimeSeriesData")


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

def archiveFig(fig, ax, name='auto', figDir = todayFigDir, figSubDir=''):
    dpi = 300
    if not os.path.exists(figDir):
        os.makedirs(figDir)
    
    saveDir = os.path.join(todayFigDir, figSubDir)
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

    
    return(infoString)
                    
                    
# %% (3) TimeSeries functions

def getCellTimeSeriesData(cellID, fromPython = True):
    if fromPython:
        allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) 
                              if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) 
                                  and f.endswith("PY.csv"))]
    else:
        allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) 
                              if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) 
                                  and f.endswith("PY.csv") and not f.endswith("PY.csv"))]
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

def plotCellTimeSeriesData(cellID):
    X = 'T'
    Y = np.array(['B', 'F', 'dx', 'dy', 'dz', 'D2', 'D3'])
    units = np.array([' (mT)', ' (pN)', ' (µm)', ' (µm)', ' (µm)', ' (µm)', ' (µm)'])
    timeSeriesDataFrame = getCellTimeSeriesData(cellID)
    print(timeSeriesDataFrame.shape)
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
        plt.gcf().show()
    else:
        print('cell not found')
        
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

def analyseTimeSeries_ctField(tsDf):
    results = {}
    results['duration'] = np.max(tsDf['T'])
    results['medianRawB'] = np.median(tsDf.B)
    results['medianThickness'] = np.median(tsDf.D3)
    results['1stDThickness'] = np.percentile(tsDf.D3, 10)
    results['9thDThickness'] = np.percentile(tsDf.D3, 90)
    results['fluctuAmpli'] = results['9thDThickness'] - results['1stDThickness']
    results['validated'] = (results['1stDThickness'] > 0)
    X, Y = tsDf['T'], tsDf['D3']
    p, residuals, rank, singular_values, rcond = np.polyfit(X, Y, deg=5, full=True)
    Y2 = np.zeros(len(X))
    for i in range(len(X)):
        deg = len(p)-1
        for k in range(deg+1):
            Y2[i] += p[k]*(X[i]**(deg-k))
    results['R2_polyFit'] = get_R2(Y, Y2)
    return(results)



def createDataDict_ctField(list_ctFieldFiles):
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
        current_tsDf = pd.read_csv(tS_DataFilePath, ';')
        current_resultDict = analyseTimeSeries_ctField(current_tsDf)
        for k in current_resultDict.keys():
            tableDict[k].append(current_resultDict[k])
    return(tableDict)



def computeGlobalTable_ctField(task = 'fromScratch', fileName = 'Global_CtFieldData', save = False):
    ctFieldTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) \
                                  if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                                      and ('thickness' in f))]
#     print(ctFieldTimeSeriesDataFiles)
    if task == 'fromScratch':
        # create a dict containing the data
        tableDict = createDataDict_ctField(ctFieldTimeSeriesDataFiles) # MAIN SUBFUNCTION
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
        new_tableDict = createDataDict_ctField(new_ctFieldTimeSeriesDataFiles) # MAIN SUBFUNCTION
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

fit_intervals = [S for S in range(0,450,50)] + [S for S in range(500,1100,100)] + [S for S in range(1500,2500,500)]
regionFitsNames = []
for ii in range(len(fit_intervals)-1):
    for jj in range(ii+1, len(fit_intervals)):
        regionFitsNames.append(str(fit_intervals[ii]) + '<s<' + str(fit_intervals[jj]))
        
fit_toPlot = ['50<s<200', '200<s<350', '350<s<500', '500<s<700', '700<s<1000', '1000<s<1500', '1500<s<2000']

#### >>> OPTION TO LIGHTEN THE COMPUTATION
regionFitsNames = fit_toPlot

mask_fitToPlot = np.array(list(map(lambda x : x in fit_toPlot, regionFitsNames)))

for rFN in regionFitsNames:
    listColumnsMeca += ['KChadwick_'+rFN, 'R2Chadwick_'+rFN, 'K2Chadwick_'+rFN, 'H0Chadwick_'+rFN, 'Npts_'+rFN, 'validatedFit_'+rFN] 
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
#     initH0, initE = initH0*(initH0>0), initE*(initE>0)
    
    initialParameters = [initE, initH0]
#     print(initialParameters)

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
        R2 = get_R2(hCompr, hPredict)
        
        Chi2 = get_Chi2(strainCompr, strainPredict, dof, err)

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
        
        R2 = get_R2(strainCompr, strainPredict)
        
        Chi2 = get_Chi2(strainCompr, strainPredict, dof, err)

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



def analyseTimeSeries_meca(f, tsDF, expDf, listColumnsMeca, PLOT, PLOT_SHOW):
    
    plotSmallElements = True
    
    #### (0) Import experimental infos
    split_f = f.split('_')
    tsDF.dx, tsDF.dy, tsDF.dz, tsDF.D2, tsDF.D3 = tsDF.dx*1000, tsDF.dy*1000, tsDF.dz*1000, tsDF.D2*1000, tsDF.D3*1000
    # thisManipID = split_f[0] + '_' + split_f[1]
    # expDf['manipID'] = expDf['date'] + '_' + expDf['manip']
    # thisExpDf = expDf.loc[expDf['manipID'] == thisManipID]
    thisManipID = findInfosInFileName(f, 'manipID')
    thisExpDf = expDf.loc[expDf['manipID'] == thisManipID]

    # Deal with the asymmetric pair case : the diameter can be for instance 4503 (float) or '4503_2691' (string)
    diameters = thisExpDf.at[thisExpDf.index.values[0], 'bead diameter']
    if len(diameters) == 2:
        DIAMETER = (int(diameters[0]) + int(diameters[1]))/2.
    else:
        DIAMETER = int(diameters[0])
    
    EXPTYPE = str(thisExpDf.at[thisExpDf.index.values[0], 'experimentType'])

    results = {}
    for c in listColumnsMeca:
        results[c] = []

    Ncomp = max(tsDF['idxAnalysis'])
    loopStruct = thisExpDf.at[thisExpDf.index.values[0], 'loop structure']
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
        compField = thisExpDf.at[thisExpDf.index.values[0], 'ramp field']
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
        # 5th plot - fig5 & ax5, gather all the stress-strain curves with local fits, and will be completed later in the code
        fig5, ax5 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
        
        if plotSmallElements:
            # 6th plot - fig6 & ax6, gather all the "small elements", and will be completed later in the code
            fig6, ax6 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))


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
            iStart = findFirst(np.abs(tsDF['idxAnalysis']), i)
            for c in colToCorrect:
                jump = tsDF[c].values[iStart+2] - tsDF[c].values[iStart-1]
                tsDF.loc[maskCompAndPrecomp, c] -= jump
                if c == 'D3':
                    D3corrected = True
                    jumpD3 = jump
        else:
            D3corrected = False
            jumpD3 = 0
            
        #### (2) Segment the compression n°i
        thisCompDf = tsDF.loc[tsDF['idxAnalysis'] == i,:]
        iStart = (findFirst(tsDF['idxAnalysis'], i))
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
                
            best_H0 = max (H0, bestH0)
            
            if not findH0_fitError:
                results['bestH0'].append(bestH0)
                
                deltaCompr = (bestH0 - hCompr)/1000
                stressCompr = fCompr / (np.pi * (DIAMETER/2000) * deltaCompr)
                strainCompr = deltaCompr / (3*(bestH0/1000))
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
                             'fitError' : [], 'validatedFit' : [], 'Npts' : []} 
            # 'E' : [], 'H0' : []
            
            #### SETTING ! Setting of the region fits
            
            
            if not findH0_fitError:
                fitConditions = []
                for ii in range(len(fit_intervals)-1):
                    for jj in range(ii+1, len(fit_intervals)):
                        fitConditions.append((stressCompr > fit_intervals[ii]) & (stressCompr < fit_intervals[jj]))
                
                #### >>> OPTION TO LIGHTEN THE COMPUTATION
                fitConditions = []
                for kk in range(len(fit_toPlot)):
                    ftP = fit_toPlot[kk]
                    lowS, highS = int(ftP.split('<s<')[0]), int(ftP.split('<s<')[1])
                    fitConditions.append((stressCompr > lowS) & (stressCompr < highS))
                
                N_Fits = len(fitConditions)
                # print(N_Fits)
                # print(regionFitsNames)
                
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
                            compressionFitChadwick_StressStrain(hCompr_region, fCompr_region, bestH0, DIAMETER)
                            
                        
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
                        K_region, R2_region  = np.nan, np.nan
                        K2_region, H0_region = np.nan, np.nan
                        # E_region  = np.nan
                        
                    # Fill dictRegionFit (reinitialized for each cell)
                    dictRegionFit['regionFitNames'].append(regionFitName)
                    dictRegionFit['Npts'].append(Npts_region)
                    dictRegionFit['K'].append(K_region)
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
                    dictRegionFit['R2'].append(np.nan)
                    dictRegionFit['K2'].append(np.nan)
                    dictRegionFit['H0'].append(np.nan)
                    # dictRegionFit['E'].append(E_region)
                    dictRegionFit['fitError'].append(True)
                    dictRegionFit['validatedFit'].append(False)
                                        
                    results['Npts_'+rFN].append(np.nan)
                    results['KChadwick_'+rFN].append(np.nan)
                    results['R2Chadwick_'+rFN].append(np.nan)
                    results['K2Chadwick_'+rFN].append(np.nan)
                    results['H0Chadwick_'+rFN].append(np.nan)
                    # results['EChadwick_'+rFN].append(np.nan)
                    results['validatedFit_'+rFN].append(False)
            
            
            #### PLOT [2/4]
            # Complete fig 1, 2, 3 with the results of the fit
            if PLOT:
                
                #### fig1
                if not fitError:
                    if validatedFit:
                        ax1.plot(thisCompDf['T'].values, thisCompDf['D3'].values-DIAMETER, color = 'chartreuse', linestyle = '-', linewidth = 1.25)
                    else:
                        ax1.plot(thisCompDf['T'].values, thisCompDf['D3'].values-DIAMETER, color = 'gold', linestyle = '-', linewidth = 1.25)
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
                    if not validatedFit:
                        titleText += '\nNON VALIDATED'
                        
                    if not findH0_fitError:
                        thisAx3.plot(stressCompr, strainCompr, 'b+')
                        thisAx3.plot([np.percentile(stressCompr,10), np.percentile(stressCompr,90)], 
                                     [np.percentile(stressCompr,10)/E, np.percentile(stressCompr,90)/E],
                                     'r--', linewidth = 1.2, label = legendText)
                        thisAx3.legend(loc = 'lower right', prop={'size': 6})
                    thisAx3.set_xlabel('sigma (Pa)')
                    thisAx3.set_ylabel('epsilon')
                else:
                    titleText += '\nFIT ERROR'
                    
                if not findH0_fitError:
                    # Computations to display the fit of the bestH0
                    min_f = np.min(fCompr)
                    low_f = np.linspace(0, min_f, 20)
                    R = DIAMETER/2
                    low_h = bestH0 - ((3*bestH0*low_f)/(np.pi*(findH0_E/1e6)*R))**0.5
                    # def inversedChadwickModel(f, E, H0):
                    #     R = DIAMETER/2
                    #     h = H0 - ((3*H0*f)/(np.pi*E*R))**0.5
                    #     return(h)
                    
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

                Npts_fitToPlot = np.array(dictRegionFit['Npts'])[mask_fitToPlot]
                K_fitToPlot = np.array(dictRegionFit['K'])[mask_fitToPlot]
                K2_fitToPlot = np.array(dictRegionFit['K2'])[mask_fitToPlot]
                H0_fitToPlot = np.array(dictRegionFit['H0'])[mask_fitToPlot]
                R2_fitToPlot = np.array(dictRegionFit['R2'])[mask_fitToPlot]
                fitError_fitToPlot = np.array(dictRegionFit['fitError'])[mask_fitToPlot]
                validatedFit_fitToPlot = np.array(dictRegionFit['validatedFit'])[mask_fitToPlot]
                if not findH0_fitError:
                    fitConditions_fitToPlot = np.array(fitConditions)[mask_fitToPlot]
                
                
                # ax2[i-1] with the 1 line plot
                if nRowsSubplot == 1:
                    thisAx4 = ax4[colSp]
                    thisAx5 = ax5[colSp]
                elif nRowsSubplot >= 1:
                    thisAx4 = ax4[rowSp,colSp]
                    thisAx5 = ax5[rowSp,colSp]
                
                main_color = 'k' # my_default_color_list[0]
                
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
                        
                        color = my_default_color_list[k+1]
                        legendText4 = ''
                    
                    
                        # hPredict_fit0 = bestH0 - ((3*bestH0*fCompr_fit)/(np.pi*(K_fit/1e6)*R))**0.5
                        # strainPredict_fit = list_strainPredict_fitToPlot[k]
                        # hPredict_fit = bestH0 * (1 - 3*strainPredict_fit)
                        # legendText4 += 'Range ' + fit + '\n'
                        # legendText4 += 'K = {:.2e}Pa'.format(K_fit)
                        # thisAx4.plot(hPredict_fit, fCompr_fit, color = color, ls = '--', linewidth = 1.8, label = legendText4)
                        
                        # legendText4 = ''
                        legendText4 += 'Range ' + fit + '\n'
                        legendText4 += 'K2 = {:.2e}Pa'.format(K2_fit)
                        thisAx4.plot(hPredict_fit2, fCompr_fit, color = color, ls = '-', linewidth = 1.8, label = legendText4)
                        thisAx4.legend(loc = 'upper right', prop={'size': 5})
                    else:
                        pass
                        # thisAx4.plot([], [], color = color, ls = '--', linewidth = 0.8, label = legendText4)
                
                
                thisAx5.set_xlabel('epsilon')
                thisAx5.set_ylabel('sigma (Pa)')
                if not fitError and not findH0_fitError:
                    thisAx5.plot(strainCompr, stressCompr, color = main_color, marker = 'o', markersize = 3, ls = '', alpha = 0.8)
                    
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
                        
                        color = my_default_color_list[k+1]
                        legendText5 = ''
                        
                        if not fitError_fit:
                            legendText5 += 'Range ' + fit + '\n'
                            legendText5 += 'K = {:.2e}Pa'.format(K_fit)
                            thisAx5.plot(strainPredict_fit, stressCompr_fit,  
                                         color = color, ls = '-', linewidth = 1.8, label = legendText5)
                            thisAx5.legend(loc = 'upper left', prop={'size': 5})
                        else:
                            pass
                            # thisAx5.plot([], [], color = color, ls = '--', linewidth = 0.8, label = legendText5)

                multiAxes = [thisAx4, thisAx5]
                
                for ax in multiAxes:
                    ax.title.set_text(titleText)
                    for item in ([ax.title, ax.xaxis.label, \
                                  ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
                        item.set_fontsize(9)
                        
                #### fig6
                if plotSmallElements:
                    
                    def def2delta(e):
                        d = 3*bestH0*e
                        return(d)
                    
                    def delta2def(d):
                        e = d/(3*bestH0)
                        return(e)

                    if nRowsSubplot == 1:
                        thisAx6 = ax6[colSp]
                    elif nRowsSubplot >= 1:
                        thisAx6 = ax6[rowSp,colSp]
                        
                    thisAx6.set_xlabel('epsilon')
                    # secax = thisAx6.secondary_xaxis('top', functions=(def2delta, delta2def))
                    # secax.set_xlabel('delta (nm)')
                    thisAx6.set_ylabel('a/h0')
                    legendText6 = ''
                    
                    if not fitError and not findH0_fitError:
                        
                        A = (deltaCompr*(DIAMETER/2000))**0.5
                        Xi = A/(bestH0/1000)
                        # legendText6 = 'bestH0 = {:.2f}nm'.format(bestH0)
                        
                        thisAx6.plot(strainCompr, Xi*(2**0.5), color = 'k', ls = '', marker = 'o', label = 'a/h0_dim', markersize = 2)
                        thisAx6.plot(strainCompr, Xi*2, color = 'k', ls = '', marker = 'd', label = 'a/h0_chad', markersize = 2)
                        thisAx6.legend(loc = 'upper left', prop={'size': 5})
                        thisAx6.set_ylim([0,5])
                        minPlot, maxPlot = thisAx6.get_xlim()
                        thisAx6.set_xlim([0,maxPlot])
                        
                        # thisAx6bis.tick_params(axis='y', labelcolor='b')
                        thisAx6bis = thisAx6.twinx()
                        color = 'firebrick'
                        thisAx6bis.set_ylabel('F (pN)', color=color)
                        thisAx6bis.plot(strainCompr, fCompr, color=color)
                        thisAx6bis.set_ylim([0,1400])
                        # thisAx6bis.tick_params(axis='y', labelcolor=color)
                        # thisAx6bis.set_yticks([0,500,1000,1500])
                        # minh = np.min(tsDF['D3'].values-DIAMETER)
                        
                        
                        lowXi = (2*Xi < 1.1)
                        if len(strainCompr[lowXi]) > 0:
                            strainLimit = np.max(strainCompr[lowXi])
                            minPlot, maxPlot = thisAx5.get_ylim()
                            thisAx5.plot([strainLimit, strainLimit], [minPlot, maxPlot], color = 'gold', ls = '--')
                            minPlot, maxPlot = thisAx6.get_ylim()
                            thisAx6.plot([strainLimit, strainLimit], [minPlot, maxPlot], color = 'gold', ls = '--')
                        
                                
                        multiAxes = [thisAx6, thisAx6bis]
                        
                        for ax in multiAxes:
                            ax.title.set_text(titleText + '\nbestH0 = {:.2f}nm'.format(bestH0))
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
            #### jumpD3
            results['jumpD3'].append(np.nan) 
            results['validatedThickness'].append(validatedThickness)
            validatedFit = False
            results['critFit'].append('Not relevant')
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
                results['H0Chadwick_'+rFN].append(np.nan)
                results['EChadwick_'+rFN].append(np.nan)
                results['R2Chadwick_'+rFN].append(np.nan)
                results['Npts_'+rFN].append(np.nan)
                results['validatedFit_'+rFN].append(False)

        if not doThisCompAnalysis:
            print('Curve not suitable for analysis !')
            print(results['cellID'][-1])
            print('Compression no ' + str(i))

    
    #### PLOT [3/4]
    
    if PLOT:
        #### (after the jump correction fix, the plot of the main curve on fig1 came here)
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
            title = thisAx3.title.get_text()
            if not 'NON VALIDATED' in title:
                if thisAx3.get_ylim()[0] < eMin:
                    eMin = thisAx3.get_ylim()[0]
                if thisAx3.get_ylim()[1] > eMax:
                    eMax = thisAx3.get_ylim()[1]
                if thisAx3.get_xlim()[0] < sMin:
                    sMin = thisAx3.get_xlim()[0]
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
            title = thisAx3.title.get_text()
            if not 'NON VALIDATED' in title:
                thisAx3.set_xlim([sMin, sMax])
                thisAx3.set_ylim([eMin, eMax])
    
    #### PLOT [4/4]
    # Save the figures
    if PLOT:
        archiveFig(fig1, ax1, name=results['cellID'][-1] + '_h(t)', figSubDir = 'MecaAnalysis_allCells')
        archiveFig(fig2, ax2, name=results['cellID'][-1] + '_F(h)', figSubDir = 'MecaAnalysis_allCells')
        archiveFig(fig3, ax3, name=results['cellID'][-1] + '_sig(eps)', figSubDir = 'MecaAnalysis_allCells')
        archiveFig(fig4, ax4, name=results['cellID'][-1] + '_F(h)_regionFits', figSubDir = 'MecaAnalysis_allCells')
        archiveFig(fig5, ax5, name=results['cellID'][-1] + '_sig(eps)_regionFits', figSubDir = 'MecaAnalysis_allCells')
        archiveFig(fig6, ax6, name=results['cellID'][-1] + '_smallElements', figSubDir = 'MecaAnalysis_allCells')
        
        
        if PLOT_SHOW:
            fig1.show()
            fig2.tight_layout()
            fig2.show()
        else:
            plt.close('all')

    return(results)



def createDataDict_meca(list_mecaFiles, listColumnsMeca, PLOT):
    """
    Subfunction of computeGlobalTable_meca
    Create the dictionnary that will be converted in a pandas table in the end.
    """
    expDf = getExperimentalConditions(experimentalDataDir)
    tableDict = {}
    Nfiles = len(list_mecaFiles)
    PLOT_SHOW = (Nfiles==0)
    if not PLOT_SHOW:
        plt.ioff()
    for c in listColumnsMeca:
        tableDict[c] = []
    for f in list_mecaFiles: #[:10]:
        tS_DataFilePath = os.path.join(timeSeriesDataDir, f)
        current_tsDF = pd.read_csv(tS_DataFilePath, sep = ';')
         # MAIN SUBFUNCTION
        current_resultDict = analyseTimeSeries_meca(f, current_tsDF, expDf, 
                                                    listColumnsMeca, PLOT, PLOT_SHOW)
        for k in current_resultDict.keys():
            tableDict[k] += current_resultDict[k]
#     for k in tableDict.keys():
#         print(k, len(tableDict[k]))
    return(tableDict)



def computeGlobalTable_meca(task = 'fromScratch', fileName = 'Global_MecaData', save = False, PLOT = False, \
                            source = 'Matlab', listColumnsMeca=listColumnsMeca):
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
                      and ('R40' in f) and not (suffixPython in f))]
        
    elif source == 'Python':
        list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir) \
                      if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv") \
                      and (('R40' in f) or ('L40' in f)) and (suffixPython in f))]
        # print(list_mecaFiles)
    
#     print(list_mecaFiles)
    
    if task == 'fromScratch':
        # create a dict containing the data
        tableDict = createDataDict_meca(list_mecaFiles, listColumnsMeca, PLOT) # MAIN SUBFUNCTION
        # create the dataframe from it
        meca_DF = pd.DataFrame(tableDict)
        
        # last step: now that the dataFrame is complete, one can use "compStartTimeThisDay" col to compute the start time of each compression relative to the first one done this day.
        allDates = list(meca_DF['date'].unique())
        for d in allDates:
            subDf = meca_DF.loc[meca_DF['date'] == d]
            experimentStartTime = np.min(subDf['compStartTimeThisDay'])
            meca_DF['compStartTimeThisDay'].loc[meca_DF['date'] == d] = meca_DF['compStartTimeThisDay'] - experimentStartTime
        
    elif task == 'updateExisting':
        # get existing table
        try:
            savePath = os.path.join(dataDir, (fileName + '.csv'))
            existing_meca_DF = pd.read_csv(savePath, sep=';')
        except:
            print('No existing table found')
            
        # find which of the time series files are new
        new_list_mecaFiles = []
        for f in list_mecaFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            if currentCellID not in existing_meca_DF.cellID.values:
                new_list_mecaFiles.append(f)
                
        # create the dict with new data
        new_tableDict = createDataDict_meca(new_list_mecaFiles, listColumnsMeca, PLOT) # MAIN SUBFUNCTION
        # create the dataframe from it
        new_meca_DF = pd.DataFrame(new_tableDict)
        # fuse the existing table with the new one
        meca_DF = pd.concat([existing_meca_DF, new_meca_DF])
        
    else: # If task is neither 'fromScratch' nor 'updateExisting'
    # Then task can be a substring that can be in some timeSeries file !
    # It will create a table with only these files, WITHOUT SAVING IT !
    # But it can plot figs from it.
        save = False
        new_list_mecaFiles = []
        for f in list_mecaFiles:
            split_f = f.split('_')
            currentCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
            if task in currentCellID:
                new_list_mecaFiles.append(f)
        # create the dict with new data
        new_tableDict = createDataDict_meca(new_list_mecaFiles, listColumnsMeca, PLOT) # MAIN SUBFUNCTION
        # create the dataframe from it
        meca_DF = pd.DataFrame(new_tableDict)
    
    for c in meca_DF.columns:
            if 'Unnamed' in c:
                meca_DF = meca_DF.drop([c], axis=1)
    
    if save:
        saveName = fileName + '.csv'
        savePath = os.path.join(dataDir, saveName)
        meca_DF.to_csv(savePath, sep=';')
    
    delta = time.time() - top
    print(delta)
    
    return(meca_DF)
            

    
def getGlobalTable_meca(fileName):
    try:
        savePath = os.path.join(dataDir, (fileName + '.csv'))
        meca_DF = pd.read_csv(savePath, sep=';')
        print('Extracted a table with ' + str(meca_DF.shape[0]) + ' lines and ' + str(meca_DF.shape[1]) + ' columns.')
    except:
        print('No existing table found')
    for c in meca_DF.columns:
            if 'Unnamed' in c:
                meca_DF = meca_DF.drop([c], axis=1)
    
    if 'ExpDay'in meca_DF.columns:
        dateExemple = meca_DF.loc[meca_DF.index[0],'ExpDay']
        if not ('manipID' in meca_DF.columns):
            meca_DF['manipID'] = meca_DF['ExpDay'] + '_' + meca_DF['CellID'].apply(lambda x: x.split('_')[0])
            
    elif 'date'in meca_DF.columns:
        dateExemple = meca_DF.loc[meca_DF.index[0],'date']
        if not ('manipID' in meca_DF.columns):
            meca_DF['manipID'] = meca_DF['date'] + '_' + meca_DF['cellName'].apply(lambda x: x.split('_')[0])
    
    if re.match(dateFormatExcel, dateExemple):
        print('bad date')
    return(meca_DF)

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
        GlobalTable_ctField = getGlobalTable_ctField()
        experimentalDataDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis//Data_Experimental"
        expDf = getExperimentalConditions(experimentalDataDir)
        fluoDf = getFluoData()
        GlobalTable_ctField = pd.merge(expDf, GlobalTable_ctField, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable_ctField = pd.merge(GlobalTable_ctField, fluoDf, how="left", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable_ctField = removeColumnsDuplicate(GlobalTable_ctField)
        print('Merged table has ' + str(GlobalTable_ctField.shape[0]) + ' lines and ' + str(GlobalTable_ctField.shape[1]) + ' columns.')
        
        # print(GlobalTable_ctField.head())
        
        return(GlobalTable_ctField)

    elif kind == 'meca_matlab':
        GlobalTable_meca_Matlab = getGlobalTable_meca('Global_MecaData')
        expDf = getExperimentalConditions(experimentalDataDir)
        fluoDf = getFluoData()
        GlobalTable_meca_Matlab = pd.merge(GlobalTable_meca_Matlab, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable_meca_Matlab = pd.merge(GlobalTable_meca_Matlab, fluoDf, how="left", left_on='CellName', right_on='cellID'
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable_meca_Matlab.shape[0]) + ' lines and ' + str(GlobalTable_meca_Matlab.shape[1]) + ' columns.')
        
        # print(GlobalTable_meca_Matlab.tail())
        GlobalTable_meca_Matlab = removeColumnsDuplicate(GlobalTable_meca_Matlab)
        return(GlobalTable_meca_Matlab)


    elif kind == 'meca_py':
        GlobalTable_meca_Py = getGlobalTable_meca('Global_MecaData_Py')
        expDf = getExperimentalConditions(experimentalDataDir)
        fluoDf = getFluoData()
        GlobalTable_meca_Py = pd.merge(GlobalTable_meca_Py, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable_meca_Py = pd.merge(GlobalTable_meca_Py, fluoDf, how="left", left_on='cellID', right_on='cellID'
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable_meca_Py.shape[0]) + ' lines and ' + str(GlobalTable_meca_Py.shape[1]) + ' columns.')
        
        # print(GlobalTable_meca_Py.tail())
        GlobalTable_meca_Py = removeColumnsDuplicate(GlobalTable_meca_Py)
        return(GlobalTable_meca_Py)


    elif kind == 'meca_py2':
        GlobalTable_meca_Py2 = getGlobalTable_meca('Global_MecaData_Py2')
        expDf = getExperimentalConditions(experimentalDataDir)
        fluoDf = getFluoData()
        GlobalTable_meca_Py2 = pd.merge(GlobalTable_meca_Py2, expDf, how="inner", on='manipID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        GlobalTable_meca_Py2 = pd.merge(GlobalTable_meca_Py2, fluoDf, how="left", left_on='cellID', right_on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        print('Merged table has ' + str(GlobalTable_meca_Py2.shape[0]) + ' lines and ' + str(GlobalTable_meca_Py2.shape[1]) + ' columns.')

        # print(GlobalTable_meca_Py2.tail())
        GlobalTable_meca_Py2 = removeColumnsDuplicate(GlobalTable_meca_Py2)
        return(GlobalTable_meca_Py2)
