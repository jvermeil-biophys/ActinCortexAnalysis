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
import re
from datetime import date
import sys
import scipy.stats as st


# Local imports
COMPUTERNAME = os.environ['COMPUTERNAME']
if COMPUTERNAME == 'ORDI-JOSEPH':
    mainDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis"
    rawDir = "D://MagneticPincherData"
    ownCloudDir = "C://Users//JosephVermeil//ownCloud//ActinCortexAnalysis"
elif COMPUTERNAME == 'LARISA':
    mainDir = "C://Users//Joseph//Desktop//ActinCortexAnalysis"
    rawDir = "F://JosephVermeil//MagneticPincherData"    
    ownCloudDir = "C://Users//Joseph//ownCloud//ActinCortexAnalysis"
elif COMPUTERNAME == 'DESKTOP-K9KOJR2':
    mainDir = "C://Users//anumi//OneDrive//Desktop//ActinCortexAnalysis"
    rawDir = "D:/Anumita/MagneticPincherData"  
elif COMPUTERNAME == '':
    mainDir = "C://Users//josep//Desktop//ActinCortexAnalysis"
    ownCloudDir = "C://Users//josep//ownCloud//ActinCortexAnalysis"

# Add the folder to path
sys.path.append(mainDir + "//Code_Python")
import utilityFunctions_JV as jvu

#%% Global constants

bead_dia = 4.503


# These regex are used to correct the stupid date conversions done by Excel
dateFormatExcel = re.compile(r'\d{2}/\d{2}/\d{4}')
dateFormatExcel2 = re.compile(r'\d{2}-\d{2}-\d{4}')
dateFormatOk = re.compile(r'\d{2}-\d{2}-\d{2}')

SCALE_100X = 15.8 # pix/µm 
NORMAL  = '\033[0m'
RED  = '\033[31m' # red
GREEN = '\033[32m' # green
ORANGE  = '\033[33m' # orange
BLUE  = '\033[36m' # blue

# %% Directories adress

experimentalDataDir = os.path.join(mainDir, "Data_Experimental_AJ")
dataDir = os.path.join(mainDir, "Data_Analysis")
timeSeriesDataDir = os.path.join(dataDir, "TimeSeriesData")

figDir = os.path.join(dataDir, "Figures")
todayFigDir = os.path.join(figDir, "Historique/" + str(date.today()))

#%% Utility functions specific to Opto experiments

def findActivation(fieldDf):
    maxZidx = fieldDf['Z'].argmax() #Finding the index of the max Z
    maxZ = fieldDf['Z'][maxZidx] #To check if the value is correct
    return maxZidx, maxZ

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


def getOptoMeta(cellID):
    date = jvu.findInfosInFileName(cellID, 'date')
    date = date.replace('-', '.')
    optoMetaDataPath = rawDir+'//Raw//'+date
    allOptoMetaDataFiles = [f for f in os.listdir(optoMetaDataPath) 
                          if (os.path.isfile(os.path.join(optoMetaDataPath, f)) 
                              and f.endswith("OptoMetadata.txt"))]
    fileFound = False
    nFile = len(allOptoMetaDataFiles)
    iFile = 0
    while (not fileFound) and (iFile < nFile):
        f = allOptoMetaDataFiles[iFile]
        if f.startswith(cellID + '_'):
            optoMetaDataPath = os.path.join(optoMetaDataPath, f)
            optoMetaDatadf = pd.read_csv(optoMetaDataPath, sep='\t')
            fileFound = True
        iFile += 1
    if not fileFound:
        optoMetaDatadf = pd.DataFrame([])
    else:
        for c in optoMetaDatadf.columns:
                if 'Unnamed' in c:
                    optoMetaDatadf = optoMetaDatadf.drop([c], axis=1)
    return(optoMetaDatadf)

#%% Statistical functions

def addStat_df(ax, data, box_pairs, param, cond, test = 'Mann-Whitney', percentHeight = 95):
    refHeight = np.percentile(data[param].values, percentHeight)
    currentHeight = refHeight
    scale = ax.get_yscale()
    xTicks = ax.get_xticklabels()
    dictXTicks = {xTicks[i].get_text() : xTicks[i].get_position()[0] for i in range(len(xTicks))}
    for bp in box_pairs:
        c1 = data[data[cond] == bp[0]][param].values
        c2 = data[data[cond] == bp[1]][param].values
        if test == 'Mann-Whitney' or test == 'Wilcox_2s' or test == 'Wilcox_greater' or test == 'Wilcox_less' or test == 't-test':
            if test=='Mann-Whitney':
                statistic, pval = st.mannwhitneyu(c1,c2)
            elif test=='Wilcox_2s':
                statistic, pval = st.wilcoxon(c1,c2, alternative = 'two-sided')
            elif test=='Wilcox_greater':
                statistic, pval = st.wilcoxon(c1,c2, alternative = 'greater')
            elif test=='Wilcox_less':
                statistic, pval = st.wilcoxon(c1,c2, alternative = 'less')
            elif test=='t-test':
                statistic, pval = st.ttest_ind(c1,c2)
            text = 'ns'
            if pval < 0.05 and pval > 0.01:
                text = '*'
            elif pval < 0.01 and pval > 0.001:
                text = '**'
            elif pval < 0.001 and pval < 0.001:
                text = '***'
            elif pval < 0.0001:
                text = '****'
            ax.plot([bp[0], bp[1]], [currentHeight, currentHeight], 'k-', lw = 1)
            XposText = (dictXTicks[bp[0]]+dictXTicks[bp[1]])/2
            if scale == 'log':
                power = 0.01* (text=='ns') + 0.000 * (text!='ns')
                YposText = currentHeight*(refHeight**power)
            else:
                factor = 0.03 * (text=='ns') + 0.000 * (text!='ns')
                YposText = currentHeight + factor*refHeight
            ax.text(XposText, YposText, text, ha = 'center', color = 'k')
    #         if text=='ns':
    #             ax.text(posText, currentHeight + 0.025*refHeight, text, ha = 'center')
    #         else:
    #             ax.text(posText, currentHeight, text, ha = 'center')
            if scale == 'log':
                currentHeight = currentHeight*(refHeight**0.05)
            else:
                currentHeight =  currentHeight + 0.15*refHeight
        ax.set_ylim([ax.get_ylim()[0], currentHeight])

        if test == 'pairwise':
            ratio = (c2/c1)
            stdError = np.nanstd(ratio)/np.sqrt(np.size(c1))
            confInt = np.nanmean(ratio) - 1.96 * stdError
            print(stdError)
            print(confInt)
            return confInt


def ctFieldThicknessIndividual(experimentalDataDir, todayFigDir, date, save = False, background = 'default'):
    try:
        os.mkdir(todayFigDir)
    except:
        pass
    
    expDf = jvu.getExperimentalConditions(experimentalDataDir, save = False, sep = ',') 
    files = os.listdir(rawDir+'/Raw/'+date)

    if background == 'dark':
        plt.style.use('dark_background')
    else:
        plt.style.use('default')

    for f in files:
        if f.endswith('.tif'):
            cellID = jvu.findInfosInFileName(f, 'cellID')
            timeSeriesDf = getCellTimeSeriesData(cellID)            
            optoMetaDataDf = getOptoMeta(cellID)
            manipID = jvu.findInfosInFileName(cellID, 'manipID')
            print(f)
            try:
               Tact = optoMetaDataDf['T_abs'].values[0]
            except:
                pass
            # Tact = timeSeriesDf['T'][timeSeriesDf['T_abs'] == Tact]
            fig = plt.figure(figsize=(20,20))
            plt.rcParams.update({'font.size': 25})
            plt.suptitle(cellID)
            t = timeSeriesDf['T'].values*1000/60
            
            plt.subplot(3, 1, 1)
            
            plt.plot(t, timeSeriesDf['D3'] - bead_dia)
            
            plt.subplot(3, 1, 2)
            plt.plot(t, timeSeriesDf['D2'] - bead_dia)
            
            plt.subplot(3, 1, 3)
            plt.plot(t, timeSeriesDf['dz'])
            
            plt.savefig(todayFigDir+'/'+cellID+'_ThicknessvTime')
            plt.show()

def ctFieldThicknessAll(experimentalDataDir, todayFigDir, date, tag = 'all', save = False, background = 'default'):
    bead_dia = 4.503
    
    try:
        os.mkdir(todayFigDir)
    except:
        pass
    
    if background == 'dark':
        plt.style.use('dark_background')
    else:
        plt.style.use('default')
    
    expDf = jvu.getExperimentalConditions(experimentalDataDir, save = False, sep = ';') 
    expDf = expDf[expDf['experimentType'] == 'optoGen']
    expDf = expDf[expDf['microscope'] == 'metamorph']
    cellConditionsDf = pd.read_csv(experimentalDataDir+'/cellConditions.csv')
    cellConditionsDf = cellConditionsDf[cellConditionsDf['excluded'] == 'no']
    allCells = cellConditionsDf['cellID'].values
    cellIDs = []
    
    if tag == 'all':
        cellIDs = allCells
    elif tag != 'all':
        selectedManipIDs = expDf['manipID'][expDf['activation type'] == tag].values
        for each in selectedManipIDs:
            cellID = (cellConditionsDf['cellID'][cellConditionsDf['cellID'].str.contains(each)])
            cellID = cellID.values
            if len(cellID) >= 1:
                cellIDs.extend(cellID)  
    cellIDs = np.asarray(cellIDs)
        
    for cellID in cellIDs:
        timeSeriesDf = getCellTimeSeriesData(cellID)            
        optoMetaDataDf = getOptoMeta(cellID)
        Tact = optoMetaDataDf['T_abs'].values[0]
        time = timeSeriesDf['T'].values*1000/60
        plt.plot(time, timeSeriesDf['D3'].values - bead_dia, label = cellID)
        plt.axvline(x = 5.0, color = 'r')
    plt.title('Thickness (um) vs Time : '+tag)
    plt.show()
    plt.legend()

def ctFieldThicknessSummary(experimentalDataDir, todayFigDir, listOfCells, background = 'default'):
    try:
        os.mkdir(todayFigDir)
    except:
        pass
    
    if background == 'dark':
        plt.style.use('dark_background')
    else:
        plt.style.use('default')
    
    
    expDf = jvu.getExperimentalConditions(experimentalDataDir, save = False, sep = ';')
    cellConditionsDf = pd.read_csv(experimentalDataDir+'/cellConditions.csv')
    summaryDict = {}
    summaryDict['cellID'] = []
    summaryDict['medianThickness'] = []
    summaryDict['medianThickness5minRange'] = []
    summaryDict['activationTag'] = []
    summaryDict['activationType'] = []
    summaryDict['fluctuations'] = []
    summaryDict['fluctuations5minRange'] = []
    summaryDict['blebCondition'] = []
    summaryDict['ratioThickness'] = []
    summaryDict['ratioFluctuations'] = []
    summaryDict['ratioFluctThick'] = []
    summaryDict['ratioFluctuations5min'] = []
    summaryDict['ratioThickness5min'] = []
    
    #### Constructing CT Field Summary Table
    for i in range(len(listOfCells)):
        
        cellID = listOfCells[i]
        print(cellID)
        timeSeriesDf = getCellTimeSeriesData(cellID)
        optoMetaDataDf = getOptoMeta(cellID)
        manipID = jvu.findInfosInFileName(cellID, 'manipID')

        Tact = optoMetaDataDf['T_abs'].values[0]
    
        summaryDict['cellID'].append(cellID)
        thicknessBefore = (timeSeriesDf['D3']-bead_dia)[(timeSeriesDf['Tabs']*1000 < Tact)]
        medianThicknessBefore = thicknessBefore.median()
        fluctBefore = np.percentile(thicknessBefore, 90) - np.percentile(thicknessBefore, 10)
        summaryDict['medianThickness'].append(medianThicknessBefore)
        summaryDict['activationTag'].append('Before')
        summaryDict['activationType'].append(expDf['activation type'][(expDf['manipID'] == manipID)].values[0])
        summaryDict['fluctuations'].append(fluctBefore)
        summaryDict['blebCondition'].append(cellConditionsDf['blebCondition'][cellConditionsDf['cellID']==cellID].values[0])
        summaryDict['ratioThickness'].append(np.nan)
        summaryDict['ratioFluctuations'].append(np.nan)
        summaryDict['ratioFluctThick'].append(np.nan)
        summaryDict['ratioFluctuations5min'].append(np.nan)
        summaryDict['ratioThickness5min'].append(np.nan)
        summaryDict['medianThickness5minRange'].append(medianThicknessBefore)
        summaryDict['fluctuations5minRange'].append(fluctBefore)
        
        summaryDict['cellID'].append(cellID)
        thicknessAfter = (timeSeriesDf['D3']-bead_dia)[(timeSeriesDf['Tabs']*1000 > Tact)]
        medianThicknessAfter = thicknessAfter.median()
        fluctAfter = np.percentile(thicknessAfter, 90) - np.percentile(thicknessAfter, 10)
        summaryDict['medianThickness'].append(medianThicknessAfter)
        summaryDict['activationTag'].append('After')
        summaryDict['activationType'].append(expDf['activation type'][(expDf['manipID'] == manipID)].values[0])
        summaryDict['fluctuations'].append(fluctAfter)
        summaryDict['blebCondition'].append(cellConditionsDf['blebCondition'][cellConditionsDf['cellID']==cellID].values[0])
        summaryDict['ratioThickness'].append(medianThicknessAfter/medianThicknessBefore)
        summaryDict['ratioFluctuations'].append(fluctAfter/fluctBefore)
        summaryDict['ratioFluctThick'].append((fluctAfter/fluctBefore)/(medianThicknessAfter/medianThicknessBefore))
        
        thicknessLast5min = thicknessAfter[timeSeriesDf['T']*1000/60 > 10.0]
        medianThicknessLast5min = thicknessLast5min.median()
        summaryDict['medianThickness5minRange'].append(medianThicknessLast5min)
        fluctuationsLast5min = np.percentile(thicknessLast5min, 90) - np.percentile(thicknessLast5min, 10)
        summaryDict['fluctuations5minRange'].append(fluctuationsLast5min)
        summaryDict['ratioFluctuations5min'].append(fluctuationsLast5min/fluctBefore)
        summaryDict['ratioThickness5min'].append(medianThicknessLast5min/medianThicknessBefore)

        summaryDf = pd.DataFrame(summaryDict)

        
    #### Subplots of each activation type, thickness
    activationType = ['global', 'at beads', 'away from beads']
    fig1, axs = plt.subplots(nrows = 1, ncols = 3)
    color = ['r', 'g', 'b']
    for j in range(len(activationType)):
        dataSpecific = summaryDf[summaryDf['activationType'] == activationType[j]]
        
        y1 = dataSpecific['medianThickness'][dataSpecific['activationTag'] == 'Before']
        y2 = dataSpecific['medianThickness'][dataSpecific['activationTag'] == 'After']
        
        axs[j].plot((np.zeros(len(y1), dtype = int), np.ones(len(y2), dtype = int)), (y1, y2), '-o', color = color[j])

        axs[j].set_title(activationType[j])
        labels = ['Before, After']
        axs[j].set_xticks = ([1,2], labels)
        axs[j].set_ylim(0, 1.5)
    plt.suptitle('Median Thickness')
    fig1.tight_layout()
    plt.savefig(todayFigDir+'/ActivationSpecificSubplots_Thickness.png')
    plt.show()
    
    #### Subplots of each activation type, fluctuations
    activationType = ['global', 'at beads', 'away from beads']
    fig2, axs = plt.subplots(nrows = 1, ncols = 3)
    color = ['r', 'g', 'b']
    for j in range(len(activationType)):
        dataSpecific = summaryDf[summaryDf['activationType'] == activationType[j]]
        
        y1 = dataSpecific['fluctuations'][dataSpecific['activationTag'] == 'Before']
        y2 = dataSpecific['fluctuations'][dataSpecific['activationTag'] == 'After']
        
        axs[j].plot((np.zeros(len(y1), dtype = int), np.ones(len(y2), dtype = int)), (y1, y2), '-o', color = color[j])
        # axs[j].plot((y1, y2), '-o', color = color[j])

        axs[j].set_title(activationType[j])
        labels = ['Before, After']
        axs[j].set_ylim(0,1.5)
        axs[j].set_xticks = ([1,2], labels)
    
    plt.suptitle('Fluctuations')
    fig2.tight_layout()
    plt.savefig(todayFigDir+'/ActivationSpecificSubplots_Fluctuations.png')
    plt.show()
    
    #### medianThickness first/last 5mins boxplots, global
    ax1a_0 = plt.figure()
    activationType = 'global'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax1a_0 = sns.boxplot(x = 'activationTag', y='medianThickness5minRange', data=data_onlyGlobal, color = 'skyblue')
    ax1a_0 = sns.swarmplot(x = 'activationTag', y='medianThickness5minRange', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax1a_0, data_onlyGlobal, [('After', 'Before')], 'medianThickness5minRange', test = 'Wilcox_greater', cond = 'activationTag')
    ax1a_0.set_ylim(0, 1)
    for patch in ax1a_0.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllThickness_Before-AfterActivation5mins_'+activationType+'.png')
    plt.show()
    
    ####  medianThickness first/last 5mins boxplots, at beads
    ax1b_0 = plt.figure()
    activationType = 'at beads'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax1b_0 = sns.boxplot(x = 'activationTag', y='medianThickness5minRange', data=data_onlyGlobal, color = 'skyblue')
    ax1b_0 = sns.swarmplot(x = 'activationTag', y='medianThickness5minRange', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax1b_0, data_onlyGlobal, [('After', 'Before')], 'medianThickness5minRange', test = 'Wilcox_greater', cond = 'activationTag')
    ax1b_0.set_ylim(0, 1)
    for patch in ax1b_0.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllThickness_Before-AfterActivation5mins_'+activationType+'.png')
    plt.show()
    
    ####  medianThickness first/last 5mins boxplots, away from beads
    ax1c_0 = plt.figure()
    activationType = 'away from beads'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax1c_0 = sns.boxplot(x = 'activationTag', y='medianThickness5minRange', data=data_onlyGlobal, color = 'skyblue')
    ax1c_0 = sns.swarmplot(x = 'activationTag', y='medianThickness5minRange', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax1c_0, data_onlyGlobal, [('After', 'Before')], 'medianThickness5minRange', test = 'Wilcox_greater', cond = 'activationTag')
    ax1c_0.set_ylim(0, 1)
    for patch in ax1c_0.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllThickness_Before-AfterActivation5mins_'+activationType+'.png')
    plt.show()
    
    #### medianThickness boxplots, global
    ax1a = plt.figure()
    activationType = 'global'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax1a = sns.boxplot(x = 'activationTag', y='medianThickness', data=data_onlyGlobal, color = 'skyblue')
    ax1a = sns.swarmplot(x = 'activationTag', y='medianThickness', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax1a, data_onlyGlobal, [('After', 'Before')], 'medianThickness', test = 'Wilcox_greater', cond = 'activationTag')
    ax1a.set_ylim(0, 1)
    for patch in ax1a.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllThickness_Before-AfterActivation_'+activationType+'.png')
    plt.show()
    
    ####  medianThickness boxplots, at beads
    ax1b = plt.figure()
    activationType = 'at beads'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax1b = sns.boxplot(x = 'activationTag', y='medianThickness', data=data_onlyGlobal, color = 'skyblue')
    ax1b = sns.swarmplot(x = 'activationTag', y='medianThickness', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax1b, data_onlyGlobal, [('After', 'Before')], 'medianThickness', test = 'Wilcox_greater', cond = 'activationTag')
    ax1b.set_ylim(0, 1)
    for patch in ax1b.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllThickness_Before-AfterActivation_'+activationType+'.png')
    plt.show()
    
    ####  medianThickness boxplots, away from beads
    ax1c = plt.figure()
    activationType = 'away from beads'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax1c = sns.boxplot(x = 'activationTag', y='medianThickness', data=data_onlyGlobal, color = 'skyblue')
    ax1c = sns.swarmplot(x = 'activationTag', y='medianThickness', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax1c, data_onlyGlobal, [('After', 'Before')], 'medianThickness', test = 'Wilcox_greater', cond = 'activationTag')
    ax1c.set_ylim(0, 1)
    for patch in ax1c.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllThickness_Before-AfterActivation_'+activationType+'.png')
    plt.show()
    
    ####  fluctuations first/last 5mins boxplots, away from beads
    ax2c_0 = plt.figure()
    activationType = 'away from beads'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax2c_0 = sns.boxplot(x = 'activationTag', y='fluctuations5minRange', data=data_onlyGlobal, color = 'skyblue')
    ax2c_0 = sns.swarmplot(x = 'activationTag', y='fluctuations5minRange', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax2c_0, data_onlyGlobal, [('After', 'Before')], 'fluctuations5minRange', test = 'Wilcox_greater', cond = 'activationTag')
    ax2c_0.set_ylim(0, 1)
    for patch in ax2c_0.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllFluctuations_Before-AfterActivation5mins_'+activationType+'.png')
    plt.show()
    
    #### fluctuations first/last 5mins boxplots, global
    ax2b_0 = plt.figure()
    activationType = 'global'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax2b_0 = sns.boxplot(x = 'activationTag', y='fluctuations5minRange', data=data_onlyGlobal, color = 'skyblue')
    ax2b_0 = sns.swarmplot(x = 'activationTag', y='fluctuations5minRange', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax2b_0, data_onlyGlobal, [('After', 'Before')], 'fluctuations5minRange', test = 'Wilcox_greater', cond = 'activationTag')
    ax2b_0.set_ylim(0, 1)
    plt.suptitle(activationType)
    for patch in ax2b_0.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.savefig(todayFigDir+'/AllFluctuations_Before-AfterActivation5mins_'+activationType+'.png')
    plt.show()
    
    #### fluctuations first/last 5mins boxplots, at beads
    ax2c_0 = plt.figure()
    activationType = 'at beads'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax2c_0 = sns.boxplot(x = 'activationTag', y='fluctuations5minRange', data=data_onlyGlobal, color = 'skyblue')
    ax2c_0 = sns.swarmplot(x = 'activationTag', y='fluctuations5minRange', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax2c_0, data_onlyGlobal, [('After', 'Before')], 'fluctuations5minRange', test = 'Wilcox_greater', cond = 'activationTag')
    ax2c_0.set_ylim(0, 1)
    for patch in ax2c_0.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllFluctuations_Before-AfterActivation5mins_'+activationType+'.png')
    plt.show()
    
    ####  fluctuations boxplots, away from beads
    ax2c = plt.figure()
    activationType = 'away from beads'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax2c = sns.boxplot(x = 'activationTag', y='fluctuations', data=data_onlyGlobal, color = 'skyblue')
    ax2c = sns.swarmplot(x = 'activationTag', y='fluctuations', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax2c, data_onlyGlobal, [('After', 'Before')], 'fluctuations', test = 'Wilcox_greater', cond = 'activationTag')
    ax2c.set_ylim(0, 1)
    for patch in ax2c.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllFluctuations_Before-AfterActivation_'+activationType+'.png')
    plt.show()
    
    #### fluctuations boxplots, global
    ax2b = plt.figure()
    activationType = 'global'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax2b = sns.boxplot(x = 'activationTag', y='fluctuations', data=data_onlyGlobal, color = 'skyblue')
    ax2b = sns.swarmplot(x = 'activationTag', y='fluctuations', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax2b, data_onlyGlobal, [('After', 'Before')], 'fluctuations', test = 'Wilcox_greater', cond = 'activationTag')
    ax2b.set_ylim(0, 1)
    plt.suptitle(activationType)
    for patch in ax2b.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.savefig(todayFigDir+'/AllFluctuations_Before-AfterActivation_'+activationType+'.png')
    plt.show()
    
    #### fluctuations boxplots, at beads
    ax2c = plt.figure()
    activationType = 'at beads'
    data_onlyGlobal = summaryDf[summaryDf['activationType'] == activationType]
    ax2c = sns.boxplot(x = 'activationTag', y='fluctuations', data=data_onlyGlobal, color = 'skyblue')
    ax2c = sns.swarmplot(x = 'activationTag', y='fluctuations', data=data_onlyGlobal, hue = 'blebCondition')
    addStat_df(ax2c, data_onlyGlobal, [('After', 'Before')], 'fluctuations', test = 'Wilcox_greater', cond = 'activationTag')
    ax2c.set_ylim(0, 1)
    for patch in ax2c.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.suptitle(activationType)
    plt.savefig(todayFigDir+'/AllFluctuations_Before-AfterActivation_'+activationType+'.png')
    plt.show()
    
    
    ####All cells thickness before/after activation
    ax2 = plt.figure()
    ax2 = sns.boxplot(x = 'activationTag', y='medianThickness', data=summaryDf, color = 'skyblue')
    ax2 = sns.swarmplot(x = 'activationTag', y='medianThickness', data=summaryDf, hue = 'activationType')
    addStat_df(ax2, summaryDf, [('After', 'Before')], 'medianThickness', test = 'Wilcox_greater', cond = 'activationTag')
    ax2.set_ylim(0, 1)
    for patch in ax2.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.savefig(todayFigDir+'/AllThickness_Before-AfterActivation.png')
    plt.show()
    
    
    ####Plot fluctuations before and after activation
    ax3 = plt.figure()
    ax3 = sns.boxplot(x = 'activationTag', y='fluctuations', data=summaryDf, color = 'skyblue')
    ax3 = sns.swarmplot(x = 'activationTag', y='fluctuations', data=summaryDf, hue = 'activationType')
    addStat_df(ax3, summaryDf, [('After', 'Before')], 'fluctuations', test = 'Wilcox_greater', cond = 'activationTag')
    ax3.set_ylim(0, 1)
    for patch in ax3.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    plt.savefig(todayFigDir+'/AllFluctuations_Before-AfterActivation.png')
    plt.show()
    
    ####Plot fluctuations vs. median thickness
    ax4 = plt.figure()
    ax4 = sns.scatterplot(x = 'medianThickness', y='fluctuations', data=summaryDf, hue = 'activationTag', style = 'activationType', s = 100)
    plt.savefig(todayFigDir+'/Summary_FluctuationsvsThickness.png')
    plt.show()
    
    #### Ratio of thicknesses first/last 5 mins after/before activation
    
    ax5_0 = plt.figure()
    ax5_0 = sns.scatterplot(x = 'cellID',  y = 'ratioThickness5min', data=summaryDf, style = 'blebCondition', hue = 'activationType', s = 100)
    ax5_0.set_xlabel('cellID')
    ax5_0.set_ylabel('After/Before Ratio Thickness (5 mins before/after)')
    # confInt = addStat_df(ax5, summaryDf, [('Before', 'After')], 'medianThickness', test = 'pairwise', cond = 'activationTag')
    # plt.suptitle("Mean - 1.96*StdError: " + str(confInt))
    ax5_0.set_xticklabels(summaryDf['cellID'].values, rotation = 90)
    ax5_0.axhline(y = 1.0, color = 'r')
    ax5_0.set_ylim(0,4)
    plt.savefig(todayFigDir+'/Summary_RatioThickness5mins.png')
    plt.show()
    
    #### Ratio of thicknesses after/before activation
    
    ax5 = plt.figure()
    ax5 = sns.scatterplot(x = 'cellID',  y = 'ratioThickness', data=summaryDf, style = 'blebCondition', hue = 'activationType', s = 100)
    ax5.set_xlabel('cellID')
    ax5.set_ylabel('After/Before Ratio Thickness')
    # confInt = addStat_df(ax5, summaryDf, [('Before', 'After')], 'medianThickness', test = 'pairwise', cond = 'activationTag')
    # plt.suptitle("Mean - 1.96*StdError: " + str(confInt))
    ax5.set_xticklabels(summaryDf['cellID'].values, rotation = 90)
    ax5.axhline(y = 1.0, color = 'r')
    ax5.set_ylim(0,4)
    plt.savefig(todayFigDir+'/Summary_RatioThickness.png')
    plt.show()
    
    #### Ratio of fluctuations first/last 5 mins after/before activation
    ax6_0 = plt.figure()
    ax6_0 = sns.scatterplot(x = 'cellID',  y = 'ratioFluctuations5min', data=summaryDf, style = 'blebCondition', hue = 'activationType', s = 100)
    ax6_0.set_xlabel('cellID')
    ax6_0.set_ylabel('After/Before Ratio Fluctuations (5 mins before/after)')
    # confInt = addStat_df(ax6, summaryDf, [('Before', 'After')], 'fluctuations', test = 'pairwise', cond = 'activationTag')
    # plt.suptitle("Mean - 1.96*StdError: " + str(confInt))
    ax6_0.set_xticklabels(summaryDf['cellID'].values, rotation = 90)
    ax6_0.axhline(y = 1.0, color = 'r')
    ax6_0.set_ylim(0,4)
    plt.savefig(todayFigDir+'/Summary_RatioFluctuations5mins.png')
    plt.show()
    
    #### Ratio of fluctuations after/before activation
    ax6 = plt.figure()
    ax6 = sns.scatterplot(x = 'cellID',  y = 'ratioFluctuations', data=summaryDf, style = 'blebCondition', hue = 'activationType', s = 100)
    ax6.set_xlabel('cellID')
    ax6.set_ylabel('After/Before Ratio Fluctuations')
    # confInt = addStat_df(ax6, summaryDf, [('Before', 'After')], 'fluctuations', test = 'pairwise', cond = 'activationTag')
    # plt.suptitle("Mean - 1.96*StdError: " + str(confInt))
    ax6.set_xticklabels(summaryDf['cellID'].values, rotation = 90)
    ax6.axhline(y = 1.0, color = 'r')
    ax6.set_ylim(0,4)
    plt.savefig(todayFigDir+'/Summary_RatioFluctuations.png')
    plt.show()
    
    #### Ratio of fluctuation/Ratio of Thickness Before/After
    ax7 = plt.figure()
    ax7 = sns.scatterplot(x = 'cellID',  y = 'ratioFluctThick', data=summaryDf, style = 'blebCondition', hue = 'activationType', s = 100)
    ax7.set_xlabel('cellID')
    ax7.set_ylabel('After/Before Ratio Fluctuations/Ratio Thickness')
    # confInt = addStat_df(ax7, summaryDf, [('Before', 'After')], 'ratioFluctThick', test = 'pairwise', cond = 'activationTag')
    # plt.suptitle("Mean - 1.96*StdError: " + str(confInt))
    ax7.set_xticklabels(summaryDf['cellID'].values, rotation = 90)
    ax7.axhline(y = 1.0, color = 'r')
    ax7.set_ylim(0,4)
    plt.savefig(todayFigDir+'/Summary_RatioFluctThickness.png')
    plt.show()
    
    #### Median Thickness After vs. Median Thickness Before 
    ax8 = plt.figure()
    x = summaryDf['medianThickness'][summaryDf['activationTag'] == 'Before'].values
    y = summaryDf['medianThickness'][summaryDf['activationTag'] == 'After'].values
    hue = summaryDf['activationType'][summaryDf['activationTag'] == 'After'].values
    style = summaryDf['blebCondition'][summaryDf['activationTag'] == 'After'].values
    ax8 = sns.scatterplot(x = x, y = y, hue = hue, style = style, s = 100)
    ax8.set_xlabel('medianThicknessBefore (um)')
    ax8.set_ylabel('medianThicknessAfter (um)')
    plt.savefig(todayFigDir+'/Summary_MedianThicknessBefore-After.png')
    plt.legend()
    plt.show()
    
     #### Fluctuations Before vs. Fluctuation After
    ax9 = plt.figure()
    x = summaryDf['fluctuations'][summaryDf['activationTag'] == 'Before'].values
    y = summaryDf['fluctuations'][summaryDf['activationTag'] == 'After'].values
    hue = summaryDf['activationType'][summaryDf['activationTag'] == 'After'].values
    style = summaryDf['blebCondition'][summaryDf['activationTag'] == 'After'].values
    ax9 = sns.scatterplot(x = x, y = y, hue = hue, style = style, s = 100)
    ax9.set_xlabel('fluctuationsBefore (um)')
    ax9.set_ylabel('fluctuationsAfter (um)')
    plt.savefig(todayFigDir+'/Summary_FluctuationsBefore-After.png')
    plt.legend()
    plt.show()
    
    return(summaryDf)


 #%% Constant field plots
# #%%% Plotting summary of thickness plots
cellDf = pd.read_csv(experimentalDataDir+'/cellConditions.csv', sep=',')
listOfCells = np.asarray(cellDf['cellID'][cellDf['excluded'] == 'no'])
summaryDf = ctFieldThicknessSummary(experimentalDataDir, todayFigDir, listOfCells)


# %% Cloe all open plots
plt.close('all')

#%%% Plotting all three plots (3D, 2D, Dz vs Time) of an experiment
# date = '22.03.01'
# ctFieldThicknessIndividual(experimentalDataDir, figDir, date, save = True, background = 'dark')


# %%
# plt.close('all')

#%%% Plotting 3D trajectories of all cells

# tag = 'global'
# ctFieldThicknessAll(experimentalDataDir, figDir, date, tag = tag, save = True, background = 'default')


# %%
# plt.close('all')


#%% shitty test plots
#%%% Plotting all three graphs (3D, 2D and Dz)

expt = '20220412_100xoil_3t3optorhoa_4.5beads_15mT_Mechanics'
folder = '22-04-12_M1_P1_C5_disc20um'
date = '22.04.12'

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

# %%% Just 3D graphs

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
# plt.xlim(0.5,8)

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
               