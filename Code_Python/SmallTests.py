# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 13:37:51 2021

@author: JosephVermeil
"""



# %% Chi2 test

# %%% Imports

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as st
import statsmodels.api as sm

import os
import time
import random
import traceback


from copy import copy
from cycler import cycler
from datetime import date
from scipy.optimize import curve_fit
# from statannot import add_stat_annotation

pd.set_option('mode.chained_assignment',None)
pd.set_option('display.max_columns', None)

import re
dateFormatExcel = re.compile('\d{2}/\d{2}/\d{4}')
dateFormatOk = re.compile('\d{2}-\d{2}-\d{2}')

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'figure.autolayout': True})

# Add the folder to path
import sys
sys.path.append("C://Users//JosephVermeil//Desktop//ActinCortexAnalysis//Code_Python")
from getExperimentalConditions import getExperimentalConditions
from PincherAnalysis_JV import *

# 5. Global constants
SCALE_100X = 15.8 # pix/µm 
NORMAL  = '\033[0m'
RED  = '\033[31m' # red
GREEN = '\033[32m' # green
ORANGE  = '\033[33m' # orange
BLUE  = '\033[36m' # blue

# %%% Directories

mainDataDir = 'D://MagneticPincherData'
rawDataDir = os.path.join(mainDataDir, 'Raw')
depthoDir = os.path.join(rawDataDir, 'EtalonnageZ')
interDataDir = os.path.join(mainDataDir, 'Intermediate')
figureDir = os.path.join(mainDataDir, 'Figures')
timeSeriesDataDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis//Data_Analysis//TimeSeriesData"
experimentalDataDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis//Data_Experimental"

# %%% Data

expDf = getExperimentalConditions(experimentalDataDir)

cellId = '21-12-16_M1_P1_C10'
df = getCellTimeSeriesData(cellId)



# %%% Functions

def segmentTimeSeries_meca(f, tsDF, expDf):
    #### (0) Import experimental infos
    split_f = f.split('_')
    tsDF.dx, tsDF.dy, tsDF.dz, tsDF.D2, tsDF.D3 = tsDF.dx*1000, tsDF.dy*1000, tsDF.dz*1000, tsDF.D2*1000, tsDF.D3*1000
    thisManipID = split_f[0] + '_' + split_f[1]
    expDf['manipID'] = expDf['date'] + '_' + expDf['manip']
    thisExpDf = expDf.loc[expDf['manipID'] == thisManipID]

    diameters = str(thisExpDf.at[thisExpDf.index.values[0], 'bead diameter'])
    diameters = diameters.split('_')
    DIAMETER = int(diameters[0])
    EXPTYPE = str(thisExpDf.at[thisExpDf.index.values[0], 'experimentType'])

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


    compField = thisExpDf.at[thisExpDf.index.values[0], 'ramp field'].split('_')
    minCompField = int(compField[0])
    maxCompField = int(compField[1])
    
    hComprList = []
    fComprList = []


    for i in range(1, Ncomp+1): #Ncomp+1):

        #### (1) Segment the compression n°i
        thisCompDf = tsDF.loc[tsDF['idxAnalysis'] == i,:]
        iStart = (findFirst(tsDF['idxAnalysis'], i))
        iStop = iStart+thisCompDf.shape[0]


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
            #### (2) Inside the compression n°i, delimit the compression and relaxation phases

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
            fCompr = (thisCompDf.F.values[jStart:jMax+1])

            # Refinement of the compression delimitation.
            # Remove the 1-2 points at the begining where there is just the viscous relaxation of the cortex
            # because of the initial decrease of B and the cortex thickness increases.
            offsetStart2 = 0
            k = 0
            while (k<len(hCompr)-10) and (hCompr[k] < np.max(hCompr[k+1:min(k+10, len(hCompr))])):
                offsetStart2 += 1
                k += 1
            # Better compressions arrays
            hCompr = np.array(hCompr[offsetStart2:])
            fCompr = np.array(fCompr[offsetStart2:])

            hComprList.append(hCompr)
            fComprList.append(fCompr)
            
    return(DIAMETER, hComprList, fComprList)



def compressionFitChadwick_V2(hCompr, fCompr, DIAMETER):
    
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
        
        residuals = hCompr-hPredict

        SSR = np.sum((residuals)**2)
        alpha = 0.975
        dof = len(fCompr)-len(params)
        q = st.t.ppf(alpha, dof) # Student coefficient
        R2 = get_R2(hCompr, hPredict)
        Chi2 = get_Chi2(hCompr, hPredict, dof)

        varE = covM[0,0]
        seE = (varE)**0.5
        E, seE = E*1e6, seE*1e6
        confIntE = [E-q*seE, E+q*seE]
        confIntEWidth = 2*q*seE

        varH0 = covM[1,1]
        seH0 = (varH0)**0.5
        confIntH0 = [H0-q*seH0, H0+q*seH0]
        confIntH0Width = 2*q*seH0
        
        
    except Exception:
        print(RED + '')
        traceback.print_exc()
        print('\n')
        error = True
        E, H0, hPredict, R2, confIntE, confIntH0 = -1, -1, np.ones(len(hCompr))*(-1), -1, [-1,-1], [-1,-1]
        print(NORMAL + '')
        
    # return(E, H0, hPredict, R2, X2, confIntE, confIntH0, error)
    return(E, H0, hPredict, R2, Chi2, residuals, error)




def get_Chi2(Ymeas, Ymodel, dof):
    residuals = Ymeas-Ymodel
    S = st.tstd(residuals)
    S = (np.sum(residuals**2)/len(residuals))**0.5
    Chi2 = np.sum((residuals/S)**2)
    Chi2_dof = Chi2/dof
    return(Chi2_dof)

def testFun(cellId):
    
    #### MAIN PART [1/2]
    expDf = getExperimentalConditions(experimentalDataDir)
    tsDF = getCellTimeSeriesData(cellId)
    DIAMETER, hComprList, fComprList = segmentTimeSeries_meca(cellId, tsDF, expDf)
    Nc = len(hComprList)

    
    
    #### PLOT [1/2]
    nColsSubplot = 5
    nRowsSubplot = ((Nc-1) // nColsSubplot) + 1
    # First plot - fig & ax, gather all the F(h) curves, and will be completed later in the code
    fig, ax = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
    # Second plot - fig2 & ax2, 
    fig2, ax2 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
    
    
    
    #### MAIN PART [2/2]
    for i in range(1,Nc+1): 
        hCompr, fCompr = hComprList[i-1], fComprList[i-1]
        E, H0, hPredict, R2, Chi2, residuals, fitError = compressionFitChadwick_V2(hCompr, fCompr, DIAMETER)
        
        NT = st.mstats.normaltest(residuals)
        # # print(NT)
        # S = st.tstd(residuals)
        print(R2, Chi2)
    
        #### PLOT [2/2]
        colSp = (i-1) % nColsSubplot
        rowSp = (i-1) // nColsSubplot

        # ax2[i-1] with the 1 line plot
        if nRowsSubplot == 1:
            thisAx = ax[colSp]
        elif nRowsSubplot >= 1:
            thisAx = ax[rowSp,colSp]
            
        if nRowsSubplot == 1:
            thisAx2 = ax2[colSp]
        elif nRowsSubplot >= 1:
            thisAx2 = ax2[rowSp,colSp]
        
        thisAx.plot(hCompr,fCompr,'b-', linewidth = 0.8)
        titleText = cellId + '__c' + str(i)
        legendText = ''
        thisAx.set_xlabel('h (nm)')
        thisAx.set_ylabel('f (pN)')
        
        legendText2 = '' + 'Normality test\nStat = {:.3f}\np-val = {:.3f}'.format(NT.statistic, NT.pvalue)
        thisAx2.plot(fCompr,residuals, 'b+', label = legendText2)
        thisAx2.legend(loc = 'upper right', prop={'size': 6})
        titleText = cellId + '__c' + str(i)
        thisAx2.set_xlabel('f (pN)')
        thisAx2.set_ylim([-15,+15])
        thisAx2.set_ylabel('residuals')
    
        if not fitError:
            legendText += 'H0 = {:.1f}nm\nE = {:.2e}Pa'.format(H0, E)
            thisAx.plot(hPredict,fCompr,'k--', linewidth = 0.8, label = legendText)
            thisAx.legend(loc = 'upper right', prop={'size': 6})

            
        else:
            titleText += '\nFIT ERROR'
    
        thisAx.title.set_text(titleText)
    
        for item in ([thisAx.title, thisAx.xaxis.label, \
                      thisAx.yaxis.label] + thisAx.get_xticklabels() + thisAx.get_yticklabels()):
            item.set_fontsize(9)
        for item in ([thisAx2.title, thisAx2.xaxis.label, \
                      thisAx2.yaxis.label] + thisAx2.get_xticklabels() + thisAx2.get_yticklabels()):
            item.set_fontsize(9)

    plt.show()

# %%% Main script

cellId = '21-12-16_M1_P1_C10' # Good fits
cellId = '21-12-16_M2_P1_C6' # Bad fits
testFun(cellId)


# %%% End

plt.close('all')






# %% Next test !


# %% Next test !


# %% Next test !






