# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 13:37:51 2021

@author: JosephVermeil
"""

# %%
import  pandas as pd
import numpy as np

dicta={'aaa' : [1,2,3,4], 'bbb' : ['p','o','p','p'], 'ccc' : ['q','q','w','w']}
df = pd.DataFrame(dicta)
N = df.loc[df['bbb'] == 'p'].shape[0]
# A = np.array(['youpi' for kk in range(N)])
df.loc[df['bbb'] == 'p', 'ccc'] = 'youpo'

# %%
import numpy as np


def myfun(a):
    a = a + [2]

a = [1]

print(a)

myfun(a)

print(a)


# %% Fastness test
import numpy as np
import time

N = int(1e4)
a = np.zeros((N,N))
b = np.ones((N,N))

T0 = time.time()
# c = (b-a)**2
c =  a.T
T1 = time.time()
# np.square(np.subtract(b,a))
T2 = time.time()

print(T1-T0)
print(T2-T1)

# %% Fastness test
import numpy as np
import time

N = int(2e4)
a = np.zeros((N,N))
b = np.ones((N,N))
A = np.repeat(np.arange(1,N+1), N, axis = 0).reshape(N,N)

T0 = time.time()
m = A / np.mean(A,axis=1)[:,None]
T1 = time.time()
m2 = (A.T / np.mean(A,axis=1)).T
T2 = time.time()
m3 = np.divide(A.T,np.mean(A,axis=1)).T
T3 = time.time()

print(T1-T0)
print(T2-T1)
print(T3-T2)

# %% Fastness test
import numpy as np
import time

def squareDistance(M, V, normalize = False): # MUCH FASTER ! **Michael Scott Voice** VERRRRY GOODE
    """
    Compute a distance between two arrays of the same size, defined as such:
    D = integral of the squared difference between the two arrays.
    It is used to compute the best fit of a slice of a bead profile on the depthograph.
    This function speed is critical for the Z computation process because it is called so many times !
    What made that function faster is the absence of 'for' loops and the use of np.repeat().
    """
    # squareDistance(self.deptho, listProfiles[i], normalize = True)
    # top = time.time()
    # n, m = M.shape[0], M.shape[1]
    # len(V) should be m
    if normalize:
        V = V/np.mean(V)
    MV = np.repeat(np.array([V]), M.shape[0], axis = 0) # Key trick for speed !
    if normalize:
        M = M / np.mean(M,axis=1)[:,None]
    R = np.sum(np.square(np.subtract(M,MV)), axis = 1)
    # print('DistanceCompTime')
    # print(time.time()-top)
    return(R)

N = int(1e4)
a = np.ones((N,N))
b = np.ones(N)*5

T0 = time.time()
R1 = squareDistance(a,b)
T1 = time.time()
R2 = squareDistance(a,b, normalize = True)
T2 = time.time()

T3 = time.time()

print(T1-T0)
print(T2-T1)
print(T3-T2)

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

cellId = '21-12-08_M1_P2_C1'
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
        residuals_h = hCompr-hPredict
        
        comprMat = np.array([hCompr, fCompr]).T
        comprMatSorted = comprMat[comprMat[:, 0].argsort()]
        hComprSorted, fComprSorted = comprMatSorted[:, 0], comprMatSorted[:, 1]
        fPredict = chadwickModel(hComprSorted, E, H0)
        residuals_f = fComprSorted-fPredict
        
        deltaCompr = (H0 - hCompr)/1000
        stressCompr = fCompr / (np.pi * (DIAMETER/2000) * deltaCompr)
        strainCompr = deltaCompr / (3*(H0/1000))
        strainPredict = stressCompr / (E*1e6)
        
        # Stats
        alpha = 0.975
        dof = len(fCompr)-len(params)
        q = st.t.ppf(alpha, dof) # Student coefficient
        R2 = get_R2(hCompr, hPredict)
        #### err
        err = 0.02
        
        # print('fCompr')
        # print(fCompr)
        # print('stressCompr')
        # print(stressCompr)
        # print('strainCompr')
        # print(strainCompr)
        # print('strainPredict')
        # print(strainPredict)

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
        
        
    except Exception:
        print(RED + '')
        traceback.print_exc()
        print('\n')
        error = True
        E, H0, hPredict, R2, confIntE, confIntH0 = -1, -1, np.ones(len(hCompr))*(-1), -1, [-1,-1], [-1,-1]
        print(NORMAL + '')
        
    # return(E, H0, hPredict, R2, X2, confIntE, confIntH0, error)
    return(E, H0, hPredict, fPredict, comprMatSorted, R2, Chi2, residuals_h, residuals_f, error)




def get_Chi2(Ymeas, Ymodel, dof, err = 10):
    residuals = Ymeas-Ymodel
    # S = st.tstd(residuals)
    # S = (np.sum(residuals**2)/len(residuals))**0.5
    Chi2 = np.sum((residuals/err)**2)
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
    # Third plot - fig3 & ax3, 
    fig3, ax3 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
    
    
    
    #### MAIN PART [2/2]
    for i in range(1,Nc+1): 
        hCompr, fCompr = hComprList[i-1], fComprList[i-1]
        E, H0, hPredict, fPredict, comprMatSorted, R2, Chi2, residuals_h, residuals_f, fitError = compressionFitChadwick_V2(hCompr, fCompr, DIAMETER)
        hComprSorted, fComprSorted = comprMatSorted[:, 0], comprMatSorted[:, 1]
        
        deltaCompr = (H0 - hCompr)/1000 # µm
        stressCompr = fCompr / (np.pi * (DIAMETER/2000) * deltaCompr)
        strainCompr = deltaCompr / (3*(H0/1000)) 
        strainPredict = stressCompr / E #((H0 - hPredict)/1000) / (3*(H0/1000))
        
        # Stats
        print(i, R2, Chi2)
    
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
            
        if nRowsSubplot == 1:
            thisAx3 = ax3[colSp]
        elif nRowsSubplot >= 1:
            thisAx3 = ax3[rowSp,colSp]
            
        titleText = cellId + '__c' + str(i)
        
        thisAx.plot(hCompr,fCompr,'bo', ls='', markersize = 2)
        # thisAx.plot(hCompr,fCompr,'b-', linewidth = 0.8)
        # thisAx.plot(hComprSorted,fComprSorted,'r-', linewidth = 0.8)
        legendText = ''
        thisAx.set_xlabel('h (nm)')
        thisAx.set_ylabel('f (pN)')
        
        thisAx2.plot(stressCompr, strainPredict - strainCompr, 'b+')
        # thisAx2.legend(loc = 'upper right', prop={'size': 6})
        thisAx2.set_xlabel('h (pN)')
        # thisAx2.set_ylim([-50,+50])
        thisAx2.set_ylabel('residuals')
        
        thisAx3.plot(stressCompr,strainCompr,'kP', ls='', markersize = 4)
        legendText3 = ''
        thisAx3.set_xlabel('sigma (Pa)')
        thisAx3.set_ylabel('epsilon')
    
        if not fitError:
            legendText += 'H0 = {:.1f}nm\nE = {:.2e}Pa'.format(H0, E)
            # thisAx.plot(hComprSorted,fPredict,'y--', linewidth = 0.8, label = legendText)
            thisAx.plot(hPredict,fCompr,'k--', linewidth = 0.8, label = legendText)
            thisAx.legend(loc = 'upper right', prop={'size': 6})
            
            legendText3 += 'H0 = {:.1f}nm\nE = {:.2e}Pa'.format(H0, E)
            thisAx3.plot(stressCompr,strainPredict,'r+', linewidth = 0.8, label = legendText3)
            thisAx3.legend(loc = 'upper right', prop={'size': 6})

            
        else:
            titleText += '\nFIT ERROR'
    
        thisAx.title.set_text(titleText)
        
        axes = [thisAx, thisAx2, thisAx3]
        for axe in axes:
            axe.title.set_text(titleText)
            for item in ([axe.title, axe.xaxis.label, axe.yaxis.label] \
                         + axe.get_xticklabels() + axe.get_yticklabels()):
                item.set_fontsize(9)

    plt.show()

# %%% Main script


# cellId = '21-12-08_M2_P1_C1' # Best fits
# cellId = '21-12-08_M1_P2_C1' # Excellent fits
# cellId = '21-12-16_M1_P1_C10' # Good fits
# cellId = '21-12-16_M2_P1_C6' # Bad fits
# cellId = '21-01-21_M1_P1_C2' # Old fits with large h
cellId = '21-12-16_M1_P1_C3' # Bad fits
# cellId = '21-12-16_M1_P1_C2' # Bad fits
# cellId = '22-01-12_M1_P1_C2' # Recent fits
cellId = '21-01-18_M1-1_P1_C1'
testFun(cellId)


# %%% End

plt.close('all')






# %% compressionFitChadwick_StressStrain

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
    minCompField = float(compField[0])
    maxCompField = float(compField[1])
    
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
        print(params)
        print(covM)
        # Previously I fitted with y=F and x=H, but it didn't work so well cause H(t) isn't monotonous:
        # params, covM = curve_fit(chadwickModel, hCompr, fCompr, initialParameters, bounds = parameterBounds)
        # Fitting with the 'inverse Chadwick model', with y=H and x=F is more convenient

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
    
    return(stressCompr, strainCompr, K, strainPredict, R2, Chi2, confIntK, error)






def testFun(cellId):
    
    def chadwickModel(h, E, H0):
        R = DIAMETER/2
        f = (np.pi*E*R*((H0-h)**2))/(3*H0)
        return(f)

    def inversedChadwickModel(f, E, H0):
        R = DIAMETER/2
        h = H0 - ((3*H0*f)/(np.pi*E*R))**0.5
        return(h)
    
    def computeStress(f, h, H0, D):
        R = D/2000
        delta = (H0 - h)/1000
        stress = f / (np.pi * R * delta)
        return(stress)
        
    def computeStrain(h, H0):
        delta = (H0 - h)/1000
        strain = delta / (3 * (H0/1000))
        return(strain)
    
    #### MAIN PART [1/2]
    expDf = getExperimentalConditions(experimentalDataDir)
    tsDF = getCellTimeSeriesData(cellId)
    DIAMETER, hComprList, fComprList = segmentTimeSeries_meca(cellId, tsDF, expDf)
    Nc = len(hComprList)    
    
    #### PLOT [1/2]
    nColsSubplot = 5
    nRowsSubplot = ((Nc-1) // nColsSubplot) + 1
    
    # 1st plot - fig & ax
    fig, ax = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
    
    # 2nd plot - fig & ax
    fig2, ax2 = plt.subplots(nRowsSubplot,nColsSubplot,figsize=(3*nColsSubplot,3*nRowsSubplot))
    
    
    #### MAIN PART [2/2]
    for i in range(1,Nc+1): 
        hCompr, fCompr = hComprList[i-1], fComprList[i-1]
        
        
        compressionStart_NbPts = 15
        hCompr_start = hCompr[:compressionStart_NbPts]
        fCompr_start = fCompr[:compressionStart_NbPts]
        
        E, H0, hPredict, R2, Chi2, confIntE, confIntH0, errorH = compressionFitChadwick(hCompr_start, fCompr_start, DIAMETER)
        
        stressCompr, strainCompr, K, strainPredict, R2, Chi2, confIntK, errorK = compressionFitChadwick_StressStrain(hCompr, fCompr, H0, DIAMETER)
        
        strainCompr = computeStrain(hCompr, H0)
        stressCompr = computeStress(fCompr, hCompr, H0, DIAMETER)
        
        #### PLOT [2/2]
        # 
        min_f = np.min(fCompr_start)
        low_f = np.linspace(0, min_f, 20)
        low_h = inversedChadwickModel(low_f, E/1e6, H0)
        
        # 
        colSp = (i-1) % nColsSubplot
        rowSp = (i-1) // nColsSubplot

        # 
        if nRowsSubplot == 1:
            thisAx = ax[colSp]
            thisAx2 = ax2[colSp]
        elif nRowsSubplot >= 1:
            thisAx = ax[rowSp,colSp]
            thisAx2 = ax2[rowSp,colSp]

            
        titleText = cellId + '__c' + str(i)
        
        
        thisAx2.plot(hCompr,fCompr,'bo', ls='', markersize = 2)
        legendText2 = ''
        thisAx2.set_xlabel('h (nm)')
        thisAx2.set_ylabel('f (pN)')
    
        if not errorH:
            legendText2 += 'H0 = {:.2f}nm'.format(H0)
            plot_startH = np.concatenate((low_h, hPredict))
            plot_startF = np.concatenate((low_f, fCompr_start))
            thisAx2.plot(plot_startH[0], plot_startF[0],'ro', markersize = 5)
            thisAx2.plot(plot_startH, plot_startF,'r--', linewidth = 0.8, label = legendText2)
            thisAx2.legend(loc = 'upper right', prop={'size': 6})
            
            
        thisAx.plot(stressCompr, strainCompr,'ko', ls='', markersize = 4)
        legendText = ''
        thisAx.set_xlabel('sigma (Pa)')
        thisAx.set_ylabel('epsilon')
            
        if not errorK:
            legendText += 'K = {:.2e}Pa'.format(K)
            thisAx.plot(stressCompr, strainPredict,'r--', linewidth = 0.8, label = legendText)
            thisAx.legend(loc = 'upper right', prop={'size': 6})

            
        else:
            titleText += '\nFIT ERROR'
    
        thisAx2.title.set_text(titleText)
        thisAx.title.set_text(titleText)
        
        axes = [thisAx, thisAx2]
        for axe in axes:
            axe.title.set_text(titleText)
            for item in ([axe.title, axe.xaxis.label, axe.yaxis.label] \
                         + axe.get_xticklabels() + axe.get_yticklabels()):
                item.set_fontsize(9)

    plt.show()
        
# %%% Main script


cellId = '22-02-09_M1_P1_C1'
testFun(cellId)


# %%% Check

cellId = '22-02-09_M1_P1_C1'
f = '22-02-09_M1_P1_C1_L40_20umdisc_PY.csv'
tsDF = getCellTimeSeriesData(cellId)

D, H, F = segmentTimeSeries_meca(f, tsDF, expDf)

# %%% End

plt.close('all')    


# %% Test plots

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

# %%% Imports

N = 3
fig, ax = plt.subplots(N, N, figsize = (3*N, 3*N))
for ii in range(N):
    for jj in range(N):
        color = colorList30[3*ii + jj]
        marker = markerList10[3*ii + jj]
        ax[ii,jj].plot([1,2,3,4], [1,2,3,4], color = color, marker = marker, mec = 'k', ms = 8, ls = '', label = str(3*ii + jj))
        ax[ii,jj].set_xlim([0,5])
        ax[ii,jj].set_ylim([0,5])
fig.legend(loc='upper center', bbox_to_anchor=(0.5,0.95), ncol = N*N)
fig.suptitle('test')
# fig.tight_layout()
fig.show()

# %% Color test !

NORMAL  = '\033[0m'
RED  = '\033[31m' # red
GREEN = '\033[32m' # green
ORANGE  = '\033[33m' # orange
BLUE  = '\033[36m' # blue
YELLOW = '\033[1;33m'

print("\033[0;37;48m Normal text\n")
print("\033[2;37;48m Underlined text\033[0;37;48m \n")
print("\033[1;37;48m Bright Colour\033[0;37;48m \n")
print("\033[3;37;48m Negative Colour\033[0;37;48m \n")
print("\033[5;37;48m Negative Colour\033[0;37;48m\n")
print("\033[1;37;40m \033[2;37:40m TextColour BlackBackground          TextColour GreyBackground                WhiteText ColouredBackground\033[0;37;40m\n")
print("\033[1;30;40m Dark Gray      \033[0m 1;30;40m            \033[0;30;47m Black      \033[0m 0;30;47m               \033[0;37;41m Black      \033[0m 0;37;41m")
print("\033[1;31;40m Bright Red     \033[0m 1;31;40m            \033[0;31;47m Red        \033[0m 0;31;47m               \033[0;37;42m Black      \033[0m 0;37;42m")
print("\033[1;32;40m Bright Green   \033[0m 1;32;40m            \033[0;32;47m Green      \033[0m 0;32;47m               \033[0;37;43m Black      \033[0m 0;37;43m")
print("\033[1;33;48m Yellow         \033[0m 1;33;48m            \033[0;33;47m Brown      \033[0m 0;33;47m               \033[0;37;44m Black      \033[0m 0;37;44m")
print("\033[1;34;40m Bright Blue    \033[0m 1;34;40m            \033[0;34;47m Blue       \033[0m 0;34;47m               \033[0;37;45m Black      \033[0m 0;37;45m")
print("\033[1;35;40m Bright Magenta \033[0m 1;35;40m            \033[0;35;47m Magenta    \033[0m 0;35;47m               \033[0;37;46m Black      \033[0m 0;37;46m")
print("\033[1;36;40m Bright Cyan    \033[0m 1;36;40m            \033[0;36;47m Cyan       \033[0m 0;36;47m               \033[0;37;47m Black      \033[0m 0;37;47m")
print("\033[1;37;40m White          \033[0m 1;37;40m            \033[0;37;40m Light Grey \033[0m 0;37;40m               \033[0;37;48m Black      \033[0m 0;37;48m")

# %% Next test !

print(YELLOW + 'yellow' + NORMAL)

