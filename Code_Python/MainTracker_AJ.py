# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 17:50:53 2021

@author: JosephVermeil
"""

# %% General imports

# 1. Imports
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import matplotlib.pyplot as plt

import os
import re
import time
import pyautogui
import matplotlib
# import cv2

# import scipy
from scipy import interpolate
from scipy import signal

# import skimage
from skimage import io, filters, exposure, measure, transform, util
from scipy.signal import find_peaks, savgol_filter
from scipy.optimize import linear_sum_assignment


# 2. Pandas settings
pd.set_option('mode.chained_assignment',None)

# 3. Plot settings
# Here we use this mode because displaying images 
# in new windows is more convenient for this code.
matplotlib.use('Qt5Agg')
# %matplotlib qt 
# To switch back to inline display, use : 
# %matplotlib widget or %matplotlib inline
# matplotlib.rcParams.update({'figure.autolayout': True})

SMALLER_SIZE = 8
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# 4. Other settings
# These regex are used to correct the stupid date conversions done by Excel
dateFormatExcel = re.compile('\d{2}/\d{2}/\d{4}')
dateFormatOk = re.compile('\d{2}-\d{2}-\d{2}')


# 5. Add the folder to path
import sys
sys.path.append("C://Users//anumi//Desktop//ActinCortexAnalysis//Code_Python")

# %% Import of the BeadTracker functions

from BeadTracker import *
import utilityFunctions_JV as jvu


# %% Setting of the directories

mainDataDir = 'D:/Anumita/MagneticPincherData'
rawDataDir = os.path.join(mainDataDir, 'Raw')
depthoDir = os.path.join(rawDataDir, 'EtalonnageZ')
interDataDir = os.path.join(mainDataDir, 'Intermediate')
figureDir = os.path.join(mainDataDir, 'Figures')
timeSeriesDataDir = "C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData"


# %% Import of the experimental conditions

experimentalDataDir = "C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Experimental_AJ"
expDf = jvu.getExperimentalConditions(experimentalDataDir, save = True, sep = ';')

# %% EXAMPLE -- 21.10.18, compressions of 3T3, M1 = M270, M2 = M450
# %%%% M1
dates = '21.10.18'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.10.18_M1_M270_100X_step20'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False)
# %%%% M2
dates = '21.10.18'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.10.18_M2_M450_100X_step20'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False)

# %% Stand alone xyz tracker: To test images from Atchoum with code
# %%%% Test run 1 with 60x objective depthographs from Atchoum

mainDataDir = 'D:/Anumita/Data'
rawDataDir = os.path.join(mainDataDir, 'Raw')
depthoDir = os.path.join(rawDataDir, 'EtalonnageZ')
interDataDir = os.path.join(mainDataDir, 'Intermediate_Py')
figureDir = os.path.join(mainDataDir, 'Figures')
timeSeriesDataDir = "C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData"

imageDir = 'D:/Anumita/Data/Raw/21.11.26'
imageName = 'Cell1_Chamber1_15mT_50ms_15s_Stack.tif'
#imageName =  '0-2-4umStack_M450only_60X.tif'
imagePath = os.path.join(imageDir, imageName)
depthoNames = '21.11.30_M450_step50_60X'

cellID = 'test'
I = io.imread(imagePath).T # trqnspose chqnge if not necessqry
print(I.shape)

manipDict = {}
manipDict['experimentType'] = 'tracking'
manipDict['scale pixel per um'] = 15.8*0.6
manipDict['optical index correction'] = 1.33/1.52
manipDict['magnetic field correction'] = 0
manipDict['beads bright spot delta'] = 0
manipDict['bead type'] = 'M450'
manipDict['bead diameter'] = 4503
manipDict['loop structure'] = '3_0_0'
manipDict['normal field multi images'] = 3
manipDict['multi image Z step'] = 500 #nm
manipDict['with fluo images'] = False

NB =  2
PTL = XYZtracking(I, cellID, NB, manipDict, depthoDir, depthoNames)

# %%%% Test run for 100X Atchoum 21/12/08 (Planes 500nm apart) - Worked shitty. The depthos were not well constructed
# because of bad Kohler illumination

mainDataDir = 'D:/Anumita/Data'
rawDataDir = os.path.join(mainDataDir, 'Raw')
depthoDir = os.path.join(rawDataDir, 'EtalonnageZ')
interDataDir = os.path.join(mainDataDir, 'Intermediate_Py')
figureDir = os.path.join(mainDataDir, 'Figures')
timeSeriesDataDir = "C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData"

imageDir = 'D:/Anumita/Data/Raw/21.12.08'
imageName = 'bead2_7-8-9frames_Stack.tif'
#imageName =  '0-2-4umStack_M450only_60X.tif'
imagePath = os.path.join(imageDir, imageName)
depthoNames = '21.12.08_M450_step50_100X'

cellID = 'test'
I = io.imread(imagePath).T # trqnspose chqnge if not necessqry
print(I.shape)

manipDict = {}
manipDict['experimentType'] = 'tracking'
manipDict['scale pixel per um'] = 15.8
manipDict['optical index correction'] = 1.33/1.52
manipDict['magnetic field correction'] = 0
manipDict['beads bright spot delta'] = 0
manipDict['bead type'] = 'M450'
manipDict['bead diameter'] = 4503
manipDict['loop structure'] = '3_0_0'
manipDict['normal field multi images'] = 3
manipDict['multi image Z step'] = 500 #nm
manipDict['with fluo images'] = False

NB =  1
PTL = XYZtracking(I, cellID, NB, manipDict, depthoDir, depthoNames)

# %%%% Test run 2 for 100X Atchoum 21/12/13 (Planes 500nm apart) - Shitty deptho.

mainDataDir = 'D:/Anumita/Data'
rawDataDir = os.path.join(mainDataDir, 'Raw')
depthoDir = os.path.join(rawDataDir, 'EtalonnageZ')
interDataDir = os.path.join(mainDataDir, 'Intermediate_Py')
figureDir = os.path.join(mainDataDir, 'Figures')
timeSeriesDataDir = "C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData"

imageDir = 'D:/Anumita/Data/Raw/21.12.13'
imageName = 'B3_3-4-5frames_Stack.tif'
#imageName =  '0-2-4umStack_M450only_60X.tif'
imagePath = os.path.join(imageDir, imageName)
depthoNames = '21.12.13_M450_step20_100X'

cellID = 'test'
I = io.imread(imagePath).T # transpose change if not necessary
print(I.shape)

manipDict = {}
manipDict['experimentType'] = 'tracking'
manipDict['scale pixel per um'] = 15.8
manipDict['optical index correction'] = 1.33/1.52
manipDict['magnetic field correction'] = 0
manipDict['beads bright spot delta'] = 0
manipDict['bead type'] = 'M450'
manipDict['bead diameter'] = 4503
manipDict['loop structure'] = '3_0_0'
manipDict['normal field multi images'] = 3
manipDict['multi image Z step'] = 500 #nm
manipDict['with fluo images'] = False

NB =  1
PTL = XYZtracking(I, cellID, NB, manipDict, depthoDir, depthoNames)


# %% OptoPincher test experiments from Atchoum
# %%%% 10/12/2021 : Shitty experiment from 21.12.10

dates = '21.12.10'
manips, wells, cells = 1, 1, 1
depthoNames = '21.12.13_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)

# %% 20/12/2021 : First experiment with the optimised illumtination

# %%%% M1
dates = '21.12.20'
manips, wells, cells = 1, 1, 1
depthoNames = '21.12.20_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% M2
dates = '21.12.20'
manips, wells, cells = 2, 3, 'all'
depthoNames = '21.12.20_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)
# %%%% M3
dates = '21.12.20'
manips, wells, cells = 3, 1, 'all'
depthoNames = '21.12.20_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)

# %% 03/02/2022 :Experiment with the re-optimised illumtination + PLL-PEG 1mg/ml
# Deptho seems a bit strange - probably have to open the secondary aperture more

# %%%% M1
dates = '22.02.03'
manips, wells, cells = 1, 1, 1
depthoNames = '22.02.03_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% M2
dates = '22.02.03'
manips, wells, cells = 2, 1, 2
depthoNames = '22.02.03_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)
# %%%% M3
dates = '22.02.03'
manips, wells, cells = 3, 1, 'all'
depthoNames = '22.02.03_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% M4
dates = '22.02.03'
manips, wells, cells = 4, 1, 'all'
depthoNames = '22.02.03_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% M5 - Riceball activation
dates = '22.02.03'
manips, wells, cells = 5, 1, 'all'
depthoNames = '22.02.03_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)

# %% 01/03/2022 :Experiment with the optimised activation parameters (1.2microWatts) and 1mg/ml PLL-PEG coated beads to prevent engulfent

# %%%% M1 : Global activation

dates = '22.03.01'
manips, wells, cells = 1, 1, 5
depthoNames = '22.03.01_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M2 : Activation away from beads

dates = '22.03.01'
manips, wells, cells = 2, 'all', 'all'
depthoNames = '22.03.01_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M3 : Activation at beads

dates = '22.03.01'
manips, wells, cells = 3, 1, 6
depthoNames = '22.03.01_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %% 22/03/2022 : Experiment with the optimised activation parameters (1.2microWatts) and 1mg/ml PLL-PEG coated beads to prevent engulfment
# To obtain as many curves as possible

# %%%% M1 : Global activation, 60s frequency

dates = '22.03.22'
manips, wells, cells = 1, 1, 1
depthoNames = '22.03.22_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M2 : Global activation, 30s frequency

dates = '22.03.22'
manips, wells, cells = 2, 'all', 'all'
depthoNames = '22.03.22_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M3 : Half activation, 30s frequency, at beads

dates = '22.03.22'
manips, wells, cells = 3, 'all', 'all'
depthoNames = '22.03.22_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M4 : Half activation, 30s frequency, away from beads

dates = '22.03.22'
manips, wells, cells = 4, 'all', 'all'
depthoNames = '22.03.22_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M5 : Global activation, 30s frequency, fixed duration

dates = '22.03.22'
manips, wells, cells = 5, 2, 4
depthoNames = '22.03.22_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %% 09/05/2022 :

# %%%% M1 : At beads - might have stopped activation after 3mins (wrong version of Metamorph state loaded by accident)

dates = '22.05.09'
manips, wells, cells = 1, 1, 1
depthoNames = '22.05.09_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M2 : Global activation, 60s frequency

dates = '22.05.09'
manips, wells, cells = 2, 2, 'all'
depthoNames = '22.05.09_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M3 : Global activation, 60s frequency

dates = '22.05.09'
manips, wells, cells = 3, 2, 'all'
depthoNames = '22.05.09_P2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M4 : Global activation, 60s frequency

dates = '22.05.09'
manips, wells, cells = 4, 2, 'all'
depthoNames = '22.05.09_P2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M5 : Global activation, 60s frequency

dates = '22.05.09'
manips, wells, cells = 5, 1, 'all'
depthoNames = '22.05.09_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M6 : Global activation, 60s frequency

dates = '22.05.09'
manips, wells, cells = 6, 2, 'all'
depthoNames = '22.05.09_P2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)


# %% 12/04/2022 : First constant field expeirments in PMMH. 0.6uW power.

# %%%% M1 : Half activation, away from beads, 50ms, 30s frequency

dates = '22.04.12'
manips, wells, cells = 1, 'all', 'all'
depthoNames = '22.04.12_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

  # %%%% M2 : Half activation, at beads, 50ms, 30s frequency

dates = '22.04.12'
manips, wells, cells = 2, 'all', 'all'
depthoNames = '22.04.12_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

  # %%%% M3 : Global activation, 50ms, 30s frequency

dates = '22.04.12'
manips, wells, cells = 3, 'all', 'all'
depthoNames = '22.04.12_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %% 09/06/2022 : 3.2 uW Power

# %%%% M1 : Half activation, away from beads, 50ms, 30s frequency

dates = '22.06.09'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.06.09_P'+str(wells)+'_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

  # %%%% M2 : Half activation, at beads, 50ms, 30s frequency

dates = '22.06.09'
manips, wells, cells = 2, 3, 'all'
depthoNames = '22.06.09_P'+str(wells)+'_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

  # %%%% M3 : Global activation, 50ms, 30s frequency

dates = '22.06.09'
manips, wells, cells = 3, 3, 'all'
depthoNames = '22.06.09_P'+str(wells)+'_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)
# %% 31/03/2022 : Experiment in PMMH Mechanics:

# %%%% M6 : Half activation, At beads, 500ms once

dates = '22.03.31'
manips, wells, cells = 6, 2, 3
depthoNames = '22.03.31_P1_M450_step20_100X'
  
timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M7 : Half activation, Away from beads, 500ms once

dates = '22.03.31'
manips, wells, cells = 7, 'all', 'all'
depthoNames = '22.03.31_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M8 : Global activation, 500ms once

dates = '22.03.31'
manips, wells, cells = 8, 2, 3
depthoNames = '22.03.31_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M9 : Global activation, 800ms once

dates = '22.03.31'
manips, wells, cells = 9, 2, 2
depthoNames = '22.03.31_P1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %% 28/04/2022 : Experiment in PMMH Mechanics:

# %%%% M1 : Half activation, away from beads, 500ms first followed by 10ms

dates = '22.04.28'
manips, wells, cells = 1, 'all', 'all'
depthoNames = '22.04.28_P1_M450_step20_100X'
  
timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %% 31/05/2022 : Experiment in PMMH Mechanics

# %%%% M1 : At beads activation, level 3 fluo intensity, filter 6, initial 500ms 
# with 50ms activation at the end of every loop

dates = '22.05.31'
manips, wells, cells = 6, 'all', 'all'
depthoNames = '22.05.31_P1_M450_step20_100X'
  
timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M4 : At beads activation, level 3 fluo intensity, filter 6, initial 500ms 
# with 50ms activation at the end of every loop

dates = '22.05.31'
manips, wells, cells = 4, 1, 'all'
depthoNames = '22.05.31_P'+str(wells)+'_M450_step20_100X'
  
timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M5 : At beads activation, level 3 fluo intensity, filter 6, initial 500ms 
# with 50ms activation at the end of every loop

dates = '22.05.31'
manips, wells, cells = 5, 1, 4
depthoNames = '22.05.31_P2_M450_step20_100X'
  
timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %%%% M7 : At beads activation, level 3 fluo intensity, filter 6, initial 500ms 
# with 50ms activation at the end of every loop

dates = '22.05.31'
manips, wells, cells = 7, 1, 8
depthoNames = '22.05.31_P'+str(wells)+'_M450_step20_100X'
  
timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)

# %% 21/06/2022 : Experiment in PMMH Mechanics

# %%%% M1 : At beads activation, level 3 fluo intensity, filter 6, initial 500ms 
# with 50ms activation at the end of every loop

dates = '22.06.21'
manips, wells, cells = 2, 2, 2
depthoNames = '22.06.21_M450_step20_100X' 
  
timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True)
