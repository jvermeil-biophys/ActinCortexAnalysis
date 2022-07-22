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
dateFormatExcel = re.compile(r'\d{2}/\d{2}/\d{4}')
dateFormatOk = re.compile(r'\d{2}-\d{2}-\d{2}')


# 5. Directories adress

COMPUTERNAME = os.environ['COMPUTERNAME']
if COMPUTERNAME == 'ORDI-JOSEPH':
    mainDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis"
    ownCloudDir = "C://Users//JosephVermeil//ownCloud//ActinCortexAnalysis"
elif COMPUTERNAME == 'LARISA':
    mainDir = "C://Users//Joseph//Desktop//ActinCortexAnalysis"
    ownCloudDir = "C://Users//Joseph//ownCloud//ActinCortexAnalysis"
elif COMPUTERNAME == '':
    mainDir = "C://Users//josep//Desktop//ActinCortexAnalysis"
    ownCloudDir = "C://Users//josep//ownCloud//ActinCortexAnalysis"
elif COMPUTERNAME == 'DATA2JHODR':
    mainDir = "C://Users//BioMecaCell//Desktop//ActinCortexAnalysis"
    ownCloudDir = ""
    

# 6. Add the folder to path
import sys
sys.path.append(mainDir + "//Code_Python")

# 7. Local imports
sys.path.append(mainDir + "//Code_Python")
from BeadTracker import mainTracker
import utilityFunctions_JV as jvu

# 8. Setting of the directories

mainDataDir = 'D://Duya//MagneticPincherData'
rawDataDir = os.path.join(mainDataDir, 'Raw')
depthoDir = os.path.join(rawDataDir, 'EtalonnageZ')
interDataDir = os.path.join(mainDataDir, 'Intermediate')

figureDir = os.path.join(mainDir, 'Figures')
timeSeriesDataDir = os.path.join(mainDir, 'Data_Analysis', 'TimeSeriesData')


try:
    ownCloud_figureDir = os.path.join(ownCloudDir, 'Figures')
    ownCloud_timeSeriesDataDir = os.path.join(ownCloudDir, 'Data_Analysis', 'TimeSeriesData')
except:
    ownCloud_figureDir = ''
    ownCloud_timeSeriesDataDir = ''

# 9. Import of the experimental conditions

experimentalDataDir = os.path.join(mainDir, 'Data_Experimental_DB')
expDf = jvu.getExperimentalConditions(experimentalDataDir, save = True, sep = ';', suffix = '_DB')

# %% Small things

#### Plot last traj

# fig, ax = plt.subplots(1,1)
# X1, Y1 = listTrajDicts[0]['X'], listTrajDicts[0]['Y']
# X2, Y2 = listTrajDicts[1]['X'], listTrajDicts[1]['Y']
# ax.plot(X1, Y1, 'r-')
# ax.plot(X2, Y2, 'b-')

# fig.show()

#### close all
plt.close('all')

# %% Experiments


# %%% 22.06.10_M1

# %%%% 22.06.10_M1; only the first cell
dates = '22.06.10'
manips, wells, cells = 1, 1, 1
depthoNames = '22.06.10_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir) 


# %%%% 22.06.10_M1; all the cells
dates = '22.06.10'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.06.10_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir) 
# %%% 22.06.16_M1

# %%%% 22.06.16_M1 ; only the first cell
dates = '22.06.16'
manips, wells, cells = 1, 1, 1
depthoNames = '22.06.16_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)



# %%%% 22.06.16_M1 ; all the cells
dates = '22.06.16'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.06.16_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)



# %%% 22.07.06

# %%%% 22.07.06_M1 ; only the first cell
dates = '22.07.06'
manips, wells, cells = 1, 1, 3
depthoNames = '22.07.06_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)

# %%%% 22.07.06_M1 ; all the cells
dates = '22.07.06'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.07.06_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)

# %%%% 22.07.06_M2 ; only the first cell
dates = '22.07.06'
manips, wells, cells = 2, 1, 6
depthoNames = '22.07.06_M2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)



# %%%% 22.07.06_M2 ; all the cells
dates = '22.07.06'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.07.06_M2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)

# %%% 22.07.12

# %%%% 22.07.12_M1 ; only the first cell
dates = '22.07.12'
manips, wells, cells = 1, 1, 8
depthoNames = '22.07.12_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)

# %%%% 22.07.12_M1 ; all the cells
dates = '22.07.12'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.07.12_M2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)


# %%%% 22.07.12_M2 ; only the first cell
dates = '22.07.12'
manips, wells, cells = 2, 1, 2
depthoNames = '22.07.12_M2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)

# %%%% 22.07.12_M1 ; all the cells
dates = '22.07.12'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.07.12_M2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)


# %% EXAMPLE FROM JV - Topic : Drugs & perturbation

# %%% Example Date : 22.03.30, compressionsLowStart of 3T3 LG +++, M450, M1 = Blebbi, M2 = LatA, M3 = Ctrl, M4 = DMSO

# %%%% Example Manip : 22.03.30_M1 C1 Seulement

dates = '22.03.30'
manips, wells, cells = 1, 1, 1
depthoNames = '22.03.30_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir)

# %%%% 22.03.30_M1
dates = '22.03.30'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.03.30_M1_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir) 

# %%%% 22.03.30_M2
dates = '22.03.30'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.03.30_M2_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir) 

# %%%% 22.03.30_M3
dates = '22.03.30'
manips, wells, cells = 3, 1, 'all'
depthoNames = '22.03.30_M3_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir) 

# %%%% 22.03.30_M4
dates = '22.03.30'
manips, wells, cells = 4, 1, 'all'
depthoNames = '22.03.30_M4_M450_step20_100X'

timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                                  sourceField = 'default',
                                  ownCloudDir = ownCloudDir, 
                                  ownCloud_figureDir = ownCloud_figureDir, 
                                  ownCloud_timeSeriesDataDir = ownCloud_timeSeriesDataDir) 











