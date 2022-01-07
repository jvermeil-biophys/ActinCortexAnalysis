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
dateFormatExcel = re.compile('\d{2}/\d{2}/\d{4}')
dateFormatOk = re.compile('\d{2}-\d{2}-\d{2}')


# 5. Add the folder to path
import sys
sys.path.append("C://Users//JosephVermeil//Desktop//ActinCortexAnalysis//Code_Python")

# %% Import of the BeadTracker functions

from BeadTracker import *


# %% Setting of the directories

mainDataDir = 'D://MagneticPincherData'
rawDataDir = os.path.join(mainDataDir, 'Raw')
depthoDir = os.path.join(rawDataDir, 'EtalonnageZ')
interDataDir = os.path.join(mainDataDir, 'Intermediate')
figureDir = os.path.join(mainDataDir, 'Figures')
timeSeriesDataDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis//Data_Analysis//TimeSeriesData"

# %% Import of the experimental conditions

experimentalDataDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis//Data_Experimental"
expDf = getExperimentalConditions(experimentalDataDir, save = True, sep = ';')

# %%%% Plot last traj

PTL.listFrames[-500].show()
for iB in range(PTL.NB):
    T = PTL.listTrajectories[iB]
    ax = plt.gca()
    T.plot(ax, iB)


# fig, ax = plt.subplots(1,1)
# X1, Y1 = listTrajDicts[0]['X'], listTrajDicts[0]['Y']
# X2, Y2 = listTrajDicts[1]['X'], listTrajDicts[1]['Y']
# ax.plot(X1, Y1, 'r-')
# ax.plot(X2, Y2, 'b-')
    
# fig.show()

# %%%%
plt.close('all')

# %%% Next manipe
# %%%% 

# %%% Next manipe
# %%%% 

# %%% 21.12.16, compressions of 3T3, M1 = M450, M2 = M270
# %%%% 21.12.16_M1 C1 Seulement
dates = '21.12.16'
manips, wells, cells = 1, 1, 9
depthoNames = '21.12.16_M1_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False)


# %%%% 21.12.16_M1
dates = '21.12.16'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.12.16_M1_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% 21.12.16_M2
dates = '21.12.16'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.12.16_M2_M270_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False, MatlabStyle = True)


# %%% 21.12.08, compressions of 3T3, M1 = M270, M2 = M450
# %%%% 21.12.08_M1 C1 Seulement
dates = '21.12.08'
manips, wells, cells = 1, 2, 4
depthoNames = '21.12.08_M1_M270_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False)


# %%%% 21.12.08_M1
dates = '21.12.08'
manips, wells, cells = 1, 2, 'all'
depthoNames = '21.12.08_M1_M270_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False)

# %%%% 21.12.08_M2
dates = '21.12.08'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.12.08_M2_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False)


# %%% 21.10.25, compressions of 3T3, M1 = M450, M2 = M270
# %%%% 21.10.25_M1 C1 Seulement
dates = '21.10.25'
manips, wells, cells = 1, 1, 1
depthoNames = '21.10.25_M1_M450_100X_step20'

listTrajDicts, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False)

# %%%% 21.10.25_M1
dates = '21.10.25'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.10.25_M1_M450_100X_step20'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False)

# %%%% 21.10.25_M2
dates = '21.10.25'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.10.25_M2_M270_100X_step20'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False)

# %%% 21.10.18, compressions of 3T3, M1 = M270, M2 = M450
# %%%% 21.10.18_M1
dates = '21.10.18'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.10.18_M1_M270_100X_step20'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.7, 
                                  redoAllSteps = False)

# %%%% 21.10.18_M2
dates = '21.10.18'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.10.18_M2_M450_100X_step20'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False)




# %% Older manipes

# %%% Next manipe
# %%%% 

# %%% 21.01.21, compressions of 3T3, M450, M1 = doxy, M2 = control, M3 = doxy

# %%%% 21.01.21_M1 - juste une cellule pour commencer
dates = '21.01.21'
manips, wells, cells = 2, 1, [9]
depthoNames = '21.01.21_M1_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False, MatlabStyle = True)



# %%%% 21.01.21_M1
dates = '21.01.21'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.01.21_M1_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% 21.01.21_M2
dates = '21.01.21'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.01.21_M2_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% 21.01.21_M3
dates = '21.01.21'
manips, wells, cells = 3, 1, 'all'
depthoNames = '21.01.21_M3_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%% 21.01.18, compressions of 3T3, M450, M1 = control, M2 = doxy, M3 = control
# %%%% 21.01.18_M1
dates = '21.01.18'
manips, wells, cells = ['1-1', '1-2'], 1, 'all'
depthoNames = '21.01.18_M1_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% 21.01.18_M2
dates = '21.01.18'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.01.18_M2_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False, MatlabStyle = True)

# %%%% 21.01.18_M3
dates = '21.01.18'
manips, wells, cells = 3, 1, 'all'
depthoNames = '21.01.18_M3_M450_step20_100X'

PTL, timeSeries_DF, dfLogF = mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, 
                                  figureDir, timeSeriesDataDir,
                                  dates, manips, wells, cells, depthoNames, 
                                  expDf, methodT = 'max_entropy', factorT = 0.8, 
                                  redoAllSteps = False, MatlabStyle = True)



















