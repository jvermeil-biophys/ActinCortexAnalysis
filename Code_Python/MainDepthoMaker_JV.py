# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 11:20:48 2021

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

from scipy import interpolate
from scipy import signal
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

# 6. Others
SCALE_100X = 15.8 # pix/µm
SCALE_63X = 9.9 # pix/µm

# %% Import of the BeadTracker functions

from BeadTracker import *

# %% 

# %% 

# %% All depthos from 21.12.16 3T3 experiments

mainDirPath = 'D://MagneticPincherData//Raw'
date = '21.12.16'

subdir = 'M1_M450'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)


subdir = 'M2_M270'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M270'
saveLabel = date + '_M2_M270_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)





# %% All depthos from 21.12.08 3T3 experiments

mainDirPath = 'D://MagneticPincherData//Raw'
date = '21.12.08'

subdir = 'M1_M270'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M270'
saveLabel = date + '_M1_M270_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)



subdir = 'M2_M450'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)


 
# %% All depthos from 21.01.21 3T3 experiments

mainDirPath = 'D://MagneticPincherData//Raw'
date = '21.01.21'

subdir = 'Deptho_M1'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)



subdir = 'Deptho_M2'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)



subdir = 'Deptho_M3'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M3_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

# %% All depthos from 21.01.18 3T3 experiments

mainDirPath = 'D://MagneticPincherData//Raw'


date = '21.01.18'
subdir = 'Deptho_M1'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)



subdir = 'Deptho_M2'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)



subdir = 'Deptho_M3'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'EtalonnageZ')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M3_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)