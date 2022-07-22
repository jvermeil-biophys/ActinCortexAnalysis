# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 21:27:55 2022

@author: JosephVermeil & AnumitaJawahar
"""

#!/usr/bin/env python
# coding: utf-8

# %% > Imports and constants

#### Main imports

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as st
import statsmodels.api as sm
import matplotlib.pyplot as plt


import os
import sys
import time
import random
import warnings
import itertools
import matplotlib

from copy import copy
from cycler import cycler
from datetime import date
from scipy.optimize import curve_fit
from matplotlib.gridspec import GridSpec
from pathlib import Path

# from statannot import add_stat_annotation
# pd.set_option('mode.chained_assignment',None)
# pd.set_option('display.max_columns', None)


#### Paths

COMPUTERNAME = os.environ['COMPUTERNAME']
if COMPUTERNAME == 'ORDI-JOSEPH':
    mainDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis"
    rawDir = "D://MagneticPincherData"
    ownCloudDir = "C://Users//JosephVermeil//ownCloud//ActinCortexAnalysis"
elif COMPUTERNAME == 'LARISA':
    mainDir = "C://Users//Joseph//Desktop//ActinCortexAnalysis"
    rawDir = "F:\JosephVermeil\MagneticPincherData"    
    ownCloudDir = "C://Users//Joseph//ownCloud//ActinCortexAnalysis"
elif COMPUTERNAME == 'DESKTOP-K9KOJR2':
    mainDir = "C://Users//anumi//OneDrive//Desktop//ActinCortexAnalysis"
    rawDir = "D:/Anumita/MagneticPincherData"    
elif COMPUTERNAME == 'DATA2JHODR':
    mainDir = "C://Users//BioMecaCell//Desktop//ActinCortexAnalysis"
    rawDir = "D:/Duya/MagneticPincherData"


experimentalDataDir = os.path.join(mainDir, "Data_Experimental_DB")
dataDir = os.path.join(mainDir, "Data_Analysis")
timeSeriesDataDir = os.path.join(dataDir, "TimeSeriesData")


figDir = os.path.join(dataDir, "Figures")
todayFigDir = os.path.join(figDir, "Historique//" + str(date.today()))


figDirLocal = os.path.join(rawDir, "Figures")
todayFigDirLocal = os.path.join(figDirLocal, "Historique//" + str(date.today()))

try:
    ownCloudFigDir = os.path.join(ownCloudDir, "Data_Analysis", "Figures")
    ownCloudTodayFigDir = os.path.join(ownCloudFigDir, "Historique//" + str(date.today()))
except:
    ownCloudFigDir, ownCloudTodayFigDir = '', ''
    
data_path=os.path.join(dataDir, "Global_MecaData_DB.csv")
file_path=os.path.join(rawDir, "Raw")
ratio_file_path=os.path.join(dataDir, "Global_ratio_file.csv")

#### Local imports
sys.path.append(mainDir + "//Code_Python")
# import PincherAnalysis_JV as jva
import PincherAnalysis_DB as dba
import utilityFunctions_JV as jvu

#### Potentially useful lines of code
# get_ipython().run_line_magic('load_ext', 'autoreload')
# get_ipython().run_line_magic('autoreload', '2')
# todayFigDirLocal

#### Pandas
pd.set_option('display.max_columns', None)
# pd.reset_option('display.max_columns')
pd.set_option('display.max_rows', None)
pd.reset_option('display.max_rows')


####  Matplotlib
matplotlib.rcParams.update({'figure.autolayout': True})
plt.style.use('default') #Dark layout

#### Fontsizes
SMALLER_SIZE = 10
SMALL_SIZE = 25
MEDIUM_SIZE = 16
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#### Bokeh
from bokeh.io import output_notebook, show
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Range1d
from bokeh.transform import factor_cmap
from bokeh.palettes import Category10
from bokeh.layouts import gridplot
output_notebook()

#### Markers
my_default_marker_list = ['o', 's', 'D', '>', '^', 'P', 'X', '<', 'v', 'p']
markerList10 = ['o', 's', 'D', '>', '^', 'P', 'X', '<', 'v', 'p']

#### Colors
# prop_cycle = plt.rcParams['axes.prop_cycle']
# colors = prop_cycle.by_key()['color']
my_default_color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                         '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
my_default_color_cycle = cycler(color=my_default_color_list)
plt.rcParams['axes.prop_cycle'] = my_default_color_cycle

pairedPalette = sns.color_palette("tab20")
pairedPalette = pairedPalette.as_hex()
pairedPalette

# clist = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a']
# sns.color_palette(clist)
colorList10 = my_default_color_list
sns.color_palette(my_default_color_list)


bigPalette1 = sns.color_palette("tab20b")
bigPalette1_hex = bigPalette1.as_hex()
bigPalette1

bigPalette2 = sns.color_palette("tab20c")
bigPalette2_hex = bigPalette2.as_hex()
bigPalette2

customPalette_hex = []
for ii in range(2, -1, -1):
    customPalette_hex.append(bigPalette2_hex[4*0 + ii]) # blue
    customPalette_hex.append(bigPalette2_hex[4*1 + ii]) # orange
    customPalette_hex.append(bigPalette2_hex[4*2 + ii]) # green
    customPalette_hex.append(bigPalette1_hex[4*3 + ii]) # red
    customPalette_hex.append(bigPalette2_hex[4*3 + ii]) # purple
    customPalette_hex.append(bigPalette1_hex[4*2 + ii]) # yellow-brown
    customPalette_hex.append(bigPalette1_hex[4*4 + ii]) # pink
    customPalette_hex.append(bigPalette1_hex[4*0 + ii]) # navy    
    customPalette_hex.append(bigPalette1_hex[4*1 + ii]) # yellow-green
    customPalette_hex.append(bigPalette2_hex[4*4 + ii]) # gray
    
# customPalette = sns.color_palette(customPalette_hex)
colorList30 = customPalette_hex

customPalette_hex = []
for ii in range(3, -1, -1):
    customPalette_hex.append(bigPalette2_hex[4*0 + ii]) # blue
    customPalette_hex.append(bigPalette2_hex[4*1 + ii]) # orange
    customPalette_hex.append(bigPalette2_hex[4*2 + ii]) # green
    customPalette_hex.append(bigPalette1_hex[4*3 + ii]) # red
    customPalette_hex.append(bigPalette2_hex[4*3 + ii]) # purple
    customPalette_hex.append(bigPalette1_hex[4*2 + ii]) # yellow-brown
    customPalette_hex.append(bigPalette1_hex[4*4 + ii]) # pink
    customPalette_hex.append(bigPalette1_hex[4*0 + ii]) # navy    
    customPalette_hex.append(bigPalette1_hex[4*1 + ii]) # yellow-green
    customPalette_hex.append(bigPalette2_hex[4*4 + ii]) # gray
    

# customPalette = sns.color_palette(customPalette_hex)
colorList40 = customPalette_hex

# TEST to get a darker list
# colorList40_darker = []
# for i in range(len(colorList40)):

#     c = colorList40[i]
#     print(c)
#     c_darker = jvu.lighten_color(c, 1.25)
#     colorList40_darker.append(c_darker)

# colorList40_darker = colorList40_darker.as_hex()
# %% Test the colors

# N = len(my_default_marker_list)
# X = np.arange(1, N+1)
# Y = np.arange(1, N+1)
# fig, ax = plt.subplots(1, 1, figsize = (3, 3))
# for i in range(N):
#     for j in range(N):
#         ax.plot([X[i]], [Y[-1-j]], color = my_default_color_list[i], marker = my_default_marker_list[j], 
#                 ls = '', markersize = 10, markeredgecolor = 'k')
# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# plt.show()

N = len(colorList40)
M = len(markerList10)
X = np.arange(1, N+1)
Y = np.arange(1, M+1)
fig, ax = plt.subplots(1, 1, figsize = (0.3*N, 0.3*M))
for i in range(N):
    for j in range(M):
        ax.plot([X[i]], [Y[-1-j]], color = colorList40[i], marker = markerList10[j], 
                ls = '', markersize = 10, markeredgecolor = 'k')
ax.set_xticks([])
ax.set_yticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.show()

# %% Reminders

# %%%  Stats
# 
# (A) How to compute confidence intervals of fitted parameters with (1-alpha) confidence:
# 
#     0) from scipy import stats
#     1) df = nb_pts - nb_parms ; se = diag(cov)**0.5
#     2) Student t coefficient : q = stat.t.ppf(1 - alpha / 2, df)
#     3) ConfInt = [params - q*se, params + q*se]

# ## 

# %% TimeSeries functions


# %%% List files
allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) \
                           if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
print(allTimeSeriesDataFiles)


# %%% Get a time series

df = dba.getCellTimeSeriesData('22-06-10_M1_P1_C1')


# # %%% Plot a time series

dba.plotCellTimeSeriesData('22-06-10_M1_P1_C8')


# %%% Plot multiple time series

allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
for f in allTimeSeriesDataFiles:
   if '22-07-06_M1' in f:
       dba.plotCellTimeSeriesData(jvu.findInfosInFileName(f, 'cellID'))

# # %%% Close all

# plt.close('all')


# %%% Plot multiple cells time on the first loop

allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
for f in allTimeSeriesDataFiles:
   if '22-07-06_M1' in f:
       dba.plotCellDifferentAdhesion(f[:-4])

# # %%% Close all

plt.close('all')


# %%% Experimental conditions

expDf = jvu.getExperimentalConditions(experimentalDataDir, save=True , sep = ';')

# =============================================================================
# %%% Constant Field

# %%%% Update the table

#dba.computeGlobalTable_ctField(task='updateExisting', fileName = '', save=False, source = 'Python')



# %%%% Refresh the whole table

# dba.computeGlobalTable_ctField(task = 'fromScratch', fileName = '', save = True, source = 'Python')


# # # %%%% Display

# df = dba.getGlobalTable_ctField().head()



# =============================================================================
# %%% Mechanics

# %%%% Update the table

dba.computeGlobalTable_meca(task = 'updateExisting', fileName = 'Global_MecaData_DB', 
                             save = True, PLOT = True, source = 'Python')

  
# %%%% Refresh the whole table

dba.computeGlobalTable_meca(task = 'fromScratch', fileName = 'Global_MecaData_DB', 
                            save = True, PLOT = True, source = 'Python')

# %%%% Specific experiments

# Task = '22-03-31_M9_P2_C2' # For instance '22-03-30 & '22-03-31'
# dba.computeGlobalTable_meca(task = Task, fileName = 'Global_MecaData_AJ', 
#                             save = True, PLOT = True, source = 'Python') # task = 'updateExisting'


# %%%% Precise dates (to plot)


date = '22-07-06 & 22-06-16 & 22-06-10 & 22-07-12'
dba.computeGlobalTable_meca(task = date, fileName = 'Global_MecaData_DB', 
                            save = True, PLOT = True, source = 'Python') # task = 'updateExisting'

# %%%% Display

df = dba.getGlobalTable_meca('Global_MecaData_Py2').tail()


# =============================================================================
# %%% Fluorescence

# %%%% Display

df = dba.getFluoData().head()




# #############################################################################
# %% > Data import & export

#### Data import

#### Display

# df1 = jvu.getExperimentalConditions().head()
# df2 = dba.getGlobalTable_ctField().head()
# df3 = dba.getGlobalTable_meca().head()
# df4 = dba.getFluoData().head()


#### GlobalTable_ctField

GlobalTable_ctField = dba.getGlobalTable(kind = 'ctField')
GlobalTable_ctField.head()


#### GlobalTable_ctField_Py

GlobalTable_ctField_Py = dba.getGlobalTable(kind = 'ctField_py')
GlobalTable_ctField_Py.head()


#### GlobalTable_meca

GlobalTable_meca = dba.getGlobalTable(kind = 'meca_matlab')
GlobalTable_meca.tail()


#### GlobalTable_meca_Py

GlobalTable_meca = dba.getGlobalTable(kind = 'Global_MecaData_DB')
GlobalTable_meca.head()


#### GlobalTable_meca_Py2

GlobalTable_meca_Py2 = dba.getGlobalTable(kind = 'meca_py2')
GlobalTable_meca_Py2.head()


#### Global_MecaData_NonLin_Py

# GlobalTable_meca_nonLin = dba.getGlobalTable(kind = 'meca_nonLin')
GlobalTable_meca_nonLin = dba.getGlobalTable(kind = 'Global_MecaData_NonLin2_Py')
GlobalTable_meca_nonLin.head()

        
# %% Plots
              
#  %%% ManipID
# manipID_to_plot='22-06-10_M1'
# x_axes=['ratio_adhesion', 'diameter', 'ratio_area','area_adhesion']

# for x_axis in x_axes:
#     dba.plot_bestH0(manipID_to_plot, x_axis)
#     dba.plot_KChadwick(manipID_to_plot, x_axis)

# %%%% All dates

# Review one more time
x_axes=['ratio_adhesion', 'diameter', 'ratio_area','area_adhesion']

for x_axis in x_axes:
    dba.plot_bestH0(x_axis)
    dba.plot_KChadwick(x_axis)



                   
        

                
                
                
            

