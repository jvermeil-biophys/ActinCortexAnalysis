# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 21:27:55 2022

@author: JosephVermeil
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

# from statannot import add_stat_annotation
# pd.set_option('mode.chained_assignment',None)
# pd.set_option('display.max_columns', None)


#### Paths

COMPUTERNAME = os.environ['COMPUTERNAME']
if COMPUTERNAME == 'ORDI-JOSEPH':
    mainDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis"
    rawDir = "D://MagneticPincherData"
    ownCloudDir = "C://Users//JosephVermeil//ownCloud//ActinCortexAnalysis"
    experimentalDataDir = os.path.join(mainDir, "Data_Experimental_JV")
elif COMPUTERNAME == 'LARISA':
    mainDir = "C://Users//Joseph//Desktop//ActinCortexAnalysis"
    rawDir = "F:\JosephVermeil\MagneticPincherData"    
    ownCloudDir = "C://Users//Joseph//ownCloud//ActinCortexAnalysis"
    experimentalDataDir = os.path.join(mainDir, "Data_Experimental_JV")
elif COMPUTERNAME == 'DESKTOP-K9KOJR2':
    mainDir = "C://Users//anumi//OneDrive//Desktop//ActinCortexAnalysis"
    experimentalDataDir = os.path.join(mainDir, "Data_Experimental_AJ")
elif COMPUTERNAME == '':
    mainDir = "C://Users//josep//Desktop//ActinCortexAnalysis"
    ownCloudDir = "C://Users//josep//ownCloud//ActinCortexAnalysis"
    experimentalDataDir = os.path.join(mainDir, "Data_Experimental_JV")


# experimentalDataDir = os.path.join(mainDir, "Data_Experimental")
dataDir = os.path.join(mainDir, "Data_Analysis")
timeSeriesDataDir = os.path.join(dataDir, "TimeSeriesData")


figDir = os.path.join(dataDir, "Figures")
todayFigDir = os.path.join(figDir, "Historique//" + str(date.today()))


figDirLocal = os.path.join(rawDir, "Figures")
todayFigDirLocal = os.path.join(figDirLocal, "Historique//" + str(date.today()))


ownCloudFigDir = os.path.join(ownCloudDir, "Data_Analysis", "Figures")
ownCloudTodayFigDir = os.path.join(ownCloudFigDir, "Historique//" + str(date.today()))

#### Local imports
sys.path.append(mainDir + "//Code_Python")
import PincherAnalysis_JV as jva
import MechanicsAnalysis_AJ as aja
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

#### Fontsizes
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

df = jva.getCellTimeSeriesData('22-02-09_M1_P1_C7')


# %%% Plot a time series

# jva.plotCellTimeSeriesData('21-02-10_M1_P1_C2')
jva.plotCellTimeSeriesData('22-03-21_M3_P1_C1_sin5-3_1Hz')

# %%% Plot a time series

jva.plotCellTimeSeriesData('22-03-21_M3_P1_C1_sin5-3_2Hz')

# %%% Variation on the plotCellTimeSeriesData function

# cellID = '22-02-09_M1_P1_C2'
# fromPython = True

# X = 'T'
# Y = np.array(['B', 'F'])
# units = np.array([' (mT)', ' (pN)'])
# timeSeriesDataFrame = jva.getCellTimeSeriesData(cellID, fromPython)
# print(timeSeriesDataFrame.shape)
# if not timeSeriesDataFrame.size == 0:
# #         plt.tight_layout()
# #         fig.show() # figsize=(20,20)
#     axes = timeSeriesDataFrame.plot(x=X, y=Y, kind='line', ax=None, subplots=True, sharex=True, sharey=False, layout=None, \
#                     figsize=(4,6), use_index=True, title = cellID + ' - f(t)', grid=None, legend=False, style=None, logx=False, logy=False, \
#                     loglog=False, xticks=None, yticks=None, xlim=None, ylim=None, rot=None, fontsize=None, colormap=None, \
#                     table=False, yerr=None, xerr=None, secondary_y=False, sort_columns=False)
#     plt.gcf().tight_layout()
#     for i in range(len(Y)):
#         axes[i].set_ylabel(Y[i] + units[i])
        
#     axes[0].set_yticks([0,2,4,6,8,10])
#     axes[0].grid(axis="y")
#     axes[1].set_yticks([0,10,20,30,40,50,100,150])
#     axes[1].grid(axis="y")
#     plt.gcf().show()
# else:
#     print('cell not found')
    
#### Nice display for B and F

# cellID = '22-01-12_M2_P1_C3'
# cellID = '22-01-12_M1_P1_C6'
cellID = '22-02-09_M1_P1_C9'
fromPython = True

X = 'T'
Y = np.array(['B', 'F'])
units = np.array([' (mT)', ' (pN)'])
tSDF = jva.getCellTimeSeriesData(cellID, fromPython)

defaultColorCycle = plt.rcParams['axes.prop_cycle']
customColorCycle = cycler(color=['purple', 'red'])
plt.rcParams['axes.prop_cycle'] = customColorCycle

if not tSDF.size == 0:
#         plt.tight_layout()
#         fig.show() # figsize=(20,20)
    axes = tSDF[tSDF['T']<=40].plot(x=X, y=Y, kind='line', 
                    subplots=True, sharex=True, sharey=False,
                    figsize=(4,4), use_index=True)
    
    plt.gcf().tight_layout()
    for i in range(len(Y)):
        axes[i].set_ylabel(Y[i] + units[i])
        
    axes[0].set_yticks([0,10,20,30,40,50])
    axes[0].tick_params(axis='y', labelrotation = 0, labelsize = 10)
    axes[0].grid(axis="y")
    axes[0].legend().set_visible(False)
    axes[1].set_yticks([k for k in range(0,1400,200)])
    axes[1].tick_params(axis='y', labelrotation = 0, labelsize = 10)
    axes[1].grid(axis="y")
    axes[1].legend().set_visible(False)
    axes[1].set_xticks([0,10,20,30,40])
    axes[1].tick_params(axis='x', labelrotation = 0, labelsize = 10)
    plt.gcf().show()
else:
    print('cell not found')
plt.rcParams['axes.prop_cycle'] = defaultColorCycle


# %%% Variation on the plotCellTimeSeriesData -> Sinus

# cellID = '22-01-12_M2_P1_C3'
# cellID = '22-01-12_M1_P1_C6'
cellID = '22-03-21_M3_P1_C1_sin5-3_1Hz'
fromPython = True

X = 'T'
Y = np.array(['B', 'F'])
units = np.array([' (mT)', ' (pN)'])
tSDF = jva.getCellTimeSeriesData(cellID, fromPython)

# defaultColorCycle = plt.rcParams['axes.prop_cycle']
# customColorCycle = cycler(color=['green', 'red'])
# plt.rcParams['axes.prop_cycle'] = customColorCycle

fig, ax = plt.subplots(1,1, figsize = (10,5))

if not tSDF.size == 0:
#         plt.tight_layout()
#         fig.show() # figsize=(20,20)
    tsDFplot = tSDF[(tSDF['T']>=3.04) & (tSDF['T']<=6.05)]
    
    ax.plot(tsDFplot['T'], 1000*(tsDFplot['D3']-4.503), color = colorList40[30], label = 'Thickness')
    ax.set_ylabel('Thickness (nm)', color = colorList40[30])
    ax.set_ylim([0, 350])
    ax.set_xlabel('Time (s)')

    axR = ax.twinx()
    axR.plot(tsDFplot['T'], tsDFplot['F'], color = colorList40[23], label = 'Force')
    axR.set_ylabel('Force (pN)', color = colorList40[23])
    axR.set_ylim([0, 150])

    # axes = tSDF[(tSDF['T']>=20) & (tSDF['T']<=60)].plot(x=X, y=Y, kind='line', 
    #                 subplots=True, sharex=True, sharey=False,
    #                 figsize=(4,4), use_index=True)
    
    # plt.gcf().tight_layout()
    # for i in range(len(Y)):
    #     axes[i].set_ylabel(Y[i] + units[i])
        
    # axes[0].set_yticks([0,10,20,30,40,50])
    # axes[0].tick_params(axis='y', labelrotation = 0, labelsize = 10)
    # axes[0].grid(axis="y")
    # axes[0].legend().set_visible(False)
    # axes[1].set_yticks([k for k in range(0,1400,200)])
    # axes[1].tick_params(axis='y', labelrotation = 0, labelsize = 10)
    # axes[1].grid(axis="y")
    # axes[1].legend().set_visible(False)
    # axes[1].set_xticks([0,10,20,30,40])
    # axes[1].tick_params(axis='x', labelrotation = 0, labelsize = 10)
    plt.gcf().show()
else:
    print('cell not found')
# plt.rcParams['axes.prop_cycle'] = defaultColorCycle

# %%% Variation on the plotCellTimeSeriesData -> Broken Ramp

# cellID = '22-01-12_M2_P1_C3'
# cellID = '22-01-12_M1_P1_C6'
cellID = '22-03-21_M4_P1_C3'
fromPython = True

X = 'T'
Y = np.array(['B', 'F'])
units = np.array([' (mT)', ' (pN)'])
tSDF = jva.getCellTimeSeriesData(cellID, fromPython)

# defaultColorCycle = plt.rcParams['axes.prop_cycle']
# customColorCycle = cycler(color=['green', 'red'])
# plt.rcParams['axes.prop_cycle'] = customColorCycle

fig, ax = plt.subplots(1,1, figsize = (6,3))

if not tSDF.size == 0:
#         plt.tight_layout()
#         fig.show() # figsize=(20,20)
    tsDFplot = tSDF[(tSDF['T']>=31) & (tSDF['T']<=34)]
    
    ax.plot(tsDFplot['T'], 1000*(tsDFplot['D3']-4.503), color = colorList40[30], label = 'Thickness')
    ax.set_ylabel('Thickness (nm)', color = colorList40[30])
    ax.set_ylim([0, 1200])
    ax.set_xlabel('Time (s)')

    axR = ax.twinx()
    axR.plot(tsDFplot['T'], tsDFplot['F'], color = colorList40[23], label = 'Force')
    axR.set_ylabel('Force (pN)', color = colorList40[23])
    axR.set_ylim([0, 1500])

    # axes = tSDF[(tSDF['T']>=20) & (tSDF['T']<=60)].plot(x=X, y=Y, kind='line', 
    #                 subplots=True, sharex=True, sharey=False,
    #                 figsize=(4,4), use_index=True)
    
    # plt.gcf().tight_layout()
    # for i in range(len(Y)):
    #     axes[i].set_ylabel(Y[i] + units[i])
        
    # axes[0].set_yticks([0,10,20,30,40,50])
    # axes[0].tick_params(axis='y', labelrotation = 0, labelsize = 10)
    # axes[0].grid(axis="y")
    # axes[0].legend().set_visible(False)
    # axes[1].set_yticks([k for k in range(0,1400,200)])
    # axes[1].tick_params(axis='y', labelrotation = 0, labelsize = 10)
    # axes[1].grid(axis="y")
    # axes[1].legend().set_visible(False)
    # axes[1].set_xticks([0,10,20,30,40])
    # axes[1].tick_params(axis='x', labelrotation = 0, labelsize = 10)
    plt.gcf().show()
else:
    print('cell not found')
# plt.rcParams['axes.prop_cycle'] = defaultColorCycle

# %%% Plot multiple time series

allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
for f in allTimeSeriesDataFiles:
    if '22-02-09_M3' in f:
        jva.plotCellTimeSeriesData(f[:-4])

allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
for f in allTimeSeriesDataFiles:
    if '22-02-09_M2' in f:
        jva.plotCellTimeSeriesData(f[:-4])


# %%% Close all

plt.close('all')


# %%% Functions acting on the trajectories

listeTraj = jva.getCellTrajData('21-12-16_M1_P1_C10', Ntraj = 2)
listeTraj[1]['pos']

def testTraj(D0):
    trajDir = os.path.join(timeSeriesDataDir, 'Trajectories')
    allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) 
                              if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) 
                                  and f.endswith(".csv"))]
    cellIDList = []
    for f in allTimeSeriesDataFiles:
        cellID = jvu.findInfosInFileName(f, 'cellID')
        if '21-12-08' in cellID:
            cellIDList.append(cellID)

    print(cellIDList)

    fig, ax = plt.subplots(1,2, figsize = (8,4))
    width = 15
    ax[0].set_title('Inside bead')
    ax[0].axis([-width,width,-width,width])
    ax[1].set_title('Outside bead')
    ax[1].axis([-width,width,-width,width])

    for C in (cellIDList):
        listeTraj = jva.getCellTrajData(C)
        iOut = (listeTraj[1]['pos'] == 'out')
        iIn = 1 - iOut
        dfIn, dfOut = listeTraj[iIn]['df'], listeTraj[iOut]['df']
        Xin = dfIn.X.values - dfIn.X.values[0]
        Yin = dfIn.Y.values - dfIn.Y.values[0]
        Xout = dfOut.X.values - dfOut.X.values[0]
        Yout = dfOut.Y.values - dfOut.Y.values[0]
        npts = len(Xin)
        for iii in range(1, npts):
            D2in = ((Xin[iii]-Xin[iii-1])**2 + (Yin[iii]-Yin[iii-1])**2) ** 0.5
            D2out = ((Xout[iii]-Xout[iii-1])**2 + (Yout[iii]-Yout[iii-1])**2) ** 0.5
            if D2in > D0 and D2out > D0:
                dxCorrIn = Xin[iii-1]-Xin[iii]
                dyCorrIn = Yin[iii-1]-Yin[iii]
                dxCorrOut = Xout[iii-1]-Xout[iii]
                dyCorrOut = Yout[iii-1]-Yout[iii]
                Xin[iii:] += dxCorrIn
                Yin[iii:] += dyCorrIn
                Xout[iii:] += dxCorrOut
                Yout[iii:] += dyCorrOut
        if np.max(np.abs(Xin)) < width and np.max(np.abs(Yin)) < width and np.max(np.abs(Xout)) < width and np.max(np.abs(Yout)) < width:
            ax[0].plot(Xin, Yin)
            ax[1].plot(Xout, Yout)

    plt.show()
    
testTraj(1)
    
testTraj(1.5)




# #############################################################################
# %% GlobalTables functions



# %%% Experimental conditions

jvu.getExperimentalConditions(experimentalDataDir, save=True , sep = ';')





# =============================================================================
# %%% Constant Field

# %%%% Update the table

# jva.computeGlobalTable_ctField(task='updateExisting', save=False)


# %%%% Refresh the whole table

# jva.computeGlobalTable_ctField(task = 'updateExisting', fileName = 'Global_CtFieldData_Py', save = True, source = 'Python') # task = 'updateExisting'


# %%%% Display

df = jva.getGlobalTable_ctField().head()






# =============================================================================
# %%% Mechanics

# %%%% Update the table

# jva.computeGlobalTable_meca(task = 'updateExisting', fileName = 'Global_MecaData_Py', 
#                             save = False, PLOT = False, source = 'Matlab') # task = 'updateExisting'


# %%%% Refresh the whole table

# jva.computeGlobalTable_meca(task = 'updateExisting', fileName = 'Global_MecaData_Py2', 
#                             save = True, PLOT = False, source = 'Python') # task = 'updateExisting'

# %%%% Drugs

drugTask = '22-03-30'
# jva.computeGlobalTable_meca(task = drugTask, fileName = 'Global_MecaData_Drugs_Py', 
#                             save = False, PLOT = True, source = 'Python') # task = 'updateExisting'


# %%%% Non-Lin

nonLinTask = '21-12-08 & 22-01-12 & 22-02-09'
# jva.computeGlobalTable_meca(task = nonLinTask, fileName = 'Global_MecaData_NonLin_Py', 
#                             save = False, PLOT = False, source = 'Python') # task = 'updateExisting'

# %%%% MCA

MCAtask = '21-01-18 & 21-01-21'
# jva.computeGlobalTable_meca(task = MCAtask, fileName = 'Global_MecaData_MCA', 
#                             save = False, PLOT = False, source = 'Python') # task = 'updateExisting'
# jva.computeGlobalTable_meca(task = MCAtask, fileName = 'Global_MecaData_MCA2', 
#                             save = True, PLOT = False, source = 'Python') # task = 'updateExisting'


# %%%% HoxB8

HoxB8task = '22-05-03' #' & 22-05-04 & 22-05-05'
jva.computeGlobalTable_meca(task = HoxB8task, fileName = 'Global_MecaData_HoxB8', 
                            save = True, PLOT = True, source = 'Python') # task = 'updateExisting'
# jva.computeGlobalTable_meca(task = MCAtask, fileName = 'Global_MecaData_MCA2', 
#                             save = True, PLOT = False, source = 'Python') # task = 'updateExisting'

# %%%% Demo for Duya

Demo = '22-06-16' #' & 22-05-04 & 22-05-05'
jva.computeGlobalTable_meca(task = Demo, fileName = 'Global_MecaData_Demo', 
                            save = True, PLOT = True, source = 'Python') # task = 'updateExisting'
# jva.computeGlobalTable_meca(task = MCAtask, fileName = 'Global_MecaData_MCA2', 
#                             save = True, PLOT = False, source = 'Python') # task = 'updateExisting'


# %%%% Precise dates (to plot)

# jva.computeGlobalTable_meca(task = '22-02-09', fileName = 'Global_MecaData_Py2', save = False, PLOT = True, source = 'Python') # task = 'updateExisting'
# jva.computeGlobalTable_meca(task = '21-01-18', fileName = 'Global_MecaData_Py2', save = False, PLOT = True, source = 'Python') # task = 'updateExisting'
# jva.computeGlobalTable_meca(task = '21-01-21', fileName = 'Global_MecaData_Py2', save = False, PLOT = True, source = 'Python') # task = 'updateExisting'
# jva.computeGlobalTable_meca(task = '22-02-09_M1', fileName = 'Global_MecaData_NonLin2_Py', 
#                             save = False, PLOT = True, source = 'Python') # task = 'updateExisting'
# jva.computeGlobalTable_meca(task = '21-01-18_M2_P1_C3', fileName = 'Global_MecaData_NonLin2_Py', 
#                             save = False, PLOT = True, source = 'Python') # task = 'updateExisting'
jva.computeGlobalTable_meca(task = '22-02-09_M1_P1_C9', fileName = 'aaa', 
                            save = False, PLOT = True, source = 'Python') # task = 'updateExisting'


# %%%% Display

df = jva.getGlobalTable_meca('Global_MecaData_Py2').tail()




# =============================================================================
# %%% Fluorescence

# %%%% Display

df = jva.getFluoData().head()




# #############################################################################
# %% > Data import & export

#### Data import

#### Display

# df1 = jvu.getExperimentalConditions().head()
# df2 = jva.getGlobalTable_ctField().head()
# df3 = jva.getGlobalTable_meca().head()
# df4 = jva.getFluoData().head()


#### GlobalTable_ctField

GlobalTable_ctField = jva.getGlobalTable(kind = 'ctField')
GlobalTable_ctField.head()


#### GlobalTable_ctField_Py

GlobalTable_ctField_Py = jva.getGlobalTable(kind = 'ctField_py')
GlobalTable_ctField_Py.head()


#### GlobalTable_meca

GlobalTable_meca = jva.getGlobalTable(kind = 'meca_matlab')
GlobalTable_meca.tail()


#### GlobalTable_meca_Py

GlobalTable_meca_Py = jva.getGlobalTable(kind = 'meca_py')
GlobalTable_meca_Py.tail()


#### GlobalTable_meca_Py2

GlobalTable_meca_Py2 = jva.getGlobalTable(kind = 'meca_py2')
GlobalTable_meca_Py2.head()


#### Global_MecaData_NonLin_Py

# GlobalTable_meca_nonLin = jva.getGlobalTable(kind = 'meca_nonLin')
GlobalTable_meca_nonLin = jva.getGlobalTable(kind = 'Global_MecaData_NonLin2_Py')
GlobalTable_meca_nonLin.head()


#### Global_MecaData_MCA

# GlobalTable_meca_MCA = jva.getGlobalTable(kind = 'meca_MCA')
GlobalTable_meca_MCA = jva.getGlobalTable(kind = 'Global_MecaData_MCA2')
GlobalTable_meca_MCA.head()

#### Global_MecaData_MCA

# GlobalTable_meca_MCA = jva.getGlobalTable(kind = 'meca_MCA')
GlobalTable_meca_MCA = jva.getGlobalTable(kind = 'Global_MecaData_MCA2')
GlobalTable_meca_MCA.head()


# %%% Custom data export

# %%%% 21-06-28 Export of E - h data for Julien

# GlobalTable_meca_Py_expJ = GlobalTable_meca_Py.loc[GlobalTable_meca_Py['validatedFit'] & GlobalTable_meca_Py['validatedThickness']]\
#                                                [['cell type', 'cell subtype', 'bead type', 'drug', 'substrate', 'compNum', \
#                                                  'EChadwick', 'H0Chadwick', 'surroundingThickness', 'ctFieldThickness']]
# GlobalTable_meca_Py_expJ = GlobalTable_meca_Py_expJ.reset_index()
# GlobalTable_meca_Py_expJ = GlobalTable_meca_Py_expJ.drop('index', axis=1)
# savePath = os.path.join(dataDir, 'mecanicsData_3T3.csv')
# GlobalTable_meca_Py_expJ.to_csv(savePath, sep=';')


# %% > Plotting Functions

# %%% Objects declaration


Styles = {''} # Project of automatic formatting according to the type of data

renameDict1 = {'SurroundingThickness':'Thickness (nm) [b&a]',
               'surroundingThickness':'Thickness (nm) [b&a]',
               'ctFieldThickness':'Thickness at low force (nm)',
               'EChadwick': 'E Chadwick (Pa)',
               'medianThickness': 'Median Thickness (nm)',               
               'fluctuAmpli': 'Fluctuations Amplitude (nm)',               
               'meanFluoPeakAmplitude' : 'Fluo Intensity (a.u.)',               
               'none':'control',               
               'doxycyclin':'expressing iMC linker',               
               'none & BSA coated glass':'control & non adherent',               
               'doxycyclin & BSA coated glass':'iMC & non adherent',               
               'none & 20um fibronectin discs':'control & adherent on fibro',               
               'doxycyclin & 20um fibronectin discs':'iMC & adherent on fibro',               
               'BSA coated glass & none':'control & non adherent',               
               'BSA coated glass & doxycyclin':'iMC & non adherent',               
               '20um fibronectin discs & none':'control & adherent on fibro',               
               '20um fibronectin discs & doxycyclin':'iMC & adherent on fibro',               
               'aSFL & none':'aSFL control',
               'aSFL & doxycyclin':'aSFL iMC',
               'aSFL-6FP & none':'aSFL-6FP control',               
               'aSFL-6FP & doxycyclin':'aSFL-6FP long-iMC',               
               'aSFL-6FP-2 & none':'aSFL-6FP-2 control',               
               'aSFL-6FP-2 & doxycyclin':'aSFL-6FP-2 long-iMC',               
               'aSFL-A8 & none':'aSFL-A8 control',               
               'aSFL-A8 & doxycyclin':'aSFL-A8 iMC',               
               'aSFL-A8-2 & none':'aSFL-A8-2 control',               
               'aSFL-A8-2 & doxycyclin':'aSFL-A8-2 iMC',               
               'dmso' : 'DMSO, no linker', 
               'smifh2' : 'SMIFH2, no linker', 
               'dmso, doxycyclin' : 'DMSO, iMC linker', 
               'smifh2, doxycyclin' : 'SMIFH2, iMC linker'}

styleDict1 =  {'none & BSA coated glass':{'color':'#ff9896','marker':'^'},               
               'doxycyclin & BSA coated glass':{'color':'#d62728','marker':'^'},               
               'none & 20um fibronectin discs':{'color':'#aec7e8','marker':'o'},               
               'doxycyclin & 20um fibronectin discs':{'color':'#1f77b4','marker':'o'},               
               'none':{'color':'#aec7e8','marker':'o'},               
               'doxycyclin':{'color':'#1f77b4','marker':'o'},               
               'BSA coated glass & none':{'color':'#ff9896','marker':'^'},               
               'BSA coated glass & doxycyclin':{'color':'#d62728','marker':'^'},               
               '20um fibronectin discs & none':{'color':'#aec7e8','marker':'o'},               
               '20um fibronectin discs & doxycyclin':{'color':'#1f77b4','marker':'o'},               
               'aSFL':{'color':'colorList40[10]','marker':'o'},               
               'aSFL-6FP':{'color':'#2ca02c','marker':'o'},               
               'aSFL-A8':{'color':'#ff7f0e','marker':'o'},                
               'aSFL & none':{'color':'colorList40[10]','marker':'o'},               
               'aSFL & doxycyclin':{'color':'colorList40[30]','marker':'o'},               
               'aSFL-6FP & none':{'color':'#98df8a','marker':'o'},               
               'aSFL-6FP & doxycyclin':{'color':'#2ca02c','marker':'o'},               
               'aSFL-A8 & none':{'color':'#ffbb78','marker':'o'},              
               'aSFL-A8 & doxycyclin':{'color':'#ff7f0e','marker':'o'},
               'DictyDB_M270':{'color':'lightskyblue','marker':'o'}, 
               'DictyDB_M450':{'color': 'maroon','marker':'o'},
               'M270':{'color':'lightskyblue','marker':'o'}, 
               'M450':{'color': 'maroon','marker':'o'}}


# These functions use matplotlib.pyplot and seaborn libraries to display 1D categorical or 2D plots

# %%% Main functions


def D1Plot(data, fig = None, ax = None, CondCol=[], Parameters=[], Filters=[], 
           Boxplot=True, AvgPerCell=False, cellID='cellID', co_order=[],
           stats=True, statMethod='Mann-Whitney', box_pairs=[], 
           figSizeFactor = 1, markersizeFactor=1, orientation = 'h', useHue = False, 
           stressBoxPlot = False, bypassLog = False, autoscale = True):
    
    data_filtered = data
    for fltr in Filters:
        data_filtered = data_filtered.loc[fltr]    
    NCond = len(CondCol)
    
    print(data_filtered.shape)
    
    if NCond == 1:
        CondCol = CondCol[0]
        
    elif NCond > 1:
        newColName = ''
        for i in range(NCond):
            newColName += CondCol[i]
            newColName += ' & '
        newColName = newColName[:-3]
        data_filtered[newColName] = ''
        for i in range(NCond):
            data_filtered[newColName] += data_filtered[CondCol[i]].astype(str)
            data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x + ' & ')
        data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x[:-3])
        CondCol = newColName
    
    if AvgPerCell:
        group = data_filtered.groupby(cellID)
        dictAggMean = getDictAggMean(data_filtered)
#         dictAggMean['EChadwick'] = 'median'
        data_filtered = group.agg(dictAggMean)
        
    data_filtered.sort_values(CondCol, axis=0, ascending=True, inplace=True)
    
    if len(co_order) > 0:
        p = getStyleLists_Sns(co_order, styleDict1)
    else:
        p = sns.color_palette()
    
    NPlots = len(Parameters)
    Conditions = list(data_filtered[CondCol].unique())
    if len(co_order) == 0:
        co_order = Conditions
    
    if fig == None:
        if orientation == 'h':
            fig, ax = plt.subplots(1, NPlots, figsize = (5*NPlots*NCond*figSizeFactor, 5))
        elif orientation == 'v':
            fig, ax = plt.subplots(NPlots, 1, figsize = (5*NCond*figSizeFactor, 5*NPlots))
    else:
        pass
    
    scaley = autoscale
    
    markersize = 5*markersizeFactor
    
    if NPlots == 1:
        ax = np.array([ax])
    
    for k in range(NPlots):

        if not bypassLog and (('EChadwick' in Parameters[k]) or ('KChadwick' in Parameters[k])):
            ax[k].set_yscale('log')


        if Boxplot:
            if stressBoxPlot:
                sns.boxplot(x=CondCol, y=Parameters[k], data=data_filtered, ax=ax[k], 
                            width = 0.5, showfliers = False, order= co_order, 
                            medianprops={"color": 'darkred', "linewidth": 2, 'alpha' : 0.8, 'zorder' : 4},
                            boxprops={"facecolor": 'None', "edgecolor": 'k',"linewidth": 2, 'alpha' : 0.7, 'zorder' : 4},
        #                   boxprops={"color": color, "linewidth": 0.5},
                            whiskerprops={"color": 'k', "linewidth": 2, 'alpha' : 0.7, 'zorder' : 4},
                            capprops={"color": 'k', "linewidth": 2, 'alpha' : 0.7, 'zorder' : 4})
                            # scaley = scaley)
            else:
                sns.boxplot(x=CondCol, y=Parameters[k], data=data_filtered, ax=ax[k], 
                            width = 0.5, showfliers = False, order= co_order, 
                            medianprops={"color": 'darkred', "linewidth": 1.5, 'alpha' : 0.8, 'zorder' : 1},
                            boxprops={"facecolor": 'None', "edgecolor": 'k',"linewidth": 1, 'alpha' : 0.7, 'zorder' : 1},
        #                   boxprops={"color": color, "linewidth": 0.5},
                            whiskerprops={"color": 'k', "linewidth": 1, 'alpha' : 0.7, 'zorder' : 1},
                            capprops={"color": 'k', "linewidth": 1, 'alpha' : 0.7, 'zorder' : 1})
                            # scaley = scaley)
            
            # data_filtered.boxplot(column=Parameters[k], by = CondCol, ax=ax[k],showfliers = False) # linewidth = 2, width = 0.5

        if stats:
            if len(box_pairs) == 0:
                box_pairs = makeBoxPairs(co_order)
            addStat_df(ax[k], data_filtered, box_pairs, Parameters[k], CondCol, test = statMethod)
            # add_stat_annotation(ax[k], x=CondCol, y=Parameters[k], data=data_filtered,box_pairs = box_pairs,test=statMethod, text_format='star',loc='inside', verbose=2)
        
        if not useHue:
            sns.swarmplot(x=CondCol, y=Parameters[k], data=data_filtered, ax=ax[k], order = co_order,
                          size=markersize, edgecolor='k', linewidth = 1*markersizeFactor, palette = p)
        else:
            sns.swarmplot(x=CondCol, y=Parameters[k], data=data_filtered, ax=ax[k], order = co_order,
                          size=markersize, edgecolor='k', linewidth = 1*markersizeFactor, palette = p, 
                          hue = 'manipID')
            legend = ax[k].legend()
            legend.remove()

        ax[k].set_xlabel('')
        ax[k].set_ylabel(Parameters[k])
        ax[k].tick_params(axis='x', labelrotation = 10)
        ax[k].yaxis.grid(True)           
    
    plt.rcParams['axes.prop_cycle'] = my_default_color_cycle
    return(fig, ax)



def D1PlotDetailed(data, CondCol=[], Parameters=[], Filters=[], Boxplot=True, cellID='cellID', 
                   co_order=[], stats=True, statMethod='Mann-Whitney', box_pairs=[],
                   figSizeFactor = 1, markersizeFactor=1, orientation = 'h', showManips = True):
    
    data_f = data
    for fltr in Filters:
        data_f = data_f.loc[fltr]    
    NCond = len(CondCol)
    
    print(min(data_f['EChadwick'].values))
    
    if NCond == 1:
        CondCol = CondCol[0]
        
    elif NCond > 1:
        newColName = ''
        for i in range(NCond):
            newColName += CondCol[i]
            newColName += ' & '
        newColName = newColName[:-3]
        data_f[newColName] = ''
        for i in range(NCond):
            data_f[newColName] += data_f[CondCol[i]].astype(str)
            data_f[newColName] = data_f[newColName].apply(lambda x : x + ' & ')
        data_f[newColName] = data_f[newColName].apply(lambda x : x[:-3])
        CondCol = newColName
        
    data_f_agg = getAggDf(data_f, cellID, CondCol, Parameters)
    
    NPlots = len(Parameters)
    Conditions = list(data_f[CondCol].unique())
        
    if orientation == 'h':
        fig, ax = plt.subplots(1, NPlots, figsize = (5*NPlots*NCond*figSizeFactor, 5))
    elif orientation == 'v':
        fig, ax = plt.subplots(NPlots, 1, figsize = (5*NCond*figSizeFactor, 5*NPlots))
    
    markersize = 5*markersizeFactor
    
    if NPlots == 1:
        ax = np.array([ax])
    
    # Colors and markers
    if len(co_order) > 0:
        Conditions = co_order
        colorList, mL = getStyleLists(co_order, styleDict1)
    else:
        co_order = Conditions
        colorList = colorList10
    markerList = markerList10
        
        
    for k in range(NPlots):
        
        Parm = Parameters[k]
        if 'EChadwick' in Parm or 'KChadwick' in Parm:
            ax[k].set_yscale('log')

        for i in range(len(Conditions)):
            c = Conditions[i]
            sub_data_f_agg = data_f_agg.loc[data_f_agg[CondCol] == c]
            Ncells = sub_data_f_agg.shape[0]
            
            color = colorList[i]
            
            if showManips:
                allManipID = list(sub_data_f_agg['manipID'].unique())
                alreadyLabeledManip = []
                dictMarker = {}
                for mi in range(len(allManipID)):
                    dictMarker[allManipID[mi]] = markerList[mi]

                
            midPos = i
            startPos = i-0.4
            stopPos = i+0.4
            
            values = sub_data_f_agg[Parm + '_mean'].values
            errors = sub_data_f_agg[Parm + '_std'].values
            
            if Boxplot:
                ax[k].boxplot(values, positions=[midPos], widths=0.4, patch_artist=True,
                            showmeans=False, showfliers=False,
                            medianprops={"color": 'k', "linewidth": 1.5, 'alpha' : 0.5},
                            boxprops={"facecolor": color, "edgecolor": 'k',
                                      "linewidth": 1, 'alpha' : 0.5},
    #                         boxprops={"color": color, "linewidth": 0.5},
                            whiskerprops={"color": 'k', "linewidth": 1.5, 'alpha' : 0.5},
                            capprops={"color": 'k', "linewidth": 1.5, 'alpha' : 0.5},
                             zorder = 1)
                

            step = 0.8/(Ncells-1)
            
            for m in range(Ncells):
                
                thisCellID = sub_data_f_agg.index.values[m]
                thisCell_data = data_f.loc[data_f[cellID] == thisCellID]
                thisCell_allValues = thisCell_data[Parm].values
                nval = len(thisCell_allValues)
                
                if showManips:
                    thisManipID = sub_data_f_agg['manipID'].values[m]
                    marker = dictMarker[thisManipID]
                else:
                    marker = 'o'
                
                cellpos = startPos + m*step
                ax[k].errorbar([cellpos], [values[m]], color = color, marker = marker, markeredgecolor = 'k', 
                               yerr=errors[m], capsize = 2, ecolor = 'k', elinewidth=0.8, barsabove = False)
                
                if showManips and thisManipID not in alreadyLabeledManip:
                    alreadyLabeledManip.append(thisManipID)
                    textLabel = thisManipID
                    ax[k].plot([cellpos], [values[m]], 
                               color = color, 
                               marker = marker, markeredgecolor = 'k', markersize = markersize,
                               label = textLabel)
                    
                ax[k].plot(cellpos*np.ones(nval), thisCell_allValues, 
                           color = color, marker = '_', ls = '', markersize = markersize)

        ax[k].set_xlabel('')
        ax[k].set_xticklabels(labels = Conditions)
        ax[k].tick_params(axis='x', labelrotation = 10)
        
        ax[k].set_ylabel(Parameters[k])
        ax[k].yaxis.grid(True)
        
        ax[k].legend(fontsize = 6)
        
        if stats:
            renameDict = {Parameters[k] + '_mean' : Parameters[k]}
            if len(box_pairs) == 0:
                box_pairs = makeBoxPairs(co_order)
            addStat_df(ax[k], data_f_agg.rename(columns = renameDict), 
                    box_pairs, Parameters[k], CondCol, test = statMethod)
        
    plt.rcParams['axes.prop_cycle'] = my_default_color_cycle
    return(fig, ax)




def D1PlotPaired(data, Parameters=[], Filters=[], Boxplot=True, cellID='cellID', 
                   co_order=[], stats=True, statMethod='Wilcox_greater', box_pairs=[],
                   figSizeFactor = 1, markersizeFactor=1, orientation = 'h', labels = []):
    
    data_f = data
    for fltr in Filters:
        data_f = data_f.loc[fltr]
    
    if orientation == 'h':
        fig, ax = plt.subplots(1, 1, figsize = (5*figSizeFactor, 5))
    elif orientation == 'v':
        fig, ax = plt.subplots(1, 1, figsize = (5*figSizeFactor, 5))
        
    markersize = 5*markersizeFactor
        
    if 'EChadwick' in Parameters[0]:
        ax.set_yscale('log')
    
    pos = 1
    start = pos - 0.4
    stop = pos + 0.4
    NParms = len(Parameters)
    width = (stop-start)/NParms
    
    if len(labels) == len(Parameters):
        tickLabels = labels
    else:
        tickLabels = Parameters
    
    posParms = [start + 0.5*width + i*width for i in range(NParms)]
    parmsValues = np.array([data_f[P].values for P in Parameters])
#     for P in Parameters:
#         print(len(data_f[P].values))
#         print(data_f[P].values)
#     print(parmsValues)
    NValues = parmsValues.shape[1]
    print('Number of values : {:.0f}'.format(NValues))
    
    for i in range(NParms):
        color = my_default_color_list[i]
        
        ax.plot([posParms[i] for k in range(NValues)], parmsValues[i,:], 
                color = color, marker = 'o', markeredgecolor = 'k', markersize = markersize,
                ls = '', zorder = 3)
        
        if Boxplot:
            ax.boxplot([parmsValues[i]], positions=[posParms[i]], widths=width*0.8, patch_artist=True,
                        showmeans=False, showfliers=False,
                        medianprops={"color": 'k', "linewidth": 2, 'alpha' : 0.75},
                        boxprops={"facecolor": color, "edgecolor": 'k',
                                  "linewidth": 1, 'alpha' : 0.5},
#                         boxprops={"color": color, "linewidth": 0.5},
                        whiskerprops={"color": 'k', "linewidth": 2, 'alpha' : 0.5},
                        capprops={"color": 'k', "linewidth": 2, 'alpha' : 0.5},
                         zorder = 1)
    
    for k in range(NValues):
        ax.plot(posParms, parmsValues[:,k], 'k', ls = '-', linewidth = 0.5, marker = '', zorder = 2)
        
    if stats:
        X = np.array(posParms)
        Y = parmsValues
        box_pairs = makeBoxPairs([i for i in range(NParms)])
        addStat_arrays(ax, X, Y, box_pairs, test = statMethod, percentHeight = 99)
        
    ax.set_xlabel('')
    ax.set_xlim([pos-0.5, pos+0.5])
    ax.set_xticks(posParms)
    ax.set_xticklabels(labels = tickLabels, fontsize = 10)
    ax.tick_params(axis='x', labelrotation = 10)
        
    return(fig, ax)



def D2Plot(data, fig = None, ax = None, XCol='', YCol='', CondCol='', Filters=[], 
           cellID='cellID', AvgPerCell=False, xscale = 'lin', yscale = 'lin', 
           figSizeFactor = 1, markers = [], markersizeFactor = 1):
    data_filtered = data
    for fltr in Filters:
        data_filtered = data_filtered.loc[fltr]
    
    NCond = len(CondCol)    
    if NCond == 1:
        CondCol = CondCol[0]
    elif NCond > 1:
        newColName = ''
        for i in range(NCond):
            newColName += CondCol[i]
            newColName += ' & '
        newColName = newColName[:-3]
        data_filtered[newColName] = ''
        for i in range(NCond):
            data_filtered[newColName] += data_filtered[CondCol[i]].astype(str)
            data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x + ' & ')
        data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x[:-3])
        CondCol = newColName
    
    if AvgPerCell:
        group = data_filtered.groupby(cellID)
        dictAggMean = getDictAggMean(data_filtered)
        data_filtered = group.agg(dictAggMean.pop(cellID)) #.reset_index(level=0, inplace=True)
        data_filtered.reset_index(level=0, inplace=True)
        
    Conditions = list(data_filtered[CondCol].unique())
    
    if fig == None:
        fig, ax = plt.subplots(1, 1, figsize = (8*figSizeFactor,5))
    else:
        pass
    
    markersize = 5 * markersizeFactor
    
    if xscale == 'log':
        ax.set_xscale('log')
    if yscale == 'log':
        ax.set_yscale('log')
    
    current_color_list = getStyleLists(Conditions, styleDict1).as_hex()
    cc = cycler(color=current_color_list)
    ax.set_prop_cycle(cc)
    
    im = 0
    for c in Conditions:
        Xraw = data_filtered[data_filtered[CondCol] == c][XCol].values
        Yraw = data_filtered[data_filtered[CondCol] == c][YCol].values
        XYraw = np.array([Xraw,Yraw]).T
        XY = XYraw[~np.isnan(XYraw).any(axis=1), :]
        X, Y = XY[:,0], XY[:,1]
        if len(markers) == 0:
            m = 'o'
        else:
            m = markers[im]
            im += 1
        
        if len(X) == 0:
            ax.plot([], [])
#             if modelFit:
#                 ax.plot([], [])
                
        elif len(X) > 0:                    
            ax.plot(X, Y, m, markersize = markersize, markeredgecolor='k', markeredgewidth = 1, label=c)
            
    ax.set_xlabel(XCol)
    ax.set_xlim([min(0,1.1*np.min(data_filtered[XCol])), 1.1*np.max(data_filtered[XCol])])
    ax.set_ylabel(YCol)
    if not yscale == 'log':
        ax.set_ylim([min(0,1.1*np.min(data_filtered[YCol])), 1.1*np.max(data_filtered[YCol])])
    ax.legend(loc='upper left')
    return(fig, ax)




def D2Plot_wFit(data, fig = None, ax = None, 
                XCol='', YCol='', CondCol='', 
                Filters=[], cellID='cellID', co_order = [],
                AvgPerCell=False, showManips = True,
                modelFit=False, modelType='y=ax+b',
                xscale = 'lin', yscale = 'lin', 
                figSizeFactor = 1, markersizeFactor = 1):
    
    data_filtered = data
    for fltr in Filters:
        data_filtered = data_filtered.loc[fltr]
    
    NCond = len(CondCol)    
    if NCond == 1:
        CondCol = CondCol[0]
    elif NCond > 1:
        newColName = ''
        for i in range(NCond):
            newColName += CondCol[i]
            newColName += ' & '
        newColName = newColName[:-3]
        data_filtered[newColName] = ''
        for i in range(NCond):
            data_filtered[newColName] += data_filtered[CondCol[i]].astype(str)
            data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x + ' & ')
        data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x[:-3])
        CondCol = newColName
    
    if AvgPerCell:
        group = data_filtered.groupby(cellID)
        dictAggMean = getDictAggMean(data_filtered)
        data_filtered = group.agg(dictAggMean.pop(cellID)) #.reset_index(level=0, inplace=True)
        data_filtered.reset_index(level=0, inplace=True)
        
    Conditions = list(data_filtered[CondCol].unique())
    
    if len(co_order) > 0:
        try:
            colorList, markerList = getStyleLists(co_order, styleDict1)
        except:
            colorList, markerList = colorList30, markerList10
    else:
        co_order = Conditions
        colorList, markerList = colorList30, markerList10
        
    colorList = [colorList40[32]]
    
    if fig == None:
        fig, ax = plt.subplots(1, 1, figsize = (8*figSizeFactor,5))
    else:
        pass
    
    markersize = 5 * markersizeFactor
    
    if xscale == 'log':
        ax.set_xscale('log')
    if yscale == 'log':
        ax.set_yscale('log')
    
#     current_color_list = getStyleLists(Conditions, styleDict1).as_hex()
#     cc = cycler(color=current_color_list)
#     ax.set_prop_cycle(cc)
    
#     if modelFit:
#         # Tweak the style cycle to plot for each condition: the points ('o') and then the fit ('-') with the same color.
# #         current_prop_cycle = plt.rcParams['axes.prop_cycle']
# #         current_color_list = prop_cycle.by_key()['color']
#         ncustom_color_list = list(np.array([current_color_list, current_color_list]).flatten(order='F'))
#         # if new_color_list was ['red', 'green', blue'], custom_color_list is now ['red', 'red', 'green', 'green', blue', blue']
#         cc = cycler(color=ncustom_color_list)
#         ax.set_prop_cycle(cc)
    
    for i in range(len(co_order)):
        c = co_order[i]
        color = colorList[i]
#         marker = my_default_marker_list[i]
        Xraw = data_filtered[data_filtered[CondCol] == c][XCol].values
        Yraw = data_filtered[data_filtered[CondCol] == c][YCol].values
        Mraw = data_filtered[data_filtered[CondCol] == c]['manipID'].values
        XYraw = np.array([Xraw,Yraw]).T
        XY = XYraw[~np.isnan(XYraw).any(axis=1), :]
        X, Y = XY[:,0], XY[:,1]
        M = Mraw[~np.isnan(XYraw).any(axis=1)]
        if len(X) == 0:
            ax.plot([], [])
            if modelFit:
                ax.plot([], [])
                
        elif len(X) > 0:
            eqnText = ''

            if modelFit:
                print('Fitting condition ' + c + ' with model ' + modelType)
                if modelType == 'y=ax+b':
                    params, results = jvu.fitLine(X, Y) # Y=a*X+b ; params[0] = b,  params[1] = a
                    pval = results.pvalues[1] # pvalue on the param 'a'
                    eqnText += " ; Y = {:.1f} X + {:.1f}".format(params[1], params[0])
                    eqnText += " ; p-val = {:.3f}".format(pval)
                    print("Y = {:.5} X + {:.5}".format(params[1], params[0]))
                    print("p-value on the 'a' coefficient: {:.4e}".format(pval))
                    # fitY = params[1]*X + params[0]
                    # imin = np.argmin(X)
                    # imax = np.argmax(X)
                    # ax.plot([X[imin],X[imax]], [fitY[imin],fitY[imax]], '--', lw = '1',
                    #         color = color, zorder = 4)
                    fitX = np.linspace(np.min(X), np.max(X), 100)
                    fitY = params[1]*fitX + params[0]
                    ax.plot(fitX, fitY, '--', lw = '1', 
                            color = color, zorder = 4)

                elif modelType == 'y=A*exp(kx)':
                    params, results = jvu.fitLine(X, np.log(Y)) # Y=a*X+b ; params[0] = b,  params[1] = a
                    pval = results.pvalues[1] # pvalue on the param 'k'
                    eqnText += " ; Y = {:.1f}*exp({:.1f}*X)".format(params[0], params[1])
                    eqnText += " ; p-val = {:.3f}".format(pval)
                    print("Y = {:.5}*exp({:.5}*X)".format(np.exp(params[0]), params[1]))
                    print("p-value on the 'k' coefficient: {:.4e}".format(pval))
                    # fitY = np.exp(params[0])*np.exp(params[1]*X)
                    # imin = np.argmin(X)
                    # imax = np.argmax(X)
                    # ax.plot([X[imin],X[imax]], [fitY[imin],fitY[imax]], '--', lw = '1',
                    #         color = color, zorder = 4)
                    fitX = np.linspace(np.min(X), np.max(X), 100)
                    fitY = np.exp(params[0])*np.exp(params[1]*fitX)
                    ax.plot(fitX, fitY, '--', lw = '1', 
                            color = color, zorder = 4)
                    
                elif modelType == 'y=k*x^a':
                    posValues = ((X > 0) & (Y > 0))
                    X, Y = X[posValues], Y[posValues]
                    params, results = jvu.fitLine(np.log(X), np.log(Y)) # Y=a*X+b ; params[0] = b,  params[1] = a
                    k = np.exp(params[0])
                    a = params[1]
                    R2 = results.rsquared
                    pval = results.pvalues[1] # pvalue on the param 'a'
                    eqnText += " ; Y = {:.1e} * X^{:.1f}".format(k, a)
                    eqnText += " ; p-val = {:.3f}".format(pval)
                    eqnText += " ; R2 = {:.2f}".format(R2)
                    print("Y = {:.4e} * X^{:.4f}".format(k, a))
                    print("p-value on the 'a' coefficient: {:.4e}".format(pval))
                    print("R2 of the fit: {:.4f}".format(R2))
                    # fitY = k * X**a
                    # imin = np.argmin(X)
                    # imax = np.argmax(X)
                    # ax.plot([X[imin],X[imax]], [fitY[imin],fitY[imax]], '--', lw = '1', 
                    #         color = color, zorder = 4)
                    fitX = np.linspace(np.min(X), np.max(X), 100)
                    fitY = k * fitX**a
                    ax.plot(fitX, fitY, '--', lw = '1', 
                            color = color, zorder = 4)
                
                print('Number of values : {:.0f}'.format(len(Y)))
                print('\n')
            
#                 if showManips:
#                     allManipID = list(data_filtered[data_filtered[CondCol] == c]['manipID'].unique())
#                     dictMarker = {}
#                     markers = []
#                     for mi in range(len(allManipID)):
#                         dictMarker[allManipID[mi]] = my_default_marker_list[mi]
#                     for k in range(len(M)):
#                         thisManipID = M[k]
#                         markers = dictMarker[thisManipID]
#                 else:
#                     markers = ['o' for i in range(len(M))]
#                 markers = np.array(markers)
                
            ax.plot(X, Y, 
                    color = color, ls = '',
                    marker = 'o', markersize = markersize, markeredgecolor='k', markeredgewidth = 1, 
                    label = c + eqnText)
            
            
    ax.set_xlabel(XCol)
    ax.set_xlim([0.9*np.min(data_filtered[XCol]), 1.1*np.max(data_filtered[XCol])])
    ax.set_ylabel(YCol)
    if not yscale == 'log':
        ax.set_ylim([0.9*np.min(data_filtered[YCol]), 1.1*np.max(data_filtered[YCol])])
    ax.legend(loc='upper left')
    
    
    return(fig, ax)


# These functions use the Bokeh library to display 1D categorical or 2D plots with interactive plots. They are less flexible but can be nice to explore the data set since you can display the cellID which is the source of each point by passing your pointer over it.



def D1PlotInteractive(data, CondCol='',Parameters=[],Filters=[],AvgPerCell=False,cellID='cellID'):
    data_filtered = data
    for fltr in Filters:
        data_filtered = data_filtered.loc[fltr]
        
#     print(data_filtered[cellID])
    if AvgPerCell:
        group = data_filtered.groupby(cellID)
        dictAggMean = getDictAggMean(data_filtered)
        data_filtered = group.agg(dictAggMean.pop(cellID)) #.reset_index(level=0, inplace=True)
        data_filtered.reset_index(level=0, inplace=True)
    
#     return(data_filtered)
    
    NPlots = len(Parameters)
    Conditions = list(data_filtered[CondCol].unique())
    if NPlots > 1:
        plots = []
        NCond = len(Conditions)
        data_filtered['X'] = 0
        data_filtered['X_jitter'] = 0.
        dictTicks = {}
        for i in range(NCond):
            mask = data_filtered[CondCol] == Conditions[i]
            data_filtered.loc[mask, 'X'] = i+1
            dictTicks[i+1] = Conditions[i]
        for i in data_filtered.index:
            data_filtered.loc[i, 'X_jitter'] = data_filtered.loc[i, 'X'] + 0.4*(np.random.rand(1)[0]-0.5)
        source = ColumnDataSource(
            data=data_filtered[[cellID]+[CondCol]+Parameters+['X','X_jitter']]
        )        
        
        for k in range(NPlots):
            hover = HoverTool(
                tooltips=[
                    ('Cell ID', "@"+cellID),
                    (Parameters[k], "@"+Parameters[k]),
                ]
            )
            index_cmap = factor_cmap(CondCol, palette=Category10[10], factors=sorted(data_filtered[CondCol].unique()), end=1)
            p = figure(plot_width=450, plot_height=500, tools=[hover], title="InteractivePlot") # 
            p.circle('X_jitter', Parameters[k], size=8, alpha = 0.6, source=source,fill_color=index_cmap,line_color='black')
            # Format
            p.x_range = Range1d(0, NCond+1)
            p.y_range = Range1d(min(0,1.1*np.min(data_filtered[Parameters[0]])), 1.1*np.max(data_filtered[Parameters[k]]))
            p.xaxis.ticker = [i for i in range(1,NCond+1)]
            p.xaxis.major_label_overrides = dictTicks
            p.xaxis.axis_label = CondCol
            p.xaxis.axis_label_text_font_size = '18pt'
            p.xaxis.major_label_text_font_size = '16pt'
            p.yaxis.axis_label = Parameters[k]
            p.yaxis.axis_label_text_font_size = '18pt'
            p.yaxis.major_label_text_font_size = '16pt'
            
            plots.append(p)
            
        p = gridplot(plots, ncols=2, toolbar_location=None)
        
        
    else:
        hover = HoverTool(
            tooltips=[
                ('Cell ID', "@"+cellID),
                (Parameters[0], "@"+Parameters[0]),
            ]
        )
        
        NCond = len(Conditions)
        data_filtered['X'] = 0
        data_filtered['X_jitter'] = 0.
        dictTicks = {}
        for i in range(NCond):
            mask = data_filtered[CondCol] == Conditions[i]
            data_filtered.loc[mask, 'X'] = i+1
            dictTicks[i+1] = Conditions[i]
        for i in data_filtered.index:
            data_filtered.loc[i, 'X_jitter'] = data_filtered.loc[i, 'X'] + 0.4*(np.random.rand(1)[0]-0.5)
        source = ColumnDataSource(
            data=data_filtered[[cellID]+[CondCol]+Parameters+['X','X_jitter']]
        )
        index_cmap = factor_cmap(CondCol, palette=Category10[10], factors=sorted(data_filtered[CondCol].unique()), end=1)
        TOOLS = "hover,pan,box_zoom,wheel_zoom,reset,save,help"
        p = figure(plot_width=500, plot_height=500, tools=TOOLS, title="InteractivePlot") # 
        p.circle('X_jitter', Parameters[0], size=8, alpha = 0.6, source=source,fill_color=index_cmap,line_color='black')
        # Format
        p.x_range = Range1d(0, NCond+1)
        p.y_range = Range1d(min(0,1.1*np.min(data_filtered[Parameters[0]])), 1.1*np.max(data_filtered[Parameters[0]]))
        p.xaxis.ticker = [i for i in range(1,NCond+1)]
        p.xaxis.major_label_overrides = dictTicks
        p.xaxis.axis_label = CondCol
        p.xaxis.axis_label_text_font_size = '18pt'
        p.xaxis.major_label_text_font_size = '16pt'
        p.yaxis.axis_label = Parameters[0]
        p.yaxis.axis_label_text_font_size = '18pt'
        p.yaxis.major_label_text_font_size = '16pt'
    return(p)



def D2PlotInteractive(data, XCol='',YCol='',CondCol='',Filters=[], cellID='cellID',AvgPerCell=False):
    
    data_filtered = data
    for fltr in Filters:
        data_filtered = data_filtered.loc[fltr]
        
    if AvgPerCell:
        group = data_filtered.groupby(cellID)
        dictAggMean = getDictAggMean(data_filtered)
        data_filtered = group.agg(dictAggMean.pop(cellID)) #.reset_index(level=0, inplace=True)
        data_filtered.reset_index(level=0, inplace=True)
    
    Conditions = list(data_filtered[CondCol].unique())

    NCond = len(Conditions)
    dictTicks = {}
    for i in range(NCond):
        dictTicks[i+1] = Conditions[i]
    
    source = ColumnDataSource(
        data=data_filtered[[cellID,CondCol,XCol,YCol]]
    )
    
    hover = HoverTool(
        tooltips=[
            ('Cell ID', "@"+cellID),
            (XCol, "@"+XCol),
            (YCol, "@"+YCol),
            (CondCol, "@"+CondCol),
        ]
    )
    
    index_cmap = factor_cmap(CondCol, palette=Category10[10], factors=sorted(data_filtered[CondCol].unique()), end=1)
    TOOLS = "pan,box_zoom,wheel_zoom,reset,save,help"
    p = figure(plot_width=900, plot_height=500, tools=TOOLS, title="InteractivePlot",toolbar_location="below") # 
    p.circle(XCol, YCol, size=8, alpha = 0.6, source=source,fill_color=index_cmap,line_color='black')
    p.add_tools(hover)
    # Format
    p.x_range = Range1d(0, 1.1*np.max(data_filtered[XCol]))
    p.y_range = Range1d(0, 1.1*np.max(data_filtered[YCol]))
    p.xaxis.axis_label = XCol
    p.xaxis.axis_label_text_font_size = '18pt'
    p.xaxis.major_label_text_font_size = '16pt'
    p.yaxis.axis_label = YCol
    p.yaxis.axis_label_text_font_size = '18pt'
    p.yaxis.major_label_text_font_size = '16pt'
    return(p)


# %%% Subfunctions



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


def getAggDf(df, cellID, CondCol, Variables):
    
    allCols = df.columns.values
    group = df.groupby(cellID)
    
    dictAggCategories = {}
    Categories = ['date', 'cellName', 'manipID', 'experimentType', 'drug', 
                  'substrate', 'cell type', 'cell subtype', 'bead type', 'bead diameter']
    if CondCol not in Categories:
        Categories += [CondCol]
    for c in Categories:
        if c in allCols:
            dictAggCategories[c] = 'first'
    dfCats = group.agg(dictAggCategories)
    
    dictAggCount = {'compNum' : 'count'}
    renameDict = {'compNum' : 'Ncomp'}
    dfCount = group.agg(dictAggCount).rename(columns=renameDict)
    
    dictAggMeans = {}
    renameDict = {}
    for v in Variables:
        dictAggMeans[v] = 'mean'
        renameDict[v] = v + '_mean'
    dfMeans = group.agg(dictAggMeans).rename(columns=renameDict)
        
    dictAggStds = {}
    renameDict = {}
    for v in Variables:
        dictAggStds[v] = 'std'
        renameDict[v] = v + '_std'
    dfStds = group.agg(dictAggStds).rename(columns=renameDict)
    
    if 'ctFieldThickness' in Variables:
        dfStds['ctFieldThickness_std'] = group.agg({'ctFieldFluctuAmpli' : 'median'}).values
        #ctFieldThickness ; ctFieldFluctuAmpli
    
#     for v in Variables:
        
    
    dfAgg = pd.merge(dfCats, dfCount, how="inner", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
    
    dfAgg = pd.merge(dfAgg, dfMeans, how="inner", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
    
    dfAgg = pd.merge(dfAgg, dfStds, how="inner", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
    
    return(dfAgg)


def getAggDf_weightedAvg(df, cellID, CondCol, Variables, WeightCols):
    
    def w_std(x, w):
        m = np.average(x, weights=w)
        v = np.average((x-m)**2, weights=w)
        std = v**0.5
        return(std)
    
    def nan2zero(x):
        if np.isnan(x):
            return(0)
        else:
            return(x)
    
    dictWeights = {}
    for i in range(len(Variables)):
        v = Variables[i]
        w = WeightCols[i]
        dictWeights[v] = w

    for v in Variables:
        df[v] = df[v].apply(nan2zero)
        df[dictWeights[v]] = df[dictWeights[v]].apply(nan2zero)
    cellIDvals = df['cellID'].unique()
    for v in Variables:
        for c in cellIDvals:
            if np.sum(df.loc[df['cellID'] == c, v].values) == 0:
                df.loc[df['cellID'] == c, dictWeights[v]].apply(lambda x: 0)
    
    
    
    allCols = df.columns.values
    group = df.groupby(cellID)
    
    dictAggCategories = {}
    Categories = ['date', 'cellName', 'manipID', 'experimentType', 'drug', 
                  'substrate', 'cell type', 'cell subtype', 'bead type', 'bead diameter']
    if CondCol not in Categories:
        Categories += [CondCol]
    for c in Categories:
        if c in allCols:
            dictAggCategories[c] = 'first'
    dfCats = group.agg(dictAggCategories)
    
    dictAggCount = {'compNum' : 'count'}
    renameDict = {'compNum' : 'Ncomp'}
    dfCount = group.agg(dictAggCount).rename(columns=renameDict) 
    
    dfAgg = pd.merge(dfCats, dfCount, how="inner", on='cellID',
    #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
    #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
    )
        
    dictAggMeans = {}
    renameDict = {}
    for i in range(len(Variables)):
        v = Variables[i]
        w = WeightCols[i]
        df_i = df[[cellID, v, w]]
        group_i = df_i.groupby(cellID)

        def WMfun(x):
            try:
                return(np.average(x, weights=df.loc[x.index, dictWeights[v]]))
            except:
                return(np.nan)
            
        def WSfun(x):
            try:
                return(w_std(x, df.loc[x.index, dictWeights[v]].values))
            except:
                return(np.nan)
        
        def Countfun(x):
            try:
                return(len(x.values[x.values != 0]))
            except:
                return(0)
        
        dictAggMeans[v] = WMfun
        renameDict[v] = v + '_Wmean'
        dfWMeans_i = group_i.agg({v : WMfun}).rename(columns={v : v + '_Wmean'})
        dfWStd_i = group_i.agg({v : WSfun}).rename(columns={v : v + '_Wstd'})
        dfWCount_i = group_i.agg({v : Countfun}).rename(columns={v : v + '_Count'})
        dfAgg = pd.merge(dfAgg, dfWMeans_i, how="inner", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        dfAgg = pd.merge(dfAgg, dfWStd_i, how="inner", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )
        
        dfAgg = pd.merge(dfAgg, dfWCount_i, how="inner", on='cellID',
        #     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
        #     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
        )

    return(dfAgg)


def makeOrder(*args):
    order = []
    listeTuple = list(itertools.product(*args, repeat=1))
    for tup in listeTuple:
        tmpText = ''
        for word in tup:
            tmpText += word
            tmpText += ' & '
        tmpText = tmpText[:-3]
        order.append(tmpText)
    return(order)


def makeBoxPairs(O):
    return(list(itertools.combinations(O, 2)))


def renameAxes(axes, rD):
    try:
        N = len(axes)
    except:
        axes = [axes]
        N = 1
    for i in range(N):
        # set xticks
        xticksTextObject = axes[i].get_xticklabels()
        xticksList = [xticksTextObject[j].get_text() for j in range(len(xticksTextObject))]
        test_hasXLabels = (len(''.join(xticksList)) > 0)
        if test_hasXLabels:
            newXticksList = [rD.get(k, k) for k in xticksList]
            axes[i].set_xticklabels(newXticksList)
        
        # set xlabel
        xlabel = axes[i].get_xlabel()
        newXlabel = rD.get(xlabel, xlabel)
        axes[i].set_xlabel(newXlabel)
        # set ylabel
        ylabel = axes[i].get_ylabel()
        newYlabel = rD.get(ylabel, ylabel)
        axes[i].set_ylabel(newYlabel)
        

def addStat_df(ax, data, box_pairs, param, cond, test = 'Mann-Whitney', percentHeight = 95):
    refHeight = np.percentile(data[param].values, percentHeight)
    currentHeight = refHeight
    scale = ax.get_yscale()
    xTicks = ax.get_xticklabels()
    dictXTicks = {xTicks[i].get_text() : xTicks[i].get_position()[0] for i in range(len(xTicks))}
    for bp in box_pairs:
        c1 = data[data[cond] == bp[0]][param].values
        c2 = data[data[cond] == bp[1]][param].values
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
    

def addStat_arrays(ax, X, Y, box_pairs, test = 'Mann-Whitney', percentHeight = 98):
    refHeight = np.percentile(Y, percentHeight)
    currentHeight = refHeight
    scale = ax.get_yscale()
    
    for bp in box_pairs:
        y1 = Y[bp[0], :]
        y2 = Y[bp[1], :]
        
        if test=='Mann-Whitney':
            statistic, pval = st.mannwhitneyu(y1, y2)
        elif test=='Wilcox_2s':
            statistic, pval = st.wilcoxon(y1, y2, alternative = 'two-sided')
        elif test=='Wilcox_greater':
            statistic, pval = st.wilcoxon(y1, y2, alternative = 'greater')
        elif test=='Wilcox_less':
            statistic, pval = st.wilcoxon(y1, y2, alternative = 'less')
        elif test=='t-test':
            statistic, pval = st.ttest_ind(y1, y2)
            
        text = 'ns'
        if pval < 0.05 and pval > 0.01:
            text = '*'
        elif pval < 0.01 and pval > 0.001:
            text = '**'
        elif pval < 0.001 and pval < 0.001:
            text = '***'
        elif pval < 0.0001:
            text = '****'
            
        ax.plot([X[bp[0]], X[bp[1]]], [currentHeight, currentHeight], 'k-', lw = 1)
        XposText = (X[bp[0]]+X[bp[1]])/2
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


def getStyleCycle(co_order, styleDict):
    colors = []
    markers = []
    linestyles = []
    linewidth = []

    for co in co_order:
        coStyle = styleDict[co]
        if 'color' in coStyle.keys():
            colors.append(coStyle['color'])
        else:
            colors.append('')
        if 'marker' in coStyle.keys():
            markers.append(coStyle['marker'])
        else:
            markers.append('')
        if 'linestyle' in coStyle.keys():
            linestyles.append(coStyle['marker'])
        else:
            linestyles.append('')
        if 'linewidth' in coStyle.keys():
            linewidth.append(coStyle['linewidth'])
        else:
            linewidth.append(1)
            
    cc = (cycler(color=colors) +
          cycler(linestyle=linestyles) +
          cycler(marker=markers) +
          cycler(linewidth=linewidth))
            
    return(cc)


def getStyleLists_Sns(co_order, styleDict):
    colors = []
    markers = []
    try:
        for co in co_order:
            coStyle = styleDict[co]
            if 'color' in coStyle.keys():
                colors.append(coStyle['color'])
            else:
                colors.append('')
            if 'marker' in coStyle.keys():
                markers.append(coStyle['marker'])
            else:
                markers.append('')
        palette = sns.color_palette(colors)
    except:
        palette = sns.color_palette()
    return(palette)

def getStyleLists(co_order, styleDict):
    colors = []
    markers = []
    for co in co_order:
        coStyle = styleDict[co]
        if 'color' in coStyle.keys():
            colors.append(coStyle['color'])
        else:
            colors.append('')
        if 'marker' in coStyle.keys():
            markers.append(coStyle['marker'])
        else:
            markers.append('')
    print(colors, markers)
    return(colors, markers)


def buildStyleDictMCA():
    # TBC
    styleDict = {}
    return(styleDict)


# %%% Tests of plotting functions


#### Test getAggDf_weightedAvg(df, cellID, CondCol, Variables, WeightCols)

# data = GlobalTable_meca_nonLin

# dates = ['22-02-09'] #['21-12-08', '22-01-12'] ['21-01-18', '21-01-21', '21-12-08']

# filterList = [(data['validatedThickness'] == True),
#               (data['cell subtype'] == 'aSFL'), 
#               (data['bead type'].apply(lambda x : x == list(['M450']))),
#               (data['substrate'] == '20um fibronectin discs'),
#               (data['date'].apply(lambda x : x in dates))]  # (data['validatedFit'] == True), 
# globalFilter = filterList[0]
# for k in range(1, len(filterList)):
#     globalFilter = globalFilter & filterList[k]

# data_f = data[globalFilter]

# fitMin = [S for S in range(25,1225,50)]
# fitMax = [S+150 for S in fitMin]
# fitCenters = np.array([S+75 for S in fitMin])
# regionFitsNames = [str(fitMin[ii]) + '<s<' + str(fitMax[ii]) for ii in range(len(fitMin))]

# listColumnsMeca = []

# KChadwick_Cols = []
# KWeight_Cols = []

# for rFN in regionFitsNames:
#     listColumnsMeca += ['KChadwick_'+rFN, 'K_CIW_'+rFN, 'R2Chadwick_'+rFN, 'K2Chadwick_'+rFN, 
#                         'H0Chadwick_'+rFN, 'Npts_'+rFN, 'validatedFit_'+rFN]
#     KChadwick_Cols += [('KChadwick_'+rFN)]

#     K_CIWidth = data_f['K_CIW_'+rFN] #.apply(lambda x : x.strip('][').split(', ')).apply(lambda x : (np.abs(float(x[0]) - float(x[1]))))
#     KWeight = (data_f['KChadwick_'+rFN]/K_CIWidth) # **2
#     data_f['K_Weight_'+rFN] = KWeight
#     data_f['K_Weight_'+rFN] *= data_f['KChadwick_'+rFN].apply(lambda x : (x<1e6))
#     data_f['K_Weight_'+rFN] *= data_f['R2Chadwick_'+rFN].apply(lambda x : (x>1e-2))
#     data_f['K_Weight_'+rFN] *= data_f['K_CIW_'+rFN].apply(lambda x : (x!=0))
#     KWeight_Cols += [('K_Weight_'+rFN)]
    

# # dictWeights = {}
# # for i in range(len(Variables)):
# #     v = Variables[i]
# #     w = WeightCols[i]
# #     dictWeights[v] = w
# # def nan2zero(x):
# #     if np.isnan(x):
# #         return(0)
# #     else:
# #         return(x)

# # for v in Variables:
# #     df[v] = df[v].apply(nan2zero)
# #     df[dictWeights[v]] = df[dictWeights[v]].apply(nan2zero)
# # cellIDvals = df['cellID'].unique()
# # for v in Variables:
# #     for c in cellIDvals:
# #         if np.sum(df.loc[df['cellID'] == c, v].values) == 0:
# #             df.loc[df['cellID'] == c, dictWeights[v]].apply(lambda x: 1)
            
# # df

# CondCol = 'date'
# Variables = KChadwick_Cols
# WeightCols = KWeight_Cols
# data_f_agg = getAggDf_weightedAvg(data_f, 'cellID', CondCol, Variables, WeightCols)
# data_f_agg


#### Test getAggDf(df, cellID, CondCol, Variables)

# data = GlobalTable_meca_Py2

# Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True)]

# data_f = data
# for fltr in Filters:
#     data_f = data_f.loc[fltr]

# # dfA = getAggDf(data_f, 'cellID', 'bead type', ['surroundingThickness', 'EChadwick'])
# dfA = getAggDf(data_f, 'cellID', 'bead type', ['ctFieldThickness', 'EChadwick'])
# dfA


#### Test the D2Plot_wFit
# data = GlobalTable_meca_Py2

# Filters = [(data['validatedFit'] == True), 
#            (data['validatedThickness'] == True), 
#            (data['substrate'] == '20um fibronectin discs'),
#            (data['date'].apply(lambda x : x in ['21-12-16','21-12-08']))]

# fig, ax = D2Plot_wFit(data, XCol='ctFieldThickness', YCol='EChadwick', CondCol = ['bead type'],
#            Filters=Filters, cellID = 'cellID', AvgPerCell=True, xscale = 'log', yscale = 'log', 
#            modelFit=True, modelType = 'y=k*x^a')

# ax.set_ylabel('EChadwick (Pa)')
# ax.set_xlabel('Thickness at low force (nm)')
# fig.suptitle('3T3aSFL: E(h)')
# ax.legend(loc = 'upper right')

# # jvu.archiveFig(fig, ax, name='E(h)_3T3aSFL_Dec21_M450-5-13_vs_M270-14-54', figDir = todayFigDirLocal + '//' + figSubDir)
# plt.show()


#### Test of the averaging per cell routine

# data = GlobalTable_meca
# CondCol='drug'
# Parameters=['SurroundingThickness','EChadwick']
# Filters = [(GlobalTable_meca['Validated'] == 1)]
# AvgPerCell=True
# cellID='CellName'

# data_filtered = data
# for fltr in Filters:
#     data_filtered = data_filtered.loc[fltr]

# group = data_filtered.groupby(cellID)
# dictAggMean = getDictAggMean(data_filtered)
# data_filtered = group.agg(dictAggMean.pop(cellID)) #.reset_index(level=0, inplace=True)
# data_filtered.reset_index(level=0, inplace=True)
# data_filtered=data_filtered[[cellID]+[CondCol]+Parameters]
# print(data_filtered)


#### Test of a routine to remove points of a list of XY positions where at least 1 of the coordinates is 'nan'

# XYraw = np.array([[np.nan, 2, 3, np.nan, 5], [10,20,30,40,50]])
# XYraw = XYraw.T
# XY = XYraw[~np.isnan(XYraw).any(axis=1), :]
# X, Y = XY[:,0], XY[:,1]
# X, Y


#### Test of a routine to double each element in a list ; example [1, 2, 3] -> [1, 1, 2, 2, 3, 3]

# newnew_color_list = np.array([new_color_list, new_color_list])
# custom_color_list = list(np.array([new_color_list, new_color_list]).flatten(order='F'))
# custom_color_list


#### Test of makeOrder function

# print(makeOrder(['none','doxycyclin'],['BSA coated glass','20um fibronectin discs']))
# print(makeOrder(['A','B']))
# print(makeOrder(['A','B'], ['C','D']))
# print(makeOrder(['A','B'], ['C','D'], ['E','F']))


#### Test of makeBoxPairs function

# O = makeOrder(['none','doxycyclin'],['BSA coated glass','20um fibronectin discs'])
# makeBoxPairs(O)


#### Test of custom cycles

# cc = (cycler(color=list('rgb')) +
#       cycler(linestyle=['', '', '']) +
#       cycler(marker=['o','*','<']) +
#      cycler(linewidth=[1,1,1]))
# cc

# plt.rcParams['axes.prop_cycle'] = cc

# fig, ax = plt.subplots()
# for i in range(10):
#     ax.plot([0,1], [0,i])
# plt.show()

















# %% Plots


# %%% Valentin's data

# %%%%% Import

rawMecaTable = jva.getGlobalTable_meca('Global_MecaData')

# rawMecaTable.head()

# rawMecaTable.loc[]
# Agréger tous les M450 WT sous un même nom


# %%%%% M450, various rates

Filters = [(rawMecaTable['Validated'] == 1), 
           ((rawMecaTable['ExpType'] == 'DictyDB_M450'))] #  | (rawMecaTable['ExpType'] == 'DictyDB_M450-Multi')
co_order = makeOrder(['02s', '05s', '1s', '2s', '4s', '7s', '12s']) # co_order = co_order,
fig, ax00 = D1Plot(rawMecaTable, CondCol=['TpsComp'],Parameters=['SurroundingThickness','EChadwick'], 
                   co_order = co_order, Filters=Filters, AvgPerCell=False, cellID='CellName',
                   stats=False, figSizeFactor = 1.8, markersizeFactor=0.5)
fig.suptitle('M450, various rates')
renameAxes(ax00,renameDict1)
plt.show()


# %%%%% M270, various rates

Filters = [(rawMecaTable['Validated'] == 1), 
           ((rawMecaTable['ExpType'] == 'DictyDB_M270'))]
co_order = makeOrder(['02s', '05s', '1s', '2s', '4s', '7s', '12s'])
fig, ax01 = D1Plot(rawMecaTable, CondCol=['TpsComp'],Parameters=['SurroundingThickness','EChadwick'],
                   co_order = co_order, Filters=Filters, AvgPerCell=False, cellID='CellName', 
                   stats=False, figSizeFactor = 1.8,markersizeFactor=0.5)
fig.suptitle('M270, various rates')
plt.show()

# %%%%% 1s compressions, dictys, sizes of beads
rawMecaTable2 = rawMecaTable

rawMecaTable2.loc[rawMecaTable2['ExpType'] == 'DictyDB_M450-Multi', 'ExpType'] = 'DictyDB_M450'


data = rawMecaTable2
ExpTypes = ['DictyDB_M270', 'DictyDB_M450']
Filters = [(data['Validated'] == 1),
           (data['SurroundingThickness'] <= 900), 
           (data['TpsComp'] == '1s'),
           (data['ExpType'].apply(lambda x : x in ExpTypes))]
co_order = makeOrder(['DictyDB_M270', 'DictyDB_M450'])

fig, ax01 = D1Plot(data, CondCol=['ExpType'],Parameters=['SurroundingThickness','EChadwick'],
                   co_order = co_order, Filters=Filters, AvgPerCell=True, cellID='CellName', 
                   stats=True, figSizeFactor = 1.0, markersizeFactor=1)

rD = {'DictyDB_M270':'M270', 'DictyDB_M450':'M450',
      'SurroundingThickness' : 'Median thickness (nm)',
      'EChadwick' : 'Elastic modulus (Pa)'}
renameAxes(ax01, rD)
fig.suptitle('Two sizes of beads\nDictyostelium cortices')
# jvu.archiveFig(fig, ax, todayFigDirLocal + '//BeadSizes', name='BeadSizeDicty', dpi = 100)
plt.show()

# %%%%% 1s compressions, dictys, sizes of beads
rawMecaTable2 = rawMecaTable

rawMecaTable2.loc[rawMecaTable2['ExpType'] == 'DictyDB_M450-Multi', 'ExpType'] = 'DictyDB_M450'


data = rawMecaTable2
ExpTypes = ['DictyDB_M270', 'DictyDB_M450']
Filters = [(data['Validated'] == 1),
           (data['SurroundingThickness'] <= 900), 
           (data['TpsComp'] == '1s'),
           (data['ExpType'].apply(lambda x : x in ExpTypes))]
co_order = makeOrder(['DictyDB_M270', 'DictyDB_M450'])

fig, ax01 = D1Plot(data, CondCol=['ExpType'],Parameters=['EChadwick'],
                   co_order = co_order, Filters=Filters, AvgPerCell=False, cellID='CellName', 
                   stats=True, figSizeFactor = 1.0, markersizeFactor=1)

rD = {'DictyDB_M270':'M270', 'DictyDB_M450':'M450',
      'SurroundingThickness' : 'Median thickness (nm)',
      'EChadwick' : 'Elastic modulus (Pa)'}
renameAxes(ax01, rD)
fig.suptitle('Two sizes of beads\nDictyostelium cortices')
# jvu.archiveFig(fig, ax, todayFigDirLocal + '//BeadSizes', name='BeadSizeDicty_2', dpi = 100)
plt.show()

# %%% MDA project - Matlab processing

# %%%%% 3T3aSFL on patterns: Ct Field


Filters = [(GlobalTable_ctField['validated'] == True), 
           (GlobalTable_ctField['medianThickness'] <= 1000)]

fig, ax = D1Plot(GlobalTable_ctField, CondCol=['drug','substrate'],
                 Parameters=['medianThickness','fluctuAmpli'],
                 Filters=Filters)

fig.suptitle('3T3aSFL on patterns: Ct Field')
plt.show()


# %%%%% 3T3aSFL on diverse substrates: Compressions


Filters = [(GlobalTable_meca['Validated'] == 1), 
           (GlobalTable_meca['cell subtype'] == 'aSFL')]

co_order = makeOrder(['none','doxycyclin'],['BSA coated glass','20um fibronectin discs'])
fig, ax = D1Plot(GlobalTable_meca, CondCol=['drug','substrate'],
                 Parameters=['SurroundingThickness','EChadwick'],
                 Filters=Filters,AvgPerCell=False,cellID='CellName', co_order=co_order,
                 markersizeFactor = 0.6)

fig.suptitle('3T3aSFL on diverse substrates: Compressions')
plt.show()

# %%%%% 3T3aSFL on patterns: Compressions


Filters = [(GlobalTable_meca['Validated'] == 1),
           (GlobalTable_meca['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca['cell subtype'] == 'aSFL-6FP')]

fig, ax = D1Plot(GlobalTable_meca, CondCol=['drug','substrate'],
                 Parameters=['SurroundingThickness','EChadwick'],
                 Filters=Filters, AvgPerCell=False, cellID='CellName',
                 markersizeFactor = 1, orientation = 'v')

fig.suptitle('3T3aSFL-6FP on patterns: Compressions')
plt.show()


# %%%%% 3T3aSFL on diverse substrates: Compressions


Filters = [(GlobalTable_meca['Validated'] == 1), 
           (GlobalTable_meca['substrate'] == '20um fibronectin discs')]

co_order = makeOrder([['aSFL','aSFL-6FP'],['none','doxycyclin']])
fig, ax = D1Plot(GlobalTable_meca, CondCol=['cell subtype','drug'],
                 Parameters=['SurroundingThickness','EChadwick'],
                 Filters=Filters,AvgPerCell=True,cellID='CellName',co_order=co_order)
fig.suptitle('3T3aSFL on diverse substrates: Compressions')
plt.show()


# %%%%% 3T3aSFL SHORT vs LONG linker: Compressions


Filters = [(GlobalTable_meca['Validated'] == 1), 
           (GlobalTable_meca['substrate'] == '20um fibronectin discs')]

co_order = makeOrder([['aSFL','aSFL-6FP'],['none','doxycyclin']])
fig, ax = D1Plot(GlobalTable_meca, CondCol=['cell subtype','drug'],
                 Parameters=['SurroundingThickness','EChadwick'],
                 Filters=Filters,AvgPerCell=True,cellID='CellName',co_order=co_order)
fig.suptitle('3T3aSFL SHORT vs LONG linker: Compressions')
plt.show()


# %%%%% 3T3aSFL - Dh = f(H)


Filters = [(GlobalTable_ctField['validated'] == True), 
           (GlobalTable_ctField['cell subtype'] == 'aSFL')]

fig, ax = D2Plot_wFit(GlobalTable_ctField, XCol='medianThickness',YCol='fluctuAmpli',
                 CondCol = ['drug'], Filters=Filters, modelFit=False)
fig.suptitle('3T3aSFL - Dh = f(H)')
# jvu.archiveFig(fig, ax, name='aSFL_Dh(h)_drug', figDir = todayFigDirLocal + '//' + 'ThicknessPlots')
plt.show()


# %%%%% 3T3aSFL - Dh = f(H) - medianThickness <= 700

# Same as above without the 2 thickest cells
Filters = [(GlobalTable_ctField['validated'] == True), 
           (GlobalTable_ctField['cell subtype'] == 'aSFL'), 
           (GlobalTable_ctField['medianThickness'] <= 700)]

fig, ax = D2Plot_wFit(GlobalTable_ctField, XCol='medianThickness',YCol='fluctuAmpli',
                 CondCol = ['drug'], Filters=Filters, modelFit=False)
fig.suptitle('3T3aSFL - Dh = f(H)')
# jvu.archiveFig(fig, ax, name='aSFL_Dh(h)_drug_wo2LastPoints', figDir = todayFigDirLocal + '//' + 'ThicknessPlots')
plt.show()


# %%%%% 3T3aSFL: E(h) - glass & fibro

Filters = [(GlobalTable_meca['Validated'] == True), 
           (GlobalTable_meca['cell subtype'] == 'aSFL')]

fig, ax = D2Plot_wFit(GlobalTable_meca, XCol='SurroundingThickness',YCol='EChadwick',
                 CondCol = ['substrate','drug'],Filters=Filters, cellID = 'CellName', 
                 AvgPerCell=True, modelFit=False, modelType='y=A*exp(kx)', yscale = 'log')

fig.suptitle('3T3aSFL: E(h)')
renameAxes(ax,renameDict1)
ax.legend(loc='upper right')
# jvu.archiveFig(fig, ax, name='aSFL_E(h)_drug&substrate_01', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% 3T3aSFL: E(h) - only fibro

Filters = [(GlobalTable_meca['Validated'] == True), 
           (GlobalTable_meca['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca['substrate'] == '20um fibronectin discs')]

fig, ax = D2Plot_wFit(GlobalTable_meca, XCol='SurroundingThickness',YCol='EChadwick',
                 CondCol = ['substrate','drug'],Filters=Filters,
                 cellID = 'CellName', AvgPerCell=False, modelFit=False, 
                 modelType='y=A*exp(kx)', xscale = 'log', yscale = 'log')

fig.suptitle('3T3aSFL: E(h)')
renameAxes(ax,renameDict1)
ax.set_xlim([90, 1500])
ax.legend(loc='upper right')
# jvu.archiveFig(fig, ax, name='aSFL_E(h)_drug&substrate_00_allComp', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% 3T3aSFL: E(h) - only fibro - several cell subtype

Filters = [(GlobalTable_meca['Validated'] == True), 
           (GlobalTable_meca['substrate'] == '20um fibronectin discs')]

fig, ax = D2Plot_wFit(GlobalTable_meca, XCol='SurroundingThickness',YCol='EChadwick',
                 CondCol = ['cell subtype','drug'], Filters=Filters, 
                 cellID = 'CellName', AvgPerCell=True, modelFit=False, 
                 modelType='y=A*exp(kx)', yscale = 'log')

fig.suptitle('3T3aSFL: E(h)')
renameAxes(ax,renameDict1)
ax.legend(loc='upper right')
# jvu.archiveFig(fig, ax, name='aSFL_E(h)_drug&substrate_02', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% 3T3aSFL: E(h) - only fibro - several cell subtype - all comp

Filters = [(GlobalTable_meca['Validated'] == True), 
           (GlobalTable_meca['substrate'] == '20um fibronectin discs')]

fig, ax = D2Plot_wFit(GlobalTable_meca, XCol='SurroundingThickness',YCol='EChadwick',
                 CondCol = ['cell subtype','drug'],Filters=Filters, 
                 cellID = 'CellName', AvgPerCell=False, modelFit=False, 
                 modelType='y=A*exp(kx)', yscale = 'log')

fig.suptitle('3T3aSFL: E(h)')
renameAxes(ax,renameDict1)
ax.legend(loc='upper right')
# jvu.archiveFig(fig, ax, name='aSFL_E(h)_drug&substrate_02_allComp', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% Some test regarding plot style

co_o = ['BSA coated glass & none', 'BSA coated glass & doxycyclin', 
        '20um fibronectin discs & doxycyclin', '20um fibronectin discs & none']

sns.color_palette(['#ff9896', '#d62728', '#1f77b4', '#aec7e8'])

# getStyleLists(co_o, styleDict1)

# xticksTextObject = ax.get_xticklabels()
# xticksList = [xticksTextObject[j].get_text() for j in range(len(xticksTextObject))]
# newXticksList = [renameDict1.get(k, k) for k in xticksList]
# newXticksList


# %%%%% 3T3aSFL - All linker types: All Compressions

Filters = [(GlobalTable_meca_Py['validatedFit'] == True), 
           (GlobalTable_meca_Py['validatedThickness'] == True), 
           (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs')]

co_order = makeOrder(['aSFL','aSFL-A8','aSFL-6FP'], ['none','doxycyclin'])

box_pairs=[('aSFL & none', 'aSFL & doxycyclin'),
 ('aSFL-A8 & none', 'aSFL-A8 & doxycyclin'),
 ('aSFL-6FP & none', 'aSFL-6FP & doxycyclin'),
 ('aSFL & none', 'aSFL-A8 & none'),
 ('aSFL & none', 'aSFL-6FP & none')]

fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype','drug'],
                 Parameters=['surroundingThickness','EChadwick'],Filters=Filters,
                 AvgPerCell=False,cellID='cellID',co_order=co_order,
                 box_pairs=box_pairs,stats=True,markersizeFactor = 0.5,
                 orientation = 'v')

renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - All linker types: All Compressions')
plt.show()

# %%% MCA project - Py1 > half python half matlab processing

figSubDir = 'MCAproject'
# print(pairedPalette)
# pairedPalette


# %%%% 3T3 aSFL (A11) on and off pattern, with or without doxy -- (january 2021)

# %%%%%


data = GlobalTable_meca # 100% matlab

Filters = [(data['Validated'] == 1), (data['cell subtype'] == 'aSFL'), 
           (data['substrate'] == 'BSA coated glass'), (data['SurroundingThickness'] > 0)]

co_order = makeOrder(['BSA coated glass'],['none','doxycyclin'])

fig, ax = D1Plot(data, CondCol=['substrate','drug'],Parameters=['SurroundingThickness','EChadwick'], 
                 Filters=Filters,AvgPerCell=False, cellID='CellName', co_order=co_order, 
                 figSizeFactor = 0.5)

renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL non-adherent: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFLonBSA_drug_SurroundingThickness&EChadwick_allComp', figSubDir = figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca # 100% matlab
dates = ['20-08-04', '20-08-05', '20-08-07', '21-01-18', '21-01-21']

Filters = [(data['Validated'] == 1), (data['cell subtype'] == 'aSFL'), 
           (data['SurroundingThickness'] > 0), (data['date'].apply(lambda x : x in dates))]

co_order = makeOrder(['BSA coated glass','20um fibronectin discs'],['none','doxycyclin'])
fig, ax = D1Plot(GlobalTable_meca, CondCol=['substrate','drug'],
                 Parameters=['SurroundingThickness','EChadwick'],
                 Filters=Filters,AvgPerCell=True, cellID='CellName', 
                 co_order=co_order, figSizeFactor = 0.8)

renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL on diverse substrates: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFL_substrate&drug_SurroundingThickness&EChadwick', figSubDir = figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py # Tracking matlab, postprocessing python
dates = ['20-08-04', '20-08-05', '20-08-07', '21-01-18', '21-01-21']

Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True), 
           (data['cell subtype'] == 'aSFL'), (data['date'].apply(lambda x : x in dates))]

co_order = makeOrder(['BSA coated glass','20um fibronectin discs'],['none','doxycyclin'])
fig, ax = D1Plot(data, CondCol=['substrate','drug'],Parameters=['surroundingThickness','EChadwick'],
                 Filters=Filters,AvgPerCell=True, cellID='cellID', co_order=co_order, 
                 figSizeFactor = 0.8,useHue = False, orientation = 'v')

renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL on diverse substrates: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFL_substrate&drug_SurroundingThickness&EChadwick_NEWTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py # Tracking matlab, postprocessing python
dates = ['20-08-04', '20-08-05', '20-08-07', '21-01-18', '21-01-21']

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['cell subtype'] == 'aSFL'),
           (data['substrate'] == '20um fibronectin discs'),
           (data['bead type'] == 'M450'),
           (data['ctFieldThickness'] <= 600),
           (data['date'].apply(lambda x : x in dates))]

co_order = makeOrder(['none','doxycyclin'])
fig, ax = D1Plot(data, CondCol=['drug'],Parameters=['ctFieldThickness','EChadwick'],Filters=Filters,                 AvgPerCell=True, cellID='cellID', co_order=co_order, figSizeFactor = 1.0,
                 useHue = False, orientation = 'h', markersizeFactor=1.2)

ax[0].set_ylim([0, ax[0].get_ylim()[1]])
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL on fibronectin discs: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFL_substrate&drug_SurroundingThickness&EChadwick_NEWTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


data = GlobalTable_ctField # Tracking matlab, postprocessing python
dates = []

Filters = [(data['validated'] == True), (data['medianThickness'] <= 1000)]

co_order = makeOrder(['20um fibronectin discs'],['none','doxycyclin'])
fig, ax = D1Plot(data, CondCol=['substrate','drug'],Parameters=['medianThickness','fluctuAmpli'],                 Filters=Filters,stats=True,co_order=co_order,figSizeFactor=0.5)

ax[0].set_ylim([0,600])
ax[1].set_ylim([0,600])
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL on patterns: Constant Field')
# jvu.archiveFig(fig, ax, name='3T3aSFL_drug_medianThickness', figSubDir = figSubDir)

plt.show()


# %%%%%


data = GlobalTable_meca # 100%  matlab
dates = ['21-01-18', '21-01-21']

Filters = [(data['Validated'] == True), 
           (data['date'].apply(lambda x : x in dates))]

# def D2Plot_wFit(data, XCol='', YCol='', CondCol='', Filters=[], cellID='cellID', 
#                 AvgPerCell=False, showManips = True,
#                 modelFit=False, modelType='y=ax+b',
#                 xscale = 'lin', yscale = 'lin', 
#                 figSizeFactor = 1):

fig, ax = D2Plot_wFit(data, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['drug'], 
                      Filters=Filters, cellID = 'CellName', AvgPerCell=True, modelFit=False, 
                      modelType = 'y=ax+b')

renameAxes(ax,renameDict1)
fig.suptitle('aSFL expressing linker: E(fluo)')
ax.set_ylim([0,80000])
ax.set_xlim([-10,1000])
yt = ax.get_yticks()
ax.set_yticklabels((yt/1000).astype(int))
ax.set_ylabel('E Chadwick (kPa)')

# jvu.archiveFig(fig, ax, name='aSFL_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


data = GlobalTable_meca # 100%  matlab
dates = ['21-01-18', '21-01-21']

#Same as above without the lonely point
Filters = [(data['Validated'] == True), 
           (data['date'].apply(lambda x : x in dates)), 
           (data['EChadwick'] <= 80000)]

fig, ax = D2Plot_wFit(data, XCol='meanFluoPeakAmplitude', YCol='EChadwick', 
                      CondCol = ['cell subtype', 'drug'], Filters=Filters, 
                      cellID = 'CellName', AvgPerCell=True, modelFit=True)

renameAxes(ax,renameDict1)
fig.suptitle('aSFL expressing linker: E(fluo)')

# jvu.archiveFig(fig, ax, name='aSFL_iMC_E(fluo)_woLastPoint', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL')]
fig, ax = D2Plot(GlobalTable_meca, XCol='meanFluoPeakAmplitude',YCol='SurroundingThickness',
                 CondCol = ['cell subtype', 'drug'], Filters=Filters, cellID = 'CellName',AvgPerCell=True)

renameAxes(ax,renameDict1)
fig.suptitle('aSFL expressing linker: H(fluo)')

# jvu.archiveFig(fig, ax, name='aSFL_iMC_H(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Filters = [(GlobalTable_ctField['validated'] == True)]
fig, ax = D2Plot(GlobalTable_ctField, XCol='meanFluoPeakAmplitude',YCol='medianThickness',
                 CondCol = ['cell subtype', 'drug'],Filters=Filters, cellID = 'cellID', AvgPerCell=True)

renameAxes(ax,renameDict1)
fig.suptitle('aSFL expressing linker: medianH(fluo)')

# jvu.archiveFig(fig, ax, name='aSFL_iMC_medianH(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%% 3T3 aSFL : standard line (A11) versus different clones (A8 - high exp ; 6FP - long linker) -- First round (april - june 2021)

# %%%%%


Filters = [(GlobalTable_meca_Py['validatedThickness'] == True), (GlobalTable_meca_Py['ctFieldThickness'] <= 1000)]
co_order = makeOrder(['aSFL','aSFL-A8','aSFL-6FP'],['none','doxycyclin'])
box_pairs=[('aSFL & none', 'aSFL & doxycyclin'),
 ('aSFL-A8 & none', 'aSFL-A8 & doxycyclin'),
 ('aSFL-6FP & none', 'aSFL-6FP & doxycyclin'),
 ('aSFL & none', 'aSFL-A8 & none'),
 ('aSFL & none', 'aSFL-6FP & none')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype','drug'],Parameters=['ctFieldThickness','ctFieldFluctuAmpli'],                 Filters=Filters,AvgPerCell=True,stats=True,co_order=co_order,box_pairs=box_pairs,figSizeFactor=1)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL on patterns: H and DH from meca expe')
# jvu.archiveFig(fig, ax, name='3T3aSFL_drug_medianThickness_fromMeca_PYTHONTABLE', figSubDir = figSubDir)

plt.show()


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == 1), (GlobalTable_meca['substrate'] == '20um fibronectin discs')]
co_order = makeOrder(['aSFL','aSFL-6FP'],['none','doxycyclin'])
fig, ax = D1Plot(GlobalTable_meca, CondCol=['cell subtype','drug'],Parameters=['SurroundingThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='CellName',co_order=co_order,stats=True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL short vs long linker: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFL_likerType&drug_SurroundingThickness&EChadwick', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['substrate'] == '20um fibronectin discs')]
co_order = makeOrder(['aSFL','aSFL-A8','aSFL-6FP'],['none','doxycyclin'])
box_pairs=[('aSFL & none', 'aSFL & doxycyclin'),
 ('aSFL-A8 & none', 'aSFL-A8 & doxycyclin'),
 ('aSFL-6FP & none', 'aSFL-6FP & doxycyclin'),
 ('aSFL & none', 'aSFL-A8 & none'),
 ('aSFL & none', 'aSFL-6FP & none')]
fig, ax = D1Plot(GlobalTable_meca, CondCol=['cell subtype','drug'],Parameters=['SurroundingThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='CellName',co_order=co_order,box_pairs=box_pairs,stats=True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - All linker types: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFL_likerType&drug_SurroundingThickness&EChadwick', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs')]
co_order = makeOrder(['aSFL','aSFL-A8','aSFL-6FP'],['none','doxycyclin'])
box_pairs=[('aSFL & none', 'aSFL & doxycyclin'),
 ('aSFL-A8 & none', 'aSFL-A8 & doxycyclin'),
 ('aSFL-6FP & none', 'aSFL-6FP & doxycyclin'),
 ('aSFL & none', 'aSFL-A8 & none'),
 ('aSFL & none', 'aSFL-6FP & none')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype','drug'],Parameters=['ctFieldThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,box_pairs=box_pairs,stats=True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - All linker types: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFL_likerType&drug_ctFThickness&EChadwick_PYTHONTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), (GlobalTable_meca_Py['drug'] == 'doxycyclin'),           (pd.isna(GlobalTable_meca_Py['meanFluoPeakAmplitude']) != True)]
co_order = makeOrder(['aSFL','aSFL-A8','aSFL-6FP'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype'],Parameters=['meanFluoPeakAmplitude'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True,figSizeFactor=1.25)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL different cell lines\nLinker expression quantif by fluo')
jvu.archiveFig(fig, ax, name='3T3aSFL_likerType&_fluoExp_PYTHONTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), ((GlobalTable_meca_Py['cell subtype'] == 'aSFL') | (GlobalTable_meca_Py['cell subtype'] == 'aSFL-A8')), (GlobalTable_meca_Py['drug'] == 'doxycyclin')]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype'],                  Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=True)
renameAxes(ax,renameDict1)
fig.suptitle('aSFL & aSFL-A8 expressing linker: E(fluo)')
# jvu.archiveFig(fig, ax, name='aSFL&A8_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


#Same as above without the lonely point
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['drug'] == 'doxycyclin'),            ((GlobalTable_meca_Py['cell subtype'] == 'aSFL') | (GlobalTable_meca_Py['cell subtype'] == 'aSFL-A8')),            (GlobalTable_meca_Py['meanFluoPeakAmplitude'] <= 1200), (GlobalTable_meca_Py['EChadwick'] >= 2000)]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype'],                  Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=True)
renameAxes(ax,renameDict1)
fig.suptitle('aSFL & aSFL-A8 expressing linker: E(fluo)')
# jvu.archiveFig(fig, ax, name='aSFL&A8_iMC_E(fluo)_fluo-1200_&_E+2000', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL-6FP')]
fig, ax = D2Plot(GlobalTable_meca, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype', 'drug'],                  Filters=Filters, cellID = 'CellName', AvgPerCell=True, modelFit=True)
renameAxes(ax,renameDict1)
fig.suptitle('aSFL-6FP expressing long linker: E(fluo)')
# jvu.archiveFig(fig, ax, name='aSFL-6FP_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


#Same as above without the lonely point
Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL-6FP'), (GlobalTable_meca['meanFluoPeakAmplitude'] <= 1000)]
fig, ax = D2Plot(GlobalTable_meca, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype', 'drug'],                  Filters=Filters, cellID = 'CellName', AvgPerCell=True, modelFit=True)
renameAxes(ax,renameDict1)
fig.suptitle('aSFL-6FP expressing long linker: E(fluo)')
# jvu.archiveFig(fig, ax, name='aSFL-6FP_iMC_E(fluo)_woLastPoint', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL-6FP')]
fig, ax = D2Plot(GlobalTable_meca, XCol='meanFluoPeakAmplitude', YCol='SurroundingThickness', CondCol = ['cell subtype', 'drug'],                 Filters=Filters, cellID = 'CellName',AvgPerCell=True)
renameAxes(ax,renameDict1)
fig.suptitle('aSFL-6FP expressing long linker: H(fluo)')
# jvu.archiveFig(fig, ax, name='aSFL-6FP_iMC_H(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%  Diverse 2D plots

# %%%%%


Filters = [(GlobalTable_ctField['validated'] == True), (GlobalTable_ctField['cell subtype'] == 'aSFL')]
fig, ax = D2Plot(GlobalTable_ctField, XCol='medianThickness',YCol='fluctuAmpli',CondCol = ['cell subtype', 'drug'],                 Filters=Filters, modelFit=True, figSizeFactor = 0.8)
fig.suptitle('3T3aSFL - Dh = f(H)')
# jvu.archiveFig(fig, ax, name='aSFL_Dh(h)_drug', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='ThicknessPlots')
plt.show()


# %%%%%


# Same as above without the 2 thickest cells
Filters = [(GlobalTable_ctField['validated'] == True), (GlobalTable_ctField['cell subtype'] == 'aSFL'),           (GlobalTable_ctField['medianThickness'] <= 700)]
fig, ax = D2Plot(GlobalTable_ctField, XCol='medianThickness',YCol='fluctuAmpli',CondCol = ['cell subtype', 'drug'],                 Filters=Filters, modelFit=True, figSizeFactor = 0.8)
fig.suptitle('3T3aSFL - Dh = f(H)')
# jvu.archiveFig(fig, ax, name='aSFL_Dh(h)_drug_wo2LastPoints', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='ThicknessPlots')
plt.show()


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL')]
fig, ax = D2Plot(GlobalTable_meca, XCol='SurroundingThickness',YCol='EChadwick',CondCol = ['substrate','drug'],           Filters=Filters, cellID = 'CellName', AvgPerCell=True, modelFit=False, modelType='y=A*exp(kx)')
fig.suptitle('3T3aSFL: E(h)')
# jvu.archiveFig(fig, ax, name='aSFL_E(h)_drug&substrate', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='')
plt.show()





# %%%% Smifh2 vs Dmso -- (Sergio's visit september 2021)

# %%%%%


pd.set_option('max_columns', None)
# pd.reset_option('max_columns')
pd.set_option('max_rows', None)
# pd.reset_option('max_rows')


# %%%%%


Outliers = ['21-09-02_M4_P1_C3', '21-09-02_M3_P1_C2', '21-09-02_M3_P1_C9']
Filters = [(GlobalTable_meca_Py['drug'].apply(lambda x : x in ['dmso', 'smifh2', 'dmso. doxycyclin', 'smifh2. doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), (GlobalTable_meca_Py['cell subtype'] == 'aSFL'),            (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
data_filtered = GlobalTable_meca_Py
for fltr in Filters:
    data_filtered = data_filtered.loc[fltr]
data_filtered


# %%%%%


Outliers = ['21-09-02_M4_P1_C3', '21-09-02_M3_P1_C2', '21-09-02_M3_P1_C9']
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['dmso', 'smifh2', 'dmso. doxycyclin', 'smifh2. doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), (GlobalTable_meca_Py['cell subtype'] == 'aSFL'),            (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
co_order = makeOrder(['dmso', 'smifh2', 'dmso. doxycyclin', 'smifh2. doxycyclin'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['drug'],                 Parameters=['ctFieldThickness','EChadwick'], Filters=Filters,                 AvgPerCell=True, cellID='cellID', co_order=co_order, figSizeFactor = 1.8, orientation = 'v')
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL smifh2')
jvu.archiveFig(fig, ax, name='3T3aSFL_smifh2&doxy_Thickness&EChadwick', figSubDir = figSubDir)
plt.show()
# 'ctFieldThickness', 'surroundingThickness', 'none', 'doxycyclin', 


# %%%% Long linker second clone (3T3 aSFL-6FP-2) -- (Sergio's visit september 2021)

# %%%%%


pd.set_option('max_columns', None)
pd.reset_option('max_columns')
pd.set_option('max_rows', None)
pd.reset_option('max_rows')


# %%%%%


# Outliers = []
# Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), \
#            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])), \
#            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL', 'aSFL-6FP', 'aSFL-6FP-2'])), \
#            (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
# co_order = makeOrder(['aSFL', 'aSFL-6FP', 'aSFL-6FP-2'], ['none', 'doxycyclin'])
# box_pairs=[('aSFL & none', 'aSFL & doxycyclin'),
#  ('aSFL-6FP & none', 'aSFL-6FP & doxycyclin'),
#  ('aSFL-6FP-2 & none', 'aSFL-6FP-2 & doxycyclin'),
#  ('aSFL & none', 'aSFL-6FP & none'),
#  ('aSFL & none', 'aSFL-6FP-2 & none'),
#  ('aSFL-6FP & none', 'aSFL-6FP-2 & none')]
# fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype', 'drug'],\
#                  Parameters=['ctFieldThickness','EChadwick'], Filters=Filters,\
#                  AvgPerCell=True, cellID='cellID', co_order=co_order, box_pairs = box_pairs, figSizeFactor = 1, orientation = 'v')
# ax[0].set_ylim([0, ax[0].get_ylim()[1]])
# renameAxes(ax,renameDict1)
# fig.suptitle('3T3aSFL 6FP-2')
# jvu.archiveFig(fig, ax, name='3T3aSFL-6FP-2_Thickness&EChadwick', figSubDir = figSubDir)
# plt.show()
# # 'ctFieldThickness', 'surroundingThickness', 'none', 'doxycyclin', 


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL', 'aSFL-6FP-2'])),            (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
co_order = makeOrder(['aSFL', 'aSFL-6FP-2'], ['none', 'doxycyclin'])
box_pairs=[('aSFL & none', 'aSFL & doxycyclin'),
 ('aSFL-6FP-2 & none', 'aSFL-6FP-2 & doxycyclin'),
 ('aSFL & none', 'aSFL-6FP-2 & none')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype', 'drug'],                 Parameters=['ctFieldThickness','EChadwick'], Filters=Filters,                 AvgPerCell=True, cellID='cellID', co_order=co_order, box_pairs = box_pairs, figSizeFactor = 1, orientation = 'v')
ax[0].set_ylim([0, ax[0].get_ylim()[1]])
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL 6FP-2')
jvu.archiveFig(fig, ax, name='3T3aSFL-6FP-2_Thickness&EChadwick', figSubDir = figSubDir)
plt.show()
# 'ctFieldThickness', 'surroundingThickness', 'none', 'doxycyclin', 


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL-6FP-2'])),            (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype', 'drug'],                  Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=True)

renameAxes(ax,renameDict1)
fig.suptitle('aSFL-6FP-2 expressing linker: E(fluo)')
yt = ax.get_yticks()
ax.set_yticklabels((yt/1000).astype(int))
ax.set_ylabel('E Chadwick (kPa)')

# jvu.archiveFig(fig, ax, name='aSFL_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL-6FP'])), \
           (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype', 'drug'],                  Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=True)

renameAxes(ax,renameDict1)
fig.suptitle('aSFL-6FP expressing linker: E(fluo)')
yt = ax.get_yticks()
ax.set_yticklabels((yt/1000).astype(int))
ax.set_ylabel('E Chadwick (kPa)')

# jvu.archiveFig(fig, ax, name='aSFL_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL'])), \
           (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype', 'drug'],                  Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=True)

renameAxes(ax,renameDict1)
fig.suptitle('aSFL expressing linker: E(fluo)')
yt = ax.get_yticks()
ax.set_yticklabels((yt/1000).astype(int))
ax.set_ylabel('E Chadwick (kPa)')

# jvu.archiveFig(fig, ax, name='aSFL_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%% High expresser experiment, second version of the clone (aSFL A8-2) -- (Sergio's visit september 2021)

# %%%%%


# Outliers = []
# Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),
#            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),
#            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
#            (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL', 'aSFL-A8', 'aSFL-A8-2'])),
#            (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
# co_order = makeOrder(['aSFL', 'aSFL-A8', 'aSFL-A8-2'], ['none', 'doxycyclin'])
# box_pairs=[('aSFL & none', 'aSFL & doxycyclin'),
#  ('aSFL-A8 & none', 'aSFL-A8 & doxycyclin'),
#  ('aSFL-A8-2 & none', 'aSFL-A8-2 & doxycyclin'),
#  ('aSFL & none', 'aSFL-A8 & none'),
#  ('aSFL & none', 'aSFL-A8-2 & none'),
#  ('aSFL-A8 & none', 'aSFL-A8-2 & none')]
# fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype', 'drug'],\
#                  Parameters=['surroundingThickness','EChadwick'], Filters=Filters,\
#                  AvgPerCell=True, cellID='cellID', co_order=co_order, box_pairs = box_pairs, 
#                  figSizeFactor = 1, orientation = 'v', useHue = True)
# ax[0].set_ylim([0, ax[0].get_ylim()[1]])
# renameAxes(ax,renameDict1)
# fig.suptitle('3T3aSFL A8-2')
# # jvu.archiveFig(fig, ax, name='3T3aSFL-A8-2_Thickness&EChadwick', figSubDir = figSubDir)
# plt.show()
# # 'ctFieldThickness', 'surroundingThickness', 'none', 'doxycyclin', 


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),
           (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),
           (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL', 'aSFL-A8-2'])),
           (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
co_order = makeOrder(['aSFL', 'aSFL-A8-2'], ['none', 'doxycyclin'])
box_pairs=[('aSFL & none', 'aSFL & doxycyclin'),
 ('aSFL-A8-2 & none', 'aSFL-A8-2 & doxycyclin'),
 ('aSFL & none', 'aSFL-A8-2 & none')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype', 'drug'],                 Parameters=['surroundingThickness','EChadwick'], Filters=Filters,                 AvgPerCell=True, cellID='cellID', co_order=co_order, box_pairs = box_pairs, 
                 figSizeFactor = 1, orientation = 'v', useHue = False)
ax[0].set_ylim([0, ax[0].get_ylim()[1]])
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL A8-2')
jvu.archiveFig(fig, ax, name='3T3aSFL-A8-2_Thickness&EChadwick', figSubDir = figSubDir)
plt.show()
# 'ctFieldThickness', 'surroundingThickness', 'none', 'doxycyclin', 


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL-A8-2'])), \
           (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype', 'drug'],                  Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=True)

renameAxes(ax,renameDict1)
fig.suptitle('aSFL-A8-2 expressing linker: E(fluo)')
yt = ax.get_yticks()
ax.set_yticklabels((yt/1000).astype(int))
ax.set_ylabel('E Chadwick (kPa)')

jvu.archiveFig(fig, ax, name='aSFL-A8-2_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL'])), \
           (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['cell subtype', 'drug'],                  Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=True)

renameAxes(ax,renameDict1)
fig.suptitle('aSFL expressing linker: E(fluo)')
yt = ax.get_yticks()
ax.set_yticklabels((yt/1000).astype(int))
ax.set_ylabel('E Chadwick (kPa)')

# jvu.archiveFig(fig, ax, name='aSFL_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),            (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),            (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL', 'aSFL-A8-2'])), \
           (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = ['drug'],                  Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=True)

renameAxes(ax,renameDict1)
fig.suptitle('{aSFL & aSFL-A8-2 mixed} expressing linker: E(fluo)')
yt = ax.get_yticks()
ax.set_yticklabels((yt/1000).astype(int))
ax.set_ylabel('E Chadwick (kPa)')

# jvu.archiveFig(fig, ax, name='aSFL_A11+A8-2_iMC_E(fluo)', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='FluoPlots')
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['doxycyclin'])), \
           (pd.isna(GlobalTable_meca_Py['meanFluoPeakAmplitude']) != True),
          (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL','aSFL-6FP','aSFL-A8','aSFL-6FP-2','aSFL-A8-2'])),]
co_order = makeOrder(['aSFL','aSFL-A8','aSFL-6FP','aSFL-A8-2','aSFL-6FP-2'])
box_pairs=[('aSFL', 'aSFL-A8'),
 ('aSFL', 'aSFL-6FP'),
 ('aSFL', 'aSFL-A8-2'),
 ('aSFL', 'aSFL-6FP-2')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype'],Parameters=['meanFluoPeakAmplitude'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,box_pairs=box_pairs,stats=True,figSizeFactor=1.75)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - different cell lines expressing the linker\nLinker expression quantification by fluo')
jvu.archiveFig(fig, ax, name='3T3aSFL_likerType&_fluoExp_NewVersion', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])), \
           (pd.isna(GlobalTable_meca_Py['meanFluoPeakAmplitude']) != True),
          (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL','aSFL-6FP','aSFL-A8','aSFL-6FP-2','aSFL-A8-2'])),]
co_order = makeOrder(['aSFL','aSFL-A8','aSFL-6FP','aSFL-A8-2','aSFL-6FP-2'], ['none', 'doxycyclin'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype', 'drug'],Parameters=['meanFluoPeakAmplitude'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=False,figSizeFactor=1.05)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL different cell lines\nLinker expression quantif by fluo')
# jvu.archiveFig(fig, ax, name='3T3aSFL_likerType&_fluoExp_NewVersion_2', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])), \
           (pd.isna(GlobalTable_meca_Py['meanFluoPeakAmplitude']) != True),
          (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL','aSFL-6FP-2','aSFL-A8-2'])),]
co_order = makeOrder(['aSFL','aSFL-A8-2','aSFL-6FP-2'], ['none', 'doxycyclin'])
box_pairs=[('aSFL & doxycyclin', 'aSFL-A8-2 & doxycyclin'),
 ('aSFL & doxycyclin', 'aSFL-6FP-2 & doxycyclin')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['cell subtype', 'drug'],Parameters=['meanFluoPeakAmplitude'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,box_pairs=box_pairs,stats=True,
                 figSizeFactor=0.95)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - different cell lines\nLinker expression quantification by fluo')
jvu.archiveFig(fig, ax, name='3T3aSFL_likerType&_fluoExp_NewVersion_2', figSubDir = figSubDir)
plt.show()


# %%% MCA project - Py2 > Full python processing

plt.close('all')
figSubDir = 'MCAproject'

# %%%% Entire table (meca_Py2)

# %%%%% Test of dataframe sorting

condition = 'drug'
co_order = ['none', 'doxycyclin']
order_dict = {}
for i in range(len(co_order)):
    order_dict[co_order[i]] = i+1

data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True),
          (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))]

data_f = data
for F in Filters:
    data_f = data_f[F]

data_f_sorted = data_f.sort_values(by=[condition], key=lambda x: x.map(order_dict))
df0 = data_f_sorted


# %%%%% 3T3 aSFL - Jan 21 - Compression experiments


data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['surroundingThickness'] <= 800),
           (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))]

fig, ax = D1Plot(data, CondCol=['drug'], Parameters=['surroundingThickness', 'EChadwick'], Filters=Filters, 
                Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], stats=True, statMethod='Mann-Whitney', 
               box_pairs=[], figSizeFactor = 0.9, markersizeFactor=1, orientation = 'h', AvgPerCell = True)
renameAxes(ax,renameDict1)
ax[0].set_ylim([0, 850])
# ax[0].legend(loc = 'upper right', fontsize = 8)
# ax[1].legend(loc = 'upper right', fontsize = 8)
fig.suptitle('3T3 aSFL - Jan 21 - Compression experiments')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_H&Echad_simple', figSubDir = figSubDir)
plt.show()


# %%%%% 3T3 aSFL - Jan 21 - Detailed compression experiments


data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['surroundingThickness'] <= 800),
           (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))]

fig, ax = D1PlotDetailed(data, CondCol=['drug'], Parameters=['bestH0', 'EChadwick'], Filters=Filters, 
                Boxplot=True, cellID='cellID', co_order = ['none', 'doxycyclin'], stats=True, statMethod='Mann-Whitney', 
                box_pairs=[], figSizeFactor = 1.8, markersizeFactor=1, orientation = 'v', showManips = True)

ax[0].legend(loc = 'upper right', fontsize = 8)
ax[1].legend(loc = 'upper right', fontsize = 8)
ax[0].set_ylim([0, 850])
ax[1].set_ylim([1, 2e5])
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_H&Echad_detailed', figSubDir = figSubDir)
plt.show()


# %%%%% 3T3aSFL on patterns: Constant Field


data = GlobalTable_ctField_Py # Tracking python, postprocessing python
dates = []

Filters = [(data['validated'] == True), (data['medianThickness'] <= 1000)]

co_order = makeOrder(['none','doxycyclin'])
fig, ax = D1Plot(data, CondCol=['drug'],Parameters=['medianThickness','fluctuAmpli'],
                 Filters=Filters,stats=True,co_order=co_order,figSizeFactor=0.9)

renameAxes(ax,renameDict1)
ax[0].set_ylim([0,600])
ax[1].set_ylim([0,600])
fig.suptitle('3T3aSFL on patterns: Constant Field')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_medianThickness', figSubDir = figSubDir)

plt.show()


# %%%%% Nonlinearity, first check

data = GlobalTable_meca_Py2

columns = data.columns.values
minS = []
maxS = []
centerS = []
colList = []

for c in columns:
    if 'KChadwick' in c and '<s<' in c:
        valS = c.split('_')[-1].split('<s<')
        valS[0], valS[1] = int(valS[0]), int(valS[1])
        minS.append(valS[0])
        maxS.append(valS[1])
        centerS.append(int(0.5*(valS[0]+valS[1])))
        colList.append(c)
#     elif 'Npts' in c and '<s<' in c:
        
        
Filters = [(data['validatedFit'] == True),
           (data['validatedThickness'] == True),
           (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))]

F = Filters[0]
for fltr in Filters[1:]:
    F = F & fltr
    
data_f = data[F]

count = []
for i in range(len(colList)):
    c = colList[i]
    S = centerS[i]
    Vals = data_f[c].values
    countValidFit = len(Vals) - np.sum(np.isnan(Vals))
    count.append(countValidFit)
    
d = {'center_stress' : centerS, 'count_valid' : count}
dfCheck = pd.DataFrame(d)


print(data_f.shape[0])
df0 = dfCheck.T
print(dfCheck.T)


# %%%%% Many tangeantial moduli of aSFL 3T3 - control vs iMC linker

data = GlobalTable_meca_Py2

listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

fig, axes = plt.subplots(1, len(listeS), figsize = (20,6))
kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in ['21-09-09']))] #, '21-12-16' , '21-01-21'

    fig, ax = D1Plot(data, fig=fig, ax=axes[kk], CondCol=['drug'], Parameters=['KChadwick_'+interval], Filters=Filters, 
                   Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], stats=True, statMethod='Mann-Whitney', 
                   AvgPerCell = False, box_pairs=[], figSizeFactor = 1, markersizeFactor=0.7, orientation = 'h', 
                   stressBoxPlot=True, bypassLog = True)

#     axes[kk].legend(loc = 'upper right', fontsize = 8)
    label = axes[kk].get_ylabel()
    axes[kk].set_title(label.split('_')[-1] + 'Pa', fontsize = 10)
    # axes[kk].set_yscale('log')
    axes[kk].set_ylim([1e2, 2e4])
    axes[kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
    if kk == 0:
        axes[kk].set_ylabel(label.split('_')[0] + ' (Pa)')
        axes[kk].tick_params(axis='x', labelsize = 10)
    else:
        axes[kk].set_ylabel(None)
        axes[kk].set_yticklabels([])
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Detailed1DPlot', figSubDir = figSubDir)
    
    kk+=1

renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
fig.tight_layout()
fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals', figSubDir = figSubDir)
plt.show()

# %%%%% Many tangeantial moduli of aSFL 3T3 - control vs iMC linker - Plot by plot


data = GlobalTable_meca_Py2

listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)
    textForFileName = str(S-75) + '-' + str(S+75) + 'Pa'

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))] #, '21-12-16'

    fig, ax = D1Plot(data, CondCol=['drug'], Parameters=['bestH0', 'KChadwick_'+interval], Filters=Filters, 
                   Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], stats=True, statMethod='Mann-Whitney', 
                   box_pairs=[], figSizeFactor = 1.0, markersizeFactor=1, orientation = 'h', AvgPerCell = True)

    # ax[0].legend(loc = 'upper right', fontsize = 6)
    # ax[1].legend(loc = 'lower right', fontsize = 6)
    ax[0].set_ylim([0,1200])
    ax[1].set_ylim([1e2, 3e4])
#     jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_H0-K_'+textForFileName, figSubDir = figSubDir)
    plt.show()
    

# %%%%% Many tangeantial moduli of aSFL 3T3 - control vs iMC linker - detailed plot by plot

data = GlobalTable_meca_Py2

listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)
    textForFileName = str(S-75) + '-' + str(S+75) + 'Pa'

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))] #, '21-12-16'

    fig, ax = D1PlotDetailed(data, CondCol=['drug'], Parameters=['bestH0', 'KChadwick_'+interval], Filters=Filters, 
                   Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], stats=True, statMethod='Mann-Whitney', 
                   box_pairs=[], figSizeFactor = 1.0, markersizeFactor=1, orientation = 'h', showManips = True)

    ax[0].legend(loc = 'upper right', fontsize = 6)
    ax[0].set_ylim([0,1200])
    ax[1].legend(loc = 'lower right', fontsize = 6)
    ax[1].set_ylim([1e2, 3e4])
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_H0-K_'+textForFileName, figSubDir = figSubDir)
    plt.show()
    


# %%%%% aSFL 3T3 - control vs iMC linker - extremal stress


# ['21-12-08', '21-12-16', '22-01-12']
data = GlobalTable_meca_Py2

dates = ['21-01-18', '21-01-21'] # ['21-01-18', '21-01-21', '21-12-08']

filterList = [(data['validatedThickness'] == True),
              (data['cell subtype'] == 'aSFL'), 
              (data['bead type'] == 'M450'),
              (data['substrate'] == '20um fibronectin discs'),
              (data['date'].apply(lambda x : x in dates))]  # (data['validatedFit'] == True), 


fig, ax = plt.subplots(1, 1, figsize = (9, 6), tight_layout=True)

globalFilter = filterList[0]
for k in range(1, len(filterList)):
    globalFilter = globalFilter & filterList[k]

data_f = data[globalFilter]
    
condition = 'drug'
conditionValues = data_f[condition].unique()
conditionFilter = {}
for val in conditionValues:
    conditionFilter[val] = (data_f[condition] == val)
    
rD = {'none' : 'Ctrl', 'doxycyclin' : 'iMC linker'}

for val in conditionValues:
    if val == 'none':
        col1, col2, col3 = 'deepskyblue', 'skyblue', 'royalblue'
    if val == 'doxycyclin':
        col1, col2, col3 = 'orangered', 'lightsalmon', 'firebrick'
    X = data_f[conditionFilter[val]]['bestH0'].values
    Ncomp = len(X)
    Smin = data_f[conditionFilter[val]]['minStress'].values
    Smax = data_f[conditionFilter[val]]['maxStress'].values

    for i in range(Ncomp):
        if i == 0:
            textLabel = rD[val]
        else:
            textLabel = None
        ax.plot([X[i], X[i]], [Smin[i], Smax[i]], 
                ls = '-', color = col1, alpha = 0.3,
                zorder = 1, label = textLabel)
        ax.plot([X[i]], [Smin[i]], 
                marker = 'o', markerfacecolor = col2, markeredgecolor = 'k', markersize = 4,
                ls = '')
        ax.plot([X[i]], [Smax[i]], 
                marker = 'o', markerfacecolor = col3, markeredgecolor = 'k', markersize = 4,
                ls = '')

ax.set_xlabel('Cortical thickness (nm)')
ax.set_xlim([0, 1200])
locator = matplotlib.ticker.MultipleLocator(100)
ax.xaxis.set_major_locator(locator)

ax.set_ylabel('Extremal stress values (Pa)')
ax.set_ylim([3e1, 5e4])
ax.set_yscale('log')
# locator = matplotlib.ticker.MultipleLocator(100)
# ax.yaxis.set_major_locator(locator)
ax.yaxis.grid(True)

ax.legend(loc = 'upper right')

ax.set_title('Jan 21 - Extremal stress values vs. thickness, all compressions, ctrl vs. iMC')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_ctrlVsLinker_StressRanges', figSubDir = figSubDir)
print(Ncomp)
plt.show()


# %%%%% aSFL 3T3 - control vs iMC linker -- E(h) -- plot  by plot

data = GlobalTable_meca_Py2

listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

# fig, axes = plt.subplots(len(listeS), 1, figsize = (9,40))
# kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)
    textForFileName = str(S-75) + '-' + str(S+75) + 'Pa'

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))] #, '21-12-16'

    fig, ax = D2Plot_wFit(data, #fig=fig, ax=axes[kk], 
                          XCol='bestH0', YCol='KChadwick_'+interval, CondCol = ['drug'], Filters=Filters, 
                          cellID = 'cellID', AvgPerCell=False, 
                          xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')
    
    # individual figs
    ax.set_xlim([8e1, 1.2e3])
    ax.set_ylim([3e2, 5e4])
    renameAxes(ax,{'none':'ctrl', 'doxycyclin':'iMC'})
    fig.set_size_inches(6,4)
    fig.tight_layout()
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_ctrlVsLinker_E(h)_' + textForFileName, figSubDir = figSubDir)
    
    
    # one big fig
#     axes[kk].set_xlim([8e1, 1.2e3])
#     axes[kk].set_ylim([3e2, 5e4])
#     kk+=1


# renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
# fig.tight_layout()
plt.show()


# %%%%% aSFL 3T3 - control vs iMC linker -- E(h)


data = GlobalTable_meca_Py2

listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

fig, axes = plt.subplots(len(listeS), 1, figsize = (9,40))
kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))] #, '21-12-16'

    fig, ax = D2Plot_wFit(data, fig=fig, ax=axes[kk], 
                          XCol='surroundingThickness', YCol='KChadwick_'+interval, CondCol = ['drug'], Filters=Filters, 
                          cellID = 'cellID', AvgPerCell=True, 
                          xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')
    
    
    
#     axes[kk].legend(loc = 'upper right', fontsize = 8)
#     label = axes[kk].get_ylabel()
#     axes[kk].set_title(label.split('_')[-1] + 'Pa', fontsize = 10)
#     axes[kk].set_yscale('log')
    axes[kk].set_xlim([8e1, 1.2e3])
#     axes[kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
#     if kk == 0:
#         axes[kk].set_ylabel(label.split('_')[0] + ' (Pa)')
#         axes[kk].tick_params(axis='x', labelsize = 10)
#     else:
#         axes[kk].set_ylabel(None)
#         axes[kk].set_yticklabels([])
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Detailed1DPlot', figSubDir = figSubDir)
    
    kk+=1

renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
fig.tight_layout()
# fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')
plt.show()


# %%%%  Table with only MCA data


# %%%%% Check the number of valid compressions for each stress interval


data = GlobalTable_meca_MCA

columns = data.columns.values
minS = []
maxS = []
centerS = []
colList = []

for c in columns:
    if 'KChadwick' in c and '<s<' in c:
        valS = c.split('_')[-1].split('<s<')
        valS[0], valS[1] = int(valS[0]), int(valS[1])
        minS.append(valS[0])
        maxS.append(valS[1])
        centerS.append(int(0.5*(valS[0]+valS[1])))
        colList.append(c)
#     elif 'Npts' in c and '<s<' in c:
        
intervalSizes = np.array(maxS) - np.array(minS)
if len(np.unique(intervalSizes)) == 1:
    intervalSize = int(intervalSizes[0])
        
Filters = [(data['validatedFit'] == True),
           (data['validatedThickness'] == True),
           (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))]

data_f = data
for F in Filters:
    data_f = data_f[F]

count = []
for i in range(len(colList)):
    c = colList[i]
    S = centerS[i]
    Vals = data_f[c].values
    countValidFit = len(Vals) - np.sum(np.isnan(Vals))
    count.append(countValidFit)
    
d = {'center_stress' : centerS, 'count_valid' : count}
dfCheck = pd.DataFrame(d)


print(data_f.shape[0])
print(dfCheck.T)


# %%%%% All stress ranges - K(stress), ctrl vs doxy

# '21-01-18', '21-01-21' > standard aSFL (A11)
# '21-04-27', '21-04-28' > long linker aSFL (6FP)
# '21-09-08' > long linker aSFL (6FP-2)
# '21-09-09' > 'high expresser' aSFL (A8)

data = GlobalTable_meca_MCA
dates = ['21-01-18', '21-01-21']
listeS = [i for i in range(100, 1100, 100)]
# listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]


fig, axes = plt.subplots(1, len(listeS), figsize = (20,6))
kk = 0

for S in listeS:
    interval = str(S-intervalSize//2) + '<s<' + str(S+intervalSize//2)

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in dates))]
    
    
    # axes[kk].set_yscale('log')
    # axes[kk].set_ylim([5e2, 5e4])
    # axes[kk].set_yscale('linear')
    axes[kk].set_ylim([-10, 3e4])

    fig, ax = D1Plot(data, fig=fig, ax=axes[kk], CondCol=['drug'], Parameters=['KChadwick_'+interval], 
                     Filters=Filters, Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], 
                     stats=True, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
                     figSizeFactor = 1, markersizeFactor=0.7, orientation = 'h', stressBoxPlot = True, 
                     bypassLog = True)

    # axes[kk].legend(loc = 'upper right', fontsize = 8)
    label = axes[kk].get_ylabel()
    axes[kk].set_title(label.split('_')[-1] + 'Pa', fontsize = 10)
    axes[kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
    axes[kk].set_ylim([-10, 3e4])
    if kk == 0:
        axes[kk].set_ylabel(label.split('_')[0] + ' (Pa)')
        axes[kk].tick_params(axis='x', labelsize = 10)
    else:
        axes[kk].set_ylabel(None)
        axes[kk].set_yticklabels([])
    
    kk+=1

renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
fig.tight_layout()
# fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')

plt.show()


# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = figSubDir)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = ownCloudTodayFigDir, figSubDir = figSubDir)

# %%%%% Only One Stress Range - K(stress), ctrl vs doxy

# '21-01-18', '21-01-21' > standard aSFL (A11)
# '21-04-27', '21-04-28' > long linker aSFL (6FP)
# '21-09-08' > long linker aSFL (6FP-2)
# '21-09-09' > 'high expresser' aSFL (A8)

data = GlobalTable_meca_MCA
dates = ['21-01-18', '21-01-21']
listeS = [300]
# listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]


fig, axes = plt.subplots(1, len(listeS), figsize = (6,5))
kk = 0

intervalSize = 200

for S in listeS:
    interval = str(S-intervalSize//2) + '<s<' + str(S+intervalSize//2)

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in dates))]
    
    
    # axes.set_yscale('log')
    # axes.set_ylim([5e2, 5e4])
    # axes.set_yscale('linear')
    axes.set_ylim([500, 3e4])

    fig, ax = D1Plot(data, fig=fig, ax=axes, CondCol=['drug'], Parameters=['KChadwick_'+interval], 
                     Filters=Filters, Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], 
                     stats=True, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
                     figSizeFactor = 1, markersizeFactor=1.1, orientation = 'h', stressBoxPlot = True, 
                     bypassLog = False)

    # axes.legend(loc = 'upper right', fontsize = 8)
    label = axes.get_ylabel()
    axes.set_title(label.split('_')[-1] + 'Pa\n', fontsize = 10)
    axes.tick_params(axis='x', labelrotation = 30, labelsize = 14)
    axes.set_ylim([-10, 3e4])
    
renameAxes(axes,{'none':'Control', 'doxycyclin':'Increased MCA', 'KChadwick_'+interval:'Tangeantial Modulus (Pa)'})
fig.tight_layout()
# fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')

plt.show()


# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = figSubDir)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = ownCloudTodayFigDir, figSubDir = figSubDir)




# %%%%% All stress ranges - K(stress), ctrl vs doxy - 200nm < H0 < 300nm

# '21-01-18', '21-01-21' > standard aSFL (A11)
# '21-04-27', '21-04-28' > long linker aSFL (6FP)
# '21-09-08' > long linker aSFL (6FP-2)
# '21-09-09' > 'high expresser' aSFL (A8)

data = GlobalTable_meca_MCA
dates = ['21-01-18', '21-01-21']
listeS = [i for i in range(200, 1100, 100)]
# listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]


fig, axes = plt.subplots(1, len(listeS), figsize = (20,6))
kk = 0

for S in listeS:
    interval = str(S-intervalSize//2) + '<s<' + str(S+intervalSize//2)

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               ((data['bestH0'] <= 300) & (data['bestH0'] >= 200)),
               (data['date'].apply(lambda x : x in dates))]
    
    fig, ax = D1Plot(data, fig=fig, ax=axes[kk], CondCol=['drug'], Parameters=['KChadwick_'+interval], 
                     Filters=Filters, Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], 
                     stats=True, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
                     figSizeFactor = 1, markersizeFactor=0.7, orientation = 'h', stressBoxPlot = True, 
                     bypassLog = True)

#     axes[kk].legend(loc = 'upper right', fontsize = 8)
    label = axes[kk].get_ylabel()
    axes[kk].set_title(label.split('_')[-1] + 'Pa', fontsize = 10)
    axes[kk].set_yscale('log')
    axes[kk].set_ylim([5e2, 5e4])
    axes[kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
    if kk == 0:
        axes[kk].set_ylabel(label.split('_')[0] + ' (Pa)')
        axes[kk].tick_params(axis='x', labelsize = 10)
    else:
        axes[kk].set_ylabel(None)
        axes[kk].set_yticklabels([])
   
    kk+=1

renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
fig.tight_layout()
fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')
plt.show()
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_200-300nm', figSubDir = figSubDir)

# %%%%% All stress ranges - K(stress), ctrl vs doxy -- With bestH0 distributions and E(h)

# '21-01-18', '21-01-21' > standard aSFL (A11)
# '21-04-27', '21-04-28' > long linker aSFL (6FP)
# '21-09-08' > long linker aSFL (6FP-2)
# '21-09-09' > 'high expresser' aSFL (A8)

data = GlobalTable_meca_MCA
dates = ['21-01-18', '21-01-21']
listeS = [i for i in range(100, 1100, 100)]
# listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]


fig, axes = plt.subplots(3, len(listeS), figsize = (20,10))
kk = 0

for S in listeS:
    interval = str(S-intervalSize//2) + '<s<' + str(S+intervalSize//2)
    print(interval)
    
    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['bestH0'] <= 1100),
               (data['date'].apply(lambda x : x in dates))]
    
    
    # axes[0, kk].set_yscale('log')
    axes[0, kk].set_ylim([5e2, 5e4])
    # axes[0, kk].set_yscale('linear')
    # axes[0, kk].set_ylim([0, 3e4])

    D1Plot(data, fig=fig, ax=axes[0, kk], CondCol=['drug'], Parameters=['KChadwick_'+interval], 
         Filters=Filters, Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], 
         stats=True, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
         figSizeFactor = 1, markersizeFactor=0.6, orientation = 'h', 
         stressBoxPlot = True)# , bypassLog = True)
    
    D1Plot(data, fig=fig, ax=axes[1, kk], CondCol=['drug'], Parameters=['bestH0'], 
         Filters=Filters, Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], 
         stats=True, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
         figSizeFactor = 1, markersizeFactor=0.6, orientation = 'h', stressBoxPlot = True)
    
    D2Plot_wFit(data, fig=fig, ax=axes[2, kk], Filters=Filters, 
          XCol='bestH0', YCol='KChadwick_'+interval, CondCol = ['drug'], co_order=['none', 'doxycyclin'], 
          cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', 
          modelFit=True, modelType='y=k*x^a', markersizeFactor = 0.6)
    
    # 0.
    # axes[0, kk].legend(loc = 'upper right', fontsize = 8)
    label0 = axes[0, kk].get_ylabel()
    axes[0, kk].set_title(label0.split('_')[-1] + 'Pa', fontsize = 10)
    axes[0, kk].set_ylim([5e2, 5e4])
    # axes[0, kk].set_ylim([1, 3e4])
    # axes[0, kk].set_xticklabels([])
    
    # 1.
    axes[1, kk].set_ylim([0, 1500])
    axes[1, kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
    
    # 2.
    # label2 = axes[2, kk].get_ylabel()
    axes[2, kk].set_xlim([8e1, 1.2e3])
    axes[2, kk].set_ylim([3e2, 5e4])
    axes[2, kk].tick_params(axis='x', labelrotation = 0, labelsize = 10)
    axes[2, kk].xaxis.label.set_size(10)
    axes[2, kk].legend().set_visible(False)
    
    for ll in range(axes.shape[0]):
        axes[ll, kk].tick_params(axis='y', labelsize = 10)
    
    if kk == 0:
        axes[0, kk].set_ylabel(label0.split('_')[0] + ' (Pa)')
        axes[1, kk].set_ylabel('Best H0 (nm)')
        axes[2, kk].set_ylabel('K_Chadwick (Pa)')
    else:
        for ll in range(axes.shape[0]):
            axes[ll, kk].set_ylabel(None)
            axes[ll, kk].set_yticklabels([])
    
    kk+=1

renameAxes(axes.flatten(),{'none':'ctrl', 'doxycyclin':'iMC', 'bestH0' : 'Best H0 (nm)'})
# fig.tight_layout()
# fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')

plt.show()


# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = figSubDir)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = ownCloudTodayFigDir, figSubDir = figSubDir)

# %%%%% All stress ranges - K(stress), ctrl vs doxy -- With bestH0 distributions and E(h) - effect of windows width


data = GlobalTable_meca_MCA
dates = ['21-01-18', '21-01-21']

fitC =  np.array([S for S in range(100, 1000, 100)])
fitW = [100, 150, 200, 250, 300]
fitCenters = np.array([[int(S) for S in fitC] for w in fitW])
fitWidth = np.array([[int(w) for S in fitC] for w in fitW])
fitMin = np.array([[int(S-(w/2)) for S in fitC] for w in fitW])
fitMax = np.array([[int(S+(w/2)) for S in fitC] for w in fitW])
# fitCenters, fitWidth = fitCenters[fitMin>0], fitWidth[fitMin>0]
# fitMin, fitMax = fitMin[fitMin>0], fitMax[fitMin>0]

for ii in range(len(fitW)):
    width = fitW[ii]
    fitMin_ii = fitMin[ii]
    N = len(fitMin_ii[fitMin_ii > 0])
    
    fig, axes = plt.subplots(3, N, figsize = (20,10))
    kk = 0
    
    for S in fitC:
        minS, maxS = int(S-width//2), int(S+width//2)
        interval = str(minS) + '<s<' + str(maxS)
        
        if minS > 0:
            print(interval)
            try:
                Filters = [(data['validatedFit_'+interval] == True),
                           (data['validatedThickness'] == True),
                           (data['bestH0'] <= 1100),
                           (data['date'].apply(lambda x : x in dates))]
                
                
                # axes[0, kk].set_yscale('log')
                axes[0, kk].set_ylim([5e2, 5e4])
                # axes[0, kk].set_yscale('linear')
                # axes[0, kk].set_ylim([0, 3e4])
            
                D1Plot(data, fig=fig, ax=axes[0, kk], CondCol=['drug'], Parameters=['KChadwick_'+interval], 
                     Filters=Filters, Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], 
                     stats=True, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
                     figSizeFactor = 1, markersizeFactor=0.6, orientation = 'h', 
                     stressBoxPlot = True)# , bypassLog = True)
                
                D1Plot(data, fig=fig, ax=axes[1, kk], CondCol=['drug'], Parameters=['bestH0'], 
                     Filters=Filters, Boxplot=True, cellID='cellID',  co_order=['none', 'doxycyclin'], 
                     stats=True, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
                     figSizeFactor = 1, markersizeFactor=0.6, orientation = 'h', stressBoxPlot = True)
                
                D2Plot_wFit(data, fig=fig, ax=axes[2, kk], Filters=Filters, 
                      XCol='bestH0', YCol='KChadwick_'+interval, CondCol = ['drug'],  co_order=['none', 'doxycyclin'], 
                      cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', 
                      modelFit=True, modelType='y=k*x^a', markersizeFactor = 0.6)
                
            except:
                pass
            
            # 0.
            # axes[0, kk].legend(loc = 'upper right', fontsize = 8)
            # label0 = axes[0, kk].get_ylabel()
            axes[0, kk].set_title('S = {:.0f} +/- {:.0f} Pa'.format(S, width//2), fontsize = 10)
            axes[0, kk].set_ylim([5e2, 5e4])
            # axes[0, kk].set_ylim([1, 3e4])
            # axes[0, kk].set_xticklabels([])
            
            # 1.
            axes[1, kk].set_ylim([0, 1500])
            axes[1, kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
            
            # 2.
            # label2 = axes[2, kk].get_ylabel()
            axes[2, kk].set_xlim([8e1, 1.2e3])
            axes[2, kk].set_ylim([3e2, 5e4])
            axes[2, kk].tick_params(axis='x', labelrotation = 0, labelsize = 10)
            axes[2, kk].xaxis.label.set_size(10)
            axes[2, kk].legend().set_visible(False)
            
            for ll in range(axes.shape[0]):
                axes[ll, kk].tick_params(axis='y', labelsize = 10)
            
            if kk == 0:
                axes[0, kk].set_ylabel(label0.split('_')[0] + ' (Pa)')
                axes[1, kk].set_ylabel('Best H0 (nm)')
                axes[2, kk].set_ylabel('K_Chadwick (Pa)')
            else:
                for ll in range(axes.shape[0]):
                    axes[ll, kk].set_ylabel(None)
                    axes[ll, kk].set_yticklabels([])
        
            kk+=1
            
        else:
            pass
    
    renameAxes(axes.flatten(),{'none':'ctrl', 'doxycyclin':'iMC', 'bestH0' : 'Best H0 (nm)'})
    # fig.tight_layout()
    # fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')
    
    fig.suptitle('Width = {:.0f}'.format(width))
    
    plt.show()


    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = figSubDir)
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = ownCloudTodayFigDir, figSubDir = figSubDir)

# %%%%% All stress ranges - K(stress), ctrl vs doxy -- effect of windows width -> 5width


data = GlobalTable_meca_MCA # Get the table in the folder "Width" or recompute it !
dates = ['21-01-18', '21-01-21']

fitC =  np.array([S for S in range(100, 1000, 100)])
fitW = [100, 150, 200, 250, 300]
fitCenters = np.array([[int(S) for S in fitC] for w in fitW])
fitWidth = np.array([[int(w) for S in fitC] for w in fitW])
fitMin = np.array([[int(S-(w/2)) for S in fitC] for w in fitW])
fitMax = np.array([[int(S+(w/2)) for S in fitC] for w in fitW])
# fitCenters, fitWidth = fitCenters[fitMin>0], fitWidth[fitMin>0]
# fitMin, fitMax = fitMin[fitMin>0], fitMax[fitMin>0]

N = len(fitC)
fig, axes = plt.subplots(5, N, figsize = (20,10))

for ii in range(len(fitW)):
    width = fitW[ii]
    fitMin_ii = fitMin[ii]
    
    kk = 0
    
    for S in fitC:
        minS, maxS = int(S-width//2), int(S+width//2)
        interval = str(minS) + '<s<' + str(maxS)
        
        if minS > 0:
            print(interval)

            Filters = [(data['validatedFit_'+interval] == True),
                       (data['validatedThickness'] == True),
                       (data['bestH0'] <= 1000),
                       (data['date'].apply(lambda x : x in dates))]
            
            
            # axes[0, kk].set_yscale('log')
            axes[0, kk].set_ylim([5e2, 5e4])
            # axes[0, kk].set_yscale('linear')
            # axes[0, kk].set_ylim([0, 3e4])
        
            D1Plot(data, fig=fig, ax=axes[ii, kk], CondCol=['drug'], Parameters=['KChadwick_'+interval], 
                 Filters=Filters, Boxplot=True, cellID='cellID', co_order=['none', 'doxycyclin'], 
                 stats=True, statMethod='Mann-Whitney', AvgPerCell = True, box_pairs=[], 
                 figSizeFactor = 1, markersizeFactor=0.6, orientation = 'h', 
                 stressBoxPlot = True)# , bypassLog = True)
            
            # D1Plot(data, fig=fig, ax=axes[1, kk], CondCol=['drug'], Parameters=['bestH0'], 
            #      Filters=Filters, Boxplot=True, cellID='cellID',  co_order=['none', 'doxycyclin'], 
            #      stats=True, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
            #      figSizeFactor = 1, markersizeFactor=0.6, orientation = 'h', stressBoxPlot = True)
            
            # D2Plot_wFit(data, fig=fig, ax=axes[2, kk], Filters=Filters, 
            #       XCol='bestH0', YCol='KChadwick_'+interval, CondCol = ['drug'],  co_order=['none', 'doxycyclin'], 
            #       cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', 
            #       modelFit=True, modelType='y=k*x^a', markersizeFactor = 0.6)
                

            
            # # 0.
            # # axes[0, kk].legend(loc = 'upper right', fontsize = 8)
            # # label0 = axes[0, kk].get_ylabel()
            # axes[0, kk].set_title('S = {:.0f} +/- {:.0f} Pa'.format(S, width//2), fontsize = 10)
            # axes[0, kk].set_ylim([5e2, 5e4])
            # # axes[0, kk].set_ylim([1, 3e4])
            # # axes[0, kk].set_xticklabels([])
            
            # # 1.
            # axes[1, kk].set_ylim([0, 1500])
            # axes[1, kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
            
            # # 2.
            # # label2 = axes[2, kk].get_ylabel()
            # axes[2, kk].set_xlim([8e1, 1.2e3])
            # axes[2, kk].set_ylim([3e2, 5e4])
            # axes[2, kk].tick_params(axis='x', labelrotation = 0, labelsize = 10)
            # axes[2, kk].xaxis.label.set_size(10)
            # axes[2, kk].legend().set_visible(False)  
        else:
            pass
        
        
        axes[ii, kk].set_title('S = {:.0f} +/- {:.0f} Pa'.format(S, width//2), fontsize = 10)
        axes[ii, kk].set_ylim([5e2, 5e4])
        axes[ii, kk].tick_params(axis='y', labelsize = 10)
        
        if kk == 0:
            axes[ii, kk].set_ylabel('K ; width={:.0f}Pa'.format(width), fontsize = 10)
        else:
            axes[ii, kk].set_ylabel(None)
            axes[ii, kk].set_yticklabels([])
            
        if ii == len(fitW)-1:
            axes[ii, kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
        else:
            axes[ii, kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
            axes[ii, kk].set_xlabel(None)
            axes[ii, kk].set_xticklabels([])
            
        kk+=1
    
    renameAxes(axes.flatten(),{'none':'ctrl', 'doxycyclin':'iMC', 'bestH0' : 'Best H0 (nm)'})
    # fig.tight_layout()
    # fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')
    
    fig.suptitle('Tangeantial moduli with different windows for the fit')
    
    plt.show()


    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = figSubDir)
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = ownCloudTodayFigDir, figSubDir = figSubDir)



# %%%%% extremalStress = f(bestH0) - whole dataset


# ['21-12-08', '21-12-16', '22-01-12']
data = GlobalTable_meca_Py2

dates = ['21-01-18', '21-01-21'] # ['21-01-18', '21-01-21', '21-12-08']

filterList = [(data['validatedThickness'] == True),
              (data['cell subtype'] == 'aSFL'), 
              (data['bead type'] == 'M450'),
              (data['substrate'] == '20um fibronectin discs'),
              (data['date'].apply(lambda x : x in dates))]  # (data['validatedFit'] == True), 


fig, ax = plt.subplots(1, 1, figsize = (9, 6), tight_layout=True)

globalFilter = filterList[0]
for k in range(1, len(filterList)):
    globalFilter = globalFilter & filterList[k]

data_f = data[globalFilter]
    
condition = 'drug'
conditionValues = data_f[condition].unique()
conditionFilter = {}
for val in conditionValues:
    conditionFilter[val] = (data_f[condition] == val)
    
rD = {'none' : 'Ctrl', 'doxycyclin' : 'iMC linker'}

for val in conditionValues:
    if val == 'none':
        col1, col2, col3 = 'deepskyblue', 'skyblue', 'royalblue'
    if val == 'doxycyclin':
        col1, col2, col3 = 'orangered', 'lightsalmon', 'firebrick'
    X = data_f[conditionFilter[val]]['bestH0'].values
    Ncomp = len(X)
    Smin = data_f[conditionFilter[val]]['minStress'].values
    Smax = data_f[conditionFilter[val]]['maxStress'].values

    for i in range(Ncomp):
        if i == 0:
            textLabel = rD[val]
        else:
            textLabel = None
        ax.plot([X[i], X[i]], [Smin[i], Smax[i]], 
                ls = '-', color = col1, alpha = 0.3,
                zorder = 1, label = textLabel)
        ax.plot([X[i]], [Smin[i]], 
                marker = 'o', markerfacecolor = col2, markeredgecolor = 'k', markersize = 4,
                ls = '')
        ax.plot([X[i]], [Smax[i]], 
                marker = 'o', markerfacecolor = col3, markeredgecolor = 'k', markersize = 4,
                ls = '')

ax.set_xlabel('Cortical thickness (nm)')
ax.set_xlim([0, 1200])
locator = matplotlib.ticker.MultipleLocator(100)
ax.xaxis.set_major_locator(locator)

ax.set_ylabel('Extremal stress values (Pa)')
ax.set_ylim([3e1, 5e4])
ax.set_yscale('log')
# locator = matplotlib.ticker.MultipleLocator(100)
# ax.yaxis.set_major_locator(locator)
ax.yaxis.grid(True)

ax.legend(loc = 'upper right')

ax.set_title('Jan 21 - Extremal stress values vs. thickness, all compressions, ctrl vs. iMC')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_ctrlVsLinker_StressRanges', figSubDir = figSubDir)
print(Ncomp)
plt.show()


# %%%%% extremalStress = f(bestH0) - just a slice


# ['21-12-08', '21-12-16', '22-01-12']
data = GlobalTable_meca_MCA

dates = ['21-01-18', '21-01-21'] # ['21-01-18', '21-01-21', '21-12-08']

filterList = [(data['validatedThickness'] == True),
              (data['cell subtype'] == 'aSFL'), 
              (data['bead type'] == 'M450'),
              ((data['bestH0'] <= 300) & (data['bestH0'] >= 200)),
              (data['substrate'] == '20um fibronectin discs'),
              (data['date'].apply(lambda x : x in dates))]  # (data['validatedFit'] == True), 


fig, ax = plt.subplots(1, 1, figsize = (9, 6), tight_layout=True)

globalFilter = filterList[0]
for k in range(1, len(filterList)):
    globalFilter = globalFilter & filterList[k]

data_f = data[globalFilter]
    
condition = 'drug'
conditionValues = data_f[condition].unique()
conditionFilter = {}
for val in conditionValues:
    conditionFilter[val] = (data_f[condition] == val)
    
rD = {'none' : 'Ctrl', 'doxycyclin' : 'iMC linker'}

for val in conditionValues:
    if val == 'none':
        col1, col2, col3 = 'deepskyblue', 'skyblue', 'royalblue'
    if val == 'doxycyclin':
        col1, col2, col3 = 'orangered', 'lightsalmon', 'firebrick'
    X = data_f[conditionFilter[val]]['bestH0'].values
    Ncomp = len(X)
    Smin = data_f[conditionFilter[val]]['minStress'].values
    Smax = data_f[conditionFilter[val]]['maxStress'].values

    for i in range(Ncomp):
        if i == 0:
            textLabel = rD[val]
        else:
            textLabel = None
        ax.plot([X[i], X[i]], [Smin[i], Smax[i]], 
                ls = '-', color = col1, alpha = 0.3,
                zorder = 1, label = textLabel)
        ax.plot([X[i]], [Smin[i]], 
                marker = 'o', markerfacecolor = col2, markeredgecolor = 'k', markersize = 4,
                ls = '')
        ax.plot([X[i]], [Smax[i]], 
                marker = 'o', markerfacecolor = col3, markeredgecolor = 'k', markersize = 4,
                ls = '')

ax.set_xlabel('Cortical thickness (nm)')
ax.set_xlim([0, 1200])
locator = matplotlib.ticker.MultipleLocator(100)
ax.xaxis.set_major_locator(locator)

ax.set_ylabel('Extremal stress values (Pa)')
ax.set_ylim([3e1, 5e4])
ax.set_yscale('log')
# locator = matplotlib.ticker.MultipleLocator(100)
# ax.yaxis.set_major_locator(locator)
ax.yaxis.grid(True)

ax.legend(loc = 'upper right')

ax.set_title('Jan 21 - Extremal stress values vs. thickness, all compressions, ctrl vs. iMC')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_ctrlVsLinker_StressRanges', figSubDir = figSubDir)
print(Ncomp)
plt.show()


# %%%%% "Three metrics" : Median H, Best H0 avg per cell, Best H0 all comp // change the dates !

# '21-01-18', '21-01-21' > standard aSFL (A11)
# '21-04-27', '21-04-28' > long linker aSFL (6FP)
# '21-09-08' > long linker aSFL (6FP-2)
# '21-09-09' > 'high expresser' aSFL (A8)


data = GlobalTable_meca_Py2 # Tracking python, postprocessing python
# data['bestH02'] = data['bestH0']**2
dates = ['21-01-18', '21-01-21']

fig, ax = plt.subplots(1, 3, figsize = (12,5))

#1
Filters = [(data['validatedThickness'] == True), 
            (data['ctFieldThickness'] <= 1000),
            (data['bead type'] == 'M450'),
            (data['substrate'] == '20um fibronectin discs'),
            (data['date'].apply(lambda x : x in dates))]

co_order = makeOrder(['none','doxycyclin'])

print(data[data['date'].apply(lambda x : x in dates)]['cell subtype'].unique())

fig, ax1 = D1Plot(data, fig=fig, ax=ax[0], CondCol=['drug'],Parameters=['ctFieldThickness'], AvgPerCell = True,
                 Filters=Filters,stats=True,co_order=co_order,figSizeFactor=1.2)
ax[0].set_title('Median thickness at nominal field B~5.5mT')
renameAxes(ax[0],renameDict1)


#2
Filters = [(data['validatedThickness'] == True),
            (data['bestH0'] <= 1000),
            (data['bead type'] == 'M450'),
            (data['substrate'] == '20um fibronectin discs'),
            (data['date'].apply(lambda x : x in dates))]

co_order = makeOrder(['none','doxycyclin'])

fig, ax2 = D1Plot(data, fig=fig, ax=ax[1], CondCol=['drug'],Parameters=['bestH0'], AvgPerCell = True,
                 Filters=Filters,stats=True,co_order=co_order,figSizeFactor=1.5)

ax[1].set_title('H0, average per cell')
renameAxes(ax[1],renameDict1)
renameAxes(ax[1],{'bestH0':'H0 (nm)'})

#3
# Filters = [(data['validatedThickness'] == True),
#             (data['bestH0'] <= 1000),
#             (data['bead type'] == 'M450'),
#             (data['substrate'] == '20um fibronectin discs'),
#             (data['date'].apply(lambda x : x in dates))]

co_order = makeOrder(['none','doxycyclin'])

fig, ax3 = D1Plot(data, fig=fig, ax=ax[2], CondCol=['drug'],Parameters=['bestH0'], AvgPerCell = False,
                 Filters=Filters,stats=True,co_order=co_order,figSizeFactor=1.5, stressBoxPlot = True, 
                 markersizeFactor = 0.5)

ax[2].set_title('H0, all compressions')
renameAxes(ax[2],renameDict1)
renameAxes(ax[2],{'bestH0':'Best H0 (nm)'})


multiAxes = ax #, thisAx7bis]
                        
for axis in multiAxes:
    axis.set_ylim([0,1050])
    for item in ([axis.title, axis.xaxis.label,axis.yaxis.label] + axis.get_xticklabels() + axis.get_yticklabels()):
        item.set_fontsize(11)

fig.suptitle('3T3aSFL on patterns - Several thickness metrics', fontsize = 14)

# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_3metricsThicknessFromMECA', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = 'ProjectAMC')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_3metricsThicknessFromMECA', figDir = ownCloudTodayFigDir, figSubDir = 'ProjectAMC')
# jvu.archiveFig(fig, ax, todayFigDirLocal + '//MCAproject', name='ThreeThickMetrics', dpi = 100)


plt.show()



# %%%%% Best H0 from fits B~5mT - All comp


data = GlobalTable_meca_MCA # Tracking python, postprocessing python
dates = []

Filters = [(data['validatedThickness'] == True), (data['bestH0'] <= 900)]

fig, ax = plt.subplots(1, 1, figsize = (6,5))

co_order = makeOrder(['none','doxycyclin'])
fig, ax = D1Plot(fig=fig, ax=ax, data=data, CondCol=['drug'],Parameters=['bestH0'], AvgPerCell = False,
                 Filters=Filters,stats=True, co_order=co_order, figSizeFactor=1.0, 
                 markersizeFactor = 0.7, stressBoxPlot = True)


ax[0].set_ylim([0, 1000])
rD = {'none':'Control', 'doxycyclin':'Increased MCA', 'bestH0' : 'Thickness (nm)'}
renameAxes(ax,rD)
# ax[0].set_ylim([0,600])
# ax[1].set_ylim([0,600])
fig.suptitle('3T3aSFL on patterns\nBest H0 from fits\nB~5mT', fontsize = 12)
fig.tight_layout()
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_medianThickness', figSubDir = figSubDir)

plt.show()


# %%%%% Best H0 from fits B~5mT - Avg per cell


data = GlobalTable_meca_MCA # Tracking python, postprocessing python
dates = []

Filters = [(data['validatedThickness'] == True), (data['bestH0'] <= 1000)]

co_order = makeOrder(['none','doxycyclin'])
fig, ax = D1Plot(data, CondCol=['drug'],Parameters=['bestH0'], AvgPerCell = True,
                 Filters=Filters,stats=True,co_order=co_order,figSizeFactor=1.5)
renameAxes(ax,renameDict1)
# ax[0].set_ylim([0,600])
# ax[1].set_ylim([0,600])
fig.suptitle('3T3aSFL on patterns\nBest H0 from fits\nB~5mT', fontsize = 12)
fig.tight_layout()
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_medianThickness', figSubDir = figSubDir)

plt.show()


# %%%%% E(h) for ranges of stress - plot by plot

data = GlobalTable_meca_Py2

listeS = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
dates = ['21-01-18', '21-01-21']
# fig, axes = plt.subplots(len(listeS), 1, figsize = (9,40))
# kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)
    textForFileName = str(S-75) + '-' + str(S+75) + 'Pa'

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in dates))] #, '21-12-16'

    fig, ax = D2Plot_wFit(data, #fig=fig, ax=axes[kk], 
                          XCol='surroundingThickness', YCol='KChadwick_'+interval, CondCol = ['drug'], Filters=Filters, 
                          cellID = 'cellID', AvgPerCell=False, 
                          xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')
    
    # individual figs
    ax.set_xlim([8e1, 1.2e3])
    ax.set_ylim([3e2, 5e4])
    renameAxes(ax,{'none':'ctrl', 'doxycyclin':'iMC'})
    fig.set_size_inches(6,4)
    fig.tight_layout()
#     jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_ctrlVsLinker_E(h)_' + textForFileName, figSubDir = figSubDir)
    
    
    # one big fig
#     axes[kk].set_xlim([8e1, 1.2e3])
#     axes[kk].set_ylim([3e2, 5e4])
#     kk+=1


# renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
# fig.tight_layout()
plt.show()
    

# %%%%% E(h) for ranges of stress - all in one


data = GlobalTable_meca_Py2

listeS = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

fig, axes = plt.subplots(len(listeS), 1, figsize = (9,40))
kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))] #, '21-12-16'

    fig, ax = D2Plot_wFit(data, fig=fig, ax=axes[kk], 
                          XCol='surroundingThickness', YCol='KChadwick_'+interval, CondCol = ['drug'], Filters=Filters, 
                          cellID = 'cellID', AvgPerCell=True, 
                          xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')
    
    
    
#     axes[kk].legend(loc = 'upper right', fontsize = 8)
#     label = axes[kk].get_ylabel()
#     axes[kk].set_title(label.split('_')[-1] + 'Pa', fontsize = 10)
#     axes[kk].set_yscale('log')
    axes[kk].set_xlim([8e1, 1.2e3])
#     axes[kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
#     if kk == 0:
#         axes[kk].set_ylabel(label.split('_')[0] + ' (Pa)')
#         axes[kk].tick_params(axis='x', labelsize = 10)
#     else:
#         axes[kk].set_ylabel(None)
#         axes[kk].set_yticklabels([])
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Detailed1DPlot', figSubDir = figSubDir)
    
    kk+=1

renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
fig.tight_layout()
# fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')
plt.show()
    

# %%%%% Effective stiffness - many ranges of s to find K


data = GlobalTable_meca_MCA # Tracking python, postprocessing python
dates = ['21-01-18', '21-01-21']
figSubDir='MCAproject'

# data['EqStiffness_200Pa'] = (2/3) * (data['bestH0']/1000)**2 * data['KChadwick_75<s<325'] / 10
# data['EqStiffness_300Pa'] = (2/3) * (data['bestH0']/1000)**2 * data['KChadwick_175<s<425'] / 10
# data['EqStiffness_400Pa'] = (2/3) * (data['bestH0']/1000)**2 * data['KChadwick_275<s<525'] / 10
# data['EqStiffness_500Pa'] = (2/3) * (data['bestH0']/1000)**2 * data['KChadwick_375<s<625'] / 10

R = 10
data['bestH0^2'] = (data['bestH0']/1000)**2
listS = [300, 400, 500, 600]
fig, ax = plt.subplots(len(listS), 3, figsize = (12,4.5*len(listS)))

for i in range(len(listS)):
    S = listS[i]
    data['EqStiffness_'+str(S)+'Pa'] = (2/3)*(data['bestH0']/1000)**2 * data['KChadwick_'+str(S-75)+'<s<'+str(S+75)] / R
    
    
    Filters = [(data['validatedThickness'] == True), (data['bestH0'] <= 1000), 
               (data['validatedFit_'+str(S-75)+'<s<'+str(S+75)] == True),
               (data['date'].apply(lambda x : x in dates))]

    co_order = makeOrder(['none','doxycyclin'])
    
    fig, ax00 = D1Plot(data, fig=fig, ax=ax[i,0], CondCol=['drug'],Parameters=['bestH0^2'], 
                       AvgPerCell = False, markersizeFactor = 0.7,
                       Filters=Filters,stats=True,co_order=co_order,figSizeFactor=1.5, stressBoxPlot = True)
    ax[i,0].set_title('bestH0^2')
    renameAxes(ax[i,0],renameDict1)
    renameAxes(ax[i,0],{'bestH0^2':'H0 squared (nm²)', 
                        'EqStiffness_'+str(S)+'Pa':'k_eff at '+str(S)+'Pa (pN/µm)'})
    
    fig, ax01 = D1Plot(data, fig=fig, ax=ax[i,1], CondCol=['drug'],Parameters=['EqStiffness_'+str(S)+'Pa'], 
                       AvgPerCell = False, markersizeFactor = 0.7,
                       Filters=Filters,stats=True,co_order=co_order,figSizeFactor=1.5, stressBoxPlot = True)
    ax[i,1].set_title('Effective spring constant at '+str(S)+'Pa, all compressions')
    renameAxes(ax[i,1],renameDict1)
    renameAxes(ax[i,1],{'bestH0':'Best H0 from fit (nm)', 
                        'EqStiffness_'+str(S)+'Pa':'k_eff at '+str(S)+'Pa (pN/µm)'})


    fig, ax02 = D1Plot(data, fig=fig, ax=ax[i,2], CondCol=['drug'],Parameters=['EqStiffness_'+str(S)+'Pa'], 
                       AvgPerCell = True, markersizeFactor = 0.7,
                       Filters=Filters,stats=True,co_order=co_order,figSizeFactor=1.5)
#     ax[i,2].set_title('Effective spring constant at '+str(S)+'Pa, average per cell')
    renameAxes(ax[i,2],renameDict1)
    renameAxes(ax[i,2],{'bestH0':'Best H0 from fit (nm)', 
                        'EqStiffness_'+str(S)+'Pa':'k_eff at '+str(S)+'Pa (pN/µm)'})



multiAxes = ax.flatten() #, thisAx7bis]
                        
for axis in multiAxes:
#     axis.set_ylim([0,1050])
    for item in ([axis.title, axis.xaxis.label,axis.yaxis.label] + axis.get_xticklabels() + axis.get_yticklabels()):
        item.set_fontsize(11)

fig.suptitle('3T3aSFL on patterns - K_eff', fontsize = 14)

# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_effectiveShellSpring', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = 'ProjectAMC')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_effectiveShellSpring', figDir = ownCloudTodayFigDir, figSubDir = 'ProjectAMC')


plt.show()

# %%%%% Effective stiffness - K(s=400Pa)


data = GlobalTable_meca_MCA # Tracking python, postprocessing python
dates = ['21-01-18', '21-01-21']
figSubDir='MCAproject'

# data['EqStiffness_200Pa'] = (2/3) * (data['bestH0']/1000)**2 * data['KChadwick_75<s<325'] / 10
# data['EqStiffness_300Pa'] = (2/3) * (data['bestH0']/1000)**2 * data['KChadwick_175<s<425'] / 10
# data['EqStiffness_400Pa'] = (2/3) * (data['bestH0']/1000)**2 * data['KChadwick_275<s<525'] / 10
# data['EqStiffness_500Pa'] = (2/3) * (data['bestH0']/1000)**2 * data['KChadwick_375<s<625'] / 10

R = 10
data['bestH0^2'] = (data['bestH0']/1000)**2
listS = [600]
width = 200

fig, ax = plt.subplots(len(listS), 2, figsize = (12,4.5*len(listS)))

for i in range(len(listS)):
    S = listS[i]
    interval = str(int(S-width//2))+'<s<'+str(int(S+width//2))
    data['EqStiffness_'+str(S)+'Pa'] = (2/3)*(data['bestH0']/1000)**2 * data['KChadwick_' + interval] / R
    
    
    Filters = [(data['validatedThickness'] == True), (data['bestH0'] <= 800), 
               (data['EqStiffness_'+str(S)+'Pa']<= 150),
               (data['validatedFit_' + interval] == True),
               (data['date'].apply(lambda x : x in dates))]

    co_order = makeOrder(['none','doxycyclin'])
    
    ax[0].set_ylim([0, 150])
    fig, ax01 = D1Plot(data, fig=fig, ax=ax[0], CondCol=['drug'],Parameters=['EqStiffness_'+str(S)+'Pa'], 
                       AvgPerCell = False, markersizeFactor = 0.7,
                       Filters=Filters,stats=True,co_order=co_order,figSizeFactor=0.8, stressBoxPlot = True)
    ax[0].set_title('Effective spring constant at '+str(S)+'Pa,\nall compressions')
    renameAxes(ax[0],renameDict1)
    renameAxes(ax[0],{'bestH0':'Best H0 from fit (nm)', 
                        'EqStiffness_'+str(S)+'Pa':'k_eff (pN/µm)'})
    ax[0].set_ylim([0, 150])
    

    ax[1].set_ylim([0, 150])
    fig, ax02 = D1Plot(data, fig=fig, ax=ax[1], CondCol=['drug'],Parameters=['EqStiffness_'+str(S)+'Pa'], 
                       AvgPerCell = True, markersizeFactor = 1.0,
                       Filters=Filters,stats=True,co_order=co_order,figSizeFactor=0.8)
    ax[1].set_title('Effective spring constant at '+str(S)+'Pa,\naverage per cell')
    renameAxes(ax[1],renameDict1)
    renameAxes(ax[1],{'bestH0':'Best H0 from fit (nm)', 
                        'EqStiffness_'+str(S)+'Pa':'k_eff (pN/µm)'})
    ax[1].set_ylim([0, 150])
    



multiAxes = ax.flatten() #, thisAx7bis]
                        
for axis in multiAxes:
#     axis.set_ylim([0,1050])
    for item in ([axis.title, axis.xaxis.label,axis.yaxis.label] + axis.get_xticklabels() + axis.get_yticklabels()):
        item.set_fontsize(11)

fig.suptitle('3T3aSFL on patterns - k_eff', fontsize = 14)

# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_effectiveShellSpring', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = 'ProjectAMC')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_effectiveShellSpring', figDir = ownCloudTodayFigDir, figSubDir = 'ProjectAMC')


plt.show()

# %%%% MCA + nonLin curves

# %%%%% K(s) for 3T3aSFL ctrl vs doxy

#### making the Dataframe 

data = GlobalTable_meca_MCA

dates = ['21-01-18', '21-01-21'] #['21-12-08', '22-01-12'] ['21-01-18', '21-01-21', '21-12-08']

filterList = [(data['validatedThickness'] == True),
              (data['cell subtype'] == 'aSFL'), 
              (data['bead type'] == 'M450'),
              (data['substrate'] == '20um fibronectin discs'),
              (data['bestH0'] <= 1000),
              (data['date'].apply(lambda x : x in dates))]  # (data['validatedFit'] == True),

globalFilter = filterList[0]
for k in range(1, len(filterList)):
    globalFilter = globalFilter & filterList[k]

data_f = data[globalFilter]

width = 200 # 200
fitCenters =  np.array([S for S in range(200, 1000, 50)])
fitMin = np.array([int(S-(width/2)) for S in fitCenters])
fitMax = np.array([int(S+(width/2)) for S in fitCenters])

regionFitsNames = [str(fitMin[ii]) + '<s<' + str(fitMax[ii]) for ii in range(len(fitMin))]

listColumnsMeca = []

KChadwick_Cols = []
KWeight_Cols = []

for rFN in regionFitsNames:
    listColumnsMeca += ['KChadwick_'+rFN, 'K_CIW_'+rFN, 'R2Chadwick_'+rFN, 'K2Chadwick_'+rFN, 
                        'H0Chadwick_'+rFN, 'Npts_'+rFN, 'validatedFit_'+rFN]
    KChadwick_Cols += [('KChadwick_'+rFN)]

    K_CIWidth = data_f['K_CIW_'+rFN] #.apply(lambda x : x.strip('][').split(', ')).apply(lambda x : (np.abs(float(x[0]) - float(x[1]))))
    KWeight = (data_f['KChadwick_'+rFN]/K_CIWidth)**2
    data_f['K_Weight_'+rFN] = KWeight
    data_f['K_Weight_'+rFN] *= data_f['KChadwick_'+rFN].apply(lambda x : (x<1e6))
    data_f['K_Weight_'+rFN] *= data_f['R2Chadwick_'+rFN].apply(lambda x : (x>1e-2))
    data_f['K_Weight_'+rFN] *= data_f['K_CIW_'+rFN].apply(lambda x : (x!=0))
    KWeight_Cols += [('K_Weight_'+rFN)]
    

#### Whole curve

def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'



regionFitsNames = [str(fitMin[ii]) + '<s<' + str(fitMax[ii]) for ii in range(len(fitMin))]

fig, ax = plt.subplots(2,1, figsize = (9,12))
conditions = ['none', 'doxycyclin']
cD = {'none':[colorList40[10], colorList40[11]], 'doxycyclin':[colorList40[30], colorList40[31]]}
oD = {'none': [-15, 1.02] , 'doxycyclin': [5, 0.97] }
lD = {'none': 'Control' , 'doxycyclin': 'iMC linker' }

for drug in conditions:
    
    Kavg = []
    Kstd = []
    D10 = []
    D90 = []
    N = []
    
    data_ff = data_f[data_f['drug'] == drug]
    for ii in range(len(fitCenters)):
        rFN = str(fitMin[ii]) + '<s<' + str(fitMax[ii])
        variable = valStr+rFN
        weight = weightStr+rFN
        
        x = data_ff[variable].apply(nan2zero).values
        w = data_ff[weight].apply(nan2zero).values
        
        if S == 250:
            d = {'x' : x, 'w' : w}
        
        m = np.average(x, weights=w)
        v = np.average((x-m)**2, weights=w)
        std = v**0.5
        
        d10, d90 = np.percentile(x[x != 0], (10, 90))
        n = len(x[x != 0])
        
        Kavg.append(m)
        Kstd.append(std)
        D10.append(d10)
        D90.append(d90)
        N.append(n)
    
    Kavg = np.array(Kavg)
    Kstd = np.array(Kstd)
    D10 = np.array(D10)
    D90 = np.array(D90)
    N = np.array(N)
    Kste = Kstd / (N**0.5)
    
    alpha = 0.975
    dof = N
    q = st.t.ppf(alpha, dof) # Student coefficient
    
    d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}
    
    
    # Weighted means -- Weighted ste 95% as error
    ax[0].errorbar(fitCenters, Kavg, yerr = q*Kste, marker = 'o', color = cD[drug][0], 
                   ecolor = 'k', elinewidth = 0.8, capsize = 3, label = lD[drug])
    ax[0].set_ylim([500,2e4])
    ax[0].set_xlim([0,1000])
    ax[0].set_title('Stress stiffening - All compressions pooled')
    
    # Weighted means -- D9-D1 as error
    ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = cD[drug][1], 
                   ecolor = 'k', elinewidth = 0.8, capsize = 3, label = lD[drug]) 
    ax[1].set_ylim([500,2e4])
    ax[1].set_xlim([0,1000])
    
    for k in range(2):
        ax[k].legend(loc = 'upper left')
        ax[k].set_yscale('log')
        ax[k].set_xlabel('Stress (Pa)')
        ax[k].set_ylabel('K (Pa)')
        for kk in range(len(N)):
            ax[k].text(x=fitCenters[kk]+oD[drug][0], y=Kavg[kk]**oD[drug][1], s='n='+str(N[kk]), fontsize = 6, color = cD[drug][k])
    
    fig.suptitle('K(s) - All compressions pooled'+'\n(fits width: {:.0f}Pa)'.format(width))    
    
    df_val = pd.DataFrame(d_val)
    dftest = pd.DataFrame(d)

plt.show()

jvu.archiveFig(fig, ax, todayFigDirLocal + '//MCAproject', name='aSFL_ctrl-vs-doxy_K(s)', dpi = 100)

#### Local zoom


Sinf, Ssup = 300, 800
extraFilters = [data_f['minStress'] <= Sinf, data_f['maxStress'] >= Ssup] # >= 800
fitCenters = fitCenters[(fitCenters>=(Sinf)) & (fitCenters<=Ssup)] # <800
fitMin = np.array([int(S-(width/2)) for S in fitCenters])
fitMax = np.array([int(S+(width/2)) for S in fitCenters])

data2_f = data_f
globalExtraFilter = extraFilters[0]
for k in range(1, len(extraFilters)):
    globalExtraFilter = globalExtraFilter & extraFilters[k]
data2_f = data2_f[globalExtraFilter]  


fig, ax = plt.subplots(2,1, figsize = (7,12))
conditions = ['none', 'doxycyclin']
# cD = {'none':[colorList30[10], colorList30[11]], 'doxycyclin':[colorList30[20], colorList30[21]]}
# oD = {'none': 1.04 , 'doxycyclin': 0.96 }
# lD = {'none': 'Control' , 'doxycyclin': 'iMC linker' }


for drug in conditions:
    
    Kavg = []
    Kstd = []
    D10 = []
    D90 = []
    N = []
    
    data_ff = data2_f[data_f['drug'] == drug]
    
    for ii in range(len(fitCenters)):
        print(fitCenters[ii])
        rFN = str(fitMin[ii]) + '<s<' + str(fitMax[ii])
        variable = valStr+rFN
        weight = weightStr+rFN
        
        x = data_ff[variable].apply(nan2zero).values
        w = data_ff[weight].apply(nan2zero).values
        
        if S == 250:
            d = {'x' : x, 'w' : w}
        
        m = np.average(x, weights=w)
        v = np.average((x-m)**2, weights=w)
        std = v**0.5
        
        d10, d90 = np.percentile(x[x != 0], (10, 90))
        n = len(x[x != 0])
        
        Kavg.append(m)
        Kstd.append(std)
        D10.append(d10)
        D90.append(d90)
        N.append(n)
    
    Kavg = np.array(Kavg)
    Kstd = np.array(Kstd)
    D10 = np.array(D10)
    D90 = np.array(D90)
    N = np.array(N)
    Kste = Kstd / (N**0.5)
    
    alpha = 0.975
    dof = N
    q = st.t.ppf(alpha, dof) # Student coefficient
    
    d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}
    
    
    # Weighted means -- Weighted ste 95% as error
    ax[0].errorbar(fitCenters, Kavg, yerr = q*Kste, marker = 'o', color = cD[drug][0], 
                   ecolor = 'k', elinewidth = 0.8, capsize = 3, label = lD[drug])
    ax[0].set_ylim([500,2e4])
    ax[0].set_xlim([200,900])
    ax[0].set_title('Stress stiffening\nOnly compressions including the [300,800]Pa range')
    
    # Weighted means -- D9-D1 as error
    ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = cD[drug][1], 
                   ecolor = 'k', elinewidth = 0.8, capsize = 3, label = lD[drug]) 
    ax[1].set_ylim([500,2e4])
    ax[1].set_xlim([200,900])
    
    for k in range(2):
        ax[k].legend(loc = 'upper left')
        ax[k].set_yscale('log')
        ax[k].set_xlabel('Stress (Pa)')
        ax[k].set_ylabel('K (Pa)')
        for kk in range(len(N)):
            ax[k].text(x=fitCenters[kk]+oD[drug][0], y=Kavg[kk]**oD[drug][1], s='n='+str(N[kk]), fontsize = 6, color = cD[drug][k])
    
    fig.suptitle('K(s) - All compressions pooled'+'\n(fits width: {:.0f}Pa)'.format(width))
    
    df_val = pd.DataFrame(d_val)
    dftest = pd.DataFrame(d)

plt.show()

jvu.archiveFig(fig, ax, todayFigDirLocal + '//MCAproject', name='aSFL_ctrl-vs-doxy_K(s)_localZoom', dpi = 100)



# %%%%%


# print(data_f.head())

# print(fitCenters)

def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in fitCenters:
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_f[variable].apply(nan2zero).values
    w = data_f[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# ax[0]
ax.errorbar(fitCenters, Kavg, yerr = q*Kste, marker = 'o', color = my_default_color_list[0], 
               ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error')
ax.set_ylim([500,2e4])

# ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# ax[1].set_ylim([500,1e6])

# for k in range(1): #2
#     ax[k].legend(loc = 'upper left')
#     ax[k].set_yscale('log')
#     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
#     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
#     for kk in range(len(N)):
#         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
ax.legend(loc = 'upper left')
ax.set_yscale('log')
ax.set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
ax.set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
for kk in range(len(N)):
    ax.text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

fig.suptitle('K(sigma)')
ax.set_title('From all compressions pooled\n22-02-09 experiment, 36 cells, 232 compression')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_K(s)globalAvg_V2', figSubDir = 'NonLin')
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


# %%%%%


# print(data_f.head())

# print(fitCenters)

extraFilters = [data_f['minStress'] <= 100, data_f['maxStress'] >= 800]
fitCenters2 = fitCenters[fitCenters<850]

data_ff = data_f
globalExtraFilter = extraFilters[0]
for k in range(1, len(extraFilters)):
    globalExtraFilter = globalExtraFilter & extraFilters[k]

data_ff = data_ff[globalExtraFilter]
data_ff
def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in fitCenters2:
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_ff[variable].apply(nan2zero).values
    w = data_ff[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# ax[0]
ax.errorbar(fitCenters2, Kavg, yerr = q*Kste, marker = 'o', color = my_default_color_list[0], 
               ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error')
ax.set_ylim([500,2e4])

# ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# ax[1].set_ylim([500,1e6])

# for k in range(1): #2
#     ax[k].legend(loc = 'upper left')
#     ax[k].set_yscale('log')
#     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
#     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
#     for kk in range(len(N)):
#         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
ax.legend(loc = 'upper left')
ax.set_yscale('log')
ax.set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
ax.set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
for kk in range(len(N)):
    ax.text(x=fitCenters2[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

fig.suptitle('K(sigma)')
ax.set_title('Only compressions including the [100, 800]Pa range')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_K(s)globalAvg_100to800', figSubDir = 'NonLin')
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


# %%%%%


# print(data_f.head())

# print(fitCenters)

extraFilters = [data_f['minStress'] <= 100, data_f['maxStress'] >= 600]
fitCenters2 = fitCenters[fitCenters<650]

data_ff = data_f
globalExtraFilter = extraFilters[0]
for k in range(1, len(extraFilters)):
    globalExtraFilter = globalExtraFilter & extraFilters[k]

data_ff = data_ff[globalExtraFilter]
data_ff
def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in fitCenters2:
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_ff[variable].apply(nan2zero).values
    w = data_ff[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# ax[0]
ax.errorbar(fitCenters2, Kavg, yerr = q*Kste, marker = 'o', color = my_default_color_list[0], 
               ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error')
ax.set_ylim([500,2e4])

# ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# ax[1].set_ylim([500,1e6])

# for k in range(1): #2
#     ax[k].legend(loc = 'upper left')
#     ax[k].set_yscale('log')
#     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
#     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
#     for kk in range(len(N)):
#         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
ax.legend(loc = 'upper left')
ax.set_yscale('log')
ax.set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
ax.set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
for kk in range(len(N)):
    ax.text(x=fitCenters2[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

fig.suptitle('K(sigma)')
ax.set_title('Only compressions including the [100, 600]Pa range')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_K(s)globalAvg_100to600', figSubDir = 'NonLin')
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


Kavg_ref, Kste_ref, N_ref, fitCenters_ref = Kavg, Kste, N, fitCenters2



# %%% Asym bead pairs -- (july 2021)

# %%%%%


rawMecaTable = jva.getGlobalTable_meca()
# rawMecaTable.head()


# %%%%%


#  GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270_M270', 'M270_M450', 'M450_M270', 'M450_M450'])


# %%%%%


Outliers = ['21-07-08_M3_P1_C10', '21-07-08_M3_P1_C5', '21-07-08_M3_P1_C6']
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), (GlobalTable_meca_Py['drug'] == 'none'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'),           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers)),            (GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270', 'M270_M450', 'M450_M270', 'M450']))]
co_order = makeOrder(['M270', 'M270_M450', 'M450_M270', 'M450'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['bead type'],                 Parameters=['surroundingThickness','EChadwick'],Filters=Filters,                 AvgPerCell=False, cellID='cellID', co_order=co_order, figSizeFactor = 2, orientation = 'v', useHue = False)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL asym bead pairs - all comps')
# jvu.archiveFig(fig, ax, name='3T3aSFL_asymBeadsAllComp_SurroundingThickness&EChadwick', figSubDir = figSubDir)
plt.show()


# %%%%%


Outliers = ['21-07-08_M3_P1_C10','21-07-08_M3_P1_C5','21-07-08_M3_P1_C6']
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), (GlobalTable_meca_Py['drug'] == 'none'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'),            (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers)),            (GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270', 'M270_M450', 'M450_M270', 'M450']))]
co_order = makeOrder(['M270', 'M270_M450', 'M450_M270', 'M450'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['bead type'],                 Parameters=['surroundingThickness','EChadwick'],Filters=Filters,                 AvgPerCell=True, cellID='cellID', co_order=co_order, figSizeFactor = 2, orientation = 'v', useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL asym bead pairs - avg per cell')
# jvu.archiveFig(fig, ax, name='3T3aSFL_asymBeadsPerCell_SurroundingThickness&EChadwick', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(rawMecaTable['Validated'] == 1), ((rawMecaTable['ExpType'] == 'DictyDB_M270') | (rawMecaTable['ExpType'] == 'DictyDB_M450')), (rawMecaTable['TpsComp'] == '1s')]
# Filters = [(rawMecaTable['Validated'] == 1), ((rawMecaTable['ExpType'] == 'DictyDB_M450')), (rawMecaTable['TpsComp'] == '1s')]
fig, ax = D1Plot(rawMecaTable, CondCol=['ExpType'],Parameters=['SurroundingThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='CellName', useHue = False)
fig.suptitle('M450 vs M270 pour compressions de 1s')
renameAxes(ax,renameDict1)
jvu.archiveFig(fig, ax, name='Dictys_beadTypes_SurroundingThickness&EChadwick', figSubDir = figSubDir)
plt.show()
# rawMecaTable[Filters[0] & Filters[1] & Filters[2]]


# %%%  Big vs Small Beads, new code

# %%%% Oct-2021

df0 = GlobalTable_meca_Py2.tail()


# %%%%% 3T3aSFL - Oct 21 experiment, bead types
data = GlobalTable_meca_Py2
Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date_x'].apply(lambda x : x in ['21-10-18', '21-10-25']))]
co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(data, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick_f<150pN'],
                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - Oct 21 experiment, bead types')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='3T3aSFL_Oct21_M450-1-10_vs_M270-5-40')
plt.show()


# %%%%% 3T3aSFL - Oct 21 experiment, bead types, 200-400nm

data = GlobalTable_meca_Py2
Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date_x'].apply(lambda x : x in ['21-10-18', '21-10-25'])),
           (data['ctFieldThickness'] < 400),
          (data['ctFieldThickness'] > 200)]
co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(data, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick_f<150pN'],
                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - Oct 21 experiment, bead types, 200-400nm')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='3T3aSFL_Oct21_M450-1-10_vs_M270-5-40_h=200-400nm')
plt.show()


# %%%%% 3T3aSFL - Oct 21 experiment, bead types, 400-600nm

data = GlobalTable_meca_Py2
Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date_x'].apply(lambda x : x in ['21-10-18', '21-10-25'])),
           (data['ctFieldThickness'] < 600),
          (data['ctFieldThickness'] > 400)]

co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(data, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick_f<150pN'],
                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - Oct 21 experiment, bead types, 400-600nm')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='3T3aSFL_Oct21_M450-1-10_vs_M270-5-40_h=400-600nm')
plt.show()


# %%%%% 3T3aSFL - Oct 21 experiment, bead types, 600-1000nm

data = GlobalTable_meca_Py2
Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date_x'].apply(lambda x : x in ['21-10-18', '21-10-25'])),
           (data['ctFieldThickness'] < 1000),
          (data['ctFieldThickness'] > 600)]

co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(data, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick_f<150pN'],
                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)

renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - Oct 21 experiment, bead types, 600-1000nm')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='3T3aSFL_Oct21_M450-1-10_vs_M270-5-40_h=600-1000nm')
plt.show()


# %%%%% 3T3aSFL: E(h)

data = GlobalTable_meca_Py2
Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date_x'].apply(lambda x : x in ['21-10-18', '21-10-25']))]

fig, ax = D2Plot(data, XCol='ctFieldThickness',YCol='EChadwick_f<150pN',CondCol = ['bead type'],
                 Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=False, 
                 modelType='y=A*exp(kx)',xscale = 'log', yscale = 'log')

fig.suptitle('3T3aSFL: E(h)')
# jvu.archiveFig(fig, ax, name='E(h)_3T3aSFL_Oct21_M450-1-10_vs_M270-5-40', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% 3T3aSFL: dz(dx) - all compressions

data = GlobalTable_meca_Py2
Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date_x'].apply(lambda x : x in ['21-10-18', '21-10-25']))]

fig, ax = D2Plot(data, XCol='surroundingDz',YCol='surroundingDx',CondCol = ['bead type'],
                 Filters=Filters, cellID = 'cellID', AvgPerCell=False, modelFit=False, 
                 xscale = 'lin', yscale = 'lin')

fig.suptitle('3T3aSFL: dz(dx) - all compressions')
# jvu.archiveFig(fig, ax, name='dz(dx)_3T3aSFL_Oct21_M450-1-10_vs_M270-5-40', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% 3T3aSFL: dz(dx) - median per cell

data = GlobalTable_meca_Py2
Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date_x'].apply(lambda x : x in ['21-10-18', '21-10-25']))]

fig, ax = D2Plot(data, XCol='ctFieldDZ',YCol='ctFieldDX',CondCol = ['bead type'],
                 Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=False, 
                 xscale = 'lin', yscale = 'lin')

fig.suptitle('3T3aSFL: dz(dx) - median per cell')
# jvu.archiveFig(fig, ax, name='DZ(DX)_3T3aSFL_Oct21_M450-1-10_vs_M270-5-40', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%% Compute nice field range matching stress-wise

# %%%%% Understanding the stress range mismatch

h0 = 500
B_M270 = 50
B_M450 = 8.5
anglefactor = 2
V_M270 = (4/3)*np.pi*(2690/2)**3 # volume [nm^3]
V_M450 = (4/3)*np.pi*(4503/2)**3 # volume [nm^3]
m_M270 = jvu.computeMag_M270computeMag_M270(B_M270) * 1e-9 * V_M270
m_M450 = jvu.computeMag_M270computeMag_M450(B_M450) * 1e-9 * V_M450
D3nm_270 = np.arange(1, h0, 1, dtype = float) + 2690
D3nm_450 = np.arange(1, h0, 1, dtype = float) + 4503
F_270 = 3e5*anglefactor*m_M270**2/D3nm_270**4
F_450 = 3e5*anglefactor*m_M450**2/D3nm_450**4

s_270 = F_270[:-20]/((2690/2e6)*(h0+2690-D3nm_270[:-20]))
s_450 = F_450[:-20]/((4503/2e6)*(h0+4503-D3nm_450[:-20]))

fig, ax = plt.subplots(2,1, figsize = (8, 8))
ax[0].plot(D3nm_270-2690, F_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[0].plot(D3nm_450-4503, F_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[0].legend()
ax[0].set_ylabel('F (pN)')
ax[1].plot(D3nm_270[:-20]-2690, s_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[1].plot(D3nm_450[:-20]-4503, s_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
# ax[1].plot(ax[1].get_xlim(), [60, 60], 'k--')
ax[1].legend()
ax[1].set_ylabel('sigma for h0=' + str(h0) + 'nm (Pa)')
fig.suptitle('Understanding the stress range mismatch\nTest', fontsize = 16)
#jvu.archiveFig(fig, ax, name='Mismatch_3T3aSFL_Test2', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% Value used for standard experiments : 5mT

h0 = 500
B_M270 = 5
B_M450 = 5
anglefactor = 2
V_M270 = (4/3)*np.pi*(2690/2)**3 # volume [nm^3]
V_M450 = (4/3)*np.pi*(4503/2)**3 # volume [nm^3]
m_M270 = jvu.computeMag_M270computeMag_M270(B_M270) * 1e-9 * V_M270
m_M450 = jvu.computeMag_M270computeMag_M450(B_M450) * 1e-9 * V_M450
D3nm_270 = np.arange(1, h0, 1, dtype = float) + 2690
D3nm_450 = np.arange(1, h0, 1, dtype = float) + 4503
F_270 = 3e5*anglefactor*m_M270**2/D3nm_270**4
F_450 = 3e5*anglefactor*m_M450**2/D3nm_450**4

s_270 = F_270[:-20]/((2690/2e6)*(h0+2690-D3nm_270[:-20]))
s_450 = F_450[:-20]/((4503/2e6)*(h0+4503-D3nm_450[:-20]))

fig, ax = plt.subplots(2,1, figsize = (8, 8))
ax[0].plot(D3nm_270-2690, F_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[0].plot(D3nm_450-4503, F_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[0].legend()
ax[0].set_ylabel('F (pN)')
ax[1].plot(D3nm_270[:-20]-2690, s_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[1].plot(D3nm_450[:-20]-4503, s_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[1].legend()
ax[1].set_ylabel('sigma for h0=' + str(h0) + 'nm (Pa)')
fig.suptitle('Understanding the stress range mismatch\nValue used for standard experiments : 5mT', fontsize = 16)
#jvu.archiveFig(fig, ax, name='Mismatch_3T3aSFL_Oct21_M450_5mT-vs-M270_5mT', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% Values used for october 2021 experiment

h0 = 500
B_M270 = 5
B_M450 = 1
anglefactor = 2
V_M270 = (4/3)*np.pi*(2690/2)**3 # volume [nm^3]
V_M450 = (4/3)*np.pi*(4503/2)**3 # volume [nm^3]
m_M270 = jvu.computeMag_M270computeMag_M270(B_M270) * 1e-9 * V_M270
m_M450 = jvu.computeMag_M270computeMag_M450(B_M450) * 1e-9 * V_M450
D3nm_270 = np.arange(1, h0, 1, dtype = float) + 2690
D3nm_450 = np.arange(1, h0, 1, dtype = float) + 4503
F_270 = 3e5*anglefactor*m_M270**2/D3nm_270**4
F_450 = 3e5*anglefactor*m_M450**2/D3nm_450**4

s_270 = F_270[:-20]/((2690/2e6)*(h0+2690-D3nm_270[:-20]))
s_450 = F_450[:-20]/((4503/2e6)*(h0+4503-D3nm_450[:-20]))

fig, ax = plt.subplots(2,1, figsize = (8, 8))
ax[0].plot(D3nm_270-2690, F_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[0].plot(D3nm_450-4503, F_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[0].legend()
ax[0].set_ylabel('F (pN)')
ax[1].plot(D3nm_270[:-20]-2690, s_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[1].plot(D3nm_450[:-20]-4503, s_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
# ax[1].plot(ax[1].get_xlim(), [60, 60], 'k--')
ax[1].legend()
ax[1].set_ylabel('sigma for h0=' + str(h0) + 'nm (Pa)')
fig.suptitle('Understanding the stress range mismatch\nValues used for october 2021 experiment', fontsize = 16)
#jvu.archiveFig(fig, ax, name='Mismatch_3T3aSFL_Oct21_M450_1mT-vs-M270_5mT', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% Proposition of values for a future experiment 1

h0 = 500
B_M270 = 14
B_M450 = 5
anglefactor = 2
V_M270 = (4/3)*np.pi*(2690/2)**3 # volume [nm^3]
V_M450 = (4/3)*np.pi*(4503/2)**3 # volume [nm^3]
m_M270 = jvu.computeMag_M270computeMag_M270(B_M270) * 1e-9 * V_M270
m_M450 = jvu.computeMag_M270computeMag_M450(B_M450) * 1e-9 * V_M450
D3nm_270 = np.arange(1, h0, 1, dtype = float) + 2690
D3nm_450 = np.arange(1, h0, 1, dtype = float) + 4503
F_270 = 3e5*anglefactor*m_M270**2/D3nm_270**4
F_450 = 3e5*anglefactor*m_M450**2/D3nm_450**4

s_270 = F_270[:-20]/((2690/2e6)*(h0+2690-D3nm_270[:-20]))
s_450 = F_450[:-20]/((4503/2e6)*(h0+4503-D3nm_450[:-20]))

fig, ax = plt.subplots(2,1, figsize = (8, 8))
ax[0].plot(D3nm_270-2690, F_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[0].plot(D3nm_450-4503, F_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[0].legend()
ax[0].set_ylabel('F (pN)')
ax[1].plot(D3nm_270[:-20]-2690, s_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[1].plot(D3nm_450[:-20]-4503, s_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[1].legend()
ax[1].set_ylabel('sigma for h0=' + str(h0) + 'nm (Pa)')
fig.suptitle('Understanding the stress range mismatch\nProposition of values for a future experiment', fontsize = 16)
# jvu.archiveFig(fig, ax, name='Mismatch_3T3aSFL_Oct21_M450_5mT-vs-M270_14mT', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% Proposition of values for a future experiment 2

h0 = 500
B_M270 = 54
B_M450 = 13
anglefactor = 2
V_M270 = (4/3)*np.pi*(2690/2)**3 # volume [nm^3]
V_M450 = (4/3)*np.pi*(4503/2)**3 # volume [nm^3]
m_M270 = jvu.computeMag_M270(B_M270) * 1e-9 * V_M270
m_M450 = jvu.computeMag_M270computeMag_M450(B_M450) * 1e-9 * V_M450
D3nm_270 = np.arange(1, h0, 1, dtype = float) + 2690
D3nm_450 = np.arange(1, h0, 1, dtype = float) + 4503
F_270 = 3e5*anglefactor*m_M270**2/D3nm_270**4
F_450 = 3e5*anglefactor*m_M450**2/D3nm_450**4

s_270 = F_270[:-20]/((2690/2e6)*(h0+2690-D3nm_270[:-20]))
s_450 = F_450[:-20]/((4503/2e6)*(h0+4503-D3nm_450[:-20]))

fig, ax = plt.subplots(2,1, figsize = (8, 8))
ax[0].plot(D3nm_270-2690, F_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[0].plot(D3nm_450-4503, F_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[0].legend()
ax[0].set_ylabel('F (pN)')
ax[1].plot(D3nm_270[:-20]-2690, s_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[1].plot(D3nm_450[:-20]-4503, s_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[1].legend()
ax[1].set_ylabel('sigma for h0=' + str(h0) + 'nm (Pa)')
fig.suptitle('Understanding the stress range mismatch\nProposition of values for a future experiment', fontsize = 16)
# jvu.archiveFig(fig, ax, name='Mismatch_3T3aSFL_Oct21_M450_13mT-vs-M270_54mT', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% Proposition of values for a future experiment 3

h0 = 500
B_M270 = 2.1
B_M450 = 1
anglefactor = 2
V_M270 = (4/3)*np.pi*(2690/2)**3 # volume [nm^3]
V_M450 = (4/3)*np.pi*(4503/2)**3 # volume [nm^3]
m_M270 = jvu.computeMag_M270computeMag_M270(B_M270) * 1e-9 * V_M270
m_M450 = jvu.computeMag_M270computeMag_M450(B_M450) * 1e-9 * V_M450
D3nm_270 = np.arange(1, h0, 1, dtype = float) + 2690
D3nm_450 = np.arange(1, h0, 1, dtype = float) + 4503
F_270 = 3e5*anglefactor*m_M270**2/D3nm_270**4
F_450 = 3e5*anglefactor*m_M450**2/D3nm_450**4

s_270 = F_270[:-20]/((2690/2e6)*(h0+2690-D3nm_270[:-20]))
s_450 = F_450[:-20]/((4503/2e6)*(h0+4503-D3nm_450[:-20]))

fig, ax = plt.subplots(2,1, figsize = (8, 8))
ax[0].plot(D3nm_270-2690, F_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[0].plot(D3nm_450-4503, F_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[0].legend()
ax[0].set_ylabel('F (pN)')
ax[1].plot(D3nm_270[:-20]-2690, s_270, 'c', label = 'M270 - B='+str(B_M270)+'mT')
ax[1].plot(D3nm_450[:-20]-4503, s_450, 'r', label = 'M450 - B='+str(B_M450)+'mT')
ax[1].legend()
ax[1].set_ylabel('sigma for h0=' + str(h0) + 'nm (Pa)')
fig.suptitle('Understanding the stress range mismatch\nProposition of values for a future experiment', fontsize = 16)
# jvu.archiveFig(fig, ax, name='Mismatch_3T3aSFL_Oct21_M450_13mT-vs-M270_54mT', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%% min>max stress for each compressions, 2 bead sizes - V1

data = GlobalTable_meca_Py2

Filters270 = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),\
           (data['cell subtype'] == 'aSFL'), 
           (data['bead type'] == 'M270'), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in ['21-10-18', '21-10-25']))]
Filter270 = Filters270[0]
for i in range(1, len(Filters270)):
    Filter270 = Filter270 & Filters270[i]
    
Filters450 = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),\
           (data['cell subtype'] == 'aSFL'), 
           (data['bead type'] == 'M450'), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in ['21-10-18', '21-10-25']))]   
Filter450 = Filters450[0]
for i in range(1, len(Filters450)):
    Filter450 = Filter450 & Filters450[i]
    

fig, ax = plt.subplots(1,2, figsize = (8, 6))
data[Filter270].plot(kind = 'scatter', ax = ax[0], x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'Min stress')
data[Filter270].plot(kind = 'scatter', ax = ax[0], x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'Max stress')
data[Filter450].plot(kind = 'scatter', ax = ax[1], x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'Min stress')
data[Filter450].plot(kind = 'scatter', ax = ax[1], x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'Max stress')
ax[0].set_xlabel('Cortical thickness (nm)')
ax[0].set_ylabel('Extremal stress values (Pa)')
ax[0].set_ylim([0, 800])
ax[0].set_xlim([0, 1200])
ax[0].set_title('M270')
ax[1].set_xlabel('Cortical thickness (nm)')
ax[1].set_ylabel('Extremal stress values (Pa)')
ax[1].set_ylim([0, 800])
ax[1].set_xlim([0, 1200])
ax[1].set_title('M450')
plt.show()


# %%%%% min>max stress for each compressions, 2 bead sizes - V2

fig, ax = plt.subplots(1,1, figsize = (8, 6))
data[Filter270].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'M270 - Min stress')
data[Filter270].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'orange', alpha = 0.3, label = 'M270 - Max stress')
data[Filter450].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'b', alpha = 0.3, label = 'M450 - Min stress')
data[Filter450].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'M450 - Max stress')
ax.set_xlabel('Cortical thickness (nm)')
ax.set_ylabel('Extremal stress values (Pa)')
ax.set_ylim([0, 1000])
ax.set_xlim([0, 1200])
plt.show()


# %%%% Dec-2021

# %%%%%

df0= GlobalTable_meca_Py2.head()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08']))]
co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(GlobalTable_meca_Py2, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick_f<150pN'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - Dec 21 experiment, bead types')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-16']))]
co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(GlobalTable_meca_Py2, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - Dec 21 experiment, bead types')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08']))]

co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(GlobalTable_meca_Py2, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL, Dec 21 - Matching ranges of stress for 2 bead types\nB{M270}=14>54mT ; B{M450}=5>13mT')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08'])),
           (GlobalTable_meca_Py2['ctFieldThickness'] < 400),
           (GlobalTable_meca_Py2['ctFieldThickness'] > 150)]

co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(GlobalTable_meca_Py2, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL, Dec 21 - Matching ranges of stress for 2 bead types\nB{M270}=14>54mT ; B{M450}=5>13mT')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08'])),
           (GlobalTable_meca_Py2['ctFieldThickness'] < 600),
           (GlobalTable_meca_Py2['ctFieldThickness'] > 400)]

co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(GlobalTable_meca_Py2, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL, Dec 21 - Matching ranges of stress for 2 bead types\nB{M270}=14>54mT ; B{M450}=5>13mT')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08']))]
fig, ax = D2Plot(GlobalTable_meca_Py2, XCol='ctFieldThickness',YCol='EChadwick',CondCol = ['bead type'],           Filters=Filters, cellID = 'cellID', AvgPerCell=True, xscale = 'log', yscale = 'log')
ax.set_ylabel('EChadwick (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-16']))]
fig, ax = D2Plot(GlobalTable_meca_Py2, XCol='ctFieldThickness',YCol='EChadwick',CondCol = ['bead type'],           Filters=Filters, cellID = 'cellID', AvgPerCell=True,xscale = 'log', yscale = 'log')
ax.set_ylabel('EChadwick (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

# plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08','21-12-16']))]

co_order = makeOrder(['M270','M450'], ['21-12-08','21-12-16'])

fig, ax = D1Plot(GlobalTable_meca_Py2, CondCol=['bead type', 'date'],Parameters=['ctFieldThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True, orientation = 'v',
                 figSizeFactor = 0.9)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL, Dec 21 - Matching ranges of stress for 2 bead types\nB{M270}=14>54mT ; B{M450}=5>13mT')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

# plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-16','21-12-08']))]
fig, ax = D2Plot(GlobalTable_meca_Py2, XCol='ctFieldThickness',YCol='EChadwick',CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=True, xscale = 'log', yscale = 'log', markers = ['o', 'o', '>', '>'])
ax.set_ylabel('EChadwick (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')
ax.legend(loc = 'upper right')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

# plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-16','21-12-08']))]
fig, ax = D2Plot_wFit(GlobalTable_meca_Py2, XCol='ctFieldThickness',YCol='EChadwick',CondCol = ['bead type'],
           Filters=Filters, cellID = 'cellID', AvgPerCell=True, xscale = 'log', yscale = 'log', 
           modelFit=True, modelType = 'y=k*x^a')
ax.set_ylabel('EChadwick (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')
ax.legend(loc = 'upper right')

# jvu.archiveFig(fig, ax, 
#                os.path.join(todayFigDirLocal, 'BigVsSmallPlots'), 
#                name='')

plt.show()


# %%%%%


Filters270 = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True),\
           (GlobalTable_meca_Py2['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py2['bead type'] == 'M270'), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08']))]
Filter270 = Filters270[0]
for i in range(1, len(Filters270)):
    Filter270 = Filter270 & Filters270[i]
    
Filters450 = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True),\
           (GlobalTable_meca_Py2['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py2['bead type'] == 'M450'), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08']))]   
Filter450 = Filters450[0]
for i in range(1, len(Filters450)):
    Filter450 = Filter450 & Filters450[i]
    

fig, ax = plt.subplots(1,2, figsize = (9, 5))
GlobalTable_meca_Py2[Filter270].plot(kind = 'scatter', ax = ax[0], x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'Min stress')
GlobalTable_meca_Py2[Filter270].plot(kind = 'scatter', ax = ax[0], x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'Max stress')
GlobalTable_meca_Py2[Filter450].plot(kind = 'scatter', ax = ax[1], x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'Min stress')
GlobalTable_meca_Py2[Filter450].plot(kind = 'scatter', ax = ax[1], x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'Max stress')
ax[0].set_xlabel('Cortical thickness (nm)')
ax[0].set_ylabel('Extremal stress values (Pa)')
ax[0].set_ylim([0, 600])
ax[0].set_xlim([0, 1200])
ax[0].set_title('M270')
ax[1].set_xlabel('Cortical thickness (nm)')
ax[1].set_ylabel('Extremal stress values (Pa)')
ax[1].set_ylim([0, 600])
ax[1].set_xlim([0, 1200])
ax[1].set_title('M450')
# plt.show()


# %%%%%


fig, ax = plt.subplots(1,1, figsize = (9, 5))
GlobalTable_meca_Py2[Filter270].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'M270 - Min stress')
GlobalTable_meca_Py2[Filter270].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'orange', alpha = 0.3, label = 'M270 - Max stress')
GlobalTable_meca_Py2[Filter450].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'b', alpha = 0.3, label = 'M450 - Min stress')
GlobalTable_meca_Py2[Filter450].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'M450 - Max stress')
ax.set_xlabel('Cortical thickness (nm)')
ax.set_ylabel('Extremal stress values (Pa)')
ax.set_title('Dec 21 - Extremal stress values vs. thickness,\nfor each compression')
# plt.show()


# %%%%%


Filters270 = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True),\
           (GlobalTable_meca_Py2['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py2['bead type'] == 'M270'), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-16']))]
Filter270 = Filters270[0]
for i in range(1, len(Filters270)):
    Filter270 = Filter270 & Filters270[i]
    
Filters450 = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True),\
           (GlobalTable_meca_Py2['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py2['bead type'] == 'M450'), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-16']))]   
Filter450 = Filters450[0]
for i in range(1, len(Filters450)):
    Filter450 = Filter450 & Filters450[i]
    

fig, ax = plt.subplots(1,2, figsize = (9, 5))
GlobalTable_meca_Py2[Filter270].plot(kind = 'scatter', ax = ax[0], x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'Min stress')
GlobalTable_meca_Py2[Filter270].plot(kind = 'scatter', ax = ax[0], x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'Max stress')
GlobalTable_meca_Py2[Filter450].plot(kind = 'scatter', ax = ax[1], x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'Min stress')
GlobalTable_meca_Py2[Filter450].plot(kind = 'scatter', ax = ax[1], x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'Max stress')
ax[0].set_xlabel('Cortical thickness (nm)')
ax[0].set_ylabel('Extremal stress values (Pa)')
ax[0].set_ylim([0, 2000])
ax[0].set_xlim([0, 1200])
ax[0].set_title('M270')
ax[1].set_xlabel('Cortical thickness (nm)')
ax[1].set_ylabel('Extremal stress values (Pa)')
ax[1].set_ylim([0, 2000])
ax[1].set_xlim([0, 1200])
ax[1].set_title('M450')
# plt.show()


# %%%%%


fig, ax = plt.subplots(1,1, figsize = (9, 5))
GlobalTable_meca_Py2[Filter270].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'c', alpha = 0.3, label = 'M270 - Min stress')
GlobalTable_meca_Py2[Filter270].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'orange', alpha = 0.3, label = 'M270 - Max stress')
GlobalTable_meca_Py2[Filter450].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'minStress', marker = 'o', color = 'b', alpha = 0.3, label = 'M450 - Min stress')
GlobalTable_meca_Py2[Filter450].plot(kind = 'scatter', ax = ax, x = 'surroundingThickness', y = 'maxStress', marker = 'o', color = 'r', alpha = 0.3, label = 'M450 - Max stress')
ax.set_xlabel('Cortical thickness (nm)')
ax.set_ylabel('Extremal stress values (Pa)')
ax.set_title('Dec 21 - Extremal stress values vs. thickness,\nfor each compression')
# plt.show()


# %%%% Jan-2022

# %%%%%


figSubDir = 'BigVsSmallPlots'


# %%%%%


GlobalTable_meca_Py2.head()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08', '22-01-12']))]
co_order = makeOrder(['M270','M450'])

fig, ax = D1Plot(GlobalTable_meca_Py2, CondCol=['bead type'],Parameters=['ctFieldThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',co_order=co_order,stats=True, useHue = True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - Different bead diameters')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Simple1DPlot', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['22-01-12']))]
fig, ax = D2Plot(GlobalTable_meca_Py2, XCol='ctFieldThickness',YCol='EChadwick',CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=True, xscale = 'log', yscale = 'log', markers = ['o', 'o'])
ax.set_ylabel('EChadwick (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')
ax.legend(loc = 'upper right')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in ['21-12-08', '21-12-16', '22-01-12']))]
fig, ax = D2Plot(GlobalTable_meca_Py2, XCol='ctFieldThickness',YCol='EChadwick',CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=True, xscale = 'log', yscale = 'log', markers = ['o', 'o','v', 'v', 's', 's'])
ax.set_ylabel('EChadwick (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')
ax.legend(loc = 'upper right')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)_all3exp', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


# Test getAggDf(df, cellID, CondCol, Variables)
data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True)]

data_f = data
for fltr in Filters:
    data_f = data_f.loc[fltr]

dfA = getAggDf(data_f, 'cellID', 'bead type', ['surroundingThickness', 'EChadwick'])
dfA.head()


# %%%%% D1PlotDetailed & NonDetailed M270 vs M450
plt.close('all')

data = GlobalTable_meca_Py2
dates = ['21-12-08', '22-01-12'] # '21-12-08', '22-01-12'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['bestH0'] <= 900),
           (data['bestH0'] >= 150),
           (data['date'].apply(lambda x : x in dates))] #, '21-12-16'
co_order  = makeOrder(['M270', 'M450'])
rD = {'bestH0' : 'H0 (nm)',
      'EChadwick' : 'Elastic modulus (Pa)'}

fig, ax = D1PlotDetailed(data, CondCol=['bead type'], Parameters=['bestH0', 'EChadwick'], 
                         Filters=Filters, Boxplot=True, cellID='cellID', 
                         co_order=co_order, stats=True, statMethod='Mann-Whitney', 
                         box_pairs=[], figSizeFactor = 1.5, markersizeFactor=1, 
                         orientation = 'v', showManips = True)
renameAxes(ax, rD)
ax[0].set_ylim([0, 1000])
ax[1].set_ylim([2e2, 2e4])

# jvu.archiveFig(fig, ax, todayFigDirLocal + '//BeadSizes', name='BeadSize_MatchDetailed', dpi = 100)
plt.show()


fig, ax = D1Plot(data, CondCol=['bead type'], Parameters=['bestH0', 'EChadwick'], 
                 Filters=Filters, AvgPerCell = False, Boxplot=True, cellID='cellID', 
                 co_order=co_order, stats=True, statMethod='Mann-Whitney', 
                 box_pairs=[], figSizeFactor = 1.2, markersizeFactor=1, 
                 stressBoxPlot = True, orientation = 'v')
renameAxes(ax, rD)
ax[0].set_ylim([0, 1000])
ax[1].set_ylim([2e2, 2e4])

# jvu.archiveFig(fig, ax, todayFigDirLocal + '//BeadSizes', name='BeadSize_MatchAllComps', dpi = 100)
plt.show()


fig, ax = D1Plot(data, CondCol=['bead type'], Parameters=['bestH0', 'EChadwick'], 
                 Filters=Filters, AvgPerCell = True, Boxplot=True, cellID='cellID', 
                 co_order=co_order, stats=True, statMethod='Mann-Whitney', 
                 box_pairs=[], figSizeFactor = 1.2, markersizeFactor=1, 
                 stressBoxPlot = False, orientation = 'v')
renameAxes(ax, rD)
ax[0].set_ylim([0, 1000])
ax[1].set_ylim([2e2, 2e4])

# jvu.archiveFig(fig, ax, todayFigDirLocal + '//BeadSizes', name='BeadSize_MatchAvg', dpi = 100)
plt.show()




# %%%%% 3 days of exp to show the anomaly


plt.close('all')
data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True),
          (data['date'].apply(lambda x : x in ['21-12-08', '21-12-16', '22-01-12']))] #, '21-12-16'

box_pairs=[('21-12-08 & M270', '21-12-08 & M450'),
             ('21-12-16 & M450', '21-12-16 & M270'),
             ('22-01-12 & M270', '22-01-12 & M450')]

fig, ax = D1PlotDetailed(data, CondCol=['date', 'bead type'], Parameters=['surroundingThickness', 'EChadwick'], Filters=Filters, 
                Boxplot=True, cellID='cellID', co_order=[], stats=True, statMethod='Mann-Whitney', 
               box_pairs=box_pairs, figSizeFactor = 0.9, markersizeFactor=1, orientation = 'v', showManips = True)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Detailed1DPlot_all3exp', figSubDir = figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True),
          (data['date'].apply(lambda x : x in ['21-12-08', '22-01-12']))] #, '21-12-16'

data_f = data
for fltr in Filters:
    data_f = data_f.loc[fltr]

rangeStress = np.arange(0, 2001)
countStress = [np.zeros(2001), np.zeros(2001)]

data_M270 = data_f.loc[data_f['bead type'] == 'M270']
data_M450 = data_f.loc[data_f['bead type'] == 'M450']

allData = [data_M270, data_M450]
text = ['M270', 'M450']

fig, ax = plt.subplots(1,1, figsize = (9, 4))


for i in range(2):
    df = allData[i]

    minS = df.minStress.values
    maxS = df.maxStress.values
    for k in range(len(minS)):
        countStress[i] += ((rangeStress > minS[k]) & (rangeStress < maxS[k]))
    
    ax.plot(rangeStress, countStress[i], label = text[i] + ' - total Ncomp = {:.0f}'.format(df.shape[0]))

ax.legend()
ax.set_xlim([0, 1000])
ax.set_xlabel('stress value (Pa)')
ax.set_ylabel('# compressions')
fig.suptitle('Number of compression which stress interval\ncontains a given stress value')
plt.show()


# %%%%%


data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True),
          (data['date'].apply(lambda x : x in ['21-12-08', '22-01-12']))] #, '21-12-16'

data_f = data
for fltr in Filters:
    data_f = data_f.loc[fltr]

binSize = 50
rangeStressLower = np.arange(0, 2000-binSize, binSize)
rangeStressUpper = np.arange(binSize, 2000, binSize)
countStress = [np.zeros(2000//binSize - 1), np.zeros(2000//binSize - 1)]

data_M270 = data_f.loc[data_f['bead type'] == 'M270']
data_M450 = data_f.loc[data_f['bead type'] == 'M450']

allData = [data_M270, data_M450]
text = ['M270', 'M450']

fig, ax = plt.subplots(1,1, figsize = (9, 6))


for i in range(2):
    df = allData[i]

    minS = df.minStress.values
    maxS = df.maxStress.values
    for k in range(1, len(minS)):
        countStress[i] += ((rangeStressLower > minS[k]) & (rangeStressUpper < maxS[k]))
    
    binPos = (rangeStressLower + rangeStressUpper)/2
    ax.plot(binPos, countStress[i], label = text[i] + ' - total Ncomp = {:.0f}'.format(df.shape[0]))

ax.legend()
ax.set_xlim([0, 800])
locator = matplotlib.ticker.MultipleLocator(binSize)
ax.xaxis.set_major_locator(locator)
ax.set_xlabel('stress value (Pa)')
ax.set_ylabel('# compressions')
locator = matplotlib.ticker.MultipleLocator(5)
ax.yaxis.set_major_locator(locator)
ax.yaxis.grid(True)
fig.suptitle('Number of compression which stress interval contains\n a given stress bin of width = {:.0f} Pa'.format(binSize))
plt.show()


# %%%%%


# ['21-12-08', '21-12-16', '22-01-12']
data = GlobalTable_meca_Py2
beadTypes = ['M270', 'M450']
dates = ['22-01-12']

fig = plt.figure(figsize = (9, 10), tight_layout=True)
gs = GridSpec(2, 2, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])
ax = [ax1, ax2, ax3]
colors = ['c', 'orange', 'b', 'r']

for i in range(len(beadTypes)):
    bT = beadTypes[i]
    filterList = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True),\
           (GlobalTable_meca_Py2['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py2['bead type'] == bT), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in dates))]
    globalFilter = filterList[0]
    
    for ii in range(1, len(filterList)):
        globalFilter = globalFilter & filterList[ii]
        
    data[globalFilter].plot(kind = 'scatter', ax = ax[i], x = 'surroundingThickness', y = 'minStress', 
                            marker = 'o', color = 'c', alpha = 0.3, label = 'Min stress')
    data[globalFilter].plot(kind = 'scatter', ax = ax[i], x = 'surroundingThickness', y = 'maxStress', 
                            marker = 'o', color = 'r', alpha = 0.3, label = 'Max stress')
    ax[i].set_xlabel('Cortical thickness (nm)')
    ax[i].set_ylabel('Extremal stress values (Pa)')
    ax[i].set_ylim([0, 2000])
    ax[i].set_xlim([0, 1200])
    ax[i].set_title(bT)
    data[globalFilter].plot(kind = 'scatter', ax = ax[2], x = 'surroundingThickness', y = 'minStress', 
                            marker = 'o', color = colors[2*i], alpha = 0.3, label = bT + ' - Min stress')
    data[globalFilter].plot(kind = 'scatter', ax = ax[2], x = 'surroundingThickness', y = 'maxStress', 
                            marker = 'o', color = colors[2*i+1], alpha = 0.3, label = bT + ' - Max stress')
    
ax[2].set_xlabel('Cortical thickness (nm)')
ax[2].set_xlim([0, 1200])

ax[2].set_ylabel('Extremal stress values (Pa)')
ax[2].set_ylim([0, 1000])
locator = matplotlib.ticker.MultipleLocator(100)
ax[2].yaxis.set_major_locator(locator)
ax[2].yaxis.grid(True)

ax[2].set_title('Jan 22 - Extremal stress values vs. thickness, for each compression')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_StressRanges', figSubDir = figSubDir)
plt.show()


# %%%%%


# ['21-12-08', '21-12-16', '22-01-12']
data = GlobalTable_meca_Py2
beadTypes = ['M270', 'M450']
dates = ['21-12-08', '22-01-12']

fig = plt.figure(figsize = (9, 10), tight_layout=True)
gs = GridSpec(2, 2, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])
ax = [ax1, ax2, ax3]
colors = ['c', 'orange', 'b', 'r']

for i in range(len(beadTypes)):
    bT = beadTypes[i]
    filterList = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True),\
           (GlobalTable_meca_Py2['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py2['bead type'] == bT), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in dates))]
    globalFilter = filterList[0]
    
    for ii in range(1, len(filterList)):
        globalFilter = globalFilter & filterList[ii]
        
    data[globalFilter].plot(kind = 'scatter', ax = ax[i], x = 'surroundingThickness', y = 'minStress', 
                            marker = 'o', color = 'c', alpha = 0.3, label = 'Min stress')
    data[globalFilter].plot(kind = 'scatter', ax = ax[i], x = 'surroundingThickness', y = 'maxStress', 
                            marker = 'o', color = 'r', alpha = 0.3, label = 'Max stress')
    ax[i].set_xlabel('Cortical thickness (nm)')
    ax[i].set_ylabel('Extremal stress values (Pa)')
    ax[i].set_ylim([0, 2000])
    ax[i].set_xlim([0, 1200])
    ax[i].set_title(bT)
    data[globalFilter].plot(kind = 'scatter', ax = ax[2], x = 'surroundingThickness', y = 'minStress', 
                            marker = 'o', color = colors[2*i], alpha = 0.3, label = bT + ' - Min stress')
    data[globalFilter].plot(kind = 'scatter', ax = ax[2], x = 'surroundingThickness', y = 'maxStress', 
                            marker = 'o', color = colors[2*i+1], alpha = 0.3, label = bT + ' - Max stress')
    
ax[2].set_xlabel('Cortical thickness (nm)')
ax[2].set_xlim([0, 1200])

ax[2].set_ylabel('Extremal stress values (Pa)')
ax[2].set_ylim([0, 1000])
locator = matplotlib.ticker.MultipleLocator(100)
ax[2].yaxis.set_major_locator(locator)
ax[2].yaxis.grid(True)

ax[2].set_title('Dec 21 - Extremal stress values vs. thickness, for each compression')

plt.show()


# %%%% Dec-2021 & Jan 2022 -> Non linearity across bead sizes

# %%%%%


# Make the dataframe

data = GlobalTable_meca_nonLin

dates = ['21-12-08', '22-01-12']  # ['22-02-09'] #['21-12-08', '22-01-12'] ['21-01-18', '21-01-21', '21-12-08']

filterList = [(data['validatedThickness'] == True),
              (data['cell subtype'] == 'aSFL'), 
#               (data['bead type'] == 'M450'),
              (data['substrate'] == '20um fibronectin discs'),
              (data['date'].apply(lambda x : x in dates))]  # (data['validatedFit'] == True), 

globalFilter = filterList[0]
for k in range(1, len(filterList)):
    globalFilter = globalFilter & filterList[k]

data_f = data[globalFilter]



fitMin = [S for S in range(25,1225,50)]
fitMax = [S+150 for S in fitMin]
fitCenters = np.array([S+75 for S in fitMin])
regionFitsNames = [str(fitMin[ii]) + '<s<' + str(fitMax[ii]) for ii in range(len(fitMin))]

listColumnsMeca = []

KChadwick_Cols = []
KWeight_Cols = []

for rFN in regionFitsNames:
    listColumnsMeca += ['KChadwick_'+rFN, 'K_CIW_'+rFN, 'R2Chadwick_'+rFN, 'K2Chadwick_'+rFN, 
                        'H0Chadwick_'+rFN, 'Npts_'+rFN, 'validatedFit_'+rFN]
    KChadwick_Cols += [('KChadwick_'+rFN)]

    K_CIWidth = data_f['K_CIW_'+rFN] #.apply(lambda x : x.strip('][').split(', ')).apply(lambda x : (np.abs(float(x[0]) - float(x[1]))))
    KWeight = (data_f['KChadwick_'+rFN]/K_CIWidth)**2
    data_f['K_Weight_'+rFN] = KWeight
    data_f['K_Weight_'+rFN] *= data_f['KChadwick_'+rFN].apply(lambda x : (x<1e6))
    data_f['K_Weight_'+rFN] *= data_f['R2Chadwick_'+rFN].apply(lambda x : (x>1e-2))
    data_f['K_Weight_'+rFN] *= data_f['K_CIW_'+rFN].apply(lambda x : (x!=0))
    KWeight_Cols += [('K_Weight_'+rFN)]
    
data_f.tail()


# %%%%%


listBT = ['M270', 'M450']
data_f2 = [data_f.loc[data_f['bead type'] == BT] for BT in listBT]
fig, ax = plt.subplots(1,1, figsize = (9,6))
color = ['skyblue', 'salmon']

for k in range(len(listBT)):
    BT = listBT[k]
    CondCol = 'date'
    Variables = KChadwick_Cols
    WeightCols = KWeight_Cols
    data_f_agg = getAggDf_weightedAvg(data_f2[k], 'cellID', CondCol, Variables, WeightCols)
    data_f_agg.T
    dictPlot = {'cellID' : [], 'Kavg' : [], 'Kstd' : [], 'Kcount' : []}
    meanCols = []
    stdCols = []
    countCols = []
    for col in data_f_agg.columns:
        if 'KChadwick' in col and 'Wmean' in col:
            meanCols.append(col)
        elif 'KChadwick' in col and '_Wstd' in col:
            stdCols.append(col)
        elif 'KChadwick' in col and '_Count' in col:
            countCols.append(col)
    meanDf = data_f_agg[meanCols]
    stdDf = data_f_agg[stdCols]
    countDf = data_f_agg[countCols]
    for c in data_f_agg.index:
        means = meanDf.T[c].values
        stds = stdDf.T[c].values
        counts = countDf.T[c].values
        dictPlot['cellID'].append(c)
        dictPlot['Kavg'].append(means)
        dictPlot['Kstd'].append(stds)
        dictPlot['Kcount'].append(counts)



    for i in range(len(dictPlot['cellID'])):
        c = dictPlot['cellID'][i]
#         color = colorList10[i%10]
    #     ax.errorbar(fitCenters, dictPlot['Kavg'][i], yerr = dictPlot['Kstd'][i], color = color)
        if i == 0:
            ax.plot(fitCenters, dictPlot['Kavg'][i], color = color[k], label = BT)
        else:
            ax.plot(fitCenters, dictPlot['Kavg'][i], color = color[k])
        low =  dictPlot['Kavg'][i] - (dictPlot['Kstd'][i] / (dictPlot['Kcount'][i]**0.5)) 
        high = dictPlot['Kavg'][i] + (dictPlot['Kstd'][i] / (dictPlot['Kcount'][i]**0.5))
        low = np.where(low < 10, 10, low)
        ax.fill_between(x=fitCenters, 
                           y1=low, 
                           y2=high,
                           color = color[k], alpha = 0.1, zorder = 1)

    ax.set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
    ax.set_ylabel('K (Pa) [tgt modulus w/ Chadwick]')
    ax.set_yscale('log')
    ax.set_ylim([1e2, 1e5])
    ax.legend(loc='lower right')
    fig.suptitle('K(sigma) - on each cell\n21-12 & 22-01 experiments')
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_K(s)allCells', figSubDir = 'NonLin')
plt.show()


# %%%%%


# print(data_f.head())

# print(fitCenters)

def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

listBT = ['M270', 'M450']
data_f2 = [data_f.loc[data_f['bead type'] == BT] for BT in listBT]
fig, ax = plt.subplots(1,1, figsize = (9,6))
color = ['skyblue', 'salmon']
ecolor = ['b', 'r']

for k in range(len(listBT)):

    Kavg = []
    Kstd = []
    D10 = []
    D90 = []
    N = []
    
    fitCenters2 = np.array([S for S in range(100,950,50)])
    
    for S in range(100,950,50):
        rFN = str(S-75) + '<s<' + str(S+75)
        variable = valStr+rFN
        weight = weightStr+rFN

        x = data_f2[k][variable].apply(nan2zero).values
        w = data_f2[k][weight].apply(nan2zero).values
#         print(S)
#         print(x)
#         print(w)

        if S == 250:
            d = {'x' : x, 'w' : w}

        m = np.average(x, weights=w)
        v = np.average((x-m)**2, weights=w)
        std = v**0.5

        d10, d90 = np.percentile(x[x != 0], (10, 90))
        n = len(x[x != 0])

        Kavg.append(m)
        Kstd.append(std)
        D10.append(d10)
        D90.append(d90)
        N.append(n)

    Kavg = np.array(Kavg)
    Kstd = np.array(Kstd)
    D10 = np.array(D10)
    D90 = np.array(D90)
    N = np.array(N)
    Kste = Kstd / (N**0.5)

    alpha = 0.975
    dof = N
    q = st.t.ppf(alpha, dof) # Student coefficient

    d_val = {'S' : fitCenters2, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

    ax.errorbar(fitCenters2, Kavg, yerr = q*Kste, marker = 'o', color = color[k], 
                   ecolor = ecolor[k], elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error\n'+listBT[k])
    ax.set_ylim([50,1e5])

#     for k in range(1):
    ax.legend(loc = 'upper left')
    ax.set_yscale('log')
    ax.set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
    ax.set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
    for kk in range(len(N)):
        ax.text(x=fitCenters2[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6, color = ecolor[k])

# try:
alpha = 0.975
dof = N_ref
q_ref = st.t.ppf(alpha, dof)
ax.errorbar(fitCenters_ref, Kavg_ref, yerr = q_ref*Kste_ref, marker = 'o', color = 'lightgreen', 
               ecolor = 'green', elinewidth = 0.8, capsize = 3, zorder = 1,
            label = 'Weighted means\nWeighted ste 95% as error\nRefExp - M450')
# for kk in range(len(N_ref)):
#         ax.text(x=fitCenters_ref[kk]+5, y=Kavg_ref[kk]**0.98, s='n='+str(N_ref[kk]), fontsize = 6, color = 'green')        
    # Kavg_ref, Kste_ref, N_ref, fitCenters_ref
# except:
#     pass

fig.suptitle('K(sigma) - 2 types of beads experiments')
ax.legend(loc = 'lower right')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_K(s)globalAvg', figSubDir = 'NonLin')
plt.show()

df_val = pd.DataFrame(d_val)
dftest = pd.DataFrame(d)

df_val


# %%%%%


# # print(data_f.head())

# # print(fitCenters)

# def w_std(x, w):
#     m = np.average(x, weights=w)
#     v = np.average((x-m)**2, weights=w)
#     std = v**0.5
#     return(std)

# def nan2zero(x):
#     if np.isnan(x):
#         return(0)
#     else:
#         return(x)

# valStr = 'KChadwick_'
# weightStr = 'K_Weight_'



# Kavg = []
# Kstd = []
# D10 = []
# D90 = []
# N = []

# for S in fitCenters:
#     rFN = str(S-75) + '<s<' + str(S+75)
#     variable = valStr+rFN
#     weight = weightStr+rFN
    
#     x = data_f[variable].apply(nan2zero).values
#     w = data_f[weight].apply(nan2zero).values
    
#     if S == 250:
#         d = {'x' : x, 'w' : w}
    
#     m = np.average(x, weights=w)
#     v = np.average((x-m)**2, weights=w)
#     std = v**0.5
    
#     d10, d90 = np.percentile(x[x != 0], (10, 90))
#     n = len(x[x != 0])
    
#     Kavg.append(m)
#     Kstd.append(std)
#     D10.append(d10)
#     D90.append(d90)
#     N.append(n)

# Kavg = np.array(Kavg)
# Kstd = np.array(Kstd)
# D10 = np.array(D10)
# D90 = np.array(D90)
# N = np.array(N)
# Kste = Kstd / (N**0.5)

# alpha = 0.975
# dof = N
# q = st.t.ppf(alpha, dof) # Student coefficient

# d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

# fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# # ax[0]
# ax.errorbar(fitCenters, Kavg, yerr = q*Kste, marker = 'o', color = my_default_color_list[0], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error')
# ax.set_ylim([500,2e4])

# # ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
# #                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# # ax[1].set_ylim([500,1e6])

# # for k in range(1): #2
# #     ax[k].legend(loc = 'upper left')
# #     ax[k].set_yscale('log')
# #     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
# #     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
# #     for kk in range(len(N)):
# #         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
# ax.legend(loc = 'upper left')
# ax.set_yscale('log')
# ax.set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
# ax.set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
# for kk in range(len(N)):
#     ax.text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

# fig.suptitle('K(sigma)')
# ax.set_title('From all compressions pooled\n22-02-09 experiment, 36 cells, 232 compression')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_K(s)globalAvg_V2', figSubDir = 'NonLin')
# plt.show()

# # df_val = pd.DataFrame(d_val)
# # dftest = pd.DataFrame(d)

# # df_val


# %%%%%


# print(data_f.head())

# print(fitCenters)

listBT = ['M270', 'M450']
data_f2 = [data_f.loc[data_f['bead type'] == BT] for BT in listBT]
fig, ax = plt.subplots(1,1, figsize = (9,6))
color = ['skyblue', 'salmon']
ecolor = ['b', 'r']

mini, maxi = 200, 450

for k in range(len(listBT)):
    
    extraFilters = [data_f2[k]['minStress'] <= mini, data_f2[k]['maxStress'] >= maxi]
    fitCenters2 = fitCenters[(fitCenters>=mini) & (fitCenters<=maxi)]

    data_ff = data_f2[k]
    globalExtraFilter = extraFilters[0]
    for kk in range(1, len(extraFilters)):
        globalExtraFilter = globalExtraFilter & extraFilters[kk]

    data_ff = data_ff[globalExtraFilter]
    data_ff
    def w_std(x, w):
        m = np.average(x, weights=w)
        v = np.average((x-m)**2, weights=w)
        std = v**0.5
        return(std)

    def nan2zero(x):
        if np.isnan(x):
            return(0)
        else:
            return(x)

    valStr = 'KChadwick_'
    weightStr = 'K_Weight_'

    Kavg = []
    Kstd = []
    D10 = []
    D90 = []
    N = []

    for S in fitCenters2:
        rFN = str(S-75) + '<s<' + str(S+75)
        variable = valStr+rFN
        weight = weightStr+rFN

        x = data_ff[variable].apply(nan2zero).values
        w = data_ff[weight].apply(nan2zero).values

        if S == 250:
            d = {'x' : x, 'w' : w}

        m = np.average(x, weights=w)
        v = np.average((x-m)**2, weights=w)
        std = v**0.5

        d10, d90 = np.percentile(x[x != 0], (10, 90))
        n = len(x[x != 0])

        Kavg.append(m)
        Kstd.append(std)
        D10.append(d10)
        D90.append(d90)
        N.append(n)

    Kavg = np.array(Kavg)
    Kstd = np.array(Kstd)
    D10 = np.array(D10)
    D90 = np.array(D90)
    N = np.array(N)
    Kste = Kstd / (N**0.5)

    alpha = 0.975
    dof = N
    q = st.t.ppf(alpha, dof) # Student coefficient

    d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

    # ax[0]
    ax.errorbar(fitCenters2, Kavg, yerr = q*Kste, marker = 'o', color = color[k], 
                   ecolor = ecolor[k], elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error\n'+listBT[k])
    ax.set_ylim([500,2e4])

    # ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
    #                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
    # ax[1].set_ylim([500,1e6])

    # for k in range(1): #2
    #     ax[k].legend(loc = 'upper left')
    #     ax[k].set_yscale('log')
    #     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
    #     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
    #     for kk in range(len(N)):
    #         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
    
    ax.set_yscale('log')
    ax.set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
    ax.set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
    for kk in range(len(N)):
        ax.text(x=fitCenters2[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6, color = ecolor[k])

alpha = 0.975
dof = N_ref
q_ref = st.t.ppf(alpha, dof)
ax.errorbar(fitCenters_ref, Kavg_ref, yerr = q_ref*Kste_ref, marker = 'o', color = 'lightgreen', 
               ecolor = 'green', elinewidth = 0.8, capsize = 3, zorder = 1,
            label = 'Weighted means\nWeighted ste 95% as error\nRefExp - M450')
        
fig.suptitle('K(sigma)')
ax.set_title('Only compressions including the [{:.0f}, {:.0f}]Pa range'.format(mini, maxi))
ax.legend(loc = 'lower right')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_K(s)globalAvg_100to800', figSubDir = 'NonLin')
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


# %%%%%


# print(data_f.head())

# print(fitCenters)

extraFilters = [data_f['minStress'] <= 100, data_f['maxStress'] >= 600]
fitCenters2 = fitCenters[fitCenters<650]

data_ff = data_f
globalExtraFilter = extraFilters[0]
for k in range(1, len(extraFilters)):
    globalExtraFilter = globalExtraFilter & extraFilters[k]

data_ff = data_ff[globalExtraFilter]
data_ff
def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in fitCenters2:
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_ff[variable].apply(nan2zero).values
    w = data_ff[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# ax[0]
ax.errorbar(fitCenters2, Kavg, yerr = q*Kste, marker = 'o', color = my_default_color_list[0], 
               ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error')
ax.set_ylim([500,2e4])

# ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# ax[1].set_ylim([500,1e6])

# for k in range(1): #2
#     ax[k].legend(loc = 'upper left')
#     ax[k].set_yscale('log')
#     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
#     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
#     for kk in range(len(N)):
#         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
ax.legend(loc = 'upper left')
ax.set_yscale('log')
ax.set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
ax.set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
for kk in range(len(N)):
    ax.text(x=fitCenters2[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

fig.suptitle('K(sigma)')
ax.set_title('Only compressions including the [100, 600]Pa range')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_K(s)globalAvg_100to600', figSubDir = 'NonLin')
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


# %%%%%


data_f


# %%%%%


warnings.filterwarnings('ignore')

data = data_f

listeS = [100, 150, 200, 250, 300, 400, 500, 600, 700, 800]

fig, axes = plt.subplots(1, len(listeS), figsize = (14,6))
kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in ['22-02-09']))] #, '21-12-16'

    fig, ax = D1Plot(data, fig=fig, ax=axes[kk], CondCol=['bead type'], Parameters=['KChadwick_'+interval], Filters=Filters, 
                   Boxplot=True, cellID='cellID', co_order=[], stats=True, statMethod='Mann-Whitney', AvgPerCell = True,
                   box_pairs=[], figSizeFactor = 1, markersizeFactor=1, orientation = 'h')

#     axes[kk].legend(loc = 'upper right', fontsize = 8)
    label = axes[kk].get_ylabel()
    axes[kk].set_title(label.split('_')[-1] + 'Pa', fontsize = 10)
    axes[kk].set_yscale('log')
    axes[kk].set_ylim([1e2, 5e4])
    axes[kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
    axes[kk].set_xticklabels([])
    if kk == 0:
        axes[kk].set_ylabel(label.split('_')[0] + ' (Pa)')
        axes[kk].tick_params(axis='x', labelsize = 10)
    else:
        axes[kk].set_ylabel(None)
        axes[kk].set_yticklabels([])
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Detailed1DPlot', figSubDir = figSubDir)
    
    kk+=1

renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
fig.tight_layout()
fig.suptitle('Tangeantial modulus of aSFL 3T3 - 0.5>50mT lowStartCompressions - avg per cell')
plt.show()
jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_nonLinK_multipleIntervals', figSubDir = 'NonLin')
    
warnings.filterwarnings('always')


# %%%%%


warnings.filterwarnings('ignore')

data = GlobalTable_meca_Py2

listeS = [100, 150, 200, 250, 300, 400, 500, 600, 700, 800]

# fig, axes = plt.subplots(len(listeS), 1, figsize = (9,40))
# kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)
    textForFileName = str(S-75) + '-' + str(S+75) + 'Pa'

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['bestH0'] <= 600),
               (data['date'].apply(lambda x : x in ['22-02-09']))] #, '21-12-16'

    fig, ax = D2Plot_wFit(data, #fig=fig, ax=axes[kk], 
                          XCol='bestH0', YCol='KChadwick_'+interval, CondCol = ['bead type'], Filters=Filters, 
                          cellID = 'cellID', AvgPerCell=False, 
                          xscale = 'linear', yscale = 'log', modelFit=True, modelType='y=k*x^a')
    
    # individual figs
    ax.set_xlim([8e1, 700])
    ax.set_ylim([1e2, 5e4])
    renameAxes(ax,{'bestH0':'best H0 (nm)', 'doxycyclin':'iMC'})
    fig.set_size_inches(6,4)
    fig.tight_layout()
#     jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_nonLin_K(h)_' + textForFileName, figSubDir = 'NonLin')
    
    
    # one big fig
#     axes[kk].set_xlim([8e1, 1.2e3])
#     axes[kk].set_ylim([3e2, 5e4])
#     kk+=1


# renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
# fig.tight_layout()
# plt.show()
    
warnings.filterwarnings('always')





# %%%  Plots that explore the impact of the delay before starting experiment

# %%%%%


figSubDir = 'ExploratoryPlots'
GlobalTable_meca_Py['compStartTimeThisDay_min'] = GlobalTable_meca_Py['compStartTimeThisDay']/60
GlobalTable_meca_Py


# %%%%%


figSubDir = 'ExploratoryPlots'
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['manipID','drug'],Parameters=['surroundingThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',stats=False)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL - Compressions - Differences with manip')
for i in range(len(ax)):
    xTicks = ax[i].get_xticklabels()
    for j in range(len(xTicks)):
        xTicks[j].set_fontsize(9)
jvu.archiveFig(fig, ax, name='3T3aSFL-Compressions-DifferencesWithManip_PYTHONTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs')]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='compStartTimeThisDay_min',YCol='EChadwick', CondCol = ['date_x','drug'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, modelFit=False)
renameAxes(ax,renameDict1)
ax.legend(loc='upper right')
ax.set_xlabel('Time since start of experiment (min)')
fig.suptitle('3T3aSFL - Compressions - Evolution with time')
jvu.archiveFig(fig, ax, name='3T3aSFL-Compressions-EvolWithTime_PYTHONTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['cell subtype'] == 'aSFL-A8'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['manipID','drug'],Parameters=['surroundingThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',stats=True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL-A8 - Compressions - Differences with manip')
for i in range(len(ax)):
    xTicks = ax[i].get_xticklabels()
    for j in range(len(xTicks)):
        xTicks[j].set_fontsize(9)
jvu.archiveFig(fig, ax, name='3T3aSFL-A8-Compressions-DifferencesWithManip_PYTHONTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['cell subtype'] == 'aSFL-A8'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs')]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='compStartTimeThisDay_min',YCol='EChadwick', CondCol = ['date_x','drug'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, modelFit=False)
renameAxes(ax,renameDict1)
ax.legend(loc='upper right')
ax.set_xlabel('Time since start of experiment (min)')
fig.suptitle('3T3aSFL-A8 - Compressions - Evolution with time')
jvu.archiveFig(fig, ax, name='3T3aSFL-A8-Compressions-EvolWithTime_PYTHONTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['cell subtype'] == 'aSFL-6FP'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs')]
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['manipID','drug'],Parameters=['surroundingThickness','EChadwick'],                 Filters=Filters,AvgPerCell=True,cellID='cellID',stats=True)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL-6FP - Compressions - Differences with manip')
for i in range(len(ax)):
    xTicks = ax[i].get_xticklabels()
    for j in range(len(xTicks)):
        xTicks[j].set_fontsize(9)
jvu.archiveFig(fig, ax, name='3T3aSFL-6FP-Compressions-DifferencesWithManip_PYTHONTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),           (GlobalTable_meca_Py['cell subtype'] == 'aSFL-6FP'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs')]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='compStartTimeThisDay_min',YCol='EChadwick', CondCol = ['date_x','drug'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, modelFit=False)
renameAxes(ax,renameDict1)
ax.legend(loc='upper right')
ax.set_xlabel('Time since start of experiment (min)')
fig.suptitle('3T3aSFL-6FP - Compressions - Evolution with time')
jvu.archiveFig(fig, ax, name='3T3aSFL-6FP-Compressions-EvolWithTime_PYTHONTABLE', figSubDir = figSubDir)
plt.show()


# %%%  Compare different data from python and matlab code

# %%%%%


Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), 
           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), GlobalTable_meca_Py['date_x'].apply(lambda x : x in ['21-01-18', '21-01-21'])]
co_order = makeOrder(['20um fibronectin discs'],['none','doxycyclin'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['substrate','drug'],                 Parameters=['surroundingThickness','EChadwick'],Filters=Filters,                 AvgPerCell = True, cellID='cellID', co_order=co_order, figSizeFactor = 0.8,
                 useHue = False, orientation = 'h')
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL on diverse substrates: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFL_substrate&drug_SurroundingThickness&EChadwick_NEWTABLE', figSubDir = figSubDir)
plt.show()


# %%%%%


Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py2['date_x'].apply(lambda x : x in ['21-01-18', '21-01-21']))]
co_order = makeOrder(['20um fibronectin discs'],['none','doxycyclin'])
fig, ax = D1Plot(GlobalTable_meca_Py2, CondCol=['substrate','drug'],                 Parameters=['surroundingThickness','EChadwick'],Filters=Filters,                 AvgPerCell = True, cellID='cellID', co_order=co_order, figSizeFactor = 0.8,
                 useHue = False, orientation = 'h')
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL on diverse substrates: Compressions')
# jvu.archiveFig(fig, ax, name='3T3aSFL_substrate&drug_SurroundingThickness&EChadwick_NEWTABLE', figSubDir = figSubDir)
plt.show()


# %%%%% TableV2


# %%%%%


FiltersPy = [(GlobalTable_meca_Py['validatedFit'] == True), 
           (GlobalTable_meca_Py['validatedThickness'] == True), 
           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py['date_x'].apply(lambda x : x in ['21-01-18', '21-01-21']))]
TableV1 = GlobalTable_meca_Py
for fltr in FiltersPy:
    TableV1 = TableV1.loc[fltr]
    
group = TableV1.groupby('cellID')
dictAggMean = getDictAggMean(TableV1)
TableV1 = group.agg(dictAggMean)

FiltersPy2 = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['cell subtype'] == 'aSFL'), 
           (GlobalTable_meca_Py2['date_x'].apply(lambda x : x in ['21-01-18', '21-01-21']))]
TableV2 = GlobalTable_meca_Py2
for fltr in FiltersPy2:
    TableV2 = TableV2.loc[fltr]
    
group = TableV2.groupby('cellID')
dictAggMean = getDictAggMean(TableV2)
TableV2 = group.agg(dictAggMean)
    
cellID = TableV2['cellID'].values

plotDict = {'cellID':[],'manipID':[],'drug':[],'E1':[],'h1':[],'E2':[],'h2':[]}
strE = 'EChadwick'
strH = 'H0Chadwick' #'surroundingThickness'
for iD in cellID[:]:
    try:
        plotDict['cellID'].append(iD)
        plotDict['manipID'] = plotDict['manipID'] + [TableV1.loc[TableV1['cellID'] == iD, 'manipID'].values[0]]
        plotDict['drug'] = plotDict['drug'] + [TableV1.loc[TableV1['cellID'] == iD, 'drug'].values[0]]
        plotDict['E1'] = plotDict['E1'] + [TableV1.loc[TableV1['cellID'] == iD, strE].values[0]]
        plotDict['h1'] = plotDict['h1'] + [TableV1.loc[TableV1['cellID'] == iD, strH].values[0]]
        plotDict['E2'] = plotDict['E2'] + [TableV2.loc[TableV2['cellID'] == iD, strE].values[0]]
        plotDict['h2'] = plotDict['h2'] + [TableV2.loc[TableV2['cellID'] == iD, strH].values[0]]
    except:
        print(iD + ' ignored')

plotDf = pd.DataFrame(plotDict)

plotDf['h2/h1'] = plotDf['h2'] / plotDf['h1']
plotDf['E2/E1'] = plotDf['E2'] / plotDf['E1']

fig, ax = plt.subplots(2,2, figsize = (9,9))
dictColor = {'doxycyclin' : 'red', 'none' : 'cyan'}

cyan_line = matplotlib.lines.Line2D([], [], color='cyan', label='Control')
red_line = matplotlib.lines.Line2D([], [], color='red', label='Doxycyclin')

for i in range(plotDf.shape[0]):
    ax[0,0].plot([1, 2], [plotDf['h1'].values[i], plotDf['h2'].values[i]], 
               marker = 'o', color = dictColor[plotDf['drug'].values[i]],
                markeredgecolor='k', linewidth = 1)
    ax[0,0].set_ylabel('h (nm)')
    
    ax[0,1].plot([1, 2], [plotDf['E1'].values[i], plotDf['E2'].values[i]], 
               marker = 'o', color = dictColor[plotDf['drug'].values[i]],
                markeredgecolor='k', linewidth = 1)
    ax[0,1].set_ylabel('E (Pa)')


sns.swarmplot(y='h2/h1', data=plotDf, ax=ax[1,0], color = 'orange', edgecolor='k', linewidth = 1)
sns.swarmplot(y='E2/E1', data=plotDf, ax=ax[1,1], color = 'green', edgecolor='k', linewidth = 1)
ax[1,0].plot(ax[1,0].get_xlim(), [1,1], 'k--')
ax[1,1].plot(ax[1,1].get_xlim(), [1,1], 'k--')
    
for k in range(2):
    ax[0,k].set_xlim([0.6,3.4])
    ax[0,k].set_xticks([1, 2])
    ax[0,k].set_xticklabels(['Matlab', 'Python'])
    ax[0,k].legend(handles=[cyan_line, red_line])
    
#     x = [1, 2]
#     labels = ['Matlab', 'Python']

    
fig.suptitle('2021-01 experiments analysed with the two algos')

plotDf





# %%%  Region fits -- First attempt (deprecated)

# %%%%%


figSubDir = 'RegionFits'


# %%%%%


Outliers = ['21-07-08_M3_P1_C10', '21-07-08_M3_P1_C5', '21-07-08_M3_P1_C6']
Filters = [(GlobalTable_meca_Py['validatedFit_f<100pN'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), 
           (GlobalTable_meca_Py['drug'] == 'none'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers)),
           (GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270', 'M270_M450', 'M450_M270', 'M450']))]
co_order = makeOrder(['M270', 'M270_M450', 'M450_M270', 'M450'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['bead type'],                 Parameters=['surroundingThickness','EChadwick_f<100pN'],Filters=Filters,                 AvgPerCell=True, cellID='cellID', co_order=co_order, figSizeFactor = 2, orientation = 'v', useHue = False)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL asym bead pairs - all comps - f<100pN')
# jvu.archiveFig(fig, ax, name='3T3aSFL_asymBeadsAllComp_SurroundingThickness&EChadwick', figSubDir = figSubDir)
plt.show()


# %%%%%


Outliers = ['21-07-08_M3_P1_C10', '21-07-08_M3_P1_C5', '21-07-08_M3_P1_C6']
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), 
           (GlobalTable_meca_Py['drug'] == 'none'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers)),
           (GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270', 'M270_M450', 'M450_M270', 'M450']))]
co_order = makeOrder(['M270', 'M270_M450', 'M450_M270', 'M450'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['bead type'],                 Parameters=['surroundingThickness','EChadwick'],Filters=Filters,                 AvgPerCell=True, cellID='cellID', co_order=co_order, figSizeFactor = 2, orientation = 'v', useHue = False)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL asym bead pairs')
jvu.archiveFig(fig, ax, name='3T3aSFL_asymBeads_SurroundingThickness&EChadwick_allDates', figSubDir = figSubDir)
plt.show()


# %%%%%


Outliers = ['21-07-08_M3_P1_C10', '21-07-08_M3_P1_C5', '21-07-08_M3_P1_C6']
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), 
           (GlobalTable_meca_Py['date_x'] == '21-07-08'),
           (GlobalTable_meca_Py['drug'] == 'none'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers)),
           (GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270', 'M270_M450', 'M450_M270', 'M450']))]
co_order = makeOrder(['M270', 'M270_M450', 'M450_M270', 'M450'])
fig, ax = D1Plot(GlobalTable_meca_Py, CondCol=['bead type'],                 Parameters=['surroundingThickness','EChadwick'],Filters=Filters,                 AvgPerCell=True, cellID='cellID', co_order=co_order, figSizeFactor = 2, orientation = 'v', useHue = False)
renameAxes(ax,renameDict1)
fig.suptitle('3T3aSFL asym bead pairs')
jvu.archiveFig(fig, ax, name='3T3aSFL_asymBeads_SurroundingThickness&EChadwick_21-07-08', figSubDir = figSubDir)
plt.show()


# %%%%%


Outliers = ['21-07-08_M3_P1_C10', '21-07-08_M3_P1_C5', '21-07-08_M3_P1_C6']
Filters = [(GlobalTable_meca_Py['validatedFit_f<100pN'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), 
           (GlobalTable_meca_Py['drug'] == 'none'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers)),
           (GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270', 'M450']))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='surroundingThickness',YCol='EChadwick_f<100pN', CondCol = ['bead type'],           Filters=Filters, cellID = 'cellID', AvgPerCell=True, modelFit=False, modelType='y=A*exp(kx)')
fig.suptitle('3T3aSFL: E(h)')
# jvu.archiveFig(fig, ax, name='aSFL_E(h)_drug&substrate', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir='')
plt.show()


# %%%%%


Outliers = ['21-07-08_M3_P1_C10', '21-07-08_M3_P1_C5', '21-07-08_M3_P1_C6']
Filters = [(GlobalTable_meca_Py['validatedFit_f<150pN'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),
           (GlobalTable_meca_Py['date_x'].apply(lambda x : x in ['21-07-08'])), 
           (GlobalTable_meca_Py['drug'] == 'none'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers)),
           (GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270', 'M450_M270', 'M270_M450', 'M450']))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='surroundingThickness',YCol='EChadwick_f<150pN', CondCol = ['bead type'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, modelFit=False, modelType='y=A*exp(kx)')
fig.suptitle('3T3aSFL: E(h) _f<150pN_21-07-08')
jvu.archiveFig(fig, ax, name='aSFL_beadTypes01_E(h)_21-07-08', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


Outliers = ['21-07-08_M3_P1_C10', '21-07-08_M3_P1_C5', '21-07-08_M3_P1_C6']
Filters = [(GlobalTable_meca_Py['validatedFit_f<150pN'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),
           (GlobalTable_meca_Py['date_x'].apply(lambda x : x in ['21-07-08'])), 
           (GlobalTable_meca_Py['drug'] == 'none'), (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py['cell subtype'] == 'aSFL'), (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers)),
           (GlobalTable_meca_Py['bead type'].apply(lambda x : x in ['M270', 'M450']))]
fig, ax = D2Plot(GlobalTable_meca_Py, XCol='surroundingThickness',YCol='EChadwick_f<150pN', CondCol = ['bead type'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, modelFit=False, modelType='y=A*exp(kx)')
fig.suptitle('3T3aSFL: E(h) _f<150pN_21-07-08')
jvu.archiveFig(fig, ax, name='aSFL_beadTypes02_E(h)_21-07-08', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%  Region fits, non linearity


# %%%% On Jan 2021 data

# %%%%%


data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['drug'] == 'none'),
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in ['21-01-18', '21-01-21']))]

data_f = data
for fltr in Filters:
    data_f = data_f.loc[fltr]
    
print(data_f[data_f['validatedFit_s<500Pa']].shape)
print(data_f[data_f['validatedFit_500<s<1000Pa']].shape)
print(data_f[data_f['validatedFit_s<500Pa'] & data_f['validatedFit_500<s<1000Pa']].shape)


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fits = ['s<500Pa', '250<s<750Pa', '500<s<1000Pa']

Filters = [(data['validatedFit'] == True),
           (data['validatedThickness'] == True),
           (data['drug'] == 'none'),
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in dates))]
Filters += [(data['validatedFit_' + f] == True) for f in fits]

Parameters = ['EChadwick_' + f for f in fits]

fig, ax = D1PlotPaired(data, Parameters=Parameters, Filters=Filters, Boxplot=True, cellID='cellID', 
                   co_order=[], stats=True, statMethod='Wilcox_less', box_pairs=[],
                   figSizeFactor = 1.5, markersizeFactor=1, orientation = 'h', labels = fits)
ax.set_ylabel('E_Chadwick (Pa)')
fig.suptitle('Jan 2021 data ; B = 3 -> 40 mT ; no drug')
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fits = ['s<500Pa', '250<s<750Pa', '500<s<1000Pa']

Filters = [(data['validatedFit'] == True),
           (data['validatedThickness'] == True),
           (data['drug'] == 'doxycyclin'),
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in dates))]
Filters += [(data['validatedFit_' + f] == True) for f in fits]

Parameters = ['EChadwick_' + f for f in fits]

fig, ax = D1PlotPaired(data, Parameters=Parameters, Filters=Filters, Boxplot=True, cellID='cellID', 
                   co_order=[], stats=True, statMethod='Wilcox_less', box_pairs=[],
                   figSizeFactor = 1.5, markersizeFactor=1, orientation = 'h', labels = fits)
ax.set_ylabel('E_Chadwick (Pa)')
fig.suptitle('Jan 2021 data ; B = 3 -> 40 mT ; doxy')
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fits = ['s<400Pa', '300<s<700Pa', '600<s<1000Pa']

Filters = [(data['validatedFit'] == True),
           (data['validatedThickness'] == True),
           (data['drug'] == 'none'),
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in dates))]
Filters += [(data['validatedFit_' + f] == True) for f in fits]

Parameters = ['EChadwick_' + f for f in fits]

fig, ax = D1PlotPaired(data, Parameters=Parameters, Filters=Filters, Boxplot=True, cellID='cellID', 
                   co_order=[], stats=True, statMethod='Wilcox_less', box_pairs=[],
                   figSizeFactor = 1.5, markersizeFactor=1, orientation = 'h', labels = fits)
ax.set_ylabel('E_Chadwick (Pa)')
fig.suptitle('Jan 2021 data ; B = 3 -> 40 mT ; no drug')
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-12-08', '22-01-12']
fits = ['s<100Pa', '100<s<200Pa', '200<s<300Pa']

Filters = [(data['validatedFit'] == True),
           (data['validatedThickness'] == True),
           (data['drug'] == 'none'),
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in dates))]
Filters += [(data['validatedFit_' + f] == True) for f in fits]

Parameters = ['EChadwick_' + f for f in fits]

fig, ax = D1PlotPaired(data, Parameters=Parameters, Filters=Filters, Boxplot=True, cellID='cellID', 
                   co_order=[], stats=True, statMethod='Wilcox_less', box_pairs=[],
                   figSizeFactor = 1.5, markersizeFactor=1, orientation = 'h', labels = fits)
ax.set_ylabel('E_Chadwick (Pa)')
fig.suptitle('Jan 2022 data ; B = 1 -> 13 mT ; no drug')
plt.show()


# %%%%%


# ['21-12-08', '21-12-16', '22-01-12']
data = GlobalTable_meca_Py2

dates = ['22-02-09'] #['21-12-08', '22-01-12'] ['21-01-18', '21-01-21', '21-12-08']

filterList = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),\
           (data['cell subtype'] == 'aSFL'), 
           (data['bead type'] == 'M450'),
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in dates))]


fig, ax = plt.subplots(1, 1, figsize = (9, 6), tight_layout=True)

globalFilter = filterList[0]
for k in range(1, len(filterList)):
    globalFilter = globalFilter & filterList[k]

X = data[globalFilter]['surroundingThickness'].values
Smin = data[globalFilter]['minStress'].values
Smax = data[globalFilter]['maxStress'].values

for i in range(len(X)):
    ax.plot([X[i], X[i]], [Smin[i], Smax[i]], 
            ls = '-', color = 'deepskyblue', alpha = 0.3,
            label = 'Stress range', zorder = 1)
    ax.plot([X[i]], [Smin[i]], 
            marker = 'o', markerfacecolor = 'skyblue', markeredgecolor = 'k', markersize = 4,
            ls = '')
    ax.plot([X[i]], [Smax[i]], 
            marker = 'o', markerfacecolor = 'royalblue', markeredgecolor = 'k', markersize = 4,
            ls = '')

ax.set_xlabel('Cortical thickness (nm)')
ax.set_xlim([0, 1200])
locator = matplotlib.ticker.MultipleLocator(100)
ax.xaxis.set_major_locator(locator)

ax.set_ylabel('Extremal stress values (Pa)')
ax.set_ylim([0, 2500])
locator = matplotlib.ticker.MultipleLocator(100)
ax.yaxis.set_major_locator(locator)
ax.yaxis.grid(True)

ax.set_title('Jan 21 & jan 22 - Extremal stress values vs. thickness, for each compression')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_StressRanges', figSubDir = figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fit = 's<200Pa_included' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['validatedFit_'+fit] == True),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp')
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fit = '100<s<300Pa_200included' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['bead type'] == 'M450'),
           (data['validatedFit_'+fit] == True),
           (data['Npts_'+fit] >= 15),
           (data['surroundingThickness'] <= 800),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)', fontsize = 12)
ax.set_xlabel('Thickness at low force (nm)', fontsize = 12)
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp', fontsize = 14)
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fit = '50<s<250Pa_150included' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['bead type'] == 'M450'),
           (data['validatedFit_'+fit] == True),
           (data['Npts_'+fit] >= 15),
           (data['surroundingThickness'] <= 800),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)', fontsize = 12)
ax.set_xlabel('Thickness at low force (nm)', fontsize = 12)
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp', fontsize = 14)
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['22-01-12']
fit = '50<s<250Pa_150included' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['bead type'] == 'M450'),
           (data['validatedFit_'+fit] == True),
           (data['Npts_'+fit] >= 15),
           (data['surroundingThickness'] <= 800),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)', fontsize = 12)
ax.set_xlabel('Thickness at low force (nm)', fontsize = 12)
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp', fontsize = 14)
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21', '22-01-12']
fit = '50<s<250Pa_150included' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['bead type'] == 'M450'),
           (data['validatedFit_'+fit] == True),
           (data['Npts_'+fit] >= 15),
           (data['surroundingThickness'] <= 800),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)', fontsize = 12)
ax.set_xlabel('Thickness at low force (nm)', fontsize = 12)
fig.suptitle('3T3aSFL: E(h) // dates = ' + 'all' + ' // ' + fit + ' // All comp', fontsize = 14)
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fit = 's<300Pa_included' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['validatedFit_'+fit] == True),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp')
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fit = '200<s<500Pa_350included' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['validatedFit_'+fit] == True),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp')
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fit = 's<200Pa' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['validatedFit_'+fit] == True),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp')
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fit = '250<s<750Pa' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['validatedFit_'+fit] == True),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp', fontsize = 17)
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-01-18', '21-01-21']
fit = '500<s<1000Pa' # 's<400Pa', '300<s<700Pa', '600<s<1000Pa'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['validatedFit_'+fit] == True),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick_' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp', fontsize = 17)
ax.legend(loc = 'upper right', fontsize = 8)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-12-08', '21-12-16', '22-01-12']
fit = '' # ['s<100Pa', '100<s<200Pa', '200<s<300Pa']

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True),
           (data['validatedFit'+fit] == True),
           (data['substrate'] == '20um fibronectin discs'),
           (data['drug'] == 'none'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick' + fit,CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
ax.legend(loc = 'upper right', fontsize = 8)
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit + ' // All comp')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['22-01-12']
fit = '_s<400Pa' # ['s<100Pa', '100<s<200Pa', '200<s<300Pa']

Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick' + fit,CondCol = ['bead type'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit[1:] +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')
ax.legend(loc = 'upper right', fontsize = 8)
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit[1:] + ' // All comp')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['22-01-12']
fit = '_s<300Pa' # ['s<100Pa', '100<s<200Pa', '200<s<300Pa']

Filters = [(GlobalTable_meca_Py2['validatedFit'] == True), 
           (GlobalTable_meca_Py2['validatedThickness'] == True), 
           (GlobalTable_meca_Py2['substrate'] == '20um fibronectin discs'),
           (GlobalTable_meca_Py2['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(data, XCol='surroundingThickness',YCol='EChadwick' + fit,CondCol = ['bead type'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit[1:] +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')
ax.legend(loc = 'upper right', fontsize = 8)
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit[1:] + ' // All comp')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-12-08', '21-12-16', '22-01-12']
fit = '' # ['s<100Pa', '100<s<200Pa', '200<s<300Pa']

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(GlobalTable_meca_Py2, XCol='ctFieldThickness',YCol='EChadwick' + fit, CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=True, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit[1:] +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')
ax.legend(loc = 'upper right', fontsize = 8)
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit[1:] + ' // All comp')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)_all3exp', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


data = GlobalTable_meca_Py2
dates = ['21-12-08', '22-01-12']
fit = '_s<400Pa' # ['s<100Pa', '100<s<200Pa', '200<s<300Pa']

# Filters = [(data['validatedFit'] == True), 
#            (data['validatedThickness'] == True),
#            (data['validatedFit_'+fit] == True),
#            (data['substrate'] == '20um fibronectin discs'),
#            (data['drug'] == 'none'),
#            (data['date'].apply(lambda x : x in dates))]

# fig, ax = D2Plot_wFit(data, XCol='surroundingThickness', YCol='EChadwick_' + fit, CondCol = ['bead type', 'date'],\
#            Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')


Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['substrate'] == '20um fibronectin discs'),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D2Plot_wFit(GlobalTable_meca_Py2, XCol='surroundingThickness',YCol='EChadwick' + fit, CondCol = ['bead type', 'date'],           Filters=Filters, cellID = 'cellID', AvgPerCell=False, xscale = 'log', yscale = 'log', modelFit=True, modelType='y=k*x^a')

ax.set_ylabel('EChadwick ['+ fit[1:] +']  (Pa)')
ax.set_xlabel('Thickness at low force (nm)')
fig.suptitle('3T3aSFL: E(h)')
ax.legend(loc = 'upper right', fontsize = 8)
fig.suptitle('3T3aSFL: E(h) // dates = ' + dates[0][:5] + ' // ' + fit[1:] + ' // All comp')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_E(h)_all3exp', figDir = todayFigDirLocal + '//' + figSubDir)
plt.show()


# %%%%%


# plt.close('all')
data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True),
          (data['date'].apply(lambda x : x in ['21-12-08', '22-01-12']))] #, '21-12-16'

fig, ax = D1PlotDetailed(data, CondCol=['bead type'], Parameters=['surroundingThickness', 'EChadwick'], Filters=Filters, 
                Boxplot=True, cellID='cellID', co_order=[], stats=True, statMethod='Mann-Whitney', 
               box_pairs=[], figSizeFactor = 1.8, markersizeFactor=1, orientation = 'v', showManips = True)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Detailed1DPlot', figSubDir = figSubDir)
plt.show()


# %%%%%


plt.close('all')
data = GlobalTable_meca_Py2

Filters = [(data['validatedFit'] == True), (data['validatedThickness'] == True),
          (data['date'].apply(lambda x : x in ['21-12-08', '21-12-16', '22-01-12']))] #, '21-12-16'

box_pairs=[('21-12-08 & M270', '21-12-08 & M450'),
             ('21-12-16 & M450', '21-12-16 & M270'),
             ('22-01-12 & M270', '22-01-12 & M450')]

fig, ax = D1PlotDetailed(data, CondCol=['date', 'bead type'], Parameters=['surroundingThickness', 'EChadwick'], Filters=Filters, 
                Boxplot=True, cellID='cellID', co_order=[], stats=True, statMethod='Mann-Whitney', 
               box_pairs=box_pairs, figSizeFactor = 0.9, markersizeFactor=1, orientation = 'v', showManips = True)
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Detailed1DPlot_all3exp', figSubDir = figSubDir)
plt.show()


# %%%% On Feb 2022 data

# %%%%%


data = GlobalTable_meca_nonLin
data


# %%%%%


# ['21-12-08', '21-12-16', '22-01-12']
data = GlobalTable_meca_nonLin

dates = ['22-02-09'] #['21-12-08', '22-01-12'] ['21-01-18', '21-01-21', '21-12-08']

filterList = [(data['validatedThickness'] == True),
              (data['cell subtype'] == 'aSFL'), 
              (data['bead type'] == 'M450'),
              (data['substrate'] == '20um fibronectin discs'),
              (data['date'].apply(lambda x : x in dates))]  # (data['validatedFit'] == True), 


fig, ax = plt.subplots(1, 1, figsize = (9, 6), tight_layout=True)

globalFilter = filterList[0]
for k in range(1, len(filterList)):
    globalFilter = globalFilter & filterList[k]


X = data[globalFilter]['bestH0'].values
Ncomp = len(X)
Smin = data[globalFilter]['minStress'].values
Smax = data[globalFilter]['maxStress'].values

for i in range(Ncomp):
    ax.plot([X[i], X[i]], [Smin[i], Smax[i]], 
            ls = '-', color = 'deepskyblue', alpha = 0.3,
            label = 'Stress range', zorder = 1)
    ax.plot([X[i]], [Smin[i]], 
            marker = 'o', markerfacecolor = 'skyblue', markeredgecolor = 'k', markersize = 4,
            ls = '')
    ax.plot([X[i]], [Smax[i]], 
            marker = 'o', markerfacecolor = 'royalblue', markeredgecolor = 'k', markersize = 4,
            ls = '')

ax.set_xlabel('Cortical thickness (nm)')
ax.set_xlim([0, 1200])
locator = matplotlib.ticker.MultipleLocator(100)
ax.xaxis.set_major_locator(locator)

ax.set_ylabel('Extremal stress values (Pa)')
ax.set_ylim([0, 2500])
locator = matplotlib.ticker.MultipleLocator(100)
ax.yaxis.set_major_locator(locator)
ax.yaxis.grid(True)

ax.set_title('Feb 22 - Extremal stress values vs. thickness, for each compression')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_CompressionsLowStart_StressRanges', figSubDir = 'NonLin')
print(Ncomp)
plt.show()


# %%%% >>> OPTION 3 - OLIVIA'S IDEA

# %%%%% Make the dataframe


# Make the dataframe

data = GlobalTable_meca_nonLin

dates = ['22-02-09'] #['21-12-08', '22-01-12'] ['21-01-18', '21-01-21', '21-12-08']

filterList = [(data['validatedThickness'] == True),
              (data['cell subtype'] == 'aSFL'), 
              (data['bead type'] == 'M450'),
              (data['substrate'] == '20um fibronectin discs'),
              (data['date'].apply(lambda x : x in dates))]  # (data['validatedFit'] == True), 
globalFilter = filterList[0]
for k in range(1, len(filterList)):
    globalFilter = globalFilter & filterList[k]

data_f = data[globalFilter]

fitMin = [S for S in range(25,1025,50)]
fitMax = [S+150 for S in fitMin]
fitCenters = np.array([S+75 for S in fitMin])
regionFitsNames = [str(fitMin[ii]) + '<s<' + str(fitMax[ii]) for ii in range(len(fitMin))]

listColumnsMeca = []

KChadwick_Cols = []
KWeight_Cols = []

for rFN in regionFitsNames:
    listColumnsMeca += ['KChadwick_'+rFN, 'K_CIW_'+rFN, 'R2Chadwick_'+rFN, 'K2Chadwick_'+rFN, 
                        'H0Chadwick_'+rFN, 'Npts_'+rFN, 'validatedFit_'+rFN]
    KChadwick_Cols += [('KChadwick_'+rFN)]

    K_CIWidth = data_f['K_CIW_'+rFN] #.apply(lambda x : x.strip('][').split(', ')).apply(lambda x : (np.abs(float(x[0]) - float(x[1]))))
    KWeight = (data_f['KChadwick_'+rFN]/K_CIWidth)**2
    data_f['K_Weight_'+rFN] = KWeight
    data_f['K_Weight_'+rFN] *= data_f['KChadwick_'+rFN].apply(lambda x : (x<1e6))
    data_f['K_Weight_'+rFN] *= data_f['R2Chadwick_'+rFN].apply(lambda x : (x>1e-2))
    data_f['K_Weight_'+rFN] *= data_f['K_CIW_'+rFN].apply(lambda x : (x!=0))
    KWeight_Cols += [('K_Weight_'+rFN)]
    
data_f.tail()


# %%%%% K(s) cell by cell


CondCol = 'date'
Variables = KChadwick_Cols
WeightCols = KWeight_Cols
data_f_agg = getAggDf_weightedAvg(data_f, 'cellID', CondCol, Variables, WeightCols)
data_f_agg.T
dictPlot = {'cellID' : [], 'Kavg' : [], 'Kstd' : [], 'Kcount' : []}
meanCols = []
stdCols = []
countCols = []
for col in data_f_agg.columns:
    if 'KChadwick' in col and 'Wmean' in col:
        meanCols.append(col)
    elif 'KChadwick' in col and '_Wstd' in col:
        stdCols.append(col)
    elif 'KChadwick' in col and '_Count' in col:
        countCols.append(col)
meanDf = data_f_agg[meanCols]
stdDf = data_f_agg[stdCols]
countDf = data_f_agg[countCols]
for c in data_f_agg.index:
    means = meanDf.T[c].values
    stds = stdDf.T[c].values
    counts = countDf.T[c].values
    dictPlot['cellID'].append(c)
    dictPlot['Kavg'].append(means)
    dictPlot['Kstd'].append(stds)
    dictPlot['Kcount'].append(counts)


fig, ax = plt.subplots(1,1, figsize = (9,6))
for i in range(len(dictPlot['cellID'])):
    c = dictPlot['cellID'][i]
    color = colorList10[i%10]
#     ax.errorbar(fitCenters, dictPlot['Kavg'][i], yerr = dictPlot['Kstd'][i], color = color)
    ax.plot(fitCenters, dictPlot['Kavg'][i], color = color)
    low =  dictPlot['Kavg'][i] - (dictPlot['Kstd'][i] / (dictPlot['Kcount'][i]**0.5)) 
    high = dictPlot['Kavg'][i] + (dictPlot['Kstd'][i] / (dictPlot['Kcount'][i]**0.5))
    low = np.where(low < 10, 10, low)
    matplotlib.pyplot.fill_between(x=fitCenters, 
                                   y1=low, 
                                   y2=high,
                                   color = color, alpha = 0.1, zorder = 1)

ax.set_xlabel('Stress (Pa)')
ax.set_ylabel('K (Pa)')
ax.set_yscale('log')
ax.set_ylim([1e2, 1e5])
fig.suptitle('Stress stiffening - On each cell') #\n1 day, 3 expts, 36 cells, 232 compression
# jvu.archiveFig(fig, ax, todayFigDirLocal + '//NonLin', name='NonLin_K(s)_cellByCell', dpi = 100)
plt.show()


# %%%%% K(s) two different error types


# print(data_f.head())

# print(fitCenters)

def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in range(100,1100,50):
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_f[variable].apply(nan2zero).values
    w = data_f[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

fig, ax = plt.subplots(2,1, figsize = (9,12))

ax[0].errorbar(fitCenters, Kavg, yerr = q*Kste, marker = 'o', color = my_default_color_list[0], 
               ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error')
ax[0].set_ylim([50,1e5])

ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
               ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
ax[1].set_ylim([500,1e6])

for k in range(2):
    ax[k].legend(loc = 'upper left')
    ax[k].set_yscale('log')
    ax[k].set_xlabel('Stress (Pa)')
    ax[k].set_ylabel('K (Pa)')
    for kk in range(len(N)):
        ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

fig.suptitle('Stress stiffening - All compressions pooled') # \n1 day, 3 expts, 36 cells, 232 compression
# jvu.archiveFig(fig, ax, todayFigDirLocal + '//NonLin', name='NonLin_K(s)_twoErrorsTypes', dpi = 100)
plt.show()

df_val = pd.DataFrame(d_val)
dftest = pd.DataFrame(d)

df_val


# %%%%% K(s) with only ste for error bars


# print(data_f.head())

# print(fitCenters)

def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in fitCenters:
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_f[variable].apply(nan2zero).values
    w = data_f[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

# FIT ?
X, Y = fitCenters, Kavg
eqnText = "Linear Fit"
params, results = jvu.fitLine(X, Y) # Y=a*X+b ; params[0] = b,  params[1] = a
pval = results.pvalues[1] # pvalue on the param 'a'
eqnText += " ; Y = {:.1f} X + {:.1f}".format(params[1], params[0])
eqnText += " ; p-val = {:.3f}".format(pval)
print("Y = {:.5} X + {:.5}".format(params[1], params[0]))
print("p-value on the 'a' coefficient: {:.4e}".format(pval))
fitY = params[1]*X + params[0]
imin = np.argmin(X)
imax = np.argmax(X)


fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# ax[0]
ax.errorbar(fitCenters, Kavg/1000, yerr = q*Kste/1000, 
            marker = 'o', color = my_default_color_list[0],
            markersize = 7.5, lw = 2,
            ecolor = 'k', elinewidth = 1.5, capsize = 5, 
            label = 'Weighted means\nWeighted ste 95% as error')
ax.set_ylim([0,1.6e4/1000])
ax.set_xlim([0,1150])

ax.plot([X[imin],X[imax]], [fitY[imin]/1000,fitY[imax]/1000], 
        ls = '--', lw = '1', color = 'darkorange', zorder = 1, label = eqnText)

# ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# ax[1].set_ylim([500,1e6])

# for k in range(1): #2
#     ax[k].legend(loc = 'upper left')
#     ax[k].set_yscale('log')
#     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
#     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
#     for kk in range(len(N)):
#         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
ax.legend(loc = 'upper left')
ax.set_yscale('linear')
ax.set_xlabel('Stress (Pa)') #' [center of a 150Pa large interval]')
ax.set_ylabel('Tangeantial Modulus (kPa)') #' [tangeant modulus w/ Chadwick]')
for kk in range(len(N)):
    # ax.text(x=fitCenters[kk]+5, y=(Kavg[kk]**0.98), s='n='+str(N[kk]), fontsize = 10)
    ax.text(x=fitCenters[kk]+5, y=(Kavg[kk]-500)/1000, s='n='+str(N[kk]), fontsize = 8)

# fig.suptitle('K(sigma)')
ax.set_title('Stress-stiffening of the actin cortex\nALL compressions pooled') #  '\n22-02-09 experiment, 36 cells, 232 compression')
jvu.archiveFig(fig, ax, todayFigDirLocal + '//NonLin', name='NonLin_K(s)_lin', dpi = 100)
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


# %%%%% 100-800Pa range


# print(data_f.head())

# print(fitCenters)

extraFilters = [data_f['minStress'] <= 100, data_f['maxStress'] >= 800]
fitCenters2 = fitCenters[fitCenters<850]

data_ff = data_f
globalExtraFilter = extraFilters[0]
for k in range(1, len(extraFilters)):
    globalExtraFilter = globalExtraFilter & extraFilters[k]

data_ff = data_ff[globalExtraFilter]
data_ff
def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in fitCenters2:
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_ff[variable].apply(nan2zero).values
    w = data_ff[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# ax[0]
ax.errorbar(fitCenters2, Kavg, yerr = q*Kste, marker = 'o', color = my_default_color_list[0], 
               ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error')
ax.set_ylim([0,1.4e4])

# ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# ax[1].set_ylim([500,1e6])

# for k in range(1): #2
#     ax[k].legend(loc = 'upper left')
#     ax[k].set_yscale('log')
#     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
#     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
#     for kk in range(len(N)):
#         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
ax.legend(loc = 'upper left')
ax.set_yscale('linear')
ax.set_xlabel('Stress (Pa)')
ax.set_ylabel('K (Pa)')
for kk in range(len(N)):
    ax.text(x=fitCenters2[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

# fig.suptitle('K(sigma)')
ax.set_title('Stress stiffening\nOnly compressions including the [100, 800]Pa range')
jvu.archiveFig(fig, ax, todayFigDirLocal + '//NonLin', name='NonLin_K(s)_100-800Pa_lin', dpi = 100)
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


# %%%%% 100-700Pa range


# print(data_f.head())

# print(fitCenters)

extraFilters = [data_f['minStress'] <= 100, data_f['maxStress'] >= 700]
fitCenters2 = fitCenters[fitCenters<=700]

data_ff = data_f
globalExtraFilter = extraFilters[0]
for k in range(1, len(extraFilters)):
    globalExtraFilter = globalExtraFilter & extraFilters[k]

data_ff = data_ff[globalExtraFilter]
data_ff
def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in fitCenters2:
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_ff[variable].apply(nan2zero).values
    w = data_ff[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# ax[0]
ax.errorbar(fitCenters2, Kavg, yerr = q*Kste, marker = 'o', color = my_default_color_list[0], 
               ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nWeighted ste 95% as error')
ax.set_ylim([0,1.4e4])

# ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# ax[1].set_ylim([500,1e6])

# for k in range(1): #2
#     ax[k].legend(loc = 'upper left')
#     ax[k].set_yscale('log')
#     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
#     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
#     for kk in range(len(N)):
#         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
ax.legend(loc = 'upper left')
ax.set_yscale('linear')
ax.set_xlabel('Stress (Pa)')
ax.set_ylabel('K (Pa)')
for kk in range(len(N)):
    ax.text(x=fitCenters2[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

# fig.suptitle('K(sigma)')
ax.set_title('Stress stiffening\nOnly compressions including the [100, 700]Pa range')
jvu.archiveFig(fig, ax, todayFigDirLocal + '//NonLin', name='NonLin_K(s)_100-700Pa_lin', dpi = 100)
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


# %%%%% 100-600Pa range


# print(data_f.head())

# print(fitCenters)

extraFilters = [data_f['minStress'] <= 100, data_f['maxStress'] >= 550]
fitCenters2 = fitCenters[fitCenters<=550]

data_ff = data_f
globalExtraFilter = extraFilters[0]
for k in range(1, len(extraFilters)):
    globalExtraFilter = globalExtraFilter & extraFilters[k]

data_ff = data_ff[globalExtraFilter]
print(len(data_ff['cellID'].unique()), data_ff['cellID'].unique())
inspectDict = {}
for cId in data_ff['cellID'].unique():
    inspectDict[cId] = data_ff[data_ff['cellID'] == cId].shape[0]
print(inspectDict)

def w_std(x, w):
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    return(std)

def nan2zero(x):
    if np.isnan(x):
        return(0)
    else:
        return(x)

valStr = 'KChadwick_'
weightStr = 'K_Weight_'

Kavg = []
Kstd = []
D10 = []
D90 = []
N = []

for S in fitCenters2:
    rFN = str(S-75) + '<s<' + str(S+75)
    variable = valStr+rFN
    weight = weightStr+rFN
    
    x = data_ff[variable].apply(nan2zero).values
    w = data_ff[weight].apply(nan2zero).values
    
    if S == 250:
        d = {'x' : x, 'w' : w}
    
    m = np.average(x, weights=w)
    v = np.average((x-m)**2, weights=w)
    std = v**0.5
    
    d10, d90 = np.percentile(x[x != 0], (10, 90))
    n = len(x[x != 0])
    
    Kavg.append(m)
    Kstd.append(std)
    D10.append(d10)
    D90.append(d90)
    N.append(n)

Kavg = np.array(Kavg)
Kstd = np.array(Kstd)
D10 = np.array(D10)
D90 = np.array(D90)
N = np.array(N)
Kste = Kstd / (N**0.5)

alpha = 0.975
dof = N
q = st.t.ppf(alpha, dof) # Student coefficient

d_val = {'S' : fitCenters, 'Kavg' : Kavg, 'Kstd' : Kstd, 'D10' : D10, 'D90' : D90, 'N' : N}

fig, ax = plt.subplots(1,1, figsize = (9,6)) # (2,1, figsize = (9,12))

# ax[0]
ax.errorbar(fitCenters2, Kavg/1000, yerr = q*Kste/1000, 
            marker = 'o', color = my_default_color_list[0],
            markersize = 10, lw = 2,
            ecolor = 'k', elinewidth = 1.5, capsize = 5, 
            label = 'Weighted means\nWeighted ste 95% as error')
ax.set_ylim([0,1.0e4/1000])

# ax[1].errorbar(fitCenters, Kavg, yerr = [D10, D90], marker = 'o', color = my_default_color_list[3], 
#                ecolor = 'k', elinewidth = 0.8, capsize = 3, label = 'Weighted means\nD9-D1 as error')
# ax[1].set_ylim([500,1e6])

# for k in range(1): #2
#     ax[k].legend(loc = 'upper left')
#     ax[k].set_yscale('log')
#     ax[k].set_xlabel('Stress (Pa) [center of a 150Pa large interval]')
#     ax[k].set_ylabel('K (Pa) [tangeant modulus w/ Chadwick]')
#     for kk in range(len(N)):
#         ax[k].text(x=fitCenters[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)
ax.legend(loc = 'upper left')
ax.set_yscale('linear')
ax.set_xlabel('Stress (Pa)') #  [center of a 150Pa large interval]
ax.set_ylabel('Tangeantial Modulus (kPa)') # [tangeant modulus w/ Chadwick]
#### Decomment to get the number of compressions again
# for kk in range(len(N)):
#     ax.text(x=fitCenters2[kk]+5, y=Kavg[kk]**0.98, s='n='+str(N[kk]), fontsize = 6)

# fig.suptitle('K(sigma)')
ax.set_title('Stress-stiffening of the actin cortex')#'\nOnly compressions including the [100, 600]Pa range')

jvu.archiveFig(fig, ax, todayFigDirLocal + '//NonLin', name='NonLin_K(s)_100-600Pa_lin', dpi = 100)
plt.show()

# df_val = pd.DataFrame(d_val)
# dftest = pd.DataFrame(d)

# df_val


Kavg_ref, Kste_ref, N_ref, fitCenters_ref = Kavg, Kste, N, fitCenters2



# %%%%% Multiple boxplot K(s)

data = data_f

listeS = [100, 150, 200, 250, 300, 400, 500, 600, 700, 800]

fig, axes = plt.subplots(1, len(listeS), figsize = (14,6))
kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['date'].apply(lambda x : x in ['22-02-09']))] #, '21-12-16'

    fig, ax = D1Plot(data, fig=fig, ax=axes[kk], CondCol=['bead type'], Parameters=['KChadwick_'+interval], Filters=Filters, 
                   Boxplot=True, cellID='cellID', co_order=[], stats=True, statMethod='Mann-Whitney', AvgPerCell = True,
                   box_pairs=[], figSizeFactor = 1, markersizeFactor=1, orientation = 'h')

#     axes[kk].legend(loc = 'upper right', fontsize = 8)
    label = axes[kk].get_ylabel()
    axes[kk].set_title(label.split('_')[-1] + 'Pa', fontsize = 10)
    axes[kk].set_yscale('log')
    axes[kk].set_ylim([1e2, 5e4])
    axes[kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
    axes[kk].set_xticklabels([])
    if kk == 0:
        axes[kk].set_ylabel(label.split('_')[0] + ' (Pa)')
        axes[kk].tick_params(axis='x', labelsize = 10)
    else:
        axes[kk].set_ylabel(None)
        axes[kk].set_yticklabels([])
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan22_M450vsM270_CompressionsLowStart_Detailed1DPlot', figSubDir = figSubDir)
    
    kk+=1

renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
fig.tight_layout()
fig.suptitle('Tangeantial modulus of aSFL 3T3 - 0.5>50mT lowStartCompressions - avg per cell')
plt.show()
jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_nonLinK_multipleIntervals', figSubDir = 'NonLin')


# %%%%% Multiple K(h)



data = GlobalTable_meca_nonLin

listeS = [100, 150, 200, 250, 300, 400, 500, 600, 700, 800]

# fig, axes = plt.subplots(len(listeS), 1, figsize = (9,40))
# kk = 0

for S in listeS:
    interval = str(S-75) + '<s<' + str(S+75)
    textForFileName = str(S-75) + '-' + str(S+75) + 'Pa'

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['bestH0'] <= 600),
               (data['date'].apply(lambda x : x in ['22-02-09']))] #, '21-12-16'

    fig, ax = D2Plot_wFit(data, #fig=fig, ax=axes[kk], 
                          XCol='bestH0', YCol='KChadwick_'+interval, CondCol = ['bead type'], Filters=Filters, 
                          cellID = 'cellID', AvgPerCell=False, 
                          xscale = 'linear', yscale = 'log', modelFit=True, modelType='y=k*x^a')
    
    # individual figs
    ax.set_xlim([8e1, 700])
    ax.set_ylim([1e2, 5e4])
    renameAxes(ax,{'bestH0':'best H0 (nm)', 'doxycyclin':'iMC'})
    fig.set_size_inches(6,4)
    fig.tight_layout()
#     jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_nonLin_K(h)_' + textForFileName, figSubDir = 'NonLin')
    
    
    # one big fig
#     axes[kk].set_xlim([8e1, 1.2e3])
#     axes[kk].set_ylim([3e2, 5e4])
#     kk+=1


# renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
# fig.tight_layout()
plt.show()


# %%%%% Multiple K(h) - 2



data = GlobalTable_meca_nonLin

listeS = [200]

# fig, axes = plt.subplots(len(listeS), 1, figsize = (9,40))
# kk = 0

for S in listeS:
    halfWidth = 125
    interval = str(S-halfWidth) + '<s<' + str(S+halfWidth)
    textForFileName = str(S-halfWidth) + '-' + str(S+halfWidth) + 'Pa'

    Filters = [(data['validatedFit_'+interval] == True),
               (data['validatedThickness'] == True),
               (data['bestH0'] <= 650),
               (data['date'].apply(lambda x : x in ['22-02-09']))] #, '21-12-16'

    fig, ax = D2Plot_wFit(data, #fig=fig, ax=axes[kk], 
                          XCol='bestH0', YCol='KChadwick_'+interval, 
                          CondCol = ['bead type'], Filters=Filters, 
                          cellID = 'cellID', AvgPerCell=False, 
                          xscale = 'log', yscale = 'log', modelFit=True, 
                          modelType='y=k*x^a', markersizeFactor = 1.2)
    
    # individual figs
    ax.set_xlim([2e2, 700])
    ax.set_ylim([5e2, 3e4])
    # ax.set_ylim([0, 8000])
    fig.suptitle(interval)
    ax.get_legend().set_visible(False)
    rD = {'bestH0':'Thickness (nm)',
          'KChadwick_'+interval : 'Tangeantial Modulus (Pa)'}
    renameAxes(ax,rD)
    fig.set_size_inches(8,7)
    fig.tight_layout()
#     jvu.archiveFig(fig, ax, name='3T3aSFL_Feb22_nonLin_K(h)_' + textForFileName, figSubDir = 'NonLin')
    
    
    # one big fig
#     axes[kk].set_xlim([8e1, 1.2e3])
#     axes[kk].set_ylim([3e2, 5e4])
#     kk+=1


# renameAxes(axes,{'none':'ctrl', 'doxycyclin':'iMC'})
# fig.tight_layout()
plt.show()

# %%%%% test


# dftest['produit'] = dftest.x.values*dftest.w.values
# wm = np.average(dftest.x.values, weights = dftest.w.values)
# print(np.sum(dftest.x.values*dftest.w.values))
# print(np.sum(dftest['produit']))
# print(np.sum(dftest.w.values))
# print((np.sum(dftest['produit']))/np.sum(dftest.w.values))
# print(wm)
# print(max(dftest.produit.values))
# dftest

# %%%% Different range width

# %%%%% All stress ranges - K(stress) -- With bestH0 distributions and E(h)


data = GlobalTable_meca_nonLin
dates = ['22-02-09']

fitC =  np.array([S for S in range(100, 1000, 100)])
fitW = [100, 150, 200, 250, 300]
fitCenters = np.array([[int(S) for S in fitC] for w in fitW])
fitWidth = np.array([[int(w) for S in fitC] for w in fitW])
fitMin = np.array([[int(S-(w/2)) for S in fitC] for w in fitW])
fitMax = np.array([[int(S+(w/2)) for S in fitC] for w in fitW])
# fitCenters, fitWidth = fitCenters[fitMin>0], fitWidth[fitMin>0]
# fitMin, fitMax = fitMin[fitMin>0], fitMax[fitMin>0]

for ii in range(len(fitW)):
    width = fitW[ii]
    fitMin_ii = fitMin[ii]
    N = len(fitMin_ii[fitMin_ii > 0])
    
    fig, axes = plt.subplots(3, N, figsize = (20,10))
    kk = 0
    
    for S in fitC:
        minS, maxS = int(S-width//2), int(S+width//2)
        interval = str(minS) + '<s<' + str(maxS)
        
        if minS > 0:
            print(interval)
            try:
                Filters = [(data['validatedFit_'+interval] == True),
                           (data['validatedThickness'] == True),
                           (data['bestH0'] <= 1100),
                           (data['date'].apply(lambda x : x in dates))]
                
                
                # axes[0, kk].set_yscale('log')
                axes[0, kk].set_ylim([5e2, 5e4])
                # axes[0, kk].set_yscale('linear')
                # axes[0, kk].set_ylim([0, 3e4])
            
                D1Plot(data, fig=fig, ax=axes[0, kk], CondCol=['bead type'], Parameters=['KChadwick_'+interval], 
                     Filters=Filters, Boxplot=True, cellID='cellID', 
                     stats=True, statMethod='Mann-Whitney', AvgPerCell = True, box_pairs=[], 
                     figSizeFactor = 1, markersizeFactor=1, orientation = 'h', 
                     stressBoxPlot = True)# , bypassLog = True)
                
                D1Plot(data, fig=fig, ax=axes[1, kk], CondCol=['bead type'], Parameters=['bestH0'], 
                     Filters=Filters, Boxplot=True, cellID='cellID', 
                     stats=True, statMethod='Mann-Whitney', AvgPerCell = True, box_pairs=[], 
                     figSizeFactor = 1, markersizeFactor=1, orientation = 'h', stressBoxPlot = True)
                
                D2Plot_wFit(data, fig=fig, ax=axes[2, kk], Filters=Filters, 
                      XCol='bestH0', YCol='KChadwick_'+interval, CondCol = ['bead type'],
                      cellID = 'cellID', AvgPerCell=True, xscale = 'log', yscale = 'log', 
                      modelFit=True, modelType='y=k*x^a', markersizeFactor = 1)
                
            except:
                pass
            
            # 0.
            # axes[0, kk].legend(loc = 'upper right', fontsize = 8)
            label0 = axes[0, kk].get_ylabel()
            axes[0, kk].set_title(label0.split('_')[-1] + 'Pa', fontsize = 10)
            axes[0, kk].set_ylim([5e2, 5e4])
            # axes[0, kk].set_ylim([1, 3e4])
            # axes[0, kk].set_xticklabels([])
            
            # 1.
            axes[1, kk].set_ylim([0, 1500])
            axes[1, kk].tick_params(axis='x', labelrotation = 50, labelsize = 10)
            
            # 2.
            # label2 = axes[2, kk].get_ylabel()
            axes[2, kk].set_xlim([8e1, 1.2e3])
            axes[2, kk].set_ylim([3e2, 5e4])
            axes[2, kk].tick_params(axis='x', labelrotation = 0, labelsize = 10)
            axes[2, kk].xaxis.label.set_size(10)
            axes[2, kk].legend().set_visible(False)
            
            for ll in range(axes.shape[0]):
                axes[ll, kk].tick_params(axis='y', labelsize = 10)
            
            if kk == 0:
                axes[0, kk].set_ylabel(label0.split('_')[0] + ' (Pa)')
                axes[1, kk].set_ylabel('Best H0 (nm)')
                axes[2, kk].set_ylabel('K_Chadwick (Pa)')
            else:
                for ll in range(axes.shape[0]):
                    axes[ll, kk].set_ylabel(None)
                    axes[ll, kk].set_yticklabels([])
        
            kk+=1
            
        else:
            pass
    
    renameAxes(axes.flatten(),{'none':'ctrl', 'doxycyclin':'iMC', 'bestH0' : 'Best H0 (nm)'})
    # fig.tight_layout()
    # fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')
    
    plt.show()


    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = todayFigDirLocal + '//' + figSubDir, figSubDir = figSubDir)
    # jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_K_multipleIntervals_250pa_lin', figDir = ownCloudTodayFigDir, figSubDir = figSubDir)


# %%%%% All stress ranges - K(stress)

#### Creation of data_new, with the different width piled on top of each other

data = GlobalTable_meca_nonLin
dates = ['22-02-09']

fitC =  np.array([S for S in range(100, 1000, 100)])
fitW = [100, 150, 200, 250, 300]
fitCenters = np.array([[int(S) for S in fitC] for w in fitW])
fitWidth = np.array([[int(w) for S in fitC] for w in fitW])
fitMin = np.array([[int(S-(w/2)) for S in fitC] for w in fitW])
fitMax = np.array([[int(S+(w/2)) for S in fitC] for w in fitW])
# fitCenters, fitWidth = fitCenters[fitMin>0], fitWidth[fitMin>0]
# fitMin, fitMax = fitMin[fitMin>0], fitMax[fitMin>0]

cols = data.columns
mainCols = []
for c in cols:
    strC = str(c)
    if '<s<' not in strC:
        mainCols.append(c)
        
data_bloc = data[mainCols]
nRows = data_bloc.shape[0]

cols_fits = np.array([['validatedFit_S='+str(S), 'KChadwick_S='+str(S), 
                       'K_CIW_S='+str(S), 'R2Chadwick_S='+str(S), 
                       'Npts_S='+str(S)] for S in fitC]).flatten()
for c in cols_fits:
    if 'validated' in c:
        data_bloc[c] = np.zeros(nRows, dtype=bool)
    else:
        data_bloc[c] = np.ones(nRows) * np.nan
        
data_bloc['fitWidth'] = np.ones(nRows)
data_new = data_bloc[:]
i = 0
for width in fitW:
    data_bloc['fitWidth'] = np.ones(nRows, dtype = int) * int(width)
    for S in fitC:
        try:
            minS, maxS = int(S-width//2), int(S+width//2)
            interval = str(minS) + '<s<' + str(maxS)
            data_bloc['validatedFit_S='+str(S)] = data['validatedFit_' + interval]
            data_bloc['KChadwick_S='+str(S)] = data['KChadwick_' + interval]
            data_bloc['K_CIW_S='+str(S)] = data['K_CIW_' + interval]
            data_bloc['R2Chadwick_S='+str(S)] = data['R2Chadwick_' + interval]
            data_bloc['Npts_S='+str(S)] = data['Npts_' + interval]
        except:
            minS, maxS = int(S-width//2), int(S+width//2)
            interval = str(minS) + '<s<' + str(maxS)
            data_bloc['ValidatedFit_S='+str(S)] = np.zeros(nRows, dtype=bool)
            data_bloc['KChadwick_S='+str(S)] = np.ones(nRows) * np.nan
            data_bloc['K_CIW_S='+str(S)] = np.ones(nRows) * np.nan
            data_bloc['R2Chadwick_S='+str(S)] = np.ones(nRows) * np.nan
            data_bloc['Npts_S='+str(S)] = np.ones(nRows) * np.nan
    if i == 0:
        data_new = data_bloc[:]
    else:
        data_new = pd.concat([data_new, data_bloc])
        
    i += 1
    
data_new = data_new.reset_index(drop=True)

#### Plot

fitCPlot = np.array([S for S in range(200, 900, 100)])
N = len(fitCPlot)
fig, axes = plt.subplots(1, N, figsize = (25,8))
kk = 0

for S in fitCPlot:
    # minS, maxS = int(S-width//2), int(S+width//2)
    # interval = str(minS) + '<s<' + str(maxS)
    Filters = [(data_new['validatedFit_S='+str(S)] == True),
               (data_new['validatedThickness'] == True),
               (data_new['bestH0'] <= 1100),
               (data_new['date'].apply(lambda x : x in dates))]
    
    
    # axes[0, kk].set_yscale('log')
    axes[kk].set_ylim([1e2, 5e4])
    # axes[0, kk].set_yscale('linear')
    # axes[0, kk].set_ylim([0, 3e4])

    D1Plot(data_new, fig=fig, ax=axes[kk], CondCol=['fitWidth'], Parameters=['KChadwick_S='+str(S)], 
         Filters=Filters, Boxplot=True, cellID='cellID', 
         stats=False, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
         figSizeFactor = 1, markersizeFactor=0.5, orientation = 'h', 
         stressBoxPlot = True)# , bypassLog = True)
    
    # D1Plot(data_new, fig=fig, ax=axes[1, kk], CondCol=['fitWidth'], Parameters=['bestH0'], 
    #      Filters=Filters, Boxplot=True, cellID='cellID', 
    #      stats=False, statMethod='Mann-Whitney', AvgPerCell = False, box_pairs=[], 
    #      figSizeFactor = 1, markersizeFactor=0.5, orientation = 'h', stressBoxPlot = True)
        
    # 0.
    # axes[0, kk].legend(loc = 'upper right', fontsize = 8)
    label0 = axes[kk].get_ylabel()
    axes[kk].set_title(label0.split('_')[-1] + 'Pa', fontsize = 10)
    axes[kk].set_ylim([4e2, 5e4])
    axes[kk].tick_params(axis='y', labelsize = 10)
    
    axes[kk].set_xlabel('Width (Pa)')
    axes[kk].tick_params(axis='x', labelrotation = 0, labelsize = 10)
    
    
    if kk == 0:
        axes[kk].set_ylabel(label0.split('_')[0] + ' (Pa)')

    else:
        axes[kk].set_ylabel(None)
        axes[kk].set_yticklabels([])

    kk+=1
    
fig.suptitle('Effect of fitting window on tangeantial modulus, data from 22-02-09')
renameAxes(axes.flatten(),{'none':'ctrl', 'doxycyclin':'iMC', 'bestH0' : 'Best H0 (nm)'})
# fig.tight_layout()
# fig.suptitle('Tangeantial modulus of aSFL 3T3 - control vs iMC linker')

plt.show()


# %%%  Drug experiments

# %%%%% 3T3 aSFL - March 22 - Compression experiments


data = GlobalTable_meca_Py2

dates = ['22-03-30']
fit = '175<s<325'

Filters = [(data['validatedFit'] == True), 
           (data['validatedThickness'] == True), 
           (data['bestH0'] <= 800),
           (data['validatedFit_' + fit] == True), 
           (data['drug'].apply(lambda x : x in ['dmso', 'blebbistatin'])),
           (data['date'].apply(lambda x : x in dates))]

fig, ax = D1Plot(data, CondCol=['drug'], Parameters=['bestH0', 'KChadwick_' + fit], Filters=Filters, 
                Boxplot=True, cellID='cellID', co_order=['dmso', 'blebbistatin'], stats=True, statMethod='Mann-Whitney', 
                box_pairs=[], figSizeFactor = 0.9, markersizeFactor=1, orientation = 'v', AvgPerCell = False)

rD = {'dmso' : 'Control', 'blebbistatin' : 'Blebbistatin',
              'bestH0' : 'Thickness (nm)', 'KChadwick_' + fit : 'Tangeantial Modulus (Pa)'}

renameAxes(ax, rD)
ax[0].set_ylim([0, 1000])
# ax[0].legend(loc = 'upper right', fontsize = 8)
# ax[1].legend(loc = 'upper right', fontsize = 8)
fig.suptitle('3T3aSFL & drugs\nPreliminary data')
# jvu.archiveFig(fig, ax, name='3T3aSFL_Jan21_drug_H&Echad_simple', figSubDir = figSubDir)
plt.show()


# %% Attempt at Interactive plots

# %%%%%


Filters = [(GlobalTable_ctField['validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL')]

p = D1PlotInteractive(GlobalTable_ctField, CondCol='drug',Parameters=['medianThickness','fluctuAmpli'],Filters=Filters)
p.children[0][0].title.text = '3T3aSFL on patterns: Ct Field'
p.children[0][0].title.text_font_size = '16pt'
p.children[1][0].title.text = ''
show(p)


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == 1), (GlobalTable_meca['cell subtype'] == 'aSFL')]
p = D1PlotInteractive(GlobalTable_meca, CondCol='drug',Parameters=['SurroundingThickness','EChadwick'],Filters=Filters,AvgPerCell=True,cellID='CellName')
p.children[0][0].title.text = '3T3aSFL on patterns: Compressions'
p.children[0][0].title.text_font_size = '14pt'
p.children[1][0].title.text = ''
show(p)


# %%%%%


Filters = [(GlobalTable_ctField['validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL')]
p = D2PlotInteractive(GlobalTable_ctField, XCol='medianThickness',YCol='fluctuAmpli',CondCol = 'drug', Filters=Filters)
show(p)


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL')]
p = D2PlotInteractive(GlobalTable_meca, XCol='SurroundingThickness',YCol='EChadwick',CondCol = 'drug', Filters=Filters, cellID = 'CellName')
show(p)


# %%%%%


Filters = [(GlobalTable_ctField['validated'] == True)]
p = D2PlotInteractive(GlobalTable_ctField, XCol='meanFluoPeakAmplitude',YCol='medianThickness',CondCol = 'drug', Filters=Filters, cellID = 'cellID',AvgPerCell=True)
show(p)


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True)]
p = D2PlotInteractive(GlobalTable_meca, XCol='meanFluoPeakAmplitude',YCol='SurroundingThickness',CondCol = 'drug', Filters=Filters, cellID = 'CellName',AvgPerCell=True)
p.title.text = '3T3aSFL expressing linker: H(fluo)'
p.title.text_font_size = '16pt'
show(p)


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL')]
p = D2PlotInteractive(GlobalTable_meca, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = 'drug', Filters=Filters, cellID = 'CellName',AvgPerCell=True)
p.title.text = 'aSFL expressing linker: E(fluo)'
p.title.text_font_size = '16pt'
show(p)


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL')]
p = D2PlotInteractive(GlobalTable_meca, XCol='meanFluoPeakAmplitude',YCol='SurroundingThickness',CondCol = 'drug', Filters=Filters, cellID = 'CellName',AvgPerCell=True)
p.title.text = 'aSFL expressing linker: H(fluo)'
p.title.text_font_size = '16pt'
show(p)


# %%%%%


Outliers = []
Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True),
           (GlobalTable_meca_Py['drug'].apply(lambda x : x in ['none', 'doxycyclin'])),
           (GlobalTable_meca_Py['substrate'] == '20um fibronectin discs'), 
           (GlobalTable_meca_Py['cell subtype'].apply(lambda x : x in ['aSFL-A8-2'])),
           (GlobalTable_meca_Py['cellID'].apply(lambda x : x not in Outliers))]
p = D2PlotInteractive(GlobalTable_meca_Py, XCol='meanFluoPeakAmplitude', YCol='EChadwick', CondCol = 'drug', 
                      Filters=Filters, cellID = 'cellID', AvgPerCell=True)
p.title.text = 'aSFL-6FP expressing long linker: E(fluo)'
p.title.text_font_size = '16pt'
show(p)


# %%%%%


Filters = [(GlobalTable_meca['Validated'] == True), (GlobalTable_meca['cell subtype'] == 'aSFL-6FP')]
p = D2PlotInteractive(GlobalTable_meca, XCol='meanFluoPeakAmplitude', YCol='SurroundingThickness', CondCol = 'drug', Filters=Filters, cellID = 'CellName',AvgPerCell=True)
p.title.text = 'aSFL-6FP expressing long linker: H(fluo)'
p.title.text_font_size = '16pt'
show(p)

# %% Utility scripts

# %%% Create a test table to test different fits

# Script !


testDir = os.path.join(dataDir, 'TestDataSet')
testFileName = 'testFitCompression.txt'
testFilePath = os.path.join(testDir, testFileName)

list_mecaFiles = [f for f in os.listdir(timeSeriesDataDir)                   if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv")                   and ('R40' in f))] # Change to allow different formats in the future
expDf = jvu.getExperimentalConditions()
tableDictTest = {}

for f in list_mecaFiles:
    print(f)
    tS_DataFilePath = os.path.join(timeSeriesDataDir, f)
    tsDF = pd.read_csv(tS_DataFilePath, ';')
    
    split_f = f.split('_')
    tsDF.dx, tsDF.dy, tsDF.dz, tsDF.D2, tsDF.D3 = tsDF.dx*1000, tsDF.dy*1000, tsDF.dz*1000, tsDF.D2*1000, tsDF.D3*1000
    thisManipID = split_f[0] + '_' + split_f[1]
    expDf['manipID'] = expDf['date'] + '_' + expDf['manip']
    thisExpDf = expDf.loc[expDf['manipID'] == thisManipID]
    DIAMETER = thisExpDf.at[thisExpDf.index.values[0], 'bead diameter']
    thisCellID = split_f[0] + '_' + split_f[1] + '_' + split_f[2] + '_' + split_f[3]
    Ncomp = max(tsDF['idxCompression'])
    
    for i in range(1, Ncomp+1): #Ncomp+1):
        thisCompHiD = thisCellID + '__' + str(i) + '__h'
        thisCompFiD = thisCellID + '__' + str(i) + '__f'
        print(thisCompHiD)
        
        thisCompDf = tsDF.loc[tsDF['idxCompression'] == i, :]
        iStart = (jvu.findFirst(tsDF['idxCompression'], i))
        iStop = iStart + thisCompDf.shape[0]
        
        # Delimit the start of the increase of B (typically the moment when the field decrease from 5 to 3)
        # and the end of its decrease (typically when it goes back from 3 to 5)
        
        listB = thisCompDf.B.values
        offsetStart, offsetStop = 0, 0
        minB, maxB = min(listB), max(listB)
        thresholdB = (maxB-minB)/50
        k = 0
        
        while (listB[k] > minB+thresholdB) or (listB[-1-k] > minB+thresholdB):
            offsetStart += int(listB[k] > minB+thresholdB)
            k += 1
        
        jStart = offsetStart
        jMax = np.argmax(thisCompDf.B)
        
        hCompr = (thisCompDf.D3.values[jStart:jMax+1] - DIAMETER)
        fCompr = (thisCompDf.F.values[jStart:jMax+1])
        
        tableDictTest[thisCompHiD] = hCompr
        tableDictTest[thisCompFiD] = fCompr
        
saveFile = open(testFilePath, 'w')
for k in tableDictTest.keys():
    saveFile.write(k)
    for i in range(len(tableDictTest[k])):
        saveFile.write(';')
        saveFile.write(str(tableDictTest[k][i]))
    saveFile.write('\n')
saveFile.close()


# %%% Create a test table to try statistical tests


# Script !


# %%%%

testDir = os.path.join(dataDir, 'TestDataSet')
GlobalTable_meca = jva.getGlobalTable_meca()
table_ExpConditions = jvu.getExperimentalConditions()
table_fluo = jva.getFluoData()
GlobalTable_meca = pd.merge(GlobalTable_meca, table_ExpConditions, how="inner", on='manipID',
#     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
#     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
)
GlobalTable_meca = pd.merge(GlobalTable_meca, table_fluo, how="left", left_on='CellName', right_on='cellID'
#     left_on=None,right_on=None,left_index=False,right_index=False,sort=True,
#     suffixes=("_x", "_y"),copy=True,indicator=False,validate=None,
)
print('Merged table has ' + str(GlobalTable_meca.shape[0]) + ' lines and ' + str(GlobalTable_meca.shape[1]) + ' columns.')


# %%%%


# Table 1

testFileName = 'testStats01_allComp.csv'
testFilePath = os.path.join(testDir, testFileName)

Filters = [(GlobalTable_meca['Validated'] == 1), (GlobalTable_meca['cell subtype'] == 'aSFL')]
co_order = makeOrder([['BSA coated glass','20um fibronectin discs'],['none','doxycyclin']])
data = GlobalTable_meca
CondCol=['substrate','drug']
Parameters=['SurroundingThickness','EChadwick']
AvgPerCell=False
cellID='CellName'

data_filtered = data
for fltr in Filters:
    data_filtered = data_filtered.loc[fltr]

NCond = len(CondCol)    
if NCond == 1:
    CondCol = CondCol[0]
elif NCond > 1:
    newColName = ''
    for i in range(NCond):
        newColName += CondCol[i]
        newColName += ' & '
    newColName = newColName[:-3]
    data_filtered[newColName] = ''
    for i in range(NCond):
        data_filtered[newColName] += data_filtered[CondCol[i]].astype(str)
        data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x + ' & ')
    data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x[:-3])
    CondCol = newColName
    
if AvgPerCell:
    group = data_filtered.groupby(cellID)
    dictAggMean = getDictAggMean(data_filtered)
    data_filtered = group.agg(dictAggMean)

data_filtered.sort_values(CondCol, axis=0, ascending=True, inplace=True)

df_output = data_filtered[[cellID, 'CompNum', newColName] + Parameters]

# df_output.to_csv(testFilePath)

df_output


# %%%%


# Table 2

testFileName = 'testStats02_avgPerCell.csv'
testFilePath = os.path.join(testDir, testFileName)

Filters = [(GlobalTable_meca['Validated'] == 1), (GlobalTable_meca['cell subtype'] == 'aSFL')]
co_order = makeOrder([['BSA coated glass','20um fibronectin discs'],['none','doxycyclin']])
data = GlobalTable_meca
CondCol=['substrate','drug']
Parameters=['SurroundingThickness','EChadwick']
AvgPerCell=True
cellID='CellName'

data_filtered = data
for fltr in Filters:
    data_filtered = data_filtered.loc[fltr]

NCond = len(CondCol)    
if NCond == 1:
    CondCol = CondCol[0]
elif NCond > 1:
    newColName = ''
    for i in range(NCond):
        newColName += CondCol[i]
        newColName += ' & '
    newColName = newColName[:-3]
    data_filtered[newColName] = ''
    for i in range(NCond):
        data_filtered[newColName] += data_filtered[CondCol[i]].astype(str)
        data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x + ' & ')
    data_filtered[newColName] = data_filtered[newColName].apply(lambda x : x[:-3])
    CondCol = newColName
    
if AvgPerCell:
    group = data_filtered.groupby(cellID)
    dictAggMean = getDictAggMean(data_filtered)
    data_filtered = group.agg(dictAggMean)

data_filtered.sort_values(CondCol, axis=0, ascending=True, inplace=True)

df_output = data_filtered[[cellID, 'CompNum', newColName] + Parameters]

# df_output.to_csv(testFilePath)

df_output


# %%%%


# Table 3
# Fake data

testDir = os.path.join(dataDir, 'TestDataSet')

testFileName = 'testStats03_FakeData.csv'
testFilePath = os.path.join(testDir, testFileName)

Npop = 10
Npoints = 30
minAvg = 100
maxAvg = 600
step = (maxAvg - minAvg)/(Npop-1)
std = 250
dictFakeData = {}
np.random.seed(11)

for i in [1, 2, 3, 10]:
    dictFakeData['Distribution_' + str(i)] = np.random.normal(loc=minAvg + step*(i-1), scale=std, size=Npoints)

dfFakeData = pd.DataFrame(dictFakeData)

# dfFakeData.to_csv(testFilePath)

dfFakeData


# %%%%

# Table 4
# Fake data 2

testDir = os.path.join(dataDir, 'TestDataSet')

testFileName = 'testStats04_FakeDataLarge.csv'
testFilePath = os.path.join(testDir, testFileName)

Npop = 10
Npoints = 300
minAvg = 100
maxAvg = 600
step = (maxAvg - minAvg)/(Npop-1)
std = 250
dictFakeData = {}
np.random.seed(11)

for i in [1, 2, 3, 10]:
    dictFakeData['Distribution_' + str(i)] = np.random.normal(loc=minAvg + step*(i-1), scale=std, size=Npoints)

dfFakeData = pd.DataFrame(dictFakeData)

# dfFakeData.to_csv(testFilePath)

dfFakeData


# %%% Comparison of stat tests

# %%%%  With the fake data

# %%%%%


testDir = os.path.join(dataDir, 'TestDataSet')

# Table 3
# Fake data
testFileName = 'testStats03_FakeData.csv'
testFilePath = os.path.join(testDir, testFileName)

dfFakeData = pd.read_csv(testFilePath)
dfFakeData = dfFakeData.drop(columns=['Unnamed: 0'])

Ncol = len(dfFakeData.columns)

refCol = dfFakeData[dfFakeData.columns[0]]
boxPlotMatrix = []

for i in range(0,Ncol):
    boxPlotMatrix.append(dfFakeData[dfFakeData.columns[i]].values)
    if i > 0:
        print('Comparison between distribution 1 and ' + str(i+1))
        tTest = st.ttest_ind(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
        print('tTest : ' + str(tTest.pvalue))
        wilcox = st.wilcoxon(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
        print('wilcox : ' + str(wilcox.pvalue))
        mannwhitneyu = st.mannwhitneyu(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
        print('mannwhitneyu : ' + str(mannwhitneyu.pvalue))
        print('')

fig, ax = plt.subplots(1,1)
ax.boxplot(boxPlotMatrix,labels=dfFakeData.columns.values)
ax.tick_params(axis='x', labelrotation = 15, labelsize = 7)
plt.show()


# %%%%%


testDir = os.path.join(dataDir, 'TestDataSet')

# Table 4
# Fake data Large
testFileName = 'testStats04_FakeDataLarge.csv'
testFilePath = os.path.join(testDir, testFileName)

dfFakeData = pd.read_csv(testFilePath)
dfFakeData = dfFakeData.drop(columns=['Unnamed: 0'])

Ncol = len(dfFakeData.columns)

refCol = dfFakeData[dfFakeData.columns[0]]
boxPlotMatrix = []

for i in range(0,Ncol):
    boxPlotMatrix.append(dfFakeData[dfFakeData.columns[i]].values)
    if i > 0:
        print('Comparison between distribution 1 and ' + str(i+1))
        tTest = st.ttest_ind(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
        print('tTest : ' + str(tTest.pvalue))
        wilcox = st.wilcoxon(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
        print('wilcox : ' + str(wilcox.pvalue))
        mannwhitneyu = st.mannwhitneyu(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
        print('mannwhitneyu : ' + str(mannwhitneyu.pvalue))
        print('')

fig, ax = plt.subplots(1,1)
ax.boxplot(boxPlotMatrix,labels=dfFakeData.columns.values)
ax.tick_params(axis='x', labelrotation = 15, labelsize = 7)
plt.show()


# %%%%%


# Summary of results: numpy.random.seed = 11
    
dictResultsFakeData = {}
dictResultsFakeData['language'] = ['python'   , 'python', 'python'      , 'R'     , 'R'          , 'Matlab' , 'Matlab' ]
dictResultsFakeData['test']     = ['ttest_ind', 'wilcox', 'mannwhitneyu', 't.test', 'wilcox.test', 'ttest2' , 'ranksum']
dictResultsFakeData['1 vs 2']   = [0.2700     , 0.2711  , 0.1353        , 0.2701  , 0.2729       , 0.2700   , 0.2707   ]
dictResultsFakeData['1 vs 3']   = [0.0714     , 0.0822  , 0.0452        , 0.0715  , 0.0906       , 0.0714   , 0.0905   ]
dictResultsFakeData['1 vs 10']  = [1.33e-11   , 4.28e-06, 2.98e-09      , 1.44e-11, 5.799e-11    , 1.33e-11 , 5.96e-09 ]
dfResultFakeData = pd.DataFrame(dictResultsFakeData)
dfResultFakeData


# %%%%%

# Summary of results: numpy.random.seed = 11

dictResultsFakeDataLarge = {}
dictResultsFakeDataLarge['language'] = ['python'   , 'python', 'python'      , 'R'     , 'R'          , 'Matlab' , 'Matlab' ]
dictResultsFakeDataLarge['test']     = ['ttest_ind', 'wilcox', 'mannwhitneyu', 't.test', 'wilcox.test', 'ttest2' , 'ranksum']
dictResultsFakeDataLarge['1 vs 2']   = [0.0082     , 0.0049  , 0.0038        , 0.0082  , 0.0077       , 0.0082   , 0.0077   ]
dictResultsFakeDataLarge['1 vs 3']   = [1.26e-06   , 9.29e-06, 1.37e-06      , 1.27e-06, 2.75e-06     , 1.26e-06 , 2.74e-06 ]
dictResultsFakeDataLarge['1 vs 10']  = [1.74e-98   , 5.64e-48, 6.39e-74      , 2.2e-16 , 2.2e-16      , 1.74e-98 , 1.27e-73 ]
dictResultsFakeDataLarge = pd.DataFrame(dictResultsFakeDataLarge)
dictResultsFakeDataLarge


# %%%%  With real data

# %%%%%


# testDir = os.path.join(dataDir, 'TestDataSet')

# Table 1
# Avg Per Cell
testFileName = 'testStats01_AvgPerCell.csv'
testFilePath = os.path.join(testDir, testFileName)

dfAvgPerCell = pd.read_csv(testFilePath)
# dfAvgPerCell = dfAvgPerCell.drop(columns=['Unnamed: 0'])
dfAvgPerCell
categories = list(dfAvgPerCell['substrate & drug'].unique())
Ncat = len(categories)
for i in range(Ncat):
    for j in range(i,Ncat):
        if j != i:
            x = dfAvgPerCell.loc[dfAvgPerCell['substrate & drug'] == categories[i], 'EChadwick'].values
            y = dfAvgPerCell.loc[dfAvgPerCell['substrate & drug'] == categories[j], 'EChadwick'].values
            print('Comparison between ' + categories[i] + ' and ' + categories[j])
            tTest = st.ttest_ind(x, y)
            print('tTest : ' + str(tTest.pvalue))
            mannwhitneyu = st.mannwhitneyu(x, y)
            print('mannwhitneyu : ' + str(mannwhitneyu.pvalue))
            print('')

# refCol = dfFakeData[dfFakeData.columns[0]]
# boxPlotMatrix = []

# for i in range(0,Ncol):
#     boxPlotMatrix.append(dfFakeData[dfFakeData.columns[i]].values)
#     if i > 0:
#         print('Comparison between distribution 1 and ' + str(i+1))
#         tTest = st.ttest_ind(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
#         print('tTest : ' + str(tTest.pvalue))
#         wilcox = st.wilcoxon(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
#         print('wilcox : ' + str(wilcox.pvalue))
#         mannwhitneyu = st.mannwhitneyu(refCol.values, dfFakeData[dfFakeData.columns[i]].values)
#         print('mannwhitneyu : ' + str(mannwhitneyu.pvalue))
#         print('')

# fig, ax = plt.subplots(1,1)
# ax.boxplot(boxPlotMatrix,labels=dfFakeData.columns.values)
# ax.tick_params(axis='x', labelrotation = 15, labelsize = 7)
# plt.show()


# %%%%%


dictResultsAvgPerCell = {}
order = ['language','test','c & NA vs iMC & A','c & NA vs c & A','iMC & NA vs iMC & A','c & NA vs iMC & NA','iMC & NA vs c & A','c & A vs iMC & A']
dictResultsAvgPerCell['language']            = ['python_statAnot'   , 'python_statAnot' ,'python'   , 'python', 'R'      , 'R'          , 'Matlab' , 'Matlab' ]
dictResultsAvgPerCell['test']                = ['ttest_ind', 'mannwhitneyu', 'ttest_ind', 'mannwhitneyu', 't.test' , 'wilcox.test', 'ttest2' , 'ranksum']
dictResultsAvgPerCell['c & NA vs iMC & NA']  = [1.000e+00  , 5.669e-01     ,0.19226082553146542  , 0.047244130659518     , 0.1503   , 0.09502       , 0.1923       , 0.0945]
dictResultsAvgPerCell['iMC & NA vs c & A']   = [9.673e-02  , 8.625e-03     ,0.016121864893694285  ,0.0007187494204219925     , 0.04082  , 0.0009726    , 0.0161     , 0.0014]
dictResultsAvgPerCell['c & A vs iMC & A']    = [2.331e-02  , 1.376e-01     ,0.00388458593467288  , 0.01146586677766893     , 0.007326 , 0.02214      , 0.0039      , 0.0229]
dictResultsAvgPerCell['c & NA vs c & A']     = [1.000e+00  , 1.000e+00     ,0.6977550928576132  , 0.2884535746840493     , 0.6948   , 0.5838       , 0.6978     , 0.5769]
dictResultsAvgPerCell['iMC & NA vs iMC & A'] = [1.000e+00  , 1.000e+00     ,0.5573451346686198  , 0.41831870120029446,  0.5031  , 0.8387       , 0.5573      , 0.8366 ]
dictResultsAvgPerCell['c & NA vs iMC & A']   = [9.726e-01  , 1.000e+00     ,0.16209530366557973  , 0.14893352365754048     , 0.04353  , 0.3043       , 0.1621       , 0.2979 ]
dfResultsAvgPerCell = pd.DataFrame(dictResultsAvgPerCell)
dfResultsAvgPerCell[order]

# dfResultsAvgPerCell = dfResultsAvgPerCell.sort_values(by='test',ascending=False)
# dfResultsAvgPerCell


# %% Plot small stuff

# %%%


SMALLER_SIZE = 6
SMALL_SIZE = 6
MEDIUM_SIZE = 10
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

F0, F1 = 500, 250
h0, h1 = 200, 50
om = 2*np.pi
delta = 1
T = np.arange(0,3,0.01)
F = F1 * np.sin(2*np.pi*T) + F0
h = h1 * np.sin(2*np.pi*T - delta) + h0

fig, ax = plt.subplots(figsize = (4,2))
color = 'firebrick'
ax.plot(T, F, color = color, lw = 1)
ax.plot([np.min(T), np.max(T)], [F0, F0], 'k--', lw = 1)
ax.set_xlabel('t (s)')
ax.set_ylabel('F (pN)', color=color)
ax.tick_params(axis='y', labelcolor=color)

axbis = ax.twinx()
color = 'blue'
axbis.set_ylabel('h (nm)', color=color)
axbis.plot(T, h, color=color, lw = 1)
axbis.tick_params(axis='y', labelcolor=color)
# axbis.set_yticks([0,500,1000,1500])
# minh = np.min(tsDF['D3'].values-DIAMETER)
# ratio = min(1/abs(minh/axM), 5)
# #             print(ratio)
#             (axmbis, axMbis) = ax1bis.get_ylim()
#             ax1bis.set_ylim([0, max(axMbis*ratio, 3*max(tsDF['F'].values))])


# %%% Experiment counter

# %%%%


# Experiment counter - Matlab table

# count ct field cells
cellID = 'cellID'
GlobalTable_ctField_CountCell = GlobalTable_ctField.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_ctField_CountCell = GlobalTable_ctField_CountCell.loc[:, [cellID]].rename(columns={cellID : 'Count cells - ctField'})

# count meca compressions
cellID = 'CellName'
GlobalTable_meca_CountComp = GlobalTable_meca.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountComp = GlobalTable_meca_CountComp.loc[:, [cellID]].rename(columns={cellID : 'Count compressions'})

# count valid meca compressions
cellID = 'CellName'
GlobalTable_meca_CountCompOK = GlobalTable_meca[GlobalTable_meca['Validated'] == True].groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountCompOK = GlobalTable_meca_CountCompOK.loc[:, [cellID]].rename(columns={cellID : 'Count OK compressions'})

# count meca cells
cellID = 'CellName'
group = GlobalTable_meca.groupby(cellID)
dictAggMean = getDictAggMean(GlobalTable_meca)
GlobalTable_meca_perCell = group.agg(dictAggMean)
GlobalTable_meca_CountCell = GlobalTable_meca_perCell.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountCell = GlobalTable_meca_CountCell.loc[:, [cellID]].rename(columns={cellID : 'Count cells - meca'})


# Fuse all the previous tables
GlobalTable_CountAll = pd.concat([GlobalTable_ctField_CountCell, GlobalTable_meca_CountCell, GlobalTable_meca_CountComp, GlobalTable_meca_CountCompOK], axis=1)
GlobalTable_CountAll = GlobalTable_CountAll.fillna(0)
GlobalTable_CountAll = GlobalTable_CountAll.loc[:,:].astype(int)
GlobalTable_CountAll


# %%%%


# Experiment counter - Python table

# count ct field cells
cellID = 'cellID'
GlobalTable_ctField_CountCell = GlobalTable_ctField.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_ctField_CountCell = GlobalTable_ctField_CountCell.loc[:, [cellID]].rename(columns={cellID : 'Count cells - ctField'})

# count meca compressions
cellID = 'cellID'
GlobalTable_meca_CountComp = GlobalTable_meca_Py.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountComp = GlobalTable_meca_CountComp.loc[:, [cellID]].rename(columns={cellID : 'Count compressions'})

# count valid meca compressions
cellID = 'cellID'
GlobalTable_meca_CountCompOK = GlobalTable_meca_Py[GlobalTable_meca_Py['validatedFit'] == True].groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountCompOK = GlobalTable_meca_CountCompOK.loc[:, [cellID]].rename(columns={cellID : 'Count OK compressions'})

# count meca cells
cellID = 'cellID'
group = GlobalTable_meca_Py.groupby(cellID)
dictAggMean = getDictAggMean(GlobalTable_meca_Py)
GlobalTable_meca_perCell = group.agg(dictAggMean)
GlobalTable_meca_CountCell = GlobalTable_meca_perCell.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountCell = GlobalTable_meca_CountCell.loc[:, [cellID]].rename(columns={cellID : 'Count cells - meca'})


# Fuse all the previous tables
GlobalTable_CountAll = pd.concat([GlobalTable_ctField_CountCell, GlobalTable_meca_CountCell, GlobalTable_meca_CountComp, GlobalTable_meca_CountCompOK], axis=1)
GlobalTable_CountAll = GlobalTable_CountAll.fillna(0)
GlobalTable_CountAll = GlobalTable_CountAll.loc[:,:].astype(int)
GlobalTable_CountAll


# %%%% Experiment counter - Python table 2


# 


# count ct field cells
cellID = 'cellID'
GlobalTable_ctField_CountCell = GlobalTable_ctField.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_ctField_CountCell = GlobalTable_ctField_CountCell.loc[:, [cellID]].rename(columns={cellID : 'Count cells - ctField'})

# count meca compressions
cellID = 'cellID'
GlobalTable_meca_CountComp = GlobalTable_meca_Py2.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountComp = GlobalTable_meca_CountComp.loc[:, [cellID]].rename(columns={cellID : 'Count compressions'})

# count valid meca compressions
cellID = 'cellID'
GlobalTable_meca_CountCompOK = GlobalTable_meca_Py2[GlobalTable_meca_Py2['validatedFit'] == True].groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountCompOK = GlobalTable_meca_CountCompOK.loc[:, [cellID]].rename(columns={cellID : 'Count OK compressions'})

# count meca cells
cellID = 'cellID'
group = GlobalTable_meca_Py2.groupby(cellID)
dictAggMean = getDictAggMean(GlobalTable_meca_Py2)
GlobalTable_meca_perCell = group.agg(dictAggMean)
GlobalTable_meca_CountCell = GlobalTable_meca_perCell.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountCell = GlobalTable_meca_CountCell.loc[:, [cellID]].rename(columns={cellID : 'Count cells - meca'})


# Fuse all the previous tables
GlobalTable_CountAll = pd.concat([GlobalTable_ctField_CountCell, GlobalTable_meca_CountCell, GlobalTable_meca_CountComp, GlobalTable_meca_CountCompOK], axis=1)
GlobalTable_CountAll = GlobalTable_CountAll.fillna(0)
GlobalTable_CountAll = GlobalTable_CountAll.loc[:,:].astype(int)
GlobalTable_CountAll

# %%%% Experiment counter


# Experiment counter - Python table 2
cellID = 'cellID'
data = GlobalTable_meca_Py2

data_f = data[data['date'].apply(lambda x : x in ['21-12-08', '22-01-12'])]

# count meca compressions

GlobalTable_meca_CountComp = data_f.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountComp = GlobalTable_meca_CountComp.loc[:, [cellID]].rename(columns={cellID : 'Count compressions'})

# count valid meca compressions

GlobalTable_meca_CountCompOK = data_f[data_f['validatedThickness'] == True].groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountCompOK = GlobalTable_meca_CountCompOK.loc[:, [cellID]].rename(columns={cellID : 'Count OK compressions'})

# count meca cells

group = data_f.groupby(cellID)
dictAggMean = getDictAggMean(data_f)
GlobalTable_meca_perCell = group.agg(dictAggMean)
GlobalTable_meca_CountCell = GlobalTable_meca_perCell.groupby(['cell type', 'cell subtype', 'bead type', 'drug', 'substrate']).count()
GlobalTable_meca_CountCell = GlobalTable_meca_CountCell.loc[:, [cellID]].rename(columns={cellID : 'Count cells - meca'})


# Fuse all the previous tables
GlobalTable_CountAll = pd.concat([GlobalTable_meca_CountCell, GlobalTable_meca_CountComp, GlobalTable_meca_CountCompOK], axis=1)
GlobalTable_CountAll = GlobalTable_CountAll.fillna(0)
GlobalTable_CountAll = GlobalTable_CountAll.loc[:,:].astype(int)
GlobalTable_CountAll


# %%%%


# Filters = [(GlobalTable_meca_Py['validatedFit'] == True), (GlobalTable_meca_Py['validatedThickness'] == True), (GlobalTable_meca_Py['cell subtype'] == 'aSFL-A8')]
# GlobalTable_meca_PyF = GlobalTable_meca_Py
# for fltr in Filters:
#         GlobalTable_meca_PyF = GlobalTable_meca_PyF.loc[fltr]
# # GlobalTable_mecaBis
# cellID = 'cellID'
# group = GlobalTable_meca_PyF.groupby(cellID)
# dictAggMean = getDictAggMean(GlobalTable_meca_PyF)
# GlobalTable_meca_PyF_perCell = group.agg(dictAggMean)
# GlobalTable_meca_PyF_perCell

