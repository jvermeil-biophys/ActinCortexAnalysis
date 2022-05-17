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
elif COMPUTERNAME == '':
    mainDir = "C://Users//josep//Desktop//ActinCortexAnalysis"
    ownCloudDir = "C://Users//josep//ownCloud//ActinCortexAnalysis"


experimentalDataDir = os.path.join(mainDir, "Data_Experimental")
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

#### Local imports
sys.path.append(mainDir + "//Code_Python")
# import PincherAnalysis_JV as jva
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
# allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) \
#                           if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
# print(allTimeSeriesDataFiles)


# %%% Get a time series

# df = aja.getCellTimeSeriesData('22-03-31_M9_P2_C2')


# # %%% Plot a time series

# aja.plotCellTimeSeriesData('22-03-31_M9_P2_C2')


# %%% Plot multiple time series

# allTimeSeriesDataFiles = [f for f in os.listdir(timeSeriesDataDir) if (os.path.isfile(os.path.join(timeSeriesDataDir, f)) and f.endswith(".csv"))]
# for f in allTimeSeriesDataFiles:
#     if '22-04-12_M2' in f:
#         aja.plotCellTimeSeriesData(f[:-4])


# # %%% Close all

# plt.close('all')


# #############################################################################
# %% GlobalTables plotting functions



# %%% Experimental conditions

expDf = jvu.getExperimentalConditions(experimentalDataDir, save=True , sep = ';')

# =============================================================================
# %%% Constant Field

# %%%% Update the table

# aja.computeGlobalTable_ctField(task='updateExisting', fileName = '', save=False, source = 'Python')



# %%%% Refresh the whole table

# aja.computeGlobalTable_ctField(task = 'fromScratch', fileName = '', save = True, source = 'Python')


# # %%%% Display

# df = aja.getGlobalTable_ctField().head()



# =============================================================================
# %%% Mechanics

# %%%% Update the table

# aja.computeGlobalTable_meca(task = 'updateExisting', fileName = '', 
#                             save = True, PLOT = True, source = 'Python')


# %%%% Refresh the whole table

aja.computeGlobalTable_meca(task = 'fromScratch', fileName = 'Global_MecaData_AJ', 
                            save = True, PLOT = True, source = 'Python')

# %%%% Specific experiments

# Task = '22-03-31_M9_P2_C2' # For instance '22-03-30 & '22-03-31'
# aja.computeGlobalTable_meca(task = Task, fileName = 'Global_MecaData_AJ', 
#                             save = True, PLOT = True, source = 'Python') # task = 'updateExisting'


# %%%% Precise dates (to plot)

# date = '22-03-31' # For instance '22-03-30 & '22-03-31'
# aja.computeGlobalTable_meca(task = date, fileName = 'Global_MecaData_AJ', 
#                             save = True, PLOT = False, source = 'Python') # task = 'updateExisting'


# %%%% Display

# df = aja.getGlobalTable_meca('Global_MecaData_Py2').tail()


# =============================================================================
# %%% Fluorescence

# %%%% Display

df = aja.getFluoData().head()




# #############################################################################
# %% > Data import & export

#### Data import

#### Display

# df1 = jvu.getExperimentalConditions().head()
# df2 = aja.getGlobalTable_ctField().head()
# df3 = aja.getGlobalTable_meca().head()
# df4 = aja.getFluoData().head()


#### GlobalTable_ctField

GlobalTable_ctField = aja.getGlobalTable(kind = 'ctField')
GlobalTable_ctField.head()


#### GlobalTable_ctField_Py

GlobalTable_ctField_Py = aja.getGlobalTable(kind = 'ctField_py')
GlobalTable_ctField_Py.head()


#### GlobalTable_meca

GlobalTable_meca = aja.getGlobalTable(kind = 'meca_matlab')
GlobalTable_meca.tail()


#### GlobalTable_meca_Py

GlobalTable_meca = aja.getGlobalTable(kind = 'Global_MecaData_AJ')
GlobalTable_meca.head()


#### GlobalTable_meca_Py2

GlobalTable_meca_Py2 = aja.getGlobalTable(kind = 'meca_py2')
GlobalTable_meca_Py2.head()


#### Global_MecaData_NonLin_Py

# GlobalTable_meca_nonLin = aja.getGlobalTable(kind = 'meca_nonLin')
GlobalTable_meca_nonLin = aja.getGlobalTable(kind = 'Global_MecaData_NonLin2_Py')
GlobalTable_meca_nonLin.head()


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

def TmodulusVsCompression(GlobalTable_meca, selectedStressRange, activationType = 'all'):
    sns.set_style('darkgrid')
    stressRanges = ['100+/-100', '150+/-100', '200+/-100', '250+/-100', \
                   '300+/-100', '350+/-100', '400+/-100', '450+/-100', \
                   '500+/-100', '550+/-100', '600+/-100', '650+/-100', \
                   '700+/-100', '750+/-100', '800+/-100', '850+/-100', \
                   '950+/-100', '1000+/-100', '1050+/-100', '1100+/-100']
    
    if not activationType == 'all':
        GlobalTable_meca = GlobalTable_meca[GlobalTable_meca['activation type'] == activationType]
    
    
    allFiles = np.unique(GlobalTable_meca['cellID'])
    # print(allFiles)
    lenSubplots = len(allFiles)
    rows= int(np.floor(np.sqrt(lenSubplots)))
    cols= int(np.ceil(lenSubplots/rows))
    fontsize = 15
    
    if selectedStressRange == 'all':
       for selectedStressRange in stressRanges:
           fig, axes = plt.subplots(nrows = rows, ncols = cols, figsize = (15,15))
           _axes = []
           for ax_array in axes:
               for ax in ax_array:
                   _axes.append(ax)
           for cellID, ax in zip(allFiles, _axes):
               print(cellID)
               GlobalTable_meca_spec = GlobalTable_meca[(GlobalTable_meca['cellID'] == cellID)]
               KChadwick = GlobalTable_meca_spec['KChadwick_S='+selectedStressRange].values
               compNum = GlobalTable_meca_spec['compNum'].values
               firstActivation = GlobalTable_meca['first activation'].values[0]
               categories = (GlobalTable_meca_spec['validatedFit_S='+selectedStressRange] == True).values
       
               categories = np.asarray(categories*1)
               colormap = np.asarray(['r', 'b'])
               labels = np.asarray(['Not valid', 'Valid'])
               ax.scatter(compNum, KChadwick, c = colormap[categories], label = labels[categories])
               ax.set_title(cellID, fontsize = 10)
               
               
               ax.set_ylim(0,10000)
               ax.set_xlim(0, len(compNum)+1)
               ax.axvline(x = firstActivation, color = 'r')
               
               plt.setp(ax.get_xticklabels(), fontsize=fontsize)
               plt.setp(ax.get_yticklabels(), fontsize=fontsize)
               fig.suptitle('KChadwick_S='+selectedStressRange+' vs. CompNum_'+activationType)
               try:
                   os.mkdir(todayFigDir)
               except:
                   pass
               
               plt.savefig(todayFigDir+'/TModulusVsComp_'+(selectedStressRange[:-6])+'_'+activationType+'.png')
           plt.show()
           plt.clf()
    else:
        fig, axes = plt.subplots(nrows = rows, ncols = cols, figsize = (15,15))
        _axes = []
        for ax_array in axes:
            for ax in ax_array:
                _axes.append(ax)
        for cellID, ax in zip(allFiles, _axes):
            print(cellID)
            GlobalTable_meca_spec = GlobalTable_meca[(GlobalTable_meca['cellID'] == cellID)]
            KChadwick = GlobalTable_meca_spec['KChadwick_S='+selectedStressRange].values
            compNum = GlobalTable_meca_spec['compNum'].values
            firstActivation = GlobalTable_meca['first activation'].values[0]
            categories = (GlobalTable_meca_spec['validatedFit_S='+selectedStressRange] == True).values
    
            categories = np.asarray(categories*1)
            colormap = np.asarray(['r', 'b'])
            labels = np.asarray(['Not valid', 'Valid'])
            ax.scatter(compNum, KChadwick, c = colormap[categories], label = labels[categories])
            ax.set_title(cellID, fontsize = 10)
            
            
            ax.set_ylim(0,10000)
            ax.set_xlim(0, len(compNum)+1)
            ax.axvline(x = firstActivation, color = 'r')
            
            plt.setp(ax.get_xticklabels(), fontsize=fontsize)
            plt.setp(ax.get_yticklabels(), fontsize=fontsize)
            fig.suptitle('KChadwick_S='+selectedStressRange+' vs. CompNum_'+activationType)
            try:
                os.mkdir(todayFigDir)
            except:
                pass
            
            plt.savefig(todayFigDir+'/TModulusVsComp_'+(selectedStressRange[:-6])+'_'+activationType+'.png')
        plt.show()
    


#%%%% Plotting TModulus vs. Compression Number
selectedStressRange = '200+/-100'
activationType = 'at beads'
TmodulusVsCompression(GlobalTable_meca, selectedStressRange, activationType)


# %%%

# H0 vs Compression Number
expt = '20220331_100xoil_3t3optorhoa_4.5beads_15mT_Mechanics'
f = '22-03-31_M8_P2_C2_disc20um_L40'
date = '22.03.31'

file = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/'+f+'_PY.csv'
tsDF = pd.read_csv(file, sep=';')


indices = GlobalTable_meca[GlobalTable_meca['cellID'] == jvu.findInfosInFileName(f, 'cellID')].index.tolist() 

plt.figure(figsize=(20,10))
plt.rcParams.update({'font.size': 35})
plt.plot(GlobalTable_meca['compNum'][indices].values, GlobalTable_meca['bestH0'][indices].values)
plt.axvline(x = 5.5, color = 'r', marker='.')
plt.xlabel('Compression Number')
plt.ylabel('bestH0 (nm)')
plt.title('bestH0 vs Compression No. | '+f)
plt.savefig('D:/Anumita/MagneticPincherData/Figures/Historique/2022-04-20/MecaAnalysis_allCells/'+f+'_H0vT.png')
plt.show()


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
                    fitY = params[1]*X + params[0]
                    imin = np.argmin(X)
                    imax = np.argmax(X)
                    ax.plot([X[imin],X[imax]], [fitY[imin],fitY[imax]], '--', lw = '1',
                            color = color, zorder = 4)

                elif modelType == 'y=A*exp(kx)':
                    params, results = jvu.fitLine(X, np.log(Y)) # Y=a*X+b ; params[0] = b,  params[1] = a
                    pval = results.pvalues[1] # pvalue on the param 'k'
                    eqnText += " ; Y = {:.1f}*exp({:.1f}*X)".format(params[0], params[1])
                    eqnText += " ; p-val = {:.3f}".format(pval)
                    print("Y = {:.5}*exp({:.5}*X)".format(np.exp(params[0]), params[1]))
                    print("p-value on the 'k' coefficient: {:.4e}".format(pval))
                    fitY = np.exp(params[0])*np.exp(params[1]*X)
                    imin = np.argmin(X)
                    imax = np.argmax(X)
                    ax.plot([X[imin],X[imax]], [fitY[imin],fitY[imax]], '--', lw = '1',
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
                    fitY = k * X**a
                    imin = np.argmin(X)
                    imax = np.argmax(X)
                    ax.plot([X[imin],X[imax]], [fitY[imin],fitY[imax]], '--', lw = '1', 
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











# %% Useful scripts

# %%% Experiment counter


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



