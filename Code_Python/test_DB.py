# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 09:43:38 2022

@author: BioMecaCell
"""


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
import glob

from copy import copy
from cycler import cycler
from datetime import date
from scipy.optimize import curve_fit
from matplotlib.gridspec import GridSpec
from pathlib import Path



#### Paths

COMPUTERNAME = os.environ['COMPUTERNAME']
  
if COMPUTERNAME == 'DATA2JHODR':
    mainDir = "C://Users//BioMecaCell//Desktop//ActinCortexAnalysis"
    rawDir = "D:/Duya/MagneticPincherData"
    ownCloudDir = ""

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
from os import listdir
from os.path import isfile, join
output_notebook()



# %% Different plots



# dataset creation
# def plotKChadwick():
#     # path and file cells I want to take are K_Chadwick, 
#     absolute_path = os.path.join(mainDir, "/Data_Analysis/Global_MecaData_DB_06_10.csv")
#     absolute_path = "C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/Global_MecaData_DB_06_10.csv"
#     tD=pd.read_csv(absolute_path, sep=';')
#     for column in tD.columns:
#         if "KChadwick_S" in column:
#             print(column)
#             mask_col = column.replace("KChadwick_S", "validatedFit_S")
#             mask = tD[mask_col]
#             df = tD[[column, "cellID"]][mask]
#             groups = df.groupby(by="cellID")
#             list_values = []
#             labels = []
#             for cell_id, data in groups:
#                 list_values.append(data[column])
#                 labels.append(cell_id[-2:])
#             if labels != []:
#                 plt.figure()
#                 plt.boxplot(list_values, labels=labels)
#                 plt.title(column)
#                 plt.tight_layout()
#                 save_file = f"D:/Duya/MagneticPincherData/Figures/22_06_10/{column.replace('/', '')}.png"
#                 plt.savefig(save_file, dpi=200)
#             else:
#                 print("Warning: Column", column, "has no values")
#     return


# def plotKChadwick():
#     # path and file cells I want to take are K_Chadwick, 
#     absolute_path = os.path.join(mainDir, "/Data_Analysis/Global_MecaData_DB_06_10.csv")
#     absolute_path = "C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/Global_MecaData_DB_06_10.csv"
#     all_files="D:/Duya/MagneticPincherData/Raw/22.06.10"
#     tD=pd.read_csv(absolute_path, sep=';')
#     cell_sizes= []
#     names=os.listdir(all_files)
#     for file in names:
#         if file.endswith(".tif"):
#             cellsize=jvu.findInfosInFileName(file, 'substrate')
#             cell_sizes.append(cellsize)
#     for column in tD.index:
#         tD = tD.replace([''],'')
        
    
#     for column in tD.columns:
#         if "KChadwick_S" in column:
#             tips = sns.load_dataset(column)
#             print(tips)
#             mask_col = column.replace("KChadwick_S", "validatedFit_S")
#             mask = tD[mask_col]
#             df = tD[[column, "cellID"]][mask]
#             groups = df.groupby(by="cellID")
#             list_values = []
#             labels = []
#             for cell_id, data in groups:
#                 list_values.append(data[column])
#                 labels.append(cell_id[-2:])
#             if labels != []:
#                 plt.figure()
#                 plt.boxplot(list_values, cell_sizes)
#                 ax = sns.swarmplot(x=cell_sizes, y=column, hue="cells", data=tips)
#                 plt.title(column)
#                 plt.tight_layout()
#                 save_file = f"D:/Duya/MagneticPincherData/Figures/22_06_10/{column.replace(('/', ''))}.png"
#                 plt.savefig(save_file, dpi=200)
#             else:
#                 print("Warning: Column", column, "has no values")
#     return


# def plotKChadwick():
#     # path and file cells I want to take are K_Chadwick, 
#     absolute_path = "C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/Global_MecaData_DB_06_10.csv"
#     tD=pd.read_csv(absolute_path, sep=';')
#     for column in tD.columns:
#         if "KChadwick_S" in column:
#             print(column)
#             mask_col = column.replace("KChadwick_S", "validatedFit_S")
#             mask = tD[mask_col]
#             df = tD[[column, "beadSize"]][mask]
#             groups = df.groupby(by="beadSize")
#             list_values = []
#             labels = []
#             for bead_size, data in groups:
#                 list_values.append(data[column])
#                 labels.append(bead_size)
#                 print(labels)
#             if labels != []:
#                 plt.figure()
#                 plt.boxplot(list_values, labels=labels)
#                 plt.title(column)
#                 plt.tight_layout()
#                 save_file = f"D:/Duya/MagneticPincherData/Figures/22-06-10_bysize/{column.replace('/', '')}.png"
#                 plt.savefig(save_file, dpi=200)
#             else:
#                 print("Warning: Column", column, "has no values")
#     return

# def plotKChadwickadhesionratio():
#     absolute_path = "C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/Global_MecaData_DB_06_10.csv"
#     tD=pd.read_csv(absolute_path, sep=';')
    
#     for column in tD.columns:
#         if "KChadwick_S" in column:
#             print(column)
#             mask_col = column.replace("KChadwick_S", "validatedFit_S")
#             mask = tD[mask_col]
#             df = tD[[column, "ratio_adhesion_cell"]][mask]
#             groups = df.groupby(by="ratio_adhesion_cell")
#             list_values = []
#             labels = []
#             for ratio, data in groups:
#                 list_values.append(data[column])
#                 labels.append(ratio)
#                 print(labels)
#             if labels != []:
#                 plt.figure()
#                 sns.boxplot(x=list_values, y=labels, data=tD)
#                 sns.swarmplot(x=list_values,y=labels, data=tD)
#                 plt.show()
#                 plt.title(column)
#                 plt.tight_layout()
#                 save_file = f"D:/Duya/MagneticPincherData/Figures/Adhesion_ratio/{column.replace('/', '')}.png"
#                 plt.savefig(save_file, dpi=200)
#             else:
#                 print("Warning: Column", column, "has no values")
#     return

def swarm_plot_KChadwick_by_diameter():
    # Load data
    data_path="C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/Global_MecaData_DB_06_10.csv"
    file_path="D:/Duya/MagneticPincherData/Raw/22.06.10"
    df = pd.read_csv(data_path, sep=";")
    # Define the cell ids to plot and remove measurement for other cell ids
    cell_ids_to_plot = df["cellID"].unique()
    mask = df["cellID"].isin(cell_ids_to_plot)
    df = df[mask]
    # Make column with adhesion diameter
    csv_files = [f.name for f in Path(file_path).glob("*")]
    cell_id_to_diam = {}
    for f in csv_files:
        if f.endswith(".tif"):
            cell_id = jvu.findInfosInFileName(f, "cellID")
            diam = jvu.findInfosInFileName(f, "substrate")[4:-2]
            cell_id_to_diam[cell_id] = diam
    df["diameter"] = df["cellID"].map(cell_id_to_diam).astype(int)
    diam_to_plot = df["diameter"].unique()
    diam_to_plot.sort()
    
    for column in df.columns:
        if "KChadwick" in column:
            col_vf = column.replace("KChadwick", "validatedFit")
            validated_fit = df[col_vf]
            df_validated = df[["diameter", column]][validated_fit]
            # Count the occurences of each diameter
            diam_counts = df_validated["diameter"].value_counts()
            xticks_labels = [str(d) for d in diam_to_plot]
            for diam, count in diam_counts.iteritems():
                idx = xticks_labels.index(str(diam))
                xticks_labels[idx] = f"{diam}\n\n{count}"
            # Add Nan values for missing cell IDS
            non_missing_diam = df_validated["diameter"].unique()
            for diam in diam_to_plot:
                if not (diam in non_missing_diam):
                    # Add NaN row
                    nan_row = pd.DataFrame({"diameter": [diam], column: [np.NaN]})
                    df_validated = pd.concat([df_validated, nan_row])
            plt.figure(dpi=200, figsize=(8, 6))
            sns.boxplot(x="diameter", y=column, data=df_validated, linewidth=0.5, boxprops=dict(alpha=.3),width=0.4)
            sns.swarmplot(x="diameter", y=column, data=df_validated)
            plt.xticks(range(len(diam_to_plot)), xticks_labels)
            plt.title(column)
            plt.tight_layout()
            save_file = f"D:/Duya/MagneticPincherData/Figures/Adhesion_size/{column.replace('/', '')}.png"
            plt.savefig(save_file, dpi=200)
            plt.close()
    return




def swarm_plot_KChadwick_by_ratio():
   
    # Load data
    data_path="C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/Global_MecaData_DB_06_10.csv"
    ratio_file_path="C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/22-06-10_Volume.csv"
    df = pd.read_csv(data_path, sep=";")
    df_ratio=pd.read_csv(ratio_file_path, sep=";")
    # Remove nan values of adhesion ratio in df_ratio
    mask = df_ratio["ratio_adhesion"].isna()
    df_ratio = df_ratio[~mask]
    # Round the values of ratio_adhesion for visibility in the plots
    df_ratio["ratio_adhesion"] = df_ratio["ratio_adhesion"].round(2)
    # Since there are several rows with the same cell ID, keep only the first one
    df_ratio = df_ratio.groupby("cellID").agg("first")
    # Merge df_ratio into df
    df_merged=df.merge(df_ratio, on='cellID')
    # Define the list of ratios to plot
    ratios_to_plot = df_merged["ratio_adhesion"].unique()
    ratios_to_plot.sort()

    # Find KChad
    for column in df_merged.columns:
        if "KChadwick" in column:
            col_vf = column.replace("KChadwick", "validatedFit")
            validated_fit = df_merged[col_vf]
            df_validated = df_merged[["ratio_adhesion", column]][validated_fit]
            # Count the occurences of each adhesion ratio
            ratio_counts = df_validated["ratio_adhesion"].value_counts()
            xticks_labels = [str(r) for r in ratios_to_plot]
            for ratio, count in ratio_counts.iteritems():
                idx = xticks_labels.index(str(ratio))
                xticks_labels[idx] = f"{ratio}\n\n{count}"
            # Add Nan values for missing adhesion ratios
            non_missing_ratio = df_validated["ratio_adhesion"].unique()
            for ratio in ratios_to_plot:
                if not (ratio in non_missing_ratio):
                    # Add NaN row
                    nan_row = pd.DataFrame({"ratio_adhesion": [ratio], column: [np.NaN]})
                    df_validated = pd.concat([df_validated, nan_row])
            plt.figure(dpi=200, figsize=(10, 8))
            sns.boxplot(x="ratio_adhesion", y=column, data=df_validated,linewidth=0.5, boxprops=dict(alpha=.3),width=0.4)
            sns.swarmplot(x="ratio_adhesion", y=column, data=df_validated)
            plt.xticks(range(len(ratios_to_plot)), xticks_labels)
            plt.title(column)
            plt.tight_layout()
            save_file = f"D:/Duya/MagneticPincherData/Figures/Adhesion_ratio/{column.replace('/', '')}.png"
            plt.savefig(save_file, dpi=200)
            plt.close()
    return

def plot_H0_by_diameter():
    data_path="C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/Global_MecaData_DB_06_10.csv"
    file_path="D:/Duya/MagneticPincherData/Raw/22.06.10"
    df = pd.read_csv(data_path, sep=";")
    # Define the cell ids to plot and remove measurement for other cell ids
    cell_ids_to_plot = df["cellID"].unique()
    mask = df["cellID"].isin(cell_ids_to_plot)
    df = df[mask]
    # Make column with adhesion diameter
    csv_files = [f.name for f in Path(file_path).glob("*")]
    cell_id_to_diam = {}
    for f in csv_files:
        if f.endswith(".tif"):
            cell_id = jvu.findInfosInFileName(f, "cellID")
            diam = jvu.findInfosInFileName(f, "substrate")[4:-2]
            cell_id_to_diam[cell_id] = diam
    df["diameter"] = df["cellID"].map(cell_id_to_diam).astype(int)
    diam_to_plot = df["diameter"].unique()
    diam_to_plot.sort()
    
    # Find_bestH0
    if "bestH0" in df.columns:
        df_validated = df[["diameter","bestH0"]]
            # Count the occurences of each diameter
        diam_counts = df_validated["diameter"].value_counts()
        xticks_labels = [str(d) for d in diam_to_plot]
        for diam, count in diam_counts.iteritems():
            idx = xticks_labels.index(str(diam))
            xticks_labels[idx] = f"{diam}\n\n{count}"
            # Add Nan values for missing cell IDS
        non_missing_diam = df_validated["diameter"].unique()
        for diam in diam_to_plot:
            if not (diam in non_missing_diam):
                    # Add NaN row
                nan_row = pd.DataFrame({"diameter": [diam], "bestH0": [np.NaN]})
                df_validated = pd.concat([df_validated, nan_row])
        plt.figure(dpi=200, figsize=(8, 6))
        sns.boxplot(x="diameter", y="bestH0", data=df_validated, linewidth=0.5, boxprops=dict(alpha=.3),width=0.4)
        sns.swarmplot(x="diameter", y="bestH0", data=df_validated)
        plt.xticks(range(len(diam_to_plot)), xticks_labels)
        plt.title("bestH0")
        plt.tight_layout()
        save_file = f"D:/Duya/MagneticPincherData/Figures/bestH0_size.png"
        plt.savefig(save_file, dpi=200)
        plt.close()
    return

def plot_bestH0_by_ratio():
   
    # Load data
    data_path="C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/Global_MecaData_DB_06_10.csv"
    ratio_file_path="C:/Users/BioMecaCell/Desktop/ActinCortexAnalysis/Data_Analysis/22-06-10_Volume.csv"
    df = pd.read_csv(data_path, sep=";")
    df_ratio=pd.read_csv(ratio_file_path, sep=";")
    # Remove nan values of adhesion ratio in df_ratio
    mask = df_ratio["ratio_adhesion"].isna()
    df_ratio = df_ratio[~mask]
    # Round the values of ratio_adhesion for visibility in the plots
    df_ratio["ratio_adhesion"] = df_ratio["ratio_adhesion"].round(2)
    # Since there are several rows with the same cell ID, keep only the first one
    df_ratio = df_ratio.groupby("cellID").agg("first")
    # Merge df_ratio into df
    df_merged=df.merge(df_ratio, on='cellID')
    # Define the list of ratios to plot
    ratios_to_plot = df_merged["ratio_adhesion"].unique()
    ratios_to_plot.sort()

    # Find bestH0
    if "bestH0" in df_merged.columns:
        df_validated = df_merged[["ratio_adhesion", "bestH0"]]
            # Count the occurences of each adhesion ratio
        ratio_counts = df_validated["ratio_adhesion"].value_counts()
        xticks_labels = [str(r) for r in ratios_to_plot]
        for ratio, count in ratio_counts.iteritems():
            idx = xticks_labels.index(str(ratio))
            xticks_labels[idx] = f"{ratio}\n\n{count}"
            # Add Nan values for missing adhesion ratios
        non_missing_ratio = df_validated["ratio_adhesion"].unique()
        for ratio in ratios_to_plot:
            if not (ratio in non_missing_ratio):
                    # Add NaN row
                nan_row = pd.DataFrame({"ratio_adhesion": [ratio], "bestH0": [np.NaN]})
                df_validated = pd.concat([df_validated, nan_row])
        plt.figure(dpi=200, figsize=(10, 8))
        sns.boxplot(x="ratio_adhesion", y="bestH0", data=df_validated,linewidth=0.5, boxprops=dict(alpha=.3),width=0.4)
        sns.swarmplot(x="ratio_adhesion", y="bestH0", data=df_validated)
        plt.xticks(range(len(ratios_to_plot)), xticks_labels)
        plt.title("bestH0")
        plt.tight_layout()
        save_file = f"D:/Duya/MagneticPincherData/Figures/bestH0_ratio.png"
        plt.savefig(save_file, dpi=200) 
    return

