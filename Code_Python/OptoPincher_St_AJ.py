# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 16:27:23 2022

@author: anumi
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
import tifffile as tf

st.title('OptoPincher Plotter')

path = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData'
folder = '21-12-20_M2_P1_C3_disc20um'

def load_data(path, folder, bead_dia = 4.503):
    data = pd.read_csv(path+'/'+folder+'_PY.csv', sep=';')
    xyz_dist = data['D3'] - bead_dia
    xy_dist = data['D2'] - bead_dia
    dz = data['dz']
    t = (data['T']*1000)/60
    return xyz_dist, xy_dist, dz, t

    
def plot(path, folder, Nan_thresh = 'nan'):
    xyz_dist, xy_dist, dz, t = load_data(path, folder)
    outliers = np.where(xyz_dist > Nan_thresh)[0]
    xyz_dist[outliers] = 'nan'
    xyz_dist[outliers] = 'nan' 
    dz[outliers] = 'nan'
    df1 = pd.DataFrame({'Time': t, 'Thickness (3D)': xyz_dist})
    df1 = df1.set_index('Time')
    st.line_chart(df1)
    
    df2 = pd.DataFrame({'Time': t, 'XY Distance': xy_dist})
    df2 = df2.set_index('Time')
    st.line_chart(df2)
    
    df3 = pd.DataFrame({'Time': t, 'Dz':dz})
    df3 = df3.set_index('Time')
    st.line_chart(df3)



# %% Tests

from PIL import Image

path = 'D:/Anumita/MagneticPincherData/Raw/21.12.20/'
folder = '21-12-20_M1_P1_C1_disc20um'

stack = Image.open(path+folder+'.tif')

for i in range(4):
    try:
        stack.seek(i)
        print(stack.getpixel( (0, 0)))
    except EOFError:
        # Not enough frames in img
        break