# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 16:27:23 2022

@author: anumi
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from PIL import Image, ImageSequence

st.title('OptoPincher Plotter')

path = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData'
folder = '21-12-20_M2_P1_C3_disc20um'
bead_pos = 'In' #'Out'

if bead_pos == 'In':
    id_pos = 1
else:
    id_pos = 0

def load_data(path, folder, bead_dia = 4.503):
    data = pd.read_csv(path+'/'+folder+'_PY.csv', sep=';')
    xyz_dist = data['D3'] - bead_dia
    xy_dist = data['D2'] - bead_dia
    dz = data['dz']
    t = (data['T']*1000)/60
    return xyz_dist, xy_dist, dz, t

def find_stack(stack, index):
    path = 'D:/Anumita/MagneticPincherData/Raw/21.12.20/'
    folder = '21-12-20_M1_P1_C1_disc20um'
    stack = Image.open(path+folder+'.tif')
    
    return img

def load_img(path, expt, folder, sliderNo):
    img = Image.open(path+'/'+expt+'/'+folder+'.tif')
    
    
data2 = 'C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Analysis/TimeSeriesData/Trajectories_raw/'+folder+'_rawTraj'+id_pos+'_PY.csv'
    
def plot(path, folder, Nan_thresh = 'nan'):
    xyz_dist, xy_dist, dz, t = load_data(path, folder)
    outlier = np.where(xyz_dist > Nan_thresh)[0]
    xyz_dist[outlier] = np.nan
    xy_dist[outlier] = np.nan
    dz[outlier] = np.nan
        
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
import matplotlib.pyplot as plt
import tifffile as tf

path = 'D:/Anumita/MagneticPincherData/Raw/21.12.20/'
folder = '21-12-20_M1_P1_C1_disc20um'
stack = tf.imread(path+folder+'.tif')

for i in range(len(stack)):
    plt.imshow(stack[i], cmap='gray')
    plt.axis('off')
    plt.show()


