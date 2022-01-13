# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 18:49:02 2022

@author: anumi
"""

import streamlit as st
import pandas as pd
import numpy as np


st.title('OptoPincher Plotter')

def load_data(file):
    data = pd.read_csv(file, sep = ';')
    return data