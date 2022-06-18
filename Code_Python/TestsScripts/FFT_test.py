# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 16:49:19 2022

@author: JosephVermeil
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq

S = pd.read_csv('C://Users//JosephVermeil//Desktop//SignalToFT.txt', header = None, names = ['D3'])
D3 =  S.D3.values
tf_D3 = fft(D3)

# Number of sample points
N = len(D3)
# sample spacing
T = 20 / N
x = np.linspace(0.0, N*T, N, endpoint=False)
xf = fftfreq(N, T)[:N//2]

plt.plot(xf, 2.0/N * np.abs(tf_D3[0:N//2]))
plt.grid()
plt.show()