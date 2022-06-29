# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 14:36:25 2022

@author: anumi
"""

#%%
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression


def gaussianWeights(x0, X, tau):
    X = X.flatten(order='C')
    return np.exp ( -((X - x0) ** 2) / tau ** 2)

Xreal = np.sort(1000 * np.random.rand(100, 1), axis=0)
Xideal = np.linspace(0,1000,10000)
x0 = 200
tau = 50
Yreal = gaussianWeights(x0, Xreal, tau)
Yideal = gaussianWeights(x0, Xideal, tau)

plt.plot(Xreal, Yreal)
plt.plot(Xideal, Yideal)

#%%
    
y = np.sort(1000 * np.random.rand(100, 1)**2, axis=0)

x = np.linspace(0, len(y), len(y))
x = x.reshape(-1, 1)
# y = x**2 + np.sin(x)
slopes = []
intercepts = []

HR = 100

target = np.arange(100,900,100)
stressCompr = y
strainCompr = x

plt.plot(x, y, color = 'black', linewidth = 4)

#%%

for i in target:
    weights = gaussianWeights(i, stressCompr, HR)
    
    regr = LinearRegression()
    regr.fit(strainCompr, stressCompr, weights)
    slope = regr.coef_
    intercept = regr.intercept_
    stressPredict = slope*x + intercept
    
    strainCompr = strainCompr.flatten()
    stressCompr_copy = np.copy(stressCompr).flatten()
    
    idx_local = np.where(np.logical_and(stressCompr_copy >= i-HR, \
                                        stressCompr_copy <= i+HR))
    
    
    stressPredict = stressPredict[idx_local]
    strainCompr = strainCompr[idx_local]
    plt.plot(strainCompr, stressPredict, label = i)

plt.legend()
plt.show()