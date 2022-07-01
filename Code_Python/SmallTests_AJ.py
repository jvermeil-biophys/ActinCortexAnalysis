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

#%%
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



#%%

fCompr = np.sort(1000 * np.random.rand(100, 1)**2, axis=0).flatten()

hCompr = np.linspace(0, len(fCompr), len(fCompr))

fitCentres_region = 300
DIAMETER = 4503
H0 = 600

HR = 100

def compressionFitChadwick_weightedLinearFit(hCompr, fCompr, fitCentres_region, H0, DIAMETER):

    error = False
    
    def computeStress(f, h, H0):
        R = DIAMETER/2000
        delta = (H0 - h)/1000
        stress = f / (np.pi * R * delta)
        return(stress)
        
    def computeStrain(h, H0):
        delta = (H0 - h)/1000
        strain = delta / (3 * (H0/1000))
        return(strain)
    
    def constitutiveRelation(strain, K, stress0):
        stress = (K * strain) + stress0
        return(stress)
    
    def gaussianWeights(fitCentres_region, stressCompr, HR):
        stressCompr = stressCompr.flatten(order='C')
        return np.exp ( -((stressCompr - fitCentres_region) ** 2) / HR ** 2)
    
    
    def weightedLinearFit(strainCompr, stressCompr, fitCentres_region, HR = 100):
        strainCompr = strainCompr.reshape(-1, 1)
        stressCompr = stressCompr.reshape(-1, 1)
        print(fitCentres_region)
        weights = gaussianWeights(fitCentres_region, stressCompr, HR)
        # print(weights)
        regr = LinearRegression()
        regr.fit(strainCompr, stressCompr, weights)
        K = regr.coef_[0]
        stress0 = regr.intercept_
        
        # print(K)
        stressCompr_copy = np.copy(stressCompr).flatten()
        stressPredict = constitutiveRelation(stressCompr_copy, K, stress0)
        
        idx_local = np.where(np.logical_and(stressCompr_copy >= fitCentres_region-HR, \
                                            stressCompr_copy <= fitCentres_region+HR))
    
        
        stressPredict = stressPredict[idx_local]
        # plt.plot(strainPredict, stressPredict, label = regionFitName)
        
        # print(stressPredict)
        # print(strainPredict)
        return(K, stressPredict)
    
    strainCompr = computeStrain(hCompr, H0)
    stressCompr = computeStress(fCompr, hCompr, H0)
    
    try:
        K, stressPredict = weightedLinearFit(strainCompr, stressCompr, fitCentres_region)
        print('here')
    except:
        print('error')
        error = True
        K, stressPredict, strainPredict = -1, np.ones(len(stressCompr))*(-1), np.ones(len(stressCompr))*(-1)
    
    # print(K)
    print(stressPredict)
    return(K, stressPredict, error)

#%%
K, stressPredict, error = compressionFitChadwick_weightedLinearFit(hCompr, fCompr, fitCentres_region, H0, DIAMETER)

#%%

def constitutiveRelation(strain, K, stress0):
    stress = (K * strain) + stress0
    return(stress)

strainCompr = strainCompr.reshape(-1, 1)
stressCompr = stressCompr.reshape(-1, 1)
print(fitCentres_region)
weights = gaussianWeights(fitCentres_region, stressCompr, HR)
# print(weights)
regr = LinearRegression()
regr.fit(strainCompr, stressCompr, weights)
K = regr.coef_[0]
stress0 = regr.intercept_
strainCompr = strainCompr.flatten()

print(K)
print(stress0)
stressCompr_copy = np.copy(stressCompr).flatten()
stressPredict = constitutiveRelation(stressCompr_copy, K, stress0)

idx_local = np.where(np.logical_and(stressCompr_copy >= fitCentres_region-HR, \
                                    stressCompr_copy <= fitCentres_region+HR))

print(idx_local)
stressPredict = stressPredict[idx_local]
strainPredict = strainCompr[idx_local]
# plt.plot(strainPredict, stressPredict, label = regionFitName)

print(stressPredict)
print(strainPredict)

#%%

def constitutiveRelation(strain, K, stress0):
    stress = (K * strain) + stress0
    return(stress)

strainCompr = strainCompr.reshape(-1, 1)
stressCompr = stressCompr.reshape(-1, 1)
print(fitCentres_region)
weights = gaussianWeights(fitCentres_region, stressCompr, HR)
# print(weights)
regr = LinearRegression()
regr.fit(strainCompr, stressCompr, weights)
K = regr.coef_[0]
stress0 = regr.intercept_
strainCompr = strainCompr.flatten()

# print(K)
stressCompr_copy = np.copy(stressCompr).flatten()
stressPredict = constitutiveRelation(stressCompr_copy, K, stress0)

idx_local = np.where(np.logical_and(stressCompr_copy >= fitCentres_region-HR, \
                                    stressCompr_copy <= fitCentres_region+HR))


stressPredict = stressPredict[idx_local]
# plt.plot(strainPredict, stressPredict, label = regionFitName)

print(stressPredict)
# print(strainPredict)
# return(K, stressPredict)