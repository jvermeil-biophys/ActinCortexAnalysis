# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 13:37:51 2021

@author: JosephVermeil
"""

# import numpy as np
# import pandas as pd

# # %%
# A = [1,2,4,6]
# x = np.searchsorted(A, 1, 'left')
# print(x)
# print(A[x])

# # %%
# d = {'a' : [0,1], 'b' : [10,100], 'c' : ['A', 'B']}
# df = pd.DataFrame(d)
# testPath = 'D://test.csv'
# df.to_csv(testPath, sep = '\t', index = False)

# # %%

# # %%

# # %%

# %%

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.ion()
# print(matplotlib.__version__)

# %%

plt.close('all')
fig, ax = plt.subplots(1,1)
# fig.show()
fig2, ax2 = plt.subplots(1,1)
# fig2.show()
# plt.show(block = False)
# plt.waitforbuttonpress()

# %%

plt.close('all')

# %%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

dirPath = 'E://21.12.08_Patterns3T3_beadSizes//M2_M450_aSFL_R40'
fileName = 'ResAcquis1 .txt'
filePath = os.path.join(dirPath, fileName)
df = pd.read_csv(filePath, sep = '\t', header = None)
df.columns = ['B', 'Bmeas', 'Zpiezo', 'T']
T0 = df['T'].values[0]
df['T'] = df['T'] - T0
df['B'] = df['B'] * 3/5
df['Bmeas'] = df['Bmeas'] * 3/5