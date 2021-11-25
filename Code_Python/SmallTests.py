# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 13:37:51 2021

@author: JosephVermeil
"""

import numpy as np
import pandas as pd

# %%
A = [1,2,4,6]
x = np.searchsorted(A, 1, 'left')
print(x)
print(A[x])

# %%
d = {'a' : [0,1], 'b' : [10,100], 'c' : ['A', 'B']}
df = pd.DataFrame(d)
testPath = 'D://test.csv'
df.to_csv(testPath, sep = '\t', index = False)

# %%

# %%

# %%
