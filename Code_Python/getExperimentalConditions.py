# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 15:00:27 2022

@author: JosephVermeil
"""

import numpy as np
import pandas as pd
import os
import re


NORMAL  = '\033[0m'
RED  = '\033[31m' # red
GREEN = '\033[32m' # green
ORANGE  = '\033[33m' # orange
BLUE  = '\033[36m' # blue


dateFormatExcel = re.compile('\d{2}/\d{2}/\d{4}')
dateFormatExcel2 = re.compile('\d{2}-\d{2}-\d{4}')
dateFormatOk = re.compile('\d{2}-\d{2}-\d{2}')
dateFormatExcel2 = re.compile('\d{2}-\d{2}-\d{4}')

def getExperimentalConditions(experimentalDataDir, save = False, sep = ';'):
    """"
    Import the table with all the conditions in a clean way.
    It is a tedious function to read because it's doing a boring job:
    Converting strings into numbers when possible
    Converting commas into dots to correct for the French decimal notation
    Converting semicolon separated values into lists when needed
    Etc
    """
    #### 0. Import the table
    experimentalDataFile = 'ExperimentalConditions.csv'
    experimentalDataFilePath = os.path.join(experimentalDataDir, experimentalDataFile)
    expConditionsDF = pd.read_csv(experimentalDataFilePath, sep=sep, header=0)
    print(BLUE + 'Importing Experimental Conditions' + NORMAL)
    print(BLUE + 'Extracted a table with ' + str(expConditionsDF.shape[0]) + ' lines and ' + str(expConditionsDF.shape[1]) + ' columns' + NORMAL)
    #### 1. Clean the table
    
    #### 1.1 Remove useless columns
    for c in expConditionsDF.columns:
        if 'Unnamed' in c:
            expConditionsDF = expConditionsDF.drop([c], axis=1)
        if '.1' in c:
            expConditionsDF = expConditionsDF.drop([c], axis=1)
    expConditionsDF = expConditionsDF.convert_dtypes()

    #### 1.2 Convert commas into dots
    listTextColumns = []
    for col in expConditionsDF.columns:
        try:
            if expConditionsDF[col].dtype == 'string':
                listTextColumns.append(col)
        except:
            pass
    expConditionsDF[listTextColumns] = expConditionsDF[listTextColumns].apply(lambda x: x.str.replace(',','.'))

    #### 1.3 Format 'scale'
    expConditionsDF['scale pixel per um'] = expConditionsDF['scale pixel per um'].astype(float)
    
    #### 1.4 Format 'optical index correction'
    try: # In case the format is 'n1/n2'
        expConditionsDF['optical index correction'] = \
                  expConditionsDF['optical index correction'].apply(lambda x: x.split('/')[0]).astype(float) \
                / expConditionsDF['optical index correction'].apply(lambda x: x.split('/')[1]).astype(float)
        print(ORANGE + 'optical index correction : format changed' + NORMAL)
    except:
        pass
    
    #### 1.5 Format 'magnetic field correction'
    expConditionsDF['magnetic field correction'] = expConditionsDF['magnetic field correction'].astype(float)
    
    #### 1.6 Format 'with fluo images'
    expConditionsDF['with fluo images'] = expConditionsDF['with fluo images'].astype(bool)

    # #### 1.7 Format 'ramp field'
    # try:
    #     print(ORANGE + 'ramp field : converted to list successfully' + NORMAL)
    #     expConditionsDF['ramp field'] = \
    #     expConditionsDF['ramp field'].apply(lambda x: [x.split(';')[0], x.split(';')[1]] if not pd.isnull(x) else [])
    # except:
    #     pass

    #### 1.8 Format 'date'
    dateExemple = expConditionsDF.loc[expConditionsDF.index[1],'date']
    if re.match(dateFormatExcel, dateExemple):
        print(ORANGE + 'dates : format corrected' + NORMAL)
        expConditionsDF.loc[:,'date'] = expConditionsDF.loc[:,'date'].apply(lambda x: x.split('/')[0] + '-' + x.split('/')[1] + '-' + x.split('/')[2][2:])        
    elif re.match(dateFormatExcel2, dateExemple):
        print(ORANGE + 'dates : format corrected' + NORMAL)
        expConditionsDF.loc[:,'date'] = expConditionsDF.loc[:,'date'].apply(lambda x: x.split('-')[0] + '-' + x.split('-')[1] + '-' + x.split('-')[2][2:])  


    #### 2. Save the table, if required
    if save:
        saveName = 'ExperimentalConditions.csv'
        savePath = os.path.join(experimentalDataDir, saveName)
        expConditionsDF.to_csv(savePath, sep=';')

    #### 3. Generate additionnal field that won't be saved
    
    def str2int(s):
        try:
            x = int(s)
        except:
            x = np.nan
        return(x)
    
    def str2float(s):
        try:
            x = float(s)
        except:
            x = np.nan
        return(x)
    
    #### 3.1 Make 'manipID'
    expConditionsDF['manipID'] = expConditionsDF['date'] + '_' + expConditionsDF['manip']
    
    # #### 3.2 Format 'bead diameter'
    # diameters = expConditionsDF.loc[:,'bead diameter'].apply(lambda x: str(x).split('_'))
    # diameters = diameters.apply(lambda x: [int(xx) for xx in x])
    # expConditionsDF.loc[:,'bead diameter'] = diameters
    # # print(ORANGE + 'ramp field : converted to list successfully' + NORMAL)
    
    # #### 3.3 Format 'bead type'
    # bt = expConditionsDF.loc[:,'bead type'].apply(lambda x: str(x).split('_'))
    # bt = bt.apply(lambda x: [str(xx) for xx in x])
    # expConditionsDF.loc[:,'bead type'] = bt
    
    # #### 3.4 Format 'ramp field'
    # rf = expConditionsDF.loc[:,'ramp field'].apply(lambda x: str(x).split('_'))
    # rf = rf.apply(lambda x: [str2float(xx) for xx in x])
    # expConditionsDF.loc[:,'ramp field'] = rf
    
    # #### 3.5 Format 'loop structure'
    # ls = expConditionsDF.loc[:,'loop structure'].apply(lambda x: str(x).split('_'))
    # ls = ls.apply(lambda x: [str2int(xx) for xx in x])
    # expConditionsDF.loc[:,'loop structure'] = ls

    #### 4. END
    return(expConditionsDF)