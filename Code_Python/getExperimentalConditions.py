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
    # Getting the table
    experimentalDataFile = 'ExperimentalConditions.csv'
    experimentalDataFilePath = os.path.join(experimentalDataDir, experimentalDataFile)
    expConditionsDF = pd.read_csv(experimentalDataFilePath, sep=sep, header=0)
    print(BLUE + 'Extracted a table with ' + str(expConditionsDF.shape[0]) + ' lines and ' + str(expConditionsDF.shape[1]) + ' columns.' + NORMAL)
    
    # Cleaning the table
#     try:
    for c in expConditionsDF.columns:
        if 'Unnamed' in c:
            expConditionsDF = expConditionsDF.drop([c], axis=1)
        if '.1' in c:
            expConditionsDF = expConditionsDF.drop([c], axis=1)
    expConditionsDF = expConditionsDF.convert_dtypes()

    listTextColumns = []
    for col in expConditionsDF.columns:
        try:
            if expConditionsDF[col].dtype == 'string':
                listTextColumns.append(col)
        except:
            aaaa=0
            #Ok

    expConditionsDF[listTextColumns] = expConditionsDF[listTextColumns].apply(lambda x: x.str.replace(',','.'))

    expConditionsDF['scale pixel per um'] = expConditionsDF['scale pixel per um'].astype(float)
    try:
        expConditionsDF['optical index correction'] = \
                  expConditionsDF['optical index correction'].apply(lambda x: x.split('/')[0]).astype(float) \
                / expConditionsDF['optical index correction'].apply(lambda x: x.split('/')[1]).astype(float)
    except:
        print('optical index correction already in ' + str(expConditionsDF['optical index correction'].dtype) + ' type.')

    expConditionsDF['magnetic field correction'] = expConditionsDF['magnetic field correction'].astype(float)
    expConditionsDF['with fluo images'] = expConditionsDF['with fluo images'].astype(bool)

    try:
        expConditionsDF['ramp field'] = \
        expConditionsDF['ramp field'].apply(lambda x: [x.split(';')[0], x.split(';')[1]] if not pd.isnull(x) else [])
    except:
        aaaa=0
        #Ok

    dateExemple = expConditionsDF.loc[expConditionsDF.index[1],'date']

    if re.match(dateFormatExcel, dateExemple):
        print('dates corrected')
        expConditionsDF.loc[1:,'date'] = expConditionsDF.loc[1:,'date'].apply(lambda x: x.split('/')[0] + '-' + x.split('/')[1] + '-' + x.split('/')[2][2:])        
    elif re.match(dateFormatExcel2, dateExemple):
        print('dates corrected')
        expConditionsDF.loc[1:,'date'] = expConditionsDF.loc[1:,'date'].apply(lambda x: x.split('-')[0] + '-' + x.split('-')[1] + '-' + x.split('-')[2][2:])        
#     except:
#         print('Unexpected bug with the cleaning step')

    if save:
        saveName = 'ExperimentalConditions.csv'
        savePath = os.path.join(experimentalDataDir, saveName)
        expConditionsDF.to_csv(savePath, sep=';')

    expConditionsDF['manipID'] = expConditionsDF['date'] + '_' + expConditionsDF['manip']
#     reorgaList = np.array([i for i in range(len(expConditionsDF.columns))])
#     reorgaList[2] = reorgaList[-1]
#     reorgaList[3:] = reorgaList[3:] - np.ones(len(reorgaList)-3)
#     expConditionsDF = expConditionsDF[expConditionsDF.columns[reorgaList]]
    
    return(expConditionsDF)