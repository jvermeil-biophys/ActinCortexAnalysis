# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:50:16 2021

@author: JosephVermeil
"""

# %% (0) Imports and settings

# 1. Imports
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import matplotlib.pyplot as plt

import os
import re
import time
import pyautogui
import matplotlib
import traceback
# import cv2

# import scipy
from scipy import interpolate
from scipy import signal

# import skimage
from skimage import io, filters, exposure, measure, transform, util, color
from scipy.signal import find_peaks, savgol_filter
from scipy.optimize import linear_sum_assignment
from matplotlib.gridspec import GridSpec

# Add the folder to path
COMPUTERNAME = os.environ['COMPUTERNAME']
if COMPUTERNAME == 'ORDI-JOSEPH':
    mainDir = "C://Users//JosephVermeil//Desktop//ActinCortexAnalysis"
elif COMPUTERNAME == 'LARISA':
    mainDir = "C://Users//Joseph//Desktop//ActinCortexAnalysis"
elif COMPUTERNAME == '':
    mainDir = "C://Users//josep//Desktop//ActinCortexAnalysis"

import sys
sys.path.append(mainDir + "//Code_Python")
from getExperimentalConditions import getExperimentalConditions

# 2. Pandas settings
pd.set_option('mode.chained_assignment', None)

# 3. Plot settings
# Here we use this mode because displaying images 
# in new windows is more convenient for this code.
# %matplotlib qt 
# matplotlib.use('Qt5Agg')
# To switch back to inline display, use : 
# %matplotlib widget or %matplotlib inline
# matplotlib.rcParams.update({'figure.autolayout': True})

SMALLER_SIZE = 8
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# 4. Other settings
# These regex are used to correct the stupid date conversions done by Excel
dateFormatExcel = re.compile('\d{2}/\d{2}/\d{4}')
dateFormatOk = re.compile('\d{2}-\d{2}-\d{2}')

# 5. Global constants
SCALE_100X = 15.8 # pix/µm 
NORMAL  = '\033[0m'
RED  = '\033[31m' # red
GREEN = '\033[32m' # green
ORANGE  = '\033[33m' # orange
BLUE  = '\033[36m' # blue



# %% (1) Utility functions

def findInfosInFileName(f, infoType):
    """
    Return a given type of info from a file name.
    Inputs : f (str), the file name.
             infoType (str), the type of info wanted.
             infoType can be equal to : 
             * 'M', 'P', 'C' -> will return the number of manip (M), well (P), or cell (C) in a cellID.
             ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'C', the function will return 8.
             * 'manipID'     -> will return the full manip ID.
             ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'manipID', the function will return '21-01-18_M2'.
             * 'cellID'     -> will return the full cell ID.
             ex : if f = '21-01-18_M2_P1_C8.tif' and infoType = 'cellID', the function will return '21-01-18_M2_P1_C8'.
    """
    if infoType in ['M', 'P', 'C']:
        acceptedChar = [str(i) for i in range(10)] + ['.', '-']
        string = '_' + infoType
        iStart = re.search(string, f).end()
        i = iStart
        infoString = '' + f[i]
        while f[i+1] in acceptedChar and i < len(f)-1:
            i += 1
            infoString += f[i]
            
    elif infoType == 'date':
        datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
        date = f[datePos.start():datePos.end()]
        infoString = date
    
    elif infoType == 'manipID':
        datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
        date = f[datePos.start():datePos.end()]
        manip = 'M' + findInfosInFileName(f, 'M')
        infoString = date + '_' + manip
        
    elif infoType == 'cellID':
        datePos = re.search(r"[\d]{1,2}-[\d]{1,2}-[\d]{2}", f)
        date = f[datePos.start():datePos.end()]
        infoString = date + '_' + 'M' + findInfosInFileName(f, 'M') + \
                            '_' + 'P' + findInfosInFileName(f, 'P') + \
                            '_' + 'C' + findInfosInFileName(f, 'C')

    
    return(infoString)

def isFileOfInterest(f, manips, wells, cells):
    """
    Determine if a file f correspond to the given criteria.
    More precisely, return a boolean saying if the manip, well and cell number are in the given range.
    f is a file name. Each of the fields 'manips', 'wells', 'cells' can be either a number, a list of numbers, or 'all'.
    Example : if f = '21-01-18_M2_P1_C8.tif'
    * manips = 'all', wells = 'all', cells = 'all' -> the function return True.
    * manips = 1, wells = 'all', cells = 'all' -> the function return False.
    * manips = [1, 2], wells = 'all', cells = 'all' -> the function return True.
    * manips = [1, 2], wells = 2, cells = 'all' -> the function return False.
    * manips = [1, 2], wells = 1, cells = [5, 6, 7, 8] -> the function return True.
    Note : if manips = 'all', the code will consider that wells = 'all', cells = 'all'.
           if wells = 'all', the code will consider that cells = 'all'.
           This means you can add filters only in this order : manips > wells > cells.
    """
    test = False
    if f.endswith(".tif"):
        if manips == 'all':
            test = True
        else:
            try:
                manips_str = [str(i) for i in manips]
            except:
                manips_str = [str(manips)]
            infoM = findInfosInFileName(f, 'M')
            if infoM in manips_str:
                if wells == 'all':
                    test = True
                else:
                    try:
                        wells_str = [str(i) for i in wells]
                    except:
                        wells_str = [str(wells)]
                    infoP = findInfosInFileName(f, 'P')
                    if infoP in wells_str:
                        if cells == 'all':
                            test = True
                        else:
                            try:
                                cells_str = [str(i) for i in cells]
                            except:
                                cells_str = [str(cells)]
                            infoC = findInfosInFileName(f, 'C')
                            if infoC in cells_str:
                                test = True
    return(test)

def compute_cost_matrix(XY1,XY2):
    """
    Compute a custom cost matrix between two arrays of XY positions.
    Here the costs are simply the squared distance between each XY positions.
    Example : M[2,1] is the sqaured distance between XY1[2] and XY2[1], 
    which is ((XY2[1,1]-XY1[2,1])**2 + (XY2[1,0]-XY1[2,0])**2)
    """
    N1, N2 = XY1.shape[0],XY2.shape[0]
    M = np.zeros((N1, N2))
    for i in range(N1):
        for j in range(N2):
            M[i,j] = (np.sum((XY2[j,:] - XY1[i,:]) ** 2))
    return(M)

def ui2array(uixy):
    """
    Translate the output of the function plt.ginput() 
    (which are lists of tuples), in an XY array with this shape:
    XY = [[x0, y0], [x1, y1], [x2, y2], ...]
    So if you need the [x, y] of 1 specific point, call XY[i]
    If you need the list of all x coordinates, call XY[:, 0]
    """
    n = len(uixy)
    XY = np.zeros((n, 2))
    for i in range(n):
        XY[i,0], XY[i,1] = uixy[i][0], uixy[i][1]
    return(XY)

def getROI(roiSize, x0, y0, nx, ny):
    """
    Return coordinates of top left (x1, y1) and bottom right (x2, y2) corner of a ROI, 
    and a boolean validROI that says if the ROI exceed the limit of the image.
    Inputs : 
    - roiSize, the width of the (square) ROI.
    - x0, y0, the position of the central pixel.
    - nx, ny, the size of the image.
    Note : the ROI is done so that the final width (= height) 
    of the ROI will always be an odd number.
    """
    roiSize += roiSize%2
    x1 = int(np.floor(x0) - roiSize*0.5) - 1
    x2 = int(np.floor(x0) + roiSize*0.5)
    y1 = int(np.floor(y0) - roiSize*0.5) - 1
    y2 = int(np.floor(y0) + roiSize*0.5)
    if min([x1,nx-x2,y1,ny-y2]) < 0:
        validROI = False
    else:
        validROI = True
    return(x1, y1, x2, y2, validROI)

def getDepthoCleanSize(D, scale):
    """
    Function that looks stupid but is quite important ! It allows to standardise 
    across all other functions the way the depthograph width is computed.
    D here is the approximative size of the bead in microns, 4.5 for M450, 2.7 for M270.
    Scale is the pixel to microns ration of the objective.
    """
    cleanSize = int(np.floor(1*D*scale))
    cleanSize += 1 + cleanSize%2
    return(cleanSize)

def squareDistance_V0(M, V, normalize = False): # MAKE FASTER !!!
    """
    DEPRECATED BECAUSE TOO SLOW
    Compute a distance between two arrays of the same size, defined as such:
    D = integral of the squared difference between the two arrays.
    It is used to compute the best fit of a slice of a bead profile on the depthograph.
    This function speed is critical for the Z computation process because it is called so many times !
    """
    top = time.time()
    n, m = M.shape[0], M.shape[1]
    # len(V) should be m
    result = np.zeros(n)
    if normalize:
        V = V/np.mean(V)
    for i in range(n):
        if normalize:
            Mi = M[i,:]/np.mean(M[i,:])
        else:
            Mi = M[i,:]
        d = np.sum((Mi-V)**2)
        result[i] = d
    print('DistanceCompTime')
    print(time.time()-top)
    return(result)

def squareDistance(M, V, normalize = False): # MUCH FASTER ! **Michael Scott Voice** VERRRRY GOODE
    """
    Compute a distance between two arrays of the same size, defined as such:
    D = integral of the squared difference between the two arrays.
    It is used to compute the best fit of a slice of a bead profile on the depthograph.
    This function speed is critical for the Z computation process because it is called so many times !
    What made that function faster is the absence of 'for' loops and the use of np.repeat().
    """
    #     top = time.time()
    n, m = M.shape[0], M.shape[1]
    # len(V) should be m
    if normalize:
        V = V/np.mean(V)
    V = np.array([V])
    MV = np.repeat(V, n, axis = 0) # Key trick for speed !
    if normalize:
        M = (M.T/np.mean(M, axis = 1).T).T
    R = np.sum((M-MV)**2, axis = 1)
#     print('DistanceCompTime')
#     print(time.time()-top)
    return(R)

def matchDists(listD, listStatus, Nup, NVox, direction):
    """
    This function transform the different distances curves computed for 
    a Nuplet of images to match their minima. By definition it is not used for singlets of images.
    In practice, it's a tedious and boring function.
    For a triplet of image, it will move the distance curve by NVox voxels to the left 
    for the first curve of a triplet, not move the second one, and move the third by NVox voxels to the right.
    The goal : align the 3 matching minima so that the sum of the three will have a clear global minimum.
    direction = 'upward'or 'downward'
    """
    N = len(listStatus)
    offsets = np.array(listStatus) - np.ones(N) * (Nup//2 + 1)
    offsets = offsets.astype(int)
    listD2 = []
    if direction == 'upward':
        for i in range(N):
            if offsets[i] < 0:
                shift = abs(offsets[i])*NVox
                D = listD[i]
                fillVal = D[-1]
                D2 = np.concatenate((D[shift:],fillVal*np.ones(shift))).astype(np.uint16)
                listD2.append(D2)
            if offsets[i] == 0:
                D = listD[i].astype(np.uint16)
                listD2.append(D)
            if offsets[i] > 0:
                shift = abs(offsets[i])*NVox
                D = listD[i]
                fillVal = D[0]
                D2 = np.concatenate((fillVal*np.ones(shift),D[:-shift])).astype(np.uint16)
                listD2.append(D2)
    elif direction == 'downward':
        for i in range(N):
            if offsets[i] > 0:
                shift = abs(offsets[i])*NVox
                D = listD[i]
                fillVal = D[-1]
                D2 = np.concatenate((D[shift:],fillVal*np.ones(shift))).astype(np.uint16)
                listD2.append(D2)
            if offsets[i] == 0:
                D = listD[i].astype(np.uint16)
                listD2.append(D)
            if offsets[i] < 0:
                shift = abs(offsets[i])*NVox
                D = listD[i]
                fillVal = D[0]
                D2 = np.concatenate((fillVal*np.ones(shift),D[:-shift])).astype(np.uint16)
                listD2.append(D2)
    return(np.array(listD2))

def uiThresholding(I, method = 'otsu', factorT = 0.8):
    """
    Interactive thresholding function to replace IJ.
    Compute an auto thresholding on a global 3D image with a method from this list:
    > 'otsu', 'max_entropy', (add the method you want ! here are the options : https://scikit-image.org/docs/stable/api/skimage.filters.html )
    Then display a figure for the user to assess the threshold fitness, and according to the user choice,
    confirm the threshold or recompute it with new parameters in a recursive way.
    """
    # 1. Compute the threshold
    nz = I.shape[0]
    if method == 'otsu':
        threshold = factorT*filters.threshold_otsu(I)
    elif method == 'max_entropy':
        bitDepth = util.dtype_limits(I)[1]+1
        I8 = util.img_as_ubyte(I)
        threshold = factorT*max_entropy_threshold(I8)*(bitDepth/2**8)
        
    # 2. Display images for the user to assess the fitness
    # New version of the plot
        nS = I.shape[0]
        loopSize = nS//4
        N = min(4, nS//loopSize)
        L_I_plot = [I[loopSize*2*k + 2] for k in range(N)]
        L_I_thresh = [I_plot > threshold for I_plot in L_I_plot]
        for i in range(N):
            I_plot = L_I_plot[i]
            I_thresh = L_I_thresh[i]
            I_plot = util.img_as_ubyte(I_plot)
            I_plot = color.gray2rgb(I_plot)
            pStart, pStop = np.percentile(I_plot, (1, 99))
            I_plot = exposure.rescale_intensity(I_plot, in_range=(pStart, pStop))
            red_multiplier = [255, 0, 0]
            I_plot[I_thresh] = red_multiplier
            L_I_plot[i] = I_plot
            
        I_thresh_all = I > threshold
        I_thresh_max = np.max(I_thresh_all, axis = 0)
        
        fig = plt.figure(tight_layout=True)
        gs = GridSpec(2, 4, figure=fig)
        ax = []
        for i in range(N):
            ax.append(fig.add_subplot(gs[i//2, i%2]))
            ax[-1].imshow(L_I_plot[i])
            ax[-1].set_title('Frame ' + str(loopSize*2*i + 2) + '/' + str(nS), fontsize = 8)
            ax[-1].axes.xaxis.set_ticks([])
            ax[-1].axes.yaxis.set_ticks([])
        ax.append(fig.add_subplot(gs[:, 2:]))
        ax[-1].imshow(I_thresh_max, cmap = 'gray')
        ax[-1].set_title('Max projection', fontsize = 10)
        ax[-1].axes.xaxis.set_ticks([])
        ax[-1].axes.yaxis.set_ticks([])
        fig.suptitle(str(threshold), fontsize = 12)
        fig.show()
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(50, 380, 1800, 650)
    
    # 3. Ask the question to the user
    QA = pyautogui.confirm(
                text='Is the threshold satisfying?',
                title='Confirm threshold', 
                buttons=['Yes', '10% Lower', '5% Lower', '1% Lower', '1% Higher', '5% Higher', '10% Higher'])
    plt.close(fig)
    
    # 4. Recall the same function with new parameters, or validate the threshold 
    # according to the user answer.
    increment = 0.1 * ('10%' in QA) + 0.05 * ('5%' in QA) + 0.01 * ('1%' in QA)
    if 'Lower' in QA:
        uiThresholding(method = method, factorT = factorT - increment)
    elif 'Higher' in QA:
        uiThresholding(method = method, factorT = factorT + increment)
    elif QA == 'Yes':
        threshold = threshold
    return(threshold)

def max_entropy(data):
    """
    Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method
    Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for Gray-Level Picture Thresholding Using the Entropy
    of the Histogram", Graphical Models and Image Processing, 29(3): 273-285
    M. Emre Celebi
    06.15.2007
    Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
    2016-04-28: Adapted for Python 2.7 by Robert Metchev from Java source of MaxEntropy() in the Autothresholder plugin
    http://rsb.info.nih.gov/ij/plugins/download/AutoThresholder.java
    :param data: Sequence representing the histogram of the image
    :return threshold: Resulting maximum entropy threshold
    """

    # calculate CDF (cumulative density function)
    cdf = data.astype(np.float).cumsum()

    # find histogram's nonzero area
    valid_idx = np.nonzero(data)[0]
    first_bin = valid_idx[0]
    last_bin = valid_idx[-1]

    # initialize search for maximum
    max_ent, threshold = 0, 0

    for it in range(first_bin, last_bin + 1):
        # Background (dark)
        hist_range = data[:it + 1]
        hist_range = hist_range[hist_range != 0] / cdf[it]  # normalize within selected range & remove all 0 elements
        tot_ent = -np.sum(hist_range * np.log(hist_range))  # background entropy

        # Foreground/Object (bright)
        hist_range = data[it + 1:]
        # normalize within selected range & remove all 0 elements
        hist_range = hist_range[hist_range != 0] / (cdf[last_bin] - cdf[it])
        tot_ent -= np.sum(hist_range * np.log(hist_range))  # accumulate object entropy

        # find max
        if tot_ent > max_ent:
            max_ent, threshold = tot_ent, it

    return(threshold)

def max_entropy_threshold(I):
    """
    Function based on the previous one that directly takes an image for argument.
    """
    H, bins = exposure.histogram(I, nbins=256, source_range='image', normalize=False)
    T = max_entropy(H)
    return(T)


# %% (2) Tracker classes
    
# experimentalDataDir = "C:/Users/anumi/OneDrive/Desktop/ActinCortexAnalysis/Data_Experimental"
# expDf = getExperimentalConditions()

    
# %%%% PincherTimeLapse
    
class PincherTimeLapse:
    """
    This class is initialised for each new .tif file analysed.
    
    It requires the following inputs :
    > I : the timelapse analysed.
    > cellID : the id of the cell currently analysed.
    > manipDict : the line of the experimental data table that concerns the current experiment.
    > NB : the number of beads of interest that will be tracked.
    
    It contains:
    * data about the 3D image I (dimensions = time, height, width),
    * a list of Frame objects listFrames, 1 per frame in the timelapse,
    * a list of Trajectory objects listTrajectories, 1 per bead of interest (Boi) in the timelapse,
    * a dictionnary dictLog, saving the status_frame of each frame (see below) 
                             and all of the user inputs (points and clicks) during the tracking,
    * a pandas DataFrame detectBeadsResult, that contains the raw output of the bead tracking,
    * metadata about the experiment (cellID, expType, loopStruct, loop_totalSize, loop_rampSize, 
                                        loop_excludedSize, nLoop, Nuplet, blackFramesPerLoop).
    
    When a PincherTimeLapse is initialised, most of these variables are initialised to zero values.
    In order to compute the different fields, the following methods should be called in this order:
    - ptl.checkIfBlackFrames() : detect if there are black images at the end of each loop in the time lapse and 
                                 classify them as not relevant by filling the appropriate fields.
    - ptl.saveFluoAside() : save the fluo images in an other folder and classify them as not relevant 
                            for the rest of the image analysis.
    - ptl.determineFramesStatus_R40() : fill the status_frame and status_nUp column of the dictLog.
                                    in the status_frame field: -1 means excluded ; 0 means ramp ; >0 means *position in* the n-uplet
                                    in the status_nUp field: -1 means excluded ; 0 means ramp ; >0 means *number of* the n-uplet
    - ptl.uiThresholding() : Compute the threshold that will be used for segmentation in an interractive way.
    - ptl.saveMetaData() : Save the computed threshold along with a few other data.
    - ptl.makeFramesList() : Initialize the frame list.
    - ptl.detectBeads() : Detect all the beads or load their positions from a pre-existing '_Results.txt' file.
    - ptl.buildTrajectories() : Do the tracking of the beads of interest, with the user help, or load pre-existing trajectories.
    [In the meantime, Z computations and neighbours detections are performed on the Trajectory objects]
    - ptl.computeForces() : when the Trajectory objects are complete (with Z and neighbours), compute the forces. 
                            Include recent corrections to the formula [October 2021].
    """
    
    def __init__(self, I, cellID, manipDict, NB = 2):
        # 1. Infos about the 3D image. The shape of the 3D image should be the following: T, Y, X !
        nS, ny, nx = I.shape[0], I.shape[1], I.shape[2]
        self.I = I
        self.nx = nx
        self.ny = ny
        self.nS = nS
        
        # 2. Infos about the experimental conditions, mainly from the DataFrame 'manipDict'.
        self.NB = NB # The number of beads of interest ! Typically 2 for a normal experiment, 4 for a triple pincher !
        self.cellID = cellID
        self.wFluo = manipDict['with fluo images']
        self.expType = manipDict['experimentType']
        self.scale = manipDict['scale pixel per um']
        self.OptCorrFactor = manipDict['optical index correction']
        self.MagCorrFactor = manipDict['magnetic field correction']
        self.Nuplet = manipDict['normal field multi images']
        self.Zstep = manipDict['multi image Z step']
        
        self.BeadsZDelta = manipDict['beads bright spot delta']
        self.BeadTypeStr = manipDict['bead type']
        self.beadTypes = [bT for bT in str(manipDict['bead type']).split('_')]
        self.beadDiameters = [int(bD) for bD in str(manipDict['bead diameter']).split('_')]
        self.dictBeadDiameters = {}
        for k in range(len(self.beadTypes)):
            self.dictBeadDiameters[self.beadTypes[k]] = self.beadDiameters[k]
         
        loopStruct = manipDict['loop structure'].split('_')
        #### Exp type dependance here
        
        # This is an ugly but necessary part of the code
        # This loopStruct field contains from 1 to 3 numbers, separated by a '_'
        # The convention for now is : 'totalLoopSize_rampSize_excludedSize'
        # totalLoopSize > compulsary. The size of an entire loop of images.
        # rampSize > compulsary only for compressions exp. The number of images belonging to the compression per loop.
        # excludedSize > optional. Indicates if some images (eg. fluorescence ones) should be systematically excluded at the end of each loop.
        if 'compressions' in self.expType or 'thickness' in self.expType:
            self.loop_totalSize = int(loopStruct[0])
            
            if self.expType == 'compressions':
                self.loop_rampSize = int(loopStruct[1])
            elif self.expType == 'compressionsLowStart':
                self.loop_rampSize = int(loopStruct[1])
            else:
                self.loop_rampSize = 0
            
            
            if len(loopStruct) == 3: # This 3rd part of the 'loopStruct' field is the nb of frames at the end
            # of each loop which are not part of the main analysis and should be excluded. Typically fluo images.
                self.loop_excludedSize = int(loopStruct[2])
            else:
                self.loop_excludedSize = 0
            self.nLoop = int(np.round(nS/self.loop_totalSize))
        
        elif 'optoGen' in self.expType:
            self.optoStep1 = int(loopStruct[0])
            self.optoStep2 = int(loopStruct[1])
            self.optoStep3 = int(loopStruct[2])
        
        # 3. Field that are just initialized for now and will be filled by calling different methods.
        self.threshold = 0
        self.listFrames = []
        self.listTrajectories = []
        
        self.dictLog = {'Slice' : np.array([i+1 for i in range(nS)]),
                        'status_frame' : np.zeros(nS, dtype = float),  # in the status_frame field: -1 means excluded ; 0 means ramp ; >0 means position in the n-uplet
                        'status_nUp' : np.zeros(nS, dtype = int), # in the status_nUp field: -1 means excluded ; 0 means ramp ; >0 means number of the n-uplet
                        'UI' : np.zeros(nS, dtype = bool),
                        'UILog' : np.array(['' for i in range(nS)], dtype = '<U16'),
                        'UIxy' : np.zeros((nS,NB,2), dtype = int)}
        
        self.detectBeadsResult = pd.DataFrame({'Area' : [], 
                                               'StdDev' : [], 
                                               'XM' : [], 
                                               'YM' : [], 
                                               'Slice' : []})
        
        self.blackFramesPerLoop = np.zeros(self.nLoop)
        self.modeNoUIactivated = False
        
        # End of the initialization !
       
    def checkIfBlackFrames(self):
        """
        Check if some images in the time lapse are completely black. 
        This happens typically when the computer is not able to save 
        properly a series of large images with a high frequency.
        To detect them, compute the checkSum = np.sum(self.I[j]).
        Then modify the 'status_frame' & 'status_nUp' fields to '-1' in the dictLog.
        """
        for i in range(self.nLoop):
            j = ((i+1)*self.loop_totalSize) - 1
            checkSum = np.sum(self.I[j])
            while checkSum == 0:
#                 self.dictLog['Black'][j] = True
                self.dictLog['status_frame'][j] = -1
                self.dictLog['status_nUp'][j] = -1
                self.blackFramesPerLoop[i] += 1
                j -= 1
                checkSum = np.sum(self.I[j])
              
    def saveFluoAside(self, fluoDirPath, f):
        """
        If wFluo = True in the expDf, modify the 'status_frame' & 'status_nUp' fields to '-1' in the dictLog.
        And if the directory for the fluo images has not be created yet,
        find and save all of the fluo images there.
        """
        if self.wFluo:
            for i in range(self.nLoop):
                j = int(((i+1)*self.loop_totalSize) - 1 - self.blackFramesPerLoop[i])
                self.dictLog['status_frame'][j] = -1
                self.dictLog['status_nUp'][j] = -1
                
            if not os.path.exists(fluoDirPath):
                os.makedirs(fluoDirPath)
                for i in range(self.nLoop):
                    j = int(((i+1)*self.loop_totalSize) - 1 - self.blackFramesPerLoop[i])
                    Ifluo = self.I[j]
                    path = os.path.join(fluoDirPath, f[:-4] + '_Fluo_' + str(j) + '.tif')
                    io.imsave(path, Ifluo)
                
    def determineFramesStatus_R40(self):
        #### Exp type dependance here
        """
        Fill the status_frame and status_nUp column of the dictLog, in the case of a compression (R40) or constant field (thickness) experiment
        > in the status_frame field: -1 means excluded ; 0 means ramp ; 10 > x > 0 means *position in* the n-uplet.
        > in the status_nUp field: -1 means excluded ; 0 means ramp ; >0 means *number of* the n-uplet.
        Not very elegant but very confortable to work with.
        """
        N0 = self.loop_totalSize
        Nramp0 = self.loop_rampSize
        Nexclu = self.loop_excludedSize
        nUp = self.Nuplet
        N = N0 - Nexclu
        Nct = N - Nramp0
        i_nUp = 1
        for i in range(self.nLoop):
            jstart = int(i*N0)
            if Nramp0 == 0:
                for j in range(N):
                    self.dictLog['status_frame'][jstart + j] = 1 + j%self.Nuplet
                    self.dictLog['status_nUp'][jstart + j] = i_nUp + j//self.Nuplet
            else:
                Nramp = Nramp0-self.blackFramesPerLoop[i]
                for j in range(Nct//2): # Ct field before ramp
                    self.dictLog['status_frame'][jstart + j] = 1 + j%self.Nuplet
                    self.dictLog['status_nUp'][jstart + j] = i_nUp + j//self.Nuplet
                i_nUp = max(self.dictLog['status_nUp']) + 1
                jstart += int(Nct//2 + Nramp) # In the ramp it self, the two status stay equal to 0
                for j in range(Nct//2): # Ct field after ramp
                    self.dictLog['status_frame'][jstart + j] = 1 + j%self.Nuplet
                    self.dictLog['status_nUp'][jstart + j] = i_nUp + j//self.Nuplet
                i_nUp = max(self.dictLog['status_nUp']) + 1
                
    def determineFramesStatus_L40(self):
        #### Exp type dependance here
        """
        Fill the status_frame and status_nUp column of the dictLog, in the case of a compression with an initial decrease (L40)
        > in the status_frame field: -1 means excluded ; 0 means ramp ; 0.1 means pre-ramp ; 10 > x > 0 means *position in* the n-uplet.
        > in the status_nUp field: -1 means excluded ; 0 means not in a n-uplet ; x > 0 means *number of* the n-uplet.
        Not very elegant but very confortable to work with.
        """
        N0 = self.loop_totalSize
        Nramp0 = self.loop_rampSize
        Nexclu = self.loop_excludedSize
        nUp = self.Nuplet
        N = N0 - Nexclu
        Nct = N - 2*Nramp0
        i_nUp = 1
        print('N0, Nramp0, Nexclu, N, Nct')
        print(N0, Nramp0, Nexclu, N, Nct)
        for i in range(self.nLoop):
            jstart = int(i*N0)
            if Nramp0 == 0:
                for j in range(N):
                    self.dictLog['status_frame'][jstart + j] = 1 + j%self.Nuplet
                    self.dictLog['status_nUp'][jstart + j] = i_nUp + j//self.Nuplet
                    
            else:
                Nramp = Nramp0-self.blackFramesPerLoop[i]
                for j in range(Nct//2): # Ct field before ramp
                    self.dictLog['status_frame'][jstart + j] = 1 + j%self.Nuplet
                    self.dictLog['status_nUp'][jstart + j] = i_nUp + j//self.Nuplet
                jstart += int(Nct//2)
                for j in range(Nramp0): # Pre-ramp
                    self.dictLog['status_frame'][jstart + j] = 0.1
                    self.dictLog['status_nUp'][jstart + j] = 0
                i_nUp = max(self.dictLog['status_nUp']) + 1
                jstart += int(Nramp0 + Nramp) # In the ramp it self, the two status stay equal to 0
                for j in range(Nct//2): # Ct field after ramp
                    self.dictLog['status_frame'][jstart + j] = 1 + j%self.Nuplet
                    self.dictLog['status_nUp'][jstart + j] = i_nUp + j//self.Nuplet
                i_nUp = max(self.dictLog['status_nUp']) + 1            
                
    def determineFramesStatus_optoGen(self):
        #### Exp type dependance here
        """
        
        """
        pass
                
                
    def saveLog(self, display = 1, save = False, path = ''):
        """
        Save the dictLog so that next time it can be directly reloaded to save time.
        """
        dL = {}
        dL['Slice'], dL['status_frame'], dL['status_nUp'] = \
            self.dictLog['Slice'], self.dictLog['status_frame'], self.dictLog['status_nUp']

        dL['UI'], dL['UILog'] = \
            self.dictLog['UI'], self.dictLog['UILog']
        for i in range(self.NB):
            dL['UIx'+str(i+1)] = self.dictLog['UIxy'][:,i,0]
            dL['UIy'+str(i+1)] = self.dictLog['UIxy'][:,i,1]
        dfLog = pd.DataFrame(dL)
        if save:
            dfLog.to_csv(path, sep='\t')
        
        if display == 1:
            print('\n\n* Initialized Log Table:\n')
            print(dfLog)
        if display == 2:
            print('\n\n* Filled Log Table:\n')
            print(dfLog[dfLog['UI']])
        
        
        
        
    def importLog(self, path):
        """
        Import the dictLog.
        """
        dfLog = pd.read_csv(path, sep='\t')
        dL = dfLog.to_dict()
        self.dictLog['Slice'], self.dictLog['status_frame'], self.dictLog['status_nUp'] = \
            dfLog['Slice'].values, dfLog['status_frame'].values, dfLog['status_nUp'].values
        self.dictLog['UI'], self.dictLog['UILog'] = \
            dfLog['UI'].values, dfLog['UILog'].values
        for i in range(self.NB):
            xkey, ykey = 'UIx'+str(i+1), 'UIy'+str(i+1)
            self.dictLog['UIxy'][:,i,0] = dfLog[xkey].values
            self.dictLog['UIxy'][:,i,1] = dfLog[ykey].values
        
    def computeThreshold(self, method, factorT):
        """
        Compute the threshold with the chosen method, multiply it with factorT, and then save it as a field.
        CURRENTLY NOT USED IN THE MAIN FUNCTION.
        """
        threshold = 0
        end_z = 2*self.loop_totalSize
        if method == 'otsu':
            threshold = factorT*filters.threshold_otsu(self.I[:end_z])
        elif method == 'max_entropy':
            bitDepth = util.dtype_limits(self.I[:end_z])[1]+1
            I8 = util.img_as_ubyte(self.I[:end_z])
            threshold = factorT*max_entropy_threshold(I8)*(bitDepth/(2**8))
        self.threshold = threshold
        
        
    def testThresholding(self):
        """
        Display the 3rd slice of the time lapse, thresholded.
        Why the 3rd you ask ? Because it will often be a 'top' image from a triplet 
        {'bottom', 'middle', 'top'} and usually it is the one with the least intensity.
        CURRENTLY NOT USED IN THE MAIN FUNCTION.
        """
        I_test = self.I[self.nS//2]
        I_thresh = I_test > self.threshold
        fig, ax = plt.subplots(1,1)
        ax.imshow(I_thresh, cmap = 'gray')
        fig.show()
        
        
    def uiThresholding(self, method, factorT):
        """
        Interactive thresholding function to replace IJ.
        Compute an auto thresholding on a global 3D image with a method from this list:
        > 'otsu', 'max_entropy', (add the method you want ! here are the options : https://scikit-image.org/docs/stable/api/skimage.filters.html )
        Then display a figure for the user to assess the threshold fitness, and according to the user choice,
        confirm the threshold or recompute it with new parameters in a recursive way.
        """
        threshold = 0
        end_z = 2*self.loop_totalSize
        if method == 'otsu':
            threshold = factorT*filters.threshold_otsu(self.I[:end_z])
        elif method == 'max_entropy':
            bitDepth = util.dtype_limits(self.I[:end_z])[1]+1
            I8 = util.img_as_ubyte(self.I[:end_z])
            threshold = factorT*max_entropy_threshold(I8)*(bitDepth/(2**8))
            
        # New version of the plot
        loopSize = self.loop_totalSize
        N = min(4, self.nS//loopSize)
        L_I_plot = [self.I[loopSize*2*k + 2] for k in range(N)]
        L_I_thresh = [I_plot > threshold for I_plot in L_I_plot]
        for i in range(N):
            I_plot = L_I_plot[i]
            I_thresh = L_I_thresh[i]
            I_plot = util.img_as_ubyte(I_plot)
            I_plot = color.gray2rgb(I_plot)
            pStart, pStop = np.percentile(I_plot, (1, 99))
            I_plot = exposure.rescale_intensity(I_plot, in_range=(pStart, pStop))
            red_multiplier = [255, 0, 0]
            I_plot[I_thresh] = red_multiplier
            L_I_plot[i] = I_plot
            
        I_thresh_all = self.I > threshold
        I_thresh_max = np.max(I_thresh_all, axis = 0)
        
        fig = plt.figure(tight_layout=True)
        gs = GridSpec(2, 4, figure=fig)
        ax = []
        for i in range(N):
            ax.append(fig.add_subplot(gs[i//2, i%2]))
            ax[-1].imshow(L_I_plot[i])
            ax[-1].set_title('Frame ' + str(loopSize*2*i + 2) + '/' + str(self.nS), fontsize = 8)
            ax[-1].axes.xaxis.set_ticks([])
            ax[-1].axes.yaxis.set_ticks([])
        ax.append(fig.add_subplot(gs[:, 2:]))
        ax[-1].imshow(I_thresh_max, cmap = 'gray')
        ax[-1].set_title('Max projection', fontsize = 10)
        ax[-1].axes.xaxis.set_ticks([])
        ax[-1].axes.yaxis.set_ticks([])
        fig.suptitle(str(threshold), fontsize = 12)
        fig.show()
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(50, 380, 1800, 650)


        QA = pyautogui.confirm(
                    text='Is the threshold satisfying?',
                    title='Confirm threshold', 
                    buttons=['Yes', '10% Lower', '5% Lower', '1% Lower', '1% Higher', '5% Higher', '10% Higher'])
        plt.close(fig)
        increment = 0.1 * ('10%' in QA) + 0.05 * ('5%' in QA) + 0.01 * ('1%' in QA)
        if 'Lower' in QA:
            self.uiThresholding(method = method, factorT = factorT - increment)
        elif 'Higher' in QA:
            self.uiThresholding(method = method, factorT = factorT + increment)
        elif QA == 'Yes':
            self.threshold = threshold
            
    def saveMetaData(self, path):
        """
        Save the computed threshold along with a few other data.
        """
        dMD = {}
        dMD['cellID'] = self.cellID
        dMD['threshold'] = self.threshold
        dMD['NB'] = self.NB
        f = open(path, 'w')
        for k in dMD.keys():
            f.write(str(k))
            f.write('\t')
            f.write(str(dMD[k]))
            f.write('\n')
        f.close()
        
    def readMetaData(self, path, infoType):
        """
        Read the metadata file, containing the computed threshold along with a few other data.
        """
        f = open(path, 'r')
        f_lines = f.readlines()
        dMD = {}
        for line in f_lines:
            splitLine = line.split('\t')
            try:
                dMD[splitLine[0]] = int(splitLine[1])
            except:
                try:
                    dMD[splitLine[0]] = float(splitLine[1])
                except:
                    dMD[splitLine[0]] = splitLine[1]
        f.close()
        return(dMD[infoType])
             
    def makeFramesList(self):
        """
        Initialize the Frame objects and add them to the PTL.listFrames list.
        """
        for i in range(self.nS):
            status_frame = self.dictLog['status_frame'][i]
            status_nUp = self.dictLog['status_nUp'][i]
            # The Nup field of a slice is = to self.Nuplet if the status_frame indicates that the frmae is part of a multi image n-uplet
            # Otherwise the image is "alone", like in a compression, and therefore Nup = 1
            Nup = (self.Nuplet * (status_nUp > 0))  +  (1 * (status_nUp <= 0))
            if self.dictLog['status_frame'][i] >= 0:
                self.listFrames.append(Frame(self.I[i], i, self.NB, self.threshold, Nup, status_frame, status_nUp, self.scale))
    
    def detectBeads(self, resFileImported, display = False):
        """
        If no '_Results.txt' file has been previously imported, ask each Frame 
        object in the listFrames to run its Frame.detectBeads() method.
        Then concatenate the small 'Frame.resDf' to the big 'PTL.detectBeadsResult' DataFrame,
        so that in the end you'll get a DataFrame that has exactly the shape of a '_Results.txt' file made from IJ.
        *
        If a '_Results.txt' file has been previously imported, just assign to each Frame 
        object in the listFrames the relevant resDf DataFrame 
        (each resDf is, as said earlier, just a fragment of the PTL.detectBeadsResult).
        """
        for frame in self.listFrames: #[:3]:
            plot = 0
            
            # TEST #################################################################################
#             listPlot = [i for i in range(1685, 1700)]
#             if frame.iS in listPlot:
#                 plot = 1
            # TEST #################################################################################
            
            if not resFileImported:
                frame.detectBeads(plot)
                self.detectBeadsResult = pd.concat([self.detectBeadsResult, frame.resDf])
                
            else:
                resDf = self.detectBeadsResult.loc[self.detectBeadsResult['Slice'] == frame.iS+1]
                frame.resDf = resDf
                
            frame.makeListBeads()
            
        if not resFileImported:
            self.detectBeadsResult = self.detectBeadsResult.convert_dtypes()
            self.detectBeadsResult.reset_index(inplace=True)
            self.detectBeadsResult.drop(['index'], axis = 1, inplace=True)
        
        if display:
            print('\n\n* Detected Beads Result:\n')
            print(self.detectBeadsResult)

    def saveBeadsDetectResult(self, path):
        """
        Save the 'PTL.detectBeadsResult' DataFrame.
        """
        self.detectBeadsResult.to_csv(path, sep='\t')
    
    def importBeadsDetectResult(self, path=''):
        """
        Import the 'PTL.detectBeadsResult' DataFrame.
        """
        df = pd.read_csv(path, sep='\t')
        for c in df.columns:
            if 'Unnamed' in c:
                df.drop([c], axis = 1, inplace=True)
        self.detectBeadsResult = df
        
    def findBestStd_V0(self):
        """
        This ugly function is my best attempt to implement sth very simple in a robust way.
        In the 'status_frame' field, -1 means excluded image, 0 means image that isn't part of a N-uplet of images, and k>0 means position in the N-uplet of images.
        For each image in the N-uplet, I want to reconsititute this N-uplet (meaning the list of Nuplet consecutive images numbered from 1 to Nuplet, minus the images eventually with no beads detected).
        Then for each N-uplet of images, i want to find the max standard deviation and report its position because it's for the max std that the X and Y detection is the most precise.
        An exemple: with these inputs:
        Nuplet = 3
        status_frame = [1,2,3,0,0,0,1,2, 3, 1, 2, 3, 1, 2]
            iS = [0,1,2,3,4,5,6,7,11,12,13,14,15,16]
           std = [1,5,9,5,5,5,1,5, 8, 2, 6, 9, 2, 6]
        The function will return bestStd, a list of boolean with the same length.
        Where status_frame = 0, bestStd = True (the image is not part of a N-uplet, thus it need to be analysed regardless of its std).
        Where satus > 0, the function will cut the lists in N-uplet of max size 3:
        status_frame -> [1,2,3] ; [1,2] ; [ 3] ; [ 1, 2, 3] ; [ 1, 2]
            iS -> [0,1,2] ; [6,7] ; [11] ; [12,13,14] ; [15,16]
           std -> [1,5,9] ; [1,5] ; [ 8] ; [ 2, 6, 9] ; [ 2, 6]
             i -> [0,1,2] ; [6,7] ; [ 8] ; [ 9,10,11] ; [12,13]
        and then will find the best std in each of those fragment and put a True at that position (list 'i') in bestStd, so in this example, at i = 2, 7, 8, 11 ,13
        So the output is: bestStd = [False, False, True, True, True, True, False, True, True, False, False, True, False, True]
        """
        
        Nuplet = self.Nuplet
        status_frame = self.listTrajectories[0].dict['status_frame']
        iS = self.listTrajectories[0].dict['iS']
        nT = self.listTrajectories[0].nT
        std = np.zeros(nT)
        for i in range(0,self.NB):
            std += np.array(self.listTrajectories[i].dict['StdDev'])
        
        bestStd = np.zeros(nT, dtype = bool)

        current_Nup_status = []
        current_Nup_iS = []
        current_Nup_std = []
        current_Nup_i = []
        for i in range(nT):
            if status_frame[i] == 0:
                bestStd[i] = True
            else:
                if len(current_Nup_i) == 0:
                    current_Nup_status = [status_frame[i]]
                    current_Nup_iS = [iS[i]]
                    current_Nup_std = [std[i]]
                    current_Nup_i = [i]
                else:
                    if status_frame[i] > current_Nup_status[-1] and (iS[i]-current_Nup_iS[-1]) < Nuplet:
                        current_Nup_status.append(status_frame[i])
                        current_Nup_iS.append(iS[i])
                        current_Nup_std.append(std[i])
                        current_Nup_i.append(i)
                    else:
    #                     print(current_Nup_status, current_Nup_iS, current_Nup_std, current_Nup_i)
                        i_bestStdInCurrentNuplet = int(np.argmax(np.array(current_Nup_std)))
    #                     print(i_bestStdInCurrentNuplet)
                        i_bestStd = int(current_Nup_i[i_bestStdInCurrentNuplet])
                        bestStd[i_bestStd] = True
                        # then:
                        current_Nup_status = [status_frame[i]]
                        current_Nup_iS = [iS[i]]
                        current_Nup_std = [std[i]]
                        current_Nup_i = [i]

        # Need to do it one last time at the end
        i_bestStdInCurrentNuplet = int(np.argmax(np.array(current_Nup_std)))
        i_bestStd = int(current_Nup_i[i_bestStdInCurrentNuplet])
        bestStd[i_bestStd] = True

        return(bestStd)
    
    def findBestStd(self):
        """
        Simpler and better than findBestStd_V0 using the status_nUp column of the dictLog.
        ---
        For each frame of the timelapse that belongs to a N-uplet, I want to reconsititute this N-uplet 
        (meaning the list of 'Nup' consecutive images numbered from 1 to Nup, 
        minus the images eventually with no beads detected).
        Then for each N-uplet of images, i want to find the max standard deviation 
        and report its position because it's for the max std that the X and Y detection is the most precise.
        ---
        This is very easy thanks to the 'status_nUp', because it contains a different number for each N-Uplet.
        """
        
        Nup = self.Nuplet
        nT = self.listTrajectories[0].nT
        status_nUp = self.listTrajectories[0].dict['status_nUp']
        std = np.zeros(nT)
        for i in range(self.NB):
            std += np.array(self.listTrajectories[i].dict['StdDev'])
        
        bestStd = np.zeros(nT, dtype = bool)
        i = 0
        while i < nT:
            if status_nUp[i] == 0:
                bestStd[i] = True
                i += 1
            elif status_nUp[i] > 0:
                s2 = status_nUp[i]
                L = [i]
                j = 0
                while i+j < nT-1 and status_nUp[i+j+1] == s2: # lazy evaluation of booleans
                    j += 1
                    L.append(i+j)
                #print(L)    
                loc_std = std[L]
                i_bestStd = i + int(np.argmax(loc_std))
                bestStd[i_bestStd] = True
                L = []
                i = i + j + 1
                
        return(bestStd)
        
    
    
    def buildTrajectories(self):
        """
        The main tracking function.
        *
        Note about the naming conventions here: 
        - 'iF': index in the list of Frames ; 
        - 'iB': index in a list of Beads or a list of Trajectories ; 
        - 'iS': index of the slice in the image I (but here python starts with 0 and IJ starts with 1);
        - 'Boi' refers to the 'Beads of interest', ie the beads that are being tracked.
        """
        
        # 1. Initialize the BoI position in the first image where they can be detect, thanks to user input.
        init_iF = 0
        init_ok = False
        while not init_ok:
            init_iS = self.listFrames[init_iF].iS
            
            if not self.dictLog['UI'][init_iS]: # Nothing in the log yet
                self.listFrames[init_iF].show()
                mngr = plt.get_current_fig_manager()
                mngr.window.setGeometry(720, 50, 1175, 1000)
                QA = pyautogui.confirm(
                    text='Can you point the beads of interest\nin the image ' + str(init_iS + 1) + '?',
                    title='Initialise tracker', 
                    buttons=['Yes', 'Next Frame', 'Quit'])
                if QA == 'Yes':
                    init_ok = True
                    ui = plt.ginput(self.NB, timeout=0)
                    uiXY = ui2array(ui)
                    self.dictLog['UI'][init_iS] = True
                    self.dictLog['UILog'][init_iS] = 'init_' + QA
                    self.dictLog['UIxy'][init_iS] = uiXY
                elif QA == 'Next Frame':
                    self.dictLog['UI'][init_iS] = True
                    self.dictLog['UILog'][init_iS] = 'init_' + QA
                    init_iF += 1
                else:
                    fig = plt.gcf()
                    plt.close(fig)
                    return('Bug')
                
                fig = plt.gcf()
                plt.close(fig)
                
            else: # Action to do already in the log
                QA = self.dictLog['UILog'][init_iS]
                if QA == 'init_Yes':
                    init_ok = True
                    uiXY = self.dictLog['UIxy'][init_iS]
                elif QA == 'init_Next Frame':
                    init_iF += 1
                else:
                    print('Strange event in the tracking init')
        
        init_BXY = self.listFrames[init_iF].beadsXYarray()
        M = compute_cost_matrix(uiXY,init_BXY)
        row_ind, col_ind = linear_sum_assignment(M) # row_ind -> clicks / col_ind -> listBeads
        
        # Sort the initial beads to have them ordered by increasing x coord.
        # SOMETHING IS WRONG HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        sortM = np.array([[init_BXY[col_ind[i],0], col_ind[i]] for i in range(len(col_ind))])
        sortM = sortM[sortM[:, 0].argsort()]
        init_iBoi = sortM[:, 1].astype(int)
        init_BoiXY = np.array([init_BXY[col_ind[i]] for i in range(len(col_ind))])
        
        # 2. Initialize the BoI position in the first image where they can be detect, thanks to user input.
        
        # Creation of the Trajectory objects
        for iB in range(self.NB):
            self.listTrajectories.append(Trajectory(self.I, self.listFrames, self.scale, self.Zstep, iB))

            self.listTrajectories[iB].dict['Bead'].append(self.listFrames[init_iF].listBeads[init_iBoi[iB]])
            self.listTrajectories[iB].dict['iF'].append(init_iF)
            self.listTrajectories[iB].dict['iS'].append(self.listFrames[init_iF].iS)
            self.listTrajectories[iB].dict['iB_inFrame'].append(init_iBoi[iB])
            self.listTrajectories[iB].dict['X'].append(init_BoiXY[iB][0])
            self.listTrajectories[iB].dict['Y'].append(init_BoiXY[iB][1])
            self.listTrajectories[iB].dict['StdDev'].append(self.listFrames[init_iF].beadsStdDevarray()[init_iBoi[iB]])
            self.listTrajectories[iB].dict['status_frame'].append(self.listFrames[init_iF].status_frame)
            self.listTrajectories[iB].dict['status_nUp'].append(self.listFrames[init_iF].status_nUp)
            
            #### Exp type dependance here (01)
            if 'compressions' in self.expType or 'thickness' in self.expType:
                self.listTrajectories[iB].dict['idxAnalysis'].append(1 * (self.listFrames[init_iF].status_frame == 0))
                
            elif 'optoGen' in self.expType:
                self.listTrajectories[iB].dict['idxAnalysis'].append(0)

        # Start the tracking
        previous_iF = init_iF
        previous_iBoi = init_iBoi
        previous_BXY = init_BXY
        previous_BoiXY = init_BoiXY
        
        for iF in range(init_iF+1, len(self.listFrames)):
            validFrame = True
            
            if self.listFrames[iF].NBdetected < self.NB:
                validFrame = False
            
            if validFrame:
                BXY = self.listFrames[iF].beadsXYarray()
                M = compute_cost_matrix(previous_BoiXY,BXY)
                row_ind, col_ind = linear_sum_assignment(M)
                
                costs = [M[row_ind[iB],col_ind[iB]] for iB in range(len(row_ind))]
          
                # Assess wether the algo should aks for user input
                askUI = False
                if (max(costs)**0.5) * (1/self.scale) > 0.5: 
                # If the distance travelled by one of the BoI is greater than 0.5 um
                    
#                     print('M')
#                     print(M)
#                     print('row_ind, col_ind')
#                     print(row_ind, col_ind)
#                     print('previous_iBoi')
#                     print(previous_iBoi)
#                     print('costs')
#                     print(costs)
                    askUI = True
                

                if not askUI: # Automatically assign the positions of the next beads
                    try:
                        iBoi = [col_ind[iB] for iB in row_ind]
                        BoiXY = np.array([BXY[iB] for iB in iBoi])
                    except:
                        print('Error for ' + str(iF))
                        askUI = True
                        print('M')
                        print(M)
                        print('row_ind, col_ind')
                        print(row_ind, col_ind)
                        print('previous_iBoi')
                        print(previous_iBoi)
                        print('costs')
                        print(costs)
                        askUI = True

                if askUI: # Ask user input to asign the positions of the next beads
            
                    iS = self.listFrames[iF].iS
                    
                    # Case 1: the UI has been previously saved in the dictLog.
                    # Then just import the previous answer from the dictLog
                    if self.dictLog['UI'][iS]:
                        QA = self.dictLog['UILog'][iS]
                        if QA == 'Yes':
                            uiXY = self.dictLog['UIxy'][iS]
                        elif QA == 'No' or QA == 'No to all':
                            validFrame = False
                            #fig = plt.gcf()
                            #plt.close(fig)
                    
                    
                    # Case 2: the UI has not been previously saved in the dictLog
                    # Then ask for UI ; and save it in the dictLog
                    elif not self.dictLog['UI'][iS]:
                        if self.modeNoUIactivated == False:
                            # Display the image, plot beads positions and current trajectories & ask the question
                            self.listFrames[iF].show()
                            for iB in range(self.NB):
                                T = self.listTrajectories[iB]
                                ax = plt.gca()
                                T.plot(ax, iB)
                            
                            mngr = plt.get_current_fig_manager()
                            mngr.window.setGeometry(720, 50, 1175, 1000)
                            QA = pyautogui.confirm(
                                text='Can you point the beads of interest\nin the image ' + str(iS + 1) + '?',
                                title='', 
                                buttons=['No', 'Yes', 'Abort!', 'No to all'])
                            
                            # According to the question's answer:
                            if QA == 'Yes':
                                ui = plt.ginput(self.NB, timeout=0)
                                uiXY = ui2array(ui)
                                self.dictLog['UI'][iS] = True
                                self.dictLog['UILog'][iS] = QA
                                self.dictLog['UIxy'][iS] = uiXY
                            elif QA == 'No':
                                validFrame = False
                                self.dictLog['UI'][iS] = True
                                self.dictLog['UILog'][iS] = QA
                            elif QA == 'Abort!':
                                validFrame = False
                                fig = plt.gcf()
                                plt.close(fig)
                                return('Bug')
                            elif QA == 'No to all':
                                validFrame = False
                                self.modeNoUIactivated = True
                                self.dictLog['UI'][iS] = True
                                self.dictLog['UILog'][iS] = QA
                            fig = plt.gcf()
                            plt.close(fig)
                            
                        elif self.modeNoUIactivated == True:
                        # This mode is in case you don't want to keep clicking 'No' for hours when
                        # you know for a fact that there is nothing else you can do with this TimeLapse.
                            iS = self.listFrames[iF].iS
                            QA = 'No'
                            validFrame = False
                            self.dictLog['UI'][iS] = True
                            self.dictLog['UILog'][iS] = QA
                
                    
            
            # If there were more than NB objects, and then the QA wasn't 'No', then the frame is valid.
            if validFrame:
                M = compute_cost_matrix(uiXY,BXY)
                row_ind, col_ind = linear_sum_assignment(M)
                sortM = np.array([[BXY[col_ind[i],0], col_ind[i]] for i in range(len(col_ind))])
                sortM = sortM[sortM[:, 0].argsort()]
                iBoi = sortM[:, 1].astype(int)
                BoiXY = np.array([BXY[iB] for iB in iBoi])
                
                #### Exp type dependance here (02)
                if 'compressions' in self.expType or 'thickness' in self.expType:
                    if self.expType == 'compressions':
                        idxAnalysis = (self.listFrames[iF].status_frame == 0) \
                            * (max(self.listTrajectories[iB].dict['idxAnalysis']) \
                               + 1 * (self.listTrajectories[iB].dict['idxAnalysis'][-1] == 0))
                                
                    elif self.expType == 'compressionsLowStart': # a pre-ramp has the same idxAnalysis than a ramp but in negative.
                        idxAnalysis = (self.listFrames[iF].status_frame == 0) \
                            * (max(self.listTrajectories[iB].dict['idxAnalysis']) + 1 * (self.listTrajectories[iB].dict['idxAnalysis'][-1] <= 0)) \
                                - (self.listFrames[iF].status_frame == 0.1) \
                            * (abs(min(self.listTrajectories[iB].dict['idxAnalysis']) - 1 * (self.listTrajectories[iB].dict['idxAnalysis'][-1] == 0)))
                                
                elif 'optoGen' in self.expType:
                    self.listTrajectories[iB].dict['idxAnalysis'].append(0)
                
                # idxAnalysis = 0 if not in a ramp, and = number of ramp else. Basically increase by 1 each time you have an interval between two ramps.
                for iB in range(self.NB):
                    #
                    self.listTrajectories[iB].dict['Bead'].append(self.listFrames[iF].listBeads[iBoi[iB]])
                    self.listTrajectories[iB].dict['iF'].append(iF)
                    self.listTrajectories[iB].dict['iS'].append(self.listFrames[iF].iS)
                    self.listTrajectories[iB].dict['iB_inFrame'].append(iBoi[iB])
                    self.listTrajectories[iB].dict['X'].append(BoiXY[iB][0])
                    self.listTrajectories[iB].dict['Y'].append(BoiXY[iB][1])
                    self.listTrajectories[iB].dict['StdDev'].append(self.listFrames[iF].beadsStdDevarray()[iBoi[iB]])
                    self.listTrajectories[iB].dict['status_frame'].append(self.listFrames[iF].status_frame)
                    self.listTrajectories[iB].dict['status_nUp'].append(self.listFrames[iF].status_nUp)
                    self.listTrajectories[iB].dict['idxAnalysis'].append(idxAnalysis)
                    
                previous_iF = iF
                previous_iBoi = iBoi
                previous_BXY = BXY
                previous_BoiXY = BoiXY
            
            else: # else, go to the next frame 
                continue
                
        for iB in range(self.NB):
            for k in self.listTrajectories[iB].dict.keys():
                self.listTrajectories[iB].dict[k] = np.array(self.listTrajectories[iB].dict[k])
                
        # Now we have a functional Trajectory object
        # it's time to refine it
        
        nT = len(self.listTrajectories[0].dict['Bead'])
        
        #### Black Images deletion in the trajectory
        
        # (a) Add the pointer to the correct line of the _Field.txt file.
        # It's just exactly the iS already saved in the dict, except if there are black images at the end of loops.
        # In that case you have to skip the X lines corresponding to the end of the ramp part, X being the nb of black images at the end of the current loop
        # This is because when black images occurs, they do because of the high frame rate during ramp parts and thus replace these last ramp images.
        
        for iB in range(self.NB):
            self.listTrajectories[iB].dict['Zr'] = np.zeros(nT)
            self.listTrajectories[iB].nT = nT
            iField = []
            for i in range(nT):
                S = self.listTrajectories[iB].dict['iS'][i]
                iLoop = (S//self.loop_totalSize)
                offset = self.blackFramesPerLoop[iLoop]
                i_lim = iLoop*self.loop_totalSize + (self.loop_totalSize - ((self.loop_totalSize-self.loop_rampSize)//2) - (self.loop_excludedSize + offset))
                # i_lim is the first index after the end of the ramp
                addOffset = (S >= i_lim)
                SField = S + int(addOffset*offset)
                iField.append(SField)
            self.listTrajectories[iB].dict['iField'] = iField
            
        # (b) Find the image with the best std within each n-uplet
            
        bestStd = self.findBestStd()
        for i in range(self.NB):
            self.listTrajectories[i].dict['bestStd'] = bestStd
            
            
            
            
            
     
    def importTrajectories(self, path, iB):
        """
        """
        self.listTrajectories.append(Trajectory(self.I, self.listFrames, self.scale, self.Zstep, iB))
        traj_df = pd.read_csv(path, sep = '\t')
        cols = traj_df.columns.values
        cols_to_remove = []
        for c in cols:
            if 'Unnamed' in c:
                cols_to_remove.append(c)
        traj_df = traj_df.drop(columns = cols_to_remove)
        self.listTrajectories[-1].dict = traj_df.to_dict(orient = 'list')
        for i in range(len(self.listTrajectories[-1].dict['iF'])):
            iBoi =  self.listTrajectories[-1].dict['iB_inFrame'][i]
            iF =  self.listTrajectories[-1].dict['iF'][i]
            self.listTrajectories[-1].dict['Bead'][i] = self.listFrames[iF].listBeads[iBoi]

        
        
        
        
        
    def computeForces(self, traj1, traj2, B0, D3, dx):
        """
        """
    
        # Magnetization functions
        def computeMag_M270(B):
            M = 0.74257*1.05*1600 * (0.001991*B**3 + 17.54*B**2 + 153.4*B) / (B**2 + 35.53*B + 158.1)
            return(M)

        def computeMag_M450(B):
            M = 1.05*1600 * (0.001991*B**3 + 17.54*B**2 + 153.4*B) / (B**2 + 35.53*B + 158.1)
            return(M)

        dictMag = {'M270' : computeMag_M270, 'M450' : computeMag_M450}
        dictBeadTypes = {2.7 : 'M270', 4.5 : 'M450'}
        
        dictLogF = {'D3' : [], 'B0' : [], 'Btot_L' : [], 'Btot_R' : [], 'F00' : [], 'F0' : [], 'dF_L' : [], 'dF_R' : [], 'Ftot' : []}

        # Correction functions
        def Bind_neighbour(B, D_BoI, neighbourType):
            if neighbourType == '' or neighbourType == 'nan':
                return(0)

            else:
                D_neighbour = self.dictBeadDiameters[neighbourType]
                f_Mag = dictMag[neighbourType] # Appropriate magnetization function
                M_neighbour = f_Mag(B) # magnetization [A.m^-1]
                V_neighbour = (4/3)*np.pi*(D_neighbour/2)**3 # volume [nm^3]
                m_neighbour = M_neighbour*V_neighbour*1e-9 # magnetic moment [A.nm^2]

                D_tot = (D_BoI + D_neighbour)/2 # Center-to-center distance [nm]
                B_ind = 2e5*m_neighbour/(D_tot**3) # Inducted mag field [mT]
                return(B_ind)

        def deltaF_neighbour(m_BoI, B, D_BoI, D_BoI2, neighbourType):
            if neighbourType == '' or neighbourType == 'nan':
                return(0)

            else:
                D_neighbour = self.dictBeadDiameters[neighbourType]
                f_Mag = dictMag[neighbourType] # Appropriate magnetization function
                M_neighbour = f_Mag(B) # magnetization [A.m^-1]
                V_neighbour = (4/3)*np.pi*(D_neighbour/2)**3 # volume [nm^3]
                m_neighbour = M_neighbour*V_neighbour*1e-9 # magnetic moment [A.nm^2]

                D_tot = D_BoI/2 + D_BoI2 + D_neighbour/2
                deltaF = 3e5*m_BoI*m_neighbour/D_tot**4 # force [pN]
                return(deltaF)

        # Let's make sure traj1 is the left bead traj and traj2 the right one.
        avgX1 = np.mean(traj1.dict['X'])
        avgX2 = np.mean(traj2.dict['X'])
        if avgX1 < avgX2:
            traj_L, traj_R = traj1, traj2
        else:
            traj_L, traj_R = traj2, traj1

        # Get useful data
        BeadType_L, BeadType_R = dictBeadTypes[traj_L.D], dictBeadTypes[traj_R.D]
        Neighbours_BL = np.concatenate(([traj_L.dict['Neighbour_L']], [traj_L.dict['Neighbour_R']]), axis = 0)
        Neighbours_BR = np.concatenate(([traj_R.dict['Neighbour_L']], [traj_R.dict['Neighbour_R']]), axis = 0)
        D_L, D_R = self.dictBeadDiameters[BeadType_L], self.dictBeadDiameters[BeadType_R]

        nT = len(B0)
        D3nm = 1000*D3
        Dxnm = 1000*dx
        F = np.zeros(nT)

        # Maybe possible to process that faster on lists themselves
        for i in range(nT):
            # Appropriate magnetization functions
            f_Mag_L = dictMag[BeadType_L]
            f_Mag_R = dictMag[BeadType_R]

            # Btot = B0 + B inducted by potential left neighbour mag + B inducted by potential right neighbour mag
            Btot_L = B0[i] + Bind_neighbour(B0[i], D_L, Neighbours_BL[0,i]) + Bind_neighbour(B0[i], D_L, Neighbours_BL[1,i])
            Btot_R = B0[i] + Bind_neighbour(B0[i], D_R, Neighbours_BR[0,i]) + Bind_neighbour(B0[i], D_R, Neighbours_BR[1,i])

            # Magnetizations
            M_L = f_Mag_L(Btot_L)
            M_R = f_Mag_R(Btot_R)

            # Volumes
            V_L = (4/3)*np.pi*(D_L/2)**3 # volume [nm^3]
            V_R = (4/3)*np.pi*(D_R/2)**3 # volume [nm^3]

            # Magnetizations
            m_L = M_L * 1e-9 * V_L
            m_R = M_R * 1e-9 * V_R

            anglefactor = abs(3*(Dxnm[i]/D3nm[i])**2 - 1)

            # Forces
            F00 = 3e5*anglefactor * (f_Mag_L(B0[i])* 1e-9*V_L) * (f_Mag_R(B0[i])*1e-9*V_R) / (D3nm[i]**4)
            F0 = 3e5*anglefactor*m_L*m_R/D3nm[i]**4
            dF_L = deltaF_neighbour(m_L, B0[i], D_L, D_R, Neighbours_BR[1,i])
            dF_R = deltaF_neighbour(m_R, B0[i], D_R, D_L, Neighbours_BL[0,i])

            # Total force = force between beads involved in the pair (F0)
            #               + small force between B_L and B_R's potential right neighbour
            #               + small force between B_R and B_L's potential left neighbour
            F[i] = F0 + dF_L + dF_R
            
            dictLogF['D3'].append(D3nm[i]-(D_L+D_R)/2)
            dictLogF['B0'].append(B0[i])
            dictLogF['Btot_L'].append(Btot_L)
            dictLogF['Btot_R'].append(Btot_R)
            dictLogF['F00'].append(F00)
            dictLogF['F0'].append(F0)
            dictLogF['dF_L'].append(dF_L)
            dictLogF['dF_R'].append(dF_R)
            dictLogF['Ftot'].append(F[i])
        
        dfLogF = pd.DataFrame(dictLogF)
        
        return(F, dfLogF)
        

# %%%% Frame
        
class Frame:
    def __init__(self, F, iS, NB, threshold, Nup, status_frame, status_nUp, scale):
        ny, nx = F.shape[0], F.shape[1]
        self.F = F # Note : Frame.F points directly to the i-th frame of the image I ! To have 2 different versions one should use np.copy(F)
        self.threshold = threshold
        self.NBoi = NB
        self.NBdetected = 0
        self.nx = nx
        self.ny = ny
        self.iS = iS
        self.listBeads = []
        self.trajPoint = []
        self.Nuplet = Nup
        self.status_frame = status_frame
        self.status_nUp = status_nUp
        self.scale = scale
        self.resDf = pd.DataFrame({'Area' : [], 'StdDev' : [], 'XM' : [], 'YM' : [], 'Slice' : []})
 
    def __str__(self):
        text = 'a'
        return(text)
    
    def show(self, strech = True):
        fig, ax = plt.subplots(1,1)
#         fig_size = plt.gcf().get_size_inches()
#         fig.set_size_inches(2 * fig_size)
        if strech:
            pStart, pStop = np.percentile(self.F, (1, 99))
            ax.imshow(self.F, cmap = 'gray', vmin = pStart, vmax = pStop)
        else:
            ax.imshow(self.F, cmap = 'gray')
        if len(self.listBeads) > 0:
            for B in self.listBeads:
                ax.plot([B.x], [B.y], c='orange', marker='+', markersize = 15)
        fig.show()
    
    def detectBeads(self, plot):
        F_bin = self.F > self.threshold
        F_lab, nObj = ndi.label(F_bin)
        props = measure.regionprops(F_lab)
        listValidLabels = []
        areas = np.zeros(nObj+1)
        
        if plot >= 1:
            A_plot = np.zeros(nObj+1)
            P_plot = np.zeros(nObj+1)
            C_plot = np.zeros(nObj+1)
            Valid_Plot = np.zeros(nObj+1)
            Plot_Labels = []
        
        for k in range(1, nObj+1):
#             try:
            bb = props[k-1].bbox
            Valid = not (min(bb) == 0 or bb[2] == self.ny or bb[3] == self.nx) # Remove objects touching the edges of the frame

#             # OPTION 1 - Compute the metrics on the filled shape ; NB: takes a lot of time
#             F_fh = ndi.binary_fill_holes((F_lab == k).astype(int))
#             tmp_props = measure.regionprops(F_fh.astype(int))
#             A = tmp_props[0].area
#             P = tmp_props[0].perimeter
#             Circ = (4 * np.pi * A) / (P * P)
#             Valid = Valid and A >= 100 and Circ >= 0.75 # Area and circularity criterion

            # OPTION 2 - Lower the criterion in circularity - NB: less selective
            A = props[k-1].area
            P = props[k-1].perimeter
            Circ = (4 * np.pi * A) / (P * P)
            if A >= 75 and Circ >= 0.7 and A < 1200:
                pass
            elif A >= 300 and A < 1200 and Circ < 0.7 and Circ > 0.3:
                F_fh = ndi.binary_fill_holes((F_lab == k).astype(int))
                tmp_props = measure.regionprops(F_fh.astype(int))
                A = tmp_props[0].area
                P = tmp_props[0].perimeter
                Circ = (4 * np.pi * A) / (P * P)

            Valid = (Valid and A >= 75 and A < 1200 and Circ >= 0.65) # Area and circularity criterion

#             # OPTION 3 - 
#             A = props[k-1].area
#             P = props[k-1].perimeter
#             Circ = (4 * np.pi * A) / (P * P)
#             if A >= 80:µ
#                 ObjFilled = props[k-1].image_filled
#                 propsObjFilled = measure.regionprops(ObjFilled)
#                 A = propsObjFilled.area
#                 P = propsObjFilled.perimeter
#                 Circ = (4 * np.pi * A) / (P * P)

#             Valid = (Valid and A >= 100 and Circ >= 0.65)

            if plot >= 1:
                A_plot[k] = A
                C_plot[k] = Circ
                P_plot[k] = P
                Valid_Plot[k] = Valid
                Plot_Labels.append(k)
                

#             except:
#                 Valid = False
                
            if Valid:
                listValidLabels.append(k)
                areas[k] = A
                
                
#         F_labValid, nObjValid = ndi.label(F_lab)    
#         fig, ax = plt.subplots(1,1)
#         ax.imshow(F_labValid)
#         fig.show()

        if plot >= 1:
            centroids_plot = np.array([ndi.center_of_mass(self.F, labels=F_lab, index=m) for m in Plot_Labels])
            color = np.where(Valid_Plot == True, 'green', 'red')
            fig, ax = plt.subplots(1,1)
            ax.imshow(F_lab)
            for m in Plot_Labels:
                ax.plot(centroids_plot[m-1,1], centroids_plot[m-1,0], marker = '+', color = color[m])
                S = 'A = {:.2f} ; P = {:.2f} \nC = {:.2f}'.format(A_plot[m], P_plot[m], C_plot[m])
                ax.text(centroids_plot[m-1,1], centroids_plot[m-1,0], S, color = color[m])
            fig.suptitle(str(self.iS))
            fig.show()
        
        centroids = np.array([ndi.center_of_mass(self.F, labels=F_lab, index=i) for i in listValidLabels])
        try:
            if centroids.ndim == 2:
                resDict = {}
                resDict['Area'] = np.array([areas[i] for i in listValidLabels]).astype(int)
                resDict['StdDev'] = np.array([ndi.standard_deviation(self.F, labels=F_lab, index=i) for i in listValidLabels])
                resDict['XM'] = centroids[:,1]
                resDict['YM'] = centroids[:,0]
                resDict['Slice'] = np.array([self.iS+1 for i in listValidLabels]).astype(int)
            elif centroids.ndim == 1 and len(centroids) >= 2:
                resDict = {}
                resDict['Area'] = np.array([areas[i] for i in listValidLabels]).astype(int)
                resDict['StdDev'] = np.array([ndi.standard_deviation(self.F, labels=F_lab, index=i) for i in listValidLabels])
                resDict['XM'] = centroids[1]
                resDict['YM'] = centroids[0]
                resDict['Slice'] = np.array([self.iS+1 for i in listValidLabels]).astype(int)
            else:
                resDict = {}
        except:
            print(centroids)
        
        self.resDf = pd.DataFrame(resDict)
        
#         print(self.resDict)

    def makeListBeads(self):
        self.NBdetected = self.resDf.shape[0]
        for i in range(self.NBdetected):
            d = {}
            for c in self.resDf.columns:
                d[c] = self.resDf[c].values[i]
            self.listBeads.append(Bead(d, self.F))
            
    def beadsXYarray(self):
        A = np.zeros((len(self.listBeads), 2))
        for i in range(len(self.listBeads)):
            b = self.listBeads[i]
            A[i,0], A[i,1] = b.x, b.y
        return(A)
    
    def beadsStdDevarray(self):
        A = np.zeros(len(self.listBeads))
        for i in range(len(self.listBeads)):
            b = self.listBeads[i]
            A[i] = b.std
        return(A)
    
    def detectDiameter(self, plot = 0):
        for i in range(len(self.listBeads)):
            B = self.listBeads[i]
            B.detectDiameter(self.nx, self.ny, self.scale, plot)
#             B = self.listBeads[i]
#             roughSize = np.floor(5*SCALE_100X)
#             roughSize += roughSize%2
#             x1_roughROI = max(int(np.floor(B.x) - roughSize*0.5) - 1, 0)
#             x2_roughROI = min(int(np.floor(B.x) + roughSize*0.5), self.nx)
#             y1_roughROI = max(int(np.floor(B.y) - roughSize*0.5) - 1, 0)
#             y2_roughROI = min(int(np.floor(B.y) + roughSize*0.5), self.ny)        
#             F_roughRoi = self.F[y1_roughROI:y2_roughROI, x1_roughROI:x2_roughROI]            
#             if plot == 1:
#                 figh, axh = plt.subplots(1,1, figsize = (4,4))
#                 axh.hist(F_roughRoi.ravel(), bins=256, histtype='step', color='black')
#             counts, binEdges = np.histogram(F_roughRoi.ravel(), bins=256)
#             peaks, peaksProp = find_peaks(counts, height=100, threshold=None, distance=None, prominence=None, \
#                                width=None, wlen=None, rel_height=0.5, plateau_size=None)       
#             peakThreshVal = 750    
#             if counts[peaks[0]] > peakThreshVal:
#                 B.D = 4.5      
#             else:
#                 B.D = 2.7   
#             if plot == 1:    
#                 print(B.D)
#                 B.show()
            pass
                
    def detectNeighbours(self, iB, beadType = 'detect'):
        BoI = self.listBeads[iB]
        listD = []
        for B in self.listBeads:
            if B.D not in listD:
                listD.append(B.D)
        Dmax = max(listD)*self.scale
        DBoI = BoI.D*self.scale
        x, y = BoI.x, BoI.y
#         print('Bead in (x,y) = (' + str(x) + ', ' + str(y) + ')')
        
        def makeRoI(x01, x02, y01, y02, nx, ny):
            x1, x2, y1, y2 = int(x01), int(x02), int(y01), int(y02)
            Lx = [x1, x2]
            Ly = [y1, y2]
            for x in Lx:
                if x < 0:
                    x = 0
                if x > nx:
                    x = nx
            for y in Ly:
                if y < 0:
                    y = 0
                if y > ny:
                    y = ny
            return([Lx, Ly])
        
        def areThereBeadsInRoI(listBeads, RoI):
            beadsFound = []
            Lx, Ly = RoI[0], RoI[1]
#             print(RoI)
            for B in listBeads:
#                 print(B.x, B.y)
                if B.x > Lx[0] and B.x < Lx[1] and B.y > Ly[0] and B.y < Ly[1]:
                    beadsFound.append(B)
            return(beadsFound)
        
        RoI_left = makeRoI(x-(0.5*DBoI+Dmax), x-(0.5*DBoI), y-0.55*Dmax, y+0.55*Dmax, self.nx, self.ny)
        RoI_right = makeRoI(x+(0.5*DBoI), x+(0.5*DBoI+Dmax), y-0.55*Dmax, y+0.55*Dmax, self.nx, self.ny)
        
        beadsFound_left = areThereBeadsInRoI(self.listBeads, RoI_left)
        beadsFound_right = areThereBeadsInRoI(self.listBeads, RoI_right)
        
        if len(beadsFound_left) >= 1:
            if beadType == 'detect':
                beadsFound_left[0].detectDiameter(self.nx, self.ny, self.scale, plot = 0)
                D = beadsFound_left[0].D
                if D == 4.5:
                    BoI.Neighbour_L = 'M450'
                elif D == 2.7:
                    BoI.Neighbour_L = 'M270'
            else:
                BoI.Neighbour_L = beadType
                
        if len(beadsFound_right) >= 1:
            if beadType == 'detect':
                beadsFound_right[0].detectDiameter(self.nx, self.ny, self.scale, plot = 0)
                D = beadsFound_right[0].D
                if D == 4.5:
                    BoI.Neighbour_R = 'M450'
                elif D == 2.7:
                    BoI.Neighbour_R = 'M270'
            else:
                BoI.Neighbour_R = beadType
        
        return(BoI.Neighbour_L, BoI.Neighbour_R)

    #
        
# %%%% Bead
        
class Bead:
    def __init__(self, d, F):
        self.x = d['XM']
        self.y = d['YM']
        self.D = 0
        self.area = d['Area']
        self.std = d['StdDev']
        self.iS = d['Slice']-1
        self.status_frame = ''
        self.Neighbour_L = ''
        self.Neighbour_R = ''
        self.F = F
        
    def detectDiameter(self, nx, ny, scale, plot = 0):
        F = self.F
        roughSize = np.floor(5*scale)
        roughSize += roughSize%2
        x1_roughROI = max(int(np.floor(self.x) - roughSize*0.5) - 1, 0)
        x2_roughROI = min(int(np.floor(self.x) + roughSize*0.5), nx)
        y1_roughROI = max(int(np.floor(self.y) - roughSize*0.5) - 1, 0)
        y2_roughROI = min(int(np.floor(self.y) + roughSize*0.5), ny)

        F_roughRoi = F[y1_roughROI:y2_roughROI, x1_roughROI:x2_roughROI]   

        if plot == 1:
            figh, axh = plt.subplots(1,1, figsize = (4,4))
            axh.hist(F_roughRoi.ravel(), bins=256, histtype='step', color='black')

        counts, binEdges = np.histogram(F_roughRoi.ravel(), bins=256)

        peaks, peaksProp = find_peaks(counts, height=100, threshold=None, distance=None, prominence=None, \
                           width=None, wlen=None, rel_height=0.5, plateau_size=None)

        peakThreshVal = 750

        if counts[peaks[0]] > peakThreshVal:
            self.D = 4.5

        else:
            self.D = 2.7

        if plot == 1:    
            print(self.D)
            self.show()

    def show(self, strech = True):
        fig, ax = plt.subplots(1,1)
        if strech:
            pStart, pStop = np.percentile(self.F, (1, 99))
            ax.imshow(self.F, cmap = 'gray', vmin = pStart, vmax = pStop)
        else:
            ax.imshow(self.F, cmap = 'gray')
        ax.plot([self.x], [self.y], c='orange', marker='o')
        fig.show()

#
        
# %%%% Trajectory
        
class Trajectory:
    def __init__(self, I, listFrames, scale, Zstep, iB):
        nS, ny, nx = I.shape[0], I.shape[1], I.shape[2]
        self.I = I
        self.listFrames = listFrames
        self.scale = scale
        self.nx = nx
        self.ny = ny
        self.nS = nS
        self.D = 0
        self.nT = 0
        self.iB = iB
        self.dict = {'X': [],'Y': [],'idxAnalysis': [],'StdDev': [], 
                     'Bead': [],'status_frame': [],'status_nUp': [],'iF': [],'iS': [],'iB_inFrame' : [],
                     'bestStd' : [], 'Zr' : [], 'Neighbour_L' : [], 'Neighbour_R' : []}
        # iF is the index in the listFrames
        # iS is the index of the slice in the raw image MINUS ONE
        self.beadInOut = ''
        self.deptho = []
        self.depthoPath = ''
        self.depthoStep = 20
        self.depthoZFocus = 200
        self.Zstep = Zstep # The step in microns between 2 consecutive frames in a multi-frame Nuplet
                
    def __str__(self):
        text = 'iS : ' + str(self.series_iS)
        text += '\n'
        text += 'XY : ' + str(self.seriesXY)
        return(text)
    
    def save(self, path):
        df = pd.DataFrame(self.dict)
        df.to_csv(path, sep = '\t', index = False)
    
    def computeZ(self, matchingDirection, plot = 0):
        if len(self.deptho) == 0:
            return('Error, no depthograph associated with this trajectory')
        
        else:
            Ddz, Ddx = self.deptho.shape[0], self.deptho.shape[1]
            iF = self.dict['iF'][0]
            previousZ = -1
            while iF <= max(self.dict['iF']):

#### Important plotting option here
# ####### Decomment these lines to enable some plots ##################
                
                # plot = 0
                # if iF >= 0 and iF <= 30:# or (iF < 190 and iF > 150):
                #     plot = 1
                
# ############################ OK?! ###################################
                
                if iF not in self.dict['iF']: # this index isn't in the trajectory list => the frame was removed for some reason.
                    iF += 1 # Let's just go to the next index
                    
                else:
                    F = self.listFrames[iF]
                    Nup = F.Nuplet
                    if Nup <= 1:
                        framesNuplet = [F]
                        iFNuplet = [iF]
                        iF += 1
                    elif Nup > 1:
                        framesNuplet = [F]
                        iFNuplet = [iF]
                        jF = 1
                        while iF+jF <= max(self.dict['iF']) and self.listFrames[iF+jF].status_nUp == F.status_nUp:
                            if iF+jF in self.dict['iF']: # One of the images of the triplet may be invalid, 
                                # and we don't want to take it. With this test we won't
                                nextF = self.listFrames[iF+jF]
                                framesNuplet.append(nextF)
                                iFNuplet.append(iF+jF)
                            jF += 1
                            
                        iF += jF
                    
                    maxDz = 40
                    Z = self.findZ_Nuplet(framesNuplet, iFNuplet, Nup, previousZ, matchingDirection, maxDz, plot)
                    previousZ = Z
                    # This Z_pix has no meaning in itself, it needs to be compared to the depthograph Z reference point,
                    # which is depthoZFocus. 
                    Zr = self.depthoZFocus - Z # If you want to find it back, Z = depthoZFocus - Zr
                    # This definition was chosen so that when Zr > 0, the plane of observation of the bead is HIGHER than the focus
                    # and accordingly when Zr < 0, the plane of observation of the bead is LOWER than the focus

                    mask = np.array([(iF in iFNuplet) for iF in self.dict['iF']])
                    self.dict['Zr'][mask] = Zr
        
                     
    def findZ_Nuplet(self, framesNuplet, iFNuplet, Nup, previousZ, matchingDirection, maxDz = 40, plot = 0):
        try:
            Nframes = len(framesNuplet)
            listStatus_1 = [F.status_frame for F in framesNuplet]
            listXY = [[self.dict['X'][np.where(self.dict['iF']==iF)][0], 
                       self.dict['Y'][np.where(self.dict['iF']==iF)][0]] for iF in iFNuplet]
            listiS = [self.dict['iS'][np.where(self.dict['iF']==iF)][0] for iF in iFNuplet]        
            cleanSize = getDepthoCleanSize(self.D, self.scale)
            hdSize = self.deptho.shape[1]
            depthoDepth = self.deptho.shape[0]
            listProfiles = np.zeros((Nframes, hdSize))
            listROI = []
            for i in range(Nframes):
                xx = np.arange(0, 5)
                yy = np.arange(0, cleanSize)
                try:
                    X, Y, iS = int(np.round(listXY[i][0])), int(np.round(listXY[i][1])), listiS[i] # > We could also try to recenter the image to keep a subpixel resolution here
                    # line that is 5 pixels wide
                    profileROI = framesNuplet[i].F[Y-cleanSize//2:Y+cleanSize//2+1, X-2:X+3]
                    f = interpolate.interp2d(xx, yy, profileROI, kind='cubic')
                    # Now use the obtained interpolation function and plot the result:
                    xxnew = xx
                    yynew = np.linspace(0, cleanSize, hdSize)
                    profileROI_hd = f(xxnew, yynew)
                
                except: # If the vertical slice doesn't work, try the horizontal one
                    print(ORANGE + 'error with the vertical slice -> trying with horizontal one')
                    print('iFNuplet')
                    print(iFNuplet)
                    print('Roi')
                    print(Y-2,Y+3, X-cleanSize//2,X+cleanSize//2+1)
                    print('' + NORMAL)
                    
                    xx, yy = yy, xx
                    X, Y, iS = int(np.round(listXY[i][0])), int(np.round(listXY[i][1])), listiS[i] # > We could also try to recenter the image to keep a subpixel resolution here
                    # line that is 5 pixels wide
                    profileROI = framesNuplet[i].F[Y-2:Y+3, X-cleanSize//2:X+cleanSize//2+1]
                    f = interpolate.interp2d(xx, yy, profileROI, kind='cubic')
                    # Now use the obtained interpolation function and plot the result:
                    xxnew = np.linspace(0, cleanSize, hdSize)
                    yynew = yy
                    profileROI_hd = f(xxnew, yynew).T
    
                listROI.append(profileROI)
    
                listProfiles[i,:] = profileROI_hd[:,5//2] * (1/5)
                for j in range(1, 1 + 5//2):
                    listProfiles[i,:] += profileROI_hd[:,5//2-j] * (1/5)
                    listProfiles[i,:] += profileROI_hd[:,5//2+j] * (1/5)       
    
            listProfiles = listProfiles.astype(np.uint16)
    
            # now use listStatus_1, listProfiles, self.deptho + data about the jump between Nuplets ! (TBA) 
            # to compute the correlation function
            nVoxels = int(np.round(self.Zstep/self.depthoStep))
    
            listDistances = np.zeros((Nframes, depthoDepth))
            listZ = np.zeros(Nframes, dtype = int)
            for i in range(Nframes):
                listDistances[i] = squareDistance(self.deptho, listProfiles[i], normalize = True) # Utility functions
                listZ[i] = np.argmin(listDistances[i])
                
            # Translate the profiles that must be translated (status_frame 1 & 3 if Nup = 3)
            # and don't move the others (status_frame 2 if Nup = 3 or the 1 profile when Nup = 1)
            if Nup > 1:
                finalDists = matchDists(listDistances, listStatus_1, Nup, nVoxels, direction = matchingDirection)
            elif Nup == 1:
                finalDists = listDistances
                
            sumFinalD = np.sum(finalDists, axis = 0)
            if previousZ == -1 or Nup > 1: # First image OR triplets => No restriction
                Z = np.argmin(sumFinalD)
            elif Nup == 1 and previousZ != -1: # Not first image AND singlet => Restriction
                limInf = max(previousZ-maxDz, 0)
                limSup = min(previousZ+maxDz, depthoDepth)
                Z = limInf + np.argmin(sumFinalD[limInf:limSup])
    
            #### Important plotting option here
            if plot >= 1:
                fig, axes = plt.subplots(5, 3, figsize = (20,10))
                fig.tight_layout()
                im = framesNuplet[0].F
                X2, Y2 = listXY[0][0], listXY[0][1]
                
                pStart, pStop = np.percentile(im, (1, 99))
                axes[0,0].imshow(im, vmin = pStart, vmax = 1.5*pStop, cmap = 'gray')
                col_im = 'cyan'
                dx, dy = 10, 60
                axes[0,0].plot([X2], [Y2], marker = '+', color = 'red')
                axes[0,0].plot([X2-dx,X2-dx], [Y2-dy,Y2+dy], ls = '--', color = col_im)
                axes[0,0].plot([X2+dx,X2+dx], [Y2-dy,Y2+dy], ls = '--', color = col_im)
                axes[0,0].plot([X2-dx,X2+dx], [Y2-dy,Y2-dy], ls = '--', color = col_im)
                axes[0,0].plot([X2-dx,X2+dx], [Y2+dy,Y2+dy], ls = '--', color = col_im)
                
                axes[0,1].imshow(self.deptho)
                # TEST !!! # -> Works well !
                XL0, YL0 = axes[0,1].get_xlim(), axes[0,1].get_ylim()
                extent = (XL0[0], YL0[0]*(5/3), YL0[0], YL0[1])
                axes[0,1].imshow(self.deptho, extent = extent)
                # TEST !!! #
    
                pixLineHD = np.arange(0, hdSize, 1)
                zPos = np.arange(0, depthoDepth, 1)
                col = ['orange', 'gold', 'green']
                for i in range(Nframes):
                    axes[1,i].imshow(listROI[i])
                    #
                    axes[2,i].plot(pixLineHD, listProfiles[i])
                    #
                    axes[3,i].plot(zPos, listDistances[i])
                    limy3 = axes[3,i].get_ylim()
                    min_i = np.argmin(listDistances[i])
                    axes[3,i].plot([min_i,min_i],limy3,ls = '--', c = col[i])
                    #
                    axes[4,i].plot(zPos, finalDists[i])
                    limy4 = axes[4,i].get_ylim()
                    min_i = np.argmin(finalDists[i])
                    axes[4,i].plot([min_i,min_i],limy4,ls = '--', c = col[i])
                    
                    axes[0,1].plot([axes[0,1].get_xlim()[0], axes[0,1].get_xlim()[1]-1], [listZ[i], listZ[i]], ls = '--', c = col[i])
                    axes[0,1].plot([axes[0,1].get_xlim()[0], axes[0,1].get_xlim()[1]-1], [Z,Z], ls = '--', c = 'red')

                    
                axes[0,2].plot(zPos, sumFinalD)
                limy0 = axes[0,2].get_ylim()
                axes[0,2].plot([Z,Z],limy0,ls = '-', c = 'red')
                axes[0,2].plot([previousZ,previousZ],limy0,ls = '--', c = 'pink')
                axes[0,2].plot([previousZ-maxDz,previousZ-maxDz],limy0,ls = '--', c = 'cyan')
                axes[0,2].plot([previousZ+maxDz,previousZ+maxDz],limy0,ls = '--', c = 'cyan')
                
                iSNuplet = [F.iS+1 for F in framesNuplet]
                fig.suptitle('Frames ' + str(iFNuplet) + ' - Slices ' + str(iSNuplet) + ' ; Z = ' + str(Z))
                Nfig = plt.gcf().number
                fig.savefig('C://Users//JosephVermeil//Desktop//TempPlot//fig'+str(Nfig)+'.png')
    
            return(Z)
        
        except Exception:
            print(RED + '')
            traceback.print_exc()
            print('\n')
            print(ORANGE + 'Error with the Z detection')
            print('iFNuplet')
            print(iFNuplet)
            print('Roi')
            print(Y-2,Y+3, X-cleanSize//2,X+cleanSize//2+1)
            print('Deptho shape')
            print(self.deptho.shape)
            print('Shapes of listDistances, finalDists, sumFinalD')
            print(listDistances.shape)
            print(finalDists.shape)
            print(sumFinalD.shape)
            print('previousZ, previousZ-maxDz, previousZ+maxDz')
            print(previousZ, previousZ-maxDz, previousZ+maxDz)
            print('' + NORMAL)
            
    
    def keepBestStdOnly(self):
        dictBestStd = {}
        bestStd = self.dict['bestStd']
        nT = int(np.sum(bestStd))
        for k in self.dict.keys():
            A = np.array(self.dict[k])
            dictBestStd[k] = A[bestStd]
        self.dict = dictBestStd
        self.nT = nT
        
    def detectNeighbours(self, frequency, beadType = 'detect'):
        init = True
        listNeighbours = []
        for i in range(len(self.dict['iF'])):
            iS = self.dict['iS'][i]
            if frequency%iS == 0 or init:
                iF = self.dict['iF'][i]
                iB = self.dict['iB_inFrame'][i]
                Frame = self.listFrames[iF]
                Neighbour_L, Neighbour_R = Frame.detectNeighbours(iB, beadType)
                listNeighbours.append([Neighbour_L, Neighbour_R])
                init = False
            else:
                self.dict['Bead'][i].Neighbour_L = Neighbour_L
                self.dict['Bead'][i].Neighbour_R = Neighbour_R
                listNeighbours.append([Neighbour_L, Neighbour_R])
        arrayNeighbours = np.array(listNeighbours)
        self.dict['Neighbour_L'] = arrayNeighbours[:,0]
        self.dict['Neighbour_R'] = arrayNeighbours[:,1]
#         print(self.dict['Neighbour_L'], self.dict['Neighbour_R'])


        
    def detectNeighbours_ui_V0(self, Nimg, frequency, beadType): # NOT VERY WELL MADE FOR NOW
        # Plots to help the user to see the neighbour of each bead
        ncols = 4
        nrows = ((Nimg-1) // ncols) + 1

        fig, ax = plt.subplots(nrows, ncols)
        for i in range(Nimg):
            pos = np.searchsorted(self.dict['iS'], i*frequency, 'left')
            iS = self.dict['iS'][pos]
            iF = self.dict['iF'][pos]
            ax[i//ncols,i%ncols].imshow(self.I[iS], cmap = 'gray')
            ax[i//ncols,i%ncols].set_title('Loop ' + str(i+1))
            ax[i//ncols,i%ncols].plot([self.dict['X'][pos]],[self.dict['Y'][pos]], 'ro')
        
        # Ask the question
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(720, 50, 1175, 1000)
        QA = pyautogui.confirm(
            text='Neighbours of the selected bead?',
            title='', 
            buttons=['Left','Right', 'Left and right'])
        
        # According to the question's answer:
        if QA == 'Left':
            Neighbour_L = beadType
            Neighbour_R = ''
        elif QA == 'Right':
            Neighbour_L = ''
            Neighbour_R = beadType
        elif QA == 'Left and right':
            Neighbour_L, Neighbour_R = beadType, beadType
        
        plt.close(fig)
        listNeighbours = []
        
        for i in range(len(self.dict['iF'])):
            self.dict['Bead'][i].Neighbour_L = Neighbour_L
            self.dict['Bead'][i].Neighbour_R = Neighbour_R
            listNeighbours.append([Neighbour_L, Neighbour_R])
            
        arrayNeighbours = np.array(listNeighbours)
        self.dict['Neighbour_L'] = arrayNeighbours[:,0]
        self.dict['Neighbour_R'] = arrayNeighbours[:,1]
        
    def detectNeighbours_ui(self, Nimg, frequency, beadType): # NOT VERY WELL MADE FOR NOW
        # Plots to help the user to see the neighbour of each bead
        ncols = 4
        nrows = ((Nimg-1) // ncols) + 1

        fig, ax = plt.subplots(nrows, ncols)
        for i in range(Nimg):
            try:
                pos = np.searchsorted(self.dict['iS'], i*frequency, 'left')
                iS = self.dict['iS'][pos]
                iF = self.dict['iF'][pos]
                pStart, pStop = np.percentile(self.I[iS], (1, 99))
                if nrows > 1:
                    ax[i//ncols,i%ncols].imshow(self.I[iS], cmap = 'gray', vmin = pStart, vmax = pStop)
                    ax[i//ncols,i%ncols].set_title('Loop ' + str(i+1))
                    ax[i//ncols,i%ncols].plot([self.dict['X'][pos]],[self.dict['Y'][pos]], 'ro')
                elif nrows == 1:
                    ax[i].imshow(self.I[iS], cmap = 'gray', vmin = pStart, vmax = pStop)
                    ax[i].set_title('Loop ' + str(i+1))
                    ax[i].plot([self.dict['X'][pos]],[self.dict['Y'][pos]], 'ro')
            except:
                print(RED  + 'ptit probleme dans le detectNeighbours_ui' + NORMAL)
        
        # Ask the question
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(720, 50, 1175, 1000)
        QA = pyautogui.confirm(
            text='Neighbours of the selected bead?',
            title='', 
            buttons=['1', '2'])
        
        
        # According to the question's answer:
        if QA == '1':
            if self.iB%2 == 0: # the bead is on the left of a pair
                Neighbour_L, Neighbour_R = '', beadType
            elif self.iB%2 == 1: # the bead is on the right of a pair
                Neighbour_L, Neighbour_R = beadType, ''
        elif QA == '2':
            Neighbour_L, Neighbour_R = beadType, beadType
        
        plt.close(fig)
        listNeighbours = []
        
        for i in range(len(self.dict['iF'])):
            self.dict['Bead'][i].Neighbour_L = Neighbour_L
            self.dict['Bead'][i].Neighbour_R = Neighbour_R
            listNeighbours.append([Neighbour_L, Neighbour_R])
            
        arrayNeighbours = np.array(listNeighbours)
        self.dict['Neighbour_L'] = arrayNeighbours[:,0]
        self.dict['Neighbour_R'] = arrayNeighbours[:,1]



    def detectInOut_ui(self, Nimg, frequency): # NOT VERY WELL MADE FOR NOW
        # Almost copy paste of detectNeighbours_ui
        ncols = 4
        nrows = ((Nimg-1) // ncols) + 1

        fig, ax = plt.subplots(nrows, ncols)
        for i in range(Nimg):
            try:
                pos = np.searchsorted(self.dict['iS'], i*frequency, 'left')
                iS = self.dict['iS'][pos]
                iF = self.dict['iF'][pos]
                pStart, pStop = np.percentile(self.I[iS], (1, 99))
                if nrows > 1:
                    ax[i//ncols,i%ncols].imshow(self.I[iS], cmap = 'gray', vmin = pStart, vmax = pStop)
                    ax[i//ncols,i%ncols].set_title('Loop ' + str(i+1))
                    ax[i//ncols,i%ncols].plot([self.dict['X'][pos]],[self.dict['Y'][pos]], 'ro')
                elif nrows == 1:
                    ax[i].imshow(self.I[iS], cmap = 'gray', vmin = pStart, vmax = pStop)
                    ax[i].set_title('Loop ' + str(i+1))
                    ax[i].plot([self.dict['X'][pos]],[self.dict['Y'][pos]], 'ro')
            except:
                print(RED  + 'ptit probleme dans le detectInOut_ui' + NORMAL)
        
        # Ask the question
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(720, 50, 1175, 1000)
        QA = pyautogui.confirm(
            text='Is it an inside or outside bead?',
            title='', 
            buttons=['In', 'Out'])
        
        self.beadInOut = QA
        plt.close(fig)
        return(QA)
        

            
    def plot(self, ax, i_color):
        colors = ['cyan', 'red', 'blue', 'orange']
        c = colors[i_color]
        ax.plot(self.dict['X'], self.dict['Y'], color=c, lw=0.5)





# %%%% Main
        
def mainTracker(mainDataDir, rawDataDir, depthoDir, interDataDir, figureDir, timeSeriesDataDir,
         dates, manips, wells, cells, depthoNames, expDf, 
         methodT, factorT, redoAllSteps = False, MatlabStyle = False):
    
    start = time.time()
    
    #### 0. Load different data sources & Preprocess : fluo, black images, sort slices (ct/ramp ; down/middle/up)
        #### 0.1 - Make list of files to analyse
    imagesToAnalyse = []
    imagesToAnalyse_Paths = []
    if not isinstance(dates, str):
        rawDirList = [os.path.join(rawDataDir, d) for d in dates]
    else:
        rawDirList = [os.path.join(rawDataDir, dates)]
    for rd in rawDirList:
        fileList = os.listdir(rd)
        for f in fileList:
            if isFileOfInterest(f, manips, wells, cells): # See Utility Functions > isFileOfInterest
                fPath = os.path.join(rd, f)
                if os.path.isfile(fPath[:-4] + '_Field.txt'):
                    imagesToAnalyse.append(f)
                    imagesToAnalyse_Paths.append(os.path.join(rd, f))    

        #### 0.2 - Begining of the Main Loop
    for i in range(len(imagesToAnalyse)): 
        f, fP = imagesToAnalyse[i], imagesToAnalyse_Paths[i]
        manipID = findInfosInFileName(f, 'manipID') # See Utility Functions > findInfosInFileName
        cellID = findInfosInFileName(f, 'cellID') # See Utility Functions > findInfosInFileName
        
        print('\n')
        print(BLUE + 'Analysis of file {:.0f}/{:.0f} : {}'.format(i+1, len(imagesToAnalyse), f))
        print('Loading image and experimental data...' + NORMAL)
        
        #### 0.3 - Load exp data
        if manipID not in expDf['manipID'].values:
            print(RED + 'Error! No experimental data found for: ' + manipID + NORMAL)
            bug
        else:
            expDf_line = expDf.loc[expDf['manipID'] == manipID]
            manipDict = {}
            for c in expDf_line.columns.values:
                manipDict[c] = expDf_line[c].values[0]
        
#         # Extract some data from exp data
#         StrBeadTypes = str(manipDict['bead type'])
#         if '_' in StrBeadTypes:
#             beadTypes = [bT for bT in StrBeadTypes.split('_')]
#         else:
#             beadTypes = [StrBeadTypes]
        
#         StrBeadDiameters = str(manipDict['bead diameter'])
#         if '_' in StrBeadDiameters:
#             beadDiameters = [int(bD) for bD in StrBeadDiameters.split('_')]
#         else:
#             beadDiameters = [int(StrBeadDiameters)]
            
#         dictBeadDiameters = {}
#         for k in range(len(beadTypes)):
#             dictBeadDiameters[beadTypes[k]] = beadDiameters[k]
    
        #### 0.4 - Load image and init PTL
        I = io.imread(fP) # Approx 0.5s per image
        PTL = PincherTimeLapse(I, cellID, manipDict, NB = 2)
    
        #### 0.5 - Load field file
        fieldFilePath = fP[:-4] + '_Field.txt'
        fieldCols = ['B_set', 'T_abs', 'B', 'Z']
        fieldDf = pd.read_csv(fieldFilePath, sep = '\t', names = fieldCols)
        
        #### 0.6 - Check if a log file exists and load it if required
        logFilePath = fP[:-4] + '_LogPY.txt'
        logFileImported = False
        if redoAllSteps:
            pass
        elif os.path.isfile(logFilePath):
            PTL.importLog(logFilePath)
            PTL.dictLog['UILog'] = PTL.dictLog['UILog'].astype(str)
            logFileImported = True
            
        print(BLUE + 'OK!')
        
        print(BLUE + 'Pretreating the image...' + NORMAL)
        
        #### 0.7 - Detect fluo & black images
        current_date = findInfosInFileName(f, 'date')
        current_date = current_date.replace("-", ".")
        fluoDirPath = os.path.join(rawDataDir, current_date + '_Fluo', f[:-4])
        
        PTL.checkIfBlackFrames()
        PTL.saveFluoAside(fluoDirPath, f)

        
        #### 0.8 - Sort slices
        #### ! Exp type dependance here !
        if not logFileImported:
            if 'R40' in f or 'thickness' in f:
                PTL.determineFramesStatus_R40()
            elif 'L40' in f:
                PTL.determineFramesStatus_L40()
            elif 'optoGen' in f:
                pass
                
        PTL.saveLog(display = False, save = (not logFileImported), path = logFilePath)
        
        #### 0.9 - Import or determine global threshold
        MDpath = fP[:-4] + '_MetaDataPY.txt'
        if MatlabStyle:
            PTL.computeThreshold(method = methodT, factorT = factorT)
        elif redoAllSteps:
            PTL.uiThresholding(method = methodT, factorT = factorT)
        else:
            try:
                PTL.threshold = PTL.readMetaData(MDpath, 'threshold')
            except:
                PTL.uiThresholding(method = methodT, factorT = factorT) # Approx 3s per image


        #### 0.10 - Save some metadata
        PTL.saveMetaData(MDpath)
        
        #### 0.11 - Create list of Frame objects
        PTL.makeFramesList()
        
        print(BLUE + 'OK!')
        
        
    #### 1. Detect beads
    
        print(BLUE + 'Detecting all the bead objects...' + NORMAL)
        Td = time.time()
        
        #### 1.1 - Check if a _Results.txt exists and import it if it's the case
        resFilePath = fP[:-4] + '_ResultsPY.txt'
        resFileImported = False
        if redoAllSteps:
            pass
        elif os.path.isfile(resFilePath):
            PTL.importBeadsDetectResult(resFilePath)
            resFileImported = True
        
        if MatlabStyle:
            try:
                resFilePath = fP[:-4] + '_ResultsPY.txt'
                PTL.importBeadsDetectResult(resFilePath)
                resFileImported = True
            except:
                resFilePath = fP[:-4] + '_Results.txt'
                PTL.importBeadsDetectResult(resFilePath)
                resFileImported = True
        
        #### 1.2 - Detect the beads
        # Detect the beads and create the BeadsDetectResult dataframe [if no file has been loaded before] 
        # OR input the results in each Frame objects [if the results have been loaded at the previous step]
        PTL.detectBeads(resFileImported, display = 0)
        
        #### 1.3 - Save the new results if necessary
        if not resFileImported:
            PTL.saveBeadsDetectResult(path=resFilePath)
            
        print(BLUE + 'OK! dT = {:.3f}'.format(time.time()-Td) + NORMAL)


    #### 2. Make trajectories for beads of interest
        # One of the main steps ! The tracking of the beads happens here !
    
        print(BLUE + 'Tracking the beads of interest...' + NORMAL)
        Tt = time.time()
        
        #### 2.1 - Check if some trajectories exist already
        trajDirRaw = os.path.join(timeSeriesDataDir, 'Trajectories_raw')
        trajFilesExist_global = False
        trajFilesImported = False
        
        if redoAllSteps:
            pass
        else:
            allTrajPaths = [os.path.join(trajDirRaw, f[:-4] + '_rawTraj' + str(iB) + '' + '_PY.csv') for iB in range(PTL.NB)]
            allTrajPaths += [os.path.join(trajDirRaw, f[:-4] + '_rawTraj' + str(iB) + '_In' + '_PY.csv') for iB in range(PTL.NB)]
            allTrajPaths += [os.path.join(trajDirRaw, f[:-4] + '_rawTraj' + str(iB) + '_Out' + '_PY.csv') for iB in range(PTL.NB)]
            allTrajPaths = np.array(allTrajPaths)
            trajFilesExist = np.array([os.path.isfile(trajPath) for trajPath in allTrajPaths])
            trajFilesExist_sum = np.sum(trajFilesExist)
        
        #### 2.2 - If yes, load them
        if trajFilesExist_sum == PTL.NB:
            trajFilesImported = True
            trajPaths = allTrajPaths[trajFilesExist]
            for iB in range(PTL.NB):
                PTL.importTrajectories(trajPaths[iB], iB)
                # print(PTL.listTrajectories[iB].dict['X'][0], PTL.listTrajectories[iB].dict['X'][1])
            print(GREEN + 'Raw traj files found and imported :)' + NORMAL)
        
        #### 2.3 - If no, compute them by tracking the beads
        if not trajFilesImported:
            issue = PTL.buildTrajectories() 
            # Main tracking function !
            if issue == 'Bug':
                continue
            else:
                pass
        
        #### 2.4 - Save the user inputs
        PTL.saveLog(display = 0, save = True, path = logFilePath)
        
        print(BLUE + 'OK! dT = {:.3f}'.format(time.time()-Tt) + NORMAL)
        
        #### 2.5 - Sort the trajectories [Maybe unnecessary]
    
    
    #### 3. Qualify - Detect boi sizes and neighbours
        
        #### 3.1 - Infer or detect Boi sizes in the first image 
        # [Detection doesn't work well !]
        if len(PTL.beadTypes) == 1:
            if 'M450' in PTL.beadTypes[0]:
                D = 4.5
            elif 'M270' in PTL.beadTypes[0]:
                D = 2.7
            first_iF = PTL.listTrajectories[0].dict['iF'][0]
            for B in PTL.listFrames[first_iF].listBeads:
                B.D = D
        else:
            PTL.listFrames[0].detectDiameter(plot = 0)
            
        # Propagate it across the trajectories
        for iB in range(PTL.NB):
            traj = PTL.listTrajectories[iB]
            B0 = traj.dict['Bead'][0]
            D = B0.D
            traj.D = D
            for B in traj.dict['Bead']:
                B.D = D
        
        ### 3.2 - Detect neighbours
        
        # Previous way, automatic detection,
        # not robust enough
                
        # for iB in range(PTL.NB):
        #     traj = PTL.listTrajectories[iB]
        #     beadType = ''
        #     if len(PTL.beadTypes) == 1:
        #         beadType = PTL.beadTypes[0] # M270 or M450
        #     else:
        #         beadType = 'detect'
        #     traj.detectNeighbours(frequency = PTL.loop_totalSize, beadType = beadType)
        
        # Current way, with user input
        if redoAllSteps or not trajFilesImported:
            for iB in range(PTL.NB):
                traj = PTL.listTrajectories[iB]
                beadType = ''
                if len(PTL.beadTypes) == 1:
                    beadType = PTL.beadTypes[0] # M270 or M450
                elif len(PTL.beadTypes) == 2 and PTL.beadTypes[0] == PTL.beadTypes[1]:
                    beadType = PTL.beadTypes[0] # M270 or M450
                else:
                    beadType = 'detect'
                traj.detectNeighbours_ui(Nimg = PTL.nLoop, frequency = PTL.loop_totalSize, beadType = beadType)
        
        
        ### 3.3 - Detect in/out bead
        
                
        if redoAllSteps or not trajFilesImported:
            for iB in range(PTL.NB):
                traj = PTL.listTrajectories[iB]
                InOut = traj.detectInOut_ui(Nimg = PTL.nLoop, frequency = PTL.loop_totalSize)

        
    #### 4. Compute dz
                
        #### 4.1 - Import depthographs
        HDZfactor = 10
        
        if len(PTL.beadTypes) == 1:
            depthoPath = os.path.join(depthoDir, depthoNames)
#             depthoExist = os.path.exists(depthoPath+'_Deptho.tif')
            deptho = io.imread(depthoPath+'_Deptho.tif')
            depthoMetadata = pd.read_csv(depthoPath+'_Metadata.csv', sep=';')
            depthoStep = depthoMetadata.loc[0,'step']
            depthoZFocus = depthoMetadata.loc[0,'focus']
            
            # increase the resolution of the deptho with interpolation
            nX, nZ = deptho.shape[1], deptho.shape[0]
            XX, ZZ = np.arange(0, nX, 1), np.arange(0, nZ, 1)
            fd = interpolate.interp2d(XX, ZZ, deptho, kind='cubic')
            ZZ_HD = np.arange(0, nZ, 1/HDZfactor)
            depthoHD = fd(XX, ZZ_HD)
            depthoStepHD = depthoStep/HDZfactor
            depthoZFocus = depthoZFocus*HDZfactor
            #
            for iB in range(PTL.NB):
                traj = PTL.listTrajectories[iB]
                traj.deptho = depthoHD
                traj.depthoPath = depthoPath
                traj.depthoStep = depthoStepHD
                traj.depthoZFocus = depthoZFocus
        
        if len(PTL.beadTypes) > 1:
            for dN in depthoNames:
                depthoPath = os.path.join(depthoDir, dN)
                deptho = io.imread(depthoPath+'_Deptho.tif')
                depthoMetadata = pd.read_csv(depthoPath+'_Metadata.csv', sep=';')
                depthoStep = depthoMetadata.loc[0,'step']
                depthoZFocus = depthoMetadata.loc[0,'focus']
                
                # increase the resolution of the deptho with interpolation
                nX, nZ = deptho.shape[1], deptho.shape[0]
                XX, ZZ = np.arange(0, nX, 1), np.arange(0, nZ, 1)
                fd = interpolate.interp2d(XX, ZZ, deptho, kind='cubic')
                ZZ_HD = np.arange(0, nZ, 1/HDZfactor)
                depthoHD = fd(XX, ZZ_HD)
                depthoStepHD = depthoStep/HDZfactor
                depthoZFocusHD = depthoZFocus*HDZfactor
                #
                if 'M450' in dN:
                    depthoD = 4.5
                elif 'M270' in dN:
                    depthoD = 2.7
                for iB in range(PTL.NB):
                    traj = PTL.listTrajectories[iB]
                    if traj.D == depthoD:
                        traj.deptho = depthoHD
                        traj.depthoPath = depthoPath
                        traj.depthoStep = depthoStepHD
                        traj.depthoZFocus = depthoZFocusHD
        
        #### 4.2 - Compute z for each traj
        
        matchingDirection = 'upward' # Change when needed !!
        
        if redoAllSteps or not trajFilesImported:
            for iB in range(PTL.NB):
                np.set_printoptions(threshold=np.inf)
                
                print(BLUE + 'Computing Z in traj  {:.0f}...'.format(iB+1) + NORMAL)
                Tz = time.time()
                traj = PTL.listTrajectories[iB]
                traj.computeZ(matchingDirection, plot = 0)
                print(BLUE + 'OK! dT = {:.3f}'.format(time.time()-Tz) + NORMAL)
        
        else:
            print(BLUE + 'Computing Z...' + NORMAL)
            print(GREEN + 'Z had been already computed :)' + NORMAL)
        
        #### 4.3 - Save the raw traj (before Std selection)
        if redoAllSteps or not trajFilesImported:
            for iB in range(PTL.NB):
                traj = PTL.listTrajectories[iB]
                traj_df = pd.DataFrame(traj.dict)
                trajPathRaw = os.path.join(timeSeriesDataDir, 'Trajectories_raw', f[:-4] + '_rawTraj' + str(iB) + '_' + traj.beadInOut + '_PY.csv')
                traj_df.to_csv(trajPathRaw, sep = '\t', index = False)
        
        #### 4.4 - Keep only the best std data in the trajectories
        for iB in range(PTL.NB):
            traj = PTL.listTrajectories[iB]
            traj.keepBestStdOnly()
        
        #### 4.5 - The trajectories won't change from now on. We can save their '.dict' field.
        if redoAllSteps or not trajFilesImported:
            for iB in range(PTL.NB):
                traj = PTL.listTrajectories[iB]
                traj_df = pd.DataFrame(traj.dict)
                trajPath = os.path.join(timeSeriesDataDir, 'Trajectories', f[:-4] + '_traj' + str(iB) + '_' + traj.beadInOut + '_PY.csv')
                traj_df.to_csv(trajPath, sep = '\t', index = False)
    
    
    #### 5. Define pairs and compute distances
        print(BLUE + 'Computing distances...' + NORMAL)
        
        #### 5.1 - In case of 1 pair of beads
        if PTL.NB == 2:
            traj1 = PTL.listTrajectories[0]
            traj2 = PTL.listTrajectories[1]
            nT = traj1.nT
            
            #### 5.1.1 - Create a dict to prepare the export of the results
            timeSeries = {
                'idxAnalysis' : np.zeros(nT),
                'T' : np.zeros(nT),
                'Tabs' : np.zeros(nT),
                'B' : np.zeros(nT),
                'F' : np.zeros(nT),
                'dx' : np.zeros(nT),
                'dy' : np.zeros(nT),
                'dz' : np.zeros(nT),
                'D2' : np.zeros(nT),
                'D3' : np.zeros(nT)
            }

            #### 5.1.2 - Input common values:
            T0 = fieldDf['T_abs'].values[0]/1000
            timeSeries['idxAnalysis'] = traj1.dict['idxAnalysis']
            timeSeries['Tabs'] = (fieldDf['T_abs'][traj1.dict['iField']])/1000
            timeSeries['T'] = timeSeries['Tabs'].values - T0*np.ones(nT)
            timeSeries['B'] = fieldDf['B_set'][traj1.dict['iField']].values
            timeSeries['B'] *= PTL.MagCorrFactor

            #### 5.1.3 - Compute distances
            timeSeries['dx'] = (traj2.dict['X'] - traj1.dict['X'])/PTL.scale
            timeSeries['dy'] = (traj2.dict['Y'] - traj1.dict['Y'])/PTL.scale
            timeSeries['D2'] = (timeSeries['dx']**2 +  timeSeries['dy']**2)**0.5
            
            timeSeries['dz'] = (traj2.dict['Zr']*traj2.depthoStep - traj1.dict['Zr']*traj1.depthoStep)/1000
            timeSeries['dz'] *= PTL.OptCorrFactor
            timeSeries['D3'] = (timeSeries['D2']**2 +  timeSeries['dz']**2)**0.5
            
            
            #print('\n\n* timeSeries:\n')
            #print(timeSeries_DF[['T','B','F','dx','dy','dz','D2','D3']])
            print(BLUE + 'OK!' + NORMAL)
            
            
    #### 6. Compute forces
        print(BLUE + 'Computing forces...' + NORMAL)
        Tf = time.time()
        if PTL.NB == 2:
            print(GREEN + '1 pair force computation' + NORMAL)
            traj1 = PTL.listTrajectories[0]
            traj2 = PTL.listTrajectories[1]
            B0 = timeSeries['B']
            D3 = timeSeries['D3']
            dx = timeSeries['dx']
            F, dfLogF = PTL.computeForces(traj1, traj2, B0, D3, dx)
            # Main force computation function
            timeSeries['F'] = F
        
        print(BLUE + 'OK! dT = {:.3f}'.format(time.time()-Tf) + NORMAL)
            
            # Magnetization [A.m^-1]
            # M270
            # M = 0.74257*1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1)
            # M450
            # M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1);    
    
          
    #### 7. Export the results
            
        #### 7.1 - Save the tables !
        if PTL.NB == 2:
            timeSeries_DF = pd.DataFrame(timeSeries)
            timeSeriesFilePath = os.path.join(timeSeriesDataDir, f[:-4] + '_PY.csv')
            timeSeries_DF.to_csv(timeSeriesFilePath, sep = ';', index=False)
    
    print(BLUE + '\nTotal time:' + NORMAL)
    print(BLUE + str(time.time()-start) + NORMAL)
    print(BLUE + '\n' + NORMAL)
    
    plt.close('all')
    
    listTrajDicts = []
    for iB in range(PTL.NB):
        listTrajDicts.append(PTL.listTrajectories[iB].dict)
    
        #### 7.2 - Return the last objects, for optional verifications
    return(PTL, timeSeries_DF, dfLogF)








# %%%% Stand-alone functions
        
# Tracker stand-Alone Functions

# 1. XYZtracking do the tracking on any image, given 
def XYZtracking(I, cellID, NB, manipDict, depthoDir, depthoNames):
    
    PTL = PincherTimeLapse(I, cellID, manipDict, NB = NB)
    PTL.determineFramesStatus_R40()
    PTL.uiThresholding(method = 'max_entropy', factorT = 1)#, increment = 0.05)
    PTL.makeFramesList()
    PTL.detectBeads(resFileImported = False, display = 1)
    issue = PTL.buildTrajectories()
    XYZtracking_assignBeadSizes(PTL)
    for iB in range(PTL.NB):
        traj = PTL.listTrajectories[iB]
        beadType = ''
        if len(PTL.beadTypes) == 1:
            beadType = PTL.beadTypes[0] # M270 or M450
        else:
            beadType = 'detect'
        traj.detectNeighbours(frequency = PTL.loop_totalSize, beadType = beadType)
    XYZtracking_computeZ(PTL, depthoDir, depthoNames, plot = 1)
    return(PTL)

# Subfuctions of XYZtracking
def XYZtracking_assignBeadSizes(PTL):
    if len(PTL.beadTypes) == 1:
            if 'M450' in PTL.beadTypes[0]:
                D = 4.5
            elif 'M270' in PTL.beadTypes[0]:
                D = 2.7
            for B in PTL.listFrames[0].listBeads:
                B.D = D
    else:
        PTL.listFrames[0].detectDiameter(plot = 0)
            
    # Propagate it across the trajectories
    for iB in range(PTL.NB):
        traj = PTL.listTrajectories[iB]
        B0 = traj.dict['Bead'][0]
        D = B0.D
        traj.D = D
        for B in traj.dict['Bead']:
            B.D = D

def XYZtracking_computeZ(PTL, depthoDir, depthoNames, plot = 0):
    HDZfactor = 10
    if len(PTL.beadTypes) == 1:
        depthoPath = os.path.join(depthoDir, depthoNames)
#             depthoExist = os.path.exists(depthoPath+'_Deptho.tif')
        deptho = io.imread(depthoPath+'_Deptho.tif')
        depthoMetadata = pd.read_csv(depthoPath+'_Metadata.csv', sep=';')
        depthoStep = depthoMetadata.loc[0,'step']
        depthoZFocus = depthoMetadata.loc[0,'focus']

        # increase the resolution of the deptho with interpolation
        nX, nZ = deptho.shape[1], deptho.shape[0]
        XX, ZZ = np.arange(0, nX, 1), np.arange(0, nZ, 1)
        fd = interpolate.interp2d(XX, ZZ, deptho, kind='cubic')
        ZZ_HD = np.arange(0, nZ, 1/HDZfactor)
        depthoHD = fd(XX, ZZ_HD)
        depthoStepHD = depthoStep/HDZfactor
        depthoZFocusHD = depthoZFocus*HDZfactor
        #
        for iB in range(PTL.NB):
            traj = PTL.listTrajectories[iB]
            traj.deptho = depthoHD
            traj.depthoPath = depthoPath
            traj.depthoStep = depthoStepHD
            traj.depthoZFocus = depthoZFocusHD
            plt.imshow(depthoHD)

    if len(PTL.beadTypes) > 1:
        for dN in depthoNames:
            depthoPath = os.path.join(depthoDir, dN)
            deptho = io.imread(depthoPath+'_Deptho.tif')
            depthoMetadata = pd.read_csv(depthoPath+'_Metadata.csv', sep=';')
            depthoStep = depthoMetadata.loc[0,'step']
            depthoZFocus = depthoMetadata.loc[0,'focus']

            # increase the resolution of the deptho with interpolation
            nX, nZ = deptho.shape[1], deptho.shape[0]
            XX, ZZ = np.arange(0, nX, 1), np.arange(0, nZ, 1)
            fd = interpolate.interp2d(XX, ZZ, deptho, kind='cubic')
            ZZ_HD = np.arange(0, nZ, 1/HDZfactor)
            depthoHD = fd(XX, ZZ_HD)
            depthoStepHD = depthoStep/HDZfactor
            depthoZFocusHD = depthoZFocus*HDZfactor
            #
            if 'M450' in dN:
                depthoD = 4.5
            elif 'M270' in dN:
                depthoD = 2.7
            for iB in range(PTL.NB):
                traj = PTL.listTrajectories[iB]
                if traj.D == depthoD:
                    traj.deptho = depthoHD
                    traj.depthoPath = depthoPath
                    traj.depthoStep = depthoStepHD
                    traj.depthoZFocus = depthoZFocusHD
                    
    # Compute z for each traj
    matchingDirection = 'downward'
    
    for iB in range(PTL.NB):
        np.set_printoptions(threshold=np.inf)
        
        print('Computing Z in traj  {:.0f}...'.format(iB+1))
        
        Tz = time.time()
        traj = PTL.listTrajectories[iB]
        traj.computeZ(matchingDirection, plot)
        print('OK! dT = {:.3f}'.format(time.time()-Tz))

    # Keep only the best std data in the trajectories
    for iB in range(PTL.NB):
        print('Picking the images with the best StdDev in traj  {:.0f}...'.format(iB+1))
        
        Tstd = time.time()
        traj = PTL.listTrajectories[iB]
        traj.keepBestStdOnly()
        print('OK! dT = {:.3f}'.format(time.time()-Tstd))
    
        
    
# %% (3) Depthograph making classes & functions
        
# %%%% BeadDeptho
        
class BeadDeptho:
    def __init__(self, I, X0, Y0, S0, bestZ, scale, beadType, fileName):
        
        nz, ny, nx = I.shape[0], I.shape[1], I.shape[2]
        
        self.I = I
        self.nz = nz
        self.ny = ny
        self.nx = nx
        self.scale = scale
        self.X0 = X0
        self.Y0 = Y0
        self.S0 = S0
        self.XYm = np.zeros((self.nz, 2))
        self.XYm[S0-1, 0] = X0
        self.XYm[S0-1, 1] = Y0
        self.fileName = fileName
        
        self.beadType = beadType
        self.D0 = 4.5 * (beadType == 'M450') + 2.7 * (beadType == 'M270')
#         self.threshold = threshold
        self.I_cleanROI = np.array([])
#         self.cleanROI = np.zeros((self.nz, 4), dtype = int)
        
        self.validBead = True
        self.iValid = -1
        
        self.bestZ = bestZ
        self.validSlice = np.zeros(nz, dtype = bool)
        self.zFirst = 0
        self.zLast = nz
        self.validDepth = nz
        
        self.valid_v = True
        self.valid_h = True
        self.depthosDict = {}
        self.profileDict = {}
        self.ZfocusDict = {}

        
    def buildCleanROI(self, plot = 0):
        # Determine if the bead is to close to the edge on the max frame
        D0 = self.D0 + 4.5*(self.D0 == 0)
        roughSize = np.floor(1.2*D0*self.scale)
        mx, Mx = np.min(self.X0 - 0.5*roughSize), np.max(self.X0 + 0.5*roughSize)
        my, My = np.min(self.Y0 - 0.5*roughSize), np.max(self.Y0 + 0.5*roughSize)
        testImageSize = mx > 0 and Mx < self.nx and my > 0 and My < self.ny
        
        # Aggregate the different validity test (for now only 1)
        validBead = testImageSize
        
        # If the bead is valid we can proceed
        if validBead:
            # Detect or infer the size of the beads we are measuring
            if self.beadType == 'detect' or self.D0 == 0:
                counts, binEdges = np.histogram(self.I[self.z_max,my:My,mx:Mx].ravel(), bins=256)
                peaks, peaksProp = find_peaks(counts, height=100, threshold=None, distance=None, prominence=None, \
                                   width=None, wlen=None, rel_height=0.5, plateau_size=None)
                peakThreshVal = 1000
                if counts[peaks[0]] > peakThreshVal:
                    self.D0 = 4.5
                    self.beadType = 'M450'
                else:
                    self.D0 = 2.7
                    self.beadType = 'M270'
        else:
            self.validBead = False
        
        if validBead:
            for z in range(self.bestZ, -1, -1):
                if not z in self.S0:
                    break
            zFirst = z
            for z in range(self.bestZ, self.nz, +1):
                if not z in self.S0:
                    break
            zLast = z-1
            
            roughSize = int(np.floor(1.15*self.D0*self.scale))
            roughSize += 1 + roughSize%2
            roughCenter = int((roughSize+1)//2)

            cleanSize = getDepthoCleanSize(self.D0, self.scale)
            
            I_cleanROI = np.zeros([self.nz, cleanSize, cleanSize])

            try:
                for i in range(zFirst, zLast):
                    xmi, ymi = self.XYm[i,0], self.XYm[i,1]
                    x1, y1, x2, y2, validBead = getROI(roughSize, xmi, ymi, self.nx, self.ny)
                    if not validBead:
                        if x1 < 0 or x2 > self.nx:
                            self.valid_h = False
                        if y1 < 0 or y2 > self.ny:
                            self.valid_v = False

        #                 fig, ax = plt.subplots(1,2)
        #                 ax[0].imshow(self.I[i])
                    xm1, ym1 = xmi-x1, ymi-y1
                    I_roughRoi = self.I[i,y1:y2,x1:x2]
        #                 ax[1].imshow(I_roughRoi)
        #                 fig.show()

                    translation = (xm1-roughCenter, ym1-roughCenter)

                    tform = transform.EuclideanTransform(rotation=0, \
                                                         translation = (xm1-roughCenter, ym1-roughCenter))

                    I_tmp = transform.warp(I_roughRoi, tform, order = 1, preserve_range = True)

                    I_cleanROI[i] = np.copy(I_tmp[roughCenter-cleanSize//2:roughCenter+cleanSize//2+1,\
                                                  roughCenter-cleanSize//2:roughCenter+cleanSize//2+1])

                if not self.valid_v and not self.valid_h:
                    self.validBead = False

                else:
                    self.zFirst = zFirst
                    self.zLast = zLast
                    self.validDepth = zLast-zFirst
                    self.I_cleanROI = I_cleanROI.astype(np.uint16)

                # VISUALISE
                if plot >= 2:
                    for i in range(zFirst, zLast, 50):
                        self.plotROI(i)

            except:
                print('Error for the file: ' + self.fileName)


    def buildDeptho(self, plot = 0):
        preferedDeptho = 'v'
        side_ROI = self.I_cleanROI.shape[1]
        mid_ROI = side_ROI//2
        nbPixToAvg = 3 # Have to be an odd number
        deptho_v = np.zeros([self.nz, side_ROI], dtype = np.float64)
        deptho_h = np.zeros([self.nz, side_ROI], dtype = np.float64)
        deptho_HD = np.zeros([self.nz, side_ROI*5], dtype = np.float64)
        
        if self.valid_v:
            for z in range(self.zFirst, self.zLast):
                templine = side_ROI
                deptho_v[z] = self.I_cleanROI[z,:,mid_ROI] * (1/nbPixToAvg)
                for i in range(1, 1 + nbPixToAvg//2):
                    deptho_v[z] += self.I_cleanROI[z,:,mid_ROI - i] * (1/nbPixToAvg) 
                    deptho_v[z] += self.I_cleanROI[z,:,mid_ROI + i] * (1/nbPixToAvg)
            deptho_v = deptho_v.astype(np.uint16)
            self.depthosDict['deptho_v'] = deptho_v
        
        if self.valid_h:
            for z in range(self.zFirst, self.zLast):
                templine = side_ROI
                deptho_h[z] = self.I_cleanROI[z,mid_ROI,:] * (1/nbPixToAvg)
                for i in range(1, 1 + nbPixToAvg//2):
                    deptho_h[z] += self.I_cleanROI[z,mid_ROI - i,:] * (1/nbPixToAvg) 
                    deptho_h[z] += self.I_cleanROI[z,mid_ROI + i,:] * (1/nbPixToAvg)
            deptho_h = deptho_h.astype(np.uint16)
            self.depthosDict['deptho_h'] = deptho_h
        
        if preferedDeptho == 'v' and not self.valid_v:
            hdDeptho = 'h'
        elif preferedDeptho == 'h' and not self.valid_h:
            hdDeptho = 'v'
        else:
            hdDeptho = preferedDeptho
            
        if hdDeptho == 'v':
            for z in range(self.zFirst, self.zLast):
                x = np.arange(mid_ROI - 2, mid_ROI + 3)
                y = np.arange(0, side_ROI)
#                 xx, yy = np.meshgrid(x, y)
                vals = self.I_cleanROI[z, :, mid_ROI-2:mid_ROI+3]
                f = interpolate.interp2d(x, y, vals, kind='cubic')
                # Now use the obtained interpolation function and plot the result:

                xnew = x
                ynew = np.arange(0, side_ROI, 0.2)
                vals_new = f(xnew, ynew)
                deptho_HD[z] = vals_new[:,5//2] * (1/nbPixToAvg)
                for i in range(1, 1 + nbPixToAvg//2):
                    deptho_HD[z] += vals_new[:,5//2-i] * (1/nbPixToAvg)
                    deptho_HD[z] += vals_new[:,5//2+i] * (1/nbPixToAvg)              
#                 if z == self.z_max:
#                     figInterp, axesInterp = plt.subplots(1,2)
#                     axesInterp[0].imshow(vals)
#                     axesInterp[0].plot([5//2, 5//2], [0, vals.shape[0]], 'r--')
#                     axesInterp[1].imshow(vals_new)
#                     axesInterp[1].plot([5//2, 5//2], [0, vals_new.shape[0]], 'r--')
#                     figInterp.show()     
            deptho_HD = deptho_HD.astype(np.uint16)    
            self.depthosDict['deptho_HD'] = deptho_HD
                
        elif hdDeptho == 'h':
            for z in range(self.zFirst, self.zLast):
                x = np.arange(0, side_ROI)
                y = np.arange(mid_ROI - 2, mid_ROI + 3)
#                 xx, yy = np.meshgrid(x, y)
                vals = self.I_cleanROI[z, mid_ROI-2:mid_ROI+3, :]
                f = interpolate.interp2d(x, y, vals, kind='cubic')
                # Now use the obtained interpolation function and plot the result:

                xnew = np.arange(0, side_ROI, 0.2)
                ynew = y
                vals_new = f(xnew, ynew)
                deptho_HD[z] = vals_new[5//2,:] * (1/nbPixToAvg)
                for i in range(1, 1 + nbPixToAvg//2):
                    deptho_HD[z] += vals_new[5//2-i,:] * (1/nbPixToAvg)
                    deptho_HD[z] += vals_new[5//2+i,:] * (1/nbPixToAvg)
#                 if z == self.z_max:
#                     figInterp, axesInterp = plt.subplots(1,2)
#                     axesInterp[0].imshow(vals)
#                     axesInterp[0].plot([0, vals.shape[1]], [5//2, 5//2], 'r--')
#                     axesInterp[1].imshow(vals_new)
#                     axesInterp[1].plot([0, vals_new.shape[1]], [5//2, 5//2], 'r--')
#                     figInterp.show()
            deptho_HD = deptho_HD.astype(np.uint16)    
            self.depthosDict['deptho_HD'] = deptho_HD
        
        # 3D caracterisation
#         I_binary = np.zeros([self.I_cleanROI.shape[0], self.I_cleanROI.shape[1], self.I_cleanROI.shape[2]])
#         I_binary[self.zFirst:self.zLast] = (self.I_cleanROI[self.zFirst:self.zLast] > self.threshold)
#         Zm3D, Ym3D, Xm3D = ndi.center_of_mass(self.I_cleanROI, labels=I_binary, index=1)
#         self.ZfocusDict['Zm3D'] = Zm3D
        
        # Raw profiles
        mid_ROI_HD = deptho_HD.shape[1]//2
        Z = np.array([z for z in range(self.I_cleanROI.shape[0])])
#         intensity_tot = np.array([np.sum(self.I_cleanROI[z][I_binary[z].astype(bool)])/(1+np.sum(I_binary[z])) for z in range(self.I_cleanROI.shape[0])]).astype(np.float64)
        intensity_v = np.array([np.sum(deptho_v[z,:])/side_ROI for z in range(deptho_v.shape[0])]).astype(np.float64)
        intensity_h = np.array([np.sum(deptho_h[z,:])/side_ROI for z in range(deptho_h.shape[0])]).astype(np.float64)
        intensity_HD = np.array([np.sum(deptho_HD[z,mid_ROI_HD-5:mid_ROI_HD+6])/11 for z in range(deptho_HD.shape[0])]).astype(np.float64)
#         
        Zm_v, Zm_h = np.argmax(intensity_v), np.argmax(intensity_h)
#         Zm_tot = np.argmax(intensity_tot)
        Zm_HD = np.argmax(intensity_HD)
        
        self.profileDict['intensity_v'] = intensity_v
        self.profileDict['intensity_h'] = intensity_h
        self.profileDict['intensity_HD'] = intensity_HD
#         self.profileDict['intensity_tot'] = intensity_tot
        self.ZfocusDict['Zm_v'] = Zm_v
        self.ZfocusDict['Zm_h'] = Zm_h
        self.ZfocusDict['Zm_HD'] = Zm_HD
#         self.ZfocusDict['Zm_tot'] = Zm_tot


        # Smoothed profiles
        Z_hd = np.arange(0, self.I_cleanROI.shape[0], 0.2)
        intensity_v_hd = np.interp(Z_hd, Z, intensity_v)
        intensity_h_hd = np.interp(Z_hd, Z, intensity_h)
        intensity_HD_hd = np.interp(Z_hd, Z, intensity_HD)
#         intensity_tot_hd = np.interp(Z_hd, Z, intensity_tot)
        
        intensity_v_smooth = savgol_filter(intensity_v_hd, 101, 5)
        intensity_h_smooth = savgol_filter(intensity_h_hd, 101, 5)
        intensity_HD_smooth = savgol_filter(intensity_HD_hd, 101, 5)
#         intensity_tot_smooth = savgol_filter(intensity_tot_hd, 101, 5)
        
        Zm_v_hd, Zm_h_hd = Z_hd[np.argmax(intensity_v_smooth)], Z_hd[np.argmax(intensity_h_smooth)]
#         Zm_tot_hd = Z_hd[np.argmax(intensity_tot_smooth)]
        Zm_HD_hd = Z_hd[np.argmax(intensity_HD_smooth)]
        
        self.profileDict['intensity_v_smooth'] = intensity_v_smooth 
        self.profileDict['intensity_h_smooth'] = intensity_h_smooth
        self.profileDict['intensity_HD_smooth'] = intensity_HD_smooth
#         self.profileDict['intensity_tot_smooth'] = intensity_tot_smooth
        self.ZfocusDict['Zm_v_hd'] = Zm_v_hd 
        self.ZfocusDict['Zm_h_hd'] = Zm_h_hd
        self.ZfocusDict['Zm_HD_hd'] = Zm_HD_hd
#         self.ZfocusDict['Zm_tot_hd'] = Zm_tot_hd
        
        # VISUALISE
        if plot >= 2:
            self.plotProfiles()            
            
            
    def saveBeadDeptho(self, path, manipID, step, bestDetphoType = 'HD', bestFocusType = 'HD_hd'):
        supDataDir = manipID + '_supData'
        supDataDirPath = os.path.join(path, supDataDir)
        if not os.path.exists(supDataDirPath):
            os.makedirs(supDataDirPath)
        
        cleanROIName = manipID + '_cleanROI.tif'
        cleanROIPath = os.path.join(path, cleanROIName)
        io.imsave(cleanROIPath, self.I_cleanROI)
        
        profilesRaw_keys = ['intensity_v', 'intensity_h', 'intensity_HD'] #, 'intensity_tot']
        profileDictRaw = {k: self.profileDict[k] for k in profilesRaw_keys}
        profileDictRaw_df = pd.DataFrame(profileDictRaw)
        profileDictRaw_df.to_csv(os.path.join(supDataDirPath, 'profiles_raw.csv'))
        
        profilesSmooth_keys = ['intensity_v_smooth', 'intensity_h_smooth', 'intensity_HD_smooth'] #, 'intensity_tot_smooth']
        profileDictSmooth = {k: self.profileDict[k] for k in profilesSmooth_keys}
        profileDictSmooth_df = pd.DataFrame(profileDictSmooth)
        profileDictSmooth_df.to_csv(os.path.join(supDataDirPath, 'profiles_smooth.csv'))
        
        ZfocusDict_df = pd.DataFrame(self.ZfocusDict, index = [1])
        ZfocusDict_df.to_csv(os.path.join(supDataDirPath, 'Zfoci.csv'))
        
        bestFocus = self.ZfocusDict['Zm_' + bestFocusType]
        metadataPath = os.path.join(path, manipID + '_Metadata.csv')
        with open(metadataPath, 'w') as f:
            f.write('step;bestFocus')
#             for k in self.ZfocusDict.keys():
#                 f.write(';')
#                 f.write(k)
            f.write('\n')
            f.write(str(step) + ';' + str(bestFocus))
#             for k in self.ZfocusDict.keys():
#                 f.write(';')
#                 f.write(str(self.ZfocusDict[k]))
                
        depthoPath = os.path.join(path, manipID + '_deptho.tif')
        bestDeptho = self.depthosDict['deptho_' + bestDetphoType]
        io.imsave(depthoPath, bestDeptho)
        
        
        
# Plot functions
        
    def plotXYm(self):
        fig, ax = plt.subplots(1,1)
        pStart, pStop = np.percentile(self.I[self.z_max], (1, 99))
        ax.imshow(self.I[self.z_max], cmap = 'gray', vmin = pStart, vmax = pStop)
        ax.plot(self.XYm[self.validSlice,0],self.XYm[self.validSlice,1],'r-')
        fig.show()
        
    def plotROI(self, i = 'auto'):
        if i == 'auto':
            i = self.z_max
        
        fig, ax = plt.subplots(1,3, figsize = (16,4))
        
        xm, ym = np.mean(self.XYm[self.validSlice,0]),  np.mean(self.XYm[self.validSlice,1])
        ROIsize_x = self.D*1.25*self.scale + (max(self.XYm[self.validSlice,0])-min(self.XYm[self.validSlice,0]))
        ROIsize_y = self.D*1.25*self.scale + (max(self.XYm[self.validSlice,1])-min(self.XYm[self.validSlice,1]))
        x1_ROI, y1_ROI, x2_ROI, y2_ROI = int(xm - ROIsize_x//2), int(ym - ROIsize_y//2), int(xm + ROIsize_x//2), int(ym + ROIsize_y//2)        

        pStart, pStop = np.percentile(self.I[i], (1, 99))
        ax[0].imshow(self.I[i], cmap = 'gray', vmin = pStart, vmax = pStop)
        ax[0].plot([x1_ROI,x1_ROI], [y1_ROI,y2_ROI], 'c--')
        ax[0].plot([x1_ROI,x2_ROI], [y2_ROI,y2_ROI], 'c--')
        ax[0].plot([x2_ROI,x2_ROI], [y1_ROI,y2_ROI], 'c--')
        ax[0].plot([x1_ROI,x2_ROI], [y1_ROI,y1_ROI], 'c--')

        I_ROI = self.I[i,y1_ROI:y2_ROI,x1_ROI:x2_ROI] 
        pStart, pStop = np.percentile(I_ROI, (1, 99))
        ax[1].imshow(I_ROI, cmap = 'gray', vmin = pStart, vmax = pStop)
        ax[1].plot(self.XYm[self.validSlice,0]-x1_ROI, self.XYm[self.validSlice,1]-y1_ROI, 'r-', lw=0.75)
        ax[1].plot(self.XYm[i,0]-x1_ROI, self.XYm[i,1]-y1_ROI, 'b+', lw=0.75)
        
        pStart, pStop = np.percentile(self.I_cleanROI[i], (1, 99))
        mid = self.I_cleanROI[i].shape[0]//2
        I_cleanROI_binary = (self.I_cleanROI[i] > self.threshold)
        y, x = ndi.center_of_mass(self.I_cleanROI[i], labels=I_cleanROI_binary, index=1)
        ax[2].imshow(self.I_cleanROI[i], cmap = 'gray', vmin = pStart, vmax = pStop)
        ax[2].plot([0,2*mid],[mid, mid], 'r--', lw = 0.5)
        ax[2].plot([mid, mid],[0,2*mid], 'r--', lw = 0.5)
        ax[2].plot([x],[y], 'b+')
        fig.show()

    def plotProfiles(self):
        Z = np.array([z for z in range(self.I_cleanROI.shape[0])])
        Z_hd = np.arange(0, self.I_cleanROI.shape[0], 0.2)
        intensity_v = self.profileDict['intensity_v']
        intensity_h = self.profileDict['intensity_h']
        intensity_HD = self.profileDict['intensity_HD']
        intensity_tot = self.profileDict['intensity_tot']
        Zm_v = self.ZfocusDict['Zm_v']
        Zm_h = self.ZfocusDict['Zm_h']
        Zm_HD = self.ZfocusDict['Zm_HD']
        Zm_tot = self.ZfocusDict['Zm_tot']
        intensity_v_smooth = self.profileDict['intensity_v_smooth']
        intensity_h_smooth = self.profileDict['intensity_h_smooth']
        intensity_HD_smooth = self.profileDict['intensity_HD_smooth']
        intensity_tot_smooth = self.profileDict['intensity_tot_smooth']
        Zm_v_hd = self.ZfocusDict['Zm_v_hd']
        Zm_h_hd = self.ZfocusDict['Zm_h_hd']
        Zm_HD_hd = self.ZfocusDict['Zm_HD_hd']
        Zm_tot_hd = self.ZfocusDict['Zm_tot_hd']
        
#             fig, ax = plt.subplots(1,3, figsize = (12, 4))
#             ax[0].plot(Z, intensity_v)
#             ax[1].plot(Z, intensity_h)
#             ax[2].plot(Z, (intensity_tot))
#             ax[0].plot([Zm_v, Zm_v], [0, ax[0].get_ylim()[1]], 'r--', lw = 0.8, label = 'Zm_v = {:.2f}'.format(Zm_v))
#             ax[1].plot([Zm_h, Zm_h], [0, ax[1].get_ylim()[1]], 'r--', lw = 0.8, label = 'Zm_h = {:.2f}'.format(Zm_h))
#             ax[2].plot([Zm_tot, Zm_tot], [0, ax[2].get_ylim()[1]], 'r--', lw = 0.8, label = 'Zm_tot = {:.2f}'.format(Zm_tot))
#             ax[0].legend(loc = 'lower right')
#             ax[1].legend(loc = 'lower right')
#             ax[2].legend(loc = 'lower right')
        
        fig, ax = plt.subplots(1,4, figsize = (16, 4))
        ax[0].plot(Z, intensity_v, 'b-')
        ax[1].plot(Z, intensity_h, 'b-')
        ax[2].plot(Z, intensity_HD, 'b-')
        ax[3].plot(Z, (intensity_tot), 'b-')
        ax[0].plot(Z_hd, intensity_v_smooth, 'k--')
        ax[1].plot(Z_hd, intensity_h_smooth, 'k--')
        ax[2].plot(Z_hd, intensity_HD_smooth, 'k--')
        ax[3].plot(Z_hd, intensity_tot_smooth, 'k--')
        ax[0].plot([Zm_v_hd, Zm_v_hd], [0, ax[0].get_ylim()[1]], 'r--', lw = 0.8, label = 'Zm_v_hd = {:.2f}'.format(Zm_v_hd))
        ax[1].plot([Zm_h_hd, Zm_h_hd], [0, ax[1].get_ylim()[1]], 'r--', lw = 0.8, label = 'Zm_h_hd = {:.2f}'.format(Zm_h_hd))
        ax[2].plot([Zm_HD_hd, Zm_HD_hd], [0, ax[2].get_ylim()[1]], 'r--', lw = 0.8, label = 'Zm_HD_hd = {:.2f}'.format(Zm_HD_hd))
        ax[3].plot([Zm_tot_hd, Zm_tot_hd], [0, ax[3].get_ylim()[1]], 'r--', lw = 0.8, label = 'Zm_tot_hd = {:.2f}'.format(Zm_tot_hd))
        ax[0].legend(loc = 'lower right')
        ax[1].legend(loc = 'lower right')
        ax[2].legend(loc = 'lower right')
        ax[3].legend(loc = 'lower right')

        #         print('Zm_v = {:.2f}, Zm_h = {:.2f}, Zm_tot = {:.2f}'\
        #               .format(Zm_v, Zm_h, Zm_tot))
        #         print('Zm_v_hd = {:.2f}, Zm_h_hd = {:.2f}, Zm_tot_hd = {:.2f}'\
        #               .format(Zm_v_hd, Zm_h_hd, Zm_tot_hd))
        
        fig.show()
        
        
    def plotDeptho(self, d = 'HD'):
        fig, ax = plt.subplots(1,1, figsize = (4, 6))
        D = self.depthosDict['deptho_' + d]
        z_focus = self.ZfocusDict['Zm_' + d + '_hd']
        ny, nx = D.shape[0], D.shape[1]
        pStart, pStop = np.percentile(D, (1, 99))
        pStop = pStop + 0.3 * (pStop-pStart)
        ax.imshow(D, cmap='plasma', vmin = pStart, vmax = pStop)
        ax.plot([0, nx], [self.zFirst, self.zFirst], 'r--')
        ax.text(nx//2, self.zFirst - 10, str(self.zFirst), c = 'r')
        ax.plot([0, nx], [self.zLast, self.zLast], 'r--')
        ax.text(nx//2, self.zLast - 10, str(self.zLast), c = 'r')
        ax.plot([nx//2], [z_focus], 'c+')
        ax.text(nx//2, z_focus - 10, str(z_focus), c = 'c')
        fig.suptitle('File ' + self.fileName + ' - Bead ' + str(self.iValid))
        fig.show()
        
# %%%% depthoMaker
        
def depthoMaker(dirPath, savePath, specif, saveLabel, scale, beadType = 'M450', step = 20, d = 'HD', plot = 0):
    rawFileList = os.listdir(dirPath)
    listFileNames = [f[:-4] for f in rawFileList if (os.path.isfile(os.path.join(dirPath, f)) and f.endswith(".tif"))]
    L = []
    
    for f in listFileNames:
        test1 = (specif in f) or (specif == 'all')
        test2 = ((f + '_Results.txt') in os.listdir(dirPath))
        valid = test1 and test2
        if valid:
            L.append(f)
    
    listFileNames = L
    
    listBD = []
#     dictBD = {}
    
#     print(listFileNames)
    for f in listFileNames:
        filePath = os.path.join(dirPath, f)
        I = io.imread(filePath + '.tif')
        resDf = pd.read_csv((filePath + '_Results.txt'), sep = '\t').drop(columns = [' '])
        # Area,StdDev,XM,YM,Slice
        X0 = resDf['XM'].values
        Y0 = resDf['YM'].values
        S0 = resDf['Slice'].values
        bestZ = S0[np.argmax(resDf['StdDev'].values)] - 1 # The index of the image with the highest Std
        # This image will be more or less the one with the brightest spot
        
        # Create the BeadDeptho object
        BD = BeadDeptho(I, X0, Y0, S0, bestZ, scale, beadType, f)
        
        # Creation of the clean ROI where the center of mass is always perfectly centered.
        BD.buildCleanROI(plot)

        # If the bead was not acceptable (for instance too close to the edge of the image)
        # then BD.validBead will be False
        if not BD.validBead:
            print(RED + 'Not acceptable file: ' + f + NORMAL)
            
        # Else, we can proceed.
        else: 
            print(BLUE + 'Job done for the file: ' + f + NORMAL)

            # Creation of the z profiles
            BD.buildDeptho(plot)
        
        listBD.append(BD)
        i = 1
        for BD in listBD:
#             BD_manipID = findInfosInFileName(BD.fileName, 'manipID')
            subFileSavePath = os.path.join(savePath, 'Intermediate_Py', saveLabel + '_step' + str(step))
            BD.saveBeadDeptho(subFileSavePath, specif +  '_' + str(i), step = step, bestDetphoType = 'HD', bestFocusType = 'HD_hd')
            i += 1

#
## If different sizes of bead at once #
#
#         for BD in listBD:
#             if BD.D0 not in dictBD.keys():
#                 dictBD[BD.D0] = [BD]
#             else:
#                 dictBD[BD.D0].append(BD)
# 
#     for size in dictBD.keys():
#         listBD = dictBD[size]
#     ... go on with the code below with an indent added !


    maxAboveZm, maxBelowZm = 0, 0
    for BD in listBD:
        Zm = int(np.round(BD.ZfocusDict['Zm_' + d + '_hd']))
        if Zm - BD.zFirst > maxAboveZm:
            maxAboveZm = Zm - BD.zFirst
        if BD.zLast - Zm > maxBelowZm:
            maxBelowZm = BD.zLast - Zm
    maxAboveZm, maxBelowZm = int(maxAboveZm), int(maxBelowZm)
    Zfocus = maxAboveZm
    depthoWidth = listBD[0].depthosDict['deptho_' + d].shape[1]
    depthoHeight = maxAboveZm + maxBelowZm
    finalDeptho = np.zeros([depthoHeight, depthoWidth], dtype = np.float64)

    for z in range(1, maxAboveZm+1):
        count = 0
        for BD in listBD:
            Zm = int(np.round(BD.ZfocusDict['Zm_' + d + '_hd']))
            currentDeptho = BD.depthosDict['deptho_' + d]
            if Zm-z >= 0 and np.sum(currentDeptho[Zm-z,:] != 0):
                count += 1
        for BD in listBD:
            Zm = int(np.round(BD.ZfocusDict['Zm_' + d + '_hd']))
            currentDeptho = BD.depthosDict['deptho_' + d]
            if Zm-z >= 0 and np.sum(currentDeptho[Zm-z,:] != 0):
                finalDeptho[Zfocus-z,:] += currentDeptho[Zm-z,:]/count

    for z in range(0, maxBelowZm):
        count = 0
        for BD in listBD:
            Zm = int(np.round(BD.ZfocusDict['Zm_' + d + '_hd']))
            currentDeptho = BD.depthosDict['deptho_' + d]
#             print(currentDeptho.shape)
            if Zm+z >= 0 and Zm+z < currentDeptho.shape[0] and np.sum(currentDeptho[Zm+z,:] != 0):
                count += 1
        for BD in listBD:
            Zm = int(np.round(BD.ZfocusDict['Zm_' + d + '_hd']))
            currentDeptho = BD.depthosDict['deptho_' + d]
            if Zm+z >= 0 and Zm+z < currentDeptho.shape[0] and np.sum(currentDeptho[Zm+z,:] != 0):
                finalDeptho[Zfocus+z,:] += currentDeptho[Zm+z,:]/count

    finalDeptho = finalDeptho.astype(np.uint16)

    fig, ax = plt.subplots(1,1)
    ax.imshow(finalDeptho)
    
    fig.suptitle(beadType)
    fig.show()

    depthoSavePath = os.path.join(savePath, saveLabel + '_Deptho.tif')
    io.imsave(depthoSavePath, finalDeptho)
    metadataPath = os.path.join(savePath, saveLabel + '_Metadata.csv')
    with open(metadataPath, 'w') as f:
        f.write('step;focus')
        f.write('\n')
        f.write(str(step) + ';' + str(Zfocus))

    print(GREEN + 'ok' + NORMAL)
               
        
# Finished !