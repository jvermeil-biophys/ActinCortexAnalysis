datafolder = 'D:\Matlab Analysis\Data_Joseph\Raw\20.06.30';
depthoDirPath = 'D:\Matlab Analysis\Data_Joseph\Raw\';
depthoname = '15-06-20_Depthograph63x';

% scale=1/15.8; % 100X
scale=1/9.7; % 63X
% En µm par pixel

isImageTriplet = 1;
filename = '30-06-20_FallingChain05bis_37Mag';
NoiseAnalysis_NoiseFromBeadsChains_J(datafolder,filename,depthoDirPath,depthoname,isImageTriplet,scale)

% isImageTriplet = 1;
% filename = '30-06-20_FallingChain06_37Mag';
% NoiseAnalysis_NoiseFromBeadsChains_J(datafolder,filename,depthoDirPath,depthoname,isImageTriplet,scale)
% 
% isImageTriplet = 0;
% filename = '30-06-20_LandedChainZstack01_37Mag';
% getNoiseFromBeadsChains_J(datafolder,filename,depthoDirPath,depthoname,isImageTriplet,scale)

% isImageTriplet = 0;
% filename = '30-06-20_LandedChainZstack02_37Mag';
% getNoiseFromBeadsChains_J(datafolder,filename,depthoDirPath,depthoname,isImageTriplet,scale)

% isImageTriplet = 0;
% filename = '30-06-20_LandedChainZstack03_37Mag';
% getNoiseFromBeadsChains_J(datafolder,filename,depthoDirPath,depthoname,isImageTriplet,scale)