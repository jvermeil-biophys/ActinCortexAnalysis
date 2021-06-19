close all;
clear all;
clc;

warning('off')

%% Paths & options
%% Paths
RawdataFolder = 'D:\MagneticPincherData\Raw';
MatfileFolder = 'D:\MagneticPincherData\MatFiles';
FigureFolder = 'D:\MagneticPincherData\Figures';
ExperimentalDataFolder = 'C:\Users\JosephVermeil\Desktop\ActinCortexAnalysis\ExperimentalData';
ExportDataFolder = 'C:\Users\JosephVermeil\Desktop\ActinCortexAnalysis\DataAnalysis';
DropboxDataFolder = '';

%% Options
set(0, 'defaultfigureposition', get(0, 'Screensize'));

set(0,'DefaultFigureWindowStyle','docked')

% Couleurs pour dictys
Cwt = [20 109 159]./255; % WT
Csev = [182 109 106]./255; % SevA- mutant
Cfim  = [255 255 109]./255; % FimA- mutant
CaA = [0 73 73]./255; % abpA- mutant
CaC  = [255 182 119]./255; % abpC-
Cwt27 = [146 109 219]./255./2; % WT M270
Cdm = [146 0 0]./255; % DMSO
Cl  = [150 0 255]./255; % Latranculin

% Couleurs pour 3T3
Cwt = [20 109 159]./255; % WT
CaSFLnodrug = 'y'; % aSFLnodrug
CaSFLdoxy = [255 255 109]./255; % aSFLdoxy
CaSFLnodrug = [255 255 109]./255; % aSFLnodrug
CaSFLdoxy = 'r'; % aSFLdoxy

% diff rates
C02 = [30 240 60]./255;
C05 = [30 220 60]./255;
C1 = [30 200 60]./255;
C2 = [30 180 60]./255;
C4 = [30 160 60]./255;
C7 = [30 140 60]./255;
C12 = [30 120 60]./255;

C02_27 = [240 30 60]./255;
C05_27 = [220 30 60]./255;
C1_27 = [200 30 60]./255;
C2_27 = [180 30 60]./255;
C4_27 = [160 30 60]./255;
C7_27 = [140 30 60]./255;
C12_27 = [120 30 60]./255;


M02 = 'o';
M05 = 's';
M1 = 'd';
M2 = '^';
M4 = 'v';
M7 = '>';
M12 = '<';

% diff tmp comp
Cdc1  = [255 127 63]./255; % diff tmp comp
Cdc2  = [127 63 255]./255; % diff tmp comp
Cdc3  = [63 255 127]./255; % diff tmp comp

% for r2v
AUTO = 0; % automode

% for v2d
NotSaved = {};
PLOTV2D = 0; % plot choice

% for d2p
PLOTDP = 0; % plot choice

% for d2p and p2s
alignlvl = 0.9; 

%for p2s
PLOTPS = 0; % plot choice

% graph choice for s2p
CORT = 1;
PEAKS = 0; 
ACTI = 1;
DIST = 1;
TIME = 0;

% for d2m
PLOTD2M = 1; % plot choice
VERBOSED2M = 0; % text display choice



% M2P

fitparams = 'Strain100-R20.9';

% Load experimental conditions file
pathExperimentalConditions = [ExperimentalDataFolder filesep 'ExperimentalConditions.csv'];
tableExperimentalConditions = loadTableExperimentalConditions(pathExperimentalConditions);


return % stop execution here

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Res2Var Comp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 DMSO/LatA
Res2Var_wFluo_multiZ_Comp('Results','20-05-20','R40','M1','DictyAx2_DMSO',1,...
    '19-05-20_DepthographM450',24,95,143,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-05-20','R40','M2','DictyAx2_DMSO',1,...
    '19-05-20_DepthographM450',24,95,143,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','22-05-20','R40','M1','DictyAx2_DMSO',1,...
    '19-05-20_DepthographM450',24,95,143,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','22-05-20','R40','M1','DictyAx2_LatA',2,...
    '19-05-20_DepthographM450',24,95,143,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','25-05-20','R40','M1','DictyAx2_DMSO',1,...
    '19-05-20_DepthographM450',24,95,143,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','25-05-20','R40','M1','DictyAx2_LatA',2,...
    '19-05-20_DepthographM450',24,95,143,AUTO,RawdataFolder,MatfileFolder)


Res2Var_wFluo_multiZ_Comp('Results','27-05-20','R40','M1','DictyAx2_DMSO',1,...
    '19-05-20_DepthographM450',24,95,143,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','05-06-20','R40','M1','DictyAx2_DMSO',1,...
    '19-05-20_DepthographM450',24,95,143,AUTO,RawdataFolder,MatfileFolder)

%% M270 Ax2 WT (Nishit)

Res2Var_wFluo_multiZ_Comp('Results','25-08-20','R90','M1','DictyAx2_M270',1,...
    '26-08-20_DepthographM270',24,95,143,AUTO,RawdataFolder,MatfileFolder)


Res2Var_wFluo_multiZ_Comp('Results','26-08-20','R90','M1','DictyAx2_M270',1,...
    '26-08-20_DepthographM270',24,95,143,AUTO,RawdataFolder,MatfileFolder)

%% M450 Ax2 WT + Mutants (DictyBase)

Res2Var_wFluo_multiZ_Comp('Results','23-09-20','R40','M1','DictyDB_WT',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','23-09-20','R40','M2','DictyDB_SevA',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','23-09-20','R40','M3','DictyDB_FimA',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)


Res2Var_wFluo_multiZ_Comp('Results','25-09-20','R40','M1','DictyDB_SevA',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','25-09-20','R40','M2','DictyDB_FimA',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','25-09-20','R40','M3','DictyDB_WT',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)



Res2Var_wFluo_multiZ_Comp('Results','06-10-20','R40','M1','DictyDB_abpC',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','06-10-20','R40','M2','DictyDB_abpA',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','06-10-20','R40','M3','DictyDB_WT',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','08-10-20','R40','M1','DictyDB_abpA',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','08-10-20','R40','M2','DictyDB_abpC',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','08-10-20','R40','M3','DictyDB_WT',1,...
    '19-05-20_DepthographM450',36,95,167,AUTO,RawdataFolder,MatfileFolder)

%% MultiTime Comp M450 Ax2 (Nishit)

Res2Var_wFluo_multiZ_Comp_R248('Results','13-02-20','R248','M1','DictyAx2-Comp',1:2,...
    '02-03-20_Depthograph100X',21,24,98,193,383,...
    AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp_R248('Results','20-02-20','R248','M1','DictyAx2-Comp',1,...
    '02-03-20_Depthograph100X',21,24,98,193,383,...
    AUTO,RawdataFolder,MatfileFolder)

%% MultiRate Comp M450 WTDB

% 20-1020
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_04s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_1s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,100,136,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_2s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_4s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_8s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_14s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,140,176,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_24s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,141,177,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_04s','M2','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_1s','M2','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,100,136,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_2s','M2','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_4s','M2','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_8s','M2','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_14s','M2','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,140,176,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_24s','M2','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,141,177,AUTO,RawdataFolder,MatfileFolder)


% 21-10-20
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_04s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_1s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,100,136,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_2s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_4s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_8s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_14s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,140,176,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_24s','M1','DictyDB_M450',1,...
    '22-10-20_DepthographM450',18,141,177,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_04s','M2','DictyDB_M450',2,...
    '22-10-20_DepthographM450',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_1s','M2','DictyDB_M450',2,...
    '22-10-20_DepthographM450',18,100,136,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_2s','M2','DictyDB_M450',2,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_4s','M2','DictyDB_M450',2,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_8s','M2','DictyDB_M450',2,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_14s','M2','DictyDB_M450',2,...
    '22-10-20_DepthographM450',18,140,176,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R50_24s','M2','DictyDB_M450',2,...
    '22-10-20_DepthographM450',18,141,177,AUTO,RawdataFolder,MatfileFolder)

% 22-10-20
Res2Var_wFluo_multiZ_Comp('Results','22-10-20','R50_04s','M1','DictyDB_M450',1:2,...
    '22-10-20_DepthographM450',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','22-10-20','R50_4s','M1','DictyDB_M450',1:2,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','22-10-20','R50_24s','M1','DictyDB_M450',1:2,...
    '22-10-20_DepthographM450',18,141,177,AUTO,RawdataFolder,MatfileFolder)

% 28-10-20
Res2Var_wFluo_multiZ_Comp_MultiRatePerCell('Results','28-10-20','R50-Multi','M1','DictyDB_M450-Multi',1:2,...
    '22-10-20_DepthographM450','RmpTimeDatas',AUTO,RawdataFolder,MatfileFolder)

% 17-11-20
Res2Var_wFluo_multiZ_Comp('Results','17-11-20','R50-04s','M2','DictyDB_M450',1,...
    '17-11-20_DepthographM450',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','17-11-20','R50-1s','M2','DictyDB_M450',1,...
    '17-11-20_DepthographM450',18,100,136,AUTO,RawdataFolder,MatfileFolder)
% low force
Res2Var_wFluo_multiZ_Comp('Results','17-11-20','R12-04s','M1','DictyDB_M45R12',1,...
    '17-11-20_DepthographM450',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','17-11-20','R12-4s','M1','DictyDB_M45R12',1,...
    '17-11-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)

% 18-11-20
Res2Var_wFluo_multiZ_Comp('Results','18-11-20','R50-8s','M1','DictyDB_M450',2,...
    '17-11-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)
% low force
Res2Var_wFluo_multiZ_Comp('Results','18-11-20','R12-2s','M1','DictyDB_M45R12',1,...
    '17-11-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)

%% MultiRate Comp M270 WTDB

% 20-1020
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R90_04s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R90_1s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,100,136,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R90_2s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R90_4s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R90_8s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R90_14s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,140,176,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R90_24s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,141,177,AUTO,RawdataFolder,MatfileFolder)

% ici le tag est R50 par erreur lors des manip sur labview, mais le champ
% va bien jusqu'a 90mt
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_04s','M2','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_1s','M2','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,100,136,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_2s','M2','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_4s','M2','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_8s','M2','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_14s','M2','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,140,176,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-10-20','R50_24s','M2','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,141,177,AUTO,RawdataFolder,MatfileFolder)


% 21-10-20
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_04s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_1s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,100,136,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_2s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_4s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_8s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_14s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,140,176,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_24s','M1','DictyDB_M270',2,...
    '27-10-20_DepthographM270',18,141,177,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_04s','M2','DictyDB_M270',1,...
    '27-10-20_DepthographM270',18,40,76,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_1s','M2','DictyDB_M270',1,...
    '27-10-20_DepthographM270',18,100,136,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_2s','M2','DictyDB_M270',1,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_4s','M2','DictyDB_M270',1,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_8s','M2','DictyDB_M270',1,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_14s','M2','DictyDB_M270',1,...
    '27-10-20_DepthographM270',18,140,176,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-10-20','R90_24s','M2','DictyDB_M270',1,...
    '27-10-20_DepthographM270',18,141,177,AUTO,RawdataFolder,MatfileFolder)

%% LongueManip comp M450/M270

Res2Var_wFluo_multiZ_Comp('Results','04-11-20','R50_2s','M1','DictyDB_M450',1:2,...
    '22-10-20_DepthographM450',18,133,169,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','04-11-20','R90_2s','M2','DictyDB_M270',1,...
    '27-10-20_DepthographM270',18,133,169,AUTO,RawdataFolder,MatfileFolder)

%% 3T3 Joseph - Internship aSFL nodrug vs. aSFL doxy

FLUO = false;
Res2Var_wFluo_multiZ_Comp('Results','20-08-04','R40','M1','3T3aSFL_BSA_nodrugs',1,...
    '20-02-20_Depthograph100x',tableExperimentalConditions,24,95,143,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-08-05','R40','M2','3T3aSFL_BSA_nodrugs',1,...
    '20-02-20_Depthograph100x',tableExperimentalConditions,24,95,143,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-08-07','R40','M2','3T3aSFL_BSA_nodrugs',1,...
    '20-02-20_Depthograph100x',tableExperimentalConditions,24,95,143,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-08-04','R40','M2','3T3aSFL_BSA_doxy',1,...
    '20-02-20_Depthograph100x',tableExperimentalConditions,24,95,143,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-08-05','R40','M1','3T3aSFL_BSA_doxy',1,...
    '20-02-20_Depthograph100x',tableExperimentalConditions,24,95,143,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','20-08-07','R40','M1','3T3aSFL_BSA_doxy',1,...
    '20-02-20_Depthograph100x',tableExperimentalConditions,24,95,143,FLUO,AUTO,RawdataFolder,MatfileFolder)

%% 3T3 Joseph - aSFL 3T3 on 20µm disc patterns of fibronectin 
% nodrug refers to the control, doxy to the activated cells expressing a membrane - cortex linker.

FLUO = false;
Res2Var_wFluo_multiZ_Comp('Results','21-01-18','R40_disc20um','M1','3T3aSFL_doxy',1,...
    '21-01-18_Deptho_M1',tableExperimentalConditions,18,133,169,FLUO,AUTO,RawdataFolder,MatfileFolder)
FLUO = true;
Res2Var_wFluo_multiZ_Comp('Results','21-01-18','R40_disc20um_wFluo','M1','3T3aSFL_doxy',1,...
    '21-01-18_Deptho_M1',tableExperimentalConditions,18,133,170,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-01-18','R40_disc20um_wFluo','M2','3T3aSFL_nodrug',1,...
    '21-01-18_Deptho_M2',tableExperimentalConditions,18,133,170,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-01-18','R40_disc20um_wFluo','M3','3T3aSFL_doxy',1,...
    '21-01-18_Deptho_M3',tableExperimentalConditions,18,133,170,FLUO,AUTO,RawdataFolder,MatfileFolder)

FLUO = true;
Res2Var_wFluo_multiZ_Comp('Results','21-01-21','R40_disc20um_wFluo','M1','3T3aSFL_nodrug',1,...
    '21-01-21_Deptho_M1',tableExperimentalConditions,18,133,170,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-01-21','R40_disc20um_wFluo','M2','3T3aSFL_doxy',1,...
    '21-01-21_Deptho_M2',tableExperimentalConditions,18,133,170,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-01-21','R40_disc20um_wFluo','M3','3T3aSFL_nodrug',1,...
    '21-01-21_Deptho_M3',tableExperimentalConditions,18,133,170,FLUO,AUTO,RawdataFolder,MatfileFolder)

% aSFL-6FP

FLUO = true;
Res2Var_wFluo_multiZ_Comp('Results','21-04-27','R40_disc20um_wFluo','M1','3T3aSFL-6FP_nodrug',1,...
    '21-04-27_Deptho_M1',tableExperimentalConditions,18,133,170,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-04-27','R40_disc20um_wFluo','M2','3T3aSFL-6FP_doxy',1,...
    '21-04-27_Deptho_M2',tableExperimentalConditions,18,133,170,FLUO,AUTO,RawdataFolder,MatfileFolder)

FLUO = false;
Res2Var_wFluo_multiZ_Comp('Results','21-04-28','R40_disc20um','M1','3T3aSFL-6FP_doxy',1,...
    '21-04-28_Deptho_M1',tableExperimentalConditions,18,133,169,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-04-28','R40_disc20um','M2','3T3aSFL-6FP_nodrug',1,...
    '21-04-28_Deptho_M2',tableExperimentalConditions,18,133,169,FLUO,AUTO,RawdataFolder,MatfileFolder)

% aSFL-A8

FLUO = false;
Res2Var_wFluo_multiZ_Comp('Results','21-06-16','R40_disc20um','M1','3T3aSFL-A8_nodrug',1,...
    '21-06-16_Deptho_M1',tableExperimentalConditions,18,133,169,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-06-16','R40_disc20um','M2','3T3aSFL-A8_doxy',1,...
    '21-06-16_Deptho_M2',tableExperimentalConditions,18,133,169,FLUO,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','21-06-17','R40_disc20um','M1','3T3aSFL-A8_doxy',1,...
    '21-06-17_Deptho_M1',tableExperimentalConditions,18,133,169,FLUO,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','21-06-17','R40_disc20um','M2','3T3aSFL-A8_nodrug',1,...
    '21-06-17_Deptho_M2',tableExperimentalConditions,18,133,169,FLUO,AUTO,RawdataFolder,MatfileFolder)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Var2Data Classic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 Ax2 DMSO/LatA (Nishit)

Var2Data_wFluo('20-05-20','R40cst',1.1,'M1','DictyAx2_DMSO',4504,MatfileFolder);
Var2Data_wFluo('20-05-20','R40cst',1.1,'M2','DictyAx2_DMSO',4504,MatfileFolder);

Var2Data_wFluo('22-05-20','R40cst',1.1,'M1','DictyAx2_DMSO',4504,MatfileFolder);
Var2Data_wFluo('22-05-20','R40cst',1.1,'M1','DictyAx2_LatA',4504,MatfileFolder);

Var2Data_wFluo('25-05-20','R40cst',1.1,'M1','DictyAx2_DMSO',4504,MatfileFolder);
Var2Data_wFluo('25-05-20','R40cst',1.1,'M1','DictyAx2_LatA',4504,MatfileFolder);

Var2Data_wFluo('27-05-20','R40cst',1.1,'M1','DictyAx2_DMSO',4504,MatfileFolder);

RemoveFromData(MatfileFolder,{'25-05-20_M1_P2_C9_R40','20-05-20_M1_P1_C6_R40'},'DictyAx2')

%% M270 Ax2 WT (Nishit)

Var2Data_wFluo('25-08-20','R90cst',0.82,'M1','DictyAx2_M270',2691,MatfileFolder);

Var2Data_wFluo('26-08-20','R90cst',0.82,'M1','DictyAx2_M270',2691,MatfileFolder);

% removing
RemoveFromData(MatfileFolder,{'25-08-20_M1_P1_C7_R90','25-08-20_M1_P1_C7_R90'},'DictyAx2')

%% M450 Ax2 WT + Mutants (DictyBase)

Var2Data_wFluo('23-09-20','R40cst',1.1,'M1','DictyDB_WT',4504,MatfileFolder);
Var2Data_wFluo('23-09-20','R40cst',1.1,'M2','DictyDB_SevA',4504,MatfileFolder);
Var2Data_wFluo('23-09-20','R40cst',1.1,'M3','DictyDB_FimA',4504,MatfileFolder);

Var2Data_wFluo('25-09-20','R40cst',1.1,'M1','DictyDB_SevA',4504,MatfileFolder);
Var2Data_wFluo('25-09-20','R40cst',1.1,'M2','DictyDB_FimA',4504,MatfileFolder);
Var2Data_wFluo('25-09-20','R40cst',1.1,'M3','DictyDB_WT',4504,MatfileFolder);

Var2Data_wFluo('06-10-20','R40cst',1.1,'M1','DictyDB_abpC',4504,MatfileFolder);
Var2Data_wFluo('06-10-20','R40cst',1.1,'M2','DictyDB_abpA',4504,MatfileFolder);
Var2Data_wFluo('06-10-20','R40cst',1.1,'M3','DictyDB_WT',4504,MatfileFolder);

Var2Data_wFluo('08-10-20','R40cst',1.1,'M1','DictyDB_abpA',4504,MatfileFolder);
Var2Data_wFluo('08-10-20','R40cst',1.1,'M2','DictyDB_abpC',4504,MatfileFolder);
Var2Data_wFluo('08-10-20','R40cst',1.1,'M3','DictyDB_WT',4504,MatfileFolder);

% removing
RemoveFromData(MatfileFolder,{'06-10-20_M3_P1_C7_R40','06-10-20_M1_P1_C1_R40','06-10-20_M1_P1_C2_R40','06-10-20_M1_P1_C4_R40',...
    '08-10-20_M3_P1_C1_R40','08-10-20_M3_P1_C4_R40','08-10-20_M3_P1_C6_R40','08-10-20_M2_P1_C2_R40',...
    '23-09-20_M1_P1_C9_R40','23-09-20_M2_P1_C3_R40','23-09-20_M2_P1_C11_R40',...
    '25-09-20_M1_P1_C12_R40','25-09-20_M1_P1_C8_R40','25-09-20_M3_P1_C1_R40',},'DictyDB')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data2Peak %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 Ax2 DMSO/LatA (Nishit)

Data2Peaks('DictyAx2_DMSO',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DictyAx2_LatA',PLOTDP,MatfileFolder,FigureFolder);

%% M270 Ax2 WT (Nishit)

Data2Peaks('DictyAx2_M270',PLOTDP,MatfileFolder,FigureFolder);

%% M450 Ax2 WT + Mutants (DictyBase)

Data2Peaks('DictyDB_WT',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DictyDB_SevA',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DictyDB_FimA',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DictyDB_abpA',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DictyDB_abpC',PLOTDP,MatfileFolder,FigureFolder);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Peaks2Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 Ax2 DMSO/LatA (Nishit)

Peaks2Stats('DictyAx2_DMSO',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
Peaks2Stats('DictyAx2_LatA',PLOTPS,alignlvl,MatfileFolder,FigureFolder);

%% M270 Ax2 WT (Nishit)
Peaks2Stats('DictyAx2_M270',PLOTPS,alignlvl,MatfileFolder,FigureFolder);

%% M450 Ax2 WT + Mutants (DictyBase)
Peaks2Stats('DictyDB_WT',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
Peaks2Stats('DictyDB_FimA',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
Peaks2Stats('DictyDB_SevA',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
Peaks2Stats('DictyDB_abpA',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
Peaks2Stats('DictyDB_abpC',PLOTPS,alignlvl,MatfileFolder,FigureFolder);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stats2Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

Stats2Plots({'DictyDB_WT','DictyDB_FimA'},{Cwt,Cfim},'DictyMutantsFimComp',CORT,DIST,PEAKS,ACTI,TIME,alignlvl,MatfileFolder,FigureFolder);

Stats2Plots({'DictyDB_WT','DictyDB_abpA','DictyDB_abpC','DictyDB_FimA'},{Cwt,CaA,CaC,Cfim},'DictyMutantsComp',CORT,DIST,PEAKS,ACTI,TIME,alignlvl,MatfileFolder,FigureFolder);

Stats2Plots({'DictyAx2_DMSO','DictyAx2_M270'},{Cdm,Cwt27},'DictyBeadsSize',CORT,DIST,PEAKS,ACTI,TIME,alignlvl,MatfileFolder,FigureFolder);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Var2Data Comp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 Ax2 DMSO/LatA (Nishit)

Var2Data_Comp('20-05-20','R40',1.1,'M1','DictyAx2_DMSO',95,143,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-05-20','R40',1.1,'M2','DictyAx2_DMSO',95,143,'1s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('22-05-20','R40',1.1,'M1','DictyAx2_DMSO',95,143,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('22-05-20','R40',1.1,'M1','DictyAx2_LatA',95,143,'1s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('25-05-20','R40',1.1,'M1','DictyAx2_DMSO',95,143,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('25-05-20','R40',1.1,'M1','DictyAx2_LatA',95,143,'1s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('27-05-20','R40',1.1,'M1','DictyAx2_DMSO',95,143,'1s',4504,MatfileFolder,FigureFolder);

%% M270 Ax2 WT (Nishit)

Var2Data_Comp('25-08-20','R90',0.82,'M1','DictyAx2_M270',95,143,'1s',2691,MatfileFolder,FigureFolder);

Var2Data_Comp('26-08-20','R90',0.82,'M1','DictyAx2_M270',95,143,'1s',2691,MatfileFolder,FigureFolder);

%% M450 Ax2 WT + Mutants (DictyBase)

Var2Data_Comp('23-09-20','R40',1.1,'M1','DictyDB_WT',95,167,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('23-09-20','R40',1.1,'M2','DictyDB_SevA',95,167,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('23-09-20','R40',1.1,'M3','DictyDB_FimA',95,167,'1s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('25-09-20','R40',1.1,'M1','DictyDB_SevA',95,167,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('25-09-20','R40',1.1,'M2','DictyDB_FimA',95,167,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('25-09-20','R40',1.1,'M3','DictyDB_WT',95,167,'1s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('06-10-20','R40',1.1,'M1','DictyDB_abpC',95,167,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('06-10-20','R40',1.1,'M2','DictyDB_abpA',95,167,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('06-10-20','R40',1.1,'M3','DictyDB_WT',95,167,'1s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('08-10-20','R40',1.1,'M1','DictyDB_abpA',95,167,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('08-10-20','R40',1.1,'M2','DictyDB_abpC',95,167,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('08-10-20','R40',1.1,'M3','DictyDB_WT',95,167,'1s',4504,MatfileFolder,FigureFolder);

%% MultiTime Comp M450 Ax2 (Nishit)

Var2Data_Comp_R248('13-02-20','R248',1.1,'M1','DictyAx2-Comp',21,24,98,193,383,4438,MatfileFolder,FigureFolder);
Var2Data_Comp_R248('20-02-20','R248',1.1,'M1','DictyAx2-Comp',21,24,98,193,383,4438,MatfileFolder,FigureFolder);

%% MultiRate Comp M450 WTDB

%20-10-20
Var2Data_Comp('20-10-20','R50_04s',0.85,'M1','DictyDB_M450',40,76,'02s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_1s',0.85,'M1','DictyDB_M450',100,136,'05s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_2s',0.85,'M1','DictyDB_M450',133,169,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_4s',0.85,'M1','DictyDB_M450',133,169,'2s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_8s',0.85,'M1','DictyDB_M450',133,169,'4s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_14s',0.85,'M1','DictyDB_M450',140,176,'7s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_24s',0.85,'M1','DictyDB_M450',141,177,'12s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('20-10-20','R50_04s',0.85,'M2','DictyDB_M450',40,76,'02s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_1s',0.85,'M2','DictyDB_M450',100,136,'05s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_2s',0.85,'M2','DictyDB_M450',133,169,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_4s',0.85,'M2','DictyDB_M450',133,169,'2s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_8s',0.85,'M2','DictyDB_M450',133,169,'4s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_14s',0.85,'M2','DictyDB_M450',140,176,'7s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_24s',0.85,'M2','DictyDB_M450',141,177,'12s',4504,MatfileFolder,FigureFolder);

% 21-10-20
Var2Data_Comp('21-10-20','R50_04s',0.85,'M1','DictyDB_M450',40,76,'02s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_1s',0.85,'M1','DictyDB_M450',100,136,'05s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_2s',0.85,'M1','DictyDB_M450',133,169,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_4s',0.85,'M1','DictyDB_M450',133,169,'2s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_8s',0.85,'M1','DictyDB_M450',133,169,'4s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_14s',0.85,'M1','DictyDB_M450',140,176,'7s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_24s',0.85,'M1','DictyDB_M450',141,177,'12s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('21-10-20','R50_04s',0.85,'M2','DictyDB_M450',40,76,'02s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_1s',0.85,'M2','DictyDB_M450',100,136,'05s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_2s',0.85,'M2','DictyDB_M450',133,169,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_4s',0.85,'M2','DictyDB_M450',133,169,'2s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_8s',0.85,'M2','DictyDB_M450',133,169,'4s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_14s',0.85,'M2','DictyDB_M450',140,176,'7s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R50_24s',0.85,'M2','DictyDB_M450',141,177,'12s',4504,MatfileFolder,FigureFolder);

% 22-10-20
Var2Data_Comp('22-10-20','R50_04s',0.85,'M1','DictyDB_M450',40,76,'02s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('22-10-20','R50_4s',0.85,'M1','DictyDB_M450',133,169,'2s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('22-10-20','R50_24s',0.85,'M1','DictyDB_M450',141,177,'12s',4504,MatfileFolder,FigureFolder);


% 28-10-20
Var2Data_Comp_MultiRatePerCell('28-10-20','R50-Multi',0.85,'M1','DictyDB_M450-Multi',...
    'RmpTimeDatas',4504,RawdataFolder,MatfileFolder,FigureFolder);

% 17-11-20
Var2Data_Comp('17-11-20','R50-04s',0.85,'M2','DictyDB_M450',40,76,'02s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('17-11-20','R50-1s',0.85,'M2','DictyDB_M450',100,136,'05s',4504,MatfileFolder,FigureFolder);
% low force
Var2Data_Comp('17-11-20','R12-04s',0.85,'M1','DictyDB_M45R12',40,76,'02s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('17-11-20','R12-4s',0.85,'M1','DictyDB_M45R12',133,169,'2s',4504,MatfileFolder,FigureFolder);

% 18-11-20
Var2Data_Comp('18-11-20','R50-8s',0.85,'M1','DictyDB_M450',133,169,'4s',4504,MatfileFolder,FigureFolder);
% low force
Var2Data_Comp('18-11-20','R12-2s',0.85,'M1','DictyDB_M45R12',133,169,'1s',4504,MatfileFolder,FigureFolder);

%% MultiRate Comp M270 WTDB

%20-10-20
Var2Data_Comp('20-10-20','R90_04s',0.85,'M1','DictyDB_M270',40,76,'02s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R90_1s',0.85,'M1','DictyDB_M270',100,136,'05s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R90_2s',0.85,'M1','DictyDB_M270',133,169,'1s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R90_4s',0.85,'M1','DictyDB_M270',133,169,'2s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R90_8s',0.85,'M1','DictyDB_M270',133,169,'4s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R90_14s',0.85,'M1','DictyDB_M270',140,176,'7s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R90_24s',0.85,'M1','DictyDB_M270',141,177,'12s',2691,MatfileFolder,FigureFolder);

Var2Data_Comp('20-10-20','R50_04s',0.85,'M2','DictyDB_M270',40,76,'02s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_1s',0.85,'M2','DictyDB_M270',100,136,'05s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_2s',0.85,'M2','DictyDB_M270',133,169,'1s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_4s',0.85,'M2','DictyDB_M270',133,169,'2s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_8s',0.85,'M2','DictyDB_M270',133,169,'4s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_14s',0.85,'M2','DictyDB_M270',140,176,'7s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('20-10-20','R50_24s',0.85,'M2','DictyDB_M270',141,177,'12s',2691,MatfileFolder,FigureFolder);

% 21-10-20
Var2Data_Comp('21-10-20','R90_04s',0.85,'M1','DictyDB_M270',40,76,'02s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_1s',0.85,'M1','DictyDB_M270',100,136,'05s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_2s',0.85,'M1','DictyDB_M270',133,169,'1s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_4s',0.85,'M1','DictyDB_M270',133,169,'2s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_8s',0.85,'M1','DictyDB_M270',133,169,'4s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_14s',0.85,'M1','DictyDB_M270',140,176,'7s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_24s',0.85,'M1','DictyDB_M270',141,177,'12s',2691,MatfileFolder,FigureFolder);

Var2Data_Comp('21-10-20','R90_04s',0.85,'M2','DictyDB_M270',40,76,'02s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_1s',0.85,'M2','DictyDB_M270',100,136,'05s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_2s',0.85,'M2','DictyDB_M270',133,169,'1s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_4s',0.85,'M2','DictyDB_M270',133,169,'2s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_8s',0.85,'M2','DictyDB_M270',133,169,'4s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_14s',0.85,'M2','DictyDB_M270',140,176,'7s',2691,MatfileFolder,FigureFolder);
Var2Data_Comp('21-10-20','R90_24s',0.85,'M2','DictyDB_M270',141,177,'12s',2691,MatfileFolder,FigureFolder);

%% LongueManip comp M450/M270 WTDB
Var2Data_Comp('04-11-20','R50_2s',0.85,'M1','DictyDB_M450',133,169,'1s',4504,MatfileFolder,FigureFolder);

Var2Data_Comp('04-11-20','R90_2s',0.85,'M2','DictyDB_M270',133,169,'1s',2691,MatfileFolder,FigureFolder);

%% 3T3 Joseph (Internship aSFL nodrug vs. aSFL doxy)
PLOTV2D = 1;
Var2Data_Comp('20-08-04','R40','M1',tableExperimentalConditions,'3T3aSFL_BSA_nodrugs',95,143,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('20-08-05','R40','M2',tableExperimentalConditions,'3T3aSFL_BSA_nodrugs',95,143,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('20-08-07','R40','M2',tableExperimentalConditions,'3T3aSFL_BSA_nodrugs',95,143,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);

Var2Data_Comp('20-08-04','R40','M2',tableExperimentalConditions,'3T3aSFL_BSA_doxy',95,143,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('20-08-05','R40','M1',tableExperimentalConditions,'3T3aSFL_BSA_doxy',95,143,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('20-08-07','R40','M1',tableExperimentalConditions,'3T3aSFL_BSA_doxy',95,143,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);

%% 3T3 Joseph - aSFL 3T3 on 20µm disc patterns of fibronectin

PLOTV2D = 1;

Var2Data_Comp('21-01-18','R40_disc20um','M1',tableExperimentalConditions,'3T3aSFL_doxy',133,169,'1s',PLOTV2D,MatfileFolder,FigureFolder, ExportDataFolder); %,ExportDataFolder);

Var2Data_Comp('21-01-18','R40_disc20um_wFluo','M1',tableExperimentalConditions,'3T3aSFL_doxy',133,170,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('21-01-18','R40_disc20um_wFluo','M2',tableExperimentalConditions,'3T3aSFL_nodrug',133,170,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('21-01-18','R40_disc20um_wFluo','M3',tableExperimentalConditions,'3T3aSFL_doxy',133,170,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);

PLOTV2D = 1;
Var2Data_Comp('21-01-21','R40_disc20um_wFluo','M1',tableExperimentalConditions,'3T3aSFL_nodrug',133,170,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('21-01-21','R40_disc20um_wFluo','M2',tableExperimentalConditions,'3T3aSFL_doxy',133,170,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('21-01-21','R40_disc20um_wFluo','M3',tableExperimentalConditions,'3T3aSFL_nodrug',133,170,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);

% aSFL-6FP

PLOTV2D = 1;
Var2Data_Comp('21-04-27','R40_disc20um_wFluo','M1',tableExperimentalConditions,'3T3aSFL-6FP_nodrug',133,170,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('21-04-27','R40_disc20um_wFluo','M2',tableExperimentalConditions,'3T3aSFL-6FP_doxy',133,170,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);

PLOTV2D = 1;
Var2Data_Comp('21-04-28','R40_disc20um','M1',tableExperimentalConditions,'3T3aSFL-6FP_doxy',133,169,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('21-04-28','R40_disc20um','M2',tableExperimentalConditions,'3T3aSFL-6FP_nodrug',133,169,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);

% aSFL-A8

PLOTV2D = 1;
Var2Data_Comp('21-06-16','R40_disc20um','M1',tableExperimentalConditions,'3T3aSFL-A8_nodrug',133,169,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('21-06-16','R40_disc20um','M2',tableExperimentalConditions,'3T3aSFL-A8_doxy',133,169,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);

PLOTV2D = 1;
Var2Data_Comp('21-06-17','R40_disc20um','M1',tableExperimentalConditions,'3T3aSFL-A8_doxy',133,169,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);
Var2Data_Comp('21-06-17','R40_disc20um','M2',tableExperimentalConditions,'3T3aSFL-A8_nodrug',133,169,'1s',PLOTV2D,MatfileFolder,FigureFolder,ExportDataFolder);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data2Meca w/ Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%



%% MultiTime Comp M450 Ax2 (Nishit)

Data2MecaLoop_BigTable('13-02-20','M1','DictyAx2-Comp','R248',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('20-02-20','M1','DictyAx2-Comp','R248',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%% M450 Ax2 DMSO/LatA (Nishit)

Data2MecaLoop_BigTable('20-05-20','M1','DictyAx2_DMSO','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-05-20','M2','DictyAx2_DMSO','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('22-05-20','M1','DictyAx2_DMSO','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('22-05-20','M1','DictyAx2_LatA','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('25-05-20','M1','DictyAx2_DMSO','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('25-05-20','M1','DictyAx2_LatA','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('27-05-20','M1','DictyAx2_DMSO','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%% M450 Ax2 WT + Mutants (DictyBase)

Data2MecaLoop_BigTable('23-09-20','M1','DictyDB_WT','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('23-09-20','M2','DictyDB_SevA','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('23-09-20','M3','DictyDB_FimA','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('25-09-20','M1','DictyDB_SevA','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('25-09-20','M2','DictyDB_FimA','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('25-09-20','M3','DictyDB_WT','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('06-10-20','M1','DictyDB_abpC','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('06-10-20','M2','DictyDB_abpA','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('06-10-20','M3','DictyDB_WT','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('08-10-20','M1','DictyDB_abpA','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('08-10-20','M2','DictyDB_abpC','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('08-10-20','M3','DictyDB_WT','R40',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%% M270 Ax2 WT (Nishit)

Data2MecaLoop_BigTable('25-08-20','M1','DictyAx2_M270','R90',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('26-08-20','M1','DictyAx2_M270','R90',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%% MultiRate Comp M450 WTDB

% 20-10-20
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M450','R50_04s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M450','R50_1s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M450','R50_2s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M450','R50_4s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M450','R50_8s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M450','R50_14s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M450','R50_24s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M450','R50_04s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M450','R50_1s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M450','R50_2s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M450','R50_4s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M450','R50_8s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M450','R50_14s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M450','R50_24s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%21-10-20
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M450','R50_04s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M450','R50_1s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M450','R50_2s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M450','R50_4s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M450','R50_8s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M450','R50_14s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M450','R50_24s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

% Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M450','R50_04s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
% Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M450','R50_1s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
% Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M450','R50_2s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
% Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M450','R50_4s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
% Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M450','R50_8s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
% Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M450','R50_14s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
% Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M450','R50_24s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%22-10-20
Data2MecaLoop_BigTable('22-10-20','M1','DictyDB_M450','R50_04s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('22-10-20','M1','DictyDB_M450','R50_4s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('22-10-20','M1','DictyDB_M450','R50_24s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%28-10-20
Data2MecaLoop_BigTable('28-10-20','M1','DictyDB_M450-Multi','R50-Multi',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

% 17-11-20
Data2MecaLoop_BigTable('17-11-20','M2','DictyDB_M450','R50-04s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('17-11-20','M2','DictyDB_M450','R50-1s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
% low force
Data2MecaLoop_BigTable('17-11-20','M1','DictyDB_M45R12','R12-04s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('17-11-20','M1','DictyDB_M45R12','R12-4s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

% 18-11-20
Data2MecaLoop_BigTable('18-11-20','M1','DictyDB_M450','R50-8s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
% low force
Data2MecaLoop_BigTable('18-11-20','M1','DictyDB_M45R12','R12-2s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%% MultiRate Comp M270 WTDB

% 20-10-20
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M270','R90_04s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M270','R90_1s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M270','R90_2s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M270','R90_4s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M270','R90_8s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M270','R90_14s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M1','DictyDB_M270','R90_24s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M270','R50_04s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M270','R50_1s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M270','R50_2s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M270','R50_4s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M270','R50_8s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M270','R50_14s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-10-20','M2','DictyDB_M270','R50_24s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%21-10-20
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M270','R90_04s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M270','R90_1s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M270','R90_2s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M270','R90_4s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M270','R90_8s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M270','R90_14s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M1','DictyDB_M270','R90_24s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M270','R90_04s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M270','R90_1s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M270','R90_2s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M270','R90_4s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M270','R90_8s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M270','R90_14s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-10-20','M2','DictyDB_M270','R90_24s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%% LongueManip comp M450/M270
Data2MecaLoop_BigTable('04-11-20','M1','DictyDB_M450','R50_2s',4504,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('04-11-20','M2','DictyDB_M270','R90_2s',2691,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

%% 3T3 Joseph 

% (Internship aSFL nodrug vs. aSFL doxy)
PLOTD2M = 0;
Data2MecaLoop_BigTable('20-08-04','M1',tableExperimentalConditions,'3T3aSFL_BSA_nodrugs','R40',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-08-05','M2',tableExperimentalConditions,'3T3aSFL_BSA_nodrugs','R40',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-08-07','M2',tableExperimentalConditions,'3T3aSFL_BSA_nodrugs','R40',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('20-08-04','M2',tableExperimentalConditions,'3T3aSFL_BSA_doxy','R40',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-08-05','M1',tableExperimentalConditions,'3T3aSFL_BSA_doxy','R40',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('20-08-07','M1',tableExperimentalConditions,'3T3aSFL_BSA_doxy','R40',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);

%% 3T3 Joseph - aSFL 3T3 on 20µm disc patterns of fibronectin
Data2MecaLoop_BigTable('21-01-18','M1',tableExperimentalConditions,'3T3aSFL_doxy','R40_disc20um',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-01-18','M1',tableExperimentalConditions,'3T3aSFL_doxy','R40_disc20um_wFluo',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-01-18','M2',tableExperimentalConditions,'3T3aSFL_nodrug','R40_disc20um_wFluo',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-01-18','M3',tableExperimentalConditions,'3T3aSFL_doxy','R40_disc20um_wFluo',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('21-01-21','M1',tableExperimentalConditions,'3T3aSFL_nodrug','R40_disc20um_wFluo',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-01-21','M2',tableExperimentalConditions,'3T3aSFL_doxy','R40_disc20um_wFluo',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-01-21','M3',tableExperimentalConditions,'3T3aSFL_nodrug','R40_disc20um_wFluo',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);

PLOTD2M = 0;
Data2MecaLoop_BigTable('21-04-27','M1',tableExperimentalConditions,'3T3aSFL-6FP_nodrug','R40_disc20um_wFluo',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-04-27','M2',tableExperimentalConditions,'3T3aSFL-6FP_doxy','R40_disc20um_wFluo',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('21-04-28','M1',tableExperimentalConditions,'3T3aSFL-6FP_doxy','R40_disc20um',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('21-04-28','M2',tableExperimentalConditions,'3T3aSFL-6FP_nodrug','R40_disc20um',MatfileFolder,FigureFolder,ExportDataFolder,PLOTD2M,VERBOSED2M);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Meca2Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% %%%%%%%%% Several R2 available : 0.7, 0.9 & 0.95

fitparams = 'Strain100-R20.9';


%% 3T3 Joseph 

% (Internship aSFL nodrug vs. aSFL doxy)

Conditions3T3aSFL = {'ExpType','3T3aSFL_BSA_nodrugs','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL_BSA_doxy','TpsComp','1s','FitParams',fitparams;...
    };

Meca2Plot_BigTable(MatfileFolder,FigureFolder,Conditions3T3aSFL, ...
    {CaSFLnodrug,CaSFLdoxy},{'o','o'},...
    {'aSFL nodrug','aSFL doxy'},'3T3_BSA_noDrugVsDoxy',1)

%% M450 mutant
ConditionsM450Multi = {'ExpType','DictyDB_WT','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_FimA','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_abpA','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_abpC','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_SevA','TpsComp','1s','FitParams',fitparams;...
    };

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450Multi, ...
{Cwt,Cfim,CaA,CaC,Csev},{'o','o','o','o','o'},...
{'WT','FimA-','abpA-','abpC-','SevA-'},'M450Mutants',1)

%% 3T3vsDicty Joseph & Val

% (Internship aSFL nodrug vs. aSFL doxy)

Conditions3T3vsDictys = {'ExpType','3T3aSFL_BSA_nodrugs','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_WT','TpsComp','1s','FitParams',fitparams;...
    };

Meca2Plot_BigTable(MatfileFolder,FigureFolder,Conditions3T3vsDictys, ...
    {CaSFLnodrug,Cwt},{'o','o'},...
    {'3T3 aSFL non activated','Dictys Ax2 WT'},'3T3_BSA_noDrugVsDictys',1)

%% 3T3 patterns aSFL Joseph

Conditions3T3aSFL = {'ExpType','3T3aSFL_nodrug','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL_doxy','TpsComp','1s','FitParams',fitparams;...
    };

Meca2Plot_BigTable(MatfileFolder,FigureFolder,Conditions3T3aSFL, ...
    {CaSFLnodrug,CaSFLdoxy},{'o','o'},...
    {'aSFL nodrug','aSFL doxy'},'3T3patterns_noDrugVsDoxy',1)

%% 3T3 patterns aSFL-6FP  Joseph

Conditions3T3aSFL6FP = {'ExpType','3T3aSFL-6FP_nodrug','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL-6FP_doxy','TpsComp','1s','FitParams',fitparams;...
    };

Meca2Plot_BigTable(MatfileFolder,FigureFolder,Conditions3T3aSFL6FP, ...
    {CaSFLnodrug,CaSFLdoxy},{'>','>'},...
    {'aSFL-6FP nodrug','aSFL-6FP doxy'},'3T3patterns_aSFL6FP_noDrugVsDoxy',1)


%% 3T3 patterns aSFL vs aSFL-6FP  Joseph

Conditions3T3aSFLVs6FP = {'ExpType','3T3aSFL_nodrug','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL_doxy','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL-6FP_nodrug','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL-6FP_doxy','TpsComp','1s','FitParams',fitparams;...
    };

Meca2Plot_BigTable(MatfileFolder,FigureFolder,Conditions3T3aSFLVs6FP, ...
    {CaSFLnodrug,CaSFLdoxy,CaSFLnodrug,CaSFLdoxy},{'^','^','>','>'},...
    {'aSFL nodrug','aSFL doxy', 'aSFL-6FP nodrug','aSFL-6FP doxy'},'3T3patterns_aSFLvsaSFL6FP_noDrugVsDoxy',1)

%% 3T3 patterns vs supsension Joseph

Conditions3T3aSFLsupsVsAdh = {'ExpType','3T3aSFL_BSA_nodrugs','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL_BSA_doxy','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL_nodrug','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL_doxy','TpsComp','1s','FitParams',fitparams;...
    };

Meca2Plot_BigTable(MatfileFolder,FigureFolder,Conditions3T3aSFLsupsVsAdh, ...
    {CaSFLnodrug,CaSFLdoxy,CaSFLnodrug,CaSFLdoxy},{'o','o','^','^'},...
    {'aSFL nodrug BSA','aSFL doxy BSA', 'aSFL nodrug fibro','aSFL doxy fibro'},'3T3patterns_supsVsAdh_noDrugVsDoxy',1)

%% Dicty NS Multitime
ConditionsM450NSMulti = {'ExpType','DictyAx2-Comp','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyAx2-Comp','TpsComp','2s','FitParams',fitparams;...
    'ExpType','DictyAx2-Comp','TpsComp','4s','FitParams',fitparams;...
    };
Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450NSMulti, ...
{C1,C2,C4},{'o','o','o'},...
{'M450-NS-1s','M450-NS-2s','M450-NS-4s'},'M450MultiRate',1)

%% M450 multirate
ConditionsM450Multi = {'ExpType','DictyDB_M450','TpsComp','02s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','05s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','2s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','4s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','7s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','12s','FitParams',fitparams;};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450Multi, ...
{C02_27,C05_27,C1,C2,C4,C7,C12},{M02,M05,M1,M2,M4,M7,M12},...
{'M450-02s','M450-05s','M450-1s','M450-2s','M450-4s','M450-7s','M450-12s'},'M450MultiRateR209',1,...
'Exclude',{'ExpDay','04-11-20'},'Exclude',{'ExpDay','17-11-20'})

    % M450 multirate 041120
ConditionsM450041120= {'ExpType','DictyDB_M450','TpsComp','1s','FitParams',fitparams,...
    'ExpDay','04-11-20';
    }

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450041120, ...
{C1},{M02},...
{'M450-1s'},'M450_1s_04-11-20',1)


    % M450 multirate 171120
ConditionsM450171120= {'ExpType','DictyDB_M450','TpsComp','02s','FitParams',fitparams,...
    'ExpDay','17-11-20';
    'ExpType','DictyDB_M450','TpsComp','05s','FitParams',fitparams,...
    'ExpDay','17-11-20';
    };

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450171120, ...
{C02_27,C05_27},{M02,M05},...
{'M450-02s','M450-05s'},'M450_17-11-20',1)

% low force  17-11
ConditionsM45R12_1711= {'ExpType','DictyDB_M450','TpsComp','02s','FitParams',fitparams,'ExpDay','17-11-20';
    'ExpType','DictyDB_M45R12','TpsComp','02s','FitParams',fitparams,'ExpDay','17-11-20';
    'ExpType','DictyDB_M45R12','TpsComp','2s','FitParams',fitparams,'ExpDay','17-11-20';
    }

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM45R12_1711, ...
{C02_27,C12_27,C4},{M02,M02,M4},...
{'M45R50-02s','M45R12-02s','M4R12-2s'},'M450_lowforce_17-11-20',1)

% low force  18-11
ConditionsM45R12_1811= {'ExpType','DictyDB_M450','TpsComp','4s','FitParams',fitparams,'ExpDay','18-11-20';
    'ExpType','DictyDB_M45R12','TpsComp','1s','FitParams',fitparams,'ExpDay','18-11-20';
    }

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM45R12_1811, ...
{C02_27,C12_27},{M02,M02},...
{'M45R50-4s','M45R12-1s'},'M450_lowforce_18-11-20',1)

% only low force
ConditionsM45R12= {'ExpType','DictyDB_M45R12','TpsComp','02s','FitParams',fitparams,'ExpDay','17-11-20'; ...
    'ExpType','DictyDB_M45R12','TpsComp','1s','FitParams',fitparams,'ExpDay','18-11-20'; ...
    'ExpType','DictyDB_M45R12','TpsComp','2s','FitParams',fitparams,'ExpDay','17-11-20'; ...
    }
Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM45R12, ...
{C02,C4,C12},{M02,M1,M2},...
{'M45R12-02s','M45R12-1s','M45R12-2s'},'M450_lowforce',1)



% M450multirate per cell
ConditionsM450MultiCell = {'ExpType','DictyDB_M450-Multi','TpsComp','02s','FitParams',fitparams;...
    'ExpType','DictyDB_M450-Multi','TpsComp','05s','FitParams',fitparams;...
    'ExpType','DictyDB_M450-Multi','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_M450-Multi','TpsComp','2s','FitParams',fitparams;...
    'ExpType','DictyDB_M450-Multi','TpsComp','4s','FitParams',fitparams;...
    'ExpType','DictyDB_M450-Multi','TpsComp','7s','FitParams',fitparams;...
    'ExpType','DictyDB_M450-Multi','TpsComp','12s','FitParams',fitparams;};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450MultiCell, ...
{C02_27,C05_27,C1,C2,C4,C7,C12},{M02,M05,M1,M2,M4,M7,M12},...
{'M450-02s','M450-05s','M450-1s','M450-2s','M450-4s','M450-7s','M450-12s'},'M450MultiRatePerCellR209',1)

%% M270 multirate

ConditionsM270Multi = {'ExpType','DictyDB_M270','TpsComp','02s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','05s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','2s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','4s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','7s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','12s','FitParams',fitparams;};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM270Multi, ...
{C02_27,C05_27,C1_27,C2_27,C4_27,C7_27,C12_27},{M02,M05,M1,M2,M4,M7,M12},...
{'M270-02s','M270-05s','M270-1s','M270-2s','M270-4s','M270-7s','M270-12s'},'M270DMultiRate',1,...
'Exclude',{'ExpDay','04-11-20'})

%% M450270 multirate
ConditionsM450270Multi = {'ExpType','DictyDB_M450','TpsComp','02s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','02s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','05s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','05s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','2s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','2s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','4s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','4s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','7s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','7s','FitParams',fitparams;...
    'ExpType','DictyDB_M450','TpsComp','12s','FitParams',fitparams;
    'ExpType','DictyDB_M270','TpsComp','12s','FitParams',fitparams;};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450270Multi, ...
{C02,C02_27,C05,C05_27,C1,C1_27,C2,C2_27,C4,C4_27,C7,C7_27,C12,C12_27},...
{M02,M02,M05,M05,M1,M1,M2,M2,M4,M4,M7,M7,M12,M12},...
{'M450-02s','M270-02s','M450-05s','M270-05s','M450-1s','M270-1s','M450-2s',...
'M270-2s','M450-4s','M270-4s','M450-7s','M270-7s','M450-12s','M270-12s'},'M450270DMultiRate',0,...
'Exclude',{'ExpDay','04-11-20'})

%% M50/M270 original
ConditionsM450270_1s = {'ExpType','DictyDB_M450','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_WT','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_M270','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyAx2_DMSO','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyAx2_M270','TpsComp','1s','FitParams',fitparams};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450270_1s, ...
{C1,C1,C1_27,Cwt,Cwt27},{'o','o','o','o','o'},...
{'DB-M450-1s','DB-WT-1s','DB-M270-1s','NS-M450-1s','NS-M270-1s'},'M45-27_1s',1,...
'Exclude',{'ExpDay','04-11-20'})



 % only NS
ConditionsM450270_1s_NS = {'ExpType','DictyAx2_DMSO','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyAx2_M270','TpsComp','1s','FitParams',fitparams};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450270_1s_NS, ...
{Cwt,Cwt27},{'o','o'},...
{'M450-1s','M270-1s'},'M45-27_1s_NS',1)




    % M450 04-11-20
ConditionsM450041120= {'ExpType','DictyDB_M450','TpsComp','1s','FitParams',fitparams,...
    'ExpDay','04-11-20'};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450041120, ...
{C1},{'o'},...
{'M450-02s'},'M450_1s_04-11-20',1)



    % M270 04-11-20
ConditionsM450041120= {'ExpType','DictyDB_M270','TpsComp','1s','FitParams',fitparams,...
    'ExpDay','04-11-20'};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450041120, ...
{C1_27},{'o'},...
{'M270-02s'},'M270_1s_04-11-20',1)
