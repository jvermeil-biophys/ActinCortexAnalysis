close all;
clear all;
clc;

warning('off')

%% Paths & options
%% Paths
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';
RawdataFolder = 'D:\Data\Raw';
MatfileFolder = 'D:\Data\MatFile';
DropboxDataFolder = '';
FigureFolder = 'D:\Data\Figures';

%% Options
set(0, 'defaultfigureposition', get(0, 'Screensize'));

set(0,'DefaultFigureWindowStyle','docked')

% Couleurs
Cwt = [20 109 159]./255; % WT

Csev = [182 109 106]./255; % SevA- mutant

Cfim  = [255 255 109]./255; % FimA- mutant

CaA = [0 73 73]./255; % abpA- mutant

CaC  = [255 182 119]./255; % abpC-

Cwt27 = [146 109 219]./255./2; % WT M270

Cdm = [146 0 0]./255; % DMSO

Cl  = [150 0 255]./255; % Latranculin

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


% diff tmp comp
Cdc1  = [255 127 63]./255; % diff tmp comp
Cdc2  = [127 63 255]./255; % diff tmp comp
Cdc3  = [63 255 127]./255; % diff tmp comp


AUTO = 0;

% pour v2d
NotSaved = {};

%pour d2p
PLOTDFCurve = 1;
PLOTDFComp = 1;

% d2m
PLOTD2M = 0;
VERBOSED2M = 0;

%p2s
PLOTPS = 1;

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Var2Data Classic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 Ax2 DMSO/LatA (Nishit)

Var2Data_wFluo('20-05-20','R40cst',1.1,'M1','DictyAx2_DMSO',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('20-05-20','R40cst',1.1,'M2','DictyAx2_DMSO',4504,RawdataFolder,MatfileFolder);

Var2Data_wFluo('22-05-20','R40cst',1.1,'M1','DictyAx2_DMSO',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('22-05-20','R40cst',1.1,'M1','DictyAx2_LatA',4504,RawdataFolder,MatfileFolder);

Var2Data_wFluo('25-05-20','R40cst',1.1,'M1','DictyAx2_DMSO',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('25-05-20','R40cst',1.1,'M1','DictyAx2_LatA',4504,RawdataFolder,MatfileFolder);

Var2Data_wFluo('27-05-20','R40cst',1.1,'M1','DictyAx2_DMSO',4504,RawdataFolder,MatfileFolder);

RemoveFromData(MatfileFolder,{'25-05-20_M1_P2_C9_R40','20-05-20_M1_P1_C6_R40'},'DictyAx2')

%% M270 Ax2 WT (Nishit)

Var2Data_wFluo('25-08-20','R90cst',0.82,'M1','DictyAx2_M270',2691,RawdataFolder,MatfileFolder);

Var2Data_wFluo('26-08-20','R90cst',0.82,'M1','DictyAx2_M270',2691,RawdataFolder,MatfileFolder);

% removing
RemoveFromData(MatfileFolder,{'25-08-20_M1_P1_C7_R90','25-08-20_M1_P1_C7_R90'},'DictyAx2')

%% M450 Ax2 WT + Mutants (DictyBase)

Var2Data_wFluo('23-09-20','R40cst',1.1,'M1','DictyDB_WT',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('23-09-20','R40cst',1.1,'M2','DictyDB_SevA',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('23-09-20','R40cst',1.1,'M13','DictyDB_FimA',4504,RawdataFolder,MatfileFolder);

Var2Data_wFluo('25-09-20','R40cst',1.1,'M1','DictyDB_SevA',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('25-09-20','R40cst',1.1,'M2','DictyDB_FimA',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('25-09-20','R40cst',1.1,'M3','DictyDB_WT',4504,RawdataFolder,MatfileFolder);

Var2Data_wFluo('06-10-20','R40cst',1.1,'M1','DictyDB_abpC',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('06-10-20','R40cst',1.1,'M2','DictyDB_abpA',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('06-10-20','R40cst',1.1,'M3','DictyDB_WT',4504,RawdataFolder,MatfileFolder);

Var2Data_wFluo('08-10-20','R40cst',1.1,'M1','DictyDB_abpA',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('08-10-20','R40cst',1.1,'M2','DictyDB_abpC',4504,RawdataFolder,MatfileFolder);
Var2Data_wFluo('08-10-20','R40cst',1.1,'M3','DictyDB_WT',4504,RawdataFolder,MatfileFolder);

% removing
RemoveFromData(MatfileFolder,{'06-10-20_M3_P1_C7_R40','06-10-20_M1_P1_C1_R40','06-10-20_M1_P1_C2_R40','06-10-20_M1_P1_C4_R40',...
    '08-10-20_M3_P1_C1_R40','08-10-20_M3_P1_C4_R40','08-10-20_M3_P1_C6_R40','08-10-20_M2_P1_C2_R40',...
    '23-09-20_M1_P1_C9_R40','23-09-20_M2_P1_C3_R40','23-09-20_M2_P1_C11_R40',...
    '25-09-20_M1_P1_C12_R40','25-09-20_M1_P1_C8_R40','25-09-20_M3_P1_C1_R40',},'DictyDB')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data2Peak %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 Ax2 DMSO/LatA (Nishit)

Data2Peaks('DictyAx2_DMSO',1,MatfileFolder,FigureFolder);
Data2Peaks('DictyAx2_LatA',1,MatfileFolder,FigureFolder);

%% M270 Ax2 WT (Nishit)

Data2Peaks('DictyAx2_M270',1,MatfileFolder,FigureFolder);

%% M450 Ax2 WT + Mutants (DictyBase)

Data2Peaks('DictyDB_WT',1,MatfileFolder,FigureFolder);
Data2Peaks('DictyDB_SevA',1,MatfileFolder,FigureFolder);
Data2Peaks('DictyDB_FimA',1,MatfileFolder,FigureFolder);
Data2Peaks('DictyDB_abpA',1,MatfileFolder,FigureFolder);
Data2Peaks('DictyDB_abpC',1,MatfileFolder,FigureFolder);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Peaks2Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% pooling

Peaks2Stats('Dicty_WT',PLOTPS,0.9,MatfileFolder,FigureFolder);

%% M450 Ax2 DMSO/LatA (Nishit)

Peaks2Stats('DictyAx2_DMSO',PLOTPS,0.9,MatfileFolder,FigureFolder);
Peaks2Stats('DictyAx2_LatA',PLOTPS,0.9,MatfileFolder,FigureFolder);

%% M270 Ax2 WT (Nishit)
Peaks2Stats('DictyAx2_M270',PLOTPS,0.9,MatfileFolder,FigureFolder);

%% M450 Ax2 WT + Mutants (DictyBase)
Peaks2Stats('DictyDB_WT',PLOTPS,0.9,MatfileFolder,FigureFolder);
Peaks2Stats('DictyDB_FimA',PLOTPS,0.9,MatfileFolder,FigureFolder);
Peaks2Stats('DictyDB_SevA',PLOTPS,0.9,MatfileFolder,FigureFolder);
Peaks2Stats('DictyDB_abpA',PLOTPS,0.9,MatfileFolder,FigureFolder);
Peaks2Stats('DictyDB_abpC',PLOTPS,0.9,MatfileFolder,FigureFolder);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stats2Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

Stats2PlotsDictys({'DictyDB_WT','DictyDB_SevA','DictyDB_FimA'},{Cwt,Csev,Cfim},'DictyMutantsComp',1,1,1,0,0.9,MatfileFolder,FigureFolder);

Stats2PlotsDictys({'Dicty_WT','DictyDB_abpA','DictyDB_abpC','DictyDB_FimA'},{Cwt,CaA,CaC,Cfim},'DictyMutantsComp',1,1,1,0,0.9,MatfileFolder,FigureFolder);

Stats2PlotsDictys({'DictyAx2_DMSO','DictyAx2_M270'},{Cdm,Cwt27},'DictyBeadsSize',1,1,1,0,0.9,MatfileFolder,FigureFolder);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Var2Data Comp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 Ax2 DMSO/LatA (Nishit)

Var2Data_Comp_V1('20-05-20','R40',1.1,'M1','DictyAx2_DMSO',95,143,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-05-20','R40',1.1,'M2','DictyAx2_DMSO',95,143,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('22-05-20','R40',1.1,'M1','DictyAx2_DMSO',95,143,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('22-05-20','R40',1.1,'M1','DictyAx2_LatA',95,143,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('25-05-20','R40',1.1,'M1','DictyAx2_DMSO',95,143,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('25-05-20','R40',1.1,'M1','DictyAx2_LatA',95,143,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('27-05-20','R40',1.1,'M1','DictyAx2_DMSO',95,143,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

%% M270 Ax2 WT (Nishit)

Var2Data_Comp_V1('25-08-20','R90',0.82,'M1','DictyAx2_M270',95,143,'1s',2691,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('26-08-20','R90',0.82,'M1','DictyAx2_M270',95,143,'1s',2691,RawdataFolder,MatfileFolder,FigureFolder);

%% M450 Ax2 WT + Mutants (DictyBase)

Var2Data_Comp_V1('23-09-20','R40',1.1,'M1','DictyDB_WT',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('23-09-20','R40',1.1,'M2','DictyDB_SevA',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('23-09-20','R40',1.1,'M3','DictyDB_FimA',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('25-09-20','R40',1.1,'M1','DictyDB_SevA',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('25-09-20','R40',1.1,'M2','DictyDB_FimA',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('25-09-20','R40',1.1,'M3','DictyDB_WT',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('06-10-20','R40',1.1,'M1','DictyDB_abpC',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('06-10-20','R40',1.1,'M2','DictyDB_abpA',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('06-10-20','R40',1.1,'M3','DictyDB_WT',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('08-10-20','R40',1.1,'M1','DictyDB_abpA',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('08-10-20','R40',1.1,'M2','DictyDB_abpC',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('08-10-20','R40',1.1,'M3','DictyDB_WT',95,167,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

%% MultiTime Comp M450 Ax2 (Nishit)

Var2Data_Comp_V1_R248('13-02-20','R248',1.1,'M1','DictyAx2-Comp',21,24,98,193,383,4438,RawdataFolder,MatfileFolder);
Var2Data_Comp_V1_R248('20-02-20','R248',1.1,'M1','DictyAx2-Comp',21,24,98,193,383,4438,RawdataFolder,MatfileFolder);

%% MultiRate Comp M450 WTDB

%20-10-20
Var2Data_Comp_V1('20-10-20','R50_04s',0.85,'M1','DictyDB_M450',40,76,'02s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_1s',0.85,'M1','DictyDB_M450',100,136,'05s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_2s',0.85,'M1','DictyDB_M450',133,169,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_4s',0.85,'M1','DictyDB_M450',133,169,'2s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_8s',0.85,'M1','DictyDB_M450',133,169,'4s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_14s',0.85,'M1','DictyDB_M450',140,176,'7s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_24s',0.85,'M1','DictyDB_M450',141,177,'12s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('20-10-20','R50_04s',0.85,'M2','DictyDB_M450',40,76,'02s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_1s',0.85,'M2','DictyDB_M450',100,136,'05s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_2s',0.85,'M2','DictyDB_M450',133,169,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_4s',0.85,'M2','DictyDB_M450',133,169,'2s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_8s',0.85,'M2','DictyDB_M450',133,169,'4s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_14s',0.85,'M2','DictyDB_M450',140,176,'7s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_24s',0.85,'M2','DictyDB_M450',141,177,'12s',4504,RawdataFolder,MatfileFolder,FigureFolder);

% 21-10-20
Var2Data_Comp_V1('21-10-20','R50_04s',0.85,'M1','DictyDB_M450',40,76,'02s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_1s',0.85,'M1','DictyDB_M450',100,136,'05s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_2s',0.85,'M1','DictyDB_M450',133,169,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_4s',0.85,'M1','DictyDB_M450',133,169,'2s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_8s',0.85,'M1','DictyDB_M450',133,169,'4s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_14s',0.85,'M1','DictyDB_M450',140,176,'7s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_24s',0.85,'M1','DictyDB_M450',141,177,'12s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('21-10-20','R50_04s',0.85,'M2','DictyDB_M450',40,76,'02s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_1s',0.85,'M2','DictyDB_M450',100,136,'05s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_2s',0.85,'M2','DictyDB_M450',133,169,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_4s',0.85,'M2','DictyDB_M450',133,169,'2s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_8s',0.85,'M2','DictyDB_M450',133,169,'4s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_14s',0.85,'M2','DictyDB_M450',140,176,'7s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R50_24s',0.85,'M2','DictyDB_M450',141,177,'12s',4504,RawdataFolder,MatfileFolder,FigureFolder);

% 22-10-20
Var2Data_Comp_V1('22-10-20','R50_04s',0.85,'M1','DictyDB_M450',40,76,'02s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('22-10-20','R50_4s',0.85,'M1','DictyDB_M450',133,169,'2s',4504,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('22-10-20','R50_24s',0.85,'M1','DictyDB_M450',141,177,'12s',4504,RawdataFolder,MatfileFolder,FigureFolder);


% 28-10-20
Var2Data_Comp_MultiRatePerCell_V1('28-10-20','R50-Multi',0.85,'M1','DictyDB_M450-Multi',...
    'RmpTimeDatas',4504,RawdataFolder,MatfileFolder,FigureFolder);



%% MultiRate Comp M270 WTDB

%20-10-20
Var2Data_Comp_V1('20-10-20','R90_04s',0.85,'M1','DictyDB_M270',40,76,'02s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R90_1s',0.85,'M1','DictyDB_M270',100,136,'05s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R90_2s',0.85,'M1','DictyDB_M270',133,169,'1s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R90_4s',0.85,'M1','DictyDB_M270',133,169,'2s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R90_8s',0.85,'M1','DictyDB_M270',133,169,'4s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R90_14s',0.85,'M1','DictyDB_M270',140,176,'7s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R90_24s',0.85,'M1','DictyDB_M270',141,177,'12s',2691,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('20-10-20','R50_04s',0.85,'M2','DictyDB_M270',40,76,'02s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_1s',0.85,'M2','DictyDB_M270',100,136,'05s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_2s',0.85,'M2','DictyDB_M270',133,169,'1s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_4s',0.85,'M2','DictyDB_M270',133,169,'2s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_8s',0.85,'M2','DictyDB_M270',133,169,'4s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_14s',0.85,'M2','DictyDB_M270',140,176,'7s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('20-10-20','R50_24s',0.85,'M2','DictyDB_M270',141,177,'12s',2691,RawdataFolder,MatfileFolder,FigureFolder);

% 21-10-20
Var2Data_Comp_V1('21-10-20','R90_04s',0.85,'M1','DictyDB_M270',40,76,'02s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_1s',0.85,'M1','DictyDB_M270',100,136,'05s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_2s',0.85,'M1','DictyDB_M270',133,169,'1s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_4s',0.85,'M1','DictyDB_M270',133,169,'2s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_8s',0.85,'M1','DictyDB_M270',133,169,'4s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_14s',0.85,'M1','DictyDB_M270',140,176,'7s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_24s',0.85,'M1','DictyDB_M270',141,177,'12s',2691,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('21-10-20','R90_04s',0.85,'M2','DictyDB_M270',40,76,'02s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_1s',0.85,'M2','DictyDB_M270',100,136,'05s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_2s',0.85,'M2','DictyDB_M270',133,169,'1s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_4s',0.85,'M2','DictyDB_M270',133,169,'2s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_8s',0.85,'M2','DictyDB_M270',133,169,'4s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_14s',0.85,'M2','DictyDB_M270',140,176,'7s',2691,RawdataFolder,MatfileFolder,FigureFolder);
Var2Data_Comp_V1('21-10-20','R90_24s',0.85,'M2','DictyDB_M270',141,177,'12s',2691,RawdataFolder,MatfileFolder,FigureFolder);

%% LongueManip comp M450/M270 WTDB
Var2Data_Comp_V1('04-11-20','R50_2s',0.85,'M1','DictyDB_M450',133,169,'1s',4504,RawdataFolder,MatfileFolder,FigureFolder);

Var2Data_Comp_V1('04-11-20','R90_2s',0.85,'M2','DictyDB_M270',133,169,'1s',2691,RawdataFolder,MatfileFolder,FigureFolder);

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Meca2Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% M450 mutant 
ConditionsM450Multi = {'ExpType','DictyDB_WT','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_FimA','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_abpA','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_abpC','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_SevA','TpsComp','1s','FitParams','Strain100-R20.7';...
    }

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450Multi, ...
{Cwt,Cfim,CaA,CaC,Csev},...
{'WT','FimA-','abpA-','abpC-','SevA-'},'M450Mutants',1)

%% Dicty NS Multitime
ConditionsM450NSMulti = {'ExpType','DictyAx2-Comp','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyAx2-Comp','TpsComp','2s','FitParams','Strain100-R20.7';...
    'ExpType','DictyAx2-Comp','TpsComp','4s','FitParams','Strain100-R20.7';...
    }
Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450NSMulti, ...
{C1,C2,C4},...
{'M450-NS-1s','M450-NS-2s','M450-NS-4s'},'M450MultiRate',1)

%% M450 multirate    
ConditionsM450Multi = {'ExpType','DictyDB_M450','TpsComp','02s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','05s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','2s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','4s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','7s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','12s','FitParams','Strain100-R20.7';}

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450Multi, ...
{C02,C05,C1,C2,C4,C7,C12},...
{'M450-02s','M450-05s','M450-1s','M450-2s','M450-4s','M450-7s','M450-12s'},'M450MultiRate',1,...
'Exclude',{'ExpDay','04-11-20'})

    % M450 multirate 041120    
ConditionsM450041120= {'ExpType','DictyDB_M450','TpsComp','1s','FitParams','Strain100-R20.7',...
    'ExpDay','04-11-20';
    }

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450041120, ...
{C1},...
{'M450-1s'},'M450_1s_04-11-20',1)
          

%% M270 multirate  

ConditionsM270Multi = {'ExpType','DictyDB_M270','TpsComp','02s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','05s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','2s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','4s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','7s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','12s','FitParams','Strain100-R20.7';}

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM270Multi, ...
{C02_27,C05_27,C1_27,C2_27,C4_27,C7_27,C12_27},...
{'M270-02s','M270-05s','M270-1s','M270-2s','M270-4s','M270-7s','M270-12s'},'M270DMultiRate',1,...
'Exclude',{'ExpDay','04-11-20'})


%% M450270 multirate  
ConditionsM450270Multi = {'ExpType','DictyDB_M450','TpsComp','02s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','02s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','05s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','05s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','2s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','2s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','4s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','4s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','7s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','7s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M450','TpsComp','12s','FitParams','Strain100-R20.7';
    'ExpType','DictyDB_M270','TpsComp','12s','FitParams','Strain100-R20.7';}

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450270Multi, ...
{C02,C02_27,C05,C05_27,C1,C1_27,C2,C2_27,C4,C4_27,C7,C7_27,C12,C12_27},...
{'M450-02s','M270-02s','M450-05s','M270-05s','M450-1s','M270-1s','M450-2s',...
'M270-2s','M450-4s','M270-4s','M450-7s','M270-7s','M450-12s','M270-12s'},'M450270DMultiRate',0,...
'Exclude',{'ExpDay','04-11-20'})



%% M50/M270 original
ConditionsM450270_1s = {'ExpType','DictyDB_M450','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_WT','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyDB_M270','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyAx2_DMSO','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyAx2_M270','TpsComp','1s','FitParams','Strain100-R20.7'}

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450270_1s, ...
{C1,C1,C1_27,Cwt,Cwt27},...
{'DB-M450-1s','DB-WT-1s','DB-M270-1s','NS-M450-1s','NS-M270-1s'},'M45-27_1s',1,...
'Exclude',{'ExpDay','04-11-20'})



 % only NS
ConditionsM450270_1s_NS = {'ExpType','DictyAx2_DMSO','TpsComp','1s','FitParams','Strain100-R20.7';...
    'ExpType','DictyAx2_M270','TpsComp','1s','FitParams','Strain100-R20.7'}

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450270_1s_NS, ...
{Cwt,Cwt27},...
{'M450-1s','M270-1s'},'M45-27_1s_NS',1)




    % M450 04-11-20    
ConditionsM450041120= {'ExpType','DictyDB_M450','TpsComp','1s','FitParams','Strain100-R20.7',...
    'ExpDay','04-11-20'};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450041120, ...
{C1},...
{'M450-02s'},'M450_1s_04-11-20',1)



    % M270 04-11-20 
ConditionsM450041120= {'ExpType','DictyDB_M270','TpsComp','1s','FitParams','Strain100-R20.7',...
    'ExpDay','04-11-20'};

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsM450041120, ...
{C1_27},...
{'M270-02s'},'M270_1s_04-11-20',1)
          

