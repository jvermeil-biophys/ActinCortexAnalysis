close all;
clear all;
clc;

%% Paths
RawdataFolder = 'D:\Data\Raw';
MatfileFolder = 'D:\Data\MatFile';
FigureFolder = 'D:\Data\Figures';

set(0, 'defaultfigureposition', get(0, 'Screensize'));

set(0,'DefaultFigureWindowStyle','docked')

%% Parameters

% Colors
Cdm = [0 109 219]./255; % control
Cnd = Cdm/2;

Cl  = [73 0 146]./255; % Latranculin

Cb  = [219 209 0]./255; % blebbi
Cy  = [255 255 109 ]./255; % y27
Ccl = [146 0 0]./255; % calyculinA 1nM
Ccl3 = Ccl/2; % calyculinA 3nM
Ccl10 = Ccl3/2; % calyculinA 10nM

Cs  = [36 255 36]./255; % smifh2

Cck = [182 109 255]./255;

Cin = [0 73 73]./255; % inside
Co  = [255 182 119]./255; % outside

% for r2v
AUTO = 0; % automatic analysis mode (skip when user input is needed)

% for v2d
NotSaved = {}; % gather a list of non saved data (because of sorting criteria)

% for d2p
PLOTDP = 0; % plot curve with cortex level
alignlvl = 0.9; % for probcum

% for p2s
PLOTPS = 0;

% for s2p
PEAKS    = 1;
ACTI  = 1;
DIST   = 1;
TIME = 1;

return % stop execution here 

%% Analysis

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Res2Var %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    %% Nouvelles manips Inhib   
    
    % Culture 1
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DCLA_CK666'    ,1:3 ,'Results',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M2','DCLA_DMSO'     ,1,'Results',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M1','DCLA_Blebbi'   ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DCLA_LatA500'  ,2,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 2
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M1','DCLA_Blebbi'   ,1:3 ,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DCLA_LatA500'  ,2,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M2','DCLA_CK666'    ,1:3 ,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 3
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M1','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M2','DCLA_Blebbi'   ,1:3 ,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DCLA_LatA500'  ,2,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M2','DCLA_CK666'  ,1:3 ,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    % Culture 5
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DCLA_SMIFH2'   ,1,'Results',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DCLA_DMSO'     ,2,'Results',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M2','DCLA_ML141'    ,1,'Results',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M2','DCLA_SMIFH2'   ,2,'Results',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 6
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DCLA_ML141'    ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DCLA_DMSO'     ,2,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M2','DCLA_ML141'    ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DCLA_DMSO'     ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DCLA_SMIFH2'   ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DCLA_ML141'    ,2,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('20-07-18', '5mT', 'M1','DCLA_Blebbi'   ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('20-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 7
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M1','DCLA_ML141'    ,1,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DCLA_Blebbi'   ,1,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DCLA_DMSO'     ,2,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M1','DCLA_LatA500'  ,1,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M1','DCLA_CK666'    ,2,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Arpin 1
    Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M1','DCLA_ArpinWT' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M2','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M3','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Arpin 2
    Res2Var_wFluo_multiZ('04-10-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('04-10-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ('05-10-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('05-10-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    
    % Culture 9
    Res2Var_wFluo_multiZ('15-11-18', '5mT', 'M1','DCLA_Y27',1:3 ,'Results',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('15-11-18', '5mT', 'M2','DCLA_NoDrug',1 ,'Results',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    
    % Culture 10
    Res2Var_wFluo_multiZ('12-12-18', '5mT', 'M1','DCLA_NoDrug',1:3 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('12-12-18', '5mT', 'M2','DCLA_Nase',1:3 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 11
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M1','DCLA_ArpC4KO',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M1','DCLA_Y27',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M2','DCLA_ArpC4WT',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M1','DCLA_Y27',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M1','DCLA_SMIFH2',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M2','DCLA_NoDrug',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M2','DCLA_Nase',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 12
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M1','DCLA_SMIFH2',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M1','DCLA_Y27',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DCLA_SMIFH2',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DCLA_LatA500',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('21-12-18', '5mT','M1','DCLA_LatNase',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    % MyoII KO 1
    Res2Var_wFluo_multiZ('11-03-19', '5mT','M1','DCLA_MyoII-WT',1:2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('11-03-19', '5mT','M2','DCLA_MyoII-KO',1:2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('12-03-19', '5mT','M1','DCLA_MyoII-KO',1:2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('12-03-19', '5mT','M2','DCLA_MyoII-WT',1:2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 13
    Res2Var_wFluo_multiZ('20-03-19', '5mT','M1','DCLA_LatTryp',1:2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('21-03-19', '5mT','M1','DCLA_LatTryp',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    % Culture 14
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M1','DCLA_CalA10',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M1','DCLA_CalA1',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M2','DCLA_CalA1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M2','DCLA_CalA10',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('02-04-19', '5mT', 'M2','DCLA_CalA1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('02-04-19', '5mT', 'M2','DCLA_CalA10',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % SiRNA1
    Res2Var_wFluo_multiZ('29-04-19', '5mT', 'M1','DCLA_SiVim',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('29-04-19', '5mT', 'M1','DCLA_SiCtrl',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('30-04-19', '5mT', 'M1','DCLA_SiVim+LatA',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('30-04-19', '5mT', 'M1','DCLA_SiCtrl+LatA',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % culture 15
    Res2Var_wFluo_multiZ('13-05-19', '5mT', 'M1','DCLA_CalA3',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('13-05-19', '5mT', 'M1','DCLA_CalA1',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('13-05-19', '5mT', 'M3','DCLA_NoDrug',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M1','DCLA_CalA3',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M1','DCLA_CalA1',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M2','DCLA_DMSO',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M2','DCLA_NoDrug',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % culture 16 avec SiRNA2
    Res2Var_wFluo_multiZ('15-05-19', '5mT', 'M1','DCLA_SiVim+LatA',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('15-05-19', '5mT', 'M1','DCLA_SiCtrl+LatA',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('15-05-19', '5mT', 'M2','DCLA_SiVim',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('15-05-19', '5mT', 'M2','DCLA_SiCtrl',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ('16-05-19', '5mT', 'M1','DCLA_SiCtrl',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('16-05-19', '5mT', 'M1','DCLA_DMSO',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    % culture 17
        
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M1','DCLA_CalA1',1 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)        
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M1','DCLA_CalA3',2 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
        Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M2','DCLA_DMSO',1 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)        
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M2','DCLA_CalA3',2 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M1','DCLA_Sep',1:2 ,'Separation',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M2','DCLA_Sep',1 ,'Separation',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    % culture 18
    Res2Var_wFluo_multiZ('02-07-19', '5mT', 'M1','DCLA_DMSO',1 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('02-07-19', '5mT', 'M2','DCLA_CalA1',1 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('02-07-19', '5mT', 'M2','DCLA_CalA3',2 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    %% endocytose et division
        Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DCLA_Endo' ,2 ,'Endo',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
        Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DCLA_Endo' ,1,'Endo',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
        Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M2','DCLA_Endo' ,1:2,'Endo',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
        Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M2','DCLA_Endo',1 ,'Endo',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
        Res2Var_wFluo_multiZ('12-03-19', '5mT', 'M1','DCLA_Endo' ,1 ,'Endo',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
        Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M2','DCLA_Endo',2 ,'Endo',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
        Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DCLA_Div',1 ,'Div',48,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)    
    
        Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M2','DCLA_Div',2 ,'Div',48,'28-08-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)    
    
        Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DCLA_Div',2 ,'Div',48,'18-07-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
          
    %% Dicty
    
     Res2Var_wFluo_multiZ('05-11-18', '5mT', 'M1','Dicty_WT',1 ,'Results',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
     Res2Var_wFluo_multiZ('05-12-18', '5mT', 'M1','Dicty_WT',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % endocytose
    
    Res2Var_wFluo_multiZ('05-11-18', '5mT', 'M1','Dicty_Endo',1 ,'Endo',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('05-12-18', '5mT', 'M1','Dicty_Endo',1 ,'Endo',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % div
    
    Res2Var_wFluo_multiZ('05-12-18', '5mT', 'M1','Dicty_Div',1 ,'Div',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Nouvelles manip dictys (Ax2)
    
    Res2Var_wFluo_multiZ('20-05-20', '5mT', 'M1','DictyAx2_DMSO',1 ,'Results',48,'19-05-20_DepthographM450' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('20-05-20', '5mT', 'M2','DictyAx2_DMSO',1 ,'Results',48,'19-05-20_DepthographM450' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ('22-05-20', '5mT', 'M1','DictyAx2_DMSO',1 ,'Results',48,'19-05-20_DepthographM450' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('22-05-20', '5mT', 'M1','DictyAx2_LatA',2 ,'Results',48,'19-05-20_DepthographM450' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('25-05-20', '5mT', 'M1','DictyAx2_DMSO',1 ,'Results',48,'19-05-20_DepthographM450' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('25-05-20', '5mT', 'M1','DictyAx2_LatA',2 ,'Results',48,'19-05-20_DepthographM450' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ('03-08-20', '5mT', 'M1','DictyAx2_BSA30',1 ,'Results',48,'19-05-20_DepthographM450' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    %% RPE1 
         Res2Var_wFluo_multiZ('04-04-19', '5mT', 'M1','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
         Res2Var_wFluo_multiZ('10-04-19', '5mT', 'M1','RPE1',[1 2] ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
         Res2Var_wFluo_multiZ('11-04-19', '5mT', 'M1','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)    
         Res2Var_wFluo_multiZ('11-04-19', '5mT', 'M2','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)    
         Res2Var_wFluo_multiZ('11-04-19', '5mT', 'M3','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
         Res2Var_wFluo_multiZ('17-04-19', '5mT', 'M1','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
         Res2Var_wFluo_multiZ('17-04-19', '5mT', 'M2','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Inside/Outside Datas %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Manip etalonnage paula et vieux 0%
    
    Res2Var_wFluo_multiZ('22-01-18', '5mT', 'M1','DC_Inside_0pc'    ,1:2,'Inside',31,'03-02-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('22-01-18', '5mT', 'M2','DC_Inside_0pc'    ,1:2,'Inside',31,'03-02-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('23-01-18', '5mT', 'M1','DC_Inside_0pc'    ,1:2,'Inside',31,'03-02-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('23-01-18', '5mT', 'M2','DC_Inside_0pc'    ,1:2,'Inside',31,'03-02-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('25-03-19', '5mT', 'M1','DC_Inside_0pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)  
    Res2Var_wFluo_multiZ('25-03-19', '5mT', 'M2','DC_Inside_5pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)  
    Res2Var_wFluo_multiZ('25-03-19', '5mT', 'M3','DC_Inside_10pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)  
    Res2Var_wFluo_multiZ('25-03-19', '5mT', 'M4','DC_Inside_12-5pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
        Res2Var_wFluo_multiZ('26-03-19', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)  
    Res2Var_wFluo_multiZ('26-03-19', '5mT', 'M2','DC_Inside_15pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)  
    Res2Var_wFluo_multiZ('26-03-19', '5mT', 'M3','DC_Inside_10pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)  
    Res2Var_wFluo_multiZ('26-03-19', '5mT', 'M4','DC_Inside_12-5pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    %% 7.5 pc dans manip
    
    % Culture 1
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M2','DC_Inside_DMSO'    ,1,'Inside',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M2','DC_Inside_7-5pc'    ,2,'Inside',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DC_Inside_DMSO'    ,1,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DC_Inside_7-5pc'    ,2,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 2
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DC_Inside_DMSO'    ,1,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,2:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DC_Inside_DMSO'    ,1,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,2:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 3
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M1','DC_Inside_DMSO'    ,1,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,2:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DC_Inside_DMSO'    ,1,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,2:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    % Culture 4
    Res2Var_wFluo_multiZ('02-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',46,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('02-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',46,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('03-07-18', '5mT','M1','DC_Inside_7-5pc'    ,1:3,'Inside',49,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('03-07-18', '5mT','M2','DC_Inside_7-5pc'    ,1:3,'Inside',49,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 5
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DC_Inside_7-5pc'    ,1,'Inside',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DC_Inside_DMSO'    ,2,'Inside',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 6
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DC_Inside_DMSO'    ,2,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DC_Inside_DMSO'    ,1,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,2,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('20-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Culture 7
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1,'Inside',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DC_Inside_DMSO'    ,2,'Inside',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Arpin 1
    Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 8
    Res2Var_wFluo_multiZ('29-10-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    % Culture 10
    Res2Var_wFluo_multiZ('12-12-18', '5mT', 'M2','DC_Nase_Inside',1:3 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 12
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M1','DC_Inside_7-5pc',2 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DC_Inside_7-5pc',1 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DC_Inside_7-5pc',2 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    
    % culture de 2019
    Res2Var_wFluo_multiZ('14-05-19', '5mT','M2','DC_Inside_DMSO',1 ,'Inside',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('16-05-19', '5mT','M1','DC_Inside_DMSO',2 ,'Inside',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ('02-07-19', '5mT','M1','DC_Inside_DMSO',1 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    %% Manip 5 - 3 - 1 mT
    
    Res2Var_wFluo_multiZ('21-01-20', '5mT','M1','DC_Inside_5mT',1 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('21-01-20', '3mT','M1','DC_Inside_3mT',1 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('21-01-20', '1mT','M1','DC_Inside_1mT',1 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
        Res2Var_wFluo_multiZ('21-01-20', '5mT','M2','DC_Inside_5mT',1:2 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('21-01-20', '3mT','M2','DC_Inside_3mT',1:2 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ('21-01-20', '1mT','M2','DC_Inside_1mT',1:2 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    %% autre manips
    
    Res2Var_wFluo_multiZ('21-12-18', '5mT','M1','DC_Nase_Inside',1 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 14
    Res2Var_wFluo_multiZ('02-04-19', '5mT', 'M1','DC_TrueOutside',1 ,'Outside',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    
    % Culture 15
    Res2Var_wFluo_multiZ('13-05-19', '5mT', 'M2','DC_TrueOutside',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder)
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Var2Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    %% Nouvelles manips Inhib
       
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18','5mT',0.56,'M1','DCLA_CK666',4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18','5mT',0.56,'M2','DCLA_DMSO',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18','5mT',0.56,'M1','DCLA_Blebbi',4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18','5mT',0.56,'M2','DCLA_DMSO',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18','5mT',0.56,'M1','DCLA_Blebbi',4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18','5mT',0.56,'M2','DCLA_DMSO',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18','5mT',0.56,'M1','DCLA_DMSO',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18','5mT',0.56,'M2','DCLA_CK666',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18','5mT',0.56,'M1','DCLA_DMSO',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18','5mT',0.56,'M2','DCLA_Blebbi',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18','5mT',0.56,'M1','DCLA_DMSO',4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18','5mT',0.56,'M2','DCLA_CK666',4415,...
        MatfileFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18','5mT',0.56,'M1','DCLA_DMSO',4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18','5mT',0.56,'M1','DCLA_SMIFH2',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18','5mT',0.56,'M2','DCLA_SMIFH2',4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18','5mT',0.56,'M1','DCLA_DMSO'  ,4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18','5mT',0.56,'M1','DCLA_DMSO'   ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18','5mT',0.56,'M2','DCLA_SMIFH2' ,4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-07-18','5mT',0.56,'M1','DCLA_Blebbi' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4415,...
        MatfileFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-11-18','5mT',0.56,'M2','DCLA_NoDrug',4438,...
        MatfileFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-12-18','5mT',0.56,'M1','DCLA_NoDrug',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-12-18','5mT',0.56,'M1','DCLA_SMIFH2',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-12-18','5mT',0.56,'M2','DCLA_NoDrug',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18','5mT',0.56,'M1','DCLA_SMIFH2',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18','5mT',0.56,'M2','DCLA_SMIFH2',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('01-04-19','5mT',0.56,'M1','DCLA_CalA1',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('01-04-19','5mT',0.56,'M2','DCLA_CalA1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-04-19','5mT',0.56,'M2','DCLA_CalA1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('13-05-19','5mT',0.56,'M1','DCLA_CalA1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('13-05-19','5mT',0.56,'M3','DCLA_NoDrug',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19','5mT',0.56,'M1','DCLA_CalA1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19','5mT',0.56,'M2','DCLA_DMSO',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('16-05-19', '5mT',0.56, 'M1','DCLA_DMSO',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18','5mT',0.56,'M2','DCLA_LatA500',4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18','5mT',0.56,'M1','DCLA_LatA500',4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18','5mT',0.56,'M1','DCLA_LatA500',4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18','5mT',0.56,'M2','DCLA_LatA500',4438,...
        MatfileFolder)};
    
    
        NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18','5mT',0.56,'M2','DCLA_ML141',4415,...
            MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18','5mT',0.56,'M1','DCLA_ML141' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18','5mT',0.56,'M2','DCLA_ML141' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18','5mT',0.56,'M2','DCLA_ML141'  ,4415,...
        MatfileFolder)};
    
    
       % Culture 7
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18','5mT',0.56,'M1','DCLA_ML141'  ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18','5mT',0.56,'M2','DCLA_Blebbi' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18','5mT',0.56,'M2','DCLA_DMSO'   ,4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('31-07-18','5mT',0.56,'M1','DCLA_LatA500' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('31-07-18','5mT',0.56,'M1','DCLA_CK666'   ,4438,...
        MatfileFolder)};
       

    NotSaved = {NotSaved{:},Var2Data_wFluo('13-05-19','5mT',0.56,'M1','DCLA_CalA3',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19','5mT',0.56,'M1','DCLA_CalA3',4438,...
        MatfileFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('01-04-19','5mT',0.56,'M1','DCLA_CalA10',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('01-04-19','5mT',0.56,'M2','DCLA_CalA10',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-04-19','5mT',0.56,'M2','DCLA_CalA10',4438,...
        MatfileFolder)};
    
        % Arpin 1
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56,'M1','DCLA_ArpinWT' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56,'M2','DCLA_ArpinKO' ,4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M1','DCLA_ArpinKO' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M2','DCLA_ArpinWT' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M3','DCLA_ArpinKO' ,4438,...
        MatfileFolder)};
    
    % Arpin 2
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-10-18', '5mT',0.56,'M1','DCLA_ArpinKO' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-10-18', '5mT',0.56,'M2','DCLA_ArpinWT' ,4438,...
        MatfileFolder)};

    
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-11-18','5mT',0.56,'M1','DCLA_Y27',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-12-18','5mT',0.56,'M2','DCLA_Nase',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-12-18','5mT',0.56,'M1','DCLA_Y27',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-12-18','5mT',0.56,'M2','DCLA_ArpC4WT',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-12-18','5mT',0.56,'M1','DCLA_ArpC4KO',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-12-18','5mT',0.56,'M1','DCLA_Y27',4438,...
        MatfileFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18','5mT',0.56,'M1','DCLA_Y27',4438,...
        MatfileFolder)};
    
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('21-12-18','5mT',0.56,'M1','DCLA_LatNase',4438,...
        MatfileFolder)};
    
    % KO myoII
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-03-19','5mT',0.56,'M1','DCLA_MyoII-WT',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-03-19','5mT',0.56,'M2','DCLA_MyoII-KO',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-03-19','5mT',0.56,'M1','DCLA_MyoII-KO',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-03-19','5mT',0.56,'M2','DCLA_MyoII-WT',4438,...
        MatfileFolder)};
    
    
    culture 14
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-03-19','5mT',0.56,'M1','DCLA_LatTryp',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('21-03-19','5mT',0.56,'M1','DCLA_LatTryp',4438,...
        MatfileFolder)};
    
    
    % SiRNA
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-04-19', '5mT',0.56, 'M1','DCLA_SiVim',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-04-19', '5mT',0.56, 'M1','DCLA_SiCtrl',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-04-19', '5mT',0.56, 'M1','DCLA_SiVim+LatA',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-04-19', '5mT',0.56, 'M1','DCLA_SiCtrl+LatA',4438,...
        MatfileFolder)};
    
    
    
     SiRNA
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-05-19', '5mT',0.56, 'M1','DCLA_SiVim+LatA',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-05-19', '5mT',0.56, 'M1','DCLA_SiCtrl+LatA',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-05-19', '5mT',0.56, 'M2','DCLA_SiVim',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-05-19', '5mT',0.56, 'M2','DCLA_SiCtrl',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('16-05-19', '5mT',0.56, 'M1','DCLA_SiCtrl',4438,...
        MatfileFolder)};
    
    
    % culture 17
        
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M1','DCLA_CalA1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M1','DCLA_CalA3',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M2','DCLA_DMSO',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M2','DCLA_CalA3',4438,...
        MatfileFolder)};
    
        NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M2','DCLA_Sep',4438,...
        MatfileFolder)};    
        NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M1','DCLA_Sep',4438,...
        MatfileFolder)};
    
    % culture 18        
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-19', '5mT',0.56, 'M1','DCLA_DMSO',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-19', '5mT',0.56, 'M2','DCLA_CalA1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-19', '5mT',0.56, 'M2','DCLA_CalA3',4438,...
        MatfileFolder)};

    %% Endocytose et division
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56, 'M1','DCLA_Endo',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18', '5mT',0.56, 'M2','DCLA_Endo',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18', '5mT',0.56, 'M2','DCLA_Endo',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-012-18', '5mT',0.56, 'M2','DCLA_Endo',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-03-19', '5mT',0.56, 'M1','DCLA_Endo',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19', '5mT',0.56, 'M2','DCLA_Endo',4438,...
        MatfileFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56, 'M1','DCLA_Div',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18', '5mT',0.56, 'M2','DCLA_Div',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56, 'M2','DCLA_Div',4438,...
        MatfileFolder)};
    
    %% Dicty
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-11-18', '5mT',0.56, 'M1','Dicty_WT',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-11-18', '5mT',0.56, 'M1','Dicty_Endo',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-12-18', '5mT',0.56, 'M1','Dicty_WT',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-12-18', '5mT',0.56, 'M1','Dicty_Endo',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-12-18', '5mT',0.56, 'M1','Dicty_Div',4438,...
        MatfileFolder)};
    
% nouvelle manips dicty
    
NotSaved = {NotSaved{:},Var2Data_wFluo('20-05-20', '5mT',1.1, 'M1','DictyAx2_DMSO',4504,...
        MatfileFolder)};
    
NotSaved = {NotSaved{:},Var2Data_wFluo('20-05-20', '5mT',1.1, 'M2','DictyAx2_DMSO',4504,...
        MatfileFolder)};
    
    
NotSaved = {NotSaved{:},Var2Data_wFluo('22-05-20', '5mT',1.1, 'M1','DictyAx2_DMSO',4504,...
        MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('22-05-20', '5mT',1.1, 'M1','DictyAx2_LatA',4504,...
        MatfileFolder)};
    
      
NotSaved = {NotSaved{:},Var2Data_wFluo('25-05-20', '5mT',1.1, 'M1','DictyAx2_DMSO',4504,...
        MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('25-05-20', '5mT',1.1, 'M1','DictyAx2_LatA',4504,...
        MatfileFolder)};

    NotSaved = {NotSaved{:},Var2Data_wFluo('03-08-20', '5mT',1.16, 'M1','DictyAx2_BSA30',4504,...
        MatfileFolder)};

    %% RPE1
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-04-19', '5mT',0.56, 'M1','RPE1',4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-04-19', '5mT',0.56, 'M1','RPE1',4438,...
        MatfileFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-04-19', '5mT',0.56, 'M1','RPE1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-04-19', '5mT',0.56, 'M2','RPE1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-04-19', '5mT',0.56, 'M3','RPE1',4438,...
        MatfileFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-04-19', '5mT',0.56, 'M1','RPE1',4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-04-19', '5mT',0.56, 'M2','RPE1',4438,...
        MatfileFolder)};
    
    %% From Paula : 
    % removefromdata : enleve les donn�es list� dans le dossier V2D
%j'enl�ve d abord les cellules mortes 100419 125 ?
RemoveFromData(MatfileFolder,{'11-04-19_M3_P1_C5_5mT','11-04-19_M3_P1_C6_5mT','11-04-19_M3_P1_C7_5mT','11-04-19_M3_P1_C8_5mT'},'RPE1')
%j'enl�ve a cause de problemes techniques et des cas bizzares
RemoveFromData(MatfileFolder,{'11-04-19_M1_P1_C3_5mT','11-04-19_M1_P1_C4_5mT','11-04-19_M1_P1_C9_5mT',...
    '11-04-19_M2_P1_C1_5mT', '11-04-19_M2_P1_C4-1_5mT', '11-04-19_M2_P1_C3_5mT','11-04-19_M2_P1_C5_5mT','10-04-19_M1_P1_C3_5mT',...
    '10-04-19_M1_P1_C3_5mT', '10-04-19_M1_P1_C5_5mT', '10-04-19_M1_P1_C6_5mT', '10-04-19_M1_P2_C4_5mT'},'RPE1')
%j'enl�ve d abord les cellules qui blebbent trop
RemoveFromData(MatfileFolder,{'17-04-19_M1_P1_C3_5mT','11-04-19_M3_P1_C1_5mT'},'RPE1')
       
    %% Inside Datas
    
    %%%% Manip etalonnage paula et vieux 0%
    
        NotSaved = {NotSaved{:},Var2Data_wFluo('22-01-18', '5mT',0.56,'M1','DC_Inside_0pc' ,4415,...
        MatfileFolder)};    
        NotSaved = {NotSaved{:},Var2Data_wFluo('22-01-18', '5mT',0.56,'M2','DC_Inside_0pc' ,4415,...
        MatfileFolder)};
    
        NotSaved = {NotSaved{:},Var2Data_wFluo('23-01-18', '5mT',0.56,'M1','DC_Inside_0pc' ,4415,...
        MatfileFolder)};    
        NotSaved = {NotSaved{:},Var2Data_wFluo('23-01-18', '5mT',0.56,'M2','DC_Inside_0pc' ,4415,...
        MatfileFolder)};
    
    
    
        NotSaved = {NotSaved{:},Var2Data_wFluo('25-03-19', '5mT',0.56,'M1','DC_Inside_0pc' ,4438,...
        MatfileFolder)};  
        NotSaved = {NotSaved{:},Var2Data_wFluo('25-03-19', '5mT',0.56,'M2','DC_Inside_5pc' ,4438,...
        MatfileFolder)};  
        NotSaved = {NotSaved{:},Var2Data_wFluo('25-03-19', '5mT',0.56,'M3','DC_Inside_10pc' ,4438,...
        MatfileFolder)};  
        NotSaved = {NotSaved{:},Var2Data_wFluo('25-03-19', '5mT',0.56,'M4','DC_Inside_12-5pc' ,4438,...
        MatfileFolder)};  
    
        NotSaved = {NotSaved{:},Var2Data_wFluo('26-03-19', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        MatfileFolder)};  
        NotSaved = {NotSaved{:},Var2Data_wFluo('26-03-19', '5mT',0.56,'M2','DC_Inside_15pc' ,4438,...
        MatfileFolder)};  
        NotSaved = {NotSaved{:},Var2Data_wFluo('26-03-19', '5mT',0.56,'M3','DC_Inside_10pc' ,4438,...
        MatfileFolder)};  
        NotSaved = {NotSaved{:},Var2Data_wFluo('26-03-19', '5mT',0.56,'M4','DC_Inside_12-5pc' ,4438,...
        MatfileFolder)};  
    
    %%%% 7.5pc
    % Culture 1
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56,'M2','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18', '5mT',0.56,'M2','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    
    % Culture 2
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18', '5mT',0.56,'M2','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18', '5mT',0.56,'M1','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    % Culture 3
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18', '5mT',0.56,'M1','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18', '5mT',0.56,'M1','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    % Culture 4
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
%     
    NotSaved = {NotSaved{:},Var2Data_wFluo('03-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('03-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    % Culture 5
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18', '5mT',0.56,'M1','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    % Culture 6
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18', '5mT',0.56,'M1','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18', '5mT',0.56,'M1','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    % Culture 7
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18', '5mT',0.56,'M2','DC_Inside_DMSO' ,4415,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('31-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('31-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        MatfileFolder)};
    
    % Arpin 1
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4438,...
        MatfileFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4438,...
        MatfileFolder)};
    
    
    % Culture 8
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        MatfileFolder)};
    
    % Culture 12
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        MatfileFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4438,...
        MatfileFolder)};
    
    % culture 2019
    
       NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19', '5mT',0.56,'M2','DC_Inside_DMSO' ,4438,...
        MatfileFolder)};
    
       NotSaved = {NotSaved{:},Var2Data_wFluo('16-05-19', '5mT',0.56,'M1','DC_Inside_DMSO' ,4438,...
        MatfileFolder)};
    
       NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-19', '5mT',0.56,'M1','DC_Inside_DMSO' ,4438,...
        MatfileFolder)};
    
      %% Manip 5 - 3 - 1 mT
    
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '5mT',1.1,'M1','DC_Inside_5mT' ,4438,...
        MatfileFolder)};
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '3mT',1.1,'M1','DC_Inside_3mT' ,4438,...
        MatfileFolder)};
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '1mT',1.1,'M1','DC_Inside_1mT' ,4438,...
        MatfileFolder)};
      
        NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '5mT',1.1,'M2','DC_Inside_5mT' ,4438,...
        MatfileFolder)};
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '3mT',1.1,'M2','DC_Inside_3mT' ,4438,...
        MatfileFolder)};
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '1mT',1.1,'M2','DC_Inside_1mT' ,4438,...
        MatfileFolder)};
       
    %% autre manips
       
    NotSaved = {NotSaved{:},Var2Data_wFluo('21-12-18', '5mT',0.56,'M1','DC_Nase_Inside' ,4438,...
        MatfileFolder)};
    
    % culture 14
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-04-19', '5mT',0.56,'M1','DC_TrueOutside' ,4438,...
        MatfileFolder)};
    
    % culture 15
    NotSaved = {NotSaved{:},Var2Data_wFluo('13-05-19', '5mT',0.56,'M2','DC_TrueOutside' ,4438,...
        MatfileFolder)};
    
    NotSaved(find(strcmp(NotSaved,'empty'))) = []
    
    
    RemoveFromData(MatfileFolder,{'17-09-18_M2_P2_C4_5mT','10-07-18_M1_P2_C1_5mT','20-12-18_M1_P2_C4_5mT'},'Inside')
    RemoveFromData(MatfileFolder,{'02-04-19_M1_P1_C1-3_5mT'},'TrueOutside')
      
    
    RemoveFromData(MatfileFolder,{'15-05-19_M1_P1_C5-1_5mT'},'Vim')
   
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data2Peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
     
    %% Nouvelle manips dictys
    
    
    Data2Peaks('DictyAx2_DMSO',PLOTDP,MatfileFolder,FigureFolder);
    Data2Peaks('DictyAx2_LatA',PLOTDP,MatfileFolder,FigureFolder);    
    
    Data2Peaks('DictyAx2_BSA30',PLOTDP,MatfileFolder,FigureFolder); 
    
    %% Data principale

Data2Peaks('RPE1',PLOTDP,MatfileFolder,FigureFolder);

Data2Peaks('DCLA_DMSO',PLOTDP,MatfileFolder,FigureFolder);
    
Data2Peaks('DCLA_NoDrug',PLOTDP,MatfileFolder,FigureFolder);
    
Data2Peaks('DCLA_LatA500',PLOTDP,MatfileFolder,FigureFolder);

Data2Peaks('DCLA_Blebbi',PLOTDP,MatfileFolder,FigureFolder);
    
Data2Peaks('DCLA_SMIFH2',PLOTDP,MatfileFolder,FigureFolder);
    
Data2Peaks('DCLA_CK666',PLOTDP,MatfileFolder,FigureFolder);
    
Data2Peaks('DC_Inside_DMSO',PLOTDP,MatfileFolder,FigureFolder);

Data2Peaks('DC_TrueOutside',PLOTDP,MatfileFolder,FigureFolder);

Data2Peaks('DCLA_CalA1',PLOTDP,MatfileFolder,FigureFolder);

Data2Peaks('DCLA_CalA3',PLOTDP,MatfileFolder,FigureFolder);

Data2Peaks('DCLA_CalA10',PLOTDP,MatfileFolder,FigureFolder);
    
    
Data2Peaks('DCLA_ML141',PLOTDP,MatfileFolder,FigureFolder);

    
    Data2Peaks('DCLA_Y27',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_Nase',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_LatNase',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_LatTryp',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_MyoII-WT',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_MyoII-KO',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_ArpinWT',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_ArpinKO',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_ArpC4WT',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_ArpC4KO',PLOTDP,MatfileFolder,FigureFolder);
 
    Data2Peaks('DC_Nase_Inside',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_SiVim',PLOTDP,MatfileFolder,FigureFolder);
    Data2Peaks('DCLA_SiCtrl',PLOTDP,MatfileFolder,FigureFolder);
    
    Data2Peaks('DCLA_SiVim+LatA',PLOTDP,MatfileFolder,FigureFolder);
    Data2Peaks('DCLA_SiCtrl+LatA',PLOTDP,MatfileFolder,FigureFolder);
    
    
    
    Data2Peaks('DCLA_Sep',PLOTDP,MatfileFolder,FigureFolder);
    
    
    Data2Peaks('Dicty_WT',1,MatfileFolder,FigureFolder);
    
    Data2Peaks('Dicty_Div',PLOTDP,MatfileFolder,FigureFolder);


Data2Peaks('DCLA_WT',PLOTDP,MatfileFolder,FigureFolder);   
    
    
Data2Endo('DCLA_Endo',PLOTDP,MatfileFolder,FigureFolder);

 %% Inside
Data2Peaks('DC_Inside_0pc',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DC_Inside_5pc',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DC_Inside_7-5pc',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DC_Inside_10pc',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DC_Inside_12-5pc',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DC_Inside_15pc',PLOTDP,MatfileFolder,FigureFolder);

    
Data2Peaks('DC_Inside_DMSO',PLOTDP,MatfileFolder,FigureFolder);

Data2Peaks('DC_Inside_1mT',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DC_Inside_3mT',PLOTDP,MatfileFolder,FigureFolder);
Data2Peaks('DC_Inside_5mT',PLOTDP,MatfileFolder,FigureFolder);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Peaks2Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    %% Nouvelles manips dicty


    Peaks2Stats('DictyAx2_DMSO',PLOTPS,alignlvl,MatfileFolder,FigureFolder); 
    
    Peaks2Stats('DictyAx2_LatA',PLOTPS,alignlvl,MatfileFolder,FigureFolder); 
    
    Peaks2Stats('DictyAx2_BSA30',PLOTPS,alignlvl,MatfileFolder,FigureFolder); 
    
%% Data principale

    Peaks2Stats('RPE1',PLOTPS,alignlvl,MatfileFolder,FigureFolder); 
    
	Peaks2Stats('DCLA_DMSO',PLOTPS,alignlvl,MatfileFolder,FigureFolder); 
    
    Peaks2Stats('DCLA_NoDrug',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_Blebbi',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_SMIFH2',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_CK666',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_CalA1',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_CalA3',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_CalA10',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_LatA500',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_ML141',PLOTPS,alignlvl,MatfileFolder,FigureFolder);    
    
    Peaks2Stats('DC_TrueOutside',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    

    
   Peaks2Stats('DCLA_ArpinWT',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_ArpinKO',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
     
   Peaks2Stats('DCLA_ArpC4WT',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_ArpC4KO',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    
    Peaks2Stats('DCLA_Y27',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_MyoII-WT',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_MyoII-KO',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_LatNase',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_LatTryp',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_Nase',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DC_Nase_Inside',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_SiVim',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_SiCtrl',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_SiVim+LatA',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DCLA_SiCtrl+LatA',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    
    Peaks2Stats('DictyAx2_DMSO',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DictyAx2_LatA',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    %% Inside     
    Peaks2Stats('DC_Inside_0pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DC_Inside_5pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DC_Inside_7-5pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DC_Inside_10pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DC_Inside_12-5pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DC_Inside_15pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DC_Inside_1mT',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DC_Inside_3mT',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    Peaks2Stats('DC_Inside_5mT',PLOTPS,alignlvl,MatfileFolder,FigureFolder);
    
    
    Peaks2Stats('DC_Inside_DMSO',PLOTPS,alignlvl,MatfileFolder,FigureFolder);

  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stats2PLots %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    
    CORT = 0;
    PEAKS    = 0;
    ACTI  = 0;
    DIST   = 0;
    TIME = 1;
    
    Stats2Plots({'DictyAx2_DMSO','DictyAx2_LatA'},{Cdm,Cb},...
        'DictyAx2LatA',CORT,DIST,PEAKS,ACTI,TIME,alignlvl,MatfileFolder,FigureFolder)
