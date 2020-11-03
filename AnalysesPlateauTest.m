close all;
clear all;
clc;

%% Paths
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';
RawdataFolder = 'D:\Data\Raw\';
MatfileFolder = 'D:\Data\MatFile';
DropboxDataFolder = 'C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Données';
FigureFolder = 'D:\Data\Figures';



%% Parameters

% Couleurs par drogues
Cb5  = [255 128 0  ]./255;
Cy5 = [200 200  0  ]./255;
Cp5  = [210 160 0  ]./255;

Cca5 = [0 128 255]./255;

Cs5  = [0   255  90]./255;

Cm5  = [200 100 80  ]./255;
Cc5  = [255 20   0  ]./255;

Cd5  = [160 160 160]./255;

Caw5  = [140 140 140]./255;
Cak5 = [0 255 255]./255;


Cl5 = [110 50 120]./255;
CL5 = [90 40 70]./255;
Cin = [0 50 150]./255;
Co = [150 50 0]./255;
% pour r2v
AUTO = 0; % automatic analysis mode (skip when user input is needed)

% pour v2d
NotSaved = {}; % gather a list of non saved data (because of sorting criteria)

% pour d2p
PLOTDP = 0; % plot curve with cortex level
cortlvl = 0.9; % cortex level

% pour p2s
PLOTPS = 0;

% pour s2p
BAR    = 1;
PROBA  = 1;
DIST   = 1;

%% Analysis

if 0 % R2V
    %% Nouvelles manips Inhib
    
    
    % Culture 1
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DCLA_CK666'    ,1:3 ,500,147,'Kimo_DCmed-7pc_170518' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M2','DCLA_Ctrl'     ,1,500,147,'Kimo_DCmed-7pc_170518' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M2','DCLA_LatA50'   ,2,500,147,'Kimo_DCmed-7pc_170518' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M1','DCLA_Blebbi'   ,1,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DCLA_Ctrl'     ,1,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DCLA_LatA500'  ,2,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 2
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M1','DCLA_Blebbi'   ,1:3 ,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DCLA_Ctrl'     ,1,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DCLA_LatA50'   ,2,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DCLA_LatA50'   ,3,0  ,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DCLA_Ctrl'     ,1,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DCLA_LatA500'  ,2,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M2','DCLA_CK666'    ,1:3 ,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 3
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M1','DCLA_Ctrl'     ,1,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M1','DCLA_LatA50'   ,2,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M2','DCLA_Blebbi'   ,1:3 ,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DCLA_Ctrl'     ,1,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DCLA_LatA500'  ,2,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M2','DCLA_CK666'  ,1:3 ,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 4
    Res2Var_wFluo_multiZ('02-07-18', '5mT', 'M1','DCLA_CTRLF'    ,1,500,46,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('03-07-18', '5mT','M1','DCLA_CTRLF'    ,2,500,49,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('03-07-18', '5mT','M2','DCLA_SMIFH2F'  ,1:3 ,500,49,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 5
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DCLA_SMIFH2'   ,1,500,48,'Kimo_DCmed-7.5pc_120718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DCLA_Ctrl'     ,2,500,48,'Kimo_DCmed-7.5pc_120718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M2','DCLA_ML141'    ,1,500,48,'Kimo_DCmed-7.5pc_120718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M2','DCLA_SMIFH2'   ,2,500,48,'Kimo_DCmed-7.5pc_120718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 6
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DCLA_ML141'    ,1,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DCLA_Ctrl'     ,2,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M2','DCLA_ML141'    ,1,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DCLA_Ctrl'     ,1,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DCLA_SMIFH2'   ,1,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DCLA_ML141'    ,2,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('20-07-18', '5mT', 'M1','DCLA_Blebbi'   ,1,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 7
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M1','DCLA_ML141'    ,1,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DCLA_Blebbi'   ,1,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DCLA_Ctrl'     ,2,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M1','DCLA_LatA500'  ,1,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M1','DCLA_CK666'    ,2,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M2','DCLA_CK20SM10' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Arpin 1
    Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M1','DCLA_ArpinWT' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M2','DCLA_ArpinKO' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M3','DCLA_ArpinKO' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Arpin 2
    Res2Var_wFluo_multiZ('04-10-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('04-10-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    Res2Var_wFluo_multiZ('05-10-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('05-10-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 8
    Res2Var_wFluo_multiZ('29-10-18', '5mT', 'M1','DCLA_CalA25',1 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('29-10-18', '5mT', 'M1','DCLA_LatA2',2 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('29-10-18', '5mT', 'M2','DCLA_Para-B',1 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('29-10-18', '5mT', 'M2','DCLA_LatA2',2 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('30-10-18', '5mT', 'M1','DCLA_Para-B',1 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-10-18', '5mT', 'M2','DCLA_CalA25',1 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 9
    Res2Var_wFluo_multiZ('15-11-18', '5mT', 'M1','DCLA_Y27',1:3 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-11-18', '5mT', 'M2','DCLA_Ctrl',1 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-11-18', '5mT', 'M2','DCLA_LatA2',2 ,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-11-18', '5mT', 'M2','DCLA_LatA2',3 ,500,49,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 10
    Res2Var_wFluo_multiZ('12-12-18', '5mT', 'M1','DCLA_Ctrl',1:3 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('12-12-18', '5mT', 'M2','DCLA_Nase',1:3 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 11
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M1','DCLA_ArpC4KO',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M1','DCLA_Y27',2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M2','DCLA_ArpC4WT',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M1','DCLA_Y27',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M1','DCLA_SMIFH2',2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M2','DCLA_Ctrl',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M2','DCLA_Nase',2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 12
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M1','DCLA_SMIFH2',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M1','DCLA_Y27',2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DCLA_SMIFH2',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DCLA_LatA500',2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('21-12-18', '5mT','M1','DCLA_LatNase',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Fibro 1
    Res2Var_wFluo_multiZ('30-01-19', '5mT','M1','MFT16_NoVim',[1 2] ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-01-19', '5mT','M2','MFT16_Vim',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Fibro 2
    Res2Var_wFluo_multiZ('15-02-19', '5mT','M1','MFT16_Vim+Lat',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-02-19', '5mT','M2','MFT16_NoVim+Lat',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-02-19', '5mT','M3','MFT16_Vim+Lat',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-02-19', '5mT','M4','MFT16_NoVim+Lat',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % MyoII KO 1
    Res2Var_wFluo_multiZ('11-03-19', '5mT','M1','DCLA_MyoII-WT',1:2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('11-03-19', '5mT','M2','DCLA_MyoII-KO',1:2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('12-03-19', '5mT','M1','DCLA_MyoII-KO',1:2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('12-03-19', '5mT','M2','DCLA_MyoII-WT',1:2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 13
    Res2Var_wFluo_multiZ('20-03-19', '5mT','M1','DCLA_LatTryp',1:2 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('21-03-19', '5mT','M1','DCLA_LatTryp',1 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 14
        Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M1','DCLA_CalA10',1 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M1','DCLA_CalA1',2 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M2','DCLA_CalA1',1 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M2','DCLA_CalA10',2 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)

    Res2Var_wFluo_multiZ('02-04-19', '5mT', 'M2','DCLA_CalA1',1 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('02-04-19', '5mT', 'M2','DCLA_CalA10',2 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    %% Inside/Outside Datas
    
    Res2Var_wFluo_multiZ_Inside('22-01-18', '5mT', 'M1','DCOld_Inside_0pc',1:2,1000,46,'Kimo_Old' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('22-01-18', '5mT', 'M2','DCOld_Inside_0pc',1:2,1000,46,'Kimo_Old' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('23-01-18', '5mT', 'M1','DCOld_Inside_0pc',1:2,1000,46,'Kimo_Old' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('23-01-18', '5mT', 'M2','DCOld_Inside_0pc',1:2,1000,46,'Kimo_Old' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Exp etalonnage
    Res2Var_wFluo_Inside('07-04-18', '0pc', 'M1','DC_Inside_0pc',147,'Kimo_DCmed-0pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_Inside('07-04-18', '5pc', 'M1','DC_Inside_5pc',147,'Kimo_DCmed-5pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_Inside('07-04-18', '7-5pc', 'M1','DC_Inside_7-5pc',147,'Kimo_DCmed-7pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    Res2Var_wFluo_Inside('24-05-18', '0pc', 'M1','DC_Inside_0pc',147,'Kimo_DCmed-0pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_Inside('24-05-18', '3pc', 'M1','DC_Inside_3pc',147,'Kimo_DCmed-3pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_Inside('24-05-18', '5pc', 'M1','DC_Inside_5pc',147,'Kimo_DCmed-5pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_Inside('24-05-18', '6pc', 'M1','DC_Inside_6pc',147,'Kimo_DCmed-6pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_Inside('24-05-18', '7pc', 'M1','DC_Inside_7pc',147,'Kimo_DCmed-7pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_Inside('24-05-18', '8pc', 'M1','DC_Inside_8pc',147,'Kimo_DCmed-8pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_Inside('24-05-18', '10pc', 'M1','DC_Inside_10pc',147,'Kimo_DCmed-0pc_170518' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % exp etalonnage 2 25-26/03/19
    Res2Var_wFluo_multiZ_Inside('25-03-19', '5mT', 'M1','DCLA_Inside_0pc',1,500,48,'Deptho_DCmed-7.5pc_150319' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('25-03-19', '5mT', 'M2','DCLA_Inside_5pc',1,500,48,'Deptho_DCmed-7.5pc_150319' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('25-03-19', '5mT', 'M3','DCLA_Inside_10pc',1,500,48,'Deptho_DCmed-7.5pc_150319' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('25-03-19', '5mT', 'M4','DCLA_Inside_12-5pc',1,500,48,'Deptho_DCmed-7.5pc_150319' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('26-03-19', '5mT', 'M1','DCLA_Inside_7-5pc',1,500,48,'Deptho_DCmed-7.5pc_150319' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('26-03-19', '5mT', 'M2','DCLA_Inside_15pc',1,500,48,'Deptho_DCmed-7.5pc_150319' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('26-03-19', '5mT', 'M3','DCLA_Inside_10pc',1,500,48,'Deptho_DCmed-7.5pc_150319' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('26-03-19', '5mT', 'M4','DCLA_Inside_12-5pc',1,500,48,'Deptho_DCmed-7.5pc_150319' ,0,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    
    % Culture 1
    Res2Var_wFluo_multiZ_Inside('28-05-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,147,'Kimo_DCmed-7pc_170518' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('28-05-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,147,'Kimo_DCmed-7pc_170518' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('29-05-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 2
    Res2Var_wFluo_multiZ_Inside('04-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('04-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('05-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('05-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 3
    Res2Var_wFluo_multiZ_Inside('07-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('07-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('08-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('08-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,45,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 4
    Res2Var_wFluo_multiZ_Inside('02-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,46,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('02-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,46,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('03-07-18', '5mT','M1','DC_Inside_7-5pc'    ,1:3,500,49,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('03-07-18', '5mT','M2','DC_Inside_7-5pc'    ,1:3,500,49,'Kimo_DCmed-7pc_170518'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 5
    Res2Var_wFluo_multiZ_Inside('10-07-18', '5mT','M1','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_120718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('10-07-18', '5mT','M2','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_120718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 6
    Res2Var_wFluo_multiZ_Inside('18-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('18-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('19-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('19-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('20-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_180718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 7
    Res2Var_wFluo_multiZ_Inside('30-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('31-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('31-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_310718',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Arpin 1
    Res2Var_wFluo_multiZ_Inside('17-09-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('17-09-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('18-09-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('18-09-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Arpin 2
    Res2Var_wFluo_multiZ_Inside('04-10-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('04-10-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,500,48,'Kimo_DCmed-7.5pc_280818',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 8
    Res2Var_wFluo_multiZ_Inside('29-10-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,500,48,'Deptho_DCmed-7.5pc_301018' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 10
    Res2Var_wFluo_multiZ_Inside('12-12-18', '5mT', 'M2','DC_Nase_Inside',1:3 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 11
    Res2Var_wFluo_multiZ_Outside('17-12-18', '5mT', 'M1','DC_Outside',1:2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Outside('17-12-18', '5mT', 'M2','DC_Outside',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Outside('18-12-18', '5mT', 'M2','DC_Outside',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 12
    Res2Var_wFluo_multiZ_Inside('20-12-18', '5mT','M1','DC_Inside_7-5pc',2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('20-12-18', '5mT','M2','DC_Inside_7-5pc',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ_Inside('20-12-18', '5mT','M2','DC_Inside_7-5pc',2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Outside('20-12-18', '5mT', 'M1','DC_Outside',1:2 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ_Inside('21-12-18', '5mT','M1','DC_Nase_Inside',1 ,500,48,'Deptho_DCmed-7.5pc_121218' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    
    % Culture 14

    Res2Var_wFluo_multiZ_Outside('02-04-19', '5mT', 'M1','DC_TrueOutside',1 ,500,48,'Deptho_DCmed-7.5pc_150319' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)

end


if 0 % V2D
    %% Nouvelles manips Inhib
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('28-05-18','5mT',0.56,'M1','DCLA_CK666',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('28-05-18','5mT',0.56,'M2','DCLA_Ctrl',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('28-05-18','5mT',0.56,'M2','DCLA_LatA50',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-05-18','5mT',0.56,'M1','DCLA_Blebbi',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-05-18','5mT',0.56,'M2','DCLA_Ctrl',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-05-18','5mT',0.56,'M2','DCLA_LatA500',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-06-18','5mT',0.56,'M1','DCLA_Blebbi',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-06-18','5mT',0.56,'M2','DCLA_Ctrl',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-06-18','5mT',0.56,'M2','DCLA_LatA50',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('05-06-18','5mT',0.56,'M1','DCLA_Ctrl',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('05-06-18','5mT',0.56,'M1','DCLA_LatA500',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('05-06-18','5mT',0.56,'M2','DCLA_CK666',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('07-06-18','5mT',0.56,'M1','DCLA_Ctrl',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('07-06-18','5mT',0.56,'M1','DCLA_LatA50',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('07-06-18','5mT',0.56,'M2','DCLA_Blebbi',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('08-06-18','5mT',0.56,'M1','DCLA_Ctrl',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('08-06-18','5mT',0.56,'M1','DCLA_LatA500',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('08-06-18','5mT',0.56,'M2','DCLA_CK666',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('02-07-18','5mT',0.56,'M1','DCLA_CtrlF',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('03-07-18','5mT',0.56,'M1','DCLA_CtrlF',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('03-07-18','5mT',0.56,'M2','DCLA_SMIFH2F',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('10-07-18','5mT',0.56,'M1','DCLA_Ctrl',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('10-07-18','5mT',0.56,'M1','DCLA_SMIFH2',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('10-07-18','5mT',0.56,'M2','DCLA_ML141',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('10-07-18','5mT',0.56,'M2','DCLA_SMIFH2',4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-07-18','5mT',0.56,'M1','DCLA_ML141' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-07-18','5mT',0.56,'M1','DCLA_Ctrl'  ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-07-18','5mT',0.56,'M2','DCLA_ML141' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('19-07-18','5mT',0.56,'M1','DCLA_Ctrl'   ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('19-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('19-07-18','5mT',0.56,'M2','DCLA_SMIFH2' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('19-07-18','5mT',0.56,'M2','DCLA_ML141'  ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-07-18','5mT',0.56,'M1','DCLA_Blebbi' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-07-18','5mT',0.56,'M1','DCLA_ML141'  ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-07-18','5mT',0.56,'M2','DCLA_Blebbi' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-07-18','5mT',0.56,'M2','DCLA_Ctrl'   ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('31-07-18','5mT',0.56,'M1','DCLA_LatA500' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('31-07-18','5mT',0.56,'M1','DCLA_CK666'   ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % nouvelles billess
%     
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-09-18','5mT',0.56,'M1','DCLA_ArpinWT' ,4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-09-18','5mT',0.56,'M2','DCLA_ArpinKO',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-09-18','5mT',0.56,'M1','DCLA_ArpinKO',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-09-18','5mT',0.56,'M2','DCLA_ArpinWT',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-09-18','5mT',0.56,'M3','DCLA_ArpinKO',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-10-18','5mT',0.56,'M1','DCLA_ArpinKO',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-10-18','5mT',0.56,'M2','DCLA_ArpinWT',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('05-10-18','5mT',0.56,'M1','DCLA_ArpinKO',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('05-10-18','5mT',0.56,'M2','DCLA_ArpinWT',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-10-18','5mT',0.56,'M1','DCLA_CalA25',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-10-18','5mT',0.56,'M1','DCLA_LatA2',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-10-18','5mT',0.56,'M2','DCLA_Para-B',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-10-18','5mT',0.56,'M2','DCLA_LatA2',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-10-18','5mT',0.56,'M1','DCLA_Para-B',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-10-18','5mT',0.56,'M2','DCLA_CalA25',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('15-11-18','5mT',0.56,'M1','DCLA_Y27',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('15-11-18','5mT',0.56,'M2','DCLA_Ctrl',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('15-11-18','5mT',0.56,'M2','DCLA_LatA2',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('12-12-18','5mT',0.56,'M1','DCLA_Ctrl',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('12-12-18','5mT',0.56,'M2','DCLA_Nase',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-12-18','5mT',0.56,'M1','DCLA_ArpC4KO',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-12-18','5mT',0.56,'M1','DCLA_Y27',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-12-18','5mT',0.56,'M2','DCLA_ArpC4WT',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-12-18','5mT',0.56,'M1','DCLA_Y27',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-12-18','5mT',0.56,'M1','DCLA_SMIFH2',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-12-18','5mT',0.56,'M2','DCLA_Ctrl',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-12-18','5mT',0.56,'M1','DCLA_Y27',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-12-18','5mT',0.56,'M1','DCLA_SMIFH2',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-12-18','5mT',0.56,'M2','DCLA_SMIFH2',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-12-18','5mT',0.56,'M2','DCLA_LatA500',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('21-12-18','5mT',0.56,'M1','DCLA_LatNase',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
%     % Fibroblast
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-01-19','5mT',0.56,'M1','MFT16_NoVim',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-01-19','5mT',0.56,'M2','MFT16_Vim',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('15-02-19','5mT',0.56,'M1','MFT16_Vim+Lat',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('15-02-19','5mT',0.56,'M2','MFT16_NoVim+Lat',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('15-02-19','5mT',0.56,'M3','MFT16_Vim+Lat',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('15-02-19','5mT',0.56,'M4','MFT16_NoVim+Lat',4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
%     RemoveFromData({'15-02-19_M1_P1_C1_5mT','15-02-19_M1_P1_C2_5mT','15-02-19_M1_P1_C9_5mT','15-02-19_M3_P1_C2_5mT','15-02-19_M3_P1_C5_5mT',...
%         '15-02-19_M2_P1_C5_5mT','15-02-19_M2_P1_C6-1_5mT','15-02-19_M2_P1_C10_5mT',...
%         '15-02-19_M4_P1_C3_5mT','15-02-19_M4_P1_C4_5mT','15-02-19_M4_P1_C5_5mT'...
%         },'MFT')
    
    % KO myoII
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('11-03-19','5mT',0.56,'M1','DCLA_MyoII-WT',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('11-03-19','5mT',0.56,'M2','DCLA_MyoII-KO',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('12-03-19','5mT',0.56,'M1','DCLA_MyoII-KO',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('12-03-19','5mT',0.56,'M2','DCLA_MyoII-WT',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    
    % culture 14
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-03-19','5mT',0.56,'M1','DCLA_LatTryp',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('21-03-19','5mT',0.56,'M1','DCLA_LatTryp',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('01-04-19','5mT',0.56,'M1','DCLA_CalA1',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('01-04-19','5mT',0.56,'M1','DCLA_CalA10',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('01-04-19','5mT',0.56,'M2','DCLA_CalA1',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('01-04-19','5mT',0.56,'M2','DCLA_CalA10',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};

    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('02-04-19','5mT',0.56,'M2','DCLA_CalA1',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('02-04-19','5mT',0.56,'M2','DCLA_CalA10',4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};


    
    %% Inside Datas
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('22-01-18', '5mT',0.56,'M1','DCOld_Inside_0pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('22-01-18', '5mT',0.56,'M2','DCOld_Inside_0pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('23-01-18', '5mT',0.56,'M1','DCOld_Inside_0pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('23-01-18', '5mT',0.56,'M2','DCOld_Inside_0pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('07-04-18', '5mT',0.56,'M1','DC_Inside_0pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
%     
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('24-05-18', '5mT',0.56,'M1','DC_Inside_0pc' ,4432,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('24-05-18', '5mT',0.56,'M1','DC_Inside_3pc' ,4432,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('24-05-18', '5mT',0.56,'M1','DC_Inside_5pc' ,4432,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('24-05-18', '5mT',0.56,'M1','DC_Inside_6pc' ,4432,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('24-05-18', '5mT',0.56,'M1','DC_Inside_7pc' ,4432,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('24-05-18', '5mT',0.56,'M1','DC_Inside_8pc' ,4432,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('24-05-18', '5mT',0.56,'M1','DC_Inside_10pc' ,4432,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    %
    % Culture 1
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('28-05-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('28-05-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-05-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % Culture 2
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('05-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('05-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % Culture 3
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('07-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('07-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('08-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('08-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % Culture 4
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('02-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('02-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('03-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('03-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % Culture 5
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('10-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('10-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % Culture 6
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('19-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('19-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % Culture 7
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('30-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('31-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('31-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4432,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % Arpin 1
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-09-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-09-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-09-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-09-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    % Arpin 2
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-10-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('04-10-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    
    % Culture 8
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('29-10-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    
%     % Culture 10
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('12-12-18', '5mT',0.56,'M2','DC_Nase_Inside' ,4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
    
%     % culture 11
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-12-18', '5mT',0.56,'M1','DC_Outside' ,4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('17-12-18', '5mT',0.56,'M2','DC_Outside' ,4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
%     NotSaved = {[[NotSaved{:}],Var2Data_wFluo('18-12-18', '5mT',0.56,'M2','DC_Outside' ,4447,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)]};
%     
    
    % Culture 12
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-12-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-12-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('20-12-18', '5mT',0.56,'M1','DC_Outside' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('21-12-18', '5mT',0.56,'M1','DC_Nase_Inside' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('25-03-19', '5mT',0.56,'M1','DCLA_Inside_0pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('25-03-19', '5mT',0.56,'M2','DCLA_Inside_5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('25-03-19', '5mT',0.56,'M3','DCLA_Inside_10pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('25-03-19', '5mT',0.56,'M4','DCLA_Inside_12-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('26-03-19', '5mT',0.56,'M1','DCLA_Inside_7-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('26-03-19', '5mT',0.56,'M2','DCLA_Inside_15pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('26-03-19', '5mT',0.56,'M3','DCLA_Inside_10pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('26-03-19', '5mT',0.56,'M4','DCLA_Inside_12-5pc' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    
    % culture 14
    NotSaved = {[[NotSaved{:}],Var2Data_wFluo('02-04-19', '5mT',0.56,'M1','DC_TrueOutside' ,4447,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)]};
    
    NotSaved(find(strcmp(NotSaved,'empty'))) = []
    
end


if 0 % D2P
    
    
    Data2Peaks({'M1','M2'},'DCLA_Ctrl',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2'},'DCLA_Blebbi',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2'},'DCLA_SMIFH2',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2'},'DCLA_CK666',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2'},'DCLA_LatA500',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
     Data2Peaks({'M1','M2'},'DCLA_CalA1',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
     
     Data2Peaks({'M1','M2'},'DCLA_CalA10',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
    
    
    Data2Peaks({'M1','M2','M3'}, 'DC_TrueOutside',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2','M3'}, 'DCOld_Inside_0pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
%     
%     Data2Peaks({'M1','M2','M3'}, 'DC_Inside_0pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Data2Peaks({'M1','M2','M3'}, 'DC_Inside_3pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Data2Peaks({'M1','M2','M3'}, 'DC_Inside_5pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Data2Peaks({'M1','M2','M3'}, 'DC_Inside_6pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Data2Peaks({'M1','M2','M3'}, 'DC_Inside_7pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Data2Peaks({'M1','M2','M3'}, 'DC_Inside_8pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Data2Peaks({'M1','M2','M3'}, 'DC_Inside_10pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
    Data2Peaks({'M1','M2','M3'}, 'DC_Inside_7-5pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
    Data2Peaks({'M1','M2','M3','M4'}, 'DCLA_Inside_0pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2','M3','M4'}, 'DCLA_Inside_5pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2','M3','M4'}, 'DCLA_Inside_7-5pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2','M3','M4'}, 'DCLA_Inside_10pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2','M3','M4'}, 'DCLA_Inside_12-5pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2','M3','M4'}, 'DCLA_Inside_15pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
    
    Data2Peaks({'M1','M2'},'DCLA_Y27',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2'},'DCLA_ML141',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
    
    
    Data2Peaks({'M1','M2','M3'},'DCLA_Nase',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1','M2','M3'}, 'DCLA_LatNase',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks({'M1'}, 'DCLA_LatTryp',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
%     
%     
%     Data2Peaks({'M1','M2'},'DCLA_CalA25',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     
%     Data2Peaks({'M1','M2'},'DCLA_MyoII-WT',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Data2Peaks({'M1','M2'},'DCLA_MyoII-KO',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Data2Peaks({'M2','M4'},'MFT16_NoVim+Lat',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Data2Peaks({'M1','M3'},'MFT16_Vim+Lat',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     
%     Data2Peaks({'M1','M2'},'MFT16_NoVim',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Data2Peaks({'M1','M2'},'MFT16_Vim',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
    
    
    %     Data2Peaks({'M1'}, 'DC_Inside_0pc',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    %     Data2Peaks({'M1','M2','M3'}, 'DC_Nase_Inside',PLOTDP,cortlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
end


if 0 % P2S
    
    Peaks2Stats('DCLA_Ctrl',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_CK666',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_LatA500',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_Blebbi',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_Y27',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_ML141',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_SMIFH2',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_CalA1',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_CalA10',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
%     Peaks2Stats('DCLA_CalA25',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
%     Peaks2Stats('DC_Outside',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DC_TrueOutside',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
%     Peaks2Stats('DCLA_MyoII-WT',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Peaks2Stats('DCLA_MyoII-KO',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_LatNase',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_LatTryp',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
    Peaks2Stats('DCOld_Inside_0pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Peaks2Stats('DC_Inside_0pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_3pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_5pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_6pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_7pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_8pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_10pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
    Peaks2Stats('DC_Inside_7-5pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
%     Peaks2Stats('DCLA_Inside_0pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DCLA_Inside_5pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DCLA_Inside_7-5pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DCLA_Inside_10pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DCLA_Inside_12-5pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DCLA_Inside_15pc',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
    
    %
    
    %     Peaks2Stats('DCLA_Nase',PLOTPS);
    
    %
%     Peaks2Stats('MFT16_NoVim+Lat',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Peaks2Stats('MFT16_Vim+Lat',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Peaks2Stats('MFT16_NoVim',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Peaks2Stats('MFT16_Vim',PLOTPS,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     
    
    
    %
    %     Peaks2Stats('DC_Nase_Inside',PLOTPS);
    
end


if 0 % Plots
    
    BAR    = 0;
    PROBA  = 1;
    DIST   = 0;
    

    
    
    Stats2Plots({'DCLA_Ctrl','DCLA_Blebbi','DCLA_CalA1','DCLA_CK666',...
        'DCLA_SMIFH2','DCLA_LatA500'},...
        'NewDCLA-Basic',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
     
    Stats2Plots({'DCLA_Ctrl','DCLA_Blebbi','DCLA_Y27','DCLA_CalA1','DCLA_CalA10'},...
        'NewDCLA-MyoFull',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_Ctrl','DCLA_Blebbi','DCLA_Y27'},...
        'NewDCLA-Myo',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_Ctrl','DCLA_Blebbi'},...
        'NewDCLA-Bleb',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_Ctrl','DCLA_CalA1'},...
        'NewDCLA-CalA1',DIST,BAR,PROBA,MatfileFolder,FigureFolder)   
    Stats2Plots({'DCLA_Ctrl','DCLA_CalA10'},...
        'NewDCLA-CalA10',DIST,BAR,PROBA,MatfileFolder,FigureFolder) 
    Stats2Plots({'DCLA_Ctrl','DCLA_CK666'},...
        'NewDCLA-Arp',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_Ctrl','DCLA_SMIFH2'},...
        'NewDCLA-SM',DIST,BAR,PROBA,MatfileFolder,FigureFolder)   
    
    
    Stats2Plots({'DCLA_Ctrl', 'DC_Inside_7-5pc','DC_TrueOutside'},...
        'NewDCLA-inout',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    
   Stats2Plots({'DCLA_Ctrl','DCLA_Blebbi','DCLA_CalA1'},...
        'NewDCLA-Fig5',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    
    
    Stats2Plots({'DCLA_Inside_0pc','DCLA_Inside_5pc','DCLA_Inside_7-5pc','DCLA_Inside_10pc',...
        'DCLA_Inside_12-5pc','DCLA_Inside_15pc'},...
        'NewDCLA-Insides',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    
    

    
    
    Stats2Plots({'DCLA_Ctrl','DCLA_Blebbi','DCLA_SMIFH2','DCLA_CK666'},...
        'NewDCLA-Fig3',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    

    

  
   
    

    
    Stats2Plots({'DCLA_MyoII-WT','DCLA_MyoII-KO'},...
        'NewDCLA-MyoKO',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    
    
    
    
    Stats2Plots({'DCLA_Ctrl','DCLA_LatA500','DCLA_LatTryp'},...
        'NewDCLA-LatA+tryp',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    

    
    Stats2Plots({'DCLA_Ctrl','DCLA_CK666'},...
        'NewDCLA-CK',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    
    
    
    Stats2Plots({'DCLA_Ctrl','DCLA_SMIFH2','DCLA_CK666'},...
        'NewDCLA-SMCK',DIST,BAR,PROBA,MatfileFolder,FigureFolder)

    
    Stats2Plots({'DCLA_Ctrl','DCLA_CK666','DCLA_Blebbi'},...
        'NewDCLA-CKB',DIST,BAR,PROBA,MatfileFolder,FigureFolder)
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data2CortexStats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %     edit FigureMaker.m
    
    FigureMaker(0,MatfileFolder,FigureFolder)
end
