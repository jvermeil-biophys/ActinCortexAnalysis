close all;
clear all;
clc;

%% Paths
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';
RawdataFolder = 'D:\Data\Raw';
MatfileFolder = 'D:\Data\MatFile';
DropboxDataFolder = 'C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Données';
FigureFolder = 'D:\Data\Figures';

set(0, 'defaultfigureposition', get(0, 'Screensize').*[20 20 0.98 0.9]);

%% Parameters

% Couleurs par drogues
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


ALL = 0;

% pour r2v
AUTO = 0; % automatic analysis mode (skip when user input is needed)

% pour v2d
NotSaved = {}; % gather a list of non saved data (because of sorting criteria)

% pour d2p
PLOTDP = 0; % plot curve with cortex level
alignlvl = 0.5; % for probcum

% pour p2s
PLOTPS = 0;

% pour s2p
PEAKS    = 1;
PROBA  = 1;
DIST   = 0;

%% Analysis

if ALL % R2V
    %% Nouvelles manips Inhib   
    
    % Culture 1
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DCLA_CK666'    ,1:3 ,'Results',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M2','DCLA_DMSO'     ,1,'Results',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M1','DCLA_Blebbi'   ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DCLA_LatA500'  ,2,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 2
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M1','DCLA_Blebbi'   ,1:3 ,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DCLA_LatA500'  ,2,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M2','DCLA_CK666'    ,1:3 ,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 3
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M1','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M2','DCLA_Blebbi'   ,1:3 ,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DCLA_DMSO'     ,1,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DCLA_LatA500'  ,2,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M2','DCLA_CK666'  ,1:3 ,'Results',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 5
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DCLA_SMIFH2'   ,1,'Results',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DCLA_DMSO'     ,2,'Results',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M2','DCLA_ML141'    ,1,'Results',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M2','DCLA_SMIFH2'   ,2,'Results',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 6
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DCLA_ML141'    ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DCLA_DMSO'     ,2,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M2','DCLA_ML141'    ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DCLA_DMSO'     ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DCLA_SMIFH2'   ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DCLA_ML141'    ,2,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('20-07-18', '5mT', 'M1','DCLA_Blebbi'   ,1,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,'Results',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 7
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M1','DCLA_ML141'    ,1,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M1','DCLA_SMIFH2'   ,2,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DCLA_Blebbi'   ,1,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DCLA_DMSO'     ,2,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M1','DCLA_LatA500'  ,1,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M1','DCLA_CK666'    ,2,'Results',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Arpin 1
%     Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M1','DCLA_ArpinWT' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M2','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     
%     Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M3','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Arpin 2
%     Res2Var_wFluo_multiZ('04-10-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     Res2Var_wFluo_multiZ('04-10-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     
%     
%     Res2Var_wFluo_multiZ('05-10-18', '5mT', 'M1','DCLA_ArpinKO' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     Res2Var_wFluo_multiZ('05-10-18', '5mT', 'M2','DCLA_ArpinWT' ,1:3 ,'Results',48,'28-08-18_Depthograph',AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    
    % Culture 9
    Res2Var_wFluo_multiZ('15-11-18', '5mT', 'M1','DCLA_Y27',1:3 ,'Results',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-11-18', '5mT', 'M2','DCLA_NoDrug',1 ,'Results',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    
    % Culture 10
    Res2Var_wFluo_multiZ('12-12-18', '5mT', 'M1','DCLA_NoDrug',1:3 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('12-12-18', '5mT', 'M2','DCLA_Nase',1:3 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 11
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M1','DCLA_ArpC4KO',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M1','DCLA_Y27',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M2','DCLA_ArpC4WT',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M1','DCLA_Y27',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M1','DCLA_SMIFH2',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M2','DCLA_NoDrug',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-12-18', '5mT','M2','DCLA_Nase',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 12
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M1','DCLA_SMIFH2',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M1','DCLA_Y27',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DCLA_SMIFH2',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DCLA_LatA500',2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('21-12-18', '5mT','M1','DCLA_LatNase',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % MyoII KO 1
    Res2Var_wFluo_multiZ('11-03-19', '5mT','M1','DCLA_MyoII-WT',1:2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('11-03-19', '5mT','M2','DCLA_MyoII-KO',1:2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('12-03-19', '5mT','M1','DCLA_MyoII-KO',1:2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('12-03-19', '5mT','M2','DCLA_MyoII-WT',1:2 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 13
    Res2Var_wFluo_multiZ('20-03-19', '5mT','M1','DCLA_LatTryp',1:2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('21-03-19', '5mT','M1','DCLA_LatTryp',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 14
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M1','DCLA_CalA10',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M1','DCLA_CalA1',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M2','DCLA_CalA1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('01-04-19', '5mT', 'M2','DCLA_CalA10',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('02-04-19', '5mT', 'M2','DCLA_CalA1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('02-04-19', '5mT', 'M2','DCLA_CalA10',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % SiRNA1
    Res2Var_wFluo_multiZ('29-04-19', '5mT', 'M1','DCLA_SiVim',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('29-04-19', '5mT', 'M1','DCLA_SiCtrl',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('30-04-19', '5mT', 'M1','DCLA_SiVim+LatA',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('30-04-19', '5mT', 'M1','DCLA_SiCtrl+LatA',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % culture 15
    Res2Var_wFluo_multiZ('13-05-19', '5mT', 'M1','DCLA_CalA3',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('13-05-19', '5mT', 'M1','DCLA_CalA1',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('13-05-19', '5mT', 'M3','DCLA_NoDrug',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M1','DCLA_CalA3',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M1','DCLA_CalA1',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M2','DCLA_DMSO',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M2','DCLA_NoDrug',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % culture 16 avec SiRNA2
    Res2Var_wFluo_multiZ('15-05-19', '5mT', 'M1','DCLA_SiVim+LatA',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-05-19', '5mT', 'M1','DCLA_SiCtrl+LatA',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-05-19', '5mT', 'M2','DCLA_SiVim',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('15-05-19', '5mT', 'M2','DCLA_SiCtrl',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    Res2Var_wFluo_multiZ('16-05-19', '5mT', 'M1','DCLA_SiCtrl',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('16-05-19', '5mT', 'M1','DCLA_DMSO',2 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % culture 17
        
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M1','DCLA_CalA1',1 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)        
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M1','DCLA_CalA3',2 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
        Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M2','DCLA_DMSO',1 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)        
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M2','DCLA_CalA3',2 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M1','DCLA_Sep',1:2 ,'Separation',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('25-06-19', '5mT', 'M2','DCLA_Sep',1 ,'Separation',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % culture 18
    Res2Var_wFluo_multiZ('02-07-19', '5mT', 'M1','DCLA_DMSO',1 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('02-07-19', '5mT', 'M2','DCLA_CalA1',1 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('02-07-19', '5mT', 'M2','DCLA_CalA3',2 ,'Results',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    %% endocytose et division
        Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DCLA_Endo' ,2 ,'Endo',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
        Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DCLA_Endo' ,1,'Endo',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
        Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M2','DCLA_Endo' ,1:2,'Endo',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
        Res2Var_wFluo_multiZ('17-12-18', '5mT', 'M2','DCLA_Endo',1 ,'Endo',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
        Res2Var_wFluo_multiZ('12-03-19', '5mT', 'M1','DCLA_Endo' ,1 ,'Endo',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
        Res2Var_wFluo_multiZ('14-05-19', '5mT', 'M2','DCLA_Endo',2 ,'Endo',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
        Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DCLA_Div',1 ,'Div',48,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)    
    
        Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M2','DCLA_Div',2 ,'Div',48,'28-08-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)    
    
        Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DCLA_Div',2 ,'Div',48,'18-07-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
       
    
    %% Dicty
    
     Res2Var_wFluo_multiZ('05-11-18', '5mT', 'M1','Dicty_WT',1 ,'Results',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
     Res2Var_wFluo_multiZ('05-12-18', '5mT', 'M1','Dicty_WT',1 ,'Results',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % endocytose
    
    Res2Var_wFluo_multiZ('05-11-18', '5mT', 'M1','Dicty_Endo',1 ,'Endo',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('05-12-18', '5mT', 'M1','Dicty_Endo',1 ,'Endo',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % div
    
    Res2Var_wFluo_multiZ('05-12-18', '5mT', 'M1','Dicty_Div',1 ,'Div',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    %% RPE1 
         Res2Var_wFluo_multiZ('04-04-19', '5mT', 'M1','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
         Res2Var_wFluo_multiZ('10-04-19', '5mT', 'M1','RPE1',[1 2] ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
         Res2Var_wFluo_multiZ('11-04-19', '5mT', 'M1','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)    
         Res2Var_wFluo_multiZ('11-04-19', '5mT', 'M2','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)    
         Res2Var_wFluo_multiZ('11-04-19', '5mT', 'M3','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
         Res2Var_wFluo_multiZ('17-04-19', '5mT', 'M1','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
         Res2Var_wFluo_multiZ('17-04-19', '5mT', 'M2','RPE1',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Inside/Outside Datas %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Manip etalonnage paula et vieux 0%
    
%     Res2Var_wFluo_multiZ('22-01-18', '5mT', 'M1','DC_Inside_0pc'    ,1:2,'Inside',31,'03-02-18_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     Res2Var_wFluo_multiZ('22-01-18', '5mT', 'M2','DC_Inside_0pc'    ,1:2,'Inside',31,'03-02-18_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     
%     Res2Var_wFluo_multiZ('23-01-18', '5mT', 'M1','DC_Inside_0pc'    ,1:2,'Inside',31,'03-02-18_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     Res2Var_wFluo_multiZ('23-01-18', '5mT', 'M2','DC_Inside_0pc'    ,1:2,'Inside',31,'03-02-18_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
    
%     Res2Var_wFluo_multiZ('25-03-19', '5mT', 'M1','DC_Inside_0pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)  
%     Res2Var_wFluo_multiZ('25-03-19', '5mT', 'M2','DC_Inside_5pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)  
%     Res2Var_wFluo_multiZ('25-03-19', '5mT', 'M3','DC_Inside_10pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)  
%     Res2Var_wFluo_multiZ('25-03-19', '5mT', 'M4','DC_Inside_12-5pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
    
        Res2Var_wFluo_multiZ('26-03-19', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)  
%     Res2Var_wFluo_multiZ('26-03-19', '5mT', 'M2','DC_Inside_15pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)  
%     Res2Var_wFluo_multiZ('26-03-19', '5mT', 'M3','DC_Inside_10pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)  
%     Res2Var_wFluo_multiZ('26-03-19', '5mT', 'M4','DC_Inside_12-5pc'    ,1:3,'Inside',45,'15-03-19_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)
%     
    
    %% 7.5 pc dans manip
    
    % Culture 1
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('28-05-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',147,'17-05-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('29-05-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 2
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('04-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('05-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 3
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('07-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('08-06-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',45,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 4
    Res2Var_wFluo_multiZ('02-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',46,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('02-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',46,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('03-07-18', '5mT','M1','DC_Inside_7-5pc'    ,1:3,'Inside',49,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('03-07-18', '5mT','M2','DC_Inside_7-5pc'    ,1:3,'Inside',49,'17-05-18_Depthograph'  ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 5
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('10-07-18', '5mT','M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'12-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 6
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('19-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('20-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'18-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 7
    Res2Var_wFluo_multiZ('30-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('31-07-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'31-07-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Arpin 1
    Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('17-09-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('18-09-18', '5mT', 'M2','DC_Inside_7-5pc'    ,1:3,'Inside',48,'28-08-18_Depthograph',AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 8
    Res2Var_wFluo_multiZ('29-10-18', '5mT', 'M1','DC_Inside_7-5pc'    ,1:3,'Inside',48,'30-10-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % Culture 10
    Res2Var_wFluo_multiZ('12-12-18', '5mT', 'M2','DC_Nase_Inside',1:3 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 12
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M1','DC_Inside_7-5pc',2 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DC_Inside_7-5pc',1 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('20-12-18', '5mT','M2','DC_Inside_7-5pc',2 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    
    % culture de 2019
    Res2Var_wFluo_multiZ('14-05-19', '5mT','M2','DC_Inside_DMSO',1 ,'Inside',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('16-05-19', '5mT','M1','DC_Inside_DMSO',2 ,'Inside',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    Res2Var_wFluo_multiZ('02-07-19', '5mT','M1','DC_Inside_DMSO',1 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    %% Manip 5 - 3 - 1 mT
    
    Res2Var_wFluo_multiZ('21-01-20', '5mT','M1','DC_Inside_5mT',1 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('21-01-20', '3mT','M1','DC_Inside_3mT',1 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('21-01-20', '1mT','M1','DC_Inside_1mT',1 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
        Res2Var_wFluo_multiZ('21-01-20', '5mT','M2','DC_Inside_5mT',1:2 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('21-01-20', '3mT','M2','DC_Inside_3mT',1:2 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    Res2Var_wFluo_multiZ('21-01-20', '1mT','M2','DC_Inside_1mT',1:2 ,'Inside',48,'26-06-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    %% autre manips
    
    Res2Var_wFluo_multiZ('21-12-18', '5mT','M1','DC_Nase_Inside',1 ,'Inside',48,'12-12-18_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 14
    Res2Var_wFluo_multiZ('02-04-19', '5mT', 'M1','DC_TrueOutside',1 ,'Outside',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
    % Culture 15
    Res2Var_wFluo_multiZ('13-05-19', '5mT', 'M2','DC_TrueOutside',1 ,'Results',48,'15-03-19_Depthograph' ,AUTO,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)
    
end


if 1 % V2D
    %% Nouvelles manips Inhib
    
    MakeClassification

    
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18','5mT',0.56,'M1','DCLA_CK666',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18','5mT',0.56,'M2','DCLA_DMSO',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18','5mT',0.56,'M1','DCLA_Blebbi',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18','5mT',0.56,'M2','DCLA_DMSO',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18','5mT',0.56,'M1','DCLA_Blebbi',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18','5mT',0.56,'M2','DCLA_DMSO',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18','5mT',0.56,'M1','DCLA_DMSO',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18','5mT',0.56,'M2','DCLA_CK666',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18','5mT',0.56,'M1','DCLA_DMSO',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18','5mT',0.56,'M2','DCLA_Blebbi',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18','5mT',0.56,'M1','DCLA_DMSO',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18','5mT',0.56,'M2','DCLA_CK666',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18','5mT',0.56,'M1','DCLA_DMSO',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18','5mT',0.56,'M1','DCLA_SMIFH2',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18','5mT',0.56,'M2','DCLA_SMIFH2',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18','5mT',0.56,'M1','DCLA_DMSO'  ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18','5mT',0.56,'M1','DCLA_DMSO'   ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18','5mT',0.56,'M2','DCLA_SMIFH2' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-07-18','5mT',0.56,'M1','DCLA_Blebbi' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-11-18','5mT',0.56,'M2','DCLA_NoDrug',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-12-18','5mT',0.56,'M1','DCLA_NoDrug',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-12-18','5mT',0.56,'M1','DCLA_SMIFH2',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-12-18','5mT',0.56,'M2','DCLA_NoDrug',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18','5mT',0.56,'M1','DCLA_SMIFH2',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18','5mT',0.56,'M2','DCLA_SMIFH2',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('01-04-19','5mT',0.56,'M1','DCLA_CalA1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('01-04-19','5mT',0.56,'M2','DCLA_CalA1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-04-19','5mT',0.56,'M2','DCLA_CalA1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('13-05-19','5mT',0.56,'M1','DCLA_CalA1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('13-05-19','5mT',0.56,'M3','DCLA_NoDrug',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19','5mT',0.56,'M1','DCLA_CalA1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19','5mT',0.56,'M2','DCLA_DMSO',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('16-05-19', '5mT',0.56, 'M1','DCLA_DMSO',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18','5mT',0.56,'M2','DCLA_LatA500',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18','5mT',0.56,'M1','DCLA_LatA500',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18','5mT',0.56,'M1','DCLA_LatA500',4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18','5mT',0.56,'M2','DCLA_LatA500',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    %     NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18','5mT',0.56,'M2','DCLA_ML141',4415,...
    %         RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18','5mT',0.56,'M1','DCLA_ML141' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18','5mT',0.56,'M2','DCLA_ML141' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18','5mT',0.56,'M2','DCLA_ML141'  ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
       % Culture 7
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18','5mT',0.56,'M1','DCLA_ML141'  ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18','5mT',0.56,'M1','DCLA_SMIFH2' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18','5mT',0.56,'M2','DCLA_Blebbi' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18','5mT',0.56,'M2','DCLA_DMSO'   ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('31-07-18','5mT',0.56,'M1','DCLA_LatA500' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('31-07-18','5mT',0.56,'M1','DCLA_CK666'   ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
       

    NotSaved = {NotSaved{:},Var2Data_wFluo('13-05-19','5mT',0.56,'M1','DCLA_CalA3',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19','5mT',0.56,'M1','DCLA_CalA3',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('01-04-19','5mT',0.56,'M1','DCLA_CalA10',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('01-04-19','5mT',0.56,'M2','DCLA_CalA10',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-04-19','5mT',0.56,'M2','DCLA_CalA10',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
        % Arpin 1
%     NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56,'M1','DCLA_ArpinWT' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
%     NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56,'M2','DCLA_ArpinKO' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
%     
%     NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M1','DCLA_ArpinKO' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
%     NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M2','DCLA_ArpinWT' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
%     NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M3','DCLA_ArpinKO' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
%     
    % Arpin 2
%     NotSaved = {NotSaved{:},Var2Data_wFluo('04-10-18', '5mT',0.56,'M1','DCLA_ArpinKO' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
%     NotSaved = {NotSaved{:},Var2Data_wFluo('04-10-18', '5mT',0.56,'M2','DCLA_ArpinWT' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
% 
%     
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-11-18','5mT',0.56,'M1','DCLA_Y27',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-12-18','5mT',0.56,'M2','DCLA_Nase',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-12-18','5mT',0.56,'M1','DCLA_Y27',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-12-18','5mT',0.56,'M2','DCLA_ArpC4WT',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-12-18','5mT',0.56,'M1','DCLA_ArpC4KO',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-12-18','5mT',0.56,'M1','DCLA_Y27',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18','5mT',0.56,'M1','DCLA_Y27',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('21-12-18','5mT',0.56,'M1','DCLA_LatNase',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % KO myoII
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-03-19','5mT',0.56,'M1','DCLA_MyoII-WT',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-03-19','5mT',0.56,'M2','DCLA_MyoII-KO',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-03-19','5mT',0.56,'M1','DCLA_MyoII-KO',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-03-19','5mT',0.56,'M2','DCLA_MyoII-WT',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    % culture 14
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-03-19','5mT',0.56,'M1','DCLA_LatTryp',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('21-03-19','5mT',0.56,'M1','DCLA_LatTryp',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    % SiRNA
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-04-19', '5mT',0.56, 'M1','DCLA_SiVim',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-04-19', '5mT',0.56, 'M1','DCLA_SiCtrl',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-04-19', '5mT',0.56, 'M1','DCLA_SiVim+LatA',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-04-19', '5mT',0.56, 'M1','DCLA_SiCtrl+LatA',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    
    %  SiRNA
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-05-19', '5mT',0.56, 'M1','DCLA_SiVim+LatA',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-05-19', '5mT',0.56, 'M1','DCLA_SiCtrl+LatA',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-05-19', '5mT',0.56, 'M2','DCLA_SiVim',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('15-05-19', '5mT',0.56, 'M2','DCLA_SiCtrl',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('16-05-19', '5mT',0.56, 'M1','DCLA_SiCtrl',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    % culture 17
        
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M1','DCLA_CalA1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M1','DCLA_CalA3',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M2','DCLA_DMSO',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M2','DCLA_CalA3',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
        NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M2','DCLA_Sep',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};    
        NotSaved = {NotSaved{:},Var2Data_wFluo('25-06-19', '5mT',0.56, 'M1','DCLA_Sep',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % culture 18        
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-19', '5mT',0.56, 'M1','DCLA_DMSO',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-19', '5mT',0.56, 'M2','DCLA_CalA1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-19', '5mT',0.56, 'M2','DCLA_CalA3',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};

    %% Endocytose et division
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56, 'M1','DCLA_Endo',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18', '5mT',0.56, 'M2','DCLA_Endo',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18', '5mT',0.56, 'M2','DCLA_Endo',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-012-18', '5mT',0.56, 'M2','DCLA_Endo',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('12-03-19', '5mT',0.56, 'M1','DCLA_Endo',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19', '5mT',0.56, 'M2','DCLA_Endo',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56, 'M1','DCLA_Div',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18', '5mT',0.56, 'M2','DCLA_Div',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56, 'M2','DCLA_Div',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    %% Dicty
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-11-18', '5mT',0.56, 'M1','Dicty_WT',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-11-18', '5mT',0.56, 'M1','Dicty_Endo',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-12-18', '5mT',0.56, 'M1','Dicty_WT',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-12-18', '5mT',0.56, 'M1','Dicty_Endo',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-12-18', '5mT',0.56, 'M1','Dicty_Div',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    %% RPE1
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-04-19', '5mT',0.56, 'M1','RPE1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-04-19', '5mT',0.56, 'M1','RPE1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-04-19', '5mT',0.56, 'M1','RPE1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-04-19', '5mT',0.56, 'M2','RPE1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('11-04-19', '5mT',0.56, 'M3','RPE1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-04-19', '5mT',0.56, 'M1','RPE1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-04-19', '5mT',0.56, 'M2','RPE1',4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    %% From Paula : 
    % removefromdata : enleve les données listé dans le dossier V2D
%j'enlève d abord les cellules mortes 100419 125 ?
RemoveFromData(MatfileFolder,{'11-04-19_M3_P1_C5_5mT','11-04-19_M3_P1_C6_5mT','11-04-19_M3_P1_C7_5mT','11-04-19_M3_P1_C8_5mT'},'RPE1')
%j'enlève a cause de problemes techniques et des cas bizzares
RemoveFromData(MatfileFolder,{'11-04-19_M1_P1_C3_5mT','11-04-19_M1_P1_C4_5mT','11-04-19_M1_P1_C9_5mT',...
    '11-04-19_M2_P1_C1_5mT', '11-04-19_M2_P1_C4-1_5mT', '11-04-19_M2_P1_C3_5mT','11-04-19_M2_P1_C5_5mT','10-04-19_M1_P1_C3_5mT',...
    '10-04-19_M1_P1_C3_5mT', '10-04-19_M1_P1_C5_5mT', '10-04-19_M1_P1_C6_5mT', '10-04-19_M1_P2_C4_5mT'},'RPE1')
%j'enlève d abord les cellules qui blebbent trop
RemoveFromData(MatfileFolder,{'17-04-19_M1_P1_C3_5mT','11-04-19_M3_P1_C1_5mT'},'RPE1')
    
    
    %% Inside Datas
    
    %% Manip etalonnage paula et vieux 0%
%     
%         NotSaved = {NotSaved{:},Var2Data_wFluo('22-01-18', '5mT',0.56,'M1','DC_Inside_0pc' ,4415,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};    
%         NotSaved = {NotSaved{:},Var2Data_wFluo('22-01-18', '5mT',0.56,'M2','DC_Inside_0pc' ,4415,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
%     
%         NotSaved = {NotSaved{:},Var2Data_wFluo('23-01-18', '5mT',0.56,'M1','DC_Inside_0pc' ,4415,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};    
%         NotSaved = {NotSaved{:},Var2Data_wFluo('23-01-18', '5mT',0.56,'M2','DC_Inside_0pc' ,4415,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
%     
%         NotSaved = {NotSaved{:},Var2Data_wFluo('25-03-19', '5mT',0.56,'M1','DC_Inside_0pc' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};  
%         NotSaved = {NotSaved{:},Var2Data_wFluo('25-03-19', '5mT',0.56,'M2','DC_Inside_5pc' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};  
%         NotSaved = {NotSaved{:},Var2Data_wFluo('25-03-19', '5mT',0.56,'M3','DC_Inside_10pc' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};  
%         NotSaved = {NotSaved{:},Var2Data_wFluo('25-03-19', '5mT',0.56,'M4','DC_Inside_12-5pc' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};  
%     
%         NotSaved = {NotSaved{:},Var2Data_wFluo('26-03-19', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};  
%         NotSaved = {NotSaved{:},Var2Data_wFluo('26-03-19', '5mT',0.56,'M2','DC_Inside_15pc' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};  
%         NotSaved = {NotSaved{:},Var2Data_wFluo('26-03-19', '5mT',0.56,'M3','DC_Inside_10pc' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};  
%         NotSaved = {NotSaved{:},Var2Data_wFluo('26-03-19', '5mT',0.56,'M4','DC_Inside_12-5pc' ,4438,...
%         RawdataFolder,MatfileFolder,DropboxDataFolder)};  
%     
    %% 7.5pc
    % Culture 1
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-05-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-05-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % Culture 2
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('04-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % Culture 3
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('07-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('08-06-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % Culture 4
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('03-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('03-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % Culture 5
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('10-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % Culture 6
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('19-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % Culture 7
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('31-07-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('31-07-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4415,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % Arpin 1
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    % Culture 8
    NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % Culture 12
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18', '5mT',0.56,'M1','DC_Inside_7-5pc' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('20-12-18', '5mT',0.56,'M2','DC_Inside_7-5pc' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    
    % culture 2019
    
       NotSaved = {NotSaved{:},Var2Data_wFluo('14-05-19', '5mT',0.56,'M2','DC_Inside_DMSO' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
       NotSaved = {NotSaved{:},Var2Data_wFluo('16-05-19', '5mT',0.56,'M1','DC_Inside_DMSO' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
       NotSaved = {NotSaved{:},Var2Data_wFluo('02-07-19', '5mT',0.56,'M1','DC_Inside_DMSO' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
      %% Manip 5 - 3 - 1 mT
    
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '5mT',1.1,'M1','DC_Inside_5mT' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '3mT',1.1,'M1','DC_Inside_3mT' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '1mT',1.1,'M1','DC_Inside_1mT' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
      
        NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '5mT',1.1,'M2','DC_Inside_5mT' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '3mT',1.1,'M2','DC_Inside_3mT' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
          NotSaved = {NotSaved{:},Var2Data_wFluo('21-01-20', '1mT',1.1,'M2','DC_Inside_1mT' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    

    
    %% autre manips
       
    NotSaved = {NotSaved{:},Var2Data_wFluo('21-12-18', '5mT',0.56,'M1','DC_Nase_Inside' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % culture 14
    NotSaved = {NotSaved{:},Var2Data_wFluo('02-04-19', '5mT',0.56,'M1','DC_TrueOutside' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    % culture 15
    NotSaved = {NotSaved{:},Var2Data_wFluo('13-05-19', '5mT',0.56,'M2','DC_TrueOutside' ,4438,...
        RawdataFolder,MatfileFolder,DropboxDataFolder)};
    
    NotSaved(find(strcmp(NotSaved,'empty'))) = []
    
    
    RemoveFromData(MatfileFolder,{'17-09-18_M2_P2_C4_5mT','10-07-18_M1_P2_C1_5mT','20-12-18_M1_P2_C4_5mT'},'Inside')
    RemoveFromData(MatfileFolder,{'02-04-19_M1_P1_C1-3_5mT'},'TrueOutside')
      
    %% data simu 
%     NicoDataV2D
    
    
    RemoveFromData(MatfileFolder,{'15-05-19_M1_P1_C5-1_5mT'},'Vim')
end


if ALL % D2P
    
    % Data principale

% Data2Peaks('NicoActiUp',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
% Data2Peaks('NicoActiDown',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
% Data2Peaks('NicoActiUpVar',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
% Data2Peaks('NicoActiDownVar',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

Data2Peaks('RPE1',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

Data2Peaks('DCLA_DMSO',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
Data2Peaks('DCLA_NoDrug',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
Data2Peaks('DCLA_LatA500',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

Data2Peaks('DCLA_Blebbi',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
Data2Peaks('DCLA_SMIFH2',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
Data2Peaks('DCLA_CK666',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
Data2Peaks('DC_Inside_DMSO',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

Data2Peaks('DC_TrueOutside',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

Data2Peaks('DCLA_CalA1',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

Data2Peaks('DCLA_CalA3',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

Data2Peaks('DCLA_CalA10',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
Data2Peaks('DCLA_ML141',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);



% Data2Peaks('DC_Inside_0pc',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
% Data2Peaks('DC_Inside_5pc',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
% Data2Peaks('DC_Inside_7-5pc',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
% Data2Peaks('DC_Inside_10pc',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
% Data2Peaks('DC_Inside_12-5pc',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
% Data2Peaks('DC_Inside_15pc',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

    
Data2Peaks('DC_Inside_DMSO',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

Data2Peaks('DC_Inside_1mT',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
Data2Peaks('DC_Inside_3mT',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
Data2Peaks('DC_Inside_5mT',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);

    
    Data2Peaks('DCLA_Y27',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_Nase',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_LatNase',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_LatTryp',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_MyoII-WT',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_MyoII-KO',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
%     Data2Peaks('DCLA_ArpinWT',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Data2Peaks('DCLA_ArpinKO',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_ArpC4WT',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_ArpC4KO',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
 
    Data2Peaks('DC_Nase_Inside',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_SiVim',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    Data2Peaks('DCLA_SiCtrl',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Data2Peaks('DCLA_SiVim+LatA',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    Data2Peaks('DCLA_SiCtrl+LatA',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
    
Data2Peaks('DCLA_Sep',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);


    Data2Peaks('Dicty_Div',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);


Data2Peaks('DCLA_WT',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);   
    
    
Data2Endo('DCLA_Endo',PLOTDP,MatfileFolder,FigureFolder,DropboxDataFolder);
end


if ALL % P2S
%     
%     Peaks2Stats('NicoActiDown',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('NicoActiUp',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('NicoActiDownVar',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('NicoActiUpVar',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);    
    

    Peaks2Stats('RPE1',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder); 
    
	Peaks2Stats('DCLA_DMSO',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder); 
    
    Peaks2Stats('DCLA_NoDrug',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_Blebbi',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_SMIFH2',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_CK666',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_CalA1',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_CalA3',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_CalA10',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_LatA500',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_ML141',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);    
    
    Peaks2Stats('DC_TrueOutside',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
%     Peaks2Stats('DC_Inside_0pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_5pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_7-5pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_10pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_12-5pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     Peaks2Stats('DC_Inside_15pc',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     
    Peaks2Stats('DC_Inside_1mT',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    Peaks2Stats('DC_Inside_3mT',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    Peaks2Stats('DC_Inside_5mT',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    
    Peaks2Stats('DC_Inside_DMSO',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
%    Peaks2Stats('DCLA_ArpinWT',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
%     Peaks2Stats('DCLA_ArpinKO',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
%     
   Peaks2Stats('DCLA_ArpC4WT',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_ArpC4KO',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    % pas encore trié
    
    
    
    
    
    Peaks2Stats('DCLA_Y27',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_MyoII-WT',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_MyoII-KO',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_LatNase',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_LatTryp',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_Nase',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DC_Nase_Inside',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_SiVim',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_SiCtrl',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_SiVim+LatA',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
    Peaks2Stats('DCLA_SiCtrl+LatA',PLOTPS,alignlvl,MatfileFolder,FigureFolder,DropboxDataFolder);
    
end


if 0 % Plots
    
    FigureMaker(1,MatfileFolder,FigureFolder)
    
    PEAKS    = 1;
    PROBA  = 1;
    DIST   = 0;
    TIME = 0;
    
    % comparaison individuelle
        Stats2Plots({'DCLA_DMSO','NicoActiDown'},...
        'NewDCLA-NICODOWN',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','NicoActiDownVar'},...
        'NewDCLA-NICODOWNVar',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','NicoActiUp'},...
        'NewDCLA-NICOUP',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','NicoActiUpVar'},...
        'NewDCLA-NICOUPVar',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
        Stats2Plots({'DC_Inside_3mT','DC_Inside_5mT'},...
        'Inside1-3-5mT',DIST,PEAKS,PROBA,TIME,alignlvl,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_DMSO','DCLA_Blebbi'},...
        'NewDCLA-Bleb',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','DCLA_CalA1','DCLA_CalA3','DCLA_CalA10'},...
        'NewDCLA-CalA_',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','DCLA_CK666'},...
        'NewDCLA-Arp',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','DCLA_SMIFH2'},...
        'NewDCLA-SM',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','DCLA_ML141'},...
        'NewDCLA-ML',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_ArpC4WT','DCLA_ArpC4KO'},...
        'NewDCLA-ArpC4',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_ArpinWT','DCLA_ArpinKO'},...
        'NewDCLA-Arpin',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_MyoII-WT','DCLA_MyoII-KO'},...
        'NewDCLA-MyoII',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
      
    Stats2Plots({'DCLA_DMSO','DCLA_LatA500'},...
        'NewDCLA-LatA',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    
    Stats2Plots({'DCLA_DMSO','Dicty_WT'},...
        'DC-Dicty',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_DMSO','RPE1'},...
        'DC-RPE1',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    % pooled figs
    
    Stats2Plots({'DCLA_DMSO','DC_Inside_7-5pc','DC_TrueOutside','DCLA_LatA500'},...
        'NewDCLA-LatA,inout',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','DC_Inside_7-5pc','DC_TrueOutside'},...
        'NewDCLA-inout',DIST,PEAKS,PROBA,0.5,MatfileFolder,FigureFolder)
        Stats2Plots({'DC_Inside_0pc','DC_Inside_5pc','DC_Inside_7-5pc','DC_Inside_10pc','DC_Inside_12-5pc','DC_Inside_15pc'},...
        'NewDCLA-insides',DIST,PEAKS,PROBA,0.5,MatfileFolder,FigureFolder)
    
        Stats2Plots({'DCLA_DMSO','DC_Inside_7-5pc','DCLA_LatA500'},...
        'NewDCLA-LatA,in',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_DMSO','DCLA_LatA500','DCLA_Blebbi','DCLA_CK666',...
        'DCLA_SMIFH2'},...
        'NewDCLA-Basic',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
        
    Stats2Plots({'DCLA_DMSO','DCLA_Blebbi','DCLA_LatA500'},...
        'NewDCLA-BBLat',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_DMSO','DCLA_Blebbi','DCLA_CalA1','DCLA_CK666',...
        'DCLA_SMIFH2'},...
        'NewDCLA-BasicnoLA',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)   

     FigureMaker(1,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_DMSO','DCLA_NoDrug'},...
        'NewDCLA-Controls',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)    
    Stats2Plots({'DCLA_NoDrug','DCLA_Y27'},...
        'NewDCLA-Y27',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)    
    Stats2Plots({'DCLA_MyoII-WT','DCLA_MyoII-KO'},...
        'NewDCLA-Myo',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)    
    Stats2Plots({'DCLA_ArpinWT','DCLA_ArpinKO'},...
        'NewDCLA-Arpin',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)    
    Stats2Plots({'DCLA_ArpC4WT','DCLA_ArpC4KO'},...
        'NewDCLA-ArpC4',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)    

    
    
    
    
    
    Stats2Plots({'DCLA_DMSO','DCLA_CK666','DCLA_ArpinWT','DCLA_ArpinKO'},...
        'NewDCLA-Arp23',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    Stats2Plots({'NicoActiDown','NicoActiDownVar','NicoActiUp','NicoActiUpVar'},...
        'NewDCLA-NICODOWN',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
            Stats2Plots({'DCLA_DMSO','NicoActiDown','NicoActiDownVar'},...
        'NewDCLA-NICODOWN+',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)

    Stats2Plots({'DCLA_DMSO','NicoActiUp','NicoActiUpVar'},...
        'NewDCLA-NICOUP+',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)

    
    
    Stats2Plots({'DCLA_DMSO','DCLA_CalA1','DCLA_CalA3','DCLA_CalA10'},...
        'NewDCLA-CalA',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
        Stats2Plots({'DCLA_DMSO','DCLA_CalA3'},...
        'NewDCLA-CalA3',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
        Stats2Plots({'DCLA_DMSO','DCLA_CalA1'},...
        'NewDCLA-CalA1',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_DMSO','DCLA_Blebbi','DCLA_CalA1','DCLA_SMIFH2','DCLA_CK666'},...
        'NewDCLA-Fig3',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    
    Stats2Plots({'DCLA_DMSO2','DCLA_SiCtrl','DCLA_SiVim','DCLA_SiCtrl+LatA','DCLA_SiVim+LatA'},...
        'NewDCLA-KDvim',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_DMSO','DCLA_Blebbi','DCLA_Y27','DCLA_CalA1'},...
        'NewDCLA-Myo',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    
    Stats2Plots({'DCLA_DMSO','DCLA_CalA1','DCLA_CalA3'},...
        'NewDCLA-DMSOCalA',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)

    
    
    Stats2Plots({'DCLA_DMSO','DCLA_Blebbi','DCLA_CalA1'},...
        'NewDCLA-BBCal1',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
        Stats2Plots({'DCLA_DMSO','DC_Inside_7-5pc','DC_TrueOutside'},...
        'NewDCLA-inout',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    Stats2Plots({'DCLA_DMSO','DCLA_LatA500','DC_Inside_7-5pc'},...
        'NewDCLA-LatAin',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    Stats2Plots({'DCLA_DMSO','DCLA_LatA500','DC_Inside_7-5pc','DC_TrueOutside'},...
        'NewDCLA-LatAinout',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    

        
 
    
        Stats2Plots({'DCLA_DMSO'},...
        'NewDCLA-DMSO',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    
    
        Stats2Plots({'DCLA_DMSO','DCLA_LatA500','DCLA_SMIFH2','DCLA_CK666'},...
        'NewDCLA-F',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
        Stats2Plots({'DCLA_DMSO','DCLA_LatA500','DCLA_Blebbi','DCLA_CalA1'},...
        'NewDCLA-Fig',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
            Stats2Plots({'DCLA_DMSO','DCLA_Blebbi','DCLA_CalA1','DCLA_CK666',...
        'DCLA_SMIFH2'},...
        'NewDCLA-BasicnoLA',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)
    
    
    
    Stats2Plots({'DCLA_DMSO','DCLA_Blebbi','DCLA_LatA500'},...
        'NewDCLA',DIST,PEAKS,PROBA,alignlvl,MatfileFolder,FigureFolder)    

end
