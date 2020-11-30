close all;
clear all;
clc;

%% Paths & options
%% Paths
RawdataFolder = 'D:\Data\Raw';
MatfileFolder = 'D:\Data\MatFile';
DropboxDataFolder = '';
FigureFolder = 'D:\Data\Figures';

set(0, 'defaultfigureposition', get(0, 'Screensize'));

set(0,'DefaultFigureWindowStyle','docked')

% Couleurs par drogues
Cct = [0 109 219]./255; % control

Cl  = [150 0 255]./255; % Latranculin

Cb  = [219 209 0]./255; % blebbi
Cy  = [255 255 109 ]./255; % y27
Ccl = [146 0 0]./255; % calyculinA 10nM
Ccl10 = [146 0 0]./255./2; % calyculinA 1nM

Cs  = [36 255 36]./255; % smifh2

Cck = [182 109 255]./255;

Cin = [0 73 73]./255; % inside
Co  = [255 182 119]./255; % outside

AUTO = 0;

% v2d
NotSaved = {};

% d2p
PLOTDFCurve = 1;
PLOTDFComp = 1;

% d2m
PLOTD2M = 0;
VERBOSED2M = 1;

% diff tmp comp
Cdc1  = [255 127 63]./255; % diff tmp comp
Cdc2  = [127 63 255]./255; % diff tmp comp
Cdc3  = [63 255 127]./255; % diff tmp comp

return % stop execution here

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Res2Var Comp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% BeadsIn

Res2Var_wFluo_multiZ_Comp('Inside','28-08-18','R80','M1','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Inside','28-08-18','R80','M2','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Inside','17-09-18','R80','M1','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Inside','17-09-18','R80','M2','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Inside','18-09-18','R80','M1','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Inside','18-09-18','R80','M2','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Inside','24-09-18','R80','M1','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Inside','24-09-18','R80','M2','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Inside','25-09-18','R80','M1','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Inside','25-09-18','R80','M2','DC-Comp_Inside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)



Res2Var_wFluo_multiZ_Comp('Inside','12-12-18','R80','M1','DC-Comp_Inside',1,...
    '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Inside','12-12-18','R80','M2','DC-Comp_NaseInside',1:3,...
    '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)


Res2Var_wFluo_multiZ_Comp('Inside','18-12-18','R80','M2','DC-Comp_NaseInside',2,...
    '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

%% Beads out
Res2Var_wFluo_multiZ_Comp('Outside','25-09-18','R80','M1','DC-Comp_Outside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Outside','25-09-18','R80','M2','DC-Comp_Outside',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Outside','30-10-18','R80','M1','DC-Comp_Outside',1:3,...
    '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Outside','30-10-18','R80','M2','DC-Comp_Outside',1:3,...
    '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%
Res2Var_wFluo_multiZ_Comp('Outside','21-01-19','R80','M1','DC-Comp_Bsa',2,...
    '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Outside','23-01-19','R80','M1','DC-Comp_SerumJ1',1,...
    '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Outside','23-01-19','R80','M1','DC-Comp_Naked',2,...
    '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Outside','28-01-19','R80','M1','DC-Comp_SerumJ6',1,...
    '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

%% InOut DC

% CultureRampe1   18s (42 + 42 img) + 2s (95img) = 179
Res2Var_wFluo_multiZ_Comp('Results','28-08-18','R80','M1','DC-Comp_Ctrl',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','28-08-18','R80','M2','DC-Comp_Ctrl',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

% Arpin 1 18s (42 + 42 img) + 2s (95img) = 179
Res2Var_wFluo_multiZ_Comp('Results','17-09-18','R80','M1','DC-Comp_ArpinWT',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','17-09-18','R80','M2','DC-Comp_ArpinKO',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

Res2Var_wFluo_multiZ_Comp('Results','18-09-18','R80','M1','DC-Comp_ArpinKO',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','18-09-18','R80','M2','DC-Comp_ArpinWT',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','18-09-18','R80','M3','DC-Comp_ArpinKO',1:3,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

% CultureRampe Inhib 1 18s (42 + 42 img) + 2s (95img) = 179
Res2Var_wFluo_multiZ_Comp('Results','24-09-18','R80','M1','DC-Comp_CK666',1,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','24-09-18','R80','M1','DC-Comp_Ctrl',2,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','24-09-18','R80','M2','DC-Comp_SMIFH2',1,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)


Res2Var_wFluo_multiZ_Comp('Results','25-09-18','R80','M1','DC-Comp_LatA500',1,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','25-09-18','R80','M1','DC-Comp_SMIFH2',2,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','25-09-18','R80','M2','DC-Comp_Ctrl',1,...
    '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','25-09-18','R80','M2','DC-Comp_CK666',2,...
%         '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

% Arpin 2 18s (42 + 42 img) + 2s (95img) = 179
%     Res2Var_wFluo_multiZ_Comp('Results','05-10-18','R80','M1','DC-Comp_ArpinKO',1:3,...
%         '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','05-10-18','R80','M2','DC-Comp_ArpinWT',1:3,...
%         '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%
% CultureRampe Inhib 2 18s (42 + 42 img) + 2s (95img) = 179
%     Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M1','DC-Comp_CalA25',1,...
%         '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M1','DC-Comp_LatA2',2,...
%         '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M2','DC-Comp_LatA2',2,...
%         '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)


%     Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M1','DC-Comp_Blebbi',2,...
%         '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M2','DC-Comp_CalA25',1,...
%         '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M2','DC-Comp_Ctrl',2,...
    '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)


% CultureRampe Inhib 3 18s (42 + 42 img) + 2s (95img) = 179
Res2Var_wFluo_multiZ_Comp('Results','12-12-18','R80','M1','DC-Comp_Ctrl',1,...
    '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','12-12-18','R80','M2','DC-Comp_Nase',1:3,...
%         '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)


% CultureRampe ArpC4 1 18s (42 + 42 img) + 2s (95img) = 179
%     Res2Var_wFluo_multiZ_Comp('Results','17-12-18','R80','M1','DC-Comp_AprC4KO',1,...
%         '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','17-12-18','R80','M2','DC-Comp_AprC4WT',1,...
%         '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%

% CultureRampe Inhib 4 18s (42 + 42 img) + 2s (95img) = 179
%     Res2Var_wFluo_multiZ_Comp('Results','18-12-18','R80','M1','DC-Comp_Y27',1,...
%         '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','18-12-18','R80','M2','DC-Comp_Nase',2,...
%         '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

% CultureRampe Inhib 5 18s (42 + 42 img) + 2s (95img) = 179
%     Res2Var_wFluo_multiZ_Comp('Results','21-12-18','R80','M1','DC-Comp_LatANase',1,...
%         '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%
%     % CultureRampe Inhib 5 18s (42 + 42 img) + 2s (95img) = 179
%     Res2Var_wFluo_multiZ_Comp('Results','30-04-19','R80','M1','DC-Comp_SiCtrl+LatA',2,...
%         '26-06-19_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
%     Res2Var_wFluo_multiZ_Comp('Results','30-04-19','R80','M1','DC-Comp_SiVim+LatA',1,...
%         '26-06-19_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)

%% Manip avec differents temps de compression
Res2Var_wFluo_multiZ_Comp_R248('Results','14-10-19','RV248','M1','DC-Comp_R248_Ctrl',1:2,...
    '15-10-19_Depthograph',21,24,98,193,383,...
    AUTO,RawdataFolder,MatfileFolder)


Res2Var_wFluo_multiZ_Comp_R248('Results','14-10-19','RV248','M2','DC-Comp_R248_LatA',1:2,...
    '15-10-19_Depthograph',21,24,98,193,383,...
    AUTO,RawdataFolder,MatfileFolder)


Res2Var_wFluo_multiZ_Comp_R248('Results','15-10-19','RV248','M1','DC-Comp_R248_Blebbi',1:2,...
    '15-10-19_Depthograph',21,24,98,193,383,...
    AUTO,RawdataFolder,MatfileFolder)


Res2Var_wFluo_multiZ_Comp_R248('Results','15-10-19','RV248','M2','DC-Comp_R248_Ctrl',1:2,...
    '15-10-19_Depthograph',21,24,98,193,383,...
    AUTO,RawdataFolder,MatfileFolder)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Var2Data Comp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% InOut DC

% manip ctrl test
Var2Data_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,'1s',4415,MatfileFolder,FigureFolder);
Var2Data_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,'1s',4415,MatfileFolder,FigureFolder);

%     % Arpin 1
%     Var2Data_Comp('17-09-18','R80',0.56,'M1','DC-Comp_ArpinWT',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('17-09-18','R80',0.56,'M2','DC-Comp_ArpinKO',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
%     Var2Data_Comp('18-09-18','R80',0.56,'M1','DC-Comp_ArpinKO',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('18-09-18','R80',0.56,'M2','DC-Comp_ArpinWT',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('18-09-18','R80',0.56,'M3','DC-Comp_ArpinKO',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
% Inhib 1

%     Var2Data_Comp('24-09-18','R80',0.56,'M1','DC-Comp_CK666',95,179,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('24-09-18','R80',0.56,'M2','DC-Comp_SMIFH2',95,179,'1s',4438,MatfileFolder,FigureFolder);

Var2Data_Comp('25-09-18','R80',0.56,'M1','DC-Comp_LatA500',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('25-09-18','R80',0.56,'M1','DC-Comp_SMIFH2',95,179,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('25-09-18','R80',0.56,'M2','DC-Comp_CK666',95,179,'1s',4438,MatfileFolder,FigureFolder);

% Arpin 2;
%     Var2Data_Comp('05-10-18','R80',0.56,'M1','DC-Comp_ArpinKO',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('05-10-18','R80',0.56,'M2','DC-Comp_ArpinWT',95,179,'1s',4438,MatfileFolder,FigureFolder);


% Inhib 2
%
%     Var2Data_Comp('29-10-18','R80',0.56,'M1','DC-Comp_CalA25',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('29-10-18','R80',0.56,'M1','DC-Comp_LatA2',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('29-10-18','R80',0.56,'M2','DC-Comp_LatA2',95,179,'1s',4438,MatfileFolder,FigureFolder);


%     Var2Data_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Blebbi',95,179,'1s',4438)
%     Var2Data_Comp('30-10-18','R80',0.56,'M1','DC-Comp_CalA25',95,179,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,'1s',4438,MatfileFolder,FigureFolder);

% inhib 3


Var2Data_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('12-12-18','R80',0.56,'M2','DC-Comp_Nase',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
%
%     % ArpC4
%
%     Var2Data_Comp('17-12-18','R80',0.56,'M1','DC-Comp_AprC4KO',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('17-12-18','R80',0.56,'M2','DC-Comp_AprC4WT',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
%     % inhib 4
%
%     Var2Data_Comp('18-12-18','R80',0.56,'M1','DC-Comp_Y27',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('18-12-18','R80',0.56,'M2','DC-Comp_Nase',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
%     % Inhib 5
%
%     Var2Data_Comp('21-12-18','R80',0.56,'M1','DC-Comp_LatANase',95,179,'1s',4438,MatfileFolder,FigureFolder);

%% 3T3 Joseph

Var2Data_Comp('04-08-20','R40_3T3aSFLnodrugs',1.15,'M1','3T3aSFL_BSA_nodrugs',95,143,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('05-08-20','R40_3T3aSFLnodrugs',1.15,'M2','3T3aSFL_BSA_nodrugs',95,143,'1s',4504,MatfileFolder,FigureFolder);
Var2Data_Comp('07-08-20','R40_3T3aSFLnodrugs',1.15,'M2','3T3aSFL_BSA_nodrugs',95,143,'1s',4504,MatfileFolder,FigureFolder);

%% billes internes 
Var2Data_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Inside',95,179,'1s',4415,MatfileFolder,FigureFolder);
Var2Data_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Inside',95,179,'1s',4415,MatfileFolder,FigureFolder);


Var2Data_Comp('17-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp('17-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,'1s',4438,MatfileFolder,FigureFolder);


Var2Data_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp('24-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,'1s',4438,MatfileFolder,FigureFolder);


Var2Data_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,'1s',4438,MatfileFolder,FigureFolder);


Var2Data_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Inside',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('12-12-18','R80',0.56,'M2','DC-Comp_NaseInside',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
%     Var2Data_Comp('18-12-18','R80',0.56,'M2','DC-Comp_NaseInside',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
%
%     Var2Data_Comp('30-04-19','R80',0.56,'M1','DC-Comp_SiCtrl+LatA',95,179,'1s',4438,MatfileFolder,FigureFolder);
%     Var2Data_Comp('30-04-19','R80',0.56,'M1','DC-Comp_SiVim+LatA',95,179,'1s',4438) ;

%% billes externes
%
Var2Data_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Outside',95,179,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Outside',95,179,'1s',4438,MatfileFolder,FigureFolder);

Var2Data_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Outside',95,179,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Outside',95,179,'1s',4438,MatfileFolder,FigureFolder);

%     Var2Data_Comp('21-01-19','R80',0.56,'M1','DC-Comp_Bsa',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
dist = []
dist = [dist Var2Data_Comp('23-01-19','R80',0.56,'M1','DC-Comp_SerumJ1',95,179,'1s',4438)];
Var2Data_Comp('23-01-19','R80',0.56,'M1','DC-Comp_Naked',95,179,'1s',4438,MatfileFolder,FigureFolder);
%
%
dist = [dist Var2Data_Comp('28-01-19','R80',0.56,'M1','DC-Comp_SerumJ6',95,179,'1s',4438)];
%

%% compression de tps diff

Var2Data_Comp_JV_R248('14-10-19','RV248',0.56,'M1','DC-Comp_R248_Ctrl',21,24,98,193,383,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp_JV_R248('14-10-19','RV248',0.56,'M2','DC-Comp_R248_LatA',21,24,98,193,383,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp_JV_R248('15-10-19','RV248',0.56,'M1','DC-Comp_R248_Blebbi',21,24,98,193,383,'1s',4438,MatfileFolder,FigureFolder);
Var2Data_Comp_JV_R248('15-10-19','RV248',0.56,'M2','DC-Comp_R248_Ctrl',21,24,98,193,383,'1s',4438,MatfileFolder,FigureFolder);

%% remove


RemoveFromData(MatfileFolder,{'24-09-18_M1_P1_C1_R80','24-09-18_M1_P1_C3_R80','24-09-18_M1_P1_C5_R80',...
    '24-09-18_M1_P1_C6_R80','24-09-18_M2_P1_C1-2_R80','24-09-18_M2_P1_C2_R80',...
    '24-09-18_M2_P1_C4_R80','25-09-18_M2_P1_C2_R80','25-09-18_M2_P1_C3-1_R80',...
    '25-09-18_M2_P2_C2_R80','25-09-18_M2_P2_C3_R80','25-09-18_M1_P1_C1_R80',...
    '25-09-18_M1_P1_C3-1_R80','25-09-18_M1_P1_C3-2_R80','25-09-18_M1_P1_C5_R80'...
    },'Comp')

%
%     % tri des insides
%     RemoveFromData({'12-12-18_M2_P1_C9_R80in','12-12-18_M2_P1_C10_R80in',...
%         '24-09-18_M2_P1_C1_R80in',...
%         '25-09-18_M1_P2_C4_R80in','28-08-18_M2_P1_C2_R80in',...
%         '25-09-18_M1_P1_C1_R80in','28-08-18_M1_P1_C1_R80in',... % Ctrl
%         '12-12-18_M2_P2_C8_R80in','18-12-18_M2_P2_C5_R80in','18-12-18_M2_P2_C12_R80in'... % Nase
%         },'Inside')
%
%
%     CutCurves({'17-09-18_M2_P1_C6_R80in#ap#150','18-09-18_M2_P2_C3_R80in#av#150',...
%         '24-09-18_M1_P2_C4_R80in#ex#80-115',...
%         '25-09-18_M2_P1_C3_R80in#ap#240','28-08-18_M2_P1_C3_R80in#ap#200',...
%         '25-09-18_M1_P1_C1_R80in#av#75','28-08-18_M1_P1_C1_R80in#av#75',...
%         });

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data2Meca %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% D2M pour these et papier


Data2Meca('28-08-18','M1','DC-Comp_Ctrl','R80',4415,PLOTD2M,MatfileFolder,FigureFolder);
Data2Meca('28-08-18','M2','DC-Comp_Ctrl','R80',4415,PLOTD2M,MatfileFolder,FigureFolder);
Data2Meca('24-09-18','M1','DC-Comp_Ctrl','R80',4438,PLOTD2M,MatfileFolder,FigureFolder);
Data2Meca('25-09-18','M2','DC-Comp_Ctrl','R80',4438,PLOTD2M,MatfileFolder,FigureFolder);
Data2Meca('30-10-18','M2','DC-Comp_Ctrl','R80',4438,PLOTD2M,MatfileFolder,FigureFolder);
Data2Meca('12-12-18','M1','DC-Comp_Ctrl','R80',4438,PLOTD2M,MatfileFolder,FigureFolder);


Meca2Plot('Comp_Ctrl',Cct,MatfileFolder,FigureFolder)

Data2Meca('24-09-18','M1','DC-Comp_CK666',4438,PLOTD2M,MatfileFolder,FigureFolder);
Data2Meca('25-09-18','M2','DC-Comp_CK666',4438,PLOTD2M,MatfileFolder,FigureFolder);

Data2Meca('24-09-18','M2','DC-Comp_SMIFH2',4438,PLOTD2M,MatfileFolder,FigureFolder) ;
Data2Meca('25-09-18','M1','DC-Comp_SMIFH2',4438,PLOTD2M,MatfileFolder,FigureFolder);


Data2Meca('30-10-18','M1','DC-Comp_Blebbi',4438,PLOTD2M,MatfileFolder,FigureFolder);


[a,b] =  Data2Meca('25-09-18','M1','DC-Comp_LatA500','r80',4438,PLOTD2M,MatfileFolder,FigureFolder);



Meca2Plot('Comp_CK666',Cck,MatfileFolder,FigureFolder)
Meca2Plot('Comp_SMIFH2',Cs,MatfileFolder,FigureFolder)
Meca2Plot('Comp_Blebbi',Cb,MatfileFolder,FigureFolder)
Meca2Plot('Comp_LatA500',Cl,MatfileFolder,FigureFolder)

%% D2M pour comp differente longueur

% DC
Data2Meca_R248('14-10-19','M1','DC-Comp_R248_Ctrl',4438,0)
%     Data2Meca_R248('14-10-19','M2','DC-Comp_R248_LatA',4438,0)
%     Data2Meca_R248('15-10-19','M1','DC-Comp_R248_Blebbi',4438,0)
Data2Meca_R248('15-10-19','M2','DC-Comp_R248_Ctrl',4438,0)





Meca2Plot_R248('Comp_R248_Ctrl',Cct,1)
Meca2Plot_R248('Comp_R248_LatA',Cl,1)
Meca2Plot_R248('Comp_R248_Blebbi',Cb,1)

Data2MecaLoop('14-10-19','M1','DC-Comp_R248_Ctrl-1s','RV248',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop('14-10-19','M1','DC-Comp_R248_Ctrl-2s','RV248',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop('14-10-19','M1','DC-Comp_R248_Ctrl-4s','RV248',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop('15-10-19','M2','DC-Comp_R248_Ctrl-1s','RV248',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop('15-10-19','M2','DC-Comp_R248_Ctrl-2s','RV248',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop('15-10-19','M2','DC-Comp_R248_Ctrl-4s','RV248',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);


Meca2Plot_Multi({'DC-Comp_R248_Ctrl-1s','DC-Comp_R248_Ctrl-2s','DC-Comp_R248_Ctrl-4s'},{Cdc1,Cdc2,Cdc3},{'DC-1s','DC-2s','DC-4s'},MatfileFolder,FigureFolder)

%% D2MLoop_BigTable
PLOTD2M = 0;
VERBOSED2M = 0;
Data2MecaLoop_BigTable('28-08-18','M1','DC-Comp_Ctrl','R80',4415,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);
Data2MecaLoop_BigTable('28-08-18','M2','DC-Comp_Ctrl','R80',4415,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);


Data2MecaLoop_BigTable('24-09-18','M1','DC-Comp_Ctrl','R80',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('25-09-18','M2','DC-Comp_Ctrl','R80',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('30-10-18','M2','DC-Comp_Ctrl','R80',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);

Data2MecaLoop_BigTable('12-12-18','M1','DC-Comp_Ctrl','R80',4438,MatfileFolder,FigureFolder,PLOTD2M,VERBOSED2M);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Meca2Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

fitparams = 'Strain100-R20.9';

% comparaison DC/Dicty/3T3
ConditionsDCDicty = {'ExpType','DC-Comp_Ctrl','TpsComp','1s','FitParams',fitparams;...
    'ExpType','DictyDB_WT','TpsComp','1s','FitParams',fitparams;...
    'ExpType','3T3aSFL_BSA_nodrugs','TpsComp','1s','FitParams',fitparams;...
    }

Meca2Plot_BigTable(MatfileFolder,FigureFolder,ConditionsDCDicty, ...
{Cct,Cb,Cck},{'o','o','o'},...
{'DC','Dicty','3T3'},'DC-Dicty-3T3',1,...
'Exclude',{'ExpDay','04-11-20'})

%% V2D D2P P2S sur cst part
% data from compression

NotSaved = {NotSaved{:},Var2Data_wFluo('28-08-18','R80cst',0.6,'M1','DC-Comp_Ctrl',4415,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('28-08-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4415,RawdataFolder,MatfileFolder)};

NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18','R80cst',0.6,'M1','DC-Comp_ArpinWT',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18','R80cst',0.6,'M2','DC-Comp_ArpinKO',4438,RawdataFolder,MatfileFolder)};

NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M1','DC-Comp_ArpinKO',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M2','DC-Comp_ArpinWT',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M3','DC-Comp_ArpinKO',4438,RawdataFolder,MatfileFolder)};

NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M2','DC-Comp_SMIFH2',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M1','DC-Comp_CK666',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M1','DC-Comp_Ctrl',4438,RawdataFolder,MatfileFolder)};

NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M1','DC-Comp_LatA500',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M1','DC-Comp_SMIFH2',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M2','DC-Comp_CK666',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4438,RawdataFolder,MatfileFolder)};

NotSaved = {NotSaved{:},Var2Data_wFluo('05-10-18','R80cst',0.6,'M1','DC-Comp_ArpinKO',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('05-10-18','R80cst',0.6,'M2','DC-Comp_ArpinWT',4438,RawdataFolder,MatfileFolder)};

% NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M1','DC-Comp_CalA25',4438,RawdataFolder,MatfileFolder)};
% NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M1','DC-Comp_LatA2',4438,RawdataFolder,MatfileFolder)};
% NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M2','DC-Comp_LatA2',4438,RawdataFolder,MatfileFolder)};

NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M1','DC-Comp_Blebbi',4438,RawdataFolder,MatfileFolder)};
% NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M2','DC-Comp_CalA25',4438,RawdataFolder,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4438,RawdataFolder,MatfileFolder)};


NotSaved = {NotSaved{:},Var2Data_wFluo('12-12-18','R80cst',0.56,'M1','DC-Comp_Ctrl',4438,RawdataFolder,MatfileFolder)};

% 3T3 joseph

NotSaved = {NotSaved{:},Var2Data_wFluo('04-08-20','R40_3T3aSFLnodrugscst',1.15,'M1','3T3aSFL_BSA_nodrugs',4504,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('05-08-20','R40_3T3aSFLnodrugscst',1.15,'M2','3T3aSFL_BSA_nodrugs',4504,MatfileFolder)};
NotSaved = {NotSaved{:},Var2Data_wFluo('07-08-20','R40_3T3aSFLnodrugscst',1.15,'M2','3T3aSFL_BSA_nodrugs',4504,MatfileFolder)};

%%%%%%%%% billes internes %%%%%%%%%
Var2Data_wFluo('28-08-18','R80cst',0.56,'M1','DC-Comp_Inside',4415,MatfileFolder)
Var2Data_wFluo('28-08-18','R80cst',0.56,'M2','DC-Comp_Inside',4415,MatfileFolder)


Var2Data_wFluo('17-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,MatfileFolder)
Var2Data_wFluo('17-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,MatfileFolder)


Var2Data_wFluo('24-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,MatfileFolder)
Var2Data_wFluo('24-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,MatfileFolder)


Var2Data_wFluo('25-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,MatfileFolder)
Var2Data_wFluo('25-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,MatfileFolder)


Var2Data_wFluo('12-12-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,MatfileFolder)
Var2Data_wFluo('12-12-18','R80cst',0.56,'M2','DC-Comp_NaseInside',4438,MatfileFolder)

Var2Data_wFluo('18-12-18','R80cst',0.56,'M2','DC-Comp_NaseInside',4438,MatfileFolder)

% billes externes

Var2Data_wFluo('25-09-18','R80cst',0.56,'M1','DC-Comp_Outside',4438,MatfileFolder)
Var2Data_wFluo('25-09-18','R80cst',0.56,'M2','DC-Comp_Outside',4438,MatfileFolder)

Var2Data_wFluo('30-10-18','R80cst',0.56,'M1','DC-Comp_Outside',4438,MatfileFolder)
Var2Data_wFluo('30-10-18','R80cst',0.56,'M2','DC-Comp_Outside',4438,MatfileFolder)

Var2Data_wFluo('21-01-19','R80cst',0.56,'M1','DC-Comp_Bsa',4438,MatfileFolder)

Var2Data_wFluo('23-01-19','R80cst',0.56,'M1','DC-Comp_Naked',4438,MatfileFolder)
Var2Data_wFluo('23-01-19','R80cst',0.56,'M1','DC-Comp_SerumJ1',4438,MatfileFolder)


Var2Data_wFluo('28-01-19','R80cst',0.56,'M1','DC-Comp_SerumJ6',4438,MatfileFolder)


% data from comp

Data2Peaks('DC-Comp_Ctrl',0,MatfileFolder,FigureFolder);

Data2Peaks({'M1','M2'},'DC-Comp_Inside',0,0.9,MatfileFolder,FigureFolder);
Data2Peaks({'M1','M2'},'DC-Comp_Outside',1,0.9,MatfileFolder,FigureFolder);

Data2Peaks({'M1','M2'},'DC-Comp_Naked',0,0.9,MatfileFolder,FigureFolder);

Data2Peaks({'M1','M2'},'DC-Comp_SerumJ6',0,0.9,MatfileFolder,FigureFolder);



Data2Peaks('3T3aSFL_BSA_nodrugs',1,MatfileFolder,FigureFolder);

% data from comp
PLOTPS = 0;
Peaks2Stats('DC-Comp_Ctrl',PLOTPS,0.9,MatfileFolder,FigureFolder);


Peaks2Stats('3T3aSFL_BSA_nodrugs',PLOTPS,0.9,MatfileFolder,FigureFolder);



Peaks2Stats('DC-Comp_Inside',PLOTPS,MatfileFolder,FigureFolder);
Peaks2Stats('DC-Comp_Outside',PLOTPS,MatfileFolder,FigureFolder);

Peaks2Stats('DC-Comp_Naked');
Peaks2Stats('DC-Comp_SerumJ6');

FigureMaker(0,MatfileFolder,FigureFolder)


Stats2PlotsDictys({'DC-Comp_Ctrl','DictyDB_WT','3T3aSFL_BSA_nodrugs'},{Cct,Cb,Cck},'DC-Dicty-3T3',1,1,1,1,0.9,MatfileFolder,FigureFolder)
