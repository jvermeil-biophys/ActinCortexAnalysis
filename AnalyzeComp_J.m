close all;
clear all;
clc;

%% Paths & options
GeneralFolder = 'D:\Matlab Analysis\';
RawdataFolder = [GeneralFolder 'Data_Joseph\Raw']; % Raw data localisation
MatfileFolder = [GeneralFolder 'Data_Joseph\Matfiles']; % for saving data in .mat
% SecondarySaveFolder = ''; % secondary spot for backup saving
FigureFolder = [GeneralFolder 'Data_Joseph\Figures']; % folder for saving figures

PLOT = 1


% Paramètres de couleurs par drogues
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

% Automatic mode
AUTO = 0;

% Autres paramètres

% pour v2d
NotSaved = {};

%pour d2p
PLOTDFCurve = 1;
PLOTDFComp = 1;

ExperimentalConditions3T3.Objectif = '100X'; % '63X'; % '100X'
ExperimentalConditions3T3.Scale = 1/15.8; % 1/9.7; % 15.8
ExperimentalConditions3T3.MagFieldCorrectionCoeff = 1.15;
ExperimentalConditions3T3.OpticalDepthCorrectionCoeff = 1.33/1.52;
ExperimentalConditions3T3.TypeCellulaire = '3T3';
ExperimentalConditions3T3.TypeBilles = 'M450';
ExperimentalConditions3T3.DiametreBilles = 4503;
ExperimentalConditions3T3.SubstrateCoating = 'BSA';
ExperimentalConditions3T3.Drug = 'No';
ExperimentalConditions3T3.NormalField = 10;
ExperimentalConditions3T3.RampField = [3 40];

ExperimentalConditionsDictysPLL.Objectif = '100X'; % '63X'; % '100X'
ExperimentalConditionsDictysPLL.Scale = 1/15.8; % 1/9.7; % 15.8
ExperimentalConditionsDictysPLL.MagFieldCorrectionCoeff = 1.15;
ExperimentalConditionsDictysPLL.OpticalDepthCorrectionCoeff = 1.33/1.52;
ExperimentalConditionsDictysPLL.TypeCellulaire = 'dictys';
ExperimentalConditionsDictysPLL.TypeBilles = 'M450';
ExperimentalConditionsDictysPLL.DiametreBilles = 4503;
ExperimentalConditionsDictysPLL.SubstrateCoating = 'BSA';
ExperimentalConditionsDictysPLL.Drug = 'No';
ExperimentalConditionsDictysPLL.NormalField = 5;
ExperimentalConditionsDictysPLL.RampField = [3 40];

ExperimentalConditionsDictysPLL2 = ExperimentalConditionsDictysPLL;
ExperimentalConditionsDictysPLL2.RampField = [5 40];

ExperimentalConditionsDictysBSA.Objectif = '100X'; % '63X'; % '100X'
ExperimentalConditionsDictysBSA.Scale = 1/15.8; % 1/9.7; % 15.8
ExperimentalConditionsDictysBSA.MagFieldCorrectionCoeff = 1.1;
ExperimentalConditionsDictysBSA.OpticalDepthCorrectionCoeff = 1.33/1.52;
ExperimentalConditionsDictysBSA.TypeCellulaire = 'dictys';
ExperimentalConditionsDictysBSA.TypeBilles = 'M450';
ExperimentalConditionsDictysBSA.DiametreBilles = 4503;
ExperimentalConditionsDictysBSA.SubstrateCoating = 'BSA';
ExperimentalConditionsDictysBSA.Drug = 'No';
ExperimentalConditionsDictysBSA.NormalField = 5;
ExperimentalConditionsDictysBSA.RampField = [3 40];

ExperimentalConditionsDictysBSA2 = ExperimentalConditionsDictysBSA;
ExperimentalConditionsDictysBSA2.NormalField = 3;
ExperimentalConditionsDictysBSA2.RampField = [1 40];

ExperimentalConditionsFloat.Objectif = '63X'; % '63X'; % '100X'
ExperimentalConditionsFloat.Scale = 1/9.7; % 1/9.7; % 15.8
ExperimentalConditionsFloat.MagFieldCorrectionCoeff = 1.15;
ExperimentalConditionsFloat.OpticalDepthCorrectionCoeff = 1.4/1.0;
ExperimentalConditionsFloat.TypeCellulaire = 'dictys';
ExperimentalConditionsFloat.TypeBilles = 'M450';
ExperimentalConditionsFloat.DiametreBilles = 4503;
ExperimentalConditionsFloat.SubstrateCoating = 'Float';
ExperimentalConditionsFloat.Drug = 'No';
ExperimentalConditionsFloat.NormalField = 5;
ExperimentalConditionsFloat.RampField = [3 40];

EXCLUDED = {'04-08-20_M2_P1_C3', '04-08-20_M2_P1_C4', '05-08-20_M1_P1_C6', '05-08-20_M2_P1_C4', '05-08-20_M2_P1_C6', '05-08-20_M2_P1_C7', '07-08-20_M1_P1_C6', '07-08-20_M1_P1_C62', '07-08-20_M2_P1_C7'};
EXCLUDED = [EXCLUDED {'06-08-20_M1_P1_C4', '06-08-20_M1_P2_C1', '06-08-20_M1_P2_C3', '06-08-20_M1_P2_C4', '06-08-20_M1_P3_C1', '06-08-20_M1_P3_C5'}];

%%
if 0
    %% Nouvelle Manip Rampe  %%% R2V
    
    
    %% Joseph
    
    % 3T3
    Res2Var_wFluo_multiZ_Comp_J('Results','07-08-20','M1','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,1,...
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','07-08-20','M2','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,1,...
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','05-08-20','M1','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,1,...
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','05-08-20','M2','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,1,...
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','04-08-20','M2','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,1,...
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','04-08-20','M1','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,1,...
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','24-07-20','M2','R403T3aSFLdoxi','3T3aSFL_Pll_doxi',ExperimentalConditions3T3,1,...
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','24-07-20','M1','R403T3aSFLnodrugs','3T3aSFL_Pll_nodrugs',ExperimentalConditions3T3,1,...
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    % Dictys
    Res2Var_wFluo_multiZ_Comp_J('Results','06-08-20','M1','R40float','DictyComp_Float-7-80',ExperimentalConditionsFloat,1:3,...
        '08-07-20_Depthograph63x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','08-07-20','M1','R40float','DictyComp_Float-30-80',ExperimentalConditionsFloat,1,...
        '08-07-20_Depthograph63x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','26-06-20','M1','R40PLL','DictyComp_PLL-30opti',ExperimentalConditionsDictysPLL,1,...
        '15-06-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','23-06-20','M1','R40float','DictyComp_Float-30-80',ExperimentalConditionsFloat,1,...
        '15-06-20_Depthograph63x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    %% Valentin
    % BeadsIn
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
    
    
    %     % Beads out
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
    Res2Var_wFluo_multiZ_Comp('Results','25-09-18','R80','M2','DC-Comp_CK666',2,...
        '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    % Arpin 2 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','05-10-18','R80','M1','DC-Comp_ArpinKO',1:3,...
        '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','05-10-18','R80','M2','DC-Comp_ArpinWT',1:3,...
        '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    % CultureRampe Inhib 2 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M1','DC-Comp_CalA25',1,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M1','DC-Comp_LatA2',2,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M2','DC-Comp_LatA2',2,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M1','DC-Comp_Blebbi',2,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M2','DC-Comp_CalA25',1,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M2','DC-Comp_Ctrl',2,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    % CultureRampe Inhib 3 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','12-12-18','R80','M1','DC-Comp_Ctrl',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','12-12-18','R80','M2','DC-Comp_Nase',1:3,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    % CultureRampe ArpC4 1 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','17-12-18','R80','M1','DC-Comp_AprC4KO',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','17-12-18','R80','M2','DC-Comp_AprC4WT',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    % CultureRampe Inhib 4 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','18-12-18','R80','M1','DC-Comp_Y27',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','18-12-18','R80','M2','DC-Comp_Nase',2,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    % CultureRampe Inhib 5 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','21-12-18','R80','M1','DC-Comp_LatANase',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    % CultureRampe Inhib 5 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','30-04-19','R80','M1','DC-Comp_SiCtrl+LatA',2,...
        '26-06-19_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','30-04-19','R80','M1','DC-Comp_SiVim+LatA',1,...
        '26-06-19_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    %% Manip avec differente temps de compression
    Res2Var_wFluo_multiZ_Comp_R248('Results','14-10-19','R248','M1','DC-Comp_R248_Ctrl',1:2,...
        '15-10-19_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp_R248('Results','14-10-19','R248','M2','DC-Comp_R248_LatA',1:2,...
        '15-10-19_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp_R248('Results','15-10-19','R248','M1','DC-Comp_R248_Blebbi',1:2,...
        '15-10-19_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp_R248('Results','15-10-19','R248','M2','DC-Comp_R248_Ctrl',1:2,...
        '15-10-19_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp_R248('Results','13-02-20','R248','M1','DictyAx2-Comp',1:2,...
        '10-02-20_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
       Res2Var_wFluo_multiZ_Comp_R248('Results','20-02-20','R248','M1','DictyAx2-Comp',1,...
        '20-02-20_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    %% V2D comp
    
    %% Joseph
    
    % 3T3
    Var2Data_wFluo_Comp_J('07-08-20','M1','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('07-08-20','M2','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('05-08-20','M1','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('05-08-20','M2','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('04-08-20','M2','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('04-08-20','M1','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('24-07-20','M2','R403T3aSFLdoxi','3T3aSFL_Pll_doxi',ExperimentalConditions3T3,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('24-07-20','M1','R403T3aSFLnodrugs','3T3aSFL_Pll_nodrugs',ExperimentalConditions3T3,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    % Dictys
    Var2Data_wFluo_Comp_J('06-08-20','M1','R40float','DictyComp_Float-7-80',ExperimentalConditionsFloat,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('08-07-20','M1','R40float','DictyComp_Float-30-80',ExperimentalConditionsFloat,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('26-06-20','M1','R40PLL','DictyComp_PLL-30opti',ExperimentalConditionsDictysPLL,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('23-06-20','M1','R40float','DictyComp_Float-30-80',ExperimentalConditionsFloat,95,143,RawdataFolder,MatfileFolder,FigureFolder)
    
    %% manip inhib
    
    % manip ctrl test
    Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4415);
    Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4415);
    
    % Arpin 1
    Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M1','DC-Comp_ArpinWT',95,179,4438);
    Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M2','DC-Comp_ArpinKO',95,179,4438);
    
    Var2Data_wFluo_Comp('18-09-18','R80',0.56,'M1','DC-Comp_ArpinKO',95,179,4438);
    Var2Data_wFluo_Comp('18-09-18','R80',0.56,'M2','DC-Comp_ArpinWT',95,179,4438);
    Var2Data_wFluo_Comp('18-09-18','R80',0.56,'M3','DC-Comp_ArpinKO',95,179,4438);
    
    % Inhib 1
    
    Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_CK666',95,179,4438);
    Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4438);
    Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M2','DC-Comp_SMIFH2',95,179,4438);
    
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_LatA500',95,179,4438);
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_SMIFH2',95,179,4438);
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4438);
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_CK666',95,179,4438);
    
    % Arpin 2;
    Var2Data_wFluo_Comp('05-10-18','R80',0.56,'M1','DC-Comp_ArpinKO',95,179,4438);
    Var2Data_wFluo_Comp('05-10-18','R80',0.56,'M2','DC-Comp_ArpinWT',95,179,4438);
    
    
    % Inhib 2
    %
    %     Var2Data_wFluo_Comp('29-10-18','R80',0.56,'M1','DC-Comp_CalA25',95,179,4438);
    %     Var2Data_wFluo_Comp('29-10-18','R80',0.56,'M1','DC-Comp_LatA2',95,179,4438);
    %     Var2Data_wFluo_Comp('29-10-18','R80',0.56,'M2','DC-Comp_LatA2',95,179,4438);
    
    
    Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Blebbi',95,179,4438)
    %     Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M1','DC-Comp_CalA25',95,179,4438);
    Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4438);
    
    % inhib 3
    
    
    Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4438);
    Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M2','DC-Comp_Nase',95,179,4438);
    %
    %
    %     % ArpC4
    %
    Var2Data_wFluo_Comp('17-12-18','R80',0.56,'M1','DC-Comp_AprC4KO',95,179,4438);
    Var2Data_wFluo_Comp('17-12-18','R80',0.56,'M2','DC-Comp_AprC4WT',95,179,4438);
    %
    %     % inhib 4
    %
    Var2Data_wFluo_Comp('18-12-18','R80',0.56,'M1','DC-Comp_Y27',95,179,4438);
    Var2Data_wFluo_Comp('18-12-18','R80',0.56,'M2','DC-Comp_Nase',95,179,4438);
    %
    %     % Inhib 5
    %
    Var2Data_wFluo_Comp('21-12-18','R80',0.56,'M1','DC-Comp_LatANase',95,179,4438);
    
    %%%%%%%%% billes internes %%%%%%%%%
    %     Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4415);
    %     Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4415);
    %
    %
    %     Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438);
    %     Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438);
    %
    %
    %     Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438);
    %     Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438);
    %
    %
    %     Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438);
    %     Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438);
    
    
    %     Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438);
    %     Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M2','DC-Comp_NaseInside',95,179,4438);
    %
    %     Var2Data_wFluo_Comp('18-12-18','R80',0.56,'M2','DC-Comp_NaseInside',95,179,4438);
    %
    %
    %     Var2Data_wFluo_Comp('30-04-19','R80',0.56,'M1','DC-Comp_SiCtrl+LatA',95,179,4438);
    %     Var2Data_wFluo_Comp('30-04-19','R80',0.56,'M1','DC-Comp_SiVim+LatA',95,179,4438) ;
    
    %     % billes externes
    %
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Outside',95,179,4438)
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Outside',95,179,4438)
    
    Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Outside',95,179,4438)
    Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Outside',95,179,4438)
    
    %     Var2Data_wFluo_Comp('21-01-19','R80',0.56,'M1','DC-Comp_Bsa',95,179,4438);
    %
    dist = []
    dist = [dist Var2Data_wFluo_Comp('23-01-19','R80',0.56,'M1','DC-Comp_SerumJ1',95,179,4438)];
    %     Var2Data_wFluo_Comp('23-01-19','R80',0.56,'M1','DC-Comp_Naked',95,179,4438);
    %
    %
    dist = [dist Var2Data_wFluo_Comp('28-01-19','R80',0.56,'M1','DC-Comp_SerumJ6',95,179,4438)];
    %
    
    %% compression de tps diff
    
    Var2Data_wFluo_Comp_R248('14-10-19','R248',0.56,'M1','DC-Comp_R248_Ctrl',21,24,98,193,383,4438);
    Var2Data_wFluo_Comp_R248('14-10-19','R248',0.56,'M2','DC-Comp_R248_LatA',21,24,98,193,383,4438);
    Var2Data_wFluo_Comp_R248('15-10-19','R248',0.56,'M1','DC-Comp_R248_Blebbi',21,24,98,193,383,4438);
    Var2Data_wFluo_Comp_R248('15-10-19','R248',0.56,'M2','DC-Comp_R248_Ctrl',21,24,98,193,383,4438);
    
    
    Var2Data_wFluo_Comp_R248('13-02-20','R248',1.1,'M1','DictyAx2-Comp',21,24,98,193,383,4438);
    Var2Data_wFluo_Comp_R248('20-02-20','R248',1.1,'M1','DictyAx2-Comp',21,24,98,193,383,4438);
    
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
    
    %% D2M Joseph
    PLOT = 1; % PLOT = 0;
    
    % 3T3
    PLOT = 0;
    Data2Meca_J('07-08-20','M1','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('07-08-20','M2','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    
    Data2Meca_J('05-08-20','M1','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('05-08-20','M2','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    
    Data2Meca_J('04-08-20','M2','R40_3T3aSFLdoxi','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('04-08-20','M1','R40_3T3aSFLnodrugs','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);

    Data2Meca_J('24-07-20','M2','R403T3aSFLdoxi','3T3aSFL_Pll_doxi',ExperimentalConditions3T3,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('24-07-20','M1','R403T3aSFLnodrugs','3T3aSFL_Pll_nodrugs',ExperimentalConditions3T3,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    
    % Dictys
    PLOT = 0;
    Data2Meca_J('27-05-20','M1','R40','DictyAx2_DMSO',ExperimentalConditionsDictysBSA2,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('25-05-20','M1','R40','DictyAx2_DMSO',ExperimentalConditionsDictysBSA,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('22-05-20','M1','R40','DictyAx2_DMSO',ExperimentalConditionsDictysBSA,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('20-05-20','M2','R40','DictyAx2_DMSO',ExperimentalConditionsDictysBSA,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('20-05-20','M1','R40','DictyAx2_DMSO',ExperimentalConditionsDictysBSA,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    
    
    Data2Meca_J('06-08-20','M1','R40float','DictyComp_Float-7-80',ExperimentalConditionsFloat,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('08-07-20','M1','R40float','DictyComp_Float-30-80',ExperimentalConditionsFloat,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('23-06-20','M1','R40float','DictyComp_Float-30-80',ExperimentalConditionsFloat,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    
    Data2Meca_J('04-08-20','M2','R40PLL','DictyAx2_PLL-7-5opti',ExperimentalConditionsDictysPLL2,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('24-07-20','M1','R40PLLJ','DictyAx2_PLLJ-7-5opti',ExperimentalConditionsDictysPLL,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    Data2Meca_J('26-06-20','M1','R40PLLJ','DictyComp_PLLJ-30opti',ExperimentalConditionsDictysPLL,PLOT,RawdataFolder,MatfileFolder,FigureFolder,EXCLUDED);
    
    
    
    %% M2P Joseph

    % 3T3
    Meca2Plot_J('3T3aSFL_BSA_doxi',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('3T3aSFL_BSA_nodrugs',Cct,MatfileFolder,FigureFolder)
   
    Meca2Plot_J('3T3aSFL_Pll_doxi',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('3T3aSFL_Pll_nodrugs',Cct,MatfileFolder,FigureFolder)
    
    % Dictys
    Meca2Plot_J('DictyAx2_DMSO',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('DictyAx2_PLL-7-5opti',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('DictyAx2_PLLJ-7-5opti',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('DictyComp_Float-7-80',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('DictyComp_Float-30-80',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('DictyComp_PLL-30opti',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('DictyComp_Float-30-80',Cct,MatfileFolder,FigureFolder)
    
%     Meca2Plot('Comp_CK666',Cck)
%     Meca2Plot('Comp_SMIFH2',Cs)
%     Meca2Plot('Comp_Blebbi',Cb)
%     Meca2Plot('Comp_LatA500',Cl)
    
    
    %% D2M pour these et papier
    PLOT = 0;
    
    [nComp(1), nCompOK(1)] = Data2Meca('28-08-18','M1','DC-Comp_Ctrl',4415,1);
    [nComp(2), nCompOK(2)] = Data2Meca('28-08-18','M2','DC-Comp_Ctrl',4415,PLOT);
    [nComp(3), nCompOK(3)] = Data2Meca('24-09-18','M1','DC-Comp_Ctrl',4438,PLOT);
    [nComp(4), nCompOK(4)] = Data2Meca('25-09-18','M2','DC-Comp_Ctrl',4438,PLOT);
    [nComp(5), nCompOK(5)] = Data2Meca('30-10-18','M2','DC-Comp_Ctrl',4438,PLOT);
    [nComp(6), nCompOK(6)] = Data2Meca('12-12-18','M1','DC-Comp_Ctrl',4438,PLOT);
    
    
    Data2Meca('24-09-18','M1','DC-Comp_CK666',4438,PLOT);
    Data2Meca('25-09-18','M2','DC-Comp_CK666',4438,PLOT);
    
    Data2Meca('24-09-18','M2','DC-Comp_SMIFH2',4438,PLOT) ;
    Data2Meca('25-09-18','M1','DC-Comp_SMIFH2',4438,PLOT);
    
    
    Data2Meca('30-10-18','M1','DC-Comp_Blebbi',4438,PLOT);
    
    
    [a,b] =  Data2Meca('25-09-18','M1','DC-Comp_LatA500',4438,PLOT);
    
    
    
    Meca2Plot('Comp_Ctrl',Cct)
    Meca2Plot('Comp_CK666',Cck)
    Meca2Plot('Comp_SMIFH2',Cs)
    Meca2Plot('Comp_Blebbi',Cb)
    Meca2Plot('Comp_LatA500',Cl)
    
    %% D2M pour comp differente longueur
    
    % DC
    Data2Meca_R248('14-10-19','M1','DC-Comp_R248_Ctrl',4438,0)
    Data2Meca_R248('14-10-19','M2','DC-Comp_R248_LatA',4438,1)
    Data2Meca_R248('15-10-19','M1','DC-Comp_R248_Blebbi',4438,1)
    Data2Meca_R248('15-10-19','M2','DC-Comp_R248_Ctrl',4438,0)
    
    Meca2Plot_R248('Comp_R248_Ctrl',Cct,1)
    Meca2Plot_R248('Comp_R248_LatA',Cl,1)
    Meca2Plot_R248('Comp_R248_Blebbi',Cb,1)
    
    % Dicty
    
    Data2Meca_R248('13-02-20','M1','DictyAx2-Comp',4438,0)    
    Data2Meca_R248('20-02-20','M1','DictyAx2-Comp',4438,0)
    
    
    Meca2Plot_R248('DictyAx2-Comp',Cct,1)
end

if 0
    %% V2D sur no comp
    % data from compression
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-08-18','R80cst',0.6,'M1','DC-Comp_Ctrl',4415)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-08-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4415)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18','R80cst',0.6,'M1','DC-Comp_ArpinWT',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18','R80cst',0.6,'M2','DC-Comp_ArpinKO',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M1','DC-Comp_ArpinKO',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M2','DC-Comp_ArpinWT',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M3','DC-Comp_ArpinKO',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M2','DC-Comp_SMIFH2',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M1','DC-Comp_CK666',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M1','DC-Comp_Ctrl',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M1','DC-Comp_LatA500',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M1','DC-Comp_SMIFH2',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M2','DC-Comp_CK666',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-10-18','R80cst',0.6,'M1','DC-Comp_ArpinKO',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-10-18','R80cst',0.6,'M2','DC-Comp_ArpinWT',4438)};
    
    % NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M1','DC-Comp_CalA25',4438)};
    % NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M1','DC-Comp_LatA2',4438)};
    % NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M2','DC-Comp_LatA2',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M1','DC-Comp_Blebbi',4438)};
    % NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M2','DC-Comp_CalA25',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4438)};
    
    
    %%%%%%%%% billes internes %%%%%%%%%
    Var2Data_wFluo('28-08-18','R80cst',0.56,'M1','DC-Comp_Inside',4415,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('28-08-18','R80cst',0.56,'M2','DC-Comp_Inside',4415,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('17-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('17-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('24-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('24-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('25-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('25-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('12-12-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('12-12-18','R80cst',0.56,'M2','DC-Comp_NaseInside',4438,...
        RawdataFolder,MatfileFolder)
    
    Var2Data_wFluo('18-12-18','R80cst',0.56,'M2','DC-Comp_NaseInside',4438,...
        RawdataFolder,MatfileFolder)
    
    % billes externes
    
    Var2Data_wFluo('25-09-18','R80cst',0.56,'M1','DC-Comp_Outside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('25-09-18','R80cst',0.56,'M2','DC-Comp_Outside',4438,...
        RawdataFolder,MatfileFolder)
    
    Var2Data_wFluo('30-10-18','R80cst',0.56,'M1','DC-Comp_Outside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('30-10-18','R80cst',0.56,'M2','DC-Comp_Outside',4438,...
        RawdataFolder,MatfileFolder)
    
    Var2Data_wFluo('21-01-19','R80cst',0.56,'M1','DC-Comp_Bsa',4438,...
        RawdataFolder,MatfileFolder)
    
    Var2Data_wFluo('23-01-19','R80cst',0.56,'M1','DC-Comp_Naked',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('23-01-19','R80cst',0.56,'M1','DC-Comp_SerumJ1',4438,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('28-01-19','R80cst',0.56,'M1','DC-Comp_SerumJ6',4438,...
        RawdataFolder,MatfileFolder)
    
    
    % data from comp
    
    Data2Peaks({'M1','M2'},'DC-Comp_Inside',0,0.9,MatfileFolder,FigureFolder);
    Data2Peaks({'M1','M2'},'DC-Comp_Outside',1,0.9,MatfileFolder,FigureFolder);
    
    Data2Peaks({'M1','M2'},'DC-Comp_Naked',0,0.9,MatfileFolder,FigureFolder);
    
    Data2Peaks({'M1','M2'},'DC-Comp_SerumJ6',0,0.9,MatfileFolder,FigureFolder);
    
    
    % data from comp
    PLOTPS = 0;
    Peaks2Stats('DC-Comp_Inside',PLOTPS,MatfileFolder,FigureFolder);
    Peaks2Stats('DC-Comp_Outside',PLOTPS,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DC-Comp_Naked');
    Peaks2Stats('DC-Comp_SerumJ6');
    
    FigureMaker(0,MatfileFolder,FigureFolder)
    
end


 %% OLD VERSION
 if 0
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

% pour v2d
NotSaved = {};

%pour d2p
PLOTDFCurve = 1;
PLOTDFComp = 1;
 end
%%
if 0
    %% Nouvelle Manip Rampe  %%% R2V
  
    %% Joseph
    %% 3T3
    Res2Var_wFluo_multiZ_Comp_J('Results','04-08-20','R40','M1','aSFLdoxi',2,... % M1_P2 = doxi
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','04-08-20','R40','M1','aSFLnodrug',1,... % M1_P1 = nodrug
        '20-02-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    %% Dictys
    Res2Var_wFluo_multiZ_Comp_J('Results','08-07-20','R40float','M1','DictyComp_Float-30-80',1,...
        '08-07-20_Depthograph63x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','26-06-20','R40PLL','M1','DictyComp_PLL-30opti',1,...
        '15-06-20_Depthograph100x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    Res2Var_wFluo_multiZ_Comp_J('Results','23-06-20','R40float','M1','DictyComp_Float-30-80',1,...
        '15-06-20_Depthograph63x',24,95,143,AUTO,RawdataFolder,MatfileFolder)
    
    %% Valentin
    % BeadsIn
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
    
    
    %     % Beads out
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
    Res2Var_wFluo_multiZ_Comp('Results','25-09-18','R80','M2','DC-Comp_CK666',2,...
        '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    % Arpin 2 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','05-10-18','R80','M1','DC-Comp_ArpinKO',1:3,...
        '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','05-10-18','R80','M2','DC-Comp_ArpinWT',1:3,...
        '28-08-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    % CultureRampe Inhib 2 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M1','DC-Comp_CalA25',1,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M1','DC-Comp_LatA2',2,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','29-10-18','R80','M2','DC-Comp_LatA2',2,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M1','DC-Comp_Blebbi',2,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M2','DC-Comp_CalA25',1,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','30-10-18','R80','M2','DC-Comp_Ctrl',2,...
        '30-10-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    % CultureRampe Inhib 3 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','12-12-18','R80','M1','DC-Comp_Ctrl',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','12-12-18','R80','M2','DC-Comp_Nase',1:3,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    % CultureRampe ArpC4 1 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','17-12-18','R80','M1','DC-Comp_AprC4KO',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','17-12-18','R80','M2','DC-Comp_AprC4WT',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    % CultureRampe Inhib 4 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','18-12-18','R80','M1','DC-Comp_Y27',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','18-12-18','R80','M2','DC-Comp_Nase',2,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    % CultureRampe Inhib 5 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','21-12-18','R80','M1','DC-Comp_LatANase',1,...
        '12-12-18_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    % CultureRampe Inhib 5 18s (42 + 42 img) + 2s (95img) = 179
    Res2Var_wFluo_multiZ_Comp('Results','30-04-19','R80','M1','DC-Comp_SiCtrl+LatA',2,...
        '26-06-19_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    Res2Var_wFluo_multiZ_Comp('Results','30-04-19','R80','M1','DC-Comp_SiVim+LatA',1,...
        '26-06-19_Depthograph',42,95,179,AUTO,RawdataFolder,MatfileFolder)
    
    
    %% Manip avec differente temps de compression
    Res2Var_wFluo_multiZ_Comp_R248('Results','14-10-19','R248','M1','DC-Comp_R248_Ctrl',1:2,...
        '15-10-19_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp_R248('Results','14-10-19','R248','M2','DC-Comp_R248_LatA',1:2,...
        '15-10-19_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp_R248('Results','15-10-19','R248','M1','DC-Comp_R248_Blebbi',1:2,...
        '15-10-19_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp_R248('Results','15-10-19','R248','M2','DC-Comp_R248_Ctrl',1:2,...
        '15-10-19_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    
    Res2Var_wFluo_multiZ_Comp_R248('Results','13-02-20','R248','M1','DictyAx2-Comp',1:2,...
        '10-02-20_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
       Res2Var_wFluo_multiZ_Comp_R248('Results','20-02-20','R248','M1','DictyAx2-Comp',1,...
        '20-02-20_Depthograph',21,24,98,193,383,...
        AUTO,RawdataFolder,MatfileFolder)
    
    %% V2D comp
    
    %% Joseph CHECK L'ECHELLE !!!!!
    Var2Data_wFluo_Comp_J('08-07-20','R40float',1.1,'M1','DictyComp_Float-30-80',95,143,4503,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('26-06-20','R40PLL',1.1,'M1','DictyComp_PLL-30opti',95,143,4503,RawdataFolder,MatfileFolder,FigureFolder)
    
    Var2Data_wFluo_Comp_J('23-06-20','R40float',1.1,'M1','DictyComp_Float-30-80',95,143,4503,RawdataFolder,MatfileFolder,FigureFolder)
    
    %% manip inhib
    
    % manip ctrl test
    Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4415);
    Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4415);
    
    % Arpin 1
    Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M1','DC-Comp_ArpinWT',95,179,4438);
    Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M2','DC-Comp_ArpinKO',95,179,4438);
    
    Var2Data_wFluo_Comp('18-09-18','R80',0.56,'M1','DC-Comp_ArpinKO',95,179,4438);
    Var2Data_wFluo_Comp('18-09-18','R80',0.56,'M2','DC-Comp_ArpinWT',95,179,4438);
    Var2Data_wFluo_Comp('18-09-18','R80',0.56,'M3','DC-Comp_ArpinKO',95,179,4438);
    
    % Inhib 1
    
    Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_CK666',95,179,4438);
    Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4438);
    Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M2','DC-Comp_SMIFH2',95,179,4438);
    
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_LatA500',95,179,4438);
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_SMIFH2',95,179,4438);
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4438);
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_CK666',95,179,4438);
    
    % Arpin 2;
    Var2Data_wFluo_Comp('05-10-18','R80',0.56,'M1','DC-Comp_ArpinKO',95,179,4438);
    Var2Data_wFluo_Comp('05-10-18','R80',0.56,'M2','DC-Comp_ArpinWT',95,179,4438);
    
    
    % Inhib 2
    %
    %     Var2Data_wFluo_Comp('29-10-18','R80',0.56,'M1','DC-Comp_CalA25',95,179,4438);
    %     Var2Data_wFluo_Comp('29-10-18','R80',0.56,'M1','DC-Comp_LatA2',95,179,4438);
    %     Var2Data_wFluo_Comp('29-10-18','R80',0.56,'M2','DC-Comp_LatA2',95,179,4438);
    
    
    Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Blebbi',95,179,4438)
    %     Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M1','DC-Comp_CalA25',95,179,4438);
    Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4438);
    
    % inhib 3
    
    
    Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4438);
    Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M2','DC-Comp_Nase',95,179,4438);
    %
    %
    %     % ArpC4
    %
    Var2Data_wFluo_Comp('17-12-18','R80',0.56,'M1','DC-Comp_AprC4KO',95,179,4438);
    Var2Data_wFluo_Comp('17-12-18','R80',0.56,'M2','DC-Comp_AprC4WT',95,179,4438);
    %
    %     % inhib 4
    %
    Var2Data_wFluo_Comp('18-12-18','R80',0.56,'M1','DC-Comp_Y27',95,179,4438);
    Var2Data_wFluo_Comp('18-12-18','R80',0.56,'M2','DC-Comp_Nase',95,179,4438);
    %
    %     % Inhib 5
    %
    Var2Data_wFluo_Comp('21-12-18','R80',0.56,'M1','DC-Comp_LatANase',95,179,4438);
    
    %%%%%%%%% billes internes %%%%%%%%%
    %     Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4415);
    %     Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4415);
    %
    %
    %     Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438);
    %     Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438);
    %
    %
    %     Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438);
    %     Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438);
    %
    %
    %     Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438);
    %     Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438);
    
    
    %     Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438);
    %     Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M2','DC-Comp_NaseInside',95,179,4438);
    %
    %     Var2Data_wFluo_Comp('18-12-18','R80',0.56,'M2','DC-Comp_NaseInside',95,179,4438);
    %
    %
    %     Var2Data_wFluo_Comp('30-04-19','R80',0.56,'M1','DC-Comp_SiCtrl+LatA',95,179,4438);
    %     Var2Data_wFluo_Comp('30-04-19','R80',0.56,'M1','DC-Comp_SiVim+LatA',95,179,4438) ;
    
    %     % billes externes
    %
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Outside',95,179,4438)
    Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Outside',95,179,4438)
    
    Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Outside',95,179,4438)
    Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Outside',95,179,4438)
    
    %     Var2Data_wFluo_Comp('21-01-19','R80',0.56,'M1','DC-Comp_Bsa',95,179,4438);
    %
    dist = []
    dist = [dist Var2Data_wFluo_Comp('23-01-19','R80',0.56,'M1','DC-Comp_SerumJ1',95,179,4438)];
    %     Var2Data_wFluo_Comp('23-01-19','R80',0.56,'M1','DC-Comp_Naked',95,179,4438);
    %
    %
    dist = [dist Var2Data_wFluo_Comp('28-01-19','R80',0.56,'M1','DC-Comp_SerumJ6',95,179,4438)];
    %
    
    %% compression de tps diff
    
    Var2Data_wFluo_Comp_R248('14-10-19','R248',0.56,'M1','DC-Comp_R248_Ctrl',21,24,98,193,383,4438);
    Var2Data_wFluo_Comp_R248('14-10-19','R248',0.56,'M2','DC-Comp_R248_LatA',21,24,98,193,383,4438);
    Var2Data_wFluo_Comp_R248('15-10-19','R248',0.56,'M1','DC-Comp_R248_Blebbi',21,24,98,193,383,4438);
    Var2Data_wFluo_Comp_R248('15-10-19','R248',0.56,'M2','DC-Comp_R248_Ctrl',21,24,98,193,383,4438);
    
    
    Var2Data_wFluo_Comp_R248('13-02-20','R248',1.1,'M1','DictyAx2-Comp',21,24,98,193,383,4438);
    Var2Data_wFluo_Comp_R248('20-02-20','R248',1.1,'M1','DictyAx2-Comp',21,24,98,193,383,4438);
    
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
    
    %% D2M Joseph
    PLOT = 1;
%     PLOT = 0;
    Data2Meca_J('08-07-20','M1','DictyComp_Float-30-80','R40float',4503,PLOT,RawdataFolder,MatfileFolder,FigureFolder);
    
    Data2Meca_J('26-06-20','M1','DictyComp_PLL-30opti','R40PLL',4503,PLOT,RawdataFolder,MatfileFolder,FigureFolder);
    
    Data2Meca_J('23-06-20','M1','DictyComp_Float-30-80','R40float',4503,PLOT,RawdataFolder,MatfileFolder,FigureFolder);
    
    %% M2P Joseph

    Meca2Plot_J('DictyComp_Float-30-80',Cct,MatfileFolder,FigureFolder)
%     Meca2Plot_V('DictyComp_Float-30-80',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('DictyComp_PLL-30opti',Cct,MatfileFolder,FigureFolder)
    
    Meca2Plot_J('DictyComp_Float-30-80',Cct,MatfileFolder,FigureFolder)
%     Meca2Plot('Comp_CK666',Cck)
%     Meca2Plot('Comp_SMIFH2',Cs)
%     Meca2Plot('Comp_Blebbi',Cb)
%     Meca2Plot('Comp_LatA500',Cl)
    
    
    %% D2M pour these et papier
    PLOT = 0;
    
    [nComp(1), nCompOK(1)] = Data2Meca('28-08-18','M1','DC-Comp_Ctrl',4415,1);
    [nComp(2), nCompOK(2)] = Data2Meca('28-08-18','M2','DC-Comp_Ctrl',4415,PLOT);
    [nComp(3), nCompOK(3)] = Data2Meca('24-09-18','M1','DC-Comp_Ctrl',4438,PLOT);
    [nComp(4), nCompOK(4)] = Data2Meca('25-09-18','M2','DC-Comp_Ctrl',4438,PLOT);
    [nComp(5), nCompOK(5)] = Data2Meca('30-10-18','M2','DC-Comp_Ctrl',4438,PLOT);
    [nComp(6), nCompOK(6)] = Data2Meca('12-12-18','M1','DC-Comp_Ctrl',4438,PLOT);
    
    
    Data2Meca('24-09-18','M1','DC-Comp_CK666',4438,PLOT);
    Data2Meca('25-09-18','M2','DC-Comp_CK666',4438,PLOT);
    
    Data2Meca('24-09-18','M2','DC-Comp_SMIFH2',4438,PLOT) ;
    Data2Meca('25-09-18','M1','DC-Comp_SMIFH2',4438,PLOT);
    
    
    Data2Meca('30-10-18','M1','DC-Comp_Blebbi',4438,PLOT);
    
    
    [a,b] =  Data2Meca('25-09-18','M1','DC-Comp_LatA500',4438,PLOT);
    
    
    
    Meca2Plot('Comp_Ctrl',Cct)
    Meca2Plot('Comp_CK666',Cck)
    Meca2Plot('Comp_SMIFH2',Cs)
    Meca2Plot('Comp_Blebbi',Cb)
    Meca2Plot('Comp_LatA500',Cl)
    
    %% D2M pour comp differente longueur
    
    % DC
    Data2Meca_R248('14-10-19','M1','DC-Comp_R248_Ctrl',4438,0)
    Data2Meca_R248('14-10-19','M2','DC-Comp_R248_LatA',4438,1)
    Data2Meca_R248('15-10-19','M1','DC-Comp_R248_Blebbi',4438,1)
    Data2Meca_R248('15-10-19','M2','DC-Comp_R248_Ctrl',4438,0)
    
    Meca2Plot_R248('Comp_R248_Ctrl',Cct,1)
    Meca2Plot_R248('Comp_R248_LatA',Cl,1)
    Meca2Plot_R248('Comp_R248_Blebbi',Cb,1)
    
    % Dicty
    
    Data2Meca_R248('13-02-20','M1','DictyAx2-Comp',4438,0)    
    Data2Meca_R248('20-02-20','M1','DictyAx2-Comp',4438,0)
    
    
    Meca2Plot_R248('DictyAx2-Comp',Cct,1)
end

if 0
    %% V2D sur no comp
    % data from compression
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-08-18','R80cst',0.6,'M1','DC-Comp_Ctrl',4415)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('28-08-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4415)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18','R80cst',0.6,'M1','DC-Comp_ArpinWT',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('17-09-18','R80cst',0.6,'M2','DC-Comp_ArpinKO',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M1','DC-Comp_ArpinKO',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M2','DC-Comp_ArpinWT',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('18-09-18','R80cst',0.6,'M3','DC-Comp_ArpinKO',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M2','DC-Comp_SMIFH2',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M1','DC-Comp_CK666',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('24-09-18','R80cst',0.6,'M1','DC-Comp_Ctrl',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M1','DC-Comp_LatA500',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M1','DC-Comp_SMIFH2',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M2','DC-Comp_CK666',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('25-09-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-10-18','R80cst',0.6,'M1','DC-Comp_ArpinKO',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('05-10-18','R80cst',0.6,'M2','DC-Comp_ArpinWT',4438)};
    
    % NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M1','DC-Comp_CalA25',4438)};
    % NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M1','DC-Comp_LatA2',4438)};
    % NotSaved = {NotSaved{:},Var2Data_wFluo('29-10-18','R80cst',0.6,'M2','DC-Comp_LatA2',4438)};
    
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M1','DC-Comp_Blebbi',4438)};
    % NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M2','DC-Comp_CalA25',4438)};
    NotSaved = {NotSaved{:},Var2Data_wFluo('30-10-18','R80cst',0.6,'M2','DC-Comp_Ctrl',4438)};
    
    
    %%%%%%%%% billes internes %%%%%%%%%
    Var2Data_wFluo('28-08-18','R80cst',0.56,'M1','DC-Comp_Inside',4415,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('28-08-18','R80cst',0.56,'M2','DC-Comp_Inside',4415,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('17-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('17-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('24-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('24-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('25-09-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('25-09-18','R80cst',0.56,'M2','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('12-12-18','R80cst',0.56,'M1','DC-Comp_Inside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('12-12-18','R80cst',0.56,'M2','DC-Comp_NaseInside',4438,...
        RawdataFolder,MatfileFolder)
    
    Var2Data_wFluo('18-12-18','R80cst',0.56,'M2','DC-Comp_NaseInside',4438,...
        RawdataFolder,MatfileFolder)
    
    % billes externes
    
    Var2Data_wFluo('25-09-18','R80cst',0.56,'M1','DC-Comp_Outside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('25-09-18','R80cst',0.56,'M2','DC-Comp_Outside',4438,...
        RawdataFolder,MatfileFolder)
    
    Var2Data_wFluo('30-10-18','R80cst',0.56,'M1','DC-Comp_Outside',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('30-10-18','R80cst',0.56,'M2','DC-Comp_Outside',4438,...
        RawdataFolder,MatfileFolder)
    
    Var2Data_wFluo('21-01-19','R80cst',0.56,'M1','DC-Comp_Bsa',4438,...
        RawdataFolder,MatfileFolder)
    
    Var2Data_wFluo('23-01-19','R80cst',0.56,'M1','DC-Comp_Naked',4438,...
        RawdataFolder,MatfileFolder)
    Var2Data_wFluo('23-01-19','R80cst',0.56,'M1','DC-Comp_SerumJ1',4438,...
        RawdataFolder,MatfileFolder)
    
    
    Var2Data_wFluo('28-01-19','R80cst',0.56,'M1','DC-Comp_SerumJ6',4438,...
        RawdataFolder,MatfileFolder)
    
    
    % data from comp
    
    Data2Peaks({'M1','M2'},'DC-Comp_Inside',0,0.9,MatfileFolder,FigureFolder);
    Data2Peaks({'M1','M2'},'DC-Comp_Outside',1,0.9,MatfileFolder,FigureFolder);
    
    Data2Peaks({'M1','M2'},'DC-Comp_Naked',0,0.9,MatfileFolder,FigureFolder);
    
    Data2Peaks({'M1','M2'},'DC-Comp_SerumJ6',0,0.9,MatfileFolder,FigureFolder);
    
    
    % data from comp
    PLOTPS = 0;
    Peaks2Stats('DC-Comp_Inside',PLOTPS,MatfileFolder,FigureFolder);
    Peaks2Stats('DC-Comp_Outside',PLOTPS,MatfileFolder,FigureFolder);
    
    Peaks2Stats('DC-Comp_Naked');
    Peaks2Stats('DC-Comp_SerumJ6');
    
    FigureMaker(0,MatfileFolder,FigureFolder)
    
end