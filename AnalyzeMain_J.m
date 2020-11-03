close all;
clear all;
clc;

%% reinitialisation du test

% GeneralFolder = 'C:\Users\BioMecaCell\Desktop\ConstantFieldAnalysis_Joseph\'; % a changer suivant le stockage des données
GeneralFolder = 'D:\Matlab Analysis\'; % a changer suivant le stockage des données


cd(GeneralFolder)
% addpath(genpath('C:\Users\BioMecaCell\Desktop\ConstantFieldAnalysis_Joseph\Code_Joseph'))
% addpath(genpath('C:\Users\BioMecaCell\Desktop\ConstantFieldAnalysis_Joseph\Code_Joseph\ZDeptho'))
% addpath(genpath('C:\Users\BioMecaCell\Desktop\ConstantFieldAnalysis_Joseph\Code_Joseph\Subfunction'))
addpath(genpath('D:\PMMH\Matlab Analysis\Code_Joseph'))
addpath(genpath('D:\PMMH\Matlab Analysis\Code_Joseph\ZDeptho'))
addpath(genpath('D:\PMMH\Matlab Analysis\Code_Joseph\Subfunction'))


% Bout de code dangereux copié du fichier exemple
% if 0
% cd('C:\Users\BioMecaCell\Desktop\ConstantFieldAnalysis_Joseph\Data_Joseph')
% 
% rmdir Figures s
% rmdir MatFiles s
% rmdir Raw\EtalonnageZ s
% 
% mkdir MatFiles
% mkdir Figures
% mkdir Raw\EtalonnageZ\Intermediate
% end

%% Paths

RawdataFolder = [GeneralFolder 'Data_Joseph\Raw']; % Raw data localisation
MatfileFolder = [GeneralFolder 'Data_Joseph\Matfiles']; % for saving data in .mat
SecondarySaveFolder = ''; % secondary spot for backup saving
FigureFolder = [GeneralFolder 'Data_Joseph\Figures']; % folder for saving figures

set(0, 'defaultfigureposition', get(0, 'Screensize'));

%% Settings

ExperimentalConditions3T3.Objectif = '100X'; % '63X'; % '100X'
ExperimentalConditions3T3.Scale = 1/15.8; % 1/9.7; % 15.8
ExperimentalConditions3T3.MagFieldMagnet = 0;
ExperimentalConditions3T3.MagFieldCorrectionCoeff = 1.15;
ExperimentalConditions3T3.OpticalDepthCorrectionCoeff = 1.33/1.52;
ExperimentalConditions3T3.TypeCellulaire = '3T3';
ExperimentalConditions3T3.TypeBilles = 'M450';
ExperimentalConditions3T3.DiametreBilles = 4503;
ExperimentalConditions3T3.SubstrateCoating = 'BSA';
ExperimentalConditions3T3.Drug = 'No';

ExperimentalConditionsDictys.Objectif = '100X'; % '63X'; % '100X'
ExperimentalConditionsDictys.Scale = 1/15.8; % 1/9.7; % 15.8
ExperimentalConditionsDictys.MagFieldMagnet = 0;
ExperimentalConditionsDictys.MagFieldCorrectionCoeff = 1.15;
ExperimentalConditionsDictys.OpticalDepthCorrectionCoeff = 1.33/1.52;
ExperimentalConditionsDictys.TypeCellulaire = 'dictys';
ExperimentalConditionsDictys.TypeBilles = 'M450';
ExperimentalConditionsDictys.DiametreBilles = 4503;
ExperimentalConditionsDictys.SubstrateCoating = 'BSA';
ExperimentalConditionsDictys.Drug = 'No';

ExperimentalConditionsFloat.Objectif = '63X'; % '63X'; % '100X'
ExperimentalConditionsFloat.Scale = 1/9.7; % 1/9.7; % 15.8
ExperimentalConditionsFloat.MagFieldMagnet = 0;
ExperimentalConditionsFloat.MagFieldCorrectionCoeff = 1.15;
ExperimentalConditionsFloat.OpticalDepthCorrectionCoeff = 1.33/1.0;
ExperimentalConditionsFloat.TypeCellulaire = 'dictys';
ExperimentalConditionsFloat.TypeBilles = 'M450';
ExperimentalConditionsFloat.DiametreBilles = 4503;
ExperimentalConditionsFloat.SubstrateCoating = 'BSA';
ExperimentalConditionsFloat.Drug = 'No';

%% Création des depthographs
if 0
    
DepthoMaker_J

Depthos2Deptho_J

end

%% Parameters

% pour r2v
AUTO = 0; % automatic analysis mode (skip when user input is needed)

% pour v2d
NotSaved = {}; % gather a list of non saved data (because of sorting criteria)

% pour d2p
PLOTDP = 1; % plot curve with cortex level

% pour p2s
PLOTPS = 1;
alignlvl = 0.5; % for cumulative probability curves



%% Analysis
if 0
%% R2V

% Example
    % Res2Var_wFluo_multiZ_J('dd-mm-yy', '5mT', 'M1','dicty'    ,1:3 ,'Results',3,'DD-MM-YY_Depthograph' ,AUTO,...
    %     RawdataFolder,MatfileFolder,SecondarySaveFolder)
 
    %(jour de manip, code manip, numero de la manip, nom de la condition  , puits a analyser ,' type du fichier .txt de resultas, nb image par loop (3 si pas d'image fluo faite), nom du depthograph a utiliser,
       % analyse automatique 0/1, RawdataFolder, MatfileFolder, SecondarySaveFolder)

% Historique
       
      Res2Var_wFluo_multiZ_J('07-08-20', 'Thickness10mT_3T3aSFLdoxi', 'M1','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,1, ...
          'Results',3,'20-02-20_Depthograph100x',AUTO,RawdataFolder,MatfileFolder,SecondarySaveFolder)
      
      Res2Var_wFluo_multiZ_J('07-08-20', 'Thickness10mT_3T3aSFLnodrugs', 'M2','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,1, ...
          'Results',3,'20-02-20_Depthograph100x',AUTO,RawdataFolder,MatfileFolder,SecondarySaveFolder)

      Res2Var_wFluo_multiZ_J('05-08-20', 'Thickness10mT_3T3aSFLdoxi', 'M1','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,1, ...
          'Results',3,'20-02-20_Depthograph100x',AUTO,RawdataFolder,MatfileFolder,SecondarySaveFolder)
      
      Res2Var_wFluo_multiZ_J('05-08-20', 'Thickness10mT_3T3aSFLnodrugs', 'M2','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,1, ...
          'Results',3,'20-02-20_Depthograph100x',AUTO,RawdataFolder,MatfileFolder,SecondarySaveFolder)

      Res2Var_wFluo_multiZ_J('04-08-20', 'Thickness10mT_3T3aSFLdoxi', 'M2','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,1, ...
          'Results',3,'20-02-20_Depthograph100x',AUTO,RawdataFolder,MatfileFolder,SecondarySaveFolder)
      
      Res2Var_wFluo_multiZ_J('04-08-20', 'Thickness10mT_3T3aSFLnodrugs', 'M1','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,1, ...
          'Results',3,'20-02-20_Depthograph100x',AUTO,RawdataFolder,MatfileFolder,SecondarySaveFolder)

%     Res2Var_wFluo_multiZ_J('20-02-20', '5mT', 'M1','dicty'    ,1:3 ,'Results',3,'20-02-20_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
%     Res2Var_wFluo_multiZ_J('21-02-20', '5mT', 'M1','test'    ,2 ,'Test',3,'20-02-20_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
%     Res2Var_wFluo_multiZ_J('10-03-20', '40mT', 'M1','3T3_BSA'    ,1 ,'Results',3,'20-02-20_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)
%     
%     Res2Var_wFluo_multiZ_J('10-03-20', '10mT', 'M1','3T3_PLL'    ,2 ,'Results',3,'20-02-20_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
%     Res2Var_wFluo_multiZ_J('25-05-20', 'XXmT', 'M1','Test_Code'    ,1:2 ,'Results',3,'20-02-20_Depthograph' ,AUTO,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Res2Var_wFluo_multiZ_J('30-06-20', '30mT_FallingChain37Mag', 'M1','30mT_FallingChain37Mag'    ,1 ,'Results',3,'15-06-20_Depthograph63x' ,AUTO,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)

    

    

%% V2D
    
    MakeClassification_J(GeneralFolder)

    Var2Data_wFluo_J('07-08-20','R40_3T3aSFLdoxicst','M1','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('07-08-20','R40_3T3aSFLnodrugscst','M2','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('05-08-20','R40_3T3aSFLdoxicst','M1','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('05-08-20','R40_3T3aSFLnodrugscst','M2','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('04-08-20','R40_3T3aSFLdoxicst','M2','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('04-08-20','R40_3T3aSFLnodrugscst','M1','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    %
    
    Var2Data_wFluo_J('07-08-20','Thickness10mT_3T3aSFLdoxi','M1','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('07-08-20','Thickness10mT_3T3aSFLnodrugs','M2','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('05-08-20','Thickness10mT_3T3aSFLdoxi','M1','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('05-08-20','Thickness10mT_3T3aSFLnodrugs','M2','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('04-08-20','Thickness10mT_3T3aSFLdoxi','M2','3T3aSFL_BSA_doxi',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
    Var2Data_wFluo_J('04-08-20','Thickness10mT_3T3aSFLnodrugs','M1','3T3aSFL_BSA_nodrugs',ExperimentalConditions3T3,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)
    
   
%     NotSaved = [NotSaved(:)',{Var2Data_wFluo_J('20-02-20','5mT',1.1,0,'M1','dicty',4415,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)}];
%     
%     NotSaved = [NotSaved(:)',{Var2Data_wFluo_J('21-02-20','5mT',1.1,0,'M1','test',4415,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)}];

%     NotSaved = [NotSaved(:)',{Var2Data_wFluo_J('10-03-20','40mT',1.19,35,'M1','3T3_BSA',4415,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)}];
%     
%     NotSaved = [NotSaved(:)',{Var2Data_wFluo_J('10-03-20','10mT',1.19,0,'M1','3T3_PLL',4415,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)}];

%  NotSaved = [NotSaved(:)',{Var2Data_wFluo_J('25-05-20','XXmT',1.19,0,'M1','Test_Code',4415,...
%         RawdataFolder,MatfileFolder,SecondarySaveFolder)}];
    
 NotSaved = [NotSaved(:)',{Var2Data_wFluo_J('30-06-20','30mT_FallingChain37Mag',1.19,0,'M1','30mT_FallingChain37Mag',4503,...
        RawdataFolder,MatfileFolder,SecondarySaveFolder)}];
   


%% D2P

    Data2Peaks_J({'04-08-20', '05-08-20', '07-08-20'},'3T3aSFL_BSA_doxi',PLOTDP,MatfileFolder,FigureFolder,SecondarySaveFolder);
    
    Data2Peaks_J({'04-08-20', '05-08-20', '07-08-20'},'3T3aSFL_BSA_nodrugs',PLOTDP,MatfileFolder,FigureFolder,SecondarySaveFolder);

%     Data2Peaks_J('dicty',PLOTDP,MatfileFolder,FigureFolder,SecondarySaveFolder);
%     
%     Data2Peaks_J('test',PLOTDP,MatfileFolder,FigureFolder,SecondarySaveFolder);
    
%     Data2Peaks_J('3T3_BSA',PLOTDP,MatfileFolder,FigureFolder,SecondarySaveFolder);
    
%     Data2Peaks_J('3T3_PLL',PLOTDP,MatfileFolder,FigureFolder,SecondarySaveFolder);
    
%     Data2Peaks_J('Test_Code',PLOTDP,MatfileFolder,FigureFolder,SecondarySaveFolder);
    
    Data2Peaks_J("30-06-20",'30mT_FallingChain37Mag',PLOTDP,MatfileFolder,FigureFolder,SecondarySaveFolder);


%% P2S
    
    Peaks2Stats_J({'04-08-20', '05-08-20', '07-08-20'},'3T3aSFL_BSA_doxi',PLOTPS,alignlvl,MatfileFolder,FigureFolder,SecondarySaveFolder);
    
    Peaks2Stats_J({'04-08-20', '05-08-20', '07-08-20'},'3T3aSFL_BSA_nodrugs',PLOTPS,alignlvl,MatfileFolder,FigureFolder,SecondarySaveFolder);

%     Peaks2Stats_J('dicty',PLOTPS,alignlvl,MatfileFolder,FigureFolder,SecondarySaveFolder);
%   
%     Peaks2Stats_J('test',PLOTPS,alignlvl,MatfileFolder,FigureFolder,SecondarySaveFolder);
    


%     Peaks2Stats_J('3T3_BSA',PLOTPS,alignlvl,MatfileFolder,FigureFolder,SecondarySaveFolder);
  
    Peaks2Stats_J('Test_Code',PLOTPS,alignlvl,MatfileFolder,FigureFolder,SecondarySaveFolder);
    
    Peaks2Stats_J('FallingChain37Mag',PLOTPS,alignlvl,MatfileFolder,FigureFolder,SecondarySaveFolder);
    

end

