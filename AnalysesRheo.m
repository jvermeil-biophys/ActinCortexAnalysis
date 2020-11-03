close all;
clear all;
clc;

%% Paths & options
% h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';
RawdataFolder = 'E:\Raw';
MatfileFolder = 'C:\Users\PMMH\Desktop\DataValentin\MatFile';
DropboxDataFolder = '';
FigureFolder = 'C:\Users\PMMH\Desktop\DataValentin\Figures';


AUTO = 0;

V2DPLOTS = 0;

D2FRPLOTS = 0;

SS2MRPLOTS = 1;

%% Rheology analysis

  if 0   
   Res2Var_Rheology('Results','09-03-20',{'RS','S02','S1','S4'},'M1','Dicty_Rheo',1,'02-03-20_Depthograph100X',...
    AUTO,RawdataFolder,MatfileFolder)


 Var2Data_Rheology('09-03-20',{'RS','S02','S1','S4'},1.19,16.4,'M1','Dicty_Rheo',4503,...
    RawdataFolder,MatfileFolder,DropboxDataFolder,V2DPLOTS)


Data2StressStrain_Rheology('09-03-20',{'M1'},'Dicty_Rheo',{'RS','S02','S1','S4'},4503,D2FRPLOTS,MatfileFolder,FigureFolder)

StressStrain2Meca_Rheology('Dicty_Rheo',{'S02','S1','S4'},SS2MRPLOTS,MatfileFolder,FigureFolder)

  end