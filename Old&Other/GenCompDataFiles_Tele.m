RawdataFolder = 'E:\Raw';
MatfileFolder = 'C:\Users\PMMH\Desktop\DataValentin\MatFile';
DropboxDataFolder = '';
FigureFolder = 'C:\Users\PMMH\Desktop\DataValentin\Figures';

% DMSO
distDMSO = [];

distDMSO = [distDMSO Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4415,RawdataFolder,MatfileFolder)];

distDMSO = [distDMSO Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4415,RawdataFolder,MatfileFolder)];
  
distDMSO = [distDMSO Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4438,RawdataFolder,MatfileFolder)];  

distDMSO = [distDMSO Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4438,RawdataFolder,MatfileFolder)];
    
distDMSO = [distDMSO Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4438,RawdataFolder,MatfileFolder)];
    
distDMSO = [distDMSO Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4438,RawdataFolder,MatfileFolder)];

MeanDist = -mean(distDMSO,2)*1000;
ErrBot = MeanDist - (std(distDMSO,[],2)/sqrt(size(distDMSO,2)))*1000;
ErrTop = MeanDist + (std(distDMSO,[],2)/sqrt(size(distDMSO,2)))*1000;

save('C:\Users\PMMH\Dropbox\TheseValentin\Papiers\MatlabFig\MatFiles\CompData_DMSO','MeanDist','ErrBot','ErrTop')

clear MeanDist ErrBot ErrTop

% Bare Beads
distNaked = [];

distNaked = [distNaked Var2Data_wFluo_Comp('23-01-19','R80',0.56,'M1','DC-Comp_Naked',95,179,4438,RawdataFolder,MatfileFolder)];
 
MeanDist = -mean(distNaked,2)*1000;
ErrBot = MeanDist - (std(distNaked,[],2)/sqrt(size(distNaked,2)))*1000;
ErrTop = MeanDist + (std(distNaked,[],2)/sqrt(size(distNaked,2)))*1000;

save('C:\Users\PMMH\Dropbox\TheseValentin\Papiers\MatlabFig\MatFiles\CompData_BareBeads','MeanDist','ErrBot','ErrTop')

clear MeanDist ErrBot ErrTop

% Latranculin A
distLatA = [];

distLatA = [distLatA Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_LatA500',95,179,4438,RawdataFolder,MatfileFolder)];
   
MeanDist = -mean(distLatA,2)*1000;
ErrBot = MeanDist - (std(distLatA,[],2)/sqrt(size(distLatA,2)))*1000;
ErrTop = MeanDist + (std(distLatA,[],2)/sqrt(size(distLatA,2)))*1000;

save('C:\Users\PMMH\Dropbox\TheseValentin\Papiers\MatlabFig\MatFiles\CompData_LatA','MeanDist','ErrBot','ErrTop')

clear MeanDist ErrBot ErrTop

% Outside Beads
distOutside = [];

distOutside = [distOutside Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Outside',95,179,4438,RawdataFolder,MatfileFolder)];

distOutside = [distOutside Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Outside',95,179,4438,RawdataFolder,MatfileFolder)];
    
distOutside = [distOutside Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Outside',95,179,4438,RawdataFolder,MatfileFolder)];
    
distOutside = [distOutside Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Outside',95,179,4438,RawdataFolder,MatfileFolder)];


MeanDist = -mean(distOutside,2)*1000;
ErrBot = MeanDist - (std(distOutside,[],2)/sqrt(size(distOutside,2)))*1000;
ErrTop = MeanDist + (std(distOutside,[],2)/sqrt(size(distOutside,2)))*1000;

save('C:\Users\PMMH\Dropbox\TheseValentin\Papiers\MatlabFig\MatFiles\CompData_Outside','MeanDist','ErrBot','ErrTop')

clear MeanDist ErrBot ErrTop

% Inside Beads
distInside = [];

distInside = [distInside Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4415,RawdataFolder,MatfileFolder)];

distInside = [distInside Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4415,RawdataFolder,MatfileFolder)];

distInside = [distInside Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438,RawdataFolder,MatfileFolder)];

distInside = [distInside Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438,RawdataFolder,MatfileFolder)];

distInside = [distInside Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438,RawdataFolder,MatfileFolder)];

distInside = [distInside Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438,RawdataFolder,MatfileFolder)];

distInside = [distInside Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438,RawdataFolder,MatfileFolder)];

distInside = [distInside Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4438,RawdataFolder,MatfileFolder)];

distInside = [distInside Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4438,RawdataFolder,MatfileFolder)];

MeanDist = -mean(distInside,2)*1000;
ErrBot = MeanDist - (std(distInside,[],2)/sqrt(size(distInside,2)))*1000;
ErrTop = MeanDist + (std(distInside,[],2)/sqrt(size(distInside,2)))*1000;

save('C:\Users\PMMH\Dropbox\TheseValentin\Papiers\MatlabFig\MatFiles\CompData_Inside','MeanDist','ErrBot','ErrTop')

clear MeanDist ErrBot ErrTop