%% 03 02 18

names = {'1' '2' '3'};

for kn = 1:length(names)
filename = ['G:\Raw\EtalonnageZ\NotUsed\New_030218\K' names{kn} '_20n_Results.txt'];
stackname = ['G:\Raw\EtalonnageZ\NotUsed\New_030218\K' names{kn} '_20n.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\03-02-18_Depthograph_' names{kn}];
    [VarName1,Area,Mean,StdDev,Min,Max,X,Y,XM,YM,Major,Minor,Angle1,Slice] = importfileOld(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 17 05 18
names = {'1' '2' '3' '4'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\18.05.17_Zkimo\DCMed_7pcOpti_Kymo20nm_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\18.05.17_Zkimo\DCMed_7pcOpti_Kymo20nm_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\17-05-18_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

names = {'1' '2' '3' '4'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\18.05.17_Zkimo\DCMed_8pcOpti_Kymo20nm_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\18.05.17_Zkimo\DCMed_8pcOpti_Kymo20nm_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\17-05-18_Depthograph_' num2str(str2double(names{kn})+4)];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 12 07 18 

names = {'1' '2-1' '2-2' '3'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\18.07.11_NewKimo\DCMed_7-5pc_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\18.07.11_NewKimo\DCMed_7-5pc_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\12-07-18_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 18 07 18

names = {'1-1' '1-2' '2-1' '2-2'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\18.07.18\180718_Kimo_M1P' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\18.07.18\180718_Kimo_M1P' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\18-07-18_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 31 07 18

names = {'1' '2' '3' '4'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\18.07.31\180731_DCmed7-5pc_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\18.07.31\180731_DCmed7-5pc_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\31-07-18_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 28 08 18

names = {'1' '2' '3'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\18.08.28\28-08-18_Kimo' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\18.08.28\28-08-18_Kimo' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\28-08-18_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 30 10 18

names = {'1' '2' '3' '4' '5'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\18.10.30\Deptho-' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\18.10.30\Deptho-' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\30-10-18_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 12 12 18

names = {'1' '2-1' '2-2' '3' '4'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\18.12.12\Deptho_18-12-12_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\18.12.12\Deptho_18-12-12_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\12-12-18_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 15 03 19

names = {'1' '2' '3' '4' '6' '7' '8' '9-1' '9-2'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\19.03.15_Deptho\Deptho' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\19.03.15_Deptho\Deptho' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\15-03-19_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 26 06 19

names = {'11' '12' '13' '21' '22' '31' '32' '41' '42' '43'};

for kn = 1:length(names)
filename = ['D:\Data\Raw\19.06.26_Deptho\Deptho-' names{kn} '_res.txt'];
stackname = ['D:\Data\Raw\19.06.26_Deptho\Deptho-' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\26-06-19_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 15 10 19

names = {'1' '2' '3' '4' '5-1' '5-2'};

for kn = 1:length(names)
filename = ['E:\Raw\19.10.15\191015_Deptho-' names{kn} '_Results.txt'];
stackname = ['E:\Raw\19.10.15\191015_Deptho-' names{kn} '.tif'];
savename = ['E:\Raw\EtalonnageZ\MultiZCorrection\15-10-19_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 04 02 2020 63X

names = {'1-1' '1-2' '1-3' '2-1' '2-2' '3-1' '3-2' '3-3' '3-4'}
names = {'1-1'}

for kn = 1:length(names)
filename = ['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Deptho_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Deptho_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\04-02-20_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 10 02 2020 63X


names = {'1' '2' '3' '4' '5' '6-1' '6-2'}


for kn = 1:length(names)
filename = ['D:\Data\Raw\20.02.10\NewTestDeptho63X_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.02.10\NewTestDeptho63X_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\10-02-20_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 02 03 2020 63X

names = {'1-1' '1-2' '2-1' '2-2' '3-1' '3-2'}


for kn = 1:length(names)
filename = ['D:\Data\Raw\20.03.02\Deptho63X_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.03.02\Deptho63X_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\02-03-20_Depthograph63X_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 02 03 2020 100X

names = {'1-1' '1-2' '1-3' '2-1' '2-2' '3-1' '3-2' '3-3' '3-4'}


for kn = 1:length(names)
filename = ['E:\Raw\20.03.02\Deptho100X_' names{kn} '_Results.txt'];
stackname = ['E:\Raw\20.03.02\Deptho100X_' names{kn} '.tif'];
savename = ['E:\Raw\EtalonnageZ\MultiZCorrection\02-03-20_Depthograph100X_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 19 05 2020 M450

names = {'1-1' '1-2' '1-3' '1-4' '2-1' '2-2' '3-1' '3-2'}
names = {'1-1'}


for kn = 1:length(names)
filename = ['E:\Raw\20.05.19\19-05-2020_Deptho_M450_HL5_' names{kn} '_Results.txt'];
stackname = ['E:\Raw\20.05.19\19-05-2020_Deptho_M450_HL5_' names{kn} '.tif'];
savename = ['E:\Raw\EtalonnageZ\MultiZCorrection\19-05-20_DepthographM450_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 29 05 2020 My1

names = {'1-1' '1-2' '2-1' '2-2' '2-3'}


for kn = 1:length(names)
filename = ['D:\Data\Raw\20.05.29\29-05-20_DepthoMy1_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.05.29\29-05-20_DepthoMy1_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\29-05-20_DepthographMy1_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2905(filename);
[K,f,stepnm] = MakeDepthoMy1(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 03 07 20 M270 100X

names = {'1-1' '1-2' '1-3' '2-1' '2-2' '2-3' '2-4' '4-1' '4-2' '4-3'}


for kn = 1:length(names)
filename = ['D:\Data\Raw\20.07.03_BeadsTest\03-07-20_M270_Deptho_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.07.03_BeadsTest\03-07-20_M270_Deptho_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\03-07-20_DepthoM270_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all


%% 26 08 20 M270 100X

names = {'1-1' '1-2' '1-3' '2-1' '2-2' '2-3' '2-4'}


for kn = 1:length(names)
filename = ['D:\Data\Raw\20.08.26\Deptho' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.08.26\Deptho' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\26-08-20_DepthoM270_' names{kn}];
     [VarName1, Area, StdDev, Min, Max, XM, YM, BX, BY, Width, Height, Slice] = importfile260820(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 22-10-20 M450 100X

names = {'11' '12' '13' '21' '22' '23' '24'}


for kn = 1:length(names)
filename = ['D:\Data\Raw\20.10.22\Depthograph_M450-' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.10.22\Depthograph_M450-' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\22-10-20_DepthoM450_' names{kn}];
     [VarName1, Area, StdDev, Min, Max, XM, YM, BX, BY, Width, Height, Slice] = importfile_221020(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 17-11-20 M450 100X

names = {'1' '2' '3' '4'}


for kn = 1:length(names)
filename = ['D:\Data\Raw\20.11.17\Deptho_M450-' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.11.17\Deptho_M450-' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\17-11-20_DepthoM450_' names{kn}];
     [VarName1, Area, StdDev, Min, Max, XM, YM, BX, BY, Width, Height, Slice] = importfile_171120(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 16-12-20 M450 100X M1

names = {'1' '2' '3' '4' '5'}


for kn = 1:length(names)
filename = ['D:\Matlab Analysis\Data_Joseph\Raw\20.12.16_Deptho\Deptho_M1\20-12-16_M1_PetriA2_aSFLnoDrug_Deptho_' names{kn} '_Results.txt'];
stackname = ['D:\Matlab Analysis\Data_Joseph\Raw\20.12.16_Deptho\Deptho_M1\20-12-16_M1_PetriA2_aSFLnoDrug_Deptho_' names{kn} '.tif'];
savename = ['D:\Matlab Analysis\Data_Joseph\Raw\EtalonnageZ\MultiZCorrection\16-12-20_Deptho_M1_' names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 16-12-20 M450 100X M2

names = {'1' '2' '3' '4'}


for kn = 1:length(names)
filename = ['D:\Matlab Analysis\Data_Joseph\Raw\20.12.16_Deptho\Deptho_M2\20-12-16_M2_PetriA1_aSFLnoDrug_Deptho_' names{kn} '_Results.txt'];
stackname = ['D:\Matlab Analysis\Data_Joseph\Raw\20.12.16_Deptho\Deptho_M2\20-12-16_M2_PetriA1_aSFLnoDrug_Deptho_' names{kn} '.tif'];
savename = ['D:\Matlab Analysis\Data_Joseph\Raw\EtalonnageZ\MultiZCorrection\16-12-20_Deptho_M2_' names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 18-01-21 M450 100X M1

names = {'1' '2' '3' '4' '5'};
folderPath = 'D:\Matlab Analysis\Data_Joseph\Raw\21.01.18_Deptho\Deptho_M1';
fileName = '21-01-18_M1_PetriE1_aSFLdoxy_Deptho_';
saveFileName = '21-01-18_Deptho_M1_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\Matlab Analysis\Data_Joseph\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 18-01-21 M450 100X M2

names = {'1' '2' '3' '4' '5'};
folderPath = 'D:\Matlab Analysis\Data_Joseph\Raw\21.01.18_Deptho\Deptho_M2';
fileName = '21-01-18_M2_PetriE3_aSFLnodrug_Deptho_';
saveFileName = '21-01-18_Deptho_M2_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\Matlab Analysis\Data_Joseph\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 18-01-21 M450 100X M3

names = {'1' '2' '3' '4' '5'};
folderPath = 'D:\Matlab Analysis\Data_Joseph\Raw\21.01.18_Deptho\Deptho_M3';
fileName = '21-01-18_M3_PetriE4_aSFLdoxy_Deptho_';
saveFileName = '21-01-18_Deptho_M3_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\Matlab Analysis\Data_Joseph\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 21-01-21 M450 100X M1

names = {'1' '2' '3' '4' '5' '6'};
folderPath = 'D:\Matlab Analysis\Data_Joseph\Raw\21.01.21_Deptho\Deptho_M1';
fileName = '21-01-21_M1_PetriF4_aSFLnodrug_Deptho_';
saveFileName = '21-01-21_Deptho_M1_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\Matlab Analysis\Data_Joseph\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 21-01-21 M450 100X M2

names = {'1' '2' '3'};
folderPath = 'D:\Matlab Analysis\Data_Joseph\Raw\21.01.21_Deptho\Deptho_M2';
fileName = '21-01-21_M2_PetriF1_aSFLdoxy_Deptho_';
saveFileName = '21-01-21_Deptho_M2_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\Matlab Analysis\Data_Joseph\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 21-01-21 M450 100X M3

names = {'1' '2' '3' '4' '5'};
folderPath = 'D:\Matlab Analysis\Data_Joseph\Raw\21.01.21_Deptho\Deptho_M3';
fileName = '21-01-21_M3_PetriF3_aSFLnodrug_Deptho_';
saveFileName = '21-01-21_Deptho_M3_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\Matlab Analysis\Data_Joseph\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 21-02-10 M450 100X M1

names = {'1' '2' '3' '4' '5' '6'};
folderPath = 'D:\ActinCortex_Data\Raw\21.02.10_Deptho\M1';
fileName = '21.02.10_Deptho_M1_';
saveFileName = '21-02-10_Deptho_M1_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\ActinCortex_Data\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 21-02-15 M450 100X M1

names = {'1' '2' '3' '4' '5' '6'};
folderPath = 'D:\ActinCortex_Data\Raw\21.02.15_Deptho\M1';
fileName = '21.02.15_Deptho_M1_';
saveFileName = '21-02-15_Deptho_M1_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\ActinCortex_Data\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 21-02-15 M450 100X M2

names = {'1' '2' '3' '4' '5' '6'};
folderPath = 'D:\ActinCortex_Data\Raw\21.02.15_Deptho\M2';
fileName = '21.02.15_Deptho_M2_';
saveFileName = '21-02-15_Deptho_M2_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\ActinCortex_Data\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all
%% 21-02-15 M450 100X M3

names = {'1' '2' '3' '4' '5' '6'};
folderPath = 'D:\ActinCortex_Data\Raw\21.02.15_Deptho\M3';
fileName = '21.02.15_Deptho_M3_';
saveFileName = '21-02-15_Deptho_M3_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\ActinCortex_Data\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all

%% 21-03-05 M450 100X M3

names = {'1' '2' '3' '4' '5' '6'};
folderPath = 'D:\ActinCortex_Data\Raw\21.03.05_Deptho\M2';
fileName = '21.03.05_Deptho_M2_';
saveFileName = '21-03-05_Deptho_M2_';

for kn = 1:length(names)
filename = [folderPath filesep fileName names{kn} '_Results.txt'];
stackname = [folderPath filesep fileName names{kn} '.tif'];
savename = ['D:\ActinCortex_Data\Raw\EtalonnageZ\MultiZCorrection' filesep saveFileName names{kn}];
     [VarName1, Area, StdDev, XM, YM, Slice] = importfile_161220(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all


