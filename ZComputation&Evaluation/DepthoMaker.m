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
filename = ['D:\Data\Raw\19.10.15\191015_Deptho-' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\19.10.15\191015_Deptho-' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\15-10-19_Depthograph_' names{kn}];
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
filename = ['D:\Data\Raw\20.03.02\Deptho100X_' names{kn} '_Results.txt'];
stackname = ['D:\Data\Raw\20.03.02\Deptho100X_' names{kn} '.tif'];
savename = ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\02-03-20_Depthograph100X_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho2(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

clear all



