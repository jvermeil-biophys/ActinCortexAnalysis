%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Permet de faire des depthographs individuels a partir de chaque
% vidéo en Z faite sur des billes collées. A lancer après l'analyse ImageJ.
% 
% A besoin de la stack d'image en .tif et du fichier de résultat du
% tracking de la bille en .txt
% 
% Attention a ce que le format du .txt soit compatible avec le fonction
% d'importation importfileDeptho.m sinon, générer une nouvelle fonction
% d'importation spécifique au format du .txt
% 
% Un bloc par jour de manip
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 08 07 20 % 63XAir
if 0
GeneralFolder = 'D:\PMMH\Matlab Code\ConstantFieldAnalysis_Joseph\'; % a changer suivant le stockage des données

datafolder = [GeneralFolder 'Data_Joseph\Raw'];
savefolder  = [GeneralFolder 'Data_Joseph\Raw\EtalonnageZ\Intermediate'];


names = {'1' '2' '3' '4'}; % liste des video en Z de billes dans le datafolder

for kn = 1:length(names)
filename = [datafolder filesep '20.07.08_Deptho\08-07-20_Deptho_P1_FallenBead_' names{kn} '_Results.txt'];
stackname = [datafolder filesep '20.07.08_Deptho\08-07-20_Deptho_P1_FallenBead_' names{kn} '.tif'];
savename = [savefolder filesep '08-07-20_Deptho_P1_FallenBead_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end
end

%% 15 06 20 % 100X
if 0
GeneralFolder = 'D:\PMMH\Matlab Code\ConstantFieldAnalysis_Joseph\'; % a changer suivant le stockage des données

datafolder = [GeneralFolder 'Data_Joseph\Raw'];
savefolder  = [GeneralFolder 'Data_Joseph\Raw\EtalonnageZ\Intermediate'];


names = {'1' '2' '3' '4' '5' '6' '7'}; % liste des video en Z de billes dans le datafolder

for kn = 1:length(names)
filename = [datafolder filesep '20.06.15_Deptho\100X_DictyMed+Optiprep40%\15-06-20_Deptho_100xHL5+40%Optiprep_' names{kn} '_Results.txt'];
stackname = [datafolder filesep '20.06.15_Deptho\100X_DictyMed+Optiprep40%\15-06-20_Deptho_100xHL5+40%Optiprep_' names{kn} '.tif'];
savename = [savefolder filesep '15-06-20_Deptho100X_DictyMed+Optiprep40%_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end
end
%% 15 06 20 % 63XAir
if 0
GeneralFolder = 'D:\PMMH\Matlab Code\ConstantFieldAnalysis_Joseph\'; % a changer suivant le stockage des données

datafolder = [GeneralFolder 'Data_Joseph\Raw'];
savefolder  = [GeneralFolder 'Data_Joseph\Raw\EtalonnageZ\Intermediate'];


names = {'1' '2' '3' '4' '5'}; % liste des video en Z de billes dans le datafolder

for kn = 1:length(names)
filename = [datafolder filesep '20.06.15_Deptho\63X_DictyMed+Optiprep40%\15-06-20_Deptho_63xHL5+40%Optiprep_' names{kn} '_Results.txt'];
stackname = [datafolder filesep '20.06.15_Deptho\63X_DictyMed+Optiprep40%\15-06-20_Deptho_63xHL5+40%Optiprep_' names{kn} '.tif'];
savename = [savefolder filesep '15-06-20_Deptho63X_DictyMed+Optiprep40%_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end
end
%% 20 02 20 %
if 0
GeneralFolder = 'C:\Users\BioMecaCell\Desktop\ConstantFieldAnalysis_Joseph\'; % a changer suivant le stockage des données

datafolder = [GeneralFolder 'Data_Joseph\Raw'];
savefolder  = [GeneralFolder 'Data_Joseph\Raw\EtalonnageZ\Intermediate'];


names = {'1' '2' '3' '3-1' '4' '5' '5-1'}; % liste des video en Z de billes dans le datafolder

for kn = 1:length(names)
filename = [datafolder filesep '20.02.20_Deptho\20-02-20_Deptho_DictyMedium_20nm_' names{kn} '_Results.txt'];
stackname = [datafolder filesep '20.02.20_Deptho\20-02-20_Deptho_DictyMedium_20nm_' names{kn} '.tif'];
savename = [savefolder filesep '20-02-20_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end
end
%% DD MM YY % date de la manip - EXEMPLE
if 0
    
GeneralFolder = 'C:\Users\BioMecaCell\Desktop\ConstantFieldAnalysis_Joseph\'; % a changer suivant le stockage des données


% le reste de la structure de dossier est a respecter

datafolder = [GeneralFolder 'Data_Joseph\Raw'];
savefolder  = [GeneralFolder 'Data_Joseph\Raw\EtalonnageZ\Intermediate'];


names = {'1' '2' '3' '4' '5' '6' '7' '8'}; % liste des video en Z de billes dans le datafolder

for kn = 1:length(names)
filename = [datafolder filesep 'YY.MM.DD_Deptho\Deptho_20nm_' names{kn} '_Results.txt'];
stackname = [datafolder filesep 'YY.MM.DD_Deptho\Deptho_20nm_' names{kn} '.tif'];
savename = [savefolder filesep 'DD-MM-YY_Depthograph_' names{kn}];
    [VarName1,Area,StdDev,XM,YM,Slice] = importfileDeptho(filename);
[K,f,stepnm] = MakeDeptho(StdDev,XM,YM,Slice,stackname,savename,20);
end

end
