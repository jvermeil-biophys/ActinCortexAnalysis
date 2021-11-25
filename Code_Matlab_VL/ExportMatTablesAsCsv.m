% Goal: translate data contained in several .mat files in .csv documents to share with collaborators.

%% time Series Data from V2D

timeSeriesDataSourceDir = 'D:\Matlab Analysis\Data_Joseph\MatFiles\V2D';
timeSeriesDataSaveDir = 'D:\Matlab Analysis\Data_Joseph\ExportFiles';
% fileNames = ["V2D_04-08-20_M2_R40_3T3aSFL_BSA_doxy","V2D_05-08-20_M1_R40_3T3aSFL_BSA_doxy","V2D_07-08-20_M1_R40_3T3aSFL_BSA_doxy",...
%     "V2D_04-08-20_M1_R40_3T3aSFL_BSA_nodrugs","V2D_05-08-20_M2_R40_3T3aSFL_BSA_nodrugs","V2D_07-08-20_M2_R40_3T3aSFL_BSA_nodrugs"];
targetFiles = [];
listFiles = dir(timeSeriesDataSourceDir);
listFiles = listFiles(3:end);
nFiles = length(listFiles);
for iFile = 1:nFiles
    fileName = listFiles(iFile).name;
    if contains(fileName,'3T3') && contains(fileName,'.mat')
        targetFiles(end+1) = fileName;
    end
end

for iFile=1:length(targetFiles)
    load([timeSeriesDataDir filesep targetFiles(iFile)])
    nCells = length(MR);
    imagePerLoop = 111;
    for i=1:nCells
        cellName = MR{i}.name;
        M = length(MR{i}.time);
        CompNum = zeros(1,M);
        rampIndices = MR{1,i}.RampData{1,2};
        for j=1:M
            if(mod(j,imagePerLoop) > 8 && mod(j,imagePerLoop) < 104)
                CompNum(j) = 1 + floor(j/imagePerLoop);
            end
        end
        T = table(CompNum', MR{i}.time, MR{i}.B, MR{i}.F, MR{i}.dx, MR{i}.dy, MR{i}.dz, MR{i}.D2, MR{i}.D3,...
                'VariableNames',{'CompNum','T','B','F','dx','dy','dz','D2','D3'});
        writetable(T,[timeSeriesDataSaveDir filesep cellName])
    %     size(CompNo)
    %     size(MR{i}.time)
    %     size(MR{i}.B)
    %     size(MR{i}.F)
    %     size(MR{i}.dz)
    %     size(MR{i}.D2)
    %     size(MR{i}.D3)
    end
end

%% MecaData3T3Table
tableDataSourceDir = 'D:\Matlab Analysis\Data_Joseph\MatFiles';
tableDataSaveDir = 'D:\Matlab Analysis\Data_Joseph\ExportFiles';
load([tableDataSourceDir filesep 'MecaDataTable'])
threeTthreeTable = BigTable(contains(BigTable.ExpType, '3T3aSFL'),:);
writetable(removevars(threeTthreeTable,{'RawDataTDFB'}),[tableDataSaveDir filesep '3T3aSFL_MecaAnalysis'])