% Goal: translate data contained in several .mat files in .csv documents to share with collaborators.

%% V2D tables

dataDir = 'D:\Matlab Analysis\Data_Joseph\MatFiles\V2D';
fileNames = {'V2D_04-08-20_M2_R40_3T3aSFL_BSA_doxy','V2D_05-08-20_M1_R40_3T3aSFL_BSA_doxy','V2D_07-08-20_M1_R40_3T3aSFL_BSA_doxy',...
    'V2D_04-08-20_M1_R40_3T3aSFL_BSA_nodrugs','V2D_05-08-20_M2_R40_3T3aSFL_BSA_nodrugs','V2D_07-08-20_M2_R40_3T3aSFL_BSA_nodrugs'};

saveDir = 'D:\Matlab Analysis\Data_Joseph\ExportFiles';

for n=1:length(fileNames)
    load([dataDir filesep fileNames{n}])
    N = length(MR);
    imagePerLoop = 111;
    for i=1:N
        cellName = MR{i}.name;
        M = length(MR{i}.time);
        CompNum = zeros(1,M);
        rampIndices = MR{1,i}.RampData{1,2};
        for j=1:M
            if(mod(j,imagePerLoop) > 8 && mod(j,imagePerLoop) < 104)
                CompNum(j) = 1 + floor(j/imagePerLoop);
            end
        end
        T = table(CompNum', MR{i}.time, MR{i}.B, MR{i}.F, MR{i}.dz, MR{i}.D2, MR{i}.D3,...
                'VariableNames',{'CompNum','T','B','F','dz','D2','D3'});
        writetable(T,[saveDir filesep cellName])
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
dataDir = 'D:\Matlab Analysis\Data_Joseph\MatFiles';
load([dataDir filesep 'MecaData3T3Table'])
saveDir = 'D:\Matlab Analysis\Data_Joseph\ExportFiles';
writetable(removevars(SubsetTable,{'RawDataTDFB'}),[saveDir filesep '3T3aSFL_MecaAnalysis'])