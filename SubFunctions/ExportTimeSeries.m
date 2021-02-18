function ExportTimeSeries(MR, ExportDataFolder)

nCells = length(MR);
for i=1:nCells
    cellName = MR{i}.name;
    T = table(MR{i}.idxRamp_all, MR{i}.time, MR{i}.B, MR{i}.F, MR{i}.dx, MR{i}.dy, MR{i}.dz, MR{i}.D2, MR{i}.D3,...
            'VariableNames',{'idxCompression','T','B','F','dx','dy','dz','D2','D3'});
    writetable(T,[ExportDataFolder filesep cellName])
end