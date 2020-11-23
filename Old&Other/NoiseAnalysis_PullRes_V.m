function NoiseAnalysis_PullRes_V(names,datafolder)

AvgD2 = [];
StdD2 = [];
AvgD3 = [];
StdD3 = [];

for k = 1:length(names)
    load([datafolder filesep 'NoiseAnalysis_' names{k} '.mat'])
    AvgD2 = [AvgD2 cellfun(@median,{resultTableNanometers.D2})];
    StdD2 = [StdD2 [resultTableNanometers.D2std]];
    AvgD3 = [AvgD3 cellfun(@median,{resultTableNanometers.D3})];
    StdD3 = [StdD3 [resultTableNanometers.D3std]];
end


figure

end