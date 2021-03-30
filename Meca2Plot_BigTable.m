function Meca2Plot_BigTable(resfolder,figurefolder,Cond,Col,Sym,Lab,savename,SigDisplay,varargin)

% Meca2Plot_BigTable(resfolder,figurefolder,conditions,colors,labels)


%% INIT

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextInterpreter','Latex')
set(0,'DefaultLineMarkerSize',8);

here = pwd;

%% params and settings

% dossier de donn�es et sauvegarde des donn�es
path = [resfolder filesep 'D2M'];
datenow = datestr(now,'yy-mm-dd');


sfig = figurefolder;
sff = [sfig filesep datenow filesep 'Meca' filesep savename];
sffc = [sff filesep 'PerCell'];
mkdir(sff)


%% get data

load([resfolder filesep 'MecaDataTable'])


%% Get fluo data
pathTableFluo = "C:\Users\JosephVermeil\Desktop\ActinCortexAnalysis\DataAnalysis\FluoQuantification.txt";
% tableFluo = readtable(pathTableFluo, 'HeaderLines',1,'Format','%s%s%u'); 'DatetimeType'
tableFluo = readtable(pathTableFluo,'DatetimeType','text','TextType','string','Delimiter',','); % , 'HeaderLines',1
%tableFluo.uniqueID = strcat(tableFluo.('date'), '_', tableFluo.('id'));
tableFluoMeca = innerjoin(tableFluo, BigTable, 'LeftKeys', 'cellID', 'RightKeys', 'CellName');
tableFluoMeca(tableFluoMeca.Validated == 0,:)=[];
% [G, FluoResults] = findgroups(tableFluoMeca.uniqueID);
[G, cellIDtable, dateTable] = findgroups(tableFluoMeca(:,'cellID'),tableFluoMeca(:,'ExpDay'));
% meanValues = splitapply(@nanmean,tableFluoMeca.EChadwick,tableFluoMeca.H0Chadwick,tableFluoMeca.SurroundingThickness,tableFluoMeca.fluoLevel,G);
FluoResults = cellIDtable;
FluoResults.date = dateTable.Variables;
FluoResults.meanEChadwick = splitapply(@nanmean,tableFluoMeca.('EChadwick'),G);
FluoResults.meanH0Chadwick = splitapply(@nanmean,tableFluoMeca.('H0Chadwick'),G);
FluoResults.meanSurroundingThickness = splitapply(@nanmean,tableFluoMeca.('SurroundingThickness'),G);
FluoResults.meanFluoLevel = splitapply(@nanmean,tableFluoMeca.('meanFluoPeakAmplitude'),G);

%% filter excluded cells
pathTableExcluded = "C:\Users\JosephVermeil\Desktop\ActinCortexAnalysis\ExperimentalData\ExcludedCells.txt";
tableExcluded = readtable(pathTableExcluded,'TextType','string','Delimiter',',','ReadVariableNames',0); % , 'HeaderLines',1
tableExcluded.Properties.VariableNames{1}='cellID';
listExcluded = tableExcluded.cellID;
BigTable.manuallyExcluded = false(height(BigTable),1);
for i=1:length(listExcluded)
    BigTable.manuallyExcluded = BigTable.manuallyExcluded | strcmp(BigTable.CellName,listExcluded(i));
end

%%
if ~(size(Cond,1) == length(Col) && size(Cond,1) == length(Lab) && length(Lab) == length(Col))
    error('Numbers of conditions, colors, and lables don''t match !')
end


% varargin handling

excludevect = zeros(height(BigTable),1);

if ~isempty(varargin)
    
    if ismember('Exclude',varargin(1:2:end))
        ptrex = find(ismember(varargin(1:2:end),'Exclude'));
 
        for kex = 1:length(ptrex)
            condex = varargin{ptrex(kex)*2};
            
            excludevecttmp = strcmp(BigTable(:,condex{1}).Variables,condex{2});
            excludevect = excludevect | excludevecttmp;
        end
    end
end


%list of valid datas and info on non valid
VALIDptr = BigTable.Validated & ~BigTable.manuallyExcluded & ~excludevect;

ptrexptypes = find(ismember(Cond(1,:),'ExpType'))+1;

exptypesList = Cond(:,ptrexptypes);

ptrfitP= find(ismember(Cond(1,:),'FitParams'))+1;

fitPList = Cond(:,ptrfitP);

DoSuccessStats(resfolder,unique(exptypesList),unique(fitPList),excludevect)

% identifying changes between conditions

CondID = [];


%% Export a nice structure
SubsetBigTable = BigTable(contains(BigTable.ExpType, '3T3'),:);
[G, cellIDTable, expTypeTable] = findgroups(SubsetBigTable(:,'CellName') ,SubsetBigTable(:,'ExpType'));
% meanValues = splitapply(@nanmean,tableFluoMeca.EChadwick,tableFluoMeca.H0Chadwick,tableFluoMeca.SurroundingThickness,tableFluoMeca.fluoLevel,G);
ExportSummaryTable = cellIDTable;
expTypeTable.ExpType(strcmp(expTypeTable.ExpType, '3T3aSFL_BSA_nodrugs')) = '3T3aSFL_nonAdhered_nodrug';
expTypeTable.ExpType(strcmp(expTypeTable.ExpType, '3T3aSFL_BSA_doxy')) = '3T3aSFL_nonAdhered_doxycyclin';
expTypeTable.ExpType(strcmp(expTypeTable.ExpType, '3T3aSFL_nodrug')) = '3T3aSFL_Adhered20umFibro_nodrug';
expTypeTable.ExpType(strcmp(expTypeTable.ExpType, '3T3aSFL_doxy')) = '3T3aSFL_Adhered20umFibro_doxycyclin';
ExportSummaryTable.ExperimentType = expTypeTable.Variables;
ExportSummaryTable.meanEfitChadwick = splitapply(@nanmean,SubsetBigTable.('EChadwick'),G);
ExportSummaryTable.meanH0fitChadwick = splitapply(@nanmean,SubsetBigTable.('H0Chadwick'),G);
ExportSummaryTable.meanH0surrounding = splitapply(@nanmean,SubsetBigTable.('SurroundingThickness'),G);
ExpTypesList = unique(ExportSummaryTable.ExperimentType,'rows');
for i = 1:length(ExpTypesList)
    %tempStruct = struct([]);
    tempTable = ExportSummaryTable(strcmp(ExportSummaryTable.ExperimentType, ExpTypesList(i)),:);
    %tempStruct = struct('ExpType', ExpTypesList(i), 'CellName', tempTable.CellName, 'meanEfitChadwick', tempTable.meanEfitChadwick, 'meanH0fitChadwick', tempTable.meanH0fitChadwick, 'meanH0surrounding', tempTable.meanH0surrounding);
    ExportSummaryStructure(i) = struct('ExpType', ExpTypesList(i), 'CellName', tempTable.CellName, 'meanEfitChadwick', tempTable.meanEfitChadwick, 'meanH0fitChadwick', tempTable.meanH0fitChadwick, 'meanH0surrounding', tempTable.meanH0surrounding);
end
save(['C:\Users\JosephVermeil\Desktop\ActinCortexAnalysis\DataAnalysis' filesep 'ExportDataForJulien'],'ExportSummaryTable', 'ExportSummaryStructure');


%%
for k =2:2:size(Cond,2)
    
    if length(unique(Cond(:,k))) >1
        
        CondID = strcat(CondID,'-',Cond(:,k));
        
    end
    
end

% per cell ?
answer = questdlg('Plot and H0 per cell ? (creates many plots)', ...
    'Individual cell plotting', ...
    'Yes','No','No');

% retrieving data for different experimental conditions
ncond = size(Cond,1);

for kcond = 1:ncond
    
    CondSpecif = Cond(kcond,:);
    
    ntag = 0;
    
    for ktag = 1:2:length(CondSpecif)
        
        ntag = ntag +1 ;
        MatPtr(:,ntag) = ismember(BigTable(:,CondSpecif(ktag)).Variables,CondSpecif(ktag+1));
        
    end
    
    ptrcond = all(MatPtr,2) & VALIDptr;
    ptrfirst = BigTable(:,'CompNum').Variables == 1;
    
    E0chad{kcond} = BigTable(ptrcond,'EChadwick').Variables/1000;
    E0chadFirst{kcond} = BigTable(ptrcond&ptrfirst,'EChadwick').Variables/1000;
    EchadCI{kcond} = BigTable(ptrcond,'CiEChadwick').Variables/1000;
    
    Hyst{kcond} = BigTable(ptrcond,'Hysteresis').Variables;
    HystFirst{kcond} = BigTable(ptrcond&ptrfirst,'Hysteresis').Variables;
    
    H0EXP{kcond} = BigTable(ptrcond,'InitialThickness').Variables;
    H0EXPFirst{kcond} = BigTable(ptrcond&ptrfirst,'InitialThickness').Variables;
    H0SUR{kcond} = BigTable(ptrcond,'SurroundingThickness').Variables;
    H0SURFirst{kcond} = BigTable(ptrcond&ptrfirst,'SurroundingThickness').Variables;
    H0BEF{kcond} = BigTable(ptrcond,'PreviousThickness').Variables;
    H0FIT{kcond} = BigTable(ptrcond,'H0Chadwick').Variables;
    H0FITFirst{kcond} = BigTable(ptrcond&ptrfirst,'H0Chadwick').Variables;
    
    PatternRadius = 20e-6;
    Eeffectif{kcond} = (E0chad{kcond}.*(H0SUR{kcond}*1e-9).^2)/PatternRadius;
    
    CellIdList{kcond} = BigTable(ptrcond,'CellName').Variables;
    UniqueCellList{kcond} = unique(CellIdList{kcond});
    
    CompNums{kcond} = BigTable(ptrcond,'CompNum').Variables;
    CompTimes{kcond} = BigTable(ptrcond,'CompTime').Variables;
    
end

FullCellList = unique(vertcat(UniqueCellList{:}));

nCell = length(FullCellList);

%% prep figs

i_figure = 1;

if ncond <=2 % figures for 1 or less datasets
    
    figure(i_figure)
    i_FigLinEchadDistrib = i_figure;
    i_figure = i_figure + 1;
    hold on
    title('Linlin distrib + median')
    L1 = legend;
    xlabel('Echad (kPa)')
    ylabel('Count')
    
    
    figure(i_figure)
    i_FigLogEchadDistrib = i_figure;
    i_figure = i_figure + 1;
    hold on
    title('Loglog distrib + weighted avg')
    L2 = legend;
    xlabel('Elastic modulus (kPa)')
    ylabel('Count')
    
    
    
end

figure(i_figure)
i_FigEvH0 = i_figure;
i_figure = i_figure + 1;
hold on
title('E v H0')
L_EvH0 = legend;
xlabel('Thickness before comp(nm)')
ylabel('E chad (kPa)')

figure(i_figure)
i_FigEvHSUR = i_figure;
i_figure = i_figure + 1;
hold on
title('E v HSUR')
L_EvHSUR = legend;
xlabel('Thickness between compression (nm)')
ylabel('E chad (kPa)')

figure(i_figure)
i_FigHFITvHSUR = i_figure;
i_figure = i_figure + 1;
hold on
title('H0FIT v H0SUR')
L_HFITvHSUR = legend;
xlabel('Thickness fitted by Chadwick (nm)')
ylabel('Thickness between compression (nm)')

figure(i_figure)
i_FigEchadDistrib = i_figure;
i_figure = i_figure + 1;
hold on
title('Echad distrib')

figure(i_figure)
i_FigEchadPerCellDistrib = i_figure;
i_figure = i_figure + 1;
hold on
title('Echad cell distrib')

figure(i_figure)
i_FigEchadStdPerCellDistrib = i_figure;
i_figure = i_figure + 1;
hold on
title('Echad variance cell distrib')

figure(i_figure)
i_FigEchadPerCellAvgVSStd = i_figure;
i_figure = i_figure + 1;
hold on
title('Echad average vs standard deviation for each cell')

figure(i_figure)
i_FigH0Distrib = i_figure;
i_figure = i_figure + 1;
hold on
title('H0 distrib')

figure(i_figure)
i_FigSurroundingThicknessDistrib = i_figure;
i_figure = i_figure + 1;
hold on
title('Surrounding Thickness distrib')

figure(i_figure)
i_FigH0SURPerCellDistrib = i_figure;
i_figure = i_figure + 1;
hold on
title('Surrounding Thickness per Cell distrib')

figure(i_figure)
i_FigFittedThicknessDistrib = i_figure;
i_figure = i_figure + 1;
hold on
title('Fitted Thickness distrib')

figure(i_figure)
i_FigH0FITPerCellDistrib = i_figure;
i_figure = i_figure + 1;
hold on
title('Fitted Thickness per Cell  distrib')

figure(i_figure)
i_FigEandH0vsTime = i_figure;
i_figure = i_figure + 1;
subplot(211)
hold on
title('Echad & H0 vs Time')

figure(i_figure)
i_FigHysteresis = i_figure;
i_figure = i_figure + 1;
hold on
title('Hysteresis')

figure(i_figure)
i_FigFluo = i_figure;
i_figure = i_figure + 1;
hold on
title('Fluo')

figure(i_figure)
i_FigEeff = i_figure;
i_figure = i_figure + 1;
hold on
title('K_AFM')

figure(i_figure)
i_FigEeffPerCell = i_figure;
i_figure = i_figure + 1;
hold on
title('K_AFM per cell')



ymax1 = 0;
ymax2 = 0;

%% Plots


for kcond = 1:ncond
    
    
    fprintf(['\n\nFor condition : ' Lab{kcond} '\n'])
    
    fprintf(['\nNumber of compressions : ' num2str(length(E0chad{kcond})) '.\nNumber of cells : ' num2str(size(UniqueCellList{kcond},1)) '.\n'])
    
    
    
    if ncond <=2
        
        %%  distribution des Echad + median
        
        [counts,edges] = histcounts(E0chad{kcond},'binwidth',2);
        
        figure(i_FigLinEchadDistrib)
        hold on
        histogram(E0chad{kcond},edges,'facecolor',Col{kcond},'facealpha',0.5)
        plot([median(E0chad{kcond}) median(E0chad{kcond})],[0 1000],'--','linewidth',1.7,'color',Col{kcond},'handlevisibility','off')
        ax = gca;
        
        ymax1tmp = max(counts)*1.05;
        
        ymax1 = max([ymax1 ymax1tmp]);
        
        ylim([0 ymax1])
        
        %% Distrib Echad avec moyenne pond�r�e en log
        
        
        EchadWeight = (E0chad{kcond}./EchadCI{kcond}).^2;
        
        EchadWeightedMean = 10^(sum(log10(E0chad{kcond}).*EchadWeight)/sum(EchadWeight));
        
        EchadWeightedStd = sqrt(sum(((E0chad{kcond}-EchadWeightedMean).^2.*EchadWeight))/sum(EchadWeight));
        
        
        [counts,edges] = histcounts(log10(E0chad{kcond}));
        
        
        figure(i_FigLogEchadDistrib)
        hold on
        histogram(E0chad{kcond},10.^edges,'facecolor',Col{kcond},'facealpha',0.5)
        
        plot([EchadWeightedMean EchadWeightedMean],[0 1000],'--','linewidth',1.7,'color',Col{kcond},'handlevisibility','off')
        ax = gca;
        ax.LineWidth = 1.5;
        ax.XScale = 'log';
        ax.XTick = [1 10 100];
        box on
        
        
        xlim([0.5 200])
        
        
        ymax3tmp = max(counts)*1.05;
        
        ymax2 = max([ymax2 ymax3tmp]);
        
        ylim([0 ymax2])
        
        %% legends & other additions
        
        L1.String{end} = [Lab{kcond}  ' - ' num2str(median(E0chad{kcond})) 'kPa'];
        
        L2.String{end} = [Lab{kcond}  ' - ' num2str(EchadWeightedMean) '+- ' num2str(EchadWeightedStd) ' kPa'];
        
    end
    
    
    %% Module vs Epaisseur depart
    
    figure(i_FigEvH0)
    hold on
    plot(H0EXP{kcond},E0chad{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
    
    if ncond == 2
        ratioEH{kcond} = (E0chad{kcond}./H0EXP{kcond})/median(E0chad{kcond});
    end
    
    ax = gca;
    ax.YScale = 'log';
    ax.XTickLabelRotation = 45;
    
    L_EvH0.String{end} = Lab{kcond};
    
    %% Module vs Epaisseur surrounding
    
    figure(i_FigEvHSUR)
    hold on
    plot(H0SUR{kcond},E0chad{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
    
    if ncond == 2
        ratioEH{kcond} = (E0chad{kcond}./H0SUR{kcond})/median(E0chad{kcond});
    end
    
    ax = gca;
    ax.YScale = 'log';
    ax.XTickLabelRotation = 45;
    
    L_EvHSUR.String{end} = Lab{kcond};
    
    %% Epaisseur fit vs Epaisseur surrounding
    
    figure(i_FigHFITvHSUR)
    hold on
    plot(H0FIT{kcond},H0SUR{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
    
    if ncond == 2
        ratioEH{kcond} = (H0SUR{kcond}./H0FIT{kcond})/median(H0SUR{kcond});
    end
    
    ax = gca;
    %ax.YScale = 'log';
    ax.XTickLabelRotation = 45;
    
    L_HFITvHSUR.String{end} = Lab{kcond};
    
    
end


%% regrouping per cells and conditions
CellVsCond_Echad = cell(ncond,nCell);
CellVsCond_H0EXP = cell(ncond,nCell);
CellVsCond_H0SUR = cell(ncond,nCell);
CellVsCond_H0FIT = cell(ncond,nCell);
CellVsCond_Eeffectif = cell(ncond,nCell);
CellVsCond_ID = cell(ncond,nCell);

for kc = 1:nCell
    for kcond = 1:ncond
        
        ptrcell = strcmp(CellIdList{kcond}, FullCellList(kc));
        
        CellVsCond_Echad{kcond,kc} = E0chad{kcond}(ptrcell);
        CellVsCond_H0EXP{kcond,kc} = H0EXP{kcond}(ptrcell);
        CellVsCond_H0SUR{kcond,kc} = H0SUR{kcond}(ptrcell);
        CellVsCond_H0FIT{kcond,kc} = H0FIT{kcond}(ptrcell);
        CellVsCond_Eeffectif{kcond,kc} = Eeffectif{kcond}(ptrcell);
        CellVsCond_ID{kcond,kc} = strcat(FullCellList(kc), CondID(kcond));
        
    end
    
    
end

%% Beeswarm plot

% all data Echad
figure(i_FigEchadDistrib)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
ax.YScale = 'log';
%set(gca, 'YScale', 'log')
box on

% all data
plotSpread_V(E0chad,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(E0chad)
    plot([ii-0.45 ii+0.45],[nanmedian(E0chad{ii})  nanmedian(E0chad{ii})],'k--','linewidth',1.3)
end

% % first comp of each cell in another color
% set(0,'DefaultLineMarkerSize',10);
% plotSpread_V(E0chadFirst,'distributionMarkers',Sym,'distributionColors','y','xNames',Lab,'spreadWidth',1)
% set(0,'DefaultLineMarkerSize',8);

ylabel('Echad (kPa)')

if SigDisplay
    distsize = prctile(vertcat(E0chad{:}),98);
    plotheight = distsize;
    for ii = 1:length(E0chad)-1
        for jj = ii+1:length(E0chad)
            sigtxt = DoParamStats(E0chad{ii},E0chad{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight * 10^(0.03), sigtxt,'HorizontalAlignment','center','fontsize',13) %+0.005*distsize
                plotheight = plotheight * 10^(0.2);
            else
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight * 10^(0.09), sigtxt,'HorizontalAlignment','center','fontsize',13) %+0.005*distsize
                plotheight = plotheight * 10^(0.2);
            end
        end
    end
    ylim([0 plotheight])
end
% 
% if SigDisplay
%     distsize = prctile(vertcat(E0chad{:}),98);
%     plotheight = distsize;
%     for ii = 1:length(E0chad)-1
%         for jj = ii+1:length(E0chad)
%             sigtxt = DoParamStats(E0chad{ii},E0chad{jj});
%             if ~strcmp(sigtxt,'NS')
%                 plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
%                 text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
%                 plotheight = plotheight + 0.05*distsize;
%             end
%         end
%     end
%     ylim([0 plotheight])
% end

ax = gca;
ax.XTickLabelRotation = 45;








% avg per cell Echad

figure(i_FigEchadPerCellDistrib)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
ax.YScale = 'log';
box on

for kcond = 1:ncond
    
    E0chadCell{kcond} = [];
    
    for kc = 1:nCell
        
        E0chadCell{kcond} = [E0chadCell{kcond} nanmean(CellVsCond_Echad{kcond,kc})];
        
    end
end

plotSpread_V(E0chadCell,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(E0chadCell)
    plot([ii-0.45 ii+0.45],[nanmedian(E0chadCell{ii})  nanmedian(E0chadCell{ii})],'k--','linewidth',1.3)
end
ylabel('Echad Average par cell (kPa)')

if SigDisplay
    distsize = prctile(horzcat(E0chadCell{:}),98);
    plotheight = distsize;
    for ii = 1:length(E0chadCell)-1
        for jj = ii+1:length(E0chadCell)
            sigtxt = DoParamStats(E0chadCell{ii},E0chadCell{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight * 10^(0.03), sigtxt,'HorizontalAlignment','center','fontsize',13) %+0.005*distsize
                plotheight = plotheight * 10^(0.2);
            else
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight * 10^(0.09), sigtxt,'HorizontalAlignment','center','fontsize',13) %+0.005*distsize
                plotheight = plotheight * 10^(0.2);
            end
        end
    end
    ylim([0 plotheight])
end
% 
% if SigDisplay
%     distsize = prctile(horzcat(E0chadCell{:}),98);
%     plotheight = distsize;
%     for ii = 1:length(E0chadCell)-1
%         for jj = ii+1:length(E0chadCell)
%             sigtxt = DoParamStats(E0chadCell{ii},E0chadCell{jj});
%             if ~strcmp(sigtxt,'NS')
%                 plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
%                 text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
%                 plotheight = plotheight + 0.05*distsize;
%             end
%         end
%     end
%     ylim([0 plotheight])
% end


ax = gca;
ax.XTickLabelRotation = 45;







% var per cell Echad

figure(i_FigEchadStdPerCellDistrib)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
ax.YScale = 'log';
box on

for kcond = 1:ncond
    
    E0chadStdCell{kcond} = [];
    
    for kc = 1:nCell
        
        E0chadStdCell{kcond} = [E0chadStdCell{kcond} nanstd(CellVsCond_Echad{kcond,kc})];
        
    end
end

plotSpread_V(E0chadStdCell,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(E0chadStdCell)
    plot([ii-0.45 ii+0.45],[nanmedian(E0chadStdCell{ii})  nanmedian(E0chadStdCell{ii})],'k--','linewidth',1.3)
end
ylabel('Echad Variance par cell (kPa)')


if SigDisplay
    distsize = prctile(horzcat(E0chadStdCell{:}),98);
    plotheight = distsize;
    for ii = 1:length(E0chadStdCell)-1
        for jj = ii+1:length(E0chadStdCell)
            sigtxt = DoParamStats(E0chadStdCell{ii},E0chadStdCell{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                plotheight = plotheight + 0.05*distsize;
            end
        end
    end
    ylim([0 plotheight])
end


ax = gca;
ax.XTickLabelRotation = 45;





% avg vs var per cell Echad

figure(i_FigEchadPerCellAvgVSStd)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

for kcond = 1:ncond
    
    E0chadCell{kcond} = [];
    E0chadStdCell{kcond} = [];
    
    for kc = 1:nCell
        
        E0chadCell{kcond} = [E0chadCell{kcond} nanmean(CellVsCond_Echad{kcond,kc})];
        E0chadStdCell{kcond} = [E0chadStdCell{kcond} nanstd(CellVsCond_Echad{kcond,kc})];
    end
    
    plot(E0chadCell{kcond},E0chadStdCell{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
end


xlabel('Echad Average for each cell (kPa)')
ylabel('Echad Standard Deviation for each cell (kPa)')

ax = gca;
ax.XTickLabelRotation = 45;



% H0

H0Beeswarm = H0EXP;
H0BeeswarmFirst = H0EXPFirst;

figure(i_FigH0Distrib)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(H0Beeswarm,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(H0Beeswarm)
    plot([ii-0.45 ii+0.45],[nanmedian(H0Beeswarm{ii}) nanmedian(H0Beeswarm{ii})],'k--','linewidth',1.3)
end

% % first comp of each cell in another color
% set(0,'DefaultLineMarkerSize',10);
% plotSpread_V(H0BeeswarmFirst,'distributionMarkers',Sym,'distributionColors','r','xNames',Lab,'spreadWidth',1)
% set(0,'DefaultLineMarkerSize',8);

ylabel('H0 init ramp')


if SigDisplay
    distsize = prctile(vertcat(H0Beeswarm{:}),95);
    plotheight = distsize;
    for ii = 1:length(E0chad)-1
        for jj = ii+1:length(H0Beeswarm)
            sigtxt = DoParamStats(H0Beeswarm{ii},H0Beeswarm{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                plotheight = plotheight+ 0.05*distsize;
            end
        end
    end
    ylim([0 plotheight]) 
end

ax = gca;
ax.XTickLabelRotation = 45;






% Surrounding Thickness

H0Beeswarm = H0SUR;
H0BeeswarmFirst = H0SURFirst;

figure(i_FigSurroundingThicknessDistrib)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(H0Beeswarm,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(H0Beeswarm)
    plot([ii-0.45 ii+0.45],[nanmedian(H0Beeswarm{ii}) nanmedian(H0Beeswarm{ii})],'k--','linewidth',1.3)
end

% % first comp of each cell in another color
% set(0,'DefaultLineMarkerSize',10);
% plotSpread_V(H0BeeswarmFirst,'distributionMarkers',Sym,'distributionColors','r','xNames',Lab,'spreadWidth',1)
% set(0,'DefaultLineMarkerSize',8);

ylabel('Median thickness around this ramp (nm)')


if SigDisplay
    distsize = prctile(vertcat(H0Beeswarm{:}),95);
    plotheight = distsize;
    for ii = 1:length(E0chad)-1
        for jj = ii+1:length(H0Beeswarm)
            sigtxt = DoParamStats(H0Beeswarm{ii},H0Beeswarm{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                plotheight = plotheight+ 0.05*distsize;
            end
        end
    end
    ylim([0 plotheight])
end

ax = gca;
ax.XTickLabelRotation = 45;


% avg per cell H0 SUR

figure(i_FigH0SURPerCellDistrib)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

for kcond = 1:ncond
    
    H0SURCell{kcond} = [];
    
    for kc = 1:nCell
        
        H0SURCell{kcond} = [H0SURCell{kcond} nanmean(CellVsCond_H0SUR{kcond,kc})];
        
    end
end

plotSpread_V(H0SURCell,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(H0SURCell)
    plot([ii-0.45 ii+0.45],[nanmedian(H0SURCell{ii})  nanmedian(H0SURCell{ii})],'k--','linewidth',1.3)
end
ylabel('H0 surrounding per cell (nm)')

if SigDisplay
    distsize = prctile(horzcat(H0SURCell{:}),98);
    plotheight = distsize;
    for ii = 1:length(H0SURCell)-1
        for jj = ii+1:length(H0SURCell)
            sigtxt = DoParamStats(H0SURCell{ii},H0SURCell{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                plotheight = plotheight + 0.05*distsize;
            else
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                plotheight = plotheight + 0.05*distsize;
            end
        end
    end
    ylim([0 plotheight])
end

ax = gca;
ax.XTickLabelRotation = 45;




% Fitted Thickness

H0Beeswarm = H0FIT;
H0BeeswarmFirst = H0FITFirst;

figure(i_FigFittedThicknessDistrib)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(H0Beeswarm,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(H0Beeswarm)
    plot([ii-0.45 ii+0.45],[nanmedian(H0Beeswarm{ii}) nanmedian(H0Beeswarm{ii})],'k--','linewidth',1.3)
end

% % first comp of each cell in another color
% set(0,'DefaultLineMarkerSize',10);
% plotSpread_V(H0BeeswarmFirst,'distributionMarkers',Sym,'distributionColors','r','xNames',Lab,'spreadWidth',1)
% set(0,'DefaultLineMarkerSize',8);

ylabel('Fitted thickness for this ramp (nm)')

if SigDisplay
    distsize = prctile(vertcat(H0Beeswarm{:}),95);
    plotheight = distsize;
    for ii = 1:length(E0chad)-1
        for jj = ii+1:length(H0Beeswarm)
            sigtxt = DoParamStats(H0Beeswarm{ii},H0Beeswarm{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                plotheight = plotheight+ 0.05*distsize;
            end
        end
    end
    ylim([0 plotheight])
end
ax = gca;
ax.XTickLabelRotation = 45;




% avg per cell H0 FIT

figure(i_FigH0FITPerCellDistrib)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

for kcond = 1:ncond
    H0FITCell{kcond} = [];
    for kc = 1:nCell
        H0FITCell{kcond} = [H0FITCell{kcond} nanmean(CellVsCond_H0FIT{kcond,kc})];
    end
end

plotSpread_V(H0FITCell,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(H0FITCell)
    plot([ii-0.45 ii+0.45],[nanmedian(H0FITCell{ii})  nanmedian(H0FITCell{ii})],'k--','linewidth',1.3)
end
ylabel('H0 fitted per cell (nm)')

if SigDisplay
    distsize = prctile(horzcat(H0FITCell{:}),98);
    plotheight = distsize;
    for ii = 1:length(H0FITCell)-1
        for jj = ii+1:length(H0FITCell)
            sigtxt = DoParamStats(H0FITCell{ii},H0FITCell{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                plotheight = plotheight + 0.05*distsize;
            else
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                plotheight = plotheight + 0.05*distsize;
            end
        end
    end
    ylim([0 plotheight])
end


ax = gca;
ax.XTickLabelRotation = 45;


% Eeffectif

figure(i_FigEeff)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
ax.YScale = 'log';
box on

plotSpread_V(Eeffectif,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(Eeffectif)
    plot([ii-0.45 ii+0.45],[nanmedian(Eeffectif{ii}) nanmedian(Eeffectif{ii})],'k--','linewidth',1.3)
end

% % first comp of each cell in another color
% set(0,'DefaultLineMarkerSize',10);
% plotSpread_V(H0BeeswarmFirst,'distributionMarkers',Sym,'distributionColors','r','xNames',Lab,'spreadWidth',1)
% set(0,'DefaultLineMarkerSize',8);

ylabel('K_AFM (N/m) [ = E*H^2/R]')

if SigDisplay
    distsize = prctile(vertcat(Eeffectif{:}),98);
    plotheight = distsize;
    for ii = 1:length(Eeffectif)-1
        for jj = ii+1:length(Eeffectif)
            sigtxt = DoParamStats(Eeffectif{ii},Eeffectif{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight * 10^(0.03), sigtxt,'HorizontalAlignment','center','fontsize',13) %+0.005*distsize
                plotheight = plotheight * 10^(0.2);
            else
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight * 10^(0.09), sigtxt,'HorizontalAlignment','center','fontsize',13) %+0.005*distsize
                plotheight = plotheight * 10^(0.2);
            end
        end
    end
    ylim([0 plotheight])
end

% if SigDisplay
%     distsize = prctile(vertcat(Eeffectif{:}),95);
%     plotheight = distsize;
%     for ii = 1:length(Eeffectif)-1
%         for jj = ii+1:length(Eeffectif)
%             sigtxt = DoParamStats(Eeffectif{ii},Eeffectif{jj});
%             if ~strcmp(sigtxt,'NS')
%                 plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
%                 text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
%                 plotheight = plotheight+ 0.05*distsize;
%             end
%         end
%     end
%     ylim([0 plotheight])
% end

ax = gca;
ax.XTickLabelRotation = 45;

% avg Eeffectif per cell

figure(i_FigEeffPerCell)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
ax.YScale = 'log';
box on

for kcond = 1:ncond
    EeffectifCell{kcond} = [];
    for kc = 1:nCell
        EeffectifCell{kcond} = [EeffectifCell{kcond} nanmean(CellVsCond_Eeffectif{kcond,kc})];
    end
end

plotSpread_V(EeffectifCell,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(EeffectifCell)
    plot([ii-0.45 ii+0.45],[nanmedian(EeffectifCell{ii})  nanmedian(EeffectifCell{ii})],'k--','linewidth',1.3)
end
ylabel('K_AFM per cell (N/m) [ = E*H^2/R]')

if SigDisplay
    distsize = prctile(horzcat(EeffectifCell{:}),98);
    plotheight = distsize;
    for ii = 1:length(EeffectifCell)-1
        for jj = ii+1:length(EeffectifCell)
            sigtxt = DoParamStats(EeffectifCell{ii},EeffectifCell{jj});
            if ~strcmp(sigtxt,'NS')
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight * 10^(0.03), sigtxt,'HorizontalAlignment','center','fontsize',13) %+0.005*distsize
                plotheight = plotheight * 10^(0.2);
            else
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight * 10^(0.09), sigtxt,'HorizontalAlignment','center','fontsize',13) %+0.005*distsize
                plotheight = plotheight * 10^(0.2);
            end
        end
    end
    ylim([0 plotheight])
end


ax = gca;
ax.XTickLabelRotation = 45;





% Hysteresis

figure(i_FigHysteresis)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(Hyst,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(Hyst)
    plot([ii-0.45 ii+0.45],[nanmedian(Hyst{ii})  nanmedian(Hyst{ii})],'k--','linewidth',1.3)
end

% % first comp of each cell in another color
% set(0,'DefaultLineMarkerSize',10);
% plotSpread_V(HystFirst,'distributionMarkers',Sym,'distributionColors','y','xNames',Lab,'spreadWidth',1)
% set(0,'DefaultLineMarkerSize',8);


ylabel('Hysteresis')


if SigDisplay
    
    distsize = prctile(vertcat(Hyst{:}),95);
    plotheight = distsize;
    
    for ii = 1:length(E0chad)-1
        for jj = ii+1:length(Hyst)
            
            sigtxt = DoParamStats(Hyst{ii},Hyst{jj});
            
            if ~strcmp(sigtxt,'NS')
                
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.015*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                
                plotheight = plotheight+ 0.05*distsize;
            end
            
        end
    end
    
    ylim([0 plotheight])
    
end
ax = gca;
ax.XTickLabelRotation = 45;

%% Time Evolution

figure(i_FigEandH0vsTime)
hold on

for kcond = 1:ncond
    
    subplot(211)
    hold on
    plot(CompTimes{kcond},E0chad{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
    ylabel('Echad (kPa)')
    subplot(212)
    hold on
    plot(CompTimes{kcond},H0EXP{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
    ylabel('H0EXP (nm)')
    
end


%% Fluo quantif
    
    figure(i_FigFluo)
    
    datesList = unique(FluoResults.date);
    nDates = length(datesList);
    plotSymbols = ['o', '>', 's', 'd', '<'];
    for i=1:nDates
        currentDateMask = strcmp(FluoResults.date, datesList(i));
        subplot(1,2,1)
        hold on
        plot(FluoResults.meanFluoLevel(currentDateMask),FluoResults.meanEChadwick(currentDateMask),['k' plotSymbols(i)],'markerfacecolor','r')
        subplot(1,2,2)
        hold on
        plot(FluoResults.meanFluoLevel(currentDateMask),FluoResults.meanSurroundingThickness(currentDateMask),['k' plotSymbols(i)],'markerfacecolor','m')
        %plot(FluoResults.meanFluoLevel(currentDateMask),FluoResults.meanH0Chadwick(currentDateMask),['k' plotSymbols(i)],'markerfacecolor','c')
    end
    
    subplot(1,2,1)
    hold on
    xlabel('Fluo Intensity (a.u.)')
    ylabel('average Echad (kPa)')
    ax = gca;
    ax.XTickLabelRotation = 45;

    subplot(1,2,2)
    hold on
    xlabel('Fluo Intensity (a.u.)')
    ylabel('average H0 (nm)')
    currentLegend = legend;
    currentLegend.String = {'H0 btw ramps'}; %  {'H0 from fit' 'H0 btw ramps'};
    ax = gca;
    ax.XTickLabelRotation = 45;
    
    






%% per cell plot vs tps comp
if strcmp(answer,'Yes')
    mkdir(sffc)
    for kc = 1:nCell
        
        figure
        hold on
        subplot(211)
        title(FullCellList(kc))
        plotSpread_V({CellVsCond_Echad{:,kc}},'distributionMarkers',Sym,'distributionColors',Col,'xNames',[])
        ylabel('Echad (kPa)')
        ax = gca;
        ax.XTickLabel = [];
        
        subplot(212)
        plotSpread_V({CellVsCond_H0EXP{:,kc}},'distributionMarkers',Sym,'distributionColors',Col,'xNames',CondID)
        ylabel('H0exp (nm)')
        
        ax = gca;
        ax.XTickLabelRotation = 45;
        ax.TickLabelInterpreter = 'none';
        
        fig = gcf;
        saveas(fig,strcat(sffc, filesep, 'Echad_',FullCellList(kc),'.png'),'png')
        saveas(fig,strcat(sffc, filesep, 'Echad_',FullCellList(kc),'.fig'),'fig')
        
        close
    end
end

%% Saving

CondToWrite = cell(size(Cond,2)+1,ncond);

CondToWrite(1,:) = Lab;
CondToWrite(2:size(Cond,2)+1,:) = Cond';

Format = [repmat(['%s '], 1, size(Cond,2)+1) '\n'];

fid = fopen([sff filesep 'ConditionsList.txt'],'w+');
fprintf(fid,Format,CondToWrite{:});
fclose(fid);

if ncond <=2
    figure(i_FigLinEchadDistrib)
    fig = gca;
    saveas(fig,[sff filesep 'EchadDist.png'],'png')
    saveas(fig,[sff filesep 'EchadDist.fig'],'fig')
    
    figure(i_FigLogEchadDistrib)
    fig = gca;
    saveas(fig,[sff filesep 'EchadDistW.png'],'png')
    saveas(fig,[sff filesep 'EchadDistW.fig'],'fig')
    
    
    
end

if ncond >=2
    
    figure(i_FigEvH0)
    fig = gca;
    saveas(fig,[sff filesep 'EchadVh0.png'],'png')
    saveas(fig,[sff filesep 'EchadVh0.fig'],'fig')
    
    figure(i_FigEvHSUR)
    fig = gca;
    saveas(fig,[sff filesep 'EchadVhSur.png'],'png')
    saveas(fig,[sff filesep 'EchadVhSur.fig'],'fig')
    
    figure(i_FigHFITvHSUR)
    fig = gca;
    saveas(fig,[sff filesep 'hFitVhSur.png'],'png')
    saveas(fig,[sff filesep 'hFitVhSur.fig'],'fig')
    
    figure(i_FigEchadDistrib)
    fig = gca;
    saveas(fig,[sff filesep 'E0chad.png'],'png')
    saveas(fig,[sff filesep 'E0chad.fig'],'fig')
    
    figure(i_FigEchadPerCellDistrib)
    fig = gca;
    saveas(fig,[sff filesep 'E0chadCell.png'],'png')
    saveas(fig,[sff filesep 'E0chadCell.fig'],'fig')
    
    figure(i_FigH0Distrib)
    fig = gca;
    saveas(fig,[sff filesep 'H0ChadBeeswarm.png'],'png')
    saveas(fig,[sff filesep 'H0ChadBeeswarm.fig'],'fig')
    
    figure(i_FigSurroundingThicknessDistrib)
    fig = gca;
    saveas(fig,[sff filesep 'SurroundingThicknessBeeswarm.png'],'png')
    saveas(fig,[sff filesep 'SurroundingThicknessBeeswarm.fig'],'fig')
    
    figure(i_FigH0SURPerCellDistrib)
    fig = gca;
    saveas(fig,[sff filesep 'SurroundingThicknessPerCellBeeswarm.png'],'png')
    saveas(fig,[sff filesep 'SurroundingThicknessPerCellBeeswarm.fig'],'fig')
    
    figure(i_FigFittedThicknessDistrib)
    fig = gca;
    saveas(fig,[sff filesep 'FitThicknessChadBeeswarm.png'],'png')
    saveas(fig,[sff filesep 'FitThicknessChadBeeswarm.fig'],'fig')
    
    figure(i_FigH0FITPerCellDistrib)
    fig = gca;
    saveas(fig,[sff filesep 'FitThicknessChadPerCellBeeswarm.png'],'png')
    saveas(fig,[sff filesep 'FitThicknessChadPerCellBeeswarm.fig'],'fig')
    
    figure(i_FigEandH0vsTime)
    fig = gca;
    saveas(fig,[sff filesep 'E+H0InTime.png'],'png')
    saveas(fig,[sff filesep 'E+H0InTime.fig'],'fig')
    
    figure(i_FigHysteresis)
    fig = gca;
    saveas(fig,[sff filesep 'Hysteresis.png'],'png')
    saveas(fig,[sff filesep 'Hysteresis.fig'],'fig')
    
    figure(i_FigFluo)
    fig = gca;
    saveas(fig,[sff filesep 'FluoIntVsEandH.png'],'png')
    saveas(fig,[sff filesep 'FluoIntVsEandH.fig'],'fig')
    
    figure(i_FigEeff)
    fig = gca;
    saveas(fig,[sff filesep 'Eeff.png'],'png')
    saveas(fig,[sff filesep 'Eeff.fig'],'fig')
    
    figure(i_FigEeffPerCell)
    fig = gca;
    saveas(fig,[sff filesep 'EeffPerCell.png'],'png')
    saveas(fig,[sff filesep 'EeffPerCell.fig'],'fig')
    
    
end


end


%% extra functions
function [txt,p] = DoParamStats(A,B)

p = ranksum(A,B);

if 5/100 > p && p > 1/100
    txt = ['*']  ;
elseif 1/100 > p && p > 1/1000
    txt = ['**'];
elseif 1/1000 > p
    txt = ['***'];
else
    txt = ['NS'];
    
end
end

