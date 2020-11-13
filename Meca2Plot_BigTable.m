function Meca2Plot_BigTable(resfolder,figurefolder,Cond,Col,Sym,Lab,savename,SigDisplay,varargin)

% Meca2Plot_BigTable(resfolder,figurefolder,conditions,colors,labels)


%% INIT

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextInterpreter','Latex')

here = pwd;

%% params and settings

% dossier de donn�es et sauvegarde des donn�es
path = [resfolder filesep 'D2M'];
datenow = datestr(now,'yy-mm-dd');


sfig = figurefolder;
sff = [sfig filesep datenow filesep 'Meca' filesep savename];
sffc = [sff filesep 'PerCell'];
mkdir(sffc)


%% get data

load([resfolder filesep 'MecaDataTable'])

if ~(size(Cond,1) == length(Col) && size(Cond,1) == length(Lab) && length(Lab) == length(Col))
    error('Numbers of conditions, colors, and lables don''t match !')
end


% varargin handling

excludevect = zeros(height(BigTable),1);

if ~isempty(varargin)
    
    if ismember('Exclude',varargin(1:2:end))
        
        
        
        ptrex = find(ismember('Exclude',varargin(1:2:end)));
        
        for kex = 1:length(ptrex)
            condex = varargin{ptrex(kex)+1};
            
            excludevecttmp = strcmp(BigTable(:,condex{1}).Variables,condex{2});
            excludevect = excludevect | excludevecttmp;
            
        end
    end
end


%list of valid datas and info on non valid
VALIDptr = BigTable.Validated & ~excludevect;

ptrexptypes = find(ismember(Cond(1,:),'ExpType'))+1;

exptypesList = Cond(:,ptrexptypes);

ptrfitP= find(ismember(Cond(1,:),'FitParams'))+1;

fitPList = Cond(:,ptrfitP);

DoSuccessStats(unique(exptypesList),unique(fitPList),excludevect)

% identifying changes between conditions

CondID = [];

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
    
    E0chad{kcond} = BigTable(ptrcond,'EChadwick').Variables/1000;
    EchadCI{kcond} = BigTable(ptrcond,'CiEChadwick').Variables/1000;
    
    Hyst{kcond} = BigTable(ptrcond,'Hysteresis').Variables;
    
    H0EXP{kcond} = BigTable(ptrcond,'InitialThickness').Variables;
    H0SUR{kcond} = BigTable(ptrcond,'SurroundingThickness').Variables;
    H0BEF{kcond} = BigTable(ptrcond,'PreviousThickness').Variables;
    H0FIT{kcond} = BigTable(ptrcond,'H0Chadwick').Variables;
    
    CellIdList{kcond} = BigTable(ptrcond,'CellName').Variables;
    UniqueCellList{kcond} = unique(CellIdList{kcond});
    
    CompNums{kcond} = BigTable(ptrcond,'CompNum').Variables;
    CompTimes{kcond} = BigTable(ptrcond,'CompTime').Variables;
    
end

FullCellList = unique(vertcat(UniqueCellList{:}));

nCell = length(FullCellList);

%% prep figs

if ncond <=2 % figures for 1 or less datasets
    
    figure(1)
    hold on
    title('Linlin distrib + median')
    L1 = legend;
    xlabel('Echad (kPa)')
    ylabel('Count')
    
    
    figure(2)
    hold on
    title('Loglog distrib + weighted avg')
    L2 = legend;
    xlabel('Elastic modulus (kPa)')
    ylabel('Count')
    
    
    
end

figure(3)
hold on
title('E v H0')
L3 = legend;
xlabel('Thickness befor comp(nm)')
ylabel('E chad (kPa)')

figure(4)
hold on
title('Echad distrib')

figure(5)
hold on
title('Echad cell distrib')

figure(6)
hold on
title('H0 distrib')

figure(7)
subplot(211)
hold on
title('Echad & H0 vs Time')

figure(8)
hold on
title('Hysteresis')

ymax1 = 0;
ymax2 = 0;

%% Plots


for kcond = 1:ncond
    
    
    fprintf(['\n\nFor condition : ' Lab{kcond} '\n'])
    
    fprintf(['\nNumber of compressions : ' num2str(length(E0chad{kcond})) '.\nNumber of cells : ' num2str(size(UniqueCellList{kcond},1)) '.\n'])
    
    
    
    if ncond <=2
        
        %%  distribution des Echad + median
        
        [counts,edges] = histcounts(E0chad{kcond},'binwidth',2);
        
        figure(1)
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
        
        
        figure(2)
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
    
    
    %% Module vs Epaisseur d�part
    
    figure(3)
    hold on
    plot(H0EXP{kcond},E0chad{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
    
    if ncond == 2
        ratioEH{kcond} = (E0chad{kcond}./H0EXP{kcond})/median(E0chad{kcond});
    end
    
    ax = gca;
    ax.YScale = 'log';
    ax.XTickLabelRotation = 45;
    
    L3.String{end} = Lab{kcond};
    
end


%% regrouping per cells and conditions
CellVsCond_Echad = cell(ncond,nCell);
CellVsCond_H0EXP = cell(ncond,nCell);
CellVsCond_ID = cell(ncond,nCell);

for kc = 1:nCell
    for kcond = 1:ncond
        
        ptrcell = strcmp(CellIdList{kcond}, FullCellList(kc));
        
        CellVsCond_Echad{kcond,kc} = E0chad{kcond}(ptrcell);
        CellVsCond_H0EXP{kcond,kc} =H0EXP{kcond}(ptrcell);
        CellVsCond_ID{kcond,kc} = strcat(FullCellList(kc), CondID(kcond));
        
    end
    
    
end

%% Beeswarm plot

% all data Echad
figure(4)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(E0chad,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(E0chad)
    plot([ii-0.45 ii+0.45],[median(E0chad{ii})  median(E0chad{ii})],'k--','linewidth',1.3)
end
ylabel('Echad (kPa)')


if SigDisplay
    
    distsize = prctile(vertcat(E0chad{:}),98);
    plotheight = distsize;
    
    for ii = 1:length(E0chad)-1
        
        for jj = ii+1:length(E0chad)
            
            sigtxt = DoParamStats(E0chad{ii},E0chad{jj});
            
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


% avg per cell Echad
figure(5)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
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
ylabel('Echad par cell (kPa)')


if SigDisplay
    
    distsize = prctile(horzcat(E0chadCell{:}),98);
    plotheight = distsize;
    
    for ii = 1:length(E0chadCell)-1
        
        for jj = ii+1:length(E0chadCell)
            
            sigtxt = DoParamStats(E0chadCell{ii},E0chadCell{jj});
            
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

% H0


H0Beeswarm = H0EXP;

figure(6)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(H0Beeswarm,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(H0Beeswarm)
    plot([ii-0.45 ii+0.45],[median(H0Beeswarm{ii})  median(H0Beeswarm{ii})],'k--','linewidth',1.3)
end
ylabel('H0 init ramp')


if SigDisplay
    
    distsize = prctile(vertcat(H0Beeswarm{:}),90);
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


% Hysteresis

figure(8)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(Hyst,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(Hyst)
    plot([ii-0.45 ii+0.45],[median(Hyst{ii})  median(Hyst{ii})],'k--','linewidth',1.3)
end
ylabel('Hysteresis')


if SigDisplay
    
    distsize = prctile(vertcat(Hyst{:}),90);
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

figure(7)
hold on

for kcond = 1:ncond
    
    subplot(211)
    hold on
    plot(CompTimes{kcond},E0chad{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
    subplot(212)
    hold on
    plot(CompTimes{kcond},H0EXP{kcond},['k' Sym{kcond}],'markerfacecolor',Col{kcond})
    
end

%% per cell plot vs tps comp
if strcmp(answer,'Yes')
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
    figure(1)
    fig = gca;
    saveas(fig,[sff filesep 'EchadDist.png'],'png')
    saveas(fig,[sff filesep 'EchadDist.fig'],'fig')
    
    figure(2)
    fig = gca;
    saveas(fig,[sff filesep 'EchadDistW.png'],'png')
    saveas(fig,[sff filesep 'EchadDistW.fig'],'fig')
    
    
    
end

if ncond >=2
    
    figure(3)
    fig = gca;
    saveas(fig,[sff filesep 'EchadVh0.png'],'png')
    saveas(fig,[sff filesep 'EchadVh0.fig'],'fig')
    
    figure(4)
    fig = gca;
    saveas(fig,[sff filesep 'E0chad.png'],'png')
    saveas(fig,[sff filesep 'E0chad.fig'],'fig')
    
    figure(5)
    fig = gca;
    saveas(fig,[sff filesep 'E0chadCell.png'],'png')
    saveas(fig,[sff filesep 'E0chadCell.fig'],'fig')
    
    figure(6)
    fig = gca;
    saveas(fig,[sff filesep 'H0ChadBeeswarm.png'],'png')
    saveas(fig,[sff filesep 'H0ChadBeeswarm.fig'],'fig')
    
    figure(7)
    fig = gca;
    saveas(fig,[sff filesep 'E+H0InTime.png'],'png')
    saveas(fig,[sff filesep 'E+H0InTime.fig'],'fig')
    
    figure(8)
    fig = gca;
    saveas(fig,[sff filesep 'Hysteresis.png'],'png')
    saveas(fig,[sff filesep 'Hysteresis.fig'],'fig')
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

