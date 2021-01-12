function Meca2Plot_BigTable(resfolder,figurefolder,Cond,Col,Sym,Lab,savename,SigDisplay,varargin)

% Meca2Plot_BigTable(resfolder,figurefolder,conditions,colors,labels)


%% INIT

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',15)
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
    ptrfirst = BigTable(:,'CompNum').Variables == 1;
    
    E0chad{kcond} = BigTable(ptrcond,'EChadwick').Variables/1000;
    E0chadFirst{kcond} = BigTable(ptrcond&ptrfirst,'EChadwick').Variables/1000;
    EchadCI{kcond} = BigTable(ptrcond,'CiEChadwick').Variables/1000;
    
    Hyst{kcond} = BigTable(ptrcond,'Hysteresis').Variables;
    HystFirst{kcond} = BigTable(ptrcond&ptrfirst,'Hysteresis').Variables;
    
    H0EXP{kcond} = BigTable(ptrcond,'InitialThickness').Variables;
    H0EXPFirst{kcond} = BigTable(ptrcond&ptrfirst,'InitialThickness').Variables;
    H0SUR{kcond} = BigTable(ptrcond,'SurroundingThickness').Variables;
    H0BEF{kcond} = BigTable(ptrcond,'PreviousThickness').Variables;
    H0FIT{kcond} = BigTable(ptrcond,'H0Chadwick').Variables;
    
    CellIdList{kcond} = BigTable(ptrcond,'CellName').Variables;
    UniqueCellList{kcond} = unique(CellIdList{kcond});
    
    CompNums{kcond} = BigTable(ptrcond,'CompNum').Variables;
    CompTimes{kcond} = BigTable(ptrcond,'CompTime').Variables;
    
    
    Deltas{kcond} = BigTable(ptrcond,'Delta').Variables;
    Forces{kcond} = BigTable(ptrcond,'Force').Variables;   
    H0el{kcond} = BigTable(ptrcond,'H0el').Variables;   
    
    
end

FullCellList = unique(vertcat(UniqueCellList{:}));

nCell = length(FullCellList);



  %% delta� vs force
    for kcond = 1:ncond
        ncomp = length(Deltas{kcond});
    for kcomp = 1:ncomp
        
        DefMax(kcomp) = max(Deltas{kcond}{kcomp})/H0el{kcond}(kcomp);
        
        InterpForce(kcomp,:) = interp1((Deltas{kcond}{kcomp}).^2./max((Deltas{kcond}{kcomp}).^2),Forces{kcond}{kcomp},0.1:0.01:1);
        
    end
    
    MeanForce = nanmean(InterpForce,1);    
    StdForce = nanstd(InterpForce,1);
    
    ptrG1 = find(DefMax>0 & DefMax<0.33);
    ptrG2 = find(DefMax>0.34 & DefMax<0.66);
    ptrG3 = find(DefMax>0.67 & DefMax<1);
    
    figure(100+kcond)
    hold on
    
    subplot(2,4,1:2)
    hold on
    plot(0.1:0.01:1,MeanForce,'linewidth',2,'color',[0.5 0.5 0.5])
    plot(0.1:0.01:1,MeanForce+StdForce,'--','linewidth',0.5,'color',[0.5 0.5 0.5])
    plot(0.1:0.01:1,MeanForce-StdForce,'--','linewidth',0.5,'color',[0.5 0.5 0.5])
    fill([0.1:0.01:1 1:-0.01:0.1 0.1],[MeanForce+StdForce MeanForce(end:-1:1)-StdForce(end:-1:1) ...
    MeanForce(1)+StdForce(1)],[0.5 0.5 0.5],'facealpha',0.2,'edgecolor','none')
    xlabel('Normalized Delta�')
    ylabel('Force (pN)')
    xlim([0.1 1])
    
    % group 1
    nG11 = round(rand(1)*length(ptrG1));
    DefMaxtmpG11 = DefMax(ptrG1(nG11));
    Delta2G11tmp = Deltas{kcond}{ptrG1(nG11)}.^2;
    ForceG11tmp = Forces{kcond}{ptrG1(nG11)};
    
    nG12 = round(rand(1)*length(ptrG1));
    DefMaxtmpG12 = DefMax(ptrG1(nG12));
    Delta2G12tmp = Deltas{kcond}{ptrG1(nG12)}.^2;
    ForceG12tmp = Forces{kcond}{ptrG1(nG12)};
    
    subplot(2,4,3)
    hold on
    xlabel('Delta� (nm�)')
%     ylabel('Force (pN)')
    title(['MaxDef = ' num2str(round(100*DefMaxtmpG11)) '%'])
    fitlin = fitlm(Delta2G11tmp,ForceG11tmp);
    plot(Delta2G11tmp,ForceG11tmp,'-ok','markerfacecolor',0.5*[1 1 0])
    plot(Delta2G11tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G11tmp,'-b','linewidth',2)
  
    subplot(2,4,4)
    hold on
    xlabel('Delta� (nm�)')
%     ylabel('Force (pN)')
    title(['MaxDef = ' num2str(round(100*DefMaxtmpG12)) '%'])
    fitlin = fitlm(Delta2G12tmp,ForceG12tmp);
    plot(Delta2G12tmp,ForceG12tmp,'-ok','markerfacecolor',0.8*[1 1 0])
    plot(Delta2G12tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G12tmp,'-b','linewidth',2)
   
    % group 2
     nG21 = round(rand(1)*length(ptrG2));
    DefMaxtmpG21 = DefMax(ptrG2(nG21));
    Delta2G21tmp = Deltas{kcond}{ptrG2(nG21)}.^2;
    ForceG21tmp = Forces{kcond}{ptrG2(nG21)};
    
    nG22 = round(rand(1)*length(ptrG2));
    DefMaxtmpG22 = DefMax(ptrG2(nG22));
    Delta2G22tmp = Deltas{kcond}{ptrG2(nG22)}.^2;
    ForceG22tmp = Forces{kcond}{ptrG2(nG22)};
    
    subplot(2,4,5)
    hold on
    xlabel('Delta� (nm�)','fontsize',15)
    ylabel('Force (pN)')
    title(['MaxDef = ' num2str(round(100*DefMaxtmpG21)) '%'])
    fitlin = fitlm(Delta2G21tmp,ForceG21tmp);
    plot(Delta2G21tmp,ForceG21tmp,'-ok','markerfacecolor',0.5*[1 0 1])
    plot(Delta2G21tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G21tmp,'-b','linewidth',2)
  
    subplot(2,4,6)
    hold on
    xlabel('Delta� (nm�)')
%     ylabel('Force (pN)')
   title(['MaxDef = ' num2str(round(100*DefMaxtmpG22)) '%'])
    fitlin = fitlm(Delta2G22tmp,ForceG22tmp);
    plot(Delta2G22tmp,ForceG22tmp,'-ok','markerfacecolor',0.8*[1 0 1])
    plot(Delta2G22tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G22tmp,'-b','linewidth',2)
   
    % group 3
     nG31 = round(rand(1)*length(ptrG3));
    DefMaxtmpG31 = DefMax(ptrG3(nG31));
    Delta2G31tmp = Deltas{kcond}{ptrG3(nG31)}.^2;
    ForceG31tmp = Forces{kcond}{ptrG3(nG31)};
    
    nG32 = round(rand(1)*length(ptrG3));
    DefMaxtmpG32 = DefMax(ptrG3(nG32));
    Delta2G32tmp = Deltas{kcond}{ptrG3(nG32)}.^2;
    ForceG32tmp = Forces{kcond}{ptrG3(nG32)};
    
    subplot(2,4,7)
    hold on
    xlabel('Delta� (nm�)')
%     ylabel('Force (pN)')
    title(['MaxDef = ' num2str(round(100*DefMaxtmpG31)) '%'])
    fitlin = fitlm(Delta2G31tmp,ForceG31tmp);
    plot(Delta2G31tmp,ForceG31tmp,'-ok','markerfacecolor',0.5*[0 1 1])
    plot(Delta2G31tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G31tmp,'-b','linewidth',2)
  
    subplot(2,4,8)
    hold on
    xlabel('Delta� (nm�)')
%     ylabel('Force (pN)')
   title(['MaxDef = ' num2str(round(100*DefMaxtmpG32)) '%'])
    fitlin = fitlm(Delta2G32tmp,ForceG32tmp);
    plot(Delta2G32tmp,ForceG32tmp,'-ok','markerfacecolor',0.8*[0 1 1])
    plot(Delta2G32tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G32tmp,'-b','linewidth',2)
   
    
    
 %% old   
    
%     for k = 1:3
%     nG1 = round(rand(1)*length(ptrG1));
%     Delta2G1tmp = Deltas{kcond}{ptrG1(nG1)}.^2;
%     ForceG1tmp = Forces{kcond}{ptrG1(nG1)};
%     col1 = (0.4+k/10)*[1 1 0];
%     col2 = (0.4+k/10)*[0 1 1];
%     col3 = (0.4+k/10)*[1 0 1];
%     
%     nG2 = round(rand(1)*length(ptrG2));
%     Delta2G2tmp = Deltas{kcond}{ptrG2(nG2)}.^2;
%     ForceG2tmp = Forces{kcond}{ptrG2(nG2)};
%     
%     nG3 = round(rand(1)*length(ptrG3));
%     Delta2G3tmp = Deltas{kcond}{ptrG3(nG3)}.^2;
%     ForceG3tmp = Forces{kcond}{ptrG3(nG3)};
%    
%     subplot(2,2,1)
%     hold on
%     xlabel('Normalized Delta�')
%     ylabel('Force (pN)')
%     title('Delta� normalized by max')
%     plot(Delta2G1tmp/max(Delta2G1tmp),ForceG1tmp,'-ok','markerfacecolor',col1)
%     plot(Delta2G2tmp/max(Delta2G2tmp),ForceG2tmp,'-ok','markerfacecolor',col2)
%     plot(Delta2G3tmp/max(Delta2G3tmp),ForceG3tmp,'-ok','markerfacecolor',col3)
%     
%     subplot(2,2,2)
%     hold on
%     xlabel('Delta� (nm�)')
%     ylabel('Force (pN)')
%     title('Max Deformation 0-33%')
%     fitlin = fitlm(Delta2G1tmp,ForceG1tmp);
%     plot(Delta2G1tmp,ForceG1tmp,'-ok','markerfacecolor',col1)
%     plot(Delta2G1tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G1tmp,'-b','linewidth',2)
%                   
%     subplot(2,2,3)
%     hold on
%     xlabel('Delta� (nm�)')
%     ylabel('Force (pN)')
%     title('Max Deformation 34-66%')
%     fitlin = fitlm(Delta2G2tmp,ForceG2tmp);
%     plot(Delta2G2tmp,ForceG2tmp,'-ok','markerfacecolor',col2)
%     plot(Delta2G2tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G2tmp,'-b','linewidth',2)
%    
%     subplot(2,2,4)
%     hold on
%     xlabel('Delta� (nm�)')
%     ylabel('Force (pN)')
%     title('Max Deformation 66-100%')
%     fitlin = fitlm(Delta2G3tmp,ForceG3tmp);
%     plot(Delta2G3tmp,ForceG3tmp,'-ok','markerfacecolor',col3)
%     plot(Delta2G3tmp,fitlin.Coefficients.Estimate(1)+fitlin.Coefficients.Estimate(2)*Delta2G3tmp,'-b','linewidth',2)
%     end             
%     
       
    end

return
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
xlabel('Thickness before comp(nm)')
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



% all data
plotSpread_V(E0chad,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(E0chad)
    plot([ii-0.45 ii+0.45],[nanmedian(E0chad{ii})  nanmedian(E0chad{ii})],'k--','linewidth',1.3)
end

% % first comp of each cell in another color
% set(0,'DefaultLineMarkerSize',10);
% plotSpread_V(E0chadFirst,'distributionMarkers',Sym,'distributionColors','y','xNames',Lab,'spreadWidth',1)
% set(0,'DefaultLineMarkerSize',8);
% 
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
H0BeeswarmFirst = H0EXPFirst;

figure(6)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(H0Beeswarm,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(H0Beeswarm)
    plot([ii-0.45 ii+0.45],[nanmedian(H0Beeswarm{ii}) nanmedian(H0Beeswarm{ii})],'k--','linewidth',1.3)
end

% first comp of each cell in another color
set(0,'DefaultLineMarkerSize',10);
plotSpread_V(H0BeeswarmFirst,'distributionMarkers',Sym,'distributionColors','y','xNames',Lab,'spreadWidth',1)
set(0,'DefaultLineMarkerSize',8);

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


% Hysteresis

figure(8)
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

plotSpread_V(Hyst,'distributionMarkers',Sym,'distributionColors',Col,'xNames',Lab,'spreadWidth',1)
for ii = 1:length(Hyst)
    plot([ii-0.45 ii+0.45],[nanmedian(Hyst{ii})  nanmedian(Hyst{ii})],'k--','linewidth',1.3)
end

% first comp of each cell in another color
set(0,'DefaultLineMarkerSize',10);
plotSpread_V(HystFirst,'distributionMarkers',Sym,'distributionColors','y','xNames',Lab,'spreadWidth',1)
set(0,'DefaultLineMarkerSize',8);


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

figure(7)
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

