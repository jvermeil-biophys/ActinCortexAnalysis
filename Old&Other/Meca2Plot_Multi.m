function Meca2Plot_Multi(specif,color,legends,resfolder,figurefolder)

%% INIT

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextInterpreter','Latex')

here = pwd;

%% params and settings

% home
% h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';


% dossier de données et sauvegarde des données
path = [resfolder filesep 'D2M'];
datenow = datestr(now,'yy-mm-dd');




EchadBeeswarm = cell(1,length(specif));
H0ChadBeeswarm = cell(1,length(specif));

%% prep figs

if length(specif) <=2 % figures for 1 or less datasets
    
    figure(1)
    hold on
    title('Linlin distrib + median')
    L1 = legend;
    xlabel('Echad (kPa)')
    ylabel('Count')
    
    
    figure(3)
    hold on
    title('Loglog distrib + weighted avg')
    L3 = legend;
    xlabel('Elastic modulus (kPa)')
    ylabel('Count')
    
    
    
end

if length(specif) >=2 % figure for 2 or more datasets
    
    
    figure(2)
    hold on
    title('E v H0')
    L2 = legend;
    xlabel('Thickness at begining of comp from fit (nm)')
    ylabel('E chad (kPa)')
    
    
    figure(5)
    hold on
    title('Echad distrib')
    
    
    figure(6)
    hold on
    title('Echad per cell ditrib')
    
end


ymax1 = 0;
ymax3 = 0;

%% get data & plot

for kcond = 1:length(specif)
    
    
    E0chad = {};
    EchadCI = {};
    
    H0ABS= {};
    H0FIT= {};
    
    StrainR = {};
    
    
    %% Retrieve datafile
    List = ls(path);
    
    for i = 3:size(List,1)
        if contains(List(i,:),specif{kcond})
            
            fprintf(['Chargement du fichier ' List(i,:) '\n\n'])
            load([path filesep List(i,:)]);
            
            E0chad = {E0chad{:} Echad{:}};
            EchadCI = {EchadCI{:} EchadConf{:}};
            
            H0ABS= {H0ABS{:} H0abs{:}};
            H0FIT = {H0FIT{:} ChadsH0{:}};
            
            StrainR = {StrainR{:} StrainRateMed{:}};
            
            
        end
    end
    
    
    
    Echadplot = vertcat(E0chad{:})/1000;
    EchadCIvect = vertcat(EchadCI{:})/1000;
    
    EchadBeeswarm{kcond} = Echadplot;
    
    H0ABSplot = vertcat(H0ABS{:});
    H0FITplot = vertcat(H0FIT{:});
    
    H0ChadBeeswarm{kcond} = H0FITplot;
    
    StrainRplot = vertcat(StrainR{:})*100;
    
    clear EchadWcell nEchad  Echadcell StdEchadCell ErrEchadWcell
    
    % per cell with 3+ measurements
    icell = 1;
    for  k = 1:length(E0chad)
        if length(E0chad{k})>0
            EchadcellWeight = (E0chad{k}./EchadCI{k}).^2;
            Echadcell(icell) = mean(E0chad{k}); % avg per cell
            EchadWcell(icell) = sum(E0chad{k}.*EchadcellWeight)/sum(EchadcellWeight); % weighted avg per cell
            nEchad(icell) = length(E0chad{k});
            ErrEchadWcell(icell) = sqrt(sum(((E0chad{k}-EchadWcell(icell)).^2.*EchadcellWeight)))/sum(EchadcellWeight);%/sqrt(nEchad(icell)); % weighted err per cell
            StdEchadCell(icell) = std(E0chad{k});
            icell = icell + 1;
        else
            EchadcellWeight = [];
            Echadcell(icell) = NaN;
            EchadWcell(icell) = NaN;
            nEchad(icell) = NaN;
            ErrEchadWcell(icell) = NaN;
            StdEchadCell(icell) = NaN;
            icell = icell + 1;
        end
    end
    
    
    EchadBeeswarmCell{kcond} = EchadWcell/1000;
    
    
    % Sig per cell vs sig between cells
    
    fprintf(['\n\nFor condition : ' specif{kcond} '\n'])
    
    fprintf(['\nNumber of compressions : ' num2str(length(Echadplot)) '.\nNumber of cells : ' num2str(length(E0chad)) '.\n'])
    
    
    
    if length(specif) <=2
        
        %%  distribution des Echad + median
        
        [counts,edges] = histcounts(Echadplot,'binwidth',2);
        
        figure(1)
        hold on
        histogram(Echadplot,edges,'facecolor',color{kcond},'facealpha',0.5)
        plot([median(Echadplot) median(Echadplot)],[0 1000],'--','linewidth',1.7,'color',color{kcond},'handlevisibility','off')
        ax = gca;
        
        ymax1tmp = max(counts)*1.05;
        
        ymax1 = max([ymax1 ymax1tmp]);
        
        ylim([0 ymax1])
        
        
        %% Distrib Echad avec moyenne pondérée en log
        
        
        EchadWeight = (Echadplot./EchadCIvect).^2;
        
        EchadWeightedMean = 10^(sum(log10(Echadplot).*EchadWeight)/sum(EchadWeight));
        
        EchadWeightedStd = sqrt(sum(((Echadplot-EchadWeightedMean).^2.*EchadWeight))/sum(EchadWeight));
        
        
        [counts,edges] = histcounts(log10(Echadplot));
        
        
        figure(3)
        hold on
        histogram(Echadplot,10.^edges,'facecolor',color{kcond},'facealpha',0.5)
        
        plot([EchadWeightedMean EchadWeightedMean],[0 1000],'--','linewidth',1.7,'color',color{kcond},'handlevisibility','off')
        ax = gca;
        ax.LineWidth = 1.5;
        ax.XScale = 'log';
        ax.XTick = [1 10 100];
        box on
        
        
        xlim([0.5 200])
        
        
        ymax3tmp = max(counts)*1.05;
        
        ymax3 = max([ymax3 ymax3tmp]);
        
        ylim([0 ymax3])
        
        
        %% legends & other additions
        
        L1.String{end} = [legends{kcond}  ' - ' num2str(median(Echadplot)) 'kPa'];
        
        L3.String{end} = [legends{kcond}  ' - ' num2str(EchadWeightedMean) '+- ' num2str(EchadWeightedStd) ' kPa'];
        
    end
    
    
    %% Module vs Epaisseur départ
    
    figure(2)
    hold on
    plot(H0FITplot,Echadplot,'ok','markerfacecolor',color{kcond})
    
    if length(specif) == 2
        ratioEH{kcond} = (Echadplot./H0ABSplot)/median(Echadplot);
    end
    
    L2.String{end} = legends{kcond};
end

if length(specif) >=2
    
    
    
    %% Beeswarm plot
    
    % all data Echad
    figure(5)
    ax = gca;
    ax.YColor = [0 0 0];
    ax.LineWidth = 1.5;
    box on
    
    plotSpread_V(EchadBeeswarm,'distributionMarkers',{'o'},'distributionColors',color,'xNames',legends,'spreadWidth',1)
    for ii = 1:length(EchadBeeswarm)
        plot([ii-0.45 ii+0.45],[median(EchadBeeswarm{ii})  median(EchadBeeswarm{ii})],'k--','linewidth',1.3)
    end
    ylabel('Echad (kPa)')
    
    
    
    distsize = prctile(horzcat(EchadBeeswarmCell{:}),98);
    plotheight = distsize;
    
    for ii = 1:length(EchadBeeswarm)-1
        
        for jj = ii+1:length(EchadBeeswarm)
            
            sigtxt = DoParamStats(EchadBeeswarm{ii},EchadBeeswarm{jj});
            
            if ~strcmp(sigtxt,'NS')
                
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.1*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                
                plotheight = plotheight + 0.2*distsize;
            end
        end
    end
    
    % per cell Echad
    
    figure(6)
    ax = gca;
    ax.YColor = [0 0 0];
    ax.LineWidth = 1.5;
    box on
    
    plotSpread_V(EchadBeeswarmCell,'distributionMarkers',{'o'},'distributionColors',color,'xNames',legends,'spreadWidth',1)
    
    ylabel('Echad - WtdAvg per cell (kPa)')
    
    distsize = prctile(horzcat(EchadBeeswarmCell{:}),90);
    plotheight = distsize;
    
    for ii = 1:length(EchadBeeswarm)-1
        for jj = ii+1:length(EchadBeeswarmCell)
            
            sigtxt = DoParamStats(EchadBeeswarmCell{ii},EchadBeeswarmCell{jj});
            
            if ~strcmp(sigtxt,'NS')
                
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.1*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                
                plotheight = plotheight+ 0.2*distsize;
            end
        end
    end
    
    
    % H0 chadwick
    
    figure(7)
    ax = gca;
    ax.YColor = [0 0 0];
    ax.LineWidth = 1.5;
    box on
    
    plotSpread_V(H0ChadBeeswarm,'distributionMarkers',{'o'},'distributionColors',color,'xNames',legends,'spreadWidth',1)
    for ii = 1:length(H0ChadBeeswarm)
        plot([ii-0.45 ii+0.45],[median(H0ChadBeeswarm{ii})  median(H0ChadBeeswarm{ii})],'k--','linewidth',1.3)
    end
    ylabel('H0chad')
    
    distsize = prctile(vertcat(H0ChadBeeswarm{:}),90);
    plotheight = distsize;
    
    for ii = 1:length(EchadBeeswarm)-1
        for jj = ii+1:length(H0ChadBeeswarm)
            
            sigtxt = DoParamStats(H0ChadBeeswarm{ii},H0ChadBeeswarm{jj});
            
            if ~strcmp(sigtxt,'NS')
                
                plot([ii+0.05 jj-0.05],[plotheight plotheight],'k-','linewidth',1.5)
                text((ii+jj)/2, plotheight + 0.1*distsize, sigtxt,'HorizontalAlignment','center','fontsize',13)
                
                plotheight = plotheight+ 0.2*distsize;
            end
            
        end
    end
    
end




%% Saving

sfig = figurefolder;
sff = [sfig filesep datenow filesep 'Meca' filesep [legends{:}] params];
mkdir(sff)

if length(specif) <=2
    figure(1)
    fig = gca;
    saveas(fig,[sff filesep 'EchadDist.png'],'png')
    saveas(fig,[sff filesep 'EchadDist.fig'],'fig')
    
    figure(3)
    fig = gca;
    saveas(fig,[sff filesep 'EchadDistW.png'],'png')
    saveas(fig,[sff filesep 'EchadDistW.fig'],'fig')
    
    
    
end

if length(specif) >=2
    
    figure(2)
    fig = gca;
    saveas(fig,[sff filesep 'EchadVh0.png'],'png')
    saveas(fig,[sff filesep 'EchadVh0.fig'],'fig')
    
    figure(5)
    fig = gca;
    saveas(fig,[sff filesep 'EchadBeeswarm.png'],'png')
    saveas(fig,[sff filesep 'EchadBeeswarm.fig'],'fig')
    
    figure(6)
    fig = gca;
    saveas(fig,[sff filesep 'EchadBeeswarmCell.png'],'png')
    saveas(fig,[sff filesep 'EchadBeeswarmCell.fig'],'fig')
    
    figure(7)
    fig = gca;
    saveas(fig,[sff filesep 'H0ChadBeeswarm.png'],'png')
    saveas(fig,[sff filesep 'H0ChadBeeswarm.fig'],'fig')
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

