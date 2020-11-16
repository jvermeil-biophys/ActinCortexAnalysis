function [CBot,MedCw,CTop] = Data2CortexStats(Cell,label,coord,col,datafolder)

% close all


%% INIT
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)

warning('off','all')

% Data location and choice

date=get_dates();

% path= 'C:\Users\Pierre\Desktop\Data\MatFile\V2D';
path= [datafolder filesep 'P2S'];






%% Data retrieval
nd = size(date,2);

  
    
    %% Data retrieval and peaks detection
    datafile=['P2S_' Cell '.mat'];
    
    
    % verification de l'existence du fichier
    E = exist([path filesep datafile]);
    
if E == 2
        
        cd(path)
        
        load(datafile)
        fprintf(['File loaded : \n' datafile '. \n\n']);

       CBot = MS.CortBots;
       CTop = MS.CortTops;
       MedCw = MS.MedCortWs;
       CTags = MS.CTags;
       
       MeanCw = MS.Cortmeans; % pour obtenir cortex thickness sans membrane et coating
       Cstd  = MS.Cortsdts;
           
end



% label2 = [label ' - ' num2str(length(CBot))];

figure(2000)

hold on
title(['Median Cortex Width'])
Spread_Julien(MedCw,coord,0.9,col,'o',75)
plot([coord-0.2 coord+0.2],[mean(MedCw) mean(MedCw)],'k','linewidth',3)
ylabel('Median Cortex Width (nm)')
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ax.XTick = 1:coord;
ax.XTickLabel{coord} = label ;
lims = ylim;

figure(2001)
hold on
title(['Median Cortex Width'])
bar(coord,mean(MedCw),'facecolor',col,'edgecolor','none')
errorbar(coord,mean(MedCw),std(CBot)/sqrt(length(MedCw)),'k')
ylabel('Median Cortex Width (nm)')
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ax.XTick = 1:coord;
ax.XTickLabel{coord} = label ;


 
figure(2002)
hold on
Edges = -150:30:800;
Counts = histcounts(MedCw,Edges,...
    'normalization','probability');
plot(Edges(1:end-1)+diff(Edges)/2,Counts,'color',col,'linewidth',2)
xlabel('Median Cortex Width (nm)')
ylabel('Percentage of cells')
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
leg = legend('location','best')
leg.String{end}= label2 ;

figure(2003)
hold on
title(['Median Cortex Width'])
Spread_Julien_bar(MedCw,coord,col,30)
plot([coord coord+0.2],[mean(MedCw) mean(MedCw)],'--k','linewidth',3)
ylabel('Median Cortex Width (nm)')
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ax.XTick = 1:coord;
ax.XTickLabel{coord} = label ;


figure(4000)
hold on
title('Median and 10-90 percentile')
errorbar(coord,mean(MedCw),std(MedCw)/sqrt(length(MedCw)),'.','color',col,'handlevisibility','off')
plot(coord,mean(MedCw),'s','markerfacecolor',col,'markeredgecolor','none','markersize',10,'handlevisibility','off')
errorbar(coord,mean(CTop),std(CTop)/sqrt(length(CTop)),'.','color',col,'handlevisibility','off')
plot(coord,mean(CTop),'v','markerfacecolor',col,'markeredgecolor','none','markersize',10,'handlevisibility','off')
errorbar(coord,mean(CBot),std(CBot)/sqrt(length(CBot)),'.','color',col,'handlevisibility','off')
plot(coord,mean(CBot),'^','markerfacecolor',col,'markeredgecolor','none','markersize',10,'handlevisibility','off')
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ax.XTick = 1:coord;
ax.XTickLabel{coord} = label ;


step = 0.9/(length(MedCw)+1);
x = (1:length(MedCw))*step-0.45;


Mat1= sortrows([MedCw' CBot' CTop']);
figure(4001)
hold on
title('Median and 10-90 percentile')
plot(coord+x,Mat1(:,1),'ks','markerfacecolor',col)
plot(coord+x,Mat1(:,3),'kv','markerfacecolor',col)
plot(coord+x,Mat1(:,2),'k^','markerfacecolor',col)
plot([coord+x; coord+x],[Mat1(:,2)'; Mat1(:,3)'],'-k'); 
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ax.XTick = 1:coord;
ax.XTickLabel{coord} = label ;



CAmp = CTop - CBot;
[corrmat p] = corrcoef(MedCw,CAmp);
figure(5000)
hold on
plot(MedCw,CAmp,'ok','markerfacecolor',col)
xlabel('Cortex Median (nm)')
ylabel('Cortex Amplitude (nm)')
ax = gca;
ax.FontSize = 20;
leg = legend;
leg.String{end} = [label ' - ' num2str(corrmat(1,2))];


[corrmat p] = corrcoef(MeanCw,Cstd);
figure(6000)
hold on
plot(MeanCw,Cstd,'ok','markerfacecolor',col)
xlabel('Cortex Mean (nm)')
ylabel('Cortex std (nm)')
ax = gca;
ax.FontSize = 20;
leg = legend;
leg.String{end} = [label ' - ' num2str(corrmat(1,2))];

end
