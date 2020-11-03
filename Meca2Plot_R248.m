function Meca2Plot_R248(specif,color,CLOSE)

%% INIT

if CLOSE
    close all
end

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextInterpreter','Latex')

here = pwd;

%% params and settings

% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';


% dossier de données et sauvegarde des données
path = 'D:\Data\MatFile\D2M';
datenow = datestr(now,'yy-mm-dd');

% sauvergarde des figures
sfig = 'D:\Data\Figures';
sff = [sfig filesep datenow filesep 'Meca' filesep specif];
mkdir(sff)

CortThick = [];
E = {[];[];[]};
FinalSig = {[];[];[]};
E016mT = {[];[];[]};
Alpha = {[];[];[]};
H0PRC = {[];[];[]};
H0ABS= {[];[];[]};
CompN = {[];[];[]};
Hyst = {[];[];[]};


%% Retrieve datafile
List = ls(path);

for i = 3:size(List,1)
    if contains(List(i,:),specif)
        
        fprintf(['Chargement du fichier ' List(i,:) '\n\n'])
        load([path filesep List(i,:)]);
        
        
        CortThick = [CortThick Cort];
        for i = 1:3
        E{i,:} = [vertcat(E{i,:}); vertcat(E0{i,:})];
        FinalSig{i,:} = [vertcat(FinalSig{i,:}); vertcat(finalSig{i,:})];
        E016mT{i,:} = [vertcat(E016mT{i,:}); vertcat(E16mT{i,:})];
        Alpha{i,:} = [vertcat(Alpha{i,:}); vertcat(AlphaNL{i,:})];
        H0PRC{i,:} = [vertcat(H0PRC{i,:}); vertcat(H0prc{i,:})];
        H0ABS{i,:}= [vertcat(H0ABS{i,:}); vertcat(H0abs{i,:})];
        CompN{i,:} = [vertcat(CompN{i,:}); vertcat(CompNum{i,:})];
        Hyst{i,:} = [vertcat(Hyst{i,:}); vertcat(Hys{i,:})];
        end
        
    end
end

Eplot2 = vertcat(E{1,:})/1000;
Eplot4 = vertcat(E{2,:})/1000;
Eplot8 = vertcat(E{3,:})/1000;

FinalSigplot2 = vertcat(FinalSig{1,:});
FinalSigplot4 = vertcat(FinalSig{2,:});
FinalSigplot8 = vertcat(FinalSig{3,:});


E16mTplot2 = vertcat(E016mT{1,:})/1000;
E16mTplot4 = vertcat(E016mT{2,:})/1000;
E16mTplot8 = vertcat(E016mT{3,:})/1000;

Hyst2 = vertcat(Hyst{1,:});
Hyst4 = vertcat(Hyst{2,:});
Hyst8 = vertcat(Hyst{3,:});



%%  distribution des E fitter sur sigma v epsi
% [~,edges] = histcounts(Eplot2);
% 
% figure(1)
% hold on
% histogram(Eplot2,'facecolor',color,'binwidth',3)
% plot([median(Eplot2) median(Eplot2)],[0 45],'r--')
% ax = gca;
% % ax.XScale = 'log';
% % ax.XTick = [0.5 1 3 10 15 20 30 40];
% xlim([0 25])
% % ylim([0 15])
% xlabel('E (kPa)  - 2s comp')
% ylabel('Count')
% title(['Elastic modulus for 2s compression - ' num2str(median(Eplot2)) 'kPa'])
% 
% legend('Data','Median')
% 
% fig = gca;
% saveas(fig,[sff filesep 'EDist2s.png'],'png')
% saveas(fig,[sff filesep 'EDist2s.fig'],'fig')
% 
% 
% 
% [~,edges] = histcounts(Eplot4);
% 
% figure(2)
% hold on
% histogram(Eplot4,'facecolor',color,'binwidth',3)
% plot([median(Eplot4) median(Eplot4)],[0 45],'r--')
% ax = gca;
% % ax.XScale = 'log';
% % ax.XTick = [0.5 1 3 10 15 20 30 40];
% xlim([0 25])
% % ylim([0 15])
% xlabel('E (kPa) - 4s comp')
% ylabel('Count')
% title(['Elastic modulus for 4s compression - ' num2str(median(Eplot4)) 'kPa'])
% 
% legend('Data','Median')
% 
% fig = gca;
% saveas(fig,[sff filesep 'EDist4s.png'],'png')
% saveas(fig,[sff filesep 'EDist4s.fig'],'fig')
% 
% 
% 
% [~,edges] = histcounts(Eplot8);
% 
% figure(3)
% hold on
% histogram(Eplot8,'facecolor',color,'binwidth',3)
% plot([median(Eplot8) median(Eplot8)],[0 45],'r--')
% ax = gca;
% % ax.XScale = 'log';
% % ax.XTick = [0.5 1 3 10 15 20 30 40];
% xlim([0 25])
% % ylim([0 15])
% xlabel('E (kPa)  - 8s comp')
% ylabel('Count')
% title(['Elastic modulus for 8s compression - ' num2str(median(Eplot8)) 'kPa'])
% 
% legend('Data','Median')
% 
% fig = gca;
% saveas(fig,[sff filesep 'EDist8s.png'],'png')
% saveas(fig,[sff filesep 'EDist8s.fig'],'fig')
% 
% 
% figure(4)
% hold on
% 
% Spread_Julien(Eplot2,1,0.9,color,'o',1)
% plot([0.57 1.45],[mean(Eplot2) mean(Eplot2)],'r--','linewidth',2,'handlevisibility','off')
% Spread_Julien(Eplot4,2,0.9,'r','o',0.5)
% plot([1.57 2.45],[mean(Eplot4) mean(Eplot4)],'r--','linewidth',2,'handlevisibility','off')
% Spread_Julien(Eplot8,3,0.9,'g','o',0.5)
% plot([2.57 3.45],[mean(Eplot8) mean(Eplot8)],'r--','linewidth',2)
% 
% 
% %stats
% DQ = DoParamStats(Eplot2,Eplot4);
% DH = DoParamStats(Eplot2,Eplot8);
% QH = DoParamStats(Eplot4,Eplot8);
% 
% plot([1 2],[15 15],'k-','linewidth',1.5)
% text(1.5, 16, DQ,'HorizontalAlignment','center','fontsize',11)
% plot([1 3],[17 17],'k-','linewidth',1.5)
% text(2, 18, DH,'HorizontalAlignment','center','fontsize',11)
% plot([2 3],[19 19],'k-','linewidth',1.5)
% text(2.5, 20, QH,'HorizontalAlignment','center','fontsize',11)
% 
% ylabel('E (kPa)')
% 
% ax = gca;
% ax.XTick = [1 2 3];
% ax.XTickLabel = {'2s Comp','4s Comp','8s Comp'}
% 
% legend('Mean')
% 
% 
% fig = gca;
% saveas(fig,[sff filesep 'EDist248.png'],'png')
% saveas(fig,[sff filesep 'EDist248.fig'],'fig')

%%  Sigma final du fit
% [~,edges] = histcounts(FinalSigplot2);
% 
% figure(11)
% hold on
% histogram(FinalSigplot2,'facecolor',color,'binwidth',300)
% plot([median(FinalSigplot2) median(FinalSigplot2)],[0 45],'r--')
% ax = gca;
% % ax.XScale = 'log';
% % ax.XTick = [0.5 1 3 10 15 20 30 40];
% xlim([0 2500])
% ylim([0 25])
% xlabel('FinalSig (Pa)  - 2s comp')
% ylabel('Count')
% title(['FinalSigfor 2s compression - ' num2str(median(FinalSigplot2)) 'Pa'])
% 
% legend('Data','Median')
% 
% fig = gca;
% saveas(fig,[sff filesep 'FinalSigDist2s.png'],'png')
% saveas(fig,[sff filesep 'FinalSigDist2s.fig'],'fig')
% 
% 
% 
% [~,edges] = histcounts(FinalSigplot4);
% 
% figure(12)
% hold on
% histogram(FinalSigplot4,'facecolor',color,'binwidth',300)
% plot([median(FinalSigplot4) median(FinalSigplot4)],[0 45],'r--')
% ax = gca;
% % ax.XScale = 'log';
% % ax.XTick = [0.5 1 3 10 15 20 30 40];
% xlim([0 4500])
% ylim([0 25])
% xlabel('FinalSig (Pa) - 4s comp')
% ylabel('Count')
% title(['FinalSig for 4s compression - ' num2str(median(FinalSigplot4)) 'Pa'])
% 
% legend('Data','Median')
% 
% fig = gca;
% saveas(fig,[sff filesep 'FinalSigDist4s.png'],'png')
% saveas(fig,[sff filesep 'FinalSigDist4s.fig'],'fig')
% 
% 
% 
% [~,edges] = histcounts(FinalSigplot8);
% 
% figure(13)
% hold on
% histogram(FinalSigplot8,'facecolor',color,'binwidth',300)
% plot([median(FinalSigplot8) median(FinalSigplot8)],[0 45],'r--')
% ax = gca;
% % ax.XScale = 'log';
% % ax.XTick = [0.5 1 3 10 15 20 30 40];
% xlim([0 2500])
% ylim([0 25])
% xlabel('FinalSig (Pa)  - 8s comp')
% ylabel('Count')
% title(['FinalSig for 8s compression - ' num2str(median(FinalSigplot8)) 'Pa'])
% 
% legend('Data','Median')
% 
% fig = gca;
% saveas(fig,[sff filesep 'FinalSigDist8s.png'],'png')
% saveas(fig,[sff filesep 'FinalSigDist8s.fig'],'fig')
% 
% 
% figure(14)
% hold on
% 
% Spread_Julien(FinalSigplot2,1,0.9,color,'o',200)
% plot([0.57 1.45],[mean(FinalSigplot2) mean(FinalSigplot2)],'r--','linewidth',2,'handlevisibility','off')
% Spread_Julien(FinalSigplot4,2,0.9,'r','o',200)
% plot([1.57 2.45],[mean(FinalSigplot4) mean(FinalSigplot4)],'r--','linewidth',2,'handlevisibility','off')
% Spread_Julien(FinalSigplot8,3,0.9,'g','o',200)
% plot([2.57 3.45],[mean(FinalSigplot8) mean(FinalSigplot8)],'r--','linewidth',2)
% 
% 
% %stats
% DQ = DoParamStats(FinalSigplot2,FinalSigplot4);
% DH = DoParamStats(FinalSigplot2,FinalSigplot8);
% QH = DoParamStats(FinalSigplot4,FinalSigplot8);
% 
% plot([1 2],[4015 4015],'k-','linewidth',1.5)
% text(1.5, 4160, DQ,'HorizontalAlignment','center','fontsize',11)
% plot([1 3],[4250 4250],'k-','linewidth',1.5)
% text(2, 4400, DH,'HorizontalAlignment','center','fontsize',11)
% plot([2 3],[4650 4650],'k-','linewidth',1.5)
% text(2.5, 4800, QH,'HorizontalAlignment','center','fontsize',11)
% 
% ylabel('FinalSig (kPa)')
% 
% ax = gca;
% ax.XTick = [1 2 3];
% ax.XTickLabel = {'2s Comp','4s Comp','8s Comp'}
% 
% legend('Mean')
% 
% 
% fig = gca;
% saveas(fig,[sff filesep 'FinalSigDist248.png'],'png')
% saveas(fig,[sff filesep 'FinalSigDist248.fig'],'fig')



%% E vs SigFinal

figure(21)
hold on
plot(FinalSigplot2,Eplot2,'o')
plot(FinalSigplot4,Eplot4,'o')
plot(FinalSigplot8,Eplot8,'o')
xlabel('Sigma max of the lin fit (Pa)')
ylabel('E from lin fit (kPa)')



%% E local 

figure
hold on
histogram(E16mTplot2,'facecolor',color,'binwidth',3)
plot([median(E16mTplot2) median(E16mTplot2)],[0 45],'r--')
ax = gca;
% ax.XScale = 'log';
% ax.XTick = [0.5 1 3 10 15 20 30 40];
xlim([0 35])
% ylim([0 15])
xlabel('E (kPa)  - 2s comp')
ylabel('Count')
title(['Local elastic modulus for 2s compression - ' num2str(median(E16mTplot2)) 'kPa'])

legend('Data','Median')

fig = gca;
saveas(fig,[sff filesep 'E31mTDist2s.png'],'png')
saveas(fig,[sff filesep 'E31mTDist2s.fig'],'fig')




figure
hold on
histogram(E16mTplot4,'facecolor',color,'binwidth',3)
plot([median(E16mTplot4) median(E16mTplot4)],[0 45],'r--')
ax = gca;
% ax.XScale = 'log';
% ax.XTick = [0.5 1 3 10 15 20 30 40];
xlim([0 35])
% ylim([0 15])
xlabel('E (kPa) - 4s comp')
ylabel('Count')
title(['Local elastic modulus for 4s compression - ' num2str(median(E16mTplot4)) 'kPa'])

legend('Data','Median')

fig = gca;
saveas(fig,[sff filesep 'E31mTDist4s.png'],'png')
saveas(fig,[sff filesep 'E31mTDist4s.fig'],'fig')




figure
hold on
histogram(E16mTplot8,'facecolor',color,'binwidth',3)
plot([median(E16mTplot8) median(E16mTplot8)],[0 45],'r--')
ax = gca;
% ax.XScale = 'log';
% ax.XTick = [0.5 1 3 10 15 20 30 40];
xlim([0 25])
% ylim([0 15])
xlabel('E (kPa)  - 8s comp')
ylabel('Count')
title(['Local elastic modulus for 8s compression - ' num2str(median(E16mTplot8)) 'kPa'])

legend('Data','Median')

fig = gca;
saveas(fig,[sff filesep 'E31mTDist8s.png'],'png')
saveas(fig,[sff filesep 'E31mTDist8s.fig'],'fig')


figure
hold on

Spread_Julien(E16mTplot2,1,0.9,color,'o',1)
plot([0.57 1.45],[median(E16mTplot2) median(E16mTplot2)],'r--','linewidth',2,'handlevisibility','off')
Spread_Julien(E16mTplot4,2,0.9,'r','o',0.5)
plot([1.57 2.45],[median(E16mTplot4) median(E16mTplot4)],'r--','linewidth',2,'handlevisibility','off')
Spread_Julien(E16mTplot8,3,0.9,'g','o',0.5)
plot([2.57 3.45],[median(E16mTplot8) median(E16mTplot8)],'r--','linewidth',2)


%stats
DQ = DoParamStats(E16mTplot2,E16mTplot4);
DH = DoParamStats(E16mTplot2,E16mTplot8);
QH = DoParamStats(E16mTplot4,E16mTplot8);

plot([1 2],[15 15],'k-','linewidth',1.5)
text(1.5, 16, DQ,'HorizontalAlignment','center','fontsize',11)
plot([1 3],[17 17],'k-','linewidth',1.5)
text(2, 18, DH,'HorizontalAlignment','center','fontsize',11)
plot([2 3],[19 19],'k-','linewidth',1.5)
text(2.5, 20, QH,'HorizontalAlignment','center','fontsize',11)

ylabel('E local (kPa)')

ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'2s Comp','4s Comp','8s Comp'}

legend('Median')


fig = gca;
saveas(fig,[sff filesep 'E31mTDist248.png'],'png')
saveas(fig,[sff filesep 'E31mTDist248.fig'],'fig')


%%  distribution des hysteresis
% [~,edges] = histcounts(Hyst2);
% 
% figure(11)
% hold on
% histogram(Hyst2,edges,'facecolor',color)
% plot([mean(Hyst2) mean(Hyst2)],[0 45],'r--')
% ax = gca;
% 
% % xlim([0 100])
% 
% xlabel('Hysteresis (%)  - 2s comp')
% ylabel('Count')
% title(['Hysteresis (%) for 2s compression - ' num2str(mean(Hyst2)) '%'])
% 
% legend('Data','Mean')
% 
% fig = gca;
% saveas(fig,[sff filesep 'Hyst2s.png'],'png')
% saveas(fig,[sff filesep 'Hyst2s.fig'],'fig')
% 
% 
% 
% [~,edges] = histcounts(Hyst4);
% 
% figure(12)
% hold on
% histogram(Hyst4,edges,'facecolor',color)
% plot([mean(Hyst4) mean(Hyst4)],[0 45],'r--')
% ax = gca;
% 
% % xlim([0 100])
% 
% xlabel('Hysteresis (%) - 4s comp')
% ylabel('Count')
% title(['Hysteresis (%) for 4s compression - ' num2str(mean(Hyst4)) '%'])
% 
% legend('Data','Mean')
% 
% fig = gca;
% saveas(fig,[sff filesep 'Hyst4s.png'],'png')
% saveas(fig,[sff filesep 'Hyst4s.fig'],'fig')
% 
% 
% 
% [~,edges] = histcounts(Hyst8);
% 
% figure(13)
% hold on
% histogram(Hyst8,edges,'facecolor',color)
% plot([mean(Hyst8) mean(Hyst8)],[0 45],'r--')
% ax = gca;
% % xlim([-50 100])
% 
% xlabel('Hysteresis (%)  - 8s comp')
% ylabel('Count')
% title(['Hysteresis (%) for 8s compression - ' num2str(mean(Hyst8)) 'kPa'])
% 
% legend('Data','Mean')
% 
% fig = gca;
% saveas(fig,[sff filesep 'Hyst8s.png'],'png')
% saveas(fig,[sff filesep 'Hyst8s.fig'],'fig')
% 
% 
% figure(14)
% hold on
% 
% Spread_Julien(Hyst2,1,0.9,color,'o',20)
% plot([0.57 1.45],[mean(Hyst2) mean(Hyst2)],'r--','linewidth',2,'handlevisibility','off')
% Spread_Julien(Hyst4,2,0.9,'r','o',20)
% plot([1.57 2.45],[mean(Hyst4) mean(Hyst4)],'r--','linewidth',2,'handlevisibility','off')
% Spread_Julien(Hyst8,3,0.9,'g','o',20)
% plot([2.57 3.45],[mean(Hyst8) mean(Hyst8)],'r--','linewidth',2)
% 
% 
% %stats
% DQ = DoParamStats(Hyst2,Hyst4);
% DH = DoParamStats(Hyst2,Hyst8);
% QH = DoParamStats(Hyst4,Hyst8);
% 
% plot([1 2],[100 100],'k-','linewidth',1.5)
% text(1.5, 110, DQ,'HorizontalAlignment','center','fontsize',11)
% plot([1 3],[120 120],'k-','linewidth',1.5)
% text(2, 130, DH,'HorizontalAlignment','center','fontsize',11)
% plot([2 3],[140 140],'k-','linewidth',1.5)
% text(2.5, 150, QH,'HorizontalAlignment','center','fontsize',11)
% 
% ylabel('Hysteresis (%)')
% 
% ax = gca;
% ax.XTick = [1 2 3];
% ax.XTickLabel = {'2s Comp','4s Comp','8s Comp'}
% 
% legend('Mean')
% 
% 
% fig = gca;
% saveas(fig,[sff filesep 'Hyst248.png'],'png')
% saveas(fig,[sff filesep 'Hyst248.fig'],'fig')

end


%% extra functions
function txt = DoParamStats(A,B)

pCM = ranksum(A,B);

if 5/100 > pCM && pCM > 1/100
    txt = ['*']  ;
elseif 1/100 > pCM && pCM > 1/1000
    txt = ['**'];
elseif 1/1000 > pCM
    txt = ['***'];
else
    txt = ['NS'];
    
end
end
