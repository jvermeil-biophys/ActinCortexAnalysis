function Meca2Plot(specif,color,resfolder,figfolder)

%% INIT

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextInterpreter','Latex')

here = pwd;

%% params and settings

% dossier de donn�es et sauvegarde des donn�es
path = [resfolder filesep 'D2M'];
datenow = datestr(now,'yy-mm-dd');

% sauvergarde des figures
sfig = figfolder;
sff = [sfig filesep datenow filesep 'Meca' filesep specif];
mkdir(sff)

CortThick = [];
E = {};
E0dim = {};
E0chad = {};
EchadCI = {};
MDef = {};
mThick = {};
Delta5 = {};
Delta20 = {};
Alpha = {};
H0DIM = {};
H0ABS= {};
H0FIT= {};
CompN = {};
Hyst = {};


%% Retrieve datafile
List = ls(path);

for i = 3:size(List,1)
    if contains(List(i,:),specif)
        
        fprintf(['Chargement du fichier ' List(i,:) '\n\n'])
        load([path filesep List(i,:)]);
        
        
        CortThick = [CortThick Cort];
        E = {E{:} E0{:}};
        E0dim = {E0dim{:} Edim{:}};
        E0chad = {E0chad{:} Echad{:}};
        EchadCI = {EchadCI{:} EchadConf{:}};
        MDef = {MDef{:} MaxDef{:}};
        mThick = {mThick{:} MinThick{:}};
        Delta5 = {Delta5{:} Delt5{:}};
        Delta20 = {Delta20{:} Delt20{:}};
        Alpha = {Alpha{:} AlphaNL{:}};
        H0DIM = {H0DIM{:} DimH0{:}};
        H0ABS= {H0ABS{:} H0abs{:}};
        H0FIT = {H0FIT{:} ChadsH0{:}};
        CompN = {CompN{:} CompNum{:}};
%         Hyst = {Hyst{:} Hys{:}};
        
    end
end

Eplot = vertcat(E{:})/1000;

Edimplot = vertcat(E0dim{:})/1000;
Echadplot = vertcat(E0chad{:})/1000;
EchadCIvect = vertcat(EchadCI{:})/1000;

Delta5plot = vertcat(Delta5{:});
Delta20plot = vertcat(Delta20{:});

MaxDefplot = vertcat(MDef{:});

Alphaplot = vertcat(Alpha{:});

length(Eplot)

CompNs = vertcat(CompN{:});

H0ABSplot = vertcat(H0ABS{:});
H0DIMplot = vertcat(H0DIM{:});
H0FITplot = vertcat(H0FIT{:});

MinThickplot = vertcat(mThick{:});


for  k = 1:length(E)
    if length(E{k})>1
        Ecell(k) = nanmean(E{k}/1000);
        ErrEcell(k) = std(E{k}/1000)/sqrt(length(E{k}));
        nE(k) = length(E{k});
    end
end

for  k = 1:length(E0dim)
    if length(E0dim{k})>1
        E0dimcell(k) = nanmean(E0dim{k}/1000);
        ErrE0dimcell(k) = std(E0dim{k}/1000)/sqrt(length(E0dim{k}));
        nE0dim(k) = length(E0dim{k});
    end
end

for  k = 1:length(E0chad)
    if length(E0chad{k})>1
        E0chadcell(k) = nanmean(E0chad{k}/1000);
        ErrE0chadcell(k) = std(E0chad{k}/1000)/sqrt(length(E0chad{k}));
        nE0chad(k) = length(E0chad{k});
    end
end


EplotNorm = [];

for k = 1:length(Ecell)
    
    EplotNorm = [EplotNorm; E{k}/(1000*Ecell(k))];

end


Alphas = vertcat(Alpha{:});

ptrnan = ~isnan(Alphas);


%% Distrib Echad 25% avec moyenne pond�r�e


EchadWeight = (Echadplot./EchadCIvect).^2;

EchadWeightedMean = sum(Echadplot.*EchadWeight)/sum(EchadWeight);

EchadWeightedStd = sqrt(sum(((Echadplot-EchadWeightedMean).^2.*EchadWeight))/sum(EchadWeight));

[~,edges] = histcounts(log10(Echadplot));

ledges = log(edges);

figure
hold on
% title(['Comp on DC. Weighter average : ' num2str(EchadWeightedMean) ' kPa'])
title(['Weighter average : ' num2str(EchadWeightedMean) '+- ' num2str(EchadWeightedStd) ' kPa'])
histogram(Echadplot,10.^edges,'facecolor',color)

plot([EchadWeightedMean EchadWeightedMean],[0 120],'r--','linewidth',1.7)
ax = gca;
ax.LineWidth = 1.5;
ax.XScale = 'log';
ax.XTick = [1 10 100];
box on


xlim([0.5 200])
% legend({'E from Chadwick fit','Weighted mean'})
xlabel('Elastic modulus (kPa)')
ylabel('Count')
% fig = gca;
% saveas(fig,[sff filesep 'EchadDist.png'],'png')
% saveas(fig,[sff filesep 'EchadDist.fig'],'fig')

ylim([0 35])

figure
hold on
xlabel('Echad (kPa)')
ylabel('Conf Int / value')
plot(Echadplot,EchadCIvect./Echadplot,'ok','markerfacecolor',color)

%% Plot Echad vs H0Chad et vs Hdepart

% % H0chad
% 
% rfit = robustfit(log(H0FITplot),log(Echadplot));
% 
% figure
% hold on
% title(['Power coeff : ' num2str(rfit(2))])
% plot(H0FITplot,Echadplot,'ok','markerfacecolor',color)
% plot(sort(H0FITplot),exp(rfit(1))*sort(H0FITplot).^rfit(2),'r-','linewidth',1.5)
% 
% xlabel('Thickness from Chadwick fit (nm)')
% ylabel('E from Chadwick fit (kPa)')
% 
% legend('Data','Robust log-log fit')
% 
% legend('location','best')
% 
% ax = gca;
% ax.XScale = 'log';
% ax.YScale = 'log';
% 
% fig = gcf;
% saveas(fig,[sff filesep 'EchadVH0chad.png'],'png')
% saveas(fig,[sff filesep 'EchadVH0chad.fig'],'fig')

% Hdepart


[rfit,stats] = robustfit(log(H0ABSplot),log(Echadplot));
[CCf,P] = corrcoef(H0ABSplot,Echadplot);

figure
hold on
% title(['Power coeff se : ' num2str(stats.se(2)) ' - Pval Corrcoeff : ' num2str(P(1,2))])
plot(H0ABSplot,Echadplot,'ok','markerfacecolor',color)
% plot(sort(H0ABSplot),exp(rfit(1))*sort(H0ABSplot).^rfit(2),'r-','linewidth',1.5)


xlabel('Thickness at the begining of compression (nm)')
ylabel('E from Chadwick fit (kPa)')

% legend('Data','Robust log-log fit')

% legend('location','best')

% ax = gca;
% ax.XScale = 'log';
% ax.YScale = 'log';

fig = gcf;
saveas(fig,[sff filesep 'EchadVH0begining.png'],'png')
saveas(fig,[sff filesep 'EchadVH0begining.fig'],'fig')

%% Plot Echad vs MinThick

% 
% figure
% hold on
% title('2s Compression on dendritic cells')
% plot(MinThickplot*1000,Echadplot,'ok','markerfacecolor',color)
% 
% xlabel('Thickness at max compression (nm)')
% ylabel('E from Chadwick fit (kPa)')
% 
% fig = gca;
% saveas(fig,[sff filesep 'EchadVminThick.png'],'png')
% saveas(fig,[sff filesep 'EchadVminThick.fig'],'fig')

%% Plot Edim vs H0dim et vs Hdepart
% 
% % H0chad
% 
% figure(4)
% hold on
% % title('2s Compression on dendritic cells')
% plot(H0DIMplot,Edimplot,'ok','markerfacecolor',color)
% 
% xlabel('Thickness from Dimitriadis fit (nm)')
% ylabel('E from Dimitriadis fit (kPa)')
% 
% fig = gca;
% saveas(fig,[sff filesep 'EdimVH0chad.png'],'png')
% saveas(fig,[sff filesep 'EdimVH0chad.fig'],'fig')
% 
% % Hdepart
% 
% figure(5)
% hold on
% % title('2s Compression on dendritic cells')
% plot(H0ABSplot,Edimplot,'ok','markerfacecolor',color)
% 
% xlabel('Thickness at the begining of compression (nm)')
% ylabel('E from Dimitriadis fit (kPa)')
% 
% 
% fig = gca;
% saveas(fig,[sff filesep 'EdimVH0begining.png'],'png')
% saveas(fig,[sff filesep 'EcdimVH0begining.fig'],'fig')
% 
% 

%% Distrib E a nous, Edim, relation E, Edim, Echad
% 
% figure(1000)
% hold on
% plot([0 25],[0 25],'k--')
% plot(Eplot,Edimplot,'ok','markerfacecolor',color)
% xlabel('E linsmalldef (kPa)')
% ylabel('E dimitriadis (kPa)')
% 
% 
% figure(2000)
% hold on
% plot([0 25],[0 25],'k--')
% plot(Eplot,Echadplot,'ok','markerfacecolor',color)
% xlabel('E linsmalldef (kPa)')
% ylabel('E chadwick (kPa)')
% 
% 
% [~,edges] = histcounts(Eplot,'binwidth',2);
% 
% 
% 
% figure(1)
% hold on
% histogram(Eplot,edges,'facecolor',color)
% plot([mean(Eplot) mean(Eplot)],[0 45],'r--')
% ax = gca;
% % ax.XScale = 'log';
% % ax.XTick = [0.5 1 3 10 15 20 30 40];
% xlim([0 25])
% % ylim([0 15])
% xlabel('E (kPa)')
% ylabel('Count')
% fig = gca;
% saveas(fig,[sff filesep 'EDist.png'],'png')
% saveas(fig,[sff filesep 'EDist.fig'],'fig')
% 
% 
% [~,edges] = histcounts(Edimplot,'binwidth',2);
% 
% 
% 
% figure(1001)
% hold on
% histogram(Edimplot,edges,'facecolor',color)
% plot([median(Edimplot) median(Edimplot)],[0 26],'r--')
% ax = gca;
% % ax.XTick = [0.01 0.05 0.1 0.5 1 3 10 15];
% % xlim([0.05 15])
% % ylim([0 25])
% xlabel('Edim (kPa)')
% ylabel('Count')
% fig = gca;
% saveas(fig,[sff filesep 'EdimDist.png'],'png')
% saveas(fig,[sff filesep 'EdimDist.fig'],'fig')
% 


% figure(3000)
% hold on
% plot(Eplot,Echadplot,'ko','markerfacecolor',color)
% xlabel('E a nous')
% ylabel('E chadwick')

%% Enfoncement a 5mT
% 
% [~,edges] = histcounts(Delta5plot,'binwidth',2);
% 
% 
% 
% figure(2002)
% hold on
% histogram(Delta5plot,edges,'facecolor',color)
% plot([nanmedian(Delta5plot) nanmedian(Delta5plot)],[0 10],'r--')
% ax = gca;
% xlim([0 160])
% xlabel('Delta 5mT (nm)')
% ylabel('Count')
% fig = gca;
% saveas(fig,[sff filesep 'Delta5Dist.png'],'png')
% saveas(fig,[sff filesep 'Delta5Dist.fig'],'fig')
% 
% figure(2003)
% hold on
% plot(H0FITplot,Delta5plot./H0FITplot,'ok','markerfacecolor',color)
% plot([100 1100],[nanmedian(Delta5plot./H0FITplot) nanmean(Delta5plot./H0FITplot)],'r--')
% xlabel('Epaisseur de d�part fitter (nm)')
% ylabel('Delta 5mT / H0 (nm)')
% 
% nanmedian(Delta5plot./H0FITplot)

% figure(2004)
% hold on
% plot(Echadplot,Delta5plot,'ok','markerfacecolor',color)
% xlabel('E chadwick (kPa)')
% ylabel('Delta 5mT (nm)')

%% Enfoncement a 20pN
% 
% [~,edges] = histcounts(Delta20plot,'binwidth',2);
% 
% 
% 
% figure(3002)
% hold on
% histogram(Delta20plot,edges,'facecolor',color)
% plot([nanmedian(Delta20plot) nanmedian(Delta20plot)],[0 10],'r--')
% ax = gca;
% xlim([0 160])
% xlabel('Delta 20pN (nm)')
% ylabel('Count')
% fig = gca;
% saveas(fig,[sff filesep 'Delta20Dist.png'],'png')
% saveas(fig,[sff filesep 'Delta20Dist.fig'],'fig')
% 
% figure(3003)
% hold on
% plot(H0FITplot,Delta20plot./H0FITplot,'ok','markerfacecolor',color)
% plot([100 1100],[nanmean(Delta20plot./H0FITplot) nanmean(Delta20plot./H0FITplot)],'r--')
% xlabel('Epaisseur de d�part fitter (nm)')
% ylabel('Delta 20pN / H0 (nm)')
% 
% 
% nanmedian(Delta20plot./H0FITplot)
% 
% figure(3004)
% hold on
% plot(Echadplot,Delta20plot,'ok','markerfacecolor',color)
% xlabel('E chadwick (kPa)')
% ylabel('Delta 20pN (nm)')

%% deformation maximum
% 
% figure
% hold on
% histogram(MaxDefplot)
% plot([nanmean(MaxDefplot) nanmean(MaxDefplot)],[0 50],'r--')
% ylabel('Count')
% xlabel('Maximum relative deformation')
% 
% figure
% hold on
% plot(MaxDefplot,Alphaplot,'o')

%% Elasticit� moyenne par cell en fonction de l'�paisseur du cortex
% 
% figure(2)
% hold on
% errorbar(CortThick(Ecell~=0),Ecell(Ecell~=0),ErrEcell(Ecell~=0),'.k')
% plot(CortThick(Ecell~=0),Ecell(Ecell~=0),'ko','markerfacecolor',color)
% xlabel('Cortex thickness (nm)')
% ylabel('Average E (kPa) per cell')
% fig = gca;
% saveas(fig,[sff filesep 'EvCort.png'],'png')
% saveas(fig,[sff filesep 'EvCort.fig'],'fig')
% 
% 
% figure(102)
% hold on
% errorbar(CortThick(E0dimcell~=0),E0dimcell(E0dimcell~=0),ErrE0dimcell(E0dimcell~=0),'.k')
% plot(CortThick(E0dimcell~=0),E0dimcell(E0dimcell~=0),'ko','markerfacecolor',color)
% xlabel('Cortex thickness (nm)')
% ylabel('Average Edim (kPa) per cell')
% fig = gca;
% saveas(fig,[sff filesep 'EdimvCort.png'],'png')
% saveas(fig,[sff filesep 'EdimvCort.fig'],'fig')

%% Elasticit� lin�aire vs Hauteur de depart 
% % 
% [~,Pmat] = corrcoef(H0ABSplot,Eplot);
% 
% figure(3)
% hold on
% % scatter(H0ABSplot,Eplot,'ko','markerfacecolor',color)
% plot(H0ABSplot,Echadplot,'*','color',color,'linewidth',2)
% % lsline
% % legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% xlabel('Thickness before compression')
% ylabel('Elastic modulus - Chadwick (kPa)')
% ylim([0 30])
% fig = gca;
% saveas(fig,[sff filesep 'EchadvThick.png'],'png')
% saveas(fig,[sff filesep 'EchadvThick.fig'],'fig')
% % % 
% % H0log = log(H0ABSplot);
% Elog = log(Eplot);
% [~,Pmat] = corrcoef(H0log,Elog);
% P = polyfit(H0log,Elog,1);
% 
% ax.XScale = 'log'
% saveas(fig,[sff filesep 'EvThickLog.png'],'png')
% saveas(fig,[sff filesep 'EvThickLog.fig'],'fig')

% 
% figure(103)
% hold on
% plot(H0ABSplot,Edimplot,'ko','markerfacecolor',color)
% % lsline
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% xlabel('Thickness before compression')
% ylabel('Elastic modulus dimitriadis (kPa)')
% ylim([0 12])
% fig = gca;
% saveas(fig,[sff filesep 'EdimvThick.png'],'png')
% saveas(fig,[sff filesep 'EdimvThick.fig'],'fig')

%% Elasticit� lin�airre norm par cell vs Hauteur de depart en %

% [~,Pmat] = corrcoef(H0PRCplot,EplotNorm);
% 
% figure(4)
% hold on
% plot(H0PRCplot,EplotNorm,'ko','markerfacecolor',color)
% lsline
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% xlabel('Thickness amplitude percentage before compression')
% ylabel('E / $E_{cell}$')
% fig = gca;
% saveas(fig,[sff filesep 'EvPRC.png'],'png')
% saveas(fig,[sff filesep 'EvPRC.fig'],'fig')

%% E en fonction du n� de la comp


% [~,Pmat] = corrcoef(CompNs,Eplot);
% 
% figure(5)
% hold on
% plot(CompNs,Eplot,'ko','markerfacecolor',color)
% lsline
% xlabel('n$^{th}$ compression on the same cell')
% ylabel('E (kPa)')
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% fig = gca;
% saveas(fig,[sff filesep 'EvNComp.png'],'png')
% saveas(fig,[sff filesep 'EvNComp.fig'],'fig')

% %% Distribution des exposant non lineaire
% figure(301)
% hold on
% histogram(Alphas,'binwidth',0.2,'facecolor',color)
% plot([mean(Alphas) mean(Alphas)],[0 45],'r--')
% xlabel('Non-linear exponent')
% ylabel('Count')
% fig = gca;
% saveas(fig,[sff filesep 'AlphaDist.png'],'png')
% saveas(fig,[sff filesep 'AlphaDist.fig'],'fig')
% 
% xlim([0 2.5])
% 
% legend('Control','Latranculin A')

%% Exposant NL par cell en fonction de l'�paisseur du cortex

% for  k = 1:length(Alpha)
%     if ~isempty(Alpha{k})
%         Acell(k) = nanmean(Alpha{k});
%         ErrAcell(k) = nanstd(Alpha{k})/sqrt(length(Alpha{k}));
%     end
% end
% 
% 
% [~,Pmat] = corrcoef(CortThick(Acell~=0),Acell(Acell~=0));
% 
% figure(102)
% hold on
% errorbar(CortThick(Acell~=0),Acell(Acell~=0),ErrAcell((Acell~=0)),'.k','handlevisibility','off')
% plot(CortThick(Acell~=0),Acell(Acell~=0),'ko','markerfacecolor',color)
% lsline
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% xlabel('Cortex thickness (nm)')
% ylabel('Average Alpha per cell')
% fig = gca;
% saveas(fig,[sff filesep 'ALPHAvCort.png'],'png')
% saveas(fig,[sff filesep 'ALPHAvCort.fig'],'fig')

%% exposant NL vs Hauteur de depart en %

% 
% [~,Pmat] = corrcoef(H0ABSplot(ptrnan),Alphas(ptrnan));
% 
% figure(103)
% hold on
% plot(H0ABSplot,Alphas,'ko','markerfacecolor',color)
% plot([0 1000],[1 1],'--')
% lsline
% xlabel('Thickness before compression')
% ylabel('Non-linear exponent')
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% fig = gca;
% saveas(fig,[sff filesep 'ALPHAvPRC.png'],'png')
% saveas(fig,[sff filesep 'ALPHAvPRC.fig'],'fig')

%% Exposant NL en fonction du n� de la comp

% [~,Pmat] = corrcoef(CompNs(ptrnan),Alphas(ptrnan));
% 
% figure(104)
% hold on
% plot(CompNs,Alphas,'ko','markerfacecolor',color)
% lsline
% xlabel('n$^{th}$ compression on the same cell')
% ylabel('Non-linear exponent')
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% fig = gca;
% saveas(fig,[sff filesep 'AlphavNComp.png'],'png')
% saveas(fig,[sff filesep 'AlphavNComp.fig'],'fig')

%% module v exposant nl
% 
% [~,Pmat] = corrcoef(Eplot(ptrnan),Alphas(ptrnan));
% 
% figure(105)
% hold on
% plot(Alphas,Eplot,'ko','markerfacecolor',color)
% lsline
% xlabel('Non-linear exponent')
% ylabel('Elastic modulus (kPa)')
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% fig = gca;
% saveas(fig,[sff filesep 'AlphavE.png'],'png')
% saveas(fig,[sff filesep 'AlphavE.fig'],'fig')

%% Distribution des hysteresis
% Hyste = vertcat(Hyst{:});
% 
% ptrhys = ~isnan(Hyste)&Hyste>0;
% 
% Hyste = Hyste(ptrhys);
% 
% figure(201)
% hold on
% histogram(Hyste,'binwidth',10,'facecolor',color)
% plot([mean(Hyste) mean(Hyste)],[0 20],'r--')
% xlabel('Hysteresis (\%)')
% ylabel('Count')
% fig = gca;
% saveas(fig,[sff filesep 'HysDist.png'],'png')
% saveas(fig,[sff filesep 'HysDist.fig'],'fig')

%% Hysteresis par cell en fonction de l'�paisseur du cortex

% for  k = 1:length(Hyst)
%     if ~isempty(Hyst{k})
%         Hcell(k) = nanmean(Hyst{k}(Hyst{k}>0));
%         ErrHcell(k) = nanstd(Hyst{k}(Hyst{k}>0))/sqrt(length(Hyst{k}(Hyst{k}>0)));
%     end
% end
% 
% 
% [~,Pmat] = corrcoef(CortThick(Hcell~=0),Hcell(Hcell~=0));
% 
% figure(202)
% hold on
% errorbar(CortThick(Hcell~=0),Hcell(Acell~=0),ErrHcell((Hcell~=0)),'.k','handlevisibility','off')
% plot(CortThick(Hcell~=0),Hcell(Hcell~=0),'ko','markerfacecolor',color)
% lsline
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% xlabel('Cortex thickness (nm)')
% ylabel('Average hysteresis per cell')
% fig = gca;
% saveas(fig,[sff filesep 'HysvCort.png'],'png')
% saveas(fig,[sff filesep 'HysvCort.fig'],'fig')

%% Hysteresis vs Hauteur de depart en %

% H0PRCplot = vertcat(H0PRC{:});
% 
% [~,Pmat] = corrcoef(H0PRCplot(ptrhys),Hyste);
% 
% figure(203)
% hold on
% plot(H0PRCplot(ptrhys),Hyste,'ko','markerfacecolor',color)
% lsline
% xlabel('Thickness amplitude percentage before compression')
% ylabel('hysteresis (\%)')
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% fig = gca;
% saveas(fig,[sff filesep 'HysvPRC.png'],'png')
% saveas(fig,[sff filesep 'HysvPRC.fig'],'fig')

%% Hysteresis en fonction du n� de la comp
% 
% [~,Pmat] = corrcoef(CompNs(ptrhys),Hyste);
% 
% figure(204)
% hold on
% plot(CompNs(ptrhys),Hyste,'ko','markerfacecolor',color)
% lsline
% xlabel('n$^{th}$ compression on the same cell')
% ylabel('Hysteresis (\%)')
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% fig = gca;
% saveas(fig,[sff filesep 'HysvNComp.png'],'png')
% saveas(fig,[sff filesep 'HysvNComp.fig'],'fig')

%% module v hysteresis

% [~,Pmat] = corrcoef(Eplot(ptrhys),Hyste);
% 
% figure(205)
% hold on
% plot(Hyste,Eplot(ptrhys),'ko','markerfacecolor',color)
% lsline
% xlabel('Hysteresis (\%)')
% ylabel('Elastic modulus (kPa)')
% legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% fig = gca;
% saveas(fig,[sff filesep 'HysvE.png'],'png')
% saveas(fig,[sff filesep 'HysvE.fig'],'fig')

%% end
% close all

end
