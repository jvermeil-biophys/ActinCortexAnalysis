%% INIT
clear all
close all

warning('off','all')

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Serif')
set(0,'DefaultAxesFontWeight','normal')


% Condition colors
Cdc = [0 109 219]./255; % DC

Cdic = [73 180 50]./255; %  Dicty

Cfib = [176 25 59]./255; % Fibro



%% Get Data
DataFolder = 'C:\Users\Valentin\Desktop\DataValentin\MatFile\P2S';
cd(DataFolder)

load('./P2S_Dicty_WT')
Dicty_MedC = MS.MedCortWs - 11; % Cortices medians for Dictys// 11 = coating
Dicty_Acti = MS.CortActis;

load('./P2S_DCLA_DMSO.mat')
DC_MedC = MS.MedCortWs - 11; % For DCs // 11 = coating
DC_Acti = MS.CortActis;

load('./P2S_3T3_PLL.mat')
F3T3_MedC = MS.MedCortWs - 11; % For Fibroblasts // 11 = coating
F3T3_Acti = MS.CortActis;


%% Graph 1 : Median thicknesses
figure
hold on


data = {DC_MedC, Dicty_MedC};

set(0,'DefaultLineMarkerSize',9);
h = plotSpread_V(data,'distributionMarkers',{'o','o'},'distributionColors',{Cdc,Cdic} ...
,'xValues',[1 2],'xNames',{'DCs' 'Dictys'},'spreadWidth',1.4)

ylabel({'Median cortex thickness (nm)'})
h{3}.XTickLabelRotation = 45;

% data
% Spread_Julien(DC_MedC,1,0.9,Cdc,'o',70)
plot([0.55 1.45],[median(DC_MedC) median(DC_MedC)],'r--','linewidth',2)
% Spread_Julien(Dicty_MedC,2,0.9,Cdic,'o',70)
plot([1.55 2.45],[median(Dicty_MedC) median(Dicty_MedC)],'r--','linewidth',2)
% Spread_Julien(F3T3_MedC,3,0.9,Cfib,'o',70)
% plot([2.55 3.45],[median(F3T3_MedC) median(F3T3_MedC)],'r--','linewidth',2)

%stats
DCDic = DoParamStats(DC_MedC,Dicty_MedC);
DCF = DoParamStats(DC_MedC,F3T3_MedC);
% FDic = DoParamStats(F3T3_MedC,Dicty_MedC);

plot([1 2],[800 800],'k-','linewidth',1.5)
text(1.5, 820, DCDic,'HorizontalAlignment','center','fontsize',11)
% plot([1 3],[850 850],'k-','linewidth',1.5)
% text(2, 870, DCF,'HorizontalAlignment','center','fontsize',11)
% plot([2 3],[920 920],'k-','linewidth',1.5)
% text(2.5, 940, FDic,'HorizontalAlignment','center','fontsize',11)

% ax1 = gca;
% ax1.XTick = [1 2 3];
% ax1.XTickLabel = {'DCs' 'Dictys' '3T3'};
% ax1.XTickLabelRotation = 45;

ax = gca;
ax.LineWidth = 1.5;
box on
ylim([0 1200])
xlim([0.45 2.55])


set(0,'DefaultAxesFontSize',18)

%% Graph 2 : Amplitude of fluctuations
figure
hold on


data = {DC_Acti, Dicty_Acti};

set(0,'DefaultLineMarkerSize',9);
h = plotSpread_V(data,'distributionMarkers',{'o','o'},'distributionColors',{Cdc,Cdic} ...
,'xValues',[1 2],'xNames',{'DCs' 'Dictys'},'spreadWidth',1.4)

h{3}.XTickLabelRotation = 45;

ylabel({'Amplitude Fluctuations (nm)'})
% data
% Spread_Julien(DC_Acti,1,0.9,Cdc,'o',70)
plot([0.55 1.45],[median(DC_Acti) median(DC_Acti)],'r--','linewidth',2)
% Spread_Julien(Dicty_Acti,2,0.9,Cdic,'o',70)
plot([1.55 2.45],[median(Dicty_Acti) median(Dicty_Acti)],'r--','linewidth',2)
% Spread_Julien(F3T3_Acti,3,0.9,Cfib,'o',70)
% plot([2.55 3.45],[median(F3T3_Acti) median(F3T3_Acti)],'r--','linewidth',2)

%stats
DCDic = DoParamStats(DC_Acti,Dicty_Acti);
% DCF = DoParamStats(DC_Acti,F3T3_Acti);
% FDic = DoParamStats(F3T3_Acti,Dicty_Acti);

plot([1 2],[800 800],'k-','linewidth',1.5)
text(1.5, 820, DCDic,'HorizontalAlignment','center','fontsize',11)
% plot([1 3],[850 850],'k-','linewidth',1.5)
% text(2, 870, DCF,'HorizontalAlignment','center','fontsize',11)
% plot([2 3],[920 920],'k-','linewidth',1.5)
% text(2.5, 940, FDic,'HorizontalAlignment','center','fontsize',11)

% ax2 = gca;
% ax2.XTick = [1 2 3];
% ax2.XTickLabel = {'DCs' 'Dictys' '3T3'};
% ax2.XTickLabelRotation = 45;

ax = gca;
ax.LineWidth = 1.5;
box on
ylim([0 1200])
xlim([0.45 2.55])


%% Graph 3 : correlation
figure(3)
subplot(121)
hold on
ylabel('Fluctuations amplitude (nm)')
xlabel('Cortex thickness (nm)')
plot(DC_MedC,DC_Acti,'o','color',Cdc,'markerfacecolor',Cdc)
ax = gca;
h =lsline(ax);
ylim([-200 2000])
xlim([0 1000])

corrmat = corrcoef(DC_MedC,DC_Acti);
slope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));

DCfit = fitlm(DC_MedC,DC_Acti);
DCci = coefCI(DCfit);

SlopeCi = mean(abs(slope-DCci(2,:)));

txt = {['Slope = ' num2str(slope,'%.2f') ' $\pm$ ' num2str(SlopeCi,'%.2f')]};

text(20,1800,txt,'Interpreter','Latex','FontSize',20)

ax = gca;
ax.LineWidth = 1.5;
box on
clear SlopeCi slope corrmat txt ax h

subplot(122)
hold on
ylabel('Fluctuations amplitude (nm)')
xlabel('Cortex thickness (nm)')
plot(Dicty_MedC,Dicty_Acti,'o','color',Cdic,'markerfacecolor',Cdic)
ax = gca;
h =lsline(ax);
ylim([-200 2000])
xlim([0 1000])

corrmat = corrcoef(Dicty_MedC,Dicty_Acti);
slope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));

Dictyfit = fitlm(Dicty_MedC,Dicty_Acti);
Dictyci = coefCI(Dictyfit);

SlopeCi = mean(abs(slope-Dictyci(2,:)));

txt = {['Slope = ' num2str(slope,'%.2f') ' $\pm$ ' num2str(SlopeCi,'%.2f')]};

ax = gca;
ax.LineWidth = 1.5;
box on
text(20,1800,txt,'Interpreter','Latex','FontSize',20)

clear SlopeCi slope corrmat txt

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


