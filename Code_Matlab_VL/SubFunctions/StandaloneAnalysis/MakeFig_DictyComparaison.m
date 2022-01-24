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
DataFolder = 'D:\Data\MatFile\P2S';
cd(DataFolder)

load('./P2S_Dicty_WT')
Ax3_MedC = MS.MedCortWs - 11; % Cortices medians for Dictys// 11 = coating
Ax3_Acti = MS.CortActis;

load('./P2S_DictyAx2_DMSO.mat')
Ax2_MedC = MS.MedCortWs - 11; % For DCs // 11 = coating
Ax2_Acti = MS.CortActis;

load('./P2S_DictyAx2_BSA30.mat')
Ax230_MedC = MS.MedCortWs - 11; % For Fibroblasts // 11 = coating
Ax230_Acti = MS.CortActis;


%% Graph 1 : Median thicknesses
figure
hold on


data = {Ax3_MedC, Ax2_MedC, Ax230_MedC};

set(0,'DefaultLineMarkerSize',9);
h = plotSpread_V(data,'distributionMarkers',{'o','o','o'},'distributionColors',{Cdc,Cdic, Cfib} ...
,'xValues',[1 2 3],'xNames',{'Ax3','Ax2','Ax230'},'spreadWidth',1.4)

ylabel({'Median cortex thickness (nm)'})
h{3}.XTickLabelRotation = 45;

% data

plot([0.55 1.45],[median(Ax2_MedC) median(Ax2_MedC)],'r--','linewidth',2)

plot([1.55 2.45],[median(Ax3_MedC) median(Ax3_MedC)],'r--','linewidth',2)

plot([2.55 3.45],[median(Ax230_MedC) median(Ax230_MedC)],'r--','linewidth',2)

%stats
DCDic = DoParamStats(Ax2_MedC,Ax3_MedC);
DCF = DoParamStats(Ax2_MedC,Ax230_MedC);
FDic = DoParamStats(Ax230_MedC,Ax3_MedC);

plot([1 2],[800 800],'k-','linewidth',1.5)
text(1.5, 820, DCDic,'HorizontalAlignment','center','fontsize',11)
plot([1 3],[850 850],'k-','linewidth',1.5)
text(2, 870, DCF,'HorizontalAlignment','center','fontsize',11)
plot([2 3],[920 920],'k-','linewidth',1.5)
text(2.5, 940, FDic,'HorizontalAlignment','center','fontsize',11)

% ax1 = gca;
% ax1.XTick = [1 2 3];
% ax1.XTickLabel = {'DCs' 'Dictys' '3T3'};
% ax1.XTickLabelRotation = 45;

ax = gca;
ax.LineWidth = 1.5;
box on
ylim([0 1200])
xlim([0.45 3.55])


set(0,'DefaultAxesFontSize',18)

%% Graph 2 : Amplitude of fluctuations
figure
hold on


data = { Ax3_Acti, Ax2_Acti, Ax230_Acti};



set(0,'DefaultLineMarkerSize',9);
h = plotSpread_V(data,'distributionMarkers',{'o','o','o'},'distributionColors',{Cdc,Cdic, Cfib} ...
,'xValues',[1 2 3],'xNames',{'Ax3','Ax2','Ax230'},'spreadWidth',1.4)

ylabel({'Median cortex thickness (nm)'})
h{3}.XTickLabelRotation = 45;

ylabel({'Amplitude Fluctuations (nm)'})
% data
% Spread_Julien(DC_Acti,1,0.9,Cdc,'o',70)
plot([0.55 1.45],[median(Ax2_Acti) median(Ax2_Acti)],'r--','linewidth',2)
% Spread_Julien(Dicty_Acti,2,0.9,Cdic,'o',70)
plot([1.55 2.45],[median(Ax3_Acti) median(Ax3_Acti)],'r--','linewidth',2)
% Spread_Julien(Ax230_Acti,3,0.9,Cfib,'o',70)
plot([2.55 3.45],[median(Ax230_Acti) median(Ax230_Acti)],'r--','linewidth',2)

%stats
DCDic = DoParamStats(Ax2_Acti,Ax3_Acti);
DCF = DoParamStats(Ax2_Acti,Ax230_Acti);
FDic = DoParamStats(Ax230_Acti,Ax3_Acti);

plot([1 2],[800 800],'k-','linewidth',1.5)
text(1.5, 820, DCDic,'HorizontalAlignment','center','fontsize',11)
plot([1 3],[850 850],'k-','linewidth',1.5)
text(2, 870, DCF,'HorizontalAlignment','center','fontsize',11)
plot([2 3],[920 920],'k-','linewidth',1.5)
text(2.5, 940, FDic,'HorizontalAlignment','center','fontsize',11)

% ax2 = gca;
% ax2.XTick = [1 2 3];
% ax2.XTickLabel = {'DCs' 'Dictys' '3T3'};
% ax2.XTickLabelRotation = 45;

ax = gca;
ax.LineWidth = 1.5;
box on
ylim([0 1600])
xlim([0.45 3.55])


%% Graph 3 : correlation
figure(3)
subplot(121)
hold on
ylabel('Fluctuations amplitude (nm)')
xlabel('Cortex thickness (nm)')
plot(Ax2_MedC,Ax2_Acti,'o','color',Cdic,'markerfacecolor',Cdic)
ax = gca;
h =lsline(ax);
ylim([-200 2000])
xlim([0 1000])

corrmat = corrcoef(Ax2_MedC,Ax2_Acti);
slope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));

DCfit = fitlm(Ax2_MedC,Ax2_Acti);
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
plot(Ax230_MedC,Ax230_Acti,'o','color',Cfib,'markerfacecolor',Cfib)
ax = gca;
h =lsline(ax);
ylim([-200 2000])
xlim([0 1000])

corrmat = corrcoef(Ax230_MedC,Ax230_Acti);
slope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));

Dictyfit = fitlm(Ax230_MedC,Ax230_Acti);
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


