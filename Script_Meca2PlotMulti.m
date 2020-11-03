%% INIT
clear all
close all

warning('off','all')

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',30)
set(0,'DefaultAxesFontName','Serif')
set(0,'DefaultAxesFontWeight','normal')
format short e

% Condition colors (3T3)
C3t3doxi = [0 109 219]./255; % 

C3t3noDrugs  = [219 209 0]./255; % 

C3t3doxiPLL  = [0 149 219]./255; % 

% Condition colors (Dicty Ax2)
Cbsa100opti7 = [0 109 219]./255; % BSA 100X 7.5% optiprep

Cpll100opti7  = [219 209 0]./255; % PLL 100X 7.5% optiprep

Cpll63opti30  = [73 0 146]./255; % PLL 63X 30% optiprep

Cflo63opti30 = [0 210 110]./255; % Flottante 63X 30% optiprep



%% Get Data
DataFolder = 'D:\Matlab Analysis\Data_Joseph\MatFiles\D2M';
cd(DataFolder)
SaveFolder = 'C:\Users\JosephVermeil\Desktop\Figures';

%% 3T3 H_median_thickness

DataFolder = 'D:\Matlab Analysis\Data_Joseph\MatFiles\P2S';
cd(DataFolder)

% 3t3doxi
Hmed3T3doxi =[];
load(['./P2S_3T3aSFL_BSA_doxi.mat'])
Hmed3T3doxi =[Hmed3T3doxi MS.MedCortWs];
Hmedmed3T3doxi = median(Hmed3T3doxi);

% 3t3 no drugs

Hmed3T3noDrugs =[];
load(['./P2S_3T3aSFL_BSA_nodrugs.mat'])
Hmed3T3noDrugs =[Hmed3T3noDrugs MS.MedCortWs];
Hmedmed3T3noDrugs = median(Hmed3T3noDrugs);

% N
display(['Cells for Hmed, doxi - N = ' num2str(length(Hmed3T3doxi))])
display(['Cells for Hmed, noDrugs - N = ' num2str(length(Hmed3T3noDrugs))])

% Plot

figure 
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

data = {Hmed3T3doxi,Hmed3T3noDrugs};
colorslist = {C3t3doxi,C3t3noDrugs};
nameslist = {'iMC', 'Control'};

set(0,'DefaultLineMarkerSize',8);
h = plotSpread_V(data,'distributionMarkers',{'o','o'},'distributionColors', ...
    colorslist,'xNames',nameslist,'spreadWidth',1);

h{3}.XTickLabelRotation = 0;

% weighted average
plot([0.55 1.45],[Hmedmed3T3doxi Hmedmed3T3doxi],'r--','linewidth',2)
plot([1.55 2.45],[Hmedmed3T3noDrugs Hmedmed3T3noDrugs],'r--','linewidth',2)

text(0.6, Hmedmed3T3doxi + 20, num2str(Hmedmed3T3doxi,'%.0f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')
text(2.4, Hmedmed3T3noDrugs + 20, num2str(Hmedmed3T3noDrugs,'%.0f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')

%stats

plotheight0 = 400;
plotheight = plotheight0;

for ii = 1:length(data)-1
    for jj = ii+1:length(data)
        
        sigtxt = DoParamStats(data{ii},data{jj});
        
        plot([ii jj],[plotheight plotheight],'k-','linewidth',1.5)
        text((ii+jj)/2, plotheight+plotheight0/50, sigtxt,'HorizontalAlignment','center','fontsize',13)
        
        plotheight = plotheight+plotheight0/20;
    end
end

ylabel({'10mT Median thickness (nm)'}, 'FontSize',28)
hold on

yl = [0 plotheight+5];
ylim(yl);

fig = gcf;
saveas(fig,[SaveFolder filesep '3T3aSFL_Hmed10mT2.png'],'png')
saveas(fig,[SaveFolder filesep '3T3aSFL_Hmed10mT.fig'],'fig')

%% 3T3 H0_Chad

DataFolder = 'D:\Matlab Analysis\Data_Joseph\MatFiles\D2M';
cd(DataFolder)

% 3t3doxi
H03T3doxi =[];

DATES3t3doxi = {'04-08-20_M2','05-08-20_M1','07-08-20_M1'};

doxiCellCount = 0;
for k = 1:length(DATES3t3doxi)
    
    load(['./D2M_' DATES3t3doxi{k} '_3T3aSFL_BSA_doxi.mat'])
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                doxiCellCount = doxiCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
    H03T3doxi =[H03T3doxi cell2mat(cellfun(@transpose,ChadsH0,'UniformOutput',false))];
    
end

H03T3doxi = H03T3doxi;
H0Median3T3doxi = median(H03T3doxi);

% 3t3 no drugs

H03T3noDrugs =[];

DATES3t3noDrugs = {'04-08-20_M1','05-08-20_M2','07-08-20_M2'};

noDrugsCellCount = 0;
for k = 1:length(DATES3t3noDrugs)
    
    load(['./D2M_' DATES3t3noDrugs{k} '_3T3aSFL_BSA_nodrugs.mat'])
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                noDrugsCellCount = noDrugsCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
    H03T3noDrugs =[H03T3noDrugs cell2mat(cellfun(@transpose,ChadsH0,'UniformOutput',false))];
    
end

H03T3noDrugs = H03T3noDrugs;
H0Median3t3noDrugs = median(H03T3noDrugs);

% N
display(['Cells for H0, doxi - N = ' num2str(doxiCellCount)])
display(['Ramps for H0, doxi - N = ' num2str(length(H03T3doxi))])
display(['Cells for H0, noDrugs - N = ' num2str(noDrugsCellCount)])
display(['Ramps for H0, noDrugs - N = ' num2str(length(H03T3noDrugs))])

% Plot

figure 
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

data = {H03T3doxi,H03T3noDrugs};
colorslist = {C3t3doxi,C3t3noDrugs};
nameslist = {'iMC','Control'};

set(0,'DefaultLineMarkerSize',8);
h = plotSpread_V(data,'distributionMarkers',{'o','o'},'distributionColors', ...
    colorslist,'xNames',nameslist,'spreadWidth',1);

h{3}.XTickLabelRotation = 0;

% weighted average
plot([0.55 1.45],[H0Median3T3doxi H0Median3T3doxi],'r--','linewidth',2)
plot([1.55 2.45],[H0Median3t3noDrugs H0Median3t3noDrugs],'r--','linewidth',2)

text(0.6, H0Median3T3doxi + 20, num2str(H0Median3T3doxi,'%.0f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')
text(2.4, H0Median3t3noDrugs + 20, num2str(H0Median3t3noDrugs,'%.0f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')

%stats

plotheight0 = 500;
plotheight = plotheight0;

for ii = 1:length(data)-1
    for jj = ii+1:length(data)
        
        sigtxt = DoParamStats(data{ii},data{jj});
        
        plot([ii jj],[plotheight plotheight],'k-','linewidth',1.5)
        text((ii+jj)/2, plotheight+plotheight0/50, sigtxt,'HorizontalAlignment','center','fontsize',13)
        
        plotheight = plotheight+plotheight0/20;
    end
end

ylabel({'Initial thickness (nm) - Chadwick'}, 'FontSize',26)
hold on

yl = [0 plotheight+5];
ylim(yl);

fig = gcf;
saveas(fig,[SaveFolder filesep '3T3aSFL_H0_Chad2.png'],'png')
saveas(fig,[SaveFolder filesep '3T3aSFL_H0_Chad.fig'],'fig')

%% 3T3 Echad

% 3t3doxi
E3t3doxi =[];
ECI3t3doxi = [];

DATES3t3doxi = {'04-08-20_M2','05-08-20_M1','07-08-20_M1'};

doxiCellCount = 0;
for k = 1:length(DATES3t3doxi)
    
    load(['./D2M_' DATES3t3doxi{k} '_3T3aSFL_BSA_doxi.mat'])
    
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                doxiCellCount = doxiCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
    E3t3doxi =[E3t3doxi cell2mat(cellfun(@transpose,Echad,'UniformOutput',false))];
    ECI3t3doxi = [ECI3t3doxi cell2mat(cellfun(@transpose,EchadConf,'UniformOutput',false))];
    
end

EchadWeight3t3doxi = (E3t3doxi./ECI3t3doxi).^2;

EchadWeightedMean3t3doxi = sum(E3t3doxi.*EchadWeight3t3doxi)/sum(EchadWeight3t3doxi);

% 3t3 no drugs

E3t3noDrugs = [];
ECI3t3noDrugs = [];

DATES3t3noDrugs = {'04-08-20_M1','05-08-20_M2','07-08-20_M2'};

noDrugsCellCount = 0;
for k = 1:length(DATES3t3noDrugs)
    
    load(['./D2M_' DATES3t3noDrugs{k} '_3T3aSFL_BSA_nodrugs.mat'])
    
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                noDrugsCellCount = noDrugsCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
    E3t3noDrugs =[E3t3noDrugs cell2mat(cellfun(@transpose,Echad,'UniformOutput',false))];
    ECI3t3noDrugs = [ECI3t3noDrugs cell2mat(cellfun(@transpose,EchadConf,'UniformOutput',false))];
    
end

EchadWeight3t3noDrugs = (E3t3noDrugs./ECI3t3noDrugs).^2;
EchadWeightedMean3t3noDrugs = sum(E3t3noDrugs.*EchadWeight3t3noDrugs)/sum(EchadWeight3t3noDrugs);

% N
display(['Cells for Echad, doxi - N = ' num2str(doxiCellCount)])
display(['Ramps for Echad, doxi - N = ' num2str(length(E3t3doxi))])
display(['Cells for Echad, noDrugs - N = ' num2str(noDrugsCellCount)])
display(['Ramps for Echad, noDrugs - N = ' num2str(length(ECI3t3noDrugs))])

% Plot

figure
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

data = {E3t3doxi/1000,E3t3noDrugs/1000};
colorslist = {C3t3doxi,C3t3noDrugs};
nameslist = {'iMC', 'Control'};

set(0,'DefaultLineMarkerSize',8);
h = plotSpread_V(data,'distributionMarkers',{'o','o'},'distributionColors', ...
    colorslist,'xNames',nameslist,'spreadWidth',1);

h{3}.XTickLabelRotation = 0;

% weighted average
plot([0.55 1.45],[EchadWeightedMean3t3doxi EchadWeightedMean3t3doxi]./1000,'r--','linewidth',2)
plot([1.55 2.45],[EchadWeightedMean3t3noDrugs EchadWeightedMean3t3noDrugs]./1000,'r--','linewidth',2)

text(0.6, EchadWeightedMean3t3doxi/1000 + 3, num2str(EchadWeightedMean3t3doxi/1000,'%.1f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')
text(2.4, EchadWeightedMean3t3noDrugs/1000 + 3, num2str(EchadWeightedMean3t3noDrugs/1000,'%.2f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')


%stats

plotheight = 100;

for ii = 1:length(data)-1
    for jj = ii+1:length(data)
        
        sigtxt = DoParamStats(data{ii},data{jj});
        
        plot([ii jj],[plotheight plotheight],'k-','linewidth',1.5)
        text((ii+jj)/2, plotheight+2, sigtxt,'HorizontalAlignment','center','fontsize',13)
        
        plotheight = plotheight+5;
    end
end

ylabel({'Elastic modulus (kPa) - Chadwick'}, "fontsize", 26)
hold on

yl = [0 plotheight+5];
ylim(yl);

fig = gcf;
saveas(fig,[SaveFolder filesep '3T3aSFL_E_Chad2.png'],'png')
saveas(fig,[SaveFolder filesep '3T3aSFL_E_Chad.fig'],'fig')


% %% Plot Echad vs H0
% figure
% hold on
% plot(H03T3doxi,E3t3doxi,'ro')
% plot(H03T3noDrugs,E3t3noDrugs,'bo')

%% 3T3 Echad with spread cells

% 3t3doxi
E3t3doxi =[];
ECI3t3doxi = [];

DATES3t3doxi = {'04-08-20_M2','05-08-20_M1','07-08-20_M1'};

doxiCellCount = 0;
for k = 1:length(DATES3t3doxi)
    
    load(['./D2M_' DATES3t3doxi{k} '_3T3aSFL_BSA_doxi.mat'])
    
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                doxiCellCount = doxiCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
    E3t3doxi =[E3t3doxi cell2mat(cellfun(@transpose,Echad,'UniformOutput',false))];
    ECI3t3doxi = [ECI3t3doxi cell2mat(cellfun(@transpose,EchadConf,'UniformOutput',false))];
    
end

EchadWeight3t3doxi = (E3t3doxi./ECI3t3doxi).^2;

EchadWeightedMean3t3doxi = sum(E3t3doxi.*EchadWeight3t3doxi)/sum(EchadWeight3t3doxi);

% 3t3 no drugs

E3t3noDrugs = [];
ECI3t3noDrugs = [];

DATES3t3noDrugs = {'04-08-20_M1','05-08-20_M2','07-08-20_M2'};

noDrugsCellCount = 0;
for k = 1:length(DATES3t3noDrugs)
    
    load(['./D2M_' DATES3t3noDrugs{k} '_3T3aSFL_BSA_nodrugs.mat'])
    
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                noDrugsCellCount = noDrugsCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
    E3t3noDrugs =[E3t3noDrugs cell2mat(cellfun(@transpose,Echad,'UniformOutput',false))];
    ECI3t3noDrugs = [ECI3t3noDrugs cell2mat(cellfun(@transpose,EchadConf,'UniformOutput',false))];
    
end

EchadWeight3t3noDrugs = (E3t3noDrugs./ECI3t3noDrugs).^2;
EchadWeightedMean3t3noDrugs = sum(E3t3noDrugs.*EchadWeight3t3noDrugs)/sum(EchadWeight3t3noDrugs);

% 3t3 doxi PLL

E3t3doxiPLL = [];
ECI3t3doxiPLL = [];

DATES3t3doxiPLL = {'24-07-20_M2'};

doxiPLLCellCount = 0;
for k = 1:length(DATES3t3doxiPLL)
    
    load(['./D2M_' DATES3t3doxiPLL{k} '_3T3aSFL_Pll_doxi.mat'])
    
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                doxiPLLCellCount = doxiPLLCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
    E3t3doxiPLL =[E3t3doxiPLL cell2mat(cellfun(@transpose,Echad,'UniformOutput',false))];
    ECI3t3doxiPLL = [ECI3t3doxiPLL cell2mat(cellfun(@transpose,EchadConf,'UniformOutput',false))];
    
end

EchadWeight3t3doxiPLL = (E3t3doxiPLL./ECI3t3doxiPLL).^2;
EchadWeightedMean3t3doxiPLL = sum(E3t3doxiPLL.*EchadWeight3t3doxiPLL)/sum(EchadWeight3t3doxiPLL);

% N
display(['Cells for Echad, doxi - N = ' num2str(doxiCellCount)])
display(['Ramps for Echad, doxi - N = ' num2str(length(E3t3doxi))])
display(['Cells for Echad, noDrugs - N = ' num2str(noDrugsCellCount)])
display(['Ramps for Echad, noDrugs - N = ' num2str(length(E3t3noDrugs))])
display(['Cells for Echad, doxiPLL - N = ' num2str(doxiPLLCellCount)])
display(['Ramps for Echad, doxiPLL - N = ' num2str(length(E3t3doxiPLL))])

% Plot

figure
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on

data = {E3t3doxi/1000,E3t3noDrugs/1000,E3t3doxiPLL/1000};
colorslist = {C3t3doxi,C3t3noDrugs,C3t3doxiPLL};
nameslist = {'iMC', 'Control', 'iMC on PLL'};

set(0,'DefaultLineMarkerSize',8);
h = plotSpread_V(data,'distributionMarkers',{'o','o','o'},'distributionColors', ...
    colorslist,'xNames',nameslist,'spreadWidth',0.9);

h{3}.XTickLabelRotation = 0;

% weighted average
plot([0.55 1.45],[EchadWeightedMean3t3doxi EchadWeightedMean3t3doxi]./1000,'r--','linewidth',2)
plot([1.55 2.45],[EchadWeightedMean3t3noDrugs EchadWeightedMean3t3noDrugs]./1000,'r--','linewidth',2)
plot([2.55 3.45],[EchadWeightedMean3t3doxiPLL EchadWeightedMean3t3doxiPLL]./1000,'r--','linewidth',2)

text(0.6, EchadWeightedMean3t3doxi/1000 + 3, num2str(EchadWeightedMean3t3doxi/1000,'%.1f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')
text(1.6, EchadWeightedMean3t3noDrugs/1000 + 3, num2str(EchadWeightedMean3t3noDrugs/1000,'%.2f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')
text(2.6, EchadWeightedMean3t3doxiPLL/1000 + 3, num2str(EchadWeightedMean3t3doxiPLL/1000,'%.2f'),'HorizontalAlignment','center','fontsize',20,'color','k','FontWeight','bold')


%stats

plotheight = 100;

for ii = 1:length(data)-1
    for jj = ii+1:length(data)
        
        sigtxt = DoParamStats(data{ii},data{jj});
        
        plot([ii jj],[plotheight plotheight],'k-','linewidth',1.5)
        text((ii+jj)/2, plotheight+2, sigtxt,'HorizontalAlignment','center','fontsize',13)
        
        plotheight = plotheight+5;
    end
end

ylabel({'Elastic modulus (kPa) - Chadwick'}, "fontsize", 26)
hold on

yl = [0 plotheight+5];
ylim(yl);

fig = gcf;
saveas(fig,[SaveFolder filesep '3T3aSFL_E_Chad_withPLL.png'],'png')
saveas(fig,[SaveFolder filesep '3T3aSFL_E_Chad_withPLL.fig'],'fig')


% %% Plot Echad vs H0
% figure
% hold on
% plot(H03T3doxi,E3t3doxi,'ro')
% plot(H03T3noDrugs,E3t3noDrugs,'bo')

%% Dictys
    
% BSA 100X 7.5% optiprep
Ebsa100opti7 =[];
ECIbsa100opti7 = [];

DATESbsa100opti7 = {'27-05-20_M1','25-05-20_M1','22-05-20_M1','20-05-20_M1','20-05-20_M2'};

bsaCellCount = 0;

for k = 1:length(DATESbsa100opti7)
    
    load(['./D2M_' DATESbsa100opti7{k} '_DictyAx2_DMSO.mat'])
    
    Ebsa100opti7 =[Ebsa100opti7 cell2mat(cellfun(@transpose,Echad,'UniformOutput',false))];
    ECIbsa100opti7 = [ECIbsa100opti7 cell2mat(cellfun(@transpose,EchadConf,'UniformOutput',false))];
    
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                bsaCellCount = bsaCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
end

EchadWeightbsa100opti7 = (Ebsa100opti7./ECIbsa100opti7).^2;

EchadWeightedMeanbsa100opti7 = sum(Ebsa100opti7.*EchadWeightbsa100opti7)/sum(EchadWeightbsa100opti7);


% PLL 100X 7.5% optiprep

Epll100opti7J =[];
ECIpll100opti7J = [];

DATESpll100opti7J = {'24-07-20_M1'};

pllCellCount = 0;
for k = 1:length(DATESpll100opti7J)
    
    load(['./D2M_' DATESpll100opti7J{k} '_DictyAx2_PLLJ-7-5opti.mat'])
    
    Epll100opti7J =[Epll100opti7J cell2mat(cellfun(@transpose,Echad,'UniformOutput',false))];
    ECIpll100opti7J = [ECIpll100opti7J cell2mat(cellfun(@transpose,EchadConf,'UniformOutput',false))];
    
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                pllCellCount = pllCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
end

EchadWeightpll100opti7J = (Epll100opti7J./ECIpll100opti7J).^2;

EchadWeightedMeanpll100opti7J = sum(Epll100opti7J.*EchadWeightpll100opti7J)/sum(EchadWeightpll100opti7J);




% Flottante 63X 7.5% optiprep

Eflo63opti7 =[];
ECIflo63opti7 = [];

DATESflo63opti7 = {'06-08-20_M1'};

floatCellCount = 0;
for k = 1:length(DATESflo63opti7)
    
    load(['./D2M_' DATESflo63opti7{k} '_DictyComp_Float-7-80.mat'])
    
    Eflo63opti7 =[Eflo63opti7 cell2mat(cellfun(@transpose,Echad,'UniformOutput',false))];
    ECIflo63opti7 = [ECIflo63opti7 cell2mat(cellfun(@transpose,EchadConf,'UniformOutput',false))];
    
    currentName = '';
    for i = 1:length(CompressionFeatures)
        if ~isempty(CompressionFeatures(i).CellName)
            if CompressionFeatures(i).SAVING_CRITERION == 1 && ~strcmp(currentName, CompressionFeatures(i).CellName)
                floatCellCount = floatCellCount + 1;
                currentName = CompressionFeatures(i).CellName;
            end
        end
    end
    
end


EchadWeightflo63opti30 = (Eflo63opti7./ECIflo63opti7).^2;

EchadWeightedMeanflo63opti30 = sum(Eflo63opti7.*EchadWeightflo63opti30)/sum(EchadWeightflo63opti30);

% N
display(['Cells for Echad, bsa - N = ' num2str(bsaCellCount)])
display(['Ramps for Echad, bsa - N = ' num2str(length(Ebsa100opti7))])
display(['Cells for Echad, pll - N = ' num2str(pllCellCount)])
display(['Ramps for Echad, pll - N = ' num2str(length(Epll100opti7J))])
display(['Cells for Echad, float - N = ' num2str(floatCellCount)])
display(['Ramps for Echad, float - N = ' num2str(length(Eflo63opti7))])

% Plot

figure 
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.5;
box on


data = {Eflo63opti7/1000,Ebsa100opti7/1000,Epll100opti7J/1000};
colorslist = {Cflo63opti30,Cbsa100opti7,Cpll100opti7};
nameslist = {'Floating','BSA','PLL'};




set(0,'DefaultLineMarkerSize',6);
h = plotSpread_V(data,'distributionMarkers',{'o','o','o'},'distributionColors', ...
    colorslist,'xNames',nameslist,'spreadWidth',0.8);

h{3}.XTickLabelRotation = 0;

% weighted average
plot([0.55 1.45],[EchadWeightedMeanflo63opti30 EchadWeightedMeanflo63opti30]./1000,'r--','linewidth',2)
plot([1.55 2.45],[EchadWeightedMeanbsa100opti7 EchadWeightedMeanbsa100opti7]./1000,'r--','linewidth',2)
plot([2.55 3.45],[EchadWeightedMeanpll100opti7J EchadWeightedMeanpll100opti7J]./1000,'r--','linewidth',2)


text(1.7, EchadWeightedMeanbsa100opti7/1000 + 5, num2str(EchadWeightedMeanbsa100opti7/1000,'%.2f'),'HorizontalAlignment','center','fontsize',22,'color','k','FontWeight','bold')
text(2.7, EchadWeightedMeanpll100opti7J/1000 + 5, num2str(EchadWeightedMeanpll100opti7J/1000,'%.1f'),'HorizontalAlignment','center','fontsize',22,'color','k','FontWeight','bold')
text(0.7, EchadWeightedMeanflo63opti30/1000 + 5, num2str(EchadWeightedMeanflo63opti30/1000,'%.2f'),'HorizontalAlignment','center','fontsize',22,'color','k','FontWeight','bold')

%stats

plotheight = 100;

for ii = 1:length(data)-1
    for jj = ii+1:length(data)
        
        sigtxt = DoParamStats(data{ii},data{jj});
        
        plot([ii jj],[plotheight plotheight],'k-','linewidth',1.5)
        text((ii+jj)/2, plotheight+4, sigtxt,'HorizontalAlignment','center','fontsize',18)
        
        plotheight = plotheight+10;
    end
end


ylabel({'Elastic modulus (kPa) - Chadwick'}, 'fontsize', 26)
hold on

yl = [0 plotheight+5];
ylim(yl);

fig = gcf;
saveas(fig,[SaveFolder filesep 'DICTYS_E_Chad2.png'],'png')
saveas(fig,[SaveFolder filesep 'DICTYS_E_Chad.fig'],'fig')


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




