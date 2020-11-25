function FigureMaker(Param,datf,figurefolder,figurefolder2)

close all

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',30)

h = 'C:\Users\Valentin\Documents\MATLAB';

% figures saving location
ff = figurefolder;

date = datestr(now,'yy-mm-dd');
ff = [ff filesep date];
if nargin == 4
    df = figurefolder2;
end
mkdir(ff)


% Couleurs par drogues
Cdm = [0 109 219]./255; % control
Cnd = Cdm/2;

Cl  = [73 0 146]./255; % Latranculin

Cb  = [219 209 0]./255; % blebbi
Cy  = [255 255 109 ]./255; % y27
Ccl = [146 0 0]./255; % calyculinA

Cs  = [36 255 36]./255; % smifh2

Cml = [109 255 255]./255; % ml141

Cck = [182 109 255]./255;

Cin = [0 73 73]./255; % inside
Co  = [255 182 119]./255; % outside

Crp = [200 0 20]./255; % RPE1



% if Param == 1
% Titlestr = '-CalA-Param';
% elseif Param == 0
% Titlestr = '-CalA-NonParam';
% end
% listctrl = [1 2 3 4];
% ncond = {4 3 2 1};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'    ,'DMSO'         ,1 ,Cdm,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_CalA1-1'  ,'CalA1-1' ,2 ,Ccl ,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_CalA1-2'  ,'CalA1-2' ,3 ,Ccl ,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DCLA_CalA1-3'  ,'CalA1-3' ,4 ,Ccl ,datf);
% [CBot{5}, MCw{5},CTop{5}] = Data2CortexStats('DCLA_CalA1-4'  ,'CalA1-4' ,5 ,Ccl ,datf);


% % 
% if Param == 1
% Titlestr = 'LAtInOut-Param';
% elseif Param == 0
% Titlestr = 'LatInOut-NonParam';
% end
% listctrl = [1];
% ncond = {3};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'   ,'DMSO'    ,1 ,Cdm,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DCLA_LatA500' ,'LatA'   ,4 ,Cl,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DC_Inside_7-5pc' ,'Inside'    ,2 ,Cin,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DC_TrueOutside' ,'Outside'   ,3 ,Co,datf);
% colorlist = {Cdm,Cl,Cin,Co};
% if Param == 1
% Titlestr = 'CKSM-Param';
% elseif Param == 0
% Titlestr = 'CKSM-NonParam';
% end
% listctrl = [1];
% ncond = {2};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'   ,'DMSO'    ,1 ,Cdm,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_CK666' ,'CK666'   ,2 ,Cck,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_SMIFH2' ,'Smifh2'    ,3 ,Cs,datf);

% if Param == 1
% Titlestr = 'BBCal-Param';
% elseif Param == 0
% Titlestr = 'BBCal-NonParam';
% end
% listctrl = [1 3];
% ncond = {3 2};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'   ,'DMSO'    ,1 ,Cdm,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_Blebbi' ,'Blebbi'   ,2 ,Cb,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_CalA1' ,'CalA - 1nM'    ,3 ,Ccl,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DCLA_CalA3' ,'CalA - 3nM'    ,4 ,Ccl,datf);
% [CBot{5}, MCw{5},CTop{5}] = Data2CortexStats('DCLA_CalA10' ,'CalA - 10nM'    ,5 ,Ccl,datf);


% if Param == 1
% Titlestr = 'DMSO-NDs-Param';
% elseif Param == 0
% Titlestr = 'DMSO-NDs-NonParam';
% end
% listctrl = [1 2];
% ncond = {2 1};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_NoDrug' ,'NoDrug'   ,1 ,Cnd,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_DMSO' ,'DMSO'   ,2 ,Cdm,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_NoDrug2' ,'NoDrug2'   ,3 ,Cnd,datf);



% if Param == 1
% Titlestr = 'Arpin-Param';
% elseif Param == 0
% Titlestr = 'Arpin-NonParam';
% end
% listctrl = [1];
% ncond = {1};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_ArpinWT' ,'ArpinWT'  ,1 ,Cin,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_ArpinKO' ,'ArpinKO'  ,2 ,Co,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_ArpinWT' ,'ArpinWT' ,3 ,Cin,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DCLA_ArpinKO' ,'ArpinKO' ,4 ,Co,datf);
% [CBot{5}, MCw{5},CTop{5}] = Data2CortexStats('DCLA_ArpC4WT' ,'ArpC4WT' ,5 ,Cin,datf);
% [CBot{6}, MCw{6},CTop{6}] = Data2CortexStats('DCLA_ArpC4KO' ,'ArpC4KO' ,6 ,Co,datf);


% if Param == 1
% Titlestr = 'LatA-Param';
% elseif Param == 0
% Titlestr = 'LatA-NonParam';
% end
% listctrl = [1 2];
% ncond = { 3 2};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'    ,'DMSO'      ,1 ,Cdm,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_LatA500' ,'LatA'      ,2 ,Cl,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_LatTryp' ,'LatA+Tryp' ,3 ,Cl/2,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DCLA_LatNase' ,'LatA+Nase' ,4 ,Cl/2,datf);
% 
% 
% if Param == 1
%     Titlestr = 'Fig4All-Param';
% elseif Param == 0
%     Titlestr = 'Fig4All-NonParam';
% end
% listctrl = [1];
% ncond = {5};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'    ,'DMSO'         ,1 ,Cdm,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_LatA500'  ,'LatA' ,2 ,Cl ,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_Blebbi'  ,'Blebbistatin' ,3 ,Cb ,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DCLA_CK666'   ,'CK666'        ,4 ,Cck,datf);
% [CBot{5}, MCw{5},CTop{5}] = Data2CortexStats('DCLA_SMIFH2'  ,'Smifh2'       ,5 ,Cs ,datf);
% [CBot{6}, MCw{6},CTop{6}] = Data2CortexStats('DCLA_ML141'   ,'ML141'       ,6 ,Cml ,datf);
% 

if Param == 1
    Titlestr = 'Dicty-Param';
elseif Param == 0
    Titlestr = 'Dicty-NonParam';
end
listctrl = [1];
ncond = {1};
[CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'    ,'Dcs'         ,1 ,Cdm,datf);
[CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_NoDrugs'  ,'ND' ,2 ,Ccl/2 ,datf);


% if Param == 1
%     Titlestr = 'RPE1-Param';
% elseif Param == 0
%     Titlestr = 'RPE1-NonParam';
% end
% listctrl = [1];
% ncond = {1};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'    ,'Dcs'         ,1 ,Cdm,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('RPE1'  ,'RPE1' ,2 ,Crp ,datf);


% if Param == 1
% Titlestr = '-ML141-Param';
% elseif Param == 0
% Titlestr = '-ML141-NonParam';
% end
% listctrl = [1];
% ncond = {1};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DCLA_DMSO'    ,'DMSO'    ,1 ,Cdm,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_ML141'  ,'ML141'    ,2 ,Cml ,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_CalA3'   ,'CalA 1.5'   ,3 ,Ccl/2,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DCLA_CalA10'   ,'CalA 5'   ,4 ,Ccl/4,datf);

%
% if Param == 1
% Titlestr = 'KDVim-Param';
% elseif Param == 0
% Titlestr = 'KDVim-NonParam';
% end
% listctrl = [1 2 4 5];
% ncond = {1 1 1 1};
% [CBot{1}, MCw{1},CTop{1}] =  Data2CortexStats('DCLA_DMSO' ,'NoSiRNA'   ,1 ,Cdm,datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DCLA_SiCtrl' ,'SiCtrl'   ,2 ,Cin,datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DCLA_SiVim'  ,'SiVim'    ,3 ,Co,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DCLA_LatA500'  ,'LatANoSi'    ,4 ,Cl,datf);
% [CBot{5}, MCw{5},CTop{5}] = Data2CortexStats('DCLA_SiCtrl+LatA' ,'SiCtrlL'   ,5 ,Cin,datf);
% [CBot{6}, MCw{6},CTop{6}] = Data2CortexStats('DCLA_SiVim+LatA' ,'SiVimL'   ,6 ,Co,datf)


% if Param == 1
%     Titlestr = 'Optiprep-Param';
% elseif Param == 0
%     Titlestr = 'Optiprep-NonParam';
% end
% listctrl = [1];
% ncond = {4};
% [CBot{1}, MCw{1},CTop{1}] = Data2CortexStats('DC_Inside_0pc'    ,'0\% Optiprep'       ,1 ,'k',datf);
% [CBot{2}, MCw{2},CTop{2}] = Data2CortexStats('DC_Inside_5pc'    ,'5\% Optiprep'       ,2 ,'b',datf);
% [CBot{3}, MCw{3},CTop{3}] = Data2CortexStats('DC_Inside_7-5pc'    ,'7.5\% Optiprep'   ,3 ,Cin,datf);
% [CBot{4}, MCw{4},CTop{4}] = Data2CortexStats('DC_Inside_10pc'    ,'10\% Optiprep'     ,4 ,'r',datf);
% [CBot{5}, MCw{5},CTop{5}] = Data2CortexStats('DC_Inside_12-5pc'    ,'12.5\% Optiprep' ,5 ,'m',datf);
% [CBot{6}, MCw{6},CTop{6}] = Data2CortexStats('DC_Inside_15pc'    ,'15\% Optiprep'     ,6 ,Cck,datf);
% 

% Param interne de plot
coord = CBot;
for i = 1:length(CBot)
    
    coord{i}(:) = i;
    
end

% Pos ref


ll = length(listctrl);

for kl = 1:ll
    
    yCmax(kl) = max([CBot{listctrl(kl):listctrl(kl)+ncond{kl}}])*1.05;
    
    yCstep(kl) = 0.05/1.1*yCmax(kl);
    
    yCmaxM(kl) = mean(MCw{listctrl(kl)})*1.05;
    yCstepM(kl) = 0.05/1.1*yCmaxM(kl);
    
    
    yCmaxT(kl) = max([CTop{listctrl(kl):listctrl(kl)+ncond{kl}}])*1.05;
    yCstepT(kl) = 0.05/1.1*yCmaxT(kl);
    
end

%% test statistique

for ki = 1:ll
    for k=1:ncond{ki}
        i = listctrl(ki);
        ii = i+k;
        
        
        if Param
            pCM = ttest2(MCw{i},MCw{ii});
            pCM = 1-pCM;
        else
            pCM = ranksum(MCw{i},MCw{ii});
        end
        
        
        if 5/100 > pCM && pCM > 1/100
            txtC = ['*']  ;
        elseif 1/100 > pCM && pCM > 1/1000
            txtC = ['**'];
        elseif 1/1000 > pCM
            txtC = ['***'];
        else
            txtC = ['NS'];
            
        end
        
        figure(2000)
        line([i+0.1 ii-0.1], [yCmaxM(ki)*1.3+(2*k-2)*yCstepM(ki) yCmaxM(ki)*1.3+(2*k-2)*yCstepM(ki)],'color','k', 'linewidth', 1.3)
        text((i+ii)/2,yCmaxM(ki)*1.3 + (2*k-1)*yCstepM(ki),txtC,'HorizontalAlignment','center','fontsize',11)
        
        xlim([0 length(CBot)+1])
        
        figure(2003)
        line([i+0.1 ii-0.1], [yCmaxM(ki)*1.3+(2*k-2)*yCstepM(ki) yCmaxM(ki)*1.3+(2*k-2)*yCstepM(ki)],'color','k', 'linewidth', 1.3)
        text((i+ii)/2,yCmaxM(ki)*1.3 + (2*k-1)*yCstepM(ki),txtC,'HorizontalAlignment','center','fontsize',11)
        
        xlim([0 length(CBot)+1])
        
        figure(2001)
        line([i+0.1 ii-0.1], [yCmaxM(ki)*1.05+(2*k-2)*yCstepM(ki) yCmaxM(ki)*1.05+(2*k-2)*yCstepM(ki)],'color','k', 'linewidth', 1.3)
        text((i+ii)/2,yCmaxM(ki)*1.05 + (2*k-1)*yCstepM(ki),txtC,'HorizontalAlignment','center','fontsize',11)
        
        xlim([0 length(CBot)+1])
        
    end
end

%% final polish and saving


figure(2001)
plot([0 length(MCw)+1],[0 0],'k--')
ylabel('Median bead separation (nm)')
ymin = min([MCw{:}])-20;

ylim([-20 max(yCmaxM)*1.05+max(yCstepM)*2*max([ncond{:}])])
% ylim([-25 35])

ax = gca;
xtic = ax.XTick;
xlab = ax.XTickLabel;

lims = ylim;
yyaxis right
ax.YColor = [0 0 0];
ylim(lims + 4441)
ylabel('Median Bead center separation (nm)')
yyaxis left

ax.XTick = [0 xtic];
ax.XTickLabel = ['Condition :\newline     n ='; xlab];

savenamepng = ['Median' Titlestr '_Bar.png'];
savenamefig = ['Median' Titlestr '_Bar.fig'];

titlebase = ax.Title.String;

title([titlebase '-' Titlestr])
fig = gcf;
saveas(fig,[ff filesep savenamepng],'png')
saveas(fig,[ff filesep savenamefig],'fig')
if nargin == 4
    saveas(fig,[df filesep savenamepng],'png')
    saveas(fig,[df filesep savenamefig],'fig')
end

figure(2000)
plot([0 length(MCw)+1],[0 0],'k--')
ylabel('Median bead surface separation (nm)')
ymin = min([MCw{:}])-20;

ax = gca;
xtic = ax.XTick;
xlab = ax.XTickLabel;

ylim([ymin max([MCw{:}])*1.05])
lims = ylim;
yyaxis right
ax.YColor = [0 0 0];
ylim(lims + 4441)
ylabel('Median Bead center separation (nm)')
yyaxis left

ax.XTick = [0 xtic];
ax.XTickLabel = ['Condition :\newline     n ='; xlab];

savenamepng = ['Median' Titlestr '.png'];
savenamefig = ['Median' Titlestr '.fig'];
titlebase = ax.Title.String;

title([titlebase '-' Titlestr])
fig = gcf;
saveas(fig,[ff filesep savenamepng],'png')
saveas(fig,[ff filesep savenamefig],'fig')
if nargin == 4
    saveas(fig,[df filesep savenamepng],'png')
    saveas(fig,[df filesep savenamefig],'fig')
end

figure(2003)
xlim([0 length(CTop)+1])



figure(4000)
hold on
plot([0 length(MCw)+1],[0 0],'k--')
ylabel(['Bead separation (nm)'])
xlim([0 length(MCw)+1])

ax = gca;
xtic = ax.XTick;
xlab = ax.XTickLabel;

ylim([ymin 600])


ax.XTick = [0 xtic];
ax.XTickLabel = ['Condition :'; xlab];

savenamepng = ['BotMedTop' Titlestr '.png'];
savenamefig = ['BotMedTop' Titlestr '.fig'];
titlebase = ax.Title.String;

title([titlebase '-' Titlestr])
fig = gcf;
saveas(fig,[ff filesep savenamepng],'png')
saveas(fig,[ff filesep savenamefig],'fig')
if nargin == 4
    saveas(fig,[df filesep savenamepng],'png')
    saveas(fig,[df filesep savenamefig],'fig')
end

figure(4001)
hold on
plot([0 length(MCw)+1],[0 0],'k--')
ylabel('Median and 10-90 percentile')
xlim([0 length(MCw)+1])

ax = gca;
xtic = ax.XTick;
xlab = ax.XTickLabel;

ax.XTick = [0 xtic];
ax.XTickLabel = ['Condition :\newline     n ='; xlab];

savenamepng = ['BotMedTopSpread' Titlestr '.png'];
savenamefig = ['BotMedTopSpread' Titlestr '.fig'];
titlebase = ax.Title.String;

title([titlebase '-' Titlestr])
fig = gcf;
saveas(fig,[ff filesep savenamepng],'png')
saveas(fig,[ff filesep savenamefig],'fig')
if nargin == 4
    saveas(fig,[df filesep savenamepng],'png')
    saveas(fig,[df filesep savenamefig],'fig')
end


figure(5000)
savenamepng = ['MedAmplCorr' Titlestr '.png'];
savenamefig = ['MedAmplCorr' Titlestr '.fig'];
fig = gcf;
saveas(fig,[ff filesep savenamepng],'png')
saveas(fig,[ff filesep savenamefig],'fig')
if nargin == 4
    saveas(fig,[df filesep savenamepng],'png')
    saveas(fig,[df filesep savenamefig],'fig')
end
cd(h)
end
