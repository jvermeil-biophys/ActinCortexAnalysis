function Stats2Plots(Type,savespecif,DIST,PEAKS,PROBA,TIME,alignlvl,datafolder,figurefolder,figurefolder2)

%% INIT

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',23)


h = 'C:\Users\Valentin\Documents\MATLAB';


%% params and settings

% data location
DataF = [datafolder filesep 'P2S'];

% figures saving location
ff = figurefolder;

date = datestr(now,'yy-mm-dd');
ff1 = [ff filesep date];
ff2 = [ff filesep date filesep savespecif];
if nargin == 10
    df = [figurefolder2 filesep savespecif];
    mkdir(df)
end
mkdir(ff2)





% loop over datas
nt = length(Type);

% magnetisation a 5mT
B = 5;
Rb = (4432+4447)/4; % rayon bille moyen mesuré sur deux batch

M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1

V = 4/3*pi*Rb^3; % en nm^3

m = M.*10^-9*V; % moment dipolaire en A.nm^2

% Prefacteur travail de la force
mu0 = 12.6*10^6; % en mT.nm/A
C = 3*mu0/2/pi*m^2;

% Couleurs par drogues
Cdm = [0 109 219]./255; % DMSO
Cnd = Cdm/2; % No drug

Cl  = [73 0 146]./255; % Latranculin

Cb  = [219 209 0]./255; % blebbi
Cy  = [255 255 109 ]./255; % y27

Ccl = [146 0 0]./255; % calyculinA 1nM
Ccl3 = Ccl/2; % calyculinA 3nM
Ccl10 = Ccl3/2; % calyculinA 10nM

Cs  = [36 255 36]./255; % smifh2

Cml = [109 255 255]./255; % ml141

Cck = [182 109 255]./255;

Cin = [0 73 73]./255; % inside et autre
C0 = 'k';
C5 = 'b';
C10 = 'r';
C12 = 'm';


Crp = [200 0 20]./255; % RPE1


Co  = [255 182 119]./255; % outside et autre
Cin2 = Cin/2; % inside et autre
Co2  = Co/2; % outside et autre

col = cell(nt,1);

leg = cell(nt,1);
label = cell(nt,1);

% choix des couleurs
name1 = cell(1,nt);
for kt = 1:nt
    
    if contains(Type{kt},'ML141')
        tag1 = 'ml';
        name1{kt} = 'ML141';
%     elseif contains(Type{kt},'Dicty')
%         tag1 = 'cl3';
%         name1{kt} = 'Dictys';
    elseif contains(Type{kt},'RPE1')
        tag1 = 'rp';
        name1{kt} = 'RPE1';
    elseif contains(Type{kt},'DMSO')
        tag1 = 'dm';
        name1{kt} = 'DCs';
    elseif contains(Type{kt},'ActiUpVar')
        tag1 = 'in2';
        name1{kt} = 'NicoUp';
    elseif contains(Type{kt},'ActiUp')
        tag1 = 'in';
        name1{kt} = 'NicoUp';
    elseif contains(Type{kt},'ActiDownVar')
        tag1 = 'o2';
        name1{kt} = 'NicoDown';
    elseif contains(Type{kt},'ActiDown')
        tag1 = 'o';
        name1{kt} = 'NicoDown';
    elseif contains(Type{kt},'NoDrug')
        tag1 = 'nd';
        name1{kt} = 'NoDrug';
    elseif contains(Type{kt},'MyoII-WT')
        tag1 = 'dm';
        name1{kt} = 'MyoII-WT';
    elseif contains(Type{kt},'LatA500')
        tag1 = 'l';
        name1{kt} = 'LatA';
    elseif contains(Type{kt},'Blebbi')
        tag1 = 'b';
        name1{kt} = 'Blebbi';
    elseif contains(Type{kt},'MyoII-KO')
        tag1 = 'b';
        name1{kt} = 'MyoII-KO';
    elseif contains(Type{kt},'Y27')
        tag1 = 'y';
        name1{kt} = 'Y27';
    elseif contains(Type{kt},'CalA10')
        tag1 = 'cl10';
        name1{kt} = 'CalA10';
    elseif contains(Type{kt},'CalA1-1')
        tag1 = 'cl';
        name1{kt} = 'CalA1-1';
    elseif contains(Type{kt},'CalA1-2')
        tag1 = 'cl3';
        name1{kt} = 'CalA1-2';
    elseif contains(Type{kt},'CalA1-3')
        tag1 = 'cl10';
        name1{kt} = 'CalA1-3';
    elseif contains(Type{kt},'CalA1-4')
        tag1 = 'ck';
        name1{kt} = 'CalA1-4';
        
    elseif contains(Type{kt},'CalA1')
        tag1 = 'cl';
        name1{kt} = 'CalA1';
    elseif contains(Type{kt},'CalA3')
        tag1 = 'cl3';
        name1{kt} = 'CalA3';
    elseif contains(Type{kt},'SMIFH2')
        tag1 = 's';
        name1{kt} = 'Smifh2';
    elseif contains(Type{kt},'CK666')
        tag1 = 'ck';
        name1{kt} = 'CK666';
    elseif contains(Type{kt},'ArpinWT')
        tag1 = 'in';
        name1{kt} = 'ArpinWT';
    elseif contains(Type{kt},'ArpinKO')
        tag1 = 'o';
        name1{kt} = 'ArpinKO';
    elseif contains(Type{kt},'ArpC4WT')
        tag1 = 'in';
        name1{kt} = 'ArpC4WT';
    elseif contains(Type{kt},'ArpC4KO')
        tag1 = 'o';
        name1{kt} = 'ArpC4KO';
    elseif contains(Type{kt},'Inside_0pc')
        tag1 = '0';
        name1{kt} = '0\%';
    elseif contains(Type{kt},'Inside_5pc')
        tag1 = '5';
        name1{kt} = '5\%';
    elseif contains(Type{kt},'Inside_7-5')
        tag1 = 'in';
        name1{kt} = '7.5\%';
    elseif contains(Type{kt},'Inside_10')
        tag1 = '10';
        name1{kt} = '10\%';
    elseif contains(Type{kt},'Inside_12-5')
        tag1 = '12';
        name1{kt} = '12.5\%';
    elseif contains(Type{kt},'Inside_1mT')
        tag1 = '5';
        name1{kt} = '1mT';
    elseif contains(Type{kt},'Inside_3mT')
        tag1 = '10';
        name1{kt} = '3mT';
    elseif contains(Type{kt},'Inside_5mT')
        tag1 = '12';
        name1{kt} = '5mT';
    elseif contains(Type{kt},'Inside_15')
        tag1 = 'ck';
        name1{kt} = '15\%';
    elseif contains(Type{kt},'Outside')
        tag1 = 'o';
        name1{kt} = 'Outside';
    elseif contains(Type{kt},'DCLA_SiVim+Lat')
        tag1 = 'o';
        name1{kt} = 'SiVim+L';
    elseif contains(Type{kt},'DCLA_SiCtrl+Lat')
        tag1 = 'ck';
        name1{kt} = 'SiCt+L';
    elseif contains(Type{kt},'DCLA_SiVim')
        tag1 = 'o';
        name1{kt} = 'SiVim';
    elseif contains(Type{kt},'DCLA_SiCtrl')
        tag1 = 'in';
        name1{kt} = 'SiCt';
    end
    
    
    
    
    eval(['col{kt} = C' tag1 ';']);
    leg{kt} = [name1{kt}];
    
end



% for data storing

P = cell(nt,1); % prominence
PN = cell(nt,1); %prominence normalized to cortex local
H = cell(nt,1); % absolute height
D = cell(nt,1); % duration
PeakCortBot = cell(nt,1); % épaisseur cortex (globale) pour chaque pic
MedCw = cell(nt,1);
nbBI = cell(nt,1);

Ttot = cell(nt,1); % temps total de manip
ProbDist= cell(nt,1);
ProbData= cell(nt,1);
ProbDataN= cell(nt,1);
Acti1090 = cell(nt,1);



%% retrieve data


for kt = 1:nt
    
    
    filename = ['P2S_' Type{kt}];
    
    
    if exist([DataF filesep filename '.mat']) ~= 0
        load([DataF filesep filename])
        fprintf(['File loaded : ' filename '\n']);
        
        Data = MS.data;
        
        
        % get data
        P{kt} = MS.proms;
        PN{kt} = MS.promsN;
        H{kt} = MS.heights;
        PeakCortBot{kt} = MS.peaksCortW;
        D{kt} = MS.widths;
        Ttot{kt} = MS.timetot;
        ProbDist{kt} = MS.D3DataPts;
        ProbData{kt} = MS.D3Data;
        ProbDataN{kt} = MS.D3DataN;
        AllD3{kt} = MS.AllD3;
        AllD3N{kt} = MS.AllD3N;
        
        Acti1090{kt} = MS.CortActis;
        Assym{kt} = MS.CortAssyms;
        
        MedCw{kt} = MS.MedCortWs;
        nbBI{kt} = MS.nbBI;
        
    end
    
    filenameTime = ['TimeDatas_' Type{kt}];
    
    
    if exist([DataF filesep filenameTime '.mat']) ~= 0
        load([DataF filesep filenameTime])
        fprintf(['File loaded : ' filenameTime '\n']);
        
        
        
        
        % get data
        AllDistances{kt} = Distances;
        AllTimes{kt} = Times;
        
        
        
    end
    
    label{kt} = [leg{kt} '\newline' num2str(length(MS.data.cell))];
    leg{kt} = [leg{kt} ' - ' num2str(length(MS.data.cell))];
    
    
    
    
end

%% Plotting

%% nb bille int vs cortex med et fluc amp

figure

% for kt = 1:nt
%     subplot(2,nt,kt)
%     
%     plot(nbBI{kt},MedCw{kt},'ko','markerfacecolor',col{kt})
%     h = lsline;
%     if kt == ceil(nt/2)
%     xlabel('Nb of beads inside cell')
%     end
%     if kt ==1
%     ylabel('Cortex width (nm)')
%     end
%     xlim([0 6.5])
%     
%     slope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
%     
%     linfit = fitlm(nbBI{kt},MedCw{kt});
%     ConfInt = coefCI(linfit);
%     
%     SlopeCi = mean(abs(slope-ConfInt(2,:)));
% 
% txt = {['Slope = ' num2str(slope,'%.2f') ' $\pm$ ' num2str(SlopeCi,'%.2f')]};
% 
% Htxt = prctile(Acti1090{kt},100);
% 
%  text(0.2,Htxt*1.2,txt,'Interpreter','Latex','FontSize',12)
% 
%  ylim([0 Htxt*1.4])
%  
%  ax = gca;
%  
%  ax.XTick = [1 2 3 4 5 6];
% 
% end





for kt = 1:nt
% %     subplot(2,nt,nt+kt)
    
    plot(nbBI{kt},Acti1090{kt},'ko','markerfacecolor',col{kt})
    h = lsline;
    if kt == ceil(nt/2)
    xlabel('Nb of beads inside cell')
    end
    if kt == 1
    ylabel('Fluctuations Amplitude (nm)')
    end
    xlim([0 6.5])
    
    slope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
    
    linfit = fitlm(nbBI{kt},Acti1090{kt});
    
    Pval = linfit.Coefficients.pValue;
    
    ConfInt = coefCI(linfit);
    
    SlopeCi = mean(abs(slope-ConfInt(2,:)));

% txt = {['Slope = ' num2str(slope,'%.2f') ' $\pm$ ' num2str(SlopeCi,'%.2f')]};

txt = {['P-value : ' num2str(Pval(2))]}

Htxt = prctile(Acti1090{kt},100);

 text(0.2,Htxt*1.2,txt,'Interpreter','Latex','FontSize',12)

 ylim([0 Htxt*1.4])
 
  ax = gca;
 
 ax.XTick = [1 2 3 4 5 6];
 
end
%% distribution de tt les points de mesures
if DIST
    
    %% Distribution of normalized data points
    figure
    hold on
    edges = (-5:0.25:20) -0.125;
    for kt = 1:nt
        
        N=histcounts(AllD3N{kt},edges,'normalization','probability');
        xDist = edges(1:end-1)-diff(edges)/2;
        
        bar(xDist,100*N,'facecolor',col{kt},'edgecolor','none',...
            'facealpha',0.5,'linewidth',2)
        
    end
    
    ylabel('% of points')
    xlabel('Normalized Height (mean normalized)')
    legend(leg)
    
    fig = gcf;
    
    saveas(fig,[ff2 filesep savespecif '_DistNorm.fig'],'fig')
    saveas(fig,[ff2 filesep savespecif '_DistNorm.png'],'png')
    if nargin == 10
        saveas(fig,[df  filesep savespecif '_DistNorm.fig'],'fig')
        saveas(fig,[df  filesep savespecif '_DistNorm.png'],'png')
    end
    close(fig)
    
    %% Distribution of data points + distrib shapes
    figure
    hold on
    edges = -500:25:1500;
    for kt = 1:nt
        
        N=histcounts(AllD3{kt},edges,'normalization','probability');
        xDist = edges(1:end-1)-diff(edges)/2;
        
        bar(xDist,100*N,'facecolor',col{kt},'edgecolor','none',...
            'facealpha',0.7,'linewidth',2)
        
        
        
    end
    
    ylabel('% of points')
    xlabel('Measured thickness (nm)')
    legend(leg)
    
    fig = gcf;
    
    saveas(fig,[ff2 filesep savespecif '_DistNorm.fig'],'fig')
    saveas(fig,[ff2 filesep savespecif '_DistNorm.png'],'png')
    if nargin == 10
        saveas(fig,[df  filesep savespecif '_DistNorm.fig'],'fig')
        saveas(fig,[df  filesep savespecif '_DistNorm.png'],'png')
    end
    close(fig)
    
    figure
    hold on
    edges = -500:25:2500;
    for kt = 1:nt
        
        N=histcounts(AllD3{kt},edges,'normalization','probability');
        xDist = edges(1:end-1)-diff(edges)/2;
        
        plot(xDist,100*N,'color',col{kt},'linewidth',2)
        
    end
    
    ylabel('% of points')
    xlabel('Measured Height (nm)')
    legend(leg)
    
    fig = gcf;
    
    saveas(fig,[ff2 filesep savespecif '_DistCurves.fig'],'fig')
    saveas(fig,[ff2 filesep savespecif '_DistCurves.png'],'png')
    if nargin == 10
        saveas(fig,[df  filesep savespecif '_DistCurves.fig'],'fig')
        saveas(fig,[df  filesep savespecif '_DistCurves.png'],'png')
    end
    close(fig)
    
end

%% Bar graph sur les stats des pics
if PEAKS
    %% Distribution de pics (height, height/cort, prom, promnorm, width)
    
    % height - minimalcortex
    figure
    hold on
    edges = 0:50:1200;
    for kt = 1:nt
        n = histcounts(H{kt}-PeakCortBot{kt},edges);
        xval = edges(1:end-1)+mean(diff(edges))/2;
        yval{kt} = n/Ttot{kt}*600;
        errval{kt} = sqrt(n)/Ttot{kt}*600;
        bar(xval,yval{kt},'facecolor',col{kt},'edgecolor','none','facealpha',0.5)
        errorbar(xval,yval{kt},errval{kt},'.','color',col{kt},'handlevisibility','off')
        
    end
    xlabel('Peaks height - cortbot (nm)')
    ylabel('Frequency (peak/10min)')
    ylim([0 6])
    legend(leg)
    fig = gcf;
    saveas(fig,[ff2 filesep 'PeaksHeightAboveDist.fig'],'fig')
    saveas(fig,[ff2 filesep 'PeaksHeightAboveDist.png'],'png')
    %     close(fig)
    clear xval yval n
    
    % prom
    figure
    hold on
    edges = [0:33:1200];
    for kt = 1:nt
        n = histcounts(P{kt},edges);
        xval = edges(1:end-1)+mean(diff(edges))/2;
        yval{kt} = n/length(P{kt})*100;%/Ttot{kt}*600;
        errval{kt} = sqrt(n)/length(P{kt})*100;%/Ttot{kt}*600;
        bar(xval,yval{kt},'facecolor',col{kt},'edgecolor','none','facealpha',0.5)
        errorbar(xval,yval{kt},errval{kt},'.','color',col{kt},'handlevisibility','off')
    end
    
    %     for kt = 2:nt
    %         plot(xval,yval{kt}./yval{1},'o-','markerfacecolor',col{kt},'color',col{kt})
    %     end
    xlabel('Peaks prominence (nm)')
    ylabel('Proportion of all peaks (%)')
    %     ylim([-0.1 2])
    %     xlim([0 300])
    
    %     plot([0 300],[1 1],'k--','handlevisibility','off')
    % plot([0 300],[0.1 0.1],'k--','handlevisibility','off')
    
    legend(leg)
    fig = gcf;
    saveas(fig,[ff2 filesep 'PeaksPromDist.fig'],'fig')
    saveas(fig,[ff2 filesep 'PeaksPromDist.png'],'png')
    %     close(fig)
    
    
    %% Width vs height and prom
    
    %       figure
    %     hold on
    %
    %     for kt = 1:nt
    %         plot(D{kt},H{kt},'.','color',col{kt})
    %     end
    %         xlabel('Peaks widths (halfprom, sec)')
    %     ylabel('Peaks height (nm)')
    %
    %     legend(leg)
    %     fig = gcf;
    %     saveas(fig,[ff2 filesep 'PeaksHeightWidth.fig'],'fig')
    %     saveas(fig,[ff2 filesep 'PeaksHeightWidth.png'],'png')
    %     close(fig)
    
end

%% Plot de proba et mesure acti 1090
if PROBA
    %% courbe proba mean sur les hauteurs
%     figure
%     hold on
%     plot([0 0 0],[0.001 50 100],'--k','handlevisibility','off')
%     plot([-500 1500],[alignlvl alignlvl]*100,'--k','handlevisibility','off')
%     
%     
%     Pref = (99.95:-0.1:0.05)/100;
%     
%     for kt= 1:nt
%         nn = length(ProbData{kt});
%         for kn = 1:nn
%             Probtmp = ProbData{kt}{kn};
%             Ptmp = (length(Probtmp):-1:1)/length(Probtmp);
%             Proba(kn,:) = interp1(Ptmp,Probtmp,Pref,'linear',max(Probtmp));
%             %             Proba(kn,:) = interp1(Ptmp,Probtmp,Pref);
%             
%         end
%         mProba = mean(Proba,1);
%         err = std(Proba,1)/sqrt(nn);
%         
%         plot(mProba,Pref*100,'*','linewidth',3,'color',col{kt})
%         
%         plot(mProba-err,Pref*100,'--','color',col{kt},'handlevisibility','off')
%         plot(mProba+err,Pref*100,'--','color',col{kt},'handlevisibility','off')
%         
%         clear Proba err
%         
%     end
%     
%     xlim([-200 800])
%     xlabel('Height above cortex bottom thickness (nm)')
%     ylabel('Fraction data points (%)')
%     legend(leg)
%     
%     fig = gcf;
%     saveas(fig,[ff2 filesep 'ProbDistLinV.fig'],'fig')
%     saveas(fig,[ff2 filesep 'ProbDistLinV.png'],'png')
%     if nargin == 10
%         saveas(fig,[df filesep 'ProbDistLinV.fig'],'fig')
%         saveas(fig,[df filesep 'ProbDistLinV.png'],'png')
%     end
%     set(gca,'YScale','log')
%     ylim([1 100])
%     xlim([-200 800])
%     
%     ax = gca;
%     ax.YTickLabels = {'1' '10' '100'};
%     
%     fig = gcf;
%     saveas(fig,[ff2 filesep 'ProbDistLogV.fig'],'fig')
%     saveas(fig,[ff2 filesep 'ProbDistLogV.png'],'png')
%     if nargin == 10
%         saveas(fig,[df filesep 'ProbDistLogV.fig'],'fig')
%         saveas(fig,[df filesep 'ProbDistLogV.png'],'png')
%     end
%     close(fig)
%     
    %% courbe proba normalisé mean sur hauteurs
%     figure
%     hold on
%     plot([0 0 0],[0.001 50 100],'--k','handlevisibility','off')
%     plot([-100 50],[alignlvl alignlvl],'--k','handlevisibility','off')
%     
%     
%     Pref = (99.95:-0.1:0.05)/100;
%     
%     for kt= 1:nt
%         nn = length(ProbDataN{kt});
%         nop = 0;
%         for kn = 1:nn
%             if ~isnan(ProbDataN{kt}{kn})
%                 Probtmp = ProbDataN{kt}{kn};
%                 Ptmp = (length(Probtmp):-1:1)/length(Probtmp);
%                 try
%                     Proba(kn,:) = interp1(Ptmp,Probtmp,Pref);
%                 catch
%                     themall
%                 end
%             else
%                 
%                 Proba(kn,:) = NaN(1,length(Pref));
%                 nop = nop+1;
%             end
%         end
%         mProba = nanmean(Proba,1);
%         err = nanstd(Proba,1)/sqrt(nn-nop);
%         plot(mProba,Pref*100,'*','linewidth',3,'color',col{kt})
%         
%         plot(mProba-err,Pref*100,'--','color',col{kt},'handlevisibility','off')
%         plot(mProba+err,Pref*100,'--','color',col{kt},'handlevisibility','off')
%         
%         clear Proba err
%         
%     end
%     
%     legend(leg)
%     xlabel('Height above minimal cortex thickness (mean normalized)')
%     ylabel('Fraction of time (%)')
%     
%     fig = gcf;
%     saveas(fig,[ff2 filesep 'ProbNDistLinV.fig'],'fig')
%     saveas(fig,[ff2 filesep 'ProbNDistLinV.png'],'png')
%     if nargin == 10
%         saveas(fig,[df filesep 'ProbNDistLinV.fig'],'fig')
%         saveas(fig,[df filesep 'ProbNDistLinV.png'],'png')
%     end
%     set(gca,'YScale','log')
%     ylim([1 100])
%     
%     ax = gca;
%     ax.YTickLabels = {'1' '10' '100'};
%     
%     fig = gcf;
%     saveas(fig,[ff2 filesep 'ProbNDistLogV.fig'],'fig')
%     saveas(fig,[ff2 filesep 'ProbNDistLogV.png'],'png')
%     if nargin == 10
%         saveas(fig,[df filesep 'ProbNDistLogV.fig'],'fig')
%         saveas(fig,[df filesep 'ProbNDistLogV.png'],'png')
%     end
%     close(fig)
%     
%     
%     cd(h)
    
    %% plot de distribution des valeur d'acti 10-90
    figure
    hold on
    for kt = 1:nt
        
        Spread_Julien(Acti1090{kt},kt,0.7,col{kt},'o',50)
        plot([kt-0.38 kt+0.38],[median(Acti1090{kt}) median(Acti1090{kt})],'k--','linewidth', 3)
        
    end
    
    ypos = prctile(Acti1090{1},95); % for *** positioning
    
    for kt = 1:nt-1 % :nt-1 pour tout comparer
        for kkt = kt+1:nt
            pC = ranksum(Acti1090{kt},Acti1090{kkt});
            
            if 5/100 > pC && pC > 1/100
                txtC = ['*']  ;
            elseif 1/100 > pC && pC > 1/1000
                txtC = ['**'];
            elseif 1/1000 > pC
                txtC = ['***'];
            else
                txtC = ['NS'];
                
            end
            
            line([kt+0.1 kkt-0.1],[ypos*(1+((kkt-kt)/5)) ypos*(1+((kkt-kt)/5))],'color','k', 'linewidth', 1.3);
            text((kt+kkt)/2,ypos*(1.05+(kkt-kt)/5),txtC,'HorizontalAlignment','center','fontsize',11)
        end
    end
    
    ax = gca;
    ax.XTick = [1:nt];
    ax.XTickLabels = leg;
    ax.XTickLabelRotation = 0;
    
    ylabel('Fluctuation amplitude (nm)','interpreter','latex','fontsize',30)
    
    fig = gcf;
    saveas(fig,[ff2 filesep 'CortActiDist.fig'],'fig')
    saveas(fig,[ff2 filesep 'CortActiDist.png'],'png')
    
    
    %% plot de distribution des valeur d'assymetrie 50-90/10-50
    figure
    hold on
    for kt = 1:nt
        
        Spread_Julien(Assym{kt},kt,0.7,col{kt},'o',0.3)
        plot([kt-0.4 kt+0.4],[median(Assym{kt}) median(Assym{kt})],'r--','linewidth',1.5)
        
    end
    
    plot([0 nt+1],[0 0],'k--')
    
    ypos = prctile([Assym{:}],99.5); % for *** positioning
    
    for kt = 1 % :nt-1 pour tout comparer
        for kkt = kt+1:nt
            pC = ranksum(Assym{kt},Assym{kkt});
            
            if 5/100 > pC && pC > 1/100
                txtC = ['*']  ;
            elseif 1/100 > pC && pC > 1/1000
                txtC = ['**'];
            elseif 1/1000 > pC
                txtC = ['***'];
            else
                txtC = ['NS'];
                
            end
            
            line([kt+0.1 kkt-0.1],[ypos*(1+((kkt-kt)/5)) ypos*(1+((kkt-kt)/5))],'color','k', 'linewidth', 1.3);
            text((kt+kkt)/2,ypos*(1.05+(kkt-kt)/5),txtC,'HorizontalAlignment','center','fontsize',11)
        end
    end
    ax = gca;
    ax.XTick = [1:nt];
    ax.XTickLabels = leg;
    ax.XTickLabelRotation = 45;
    
    ylabel('Fluctuation assymetry')
    
    fig = gcf;
    saveas(fig,[ff2 filesep 'CortActiAssym.fig'],'fig')
    saveas(fig,[ff2 filesep 'CortActiAssym.png'],'png')
    
    
    ax = gca;
    ax.XTick = [1:5];
    ax.XTickLabels = {'DMSO','SMIFH2','CK666','Blebbi','LatA'};
    ax.XTickLabelRotation = 0;
    
end % end function


%% Etude temporelle
if TIME
    
    MinTime = 300;
    
    for kt = 1:nt
        
        Distances = AllDistances{kt};
        Times = AllTimes{kt};
        
        AllC{kt} = [];
        AllTau{kt} = [];
        
        Nok = 0;
        
%         figure
%         hold on
        
        for k = 1:length(Distances)
            T = Times{k};
            D3 = Distances{k};
            
            
            if (T(end) - T(1) > MinTime) &&  ((T(end)-T(1))/(length(T)-1) <= 1)
                Nok = Nok +1;
                
                [C,L] = AUTOCORR(T,D3,2);
           
                expfitfun = fittype('exp(-x/T)');

                expfit = fit(L',C',expfitfun);
                
                Tau = coeffvalues(expfit);
                
%                 figure
%                 plot(L,C)
%                 hold on
%                 plot(L,expfit(L))
               
                
                AllC{kt} = [AllC{kt}; C(1:(MinTime*2+1))];
                AllTau{kt} = [AllTau{kt} Tau];
                
            end
        end
        
        MedCorr{kt} = mean(AllC{kt});
        MedLag{kt} =  0:0.5:MinTime;
        ErrCorr{kt} = std(AllC{kt})/sqrt(Nok);
        
        
    end
    
    figure
    hold on
    
    for kt = 1:nt
plot(MedLag{kt},MedCorr{kt},'color',col{kt},'linewidth',2)
% plot(MedLag{kt},MedCorr{kt}+ErrCorr{kt},'--','color',col{kt},'linewidth',1,'handlevisibility','off')
% plot(MedLag{kt},MedCorr{kt}-ErrCorr{kt},'--','color',col{kt},'linewidth',1,'handlevisibility','off')
    end

xlabel('Lag (sec)')
ylabel('AutoCorrelation')

xlim([0 350])
ylim([-0.3 1])
legend(leg)

figure
hold on
for kt = 1:nt

Spread_Julien(AllTau{kt},kt,0.9,col{kt},'o',5)
plot([kt-0.45 kt+0.45],[median(AllTau{kt}) median(AllTau{kt})],'r--','linewidth',2)

end

ypos = prctile([AllTau{:}],99.5); % for *** positioning

    for kt = 2:nt % :nt-1 pour tout comparer
            pC = ranksum(AllTau{1},AllTau{kt});
            
            if 5/100 > pC && pC > 1/100
                txtC = ['*']  ;
            elseif 1/100 > pC && pC > 1/1000
                txtC = ['**'];
            elseif 1/1000 > pC
                txtC = ['***'];
            else
                txtC = ['NS'];
                
            end
            
            line([1+0.1 kt-0.1],[ypos*(1+((kt-1)/5)) ypos*(1+((kt-1)/5))],'color','k', 'linewidth', 1.3);
            text((1+kt)/2,ypos*(1.05+(kt-1)/5),txtC,'HorizontalAlignment','center','fontsize',11)
    end

ax = gca;
ax.XTick = 1:nt;
ax.XTickLabels = leg;

ylabel('Decorrelation time (s)')

    fig = gcf;
    saveas(fig,[ff2 filesep 'AutoCorrFit.fig'],'fig')
    saveas(fig,[ff2 filesep 'AutoCorrFit.png'],'png')

end
