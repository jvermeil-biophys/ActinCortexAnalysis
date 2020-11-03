function Stats2PlotsDictys(Type,savespecif,DIST,PEAKS,PROBA,TIME,alignlvl,datafolder,figurefolder,figurefolder2)

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

Cl  = [73 0 146]./255; % Latranculin


col = cell(nt,1);

leg = cell(nt,1);
label = cell(nt,1);

% choix des couleurs
name1 = cell(1,nt);
for kt = 1:nt
    
    if contains(Type{kt},'DMSO')
        tag1 = 'dm';
        name1{kt} = 'DMSO'
    elseif contains(Type{kt},'LatA')
        tag1 = 'l';
        name1{kt} = 'LatA';
    
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
        CortThick{kt} = MS.MedCortWs;
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

%% Cortex thickness

    figure
    hold on
    for kt = 1:nt
        
        Spread_Julien(CortThick{kt},kt,0.7,col{kt},'o',50)
        plot([kt-0.38 kt+0.38],[median(CortThick{kt}) median(CortThick{kt})],'k--','linewidth', 3)
        
    end
    
    ypos = prctile(CortThick{1},95); % for *** positioning
    
    for kt = 1:nt-1 % :nt-1 pour tout comparer
        for kkt = kt+1:nt
            pC = ranksum(CortThick{kt},CortThick{kkt});
            
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
    
    ylabel('Cortex thickness (nm)','interpreter','latex','fontsize',30)
    
    fig = gcf;
    saveas(fig,[ff2 filesep 'CortThick.fig'],'fig')
    saveas(fig,[ff2 filesep 'CortThick.png'],'png')
   
%% plot de distribution des valeur d'acti 10-90
    figure
    hold on
    for kt = 1:nt
        
        Spread_Julien(Acti1090{kt},kt,0.7,col{kt},'o')
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


%% distribution de tt les points de mesures
if DIST
    
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
%     close(fig)
    
%     figure
%     hold on
%     edges = -500:25:2500;
%     for kt = 1:nt
%         
%         N=histcounts(AllD3{kt},edges,'normalization','probability');
%         xDist = edges(1:end-1)-diff(edges)/2;
%         
%         plot(xDist,100*N,'color',col{kt},'linewidth',2)
%         
%     end
%     
%     ylabel('% of points')
%     xlabel('Measured Height (nm)')
%     legend(leg)
%     
%     fig = gcf;
%     
%     saveas(fig,[ff2 filesep savespecif '_DistCurves.fig'],'fig')
%     saveas(fig,[ff2 filesep savespecif '_DistCurves.png'],'png')
%     if nargin == 10
%         saveas(fig,[df  filesep savespecif '_DistCurves.fig'],'fig')
%         saveas(fig,[df  filesep savespecif '_DistCurves.png'],'png')
%     end
% %     close(fig)
    
end

%% Bar graph sur les stats des pics
if PEAKS
    %% Distribution de pics
    
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
      
end

%% Plot de proba 
if PROBA
   
    %% courbe proba mean sur les hauteurs
    figure
    hold on
    plot([0 0 0],[0.001 50 100],'--k','handlevisibility','off')
    plot([-500 1500],[alignlvl alignlvl]*100,'--k','handlevisibility','off')
    
    
    Pref = (99.95:-0.1:0.05)/100;
    
    for kt= 1:nt
        nn = length(ProbData{kt});
        for kn = 1:nn
            Probtmp = ProbData{kt}{kn};
            Ptmp = (length(Probtmp):-1:1)/length(Probtmp);
            Proba(kn,:) = interp1(Ptmp,Probtmp,Pref,'linear',max(Probtmp));
            %             Proba(kn,:) = interp1(Ptmp,Probtmp,Pref);
            
        end
        mProba = mean(Proba,1);
        err = std(Proba,1)/sqrt(nn);
        
        plot(mProba,Pref*100,'*','linewidth',3,'color',col{kt})
        
        plot(mProba-err,Pref*100,'--','color',col{kt},'handlevisibility','off')
        plot(mProba+err,Pref*100,'--','color',col{kt},'handlevisibility','off')
        
        clear Proba err
        
    end
    
    xlim([-200 800])
    xlabel('Height above cortex bottom thickness (nm)')
    ylabel('Fraction data points (%)')
    legend(leg)
    
    fig = gcf;
    saveas(fig,[ff2 filesep 'ProbDistLinV.fig'],'fig')
    saveas(fig,[ff2 filesep 'ProbDistLinV.png'],'png')
    if nargin == 10
        saveas(fig,[df filesep 'ProbDistLinV.fig'],'fig')
        saveas(fig,[df filesep 'ProbDistLinV.png'],'png')
    end

   
    
    
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
