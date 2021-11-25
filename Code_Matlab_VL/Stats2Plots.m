function Stats2Plots(Cond,col,savespecif,CORT,DIST,PEAKS,ACTI,TIME,alignlvl,resfolder,figurefolder)

%
% Stats2Plots is used for analyzing constant field curves. It load the data
% from Peaks2Stats and for different experimental conditions. It makes a
% variety of plots that can be selected through the CORT, DIST, PEAKS, PROBA and
% TIME variables. 
%
% Stats2Plots(Cond,col,savespecif,DIST,PEAKS,ACTI,TIME,alignlvl,resfolder,figurefolder)
%
% Cond : list of experimental conditions (in bracket) to plot (ex :
% {'DictyAx2_DMSO','DictyAx2_LatA'})
% col : list of colors corresponding to each condition (in brackets, ex :
% {Cdm, Cb} or {'r', 'g'} or {[0.2 0.5 0.3], [ 0.4 0.58 0.14]} )
% savename for the folder containing figures
% Cort : 1 = plot cortex thickness, fluctuations and correlation between the two, 0 = no plot
% DIST : 1 = plot all points distributions, 0 = no plot
% PEAKS : 1 = plot peaks analysis, 0 = no plot
% ACTI : 1 = plot activity curves, 0 = no plot
% TIME : 1 = plot time analysis, 0 = no plot
% resfolder : path to the matlab data, should be MatFileFolder in main
% figurefolder : path to figure saving folder, should be FigureFolder
% in main
%

% figure options
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)
set(0, 'defaultAxesTickLabelInterpreter','none'); 
set(0, 'defaultLegendInterpreter','none');
warning('off','all')
close all


%% params and settings

% data location
DataF = [resfolder filesep 'P2S'];

% figures saving location
ff = figurefolder;

date = datestr(now,'yy-mm-dd');
ff2 = [ff filesep date filesep savespecif];
mkdir(ff2)


ncond = length(Cond); % number of conditions

leg = cell(ncond,1); % legend initialisation
label = cell(ncond,1); % labels initialisation

% preparing legends
for kC = 1:ncond
    
    pos_ = strfind(Cond{kC},'_');
    leg{kC} = Cond{kC};
    label{kC} = [Cond{kC}(1:pos_-1) '\newline' Cond{kC}(pos_+1:end)];
    
end


%% retrieve data
for kC = 1:ncond

    % loading stats
    filename = ['P2S_' Cond{kC}];    
    
    if exist([DataF filesep filename '.mat']) ~= 0 
        
        load([DataF filesep filename])
        fprintf(['File loaded : ' filename '\n']);       
      
        
        % get data
        CortThick{kC} = MS.MedCortWs; % cortex thicknesses
        P{kC} = MS.proms; % peaks prominence
        H{kC} = MS.heights; % peaks height
        PeakCortBot{kC} = MS.peaksCortW; % cell cortex thickness for each peak
        D{kC} = MS.widths; % peaks duration
        Ttot{kC} = MS.timetot; % total time of experiment
        ProbDist{kC} = MS.D3DataPts; % all datapoints centered on alignlvl
        AllD3{kC} = MS.AllD3; % all datapoints pulled together
        
        Acti1090{kC} = MS.CortActis; % fluctuations amplitude
        Assym{kC} = MS.CortAssyms; % cortex assymetry
        
        CortBot{kC} = MS.CortBots; % cortex bottom level (10%)
        CortTop{kC} = MS.CortTops; % cortex top level (90%)
        
        MedCw{kC} = MS.MedCortWs; % cortex median thickness
        
        ActivityCurves{kC} = MS.Acticurves; % activity curves
        
        Lags{kC} = MS.lags; % lags from autocorrelation
        Corrs{kC} = MS.corrs; % correlation values from autocorrelation
        
        MinTime = MS.mintimecorr; % minimum length of curve for autocorr
        Fs = MS.freqsamp; % autocorr sampling frequency
        
    end
    
    % loading time data
    filenameTime = ['TimeDatas_' Cond{kC}];
    
    if exist([DataF filesep filenameTime '.mat']) ~= 0
        
        load([DataF filesep filenameTime])
        fprintf(['File loaded : ' filenameTime '\n']);      
        
        % get data
        AllDistances{kC} = Distances;
        AllTimes{kC} = Times;
               
    end
    
    % finalizing labels and legens with numbers of cells
    label{kC} = [label{kC} '\newline' num2str(length(CortThick{kC}))];
    leg{kC} = [leg{kC} ' - ' num2str(length(CortThick{kC}))];     
    
end


%% Plotting

if CORT
%% Cortex thickness

    figure
    hold on
    for kC = 1:ncond
        
        Spread_Julien(CortThick{kC},kC,0.7,col{kC},'o',50)
        plot([kC-0.38 kC+0.38],[median(CortThick{kC}) median(CortThick{kC})],'k--','linewidth', 3)
        
    end
    
    ypos = prctile(CortThick{1},95); % for *** positioning
    
    % statistic comparison
    for kC = 1:ncond-1 % kC = 1:1 for comparing everything to the first condition only
        for kkt = kC+1:ncond
            pC = ranksum(CortThick{kC},CortThick{kkt});
            
                    
            txtsup = ['p = ' num2str(pC,'%.3f')];
            
            if 5/100 > pC && pC > 1/100
                txtC = ['*  ' txtsup]  ;
            elseif 1/100 > pC && pC > 1/1000
                txtC = ['**  ' txtsup];
            elseif 1/1000 > pC
                txtC = ['***  ' txtsup];
            else
                txtC = ['NS  ' txtsup];
                
            end
            
            line([kC+0.1 kkt-0.1],[ypos ypos],'color','k', 'linewidth', 1.3);
            text((kC+kkt)/2,ypos + 20,txtC,'HorizontalAlignment','center','fontsize',11)
            
            ypos = ypos + 30;
            
        end
    end
    
    % put labels
    ax = gca;
    ax.XTick = 1:ncond;
    ax.XTickLabels = leg;
    ax.XTickLabelRotation = 0;
    ax.TickLabelInterpreter = 'none';
    
    ylim([0 ypos+30])
    ylabel('Cortex thickness (nm)','interpreter','latex','fontsize',30)
    
    % save figure
    fig = gcf;
    saveas(fig,[ff2 filesep 'CortThick.fig'],'fig')
    saveas(fig,[ff2 filesep 'CortThick.png'],'png')
   
%% Fluctuations amplitude
    figure
    hold on
    for kC = 1:ncond
        
        Spread_Julien(Acti1090{kC},kC,0.7,col{kC},'o')
        plot([kC-0.38 kC+0.38],[median(Acti1090{kC}) median(Acti1090{kC})],'k--','linewidth', 3)
        
    end
    
    ypos = prctile(Acti1090{1},95); % for *** positioning
    
    % statistic comparison
    for kC = 1:ncond-1 % kC = 1:1 for comparing everything to the first condition only
        for kkt = kC+1:ncond
            pC = ranksum(Acti1090{kC},Acti1090{kkt});
            
            
            
            txtsup = ['p = ' num2str(pC,'%.3f')];
            
            if 5/100 > pC && pC > 1/100
                txtC = ['*  ' txtsup]  ;
            elseif 1/100 > pC && pC > 1/1000
                txtC = ['**  ' txtsup];
            elseif 1/1000 > pC
                txtC = ['***  ' txtsup];
            else
                txtC = ['NS  ' txtsup];
                
            end
            
          line([kC+0.1 kkt-0.1],[ypos ypos],'color','k', 'linewidth', 1.3);
            text((kC+kkt)/2,ypos + 20,txtC,'HorizontalAlignment','center','fontsize',11)
            
            ypos = ypos + 30;
        end
    end
    
    % put labels
    ax = gca;
    ax.XTick = 1:ncond;
    ax.XTickLabels = leg;
    ax.XTickLabelRotation = 0;
    ax.TickLabelInterpreter = 'none';
    
    ylim([0 ypos+30])
    ylabel('Fluctuation amplitude (nm)','interpreter','latex','fontsize',30)
    
    % save figure
    fig = gcf;
    saveas(fig,[ff2 filesep 'CortActiDist.fig'],'fig')
    saveas(fig,[ff2 filesep 'CortActiDist.png'],'png')

%% Median + 10-90 level plot
   
for kC = 1:ncond
    
    MedCw = CortThick{kC};
    P10 = CortBot{kC};
    P90 = CortTop{kC};
    
    figure
    hold on
    PlotMedP1090(MedCw,P10,P90,kC,col{kC})
    ylabel('Bead separation (nm)')
    xlabel(leg{kC})
    
end

%% Correlation Cortex/Activity
 for kC = 1:ncond

figure
hold on
ylabel('Fluctuation amplitude (nm)')
xlabel('Median cortex thickness (nm)')
plot(CortThick{kC},Acti1090{kC},'o','color','k','markerfacecolor',col{kC},'markersize',8)

Correlfit = fitlm(CortThick{kC},Acti1090{kC});

Pval = Correlfit.Coefficients.pValue

Xnew = linspace(0,1000,1000);
[y,~] = predict(Correlfit,Xnew','Alpha',0.05,'Simultaneous',true);


plot(Xnew,y,'linewidth',2,'color',col{kC})

ylim([-200 800])
xlim([0 500])
box on

legend(leg{kC},['Slope = ' num2str(Correlfit.Coefficients{2,1}) '. Pvalue = ' num2str(Pval(2))],'interpreter','none')
legend('location','northwest')
 end
 
end

if DIST
%% Distribution of all datapoints

    
    figure
    hold on
    title('Distributions of all data points per condition')
    edges = -500:25:1500;
    
    for kC = 1:ncond
        
        N=histcounts(AllD3{kC},edges,'normalization','probability');
        xDist = edges(1:end-1)-diff(edges)/2;
        
        bar(xDist,100*N,'facecolor',col{kC},'edgecolor','none',...
            'facealpha',0.7,'linewidth',2)
               
    end
    
    ylabel('% of points')
    xlabel('Measured thickness (nm)')
    legend(leg,'interpreter','none')
    
    fig = gcf;
    
    saveas(fig,[ff2 filesep savespecif '_AllPointsDisttribution.fig'],'fig')
    saveas(fig,[ff2 filesep savespecif '_AllPointsDistribution.png'],'png')   
    
end

if PEAKS
%% Peaks statistics
    
    % prominence distribution
    figure
    hold on
    edges = [0:33:1200];
    for kC = 1:ncond
        n = histcounts(P{kC},edges);
        xval = edges(1:end-1)+mean(diff(edges))/2;
        yval{kC} = n/length(P{kC})*100;%/Ttot{kt}*600;
        errval{kC} = sqrt(n)/length(P{kC})*100;%/Ttot{kt}*600;
        bar(xval,yval{kC},'facecolor',col{kC},'edgecolor','none','facealpha',0.5)
        errorbar(xval,yval{kC},errval{kC},'.','color',col{kC},'handlevisibility','off')
    end
    
    xlabel('Peaks prominence (nm)')
    ylabel('Proportion of all peaks (%)')

    legend(leg,'interpreter','none')
    fig = gcf;
    saveas(fig,[ff2 filesep 'PeaksPromDist.fig'],'fig')
    saveas(fig,[ff2 filesep 'PeaksPromDist.png'],'png')
    
    % heigth distribution
    figure
    hold on
    edges = [0:33:1200];
    for kC = 1:ncond
        n = histcounts(H{kC},edges);
        xval = edges(1:end-1)+mean(diff(edges))/2;
        yval{kC} = n/length(H{kC})*100;%/Ttot{kt}*600;
        errval{kC} = sqrt(n)/length(H{kC})*100;%/Ttot{kt}*600;
        bar(xval,yval{kC},'facecolor',col{kC},'edgecolor','none','facealpha',0.5)
        errorbar(xval,yval{kC},errval{kC},'.','color',col{kC},'handlevisibility','off')
    end
    
    xlabel('Peaks prominence (nm)')
    ylabel('Proportion of all peaks (%)')

    legend(leg,'interpreter','none')
    fig = gcf;
    saveas(fig,[ff2 filesep 'PeaksHeigthDist.fig'],'fig')
    saveas(fig,[ff2 filesep 'PeaksHeigthDist.png'],'png')
    
    % duration distribution
    figure
    hold on
    edges = [0:2:50];
    for kC = 1:ncond
        n = histcounts(D{kC},edges);
        xval = edges(1:end-1)+mean(diff(edges))/2;
        yval{kC} = n/length(D{kC})*100;%/Ttot{kt}*600;
        errval{kC} = sqrt(n)/length(D{kC})*100;%/Ttot{kt}*600;
        bar(xval,yval{kC},'facecolor',col{kC},'edgecolor','none','facealpha',0.5)
        errorbar(xval,yval{kC},errval{kC},'.','color',col{kC},'handlevisibility','off')
    end
    
    xlabel('Peaks duration (s)')
    ylabel('Proportion of all peaks (%)')

    legend(leg,'interpreter','none')
    fig = gcf;
    saveas(fig,[ff2 filesep 'PeaksDurDist.fig'],'fig')
    saveas(fig,[ff2 filesep 'PeaksDurDist.png'],'png')
     
end

if ACTI
%% Activity curves 

   
    figure
    hold on
    plot([0 0 0],[0.001 50 100],'--k','handlevisibility','off')
    plot([-500 1500],[alignlvl alignlvl]*100,'--k','handlevisibility','off')
    
    ActivityXglobal = 0:0.05:100;
    
    for kC= 1:ncond
        nn = length(ActivityCurves{kC});
        for kn = 1:nn
            ActivityY = ActivityCurves{kC}{kn}{1};
            ActivityX = ActivityCurves{kC}{kn}{2};
            Proba(kn,:) = interp1(ActivityX,ActivityY,ActivityXglobal,'linear',max(ActivityY));
            
        end
        
        mProba = mean(Proba,1);
        err = std(Proba,1)/sqrt(nn);
        
        plot(mProba,ActivityXglobal,'*','linewidth',3,'color',col{kC})
        
        plot(mProba-err,ActivityXglobal,'--','color',col{kC},'handlevisibility','off')
        plot(mProba+err,ActivityXglobal,'--','color',col{kC},'handlevisibility','off')
        
        clear Proba err
        
    end
    
    xlim([-200 800])
    xlabel('Height above cortex bottom thickness (nm)')
    ylabel('Fraction data points (%)')
    legend(leg)
    
    fig = gcf;
    saveas(fig,[ff2 filesep 'ProbDistLinV.fig'],'fig')
    saveas(fig,[ff2 filesep 'ProbDistLinV.png'],'png')    
    
end 

if TIME
%% Time analysis

    for kC = 1:ncond
        
        AllC{kC} = []; % for all correlation curves
        AllTau{kC} = []; % for all decorelation times
        
        for k = 1:length(Lags{kC})

           
                expfitfun = fittype('exp(-x/T)'); % exponential fit function

                expfit = fit(Lags{kC}{k},Corrs{kC}{k},expfitfun);
                
                Tau = coeffvalues(expfit);
                
                AllC{kC} = [AllC{kC}; Corrs{kC}{k}(1:Fs*MinTime+1)'];
                AllTau{kC} = [AllTau{kC} Tau];
                
            
        end
        
        MedCorr{kC} = mean([AllC{kC}]);
        MedLag{kC} =  0:(1/Fs):MinTime;
        ErrCorr{kC} = std(AllC{kC})/sqrt(length(Lags{kC}));
        
        
    end
    
    figure
    hold on
    
    for kC = 1:ncond
plot(MedLag{kC},MedCorr{kC},'color',col{kC},'linewidth',2)
plot(MedLag{kC},MedCorr{kC}-ErrCorr{kC},'--','linewidth',0.5,'color',col{kC},'handlevisibility','off')
plot(MedLag{kC},MedCorr{kC}+ErrCorr{kC},'--','linewidth',0.5,'color',col{kC},'handlevisibility','off')
    end

xlabel('Lag (sec)')
ylabel('AutoCorrelation')

xlim([0 350])
ylim([-0.3 1])
legend(leg)

    fig = gcf;
    saveas(fig,[ff2 filesep 'AutoCorr.fig'],'fig')
    saveas(fig,[ff2 filesep 'AutoCorr.png'],'png')

    
figure
hold on
for kC = 1:ncond

Spread_Julien(AllTau{kC},kC,0.9,col{kC},'o',5)
plot([kC-0.45 kC+0.45],[median(AllTau{kC}) median(AllTau{kC})],'r--','linewidth',2)

end

ypos = prctile([AllTau{:}],99.5); % for *** positioning
for kC = 1:ncond-1 % kC = 1:1 for comparing everything to the first condition only
        for kkt = kC+1:ncond
            pC = ranksum(Acti1090{kC},Acti1090{kkt});
            
            
            
            txtsup = ['p = ' num2str(pC,'%.3f')];
            
            if 5/100 > pC && pC > 1/100
                txtC = ['*  ' txtsup]  ;
            elseif 1/100 > pC && pC > 1/1000
                txtC = ['**  ' txtsup];
            elseif 1/1000 > pC
                txtC = ['***  ' txtsup];
            else
                txtC = ['NS  ' txtsup];
                
            end

             line([kC+0.1 kkt-0.1],[ypos ypos],'color','k', 'linewidth', 1.3);
            text((kC+kkt)/2,ypos + 2,txtC,'HorizontalAlignment','center','fontsize',11)
            
            ypos = ypos + 3;
    end
end
ax = gca;
ax.XTick = 1:ncond;
ax.XTickLabels = leg;

ylim([0 ypos+3])
ylabel('Decorrelation time (s)')

    fig = gcf;
    saveas(fig,[ff2 filesep 'DecorrelationTimes.fig'],'fig')
    saveas(fig,[ff2 filesep 'DecorrelationTimes.png'],'png')

end

end

function PlotMedP1090(MedCw,P10,P90,coord,col) % function for plotting all cells cortex thickness and fluctuation per cell

step = 0.9/(length(MedCw)+1);
x = (1:length(MedCw))*step-0.45;

Mat1= sortrows([MedCw' P10' P90']);
plot([coord+x; coord+x],[Mat1(:,2)'; Mat1(:,3)'],'-k'); 
plot(coord+x,Mat1(:,1),'ks','markerfacecolor',col)
plot(coord+x,Mat1(:,3),'kv','markerfacecolor',col)
plot(coord+x,Mat1(:,2),'k^','markerfacecolor',col)



end
