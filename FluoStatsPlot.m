function FluoStatsPlot(specifs,BAR,PROB,DIST)

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)

%data folder
datf   = 'D:\Data\MatFile\MA';

nspe = length(specifs);

Ttot = zeros(1,nspe);
HInt = cell(1,nspe);
HDen = cell(1,nspe);
HSiz = cell(1,nspe);

% Couleurs par drogues
Cy = [200 200  0  ]./255;

Cs  = [0   255  90]./255;

Cm = [200 100 80  ]./255;
Cc = [255 20   0  ]./255;

Cd  = [160 160 160]./255;

Cl = [110 50 120]./255;


% figures saving location
ff = 'D:\Data\Figures';

date = datestr(now,'yy-mm-dd');
ff = [ff filesep date filesep 'FluoPeaks'];
df = ['C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Figure' filesep 'FluoPeaks'];
mkdir(ff)
mkdir(df)

col = cell(nspe);

for ks = 1:nspe
    
    if strcmp('Ctrl',specifs{ks})
        tag1 = 'd';
    elseif strcmp('CK666',specifs{ks})
        tag1 = 'c';
    elseif strcmp('Y27',specifs{ks})
        tag1 = 'y';
    elseif strcmp('Smifh2',specifs{ks})
        tag1 = 's';
    elseif strcmp('LatA',specifs{ks})
        tag1 = 'l';
    end
    eval(['col{ks} = C' tag1 ';']);
end

% récupération des données
for ks = 1:nspe
    load([datf filesep 'MA_' specifs{ks}])
    ProbInts{ks} = sort(Prob.allints);
    ProbDens{ks} = sort(Prob.alldens);
    ProbSizs{ks} = sort(Prob.allsizs);
    
    P{ks} = (length(ProbInts{ks}):-1:1)/length(ProbInts{ks});
    
    for km = 1:length(MA)
        name = MA{km}.name;
        for kp = 1:length(MA{km}.pos)
            
            Ttot(ks) = Ttot(ks) + MA{km}.time{kp}(end);
            
            HInt{ks} = [HInt{ks} MA{km}.peaksInt{kp}(:,4)'];
            HDen{ks} = [HDen{ks} MA{km}.peaksDen{kp}(:,4)'];
            HSiz{ks} = [HSiz{ks} MA{km}.peaksSiz{kp}(:,1)'];
        end
    end
end

if BAR
    % %% Bar graph Int
    % edges = [18 62 180 Inf];
    %
    % edgesloc = 1:length(edges);
    %
    % Bw = 0.9/nspe;
    %
    % figure
    % hold on
    % title('Frequencies of events per intensity')
    % for ks = 1:nspe
    %
    %
    %
    %     for i = 1:length(edges)-1
    %         counts(i) = sum((HInt{ks}>edges(i))&(HInt{ks}<edges(i+1)));
    %     end
    %
    %     str = strtrim(cellstr(num2str(round(counts)'))');
    %
    %     err = sqrt(counts);
    %     counts = counts/Ttot(ks)*60;
    %     err = err/Ttot(ks)*60;
    %
    %
    %     for kpl = 1:length(edges)-1
    %         subplot(1,length(edges)-1,kpl)
    %         hold on
    %         if ks == 1
    %             Scale(kpl) = counts(kpl)*1.8;
    %         end
    %         bar(1+(ks-1)*Bw,counts(kpl),1/(nspe+1),'facecolor',col{ks},'edgecolor','none')
    %         errorbar(1+(ks-1)*Bw,counts(kpl),err(kpl),'k.','handlevisibility','off')
    %         xlim([0.8 2])
    %         ax = gca;
    %         ax.XTick = 1.4;
    %         ax.XTickLabels = [num2str(edges(kpl)') '-' num2str(edges(kpl+1))];
    %         clear ax
    %     end
    %
    %
    %     clear err counts
    %
    %
    %
    %
    % end
    %
    % subplot(1,length(edges)-1,1)
    % ylabel('F (pic/min)')
    % ylim([0 Scale(1)])
    %
    %
    % subplot(1,length(edges)-1,round((length(edges)-1)/2))
    % xlabel(['Intensity (% of mean)'])
    % ylim([0 Scale(2)])
    %
    % subplot(1,length(edges)-1,3)
    % ylim([0 Scale(3)])
    % errorbar(-20,0,0,'k.')
    % leg = specifs;
    % leg{end+1} = 'errors are sqrt(N)';
    % legpos = [0.5 0.975 0 0];
    % legend(leg,'position',legpos,'Orientation','horizontal','fontsize',18)
    % fig = gcf;
    %
    %
    %
    % saveas(fig,[ff filesep 'BarFluoInt.fig'],'fig')
    % saveas(fig,[ff filesep 'BarFluoInt.png'],'png')
    % saveas(fig,[df filesep 'BarFluoInt.fig'],'fig')
    % saveas(fig,[df filesep 'BarFluoInt.png'],'png')
    %
    %
    %
    % close(fig)
    %
    % %% Bar graph Den
    % edges = [4 13 56 Inf];
    %
    % edgesloc = 1:length(edges);
    %
    % Bw = 0.9/nspe;
    %
    % figure
    % hold on
    % title('Frequencies of events per density')
    % for ks = 1:nspe
    %
    %
    %
    %     for i = 1:length(edges)-1
    %         counts(i) = sum((HDen{ks}>edges(i))&(HDen{ks}<edges(i+1)));
    %     end
    %
    %     str = strtrim(cellstr(num2str(round(counts)'))');
    %
    %     err = sqrt(counts);
    %     counts = counts/Ttot(ks)*60;
    %     err = err/Ttot(ks)*60;
    %
    %
    %     for kpl = 1:length(edges)-1
    %         subplot(1,length(edges)-1,kpl)
    %         hold on
    %         if ks == 1
    %             Scale(kpl) = counts(kpl)*1.8;
    %         end
    %         bar(1+(ks-1)*Bw,counts(kpl),1/(nspe+1),'facecolor',col{ks},'edgecolor','none')
    %         errorbar(1+(ks-1)*Bw,counts(kpl),err(kpl),'k.','handlevisibility','off')
    %         xlim([0.8 2])
    %         ax = gca;
    %         ax.XTick = 1.4;
    %         ax.XTickLabels = [num2str(edges(kpl)') '-' num2str(edges(kpl+1))];
    %         clear ax
    %     end
    %
    %
    %     clear err counts
    %
    %
    %
    %
    % end
    %
    % for kpl = 1:length(edges)-1
    %     subplot(1,length(edges)-1,kpl)
    %     ylim([0 Scale(kpl)])
    % end
    %
    % subplot(1,length(edges)-1,1)
    % ylabel('F (pic/min)')
    %
    %
    % subplot(1,length(edges)-1,round((length(edges)-1)/2))
    % xlabel(['Density'])
    %
    % subplot(1,length(edges)-1,3)
    % errorbar(-20,0,0,'k.')
    % leg = specifs;
    % leg{end+1} = 'errors are sqrt(N)';
    % legpos = [0.5 0.975 0 0];
    % legend(leg,'position',legpos,'Orientation','horizontal','fontsize',18)
    % fig = gcf;
    %
    %
    % saveas(fig,[ff filesep 'BarFluoDen.fig'],'fig')
    % saveas(fig,[ff filesep 'BarFluoDen.png'],'png')
    % saveas(fig,[df filesep 'BarFluoDen.fig'],'fig')
    % saveas(fig,[df filesep 'BarFluoDen.png'],'png')
    %
    %
    %
    % close(fig)
    
    %% Bar graph Siz
    edges = [0 1500 3000 3800 Inf];
    
    edgesloc = 1:length(edges);
    
    Bw = 0.9/nspe;
    
    figure
    hold on
    title('Frequencies of events per Size')
    for ks = 1:nspe
        
        
        
        for i = 1:length(edges)-1
            counts(i) = sum((HSiz{ks}>edges(i))&(HSiz{ks}<edges(i+1)));
        end
        
        str = strtrim(cellstr(num2str(round(counts)'))');
        
        err = sqrt(counts);
        counts = counts/Ttot(ks)*60;
        err = err/Ttot(ks)*60;
        
        
        for kpl = 1:length(edges)-1
            subplot(1,length(edges)-1,kpl)
            hold on
            if ks == 1
                Scale(kpl) = counts(kpl)*1.8;
            end
            bar(1+(ks-1)*Bw,counts(kpl),1/(nspe+1),'facecolor',col{ks},'edgecolor','none')
            errorbar(1+(ks-1)*Bw,counts(kpl),err(kpl),'k.','handlevisibility','off')
            xlim([0.8 2])
            ax = gca;
            ax.XTick = 1.4;
            ax.XTickLabels = [num2str(edges(kpl)') '-' num2str(edges(kpl+1))];
            clear ax
        end
        
        
        clear err counts
        
        
        
        
    end
    
    for kpl = 1:length(edges)-1
        subplot(1,length(edges)-1,kpl)
        ylim([0 Scale(kpl)])
    end
    
    subplot(1,length(edges)-1,1)
    ylabel('F (pic/min)')
    
    subplot(1,length(edges)-1,floor(length(edges)/2))
    xlabel(['Size (nm)'])
    
    subplot(1,length(edges)-1,length(edges)-1)
    errorbar(-20,0,0,'k.')
    leg = specifs;
    leg{end+1} = 'errors are sqrt(N)';
    legpos = [0.5 0.975 0 0];
    legend(leg,'position',legpos,'Orientation','horizontal','fontsize',18)
    fig = gcf;
    
    
    saveas(fig,[ff filesep 'BarFluoSiz.fig'],'fig')
    saveas(fig,[ff filesep 'BarFluoSiz.png'],'png')
    saveas(fig,[df filesep 'BarFluoSiz.fig'],'fig')
    saveas(fig,[df filesep 'BarFluoSiz.png'],'png')
    
    
    
    % close(fig)
    
end

if PROB
    %% Prob Curve Int
    
    figure
    hold on
    
    
    for ks = 1:nspe
        plot(ProbInts{ks},P{ks}...
            ,'*','linewidth',2,'color',col{ks})
    end
    
    xlabel('Intensity above mean')
    ylabel('P(h)')
    legend(specifs,'location','best')
    
    set(gca,'YScale','log')
    
    xl = xlim;
    yl = ylim;
    
    xlim([-200 1000])
    ylim([0.0001 yl(2)])
    
    xpatch = [xl(1) xl(2) xl(2) xl(1)];
    ypatch = [yl(1) yl(1) 0.01 0.01];
    patch(xpatch,ypatch,[0.3 0.3 0.3],'facealpha','interp','facevertexalphadata',[0.5 0.5 0 0]',...
        'edgecolor','none','handlevisibility','off')
    
    
    ypatch = [yl(1) yl(1) 0.001 0.001];
    patch(xpatch,ypatch,[0.3 0.3 0.3],'facealpha','interp','facevertexalphadata',[0.5 0.5 0 0]',...
        'edgecolor','none','handlevisibility','off')
    
    fig = gcf;
    
    
    saveas(fig,[ff filesep 'ProbaFluoInt.fig'],'fig')
    saveas(fig,[ff filesep 'ProbaFluoInt.png'],'png')
    saveas(fig,[df filesep 'ProbaFluoInt.fig'],'fig')
    saveas(fig,[df filesep 'ProbaFluoInt.png'],'png')
    
    close(fig)
    %% Prob Curve Den
    
    figure
    hold on
    
    
    for ks = 1:nspe
        plot(ProbDens{ks},P{ks}...
            ,'*','linewidth',2,'color',col{ks})
    end
    
    xlabel('Density above threshold')
    ylabel('P(h)')
    legend(specifs,'location','best')
    
    set(gca,'YScale','log')
    
    xl = xlim;
    yl = ylim;
    
    xlim([-10 100])
    ylim([0.0001 yl(2)])
    
    xpatch = [xl(1) xl(2) xl(2) xl(1)];
    ypatch = [yl(1) yl(1) 0.01 0.01];
    patch(xpatch,ypatch,[0.3 0.3 0.3],'facealpha','interp','facevertexalphadata',[0.5 0.5 0 0]',...
        'edgecolor','none','handlevisibility','off')
    
    
    ypatch = [yl(1) yl(1) 0.001 0.001];
    patch(xpatch,ypatch,[0.3 0.3 0.3],'facealpha','interp','facevertexalphadata',[0.5 0.5 0 0]',...
        'edgecolor','none','handlevisibility','off')
    
    fig = gcf;
    
    
    saveas(fig,[ff filesep 'ProbaFluoDen.fig'],'fig')
    saveas(fig,[ff filesep 'ProbaFluoDen.png'],'png')
    saveas(fig,[df filesep 'ProbaFluoDen.fig'],'fig')
    saveas(fig,[df filesep 'ProbaFluoDen.png'],'png')
    
    close(fig)
    %% Prob Curve Siz
    
    figure
    hold on
    
    
    for ks = 1:nspe
        plot(ProbSizs{ks},P{ks}...
            ,'*','linewidth',2,'color',col{ks})
    end
    
    xlabel('Size of protrusion (px)')
    ylabel('P(h)')
    legend(specifs,'location','best')
    
    set(gca,'YScale','log')
    
    xl = xlim;
    yl = ylim;
    
    xlim([0 700])
    ylim([0.0001 yl(2)])
    
    xpatch = [xl(1) xl(2) xl(2) xl(1)];
    ypatch = [yl(1) yl(1) 0.01 0.01];
    patch(xpatch,ypatch,[0 0 0],'facealpha','interp','facevertexalphadata',[0.5 0.5 0 0]',...
        'edgecolor','none','handlevisibility','off')
    
    
    ypatch = [yl(1) yl(1) 0.001 0.001];
    patch(xpatch,ypatch,[0 0 0],'facealpha','interp','facevertexalphadata',[0.5 0.5 0 0]',...
        'edgecolor','none','handlevisibility','off')
    
    fig = gcf;
    
    
    saveas(fig,[ff filesep 'ProbaFluoSiz.fig'],'fig')
    saveas(fig,[ff filesep 'ProbaFluoSiz.png'],'png')
    saveas(fig,[df filesep 'ProbaFluoSiz.fig'],'fig')
    saveas(fig,[df filesep 'ProbaFluoSiz.png'],'png')
    
    close(fig)
end

if DIST
    
    figure
    hold on
    edges =([0:500:5000 10000]);
    for ks = 1:nspe
        n = histcounts(HSiz{ks},edges);
        xval = edges(1:end-1)+diff(edges)/2;
        yval{ks} = n/Ttot(ks)*600;
        bar(xval,yval{ks},'facecolor',col{ks},'edgecolor','none','facealpha',0.5)
    end
    xlabel('Peak size (nm)')
    ylabel('Frequency (peak/10minutes)')
    ylim([0 12])
    fig=gcf;
    saveas(fig,[ff filesep 'SizePeaksDist' [specifs{:}] '.png'],'png')
    saveas(fig,[ff filesep 'SizePeaksDist' [specifs{:}] '.fig'],'fig')
    
end

end

