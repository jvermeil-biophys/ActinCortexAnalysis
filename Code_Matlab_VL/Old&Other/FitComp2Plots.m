function EPlot = FitComp2Plots(specif,PLOT,DISTCOMP,INDENT,col)

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultAxesFontSize',23)
% data folder
dataf   = 'D:\Data\MatFile\D2F';

list = ls(dataf);
list = list(3:end,:);

nfiles = size(list,1);


% sauvergarde des figures
sf = 'D:\Data\Figures';

datenow = datestr(now,'yy-mm-dd');

figfolder = [sf filesep datenow filesep 'Meca' filesep specif];

mkdir(figfolder)

ptr = find(ismember(specif,'_'));
savespecif = specif(ptr(end)+1:end);

% global var
GE = [];
GR2 = [];
GEsmall = [];
GEcomp = [];
GSStif = [];
GEpeak = [];
GHys = [];
GHysN = [];
Gddeb = [];
Gh0 = [];
Ghexp = [];
Gh0cort = [];
Ghexpcort = [];
Gcort = [];
Gdeltamax = [];
Gbasenoise = [];
Gforcenoise = [];
AllRamps = [];


Pos_ = find(specif == '_');
SpecifLeg = specif(Pos_(end) + 1:end);

for kf = 1:nfiles
    if contains(list(kf,:),specif)
        
        load([dataf filesep list(kf,:)])
        
        E = GlobalVar.E;
        R2 = GlobalVar.R2;
        Esmall = GlobalVar.Esmall;
        Ecomp = GlobalVar.Ecomp;
        StrStif = GlobalVar.StrStif;
        h0 = GlobalVar.h0;
        hexp = GlobalVar.hexp;
        h0cort = GlobalVar.h0cort;
        hexpcort = GlobalVar.hexpcort;
        Epeak = GlobalVar.Epeak;
        Hys = GlobalVar.Hysteresis;
        HysN = GlobalVar.HysteresisNorm;
        ddeb = GlobalVar.HeightBegining;
        cort = GlobalVar.Cortex;
        Deltamax = GlobalVar.Deltamax;
        ForceNoise = GlobalVar.ForceNoise;
        BaseNoise = GlobalVar.BaseNoise;
        
%         AllRamps = [AllRamps; AllRamp];
        

        
        %% Data globales
        

        GEsmall = [GEsmall Esmall(:)'];
        GE = [GE E(:)'];
        GR2 = [GR2 R2(:)'];
        GEcomp = [GEcomp Ecomp(:)'];
        GSStif = [GSStif StrStif(:)'];
        Gh0 = [Gh0 h0(:)'];
        Ghexp = [Ghexp hexp(:)'];
        Gh0cort = [Gh0cort h0cort(:)'];
        Ghexpcort = [Ghexpcort hexpcort(:)'];
        GEpeak = [GEpeak Epeak(:)'];
        GHys = [GHys Hys(:)'];
        GHysN = [GHysN HysN(:)'];
        Gddeb = [Gddeb ddeb(:)'];
        Gcort = [Gcort cort(:)'];
        Gdeltamax = [Gdeltamax 1000*Deltamax(:)'];  
        Gforcenoise = [Gforcenoise 1000*ForceNoise(:)'];  
        Gbasenoise = [Gbasenoise 1000*BaseNoise(:)'];   
        
        clear Ec Epeak EcC  EpeakC Hys HysN ddeb ptrnan
        
    end
end

ptrPlot = ~isnan(GE)&GE>0&GE<100000&GR2>0.85;

H0Plot = Ghexp(ptrPlot);
H0CortPlot = Ghexpcort(ptrPlot);
EPlot = GE(ptrPlot)/1000;
R2Plot = GR2(ptrPlot);





%% plot E vs H0
if PLOT   
 %% H0 exp
    [P,~] = polyfit(H0Plot,EPlot,1);
    SupFit = min(H0Plot):max(H0Plot);
    Efit = P(2)+SupFit*P(1);

    
    SSres = sum((polyval(P,H0Plot)-EPlot).^2);
    SStot = sum((EPlot-mean(EPlot)).^2);
    
    R2fit = 1-(SSres/SStot);
    
    figure
    hold on
    plot(H0Plot,EPlot,'*','markersize',5,'linewidth',1,'color',col)
    plot(SupFit,Efit,'linewidth',2,'color',col)
    xlabel('Epaisseur de départ (nm)')
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    ylabel('E (kPa)')
 
    title('Elinfit vs. H0exp')
    legend([SpecifLeg '- Esmall'],[SpecifLeg '- R2 : ' num2str(R2fit)]);

    saveas(gcf,[figfolder filesep 'EvH_' savespecif '.png'])
    saveas(gcf,[figfolder filesep 'EvH_' savespecif '.fig'])
    close

    figure(250)
    hold on
    LEG250 = legend;
    plot(SupFit,Efit,'linewidth',2,'color',col)
    LEG250.String{end} = [num2str(-P(1)) '-' num2str(R2fit)];
    xlabel('Epaisseur de départ (nm)')
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    ylabel('E (kPa)')
    
    title('Elinfit vs. H0exp')
    
    clear P SSres SStot R2fit
    
  %% H0 exp normalized by cort

  [P,~] = polyfit(H0CortPlot,EPlot,1);
    SupFit = min(H0CortPlot):max(H0CortPlot);
    Efit = P(2)+SupFit*P(1);

    
    SSres = sum((polyval(P,H0CortPlot)-EPlot).^2);
    SStot = sum((EPlot-mean(EPlot)).^2);
    
    R2fit = 1-(SSres/SStot);
    
    figure
    hold on
    plot(H0CortPlot,EPlot,'*','markersize',5,'linewidth',1,'color',col)
    plot(SupFit,Efit,'linewidth',2,'color',col)
    xlabel('Epaisseur de départ Norm (CU)')
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    ylabel('E (kPa)')
 
    title('Elinfit vs. H0expN')
    legend([SpecifLeg '- Esmall'],[SpecifLeg '- R2 : ' num2str(R2fit)]);
  saveas(gcf,[figfolder filesep 'ENvH_' savespecif '.png'])
  saveas(gcf,[figfolder filesep 'ENvH_' savespecif '.fig'])
    close

    figure(260)
    hold on
    LEG250 = legend;
    plot(SupFit,Efit,'linewidth',2,'color',col)
    LEG250.String{end} = [num2str(-P(1)) '-' num2str(R2fit)];
    xlabel('Epaisseur de départ Norm (CU)')
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    ylabel('E (kPa)')
    
    title('Elinfit vs. H0expN')
  
    
    clear P SSres SStot R2fit
end

%% Module distributions
if DISTCOMP
    
    figure(3000)
    hold on  
    LEG3 = legend;
    n = length(LEG3.String)+1;
    Spread_Julien(EPlot,n,0.7,col,'o')
    plot([n-0.35 n+0.35],...
        [nanmedian(EPlot) nanmedian(EPlot)],...
        'k','linewidth',2,'handlevisibility','off')
    title('E distributions')
    ylabel('E (kPa)')
    LEG3.String{end} = SpecifLeg;
    legend('location','best')
    

    
end

if INDENT
%% courbe moyenne des compressions
% 
% 
% Y = mean(AllRamps*1000,1);
% err = std(AllRamps*1000,1)/sqrt(size(AllRamps,1));
% 
% polygonX = [1:length(Y) length(Y):-1:1];
% polygonY = [Y+err Y(end:-1:1)-err(end:-1:1)];
% 
% figure(9874)
% hold on
% plot(Y,'linewidth',3,'color',col)
% patch(polygonX,polygonY,col,'facealpha',0.25,'edgecolor','none','handlevisibility','off')

%% enfoncement maximum par compresion et bruits de mesure


    figure(4000)
    hold on
    
    LEG4 = legend;
    n = length(LEG4.String)+1;
    
    spreadval = 0.7;
    
    Spread_Julien(Gdeltamax,n,2,spreadval,col)
    plot([n-spreadval/2 n+spreadval/2],[nanmedian(Gdeltamax) nanmedian(Gdeltamax)],...
        'k','linewidth',2,'handlevisibility','off')
%     plot([n-spreadval/2 n+spreadval/2],[nanmean(Gdeltamax) nanmean(Gdeltamax)],...
%         'r','linewidth',2,'handlevisibility','off')
    
    title('Max indentation')
    ylabel('$\delta_{max}$ (nm)')
    LEG4.String{end} = SpecifLeg;
    
    
        figure(5000)
    hold on
    
    LEG5 = legend;
    n = length(LEG5.String)+1;
   
    Spread_Julien(Gforcenoise,n,0.1,spreadval,col)
    plot([n-spreadval/2 n+spreadval/2],[nanmedian(Gforcenoise) nanmedian(Gforcenoise)],...
        'k','linewidth',2,'handlevisibility','off')
%     plot([n-spreadval/2 n+spreadval/2],[nanmean(Gforcenoise) nanmean(Gforcenoise)],...
%         'r','linewidth',2,'handlevisibility','off')
    
    title('Force Noise')
    ylabel('$Force noise$ (nm)')
    LEG5.String{end} = SpecifLeg;
    
    
        figure(6000)
    hold on
    
    LEG6 = legend;
    n = length(LEG6.String)+1;
    
    spreadval = 0.7;
    
    Spread_Julien(Gbasenoise,n,0.1,spreadval,col)
    plot([n-spreadval/2 n+spreadval/2],[nanmedian(Gbasenoise) nanmedian(Gbasenoise)],...
        'k','linewidth',2,'handlevisibility','off')
%     plot([n-spreadval/2 n+spreadval/2],[nanmean(Gbasenoise) nanmean(Gbasenoise)],...
%         'r','linewidth',2,'handlevisibility','off')
    
    title('Base Noise')
    ylabel('$Base noise$ (nm)')
    LEG6.String{end} = SpecifLeg;
    
end


end
