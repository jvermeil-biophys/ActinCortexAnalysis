function StressStrain2Meca_Rheology(specif,tags,PLOT,resfolder,figfolder)
%% INIT

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultAxesFontSize',17)

symbollist = {'o' 'd' 's' '>' '<' 'v'};
colorlist = {'r' 'b' 'g' 'm' 'c' 'k'};

%% params and settings
% dossier de données et sauvegarde des données
path = [resfolder filesep 'D2SS'];
sf   = [resfolder filesep 'SS2M'];

datenow = datestr(now,'yy-mm-dd');

mkdir(sf)

savename=['SS2M_' [tags{1:end}] '_' specif];
savenamepart = [savename '_Part'];


if PLOT
    % sauvergarde des figures
    sfig = figfolder;
    sff = [sfig filesep datenow filesep 'Meca' filesep specif filesep 'Res'];
    mkdir(sff)
end


%% Retrieve datafile and analyze
listm = ls(path);

for ks = 1:size(listm,1)
    list{ks} = listm(ks,:);
end

loadnamepart=[tags{:} '_' specif '.mat'];

ptrlist = find(contains(list,loadnamepart));

cd(path)

LP = length(ptrlist);

for kl = LP
    
    loadname = list{ptrlist(kl)};
    
    fprintf(['Chargement du fichier ' loadname '\n\n'])
    load(loadname);
    
    LM = length(MF);
    
    MeanTimeCos = cell(LP,LM);
    MeanTimeFFT = cell(LP,LM);
    StorageE = cell(LP,LM);
    LossE = cell(LP,LM);
    ThickLoc = cell(LP,LM);
    DephasFFT = cell(LP,LM);
    DephasCos = cell(LP,LM);
    ModuleCos = cell(LP,LM);
    
    
    FreqExp = cell(LP,LM);
    
    NameList = {};
    
    PlotColor = cell(LP,LM);
    
    for kc = 1:LM
        
        %% récuperation des données
        
        Name = MF{kc}.CellName;
        
        if ~ismember(Name,NameList)
            NameList = {NameList{:} Name};
            PlotColor{kl,kc} = colorlist{1};
            colorlist = {colorlist{2:end}};
        else
            NameList = {NameList{:} Name};
            ptr = find(strcmp(NameList,Name),1);
            PlotColor{kl,kc} = PlotColor{ptr};
        end
        
        FreqExp{kl,kc} = MF{kc}.ExpFreq;
        FreqEch = MF{kc}.FreqEch;
        
        Champ = MF{kc}.field;
        
        ExpType = MF{kc}.Exp;        
        
        Sigma = MF{kc}.Sigma;
        Epsilon = MF{kc}.Epsilon;
        

        BaseDist = MF{kc}.BaseSignal;
        
        T = MF{kc}.time;
        
        
   if FreqExp{kl,kc} > 0.5 && FreqExp{kl,kc} < 2
        
        %% Calcul des zones ~constante en epaisseur

        
        k1 = 1;
        k2 = 1;
        partsList ={};
        
        while k1 < length(BaseDist)
       
            stdtmp = 0;
            meantmp = 1;
            
            while stdtmp < 50 && k2 < length(BaseDist)
                
                k2 = k2+1;
                
                meantmp = mean(BaseDist(k1:k2));
                stdtmp = max(BaseDist(k1:k2))-min(BaseDist(k1:k2));  
                
            end
            
            partsList{end+1} = k1:k2-1;
            
            k1 = k2;
            
        end
%         
%         figure
%         subplot(311)
%         hold on
%         title('Oscilation amplitude < 50 nm')
%         ylabel('Base Thick (nm)')
%         for k = 1:length(partsList)
%             
%             plot(T(partsList{k}),BaseDist(partsList{k}),'linewidth',2)
%             
%         end
%         
%         subplot(312)
%         hold on
%         ylabel('Stress (Pa)')
%         for k = 1:length(partsList)
%             
%             plot(T(partsList{k}),SigmaFilt(partsList{k}),'linewidth',2)
%             
%         end
%         
%         subplot(313)
%         hold on
%         ylabel('Strain (\%)')
%         xlabel('Time (s)')
%         for k = 1:length(partsList)
%             
%             plot(T(partsList{k}),EpsilonFilt(partsList{k}),'linewidth',2)
%             
%         end
%         


        %% Calcul module complexe et déphasage  sur 3 periodes avec Fourier
        % avec fourier         
        
        % calcul sur 3 periode glissante
        
        Texp = 1/FreqExp{kl,kc}; % periode en s
        
        Ncycle = T(end)/Texp;
        
        if Ncycle > 5
            
            FitLength = floor(3*Texp*FreqEch);
            FitShift = floor(Texp*FreqEch*3);
            Nfit = floor(Ncycle/3);
            
            MeanTimeFFT{kl,kc} = zeros(1,Nfit);
            
            ThickLoc{kl,kc} = zeros(1,Nfit);
            
            StorageE{kl,kc} = zeros(1,Nfit);
            LossE{kl,kc} = zeros(1,Nfit);
            DephasFFT{kl,kc} = zeros(1,Nfit);
            
            for k = 1:Nfit
                
                ChampPart = Champ((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                Tpart = T((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                
                SigPartToFT = Sigma((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                EpsPartToFT = Epsilon((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
             
                
                FtSigPart = fft(SigPartToFT);
                FtEpsPart = fft(EpsPartToFT);
                FtChampPart = fft(ChampPart-mean(ChampPart));
                
                FreqsPart = (0:length(FtSigPart)-1)*FreqEch/length(FtSigPart);
                
                
                FreqExpLoc = FreqsPart(find(abs(FtChampPart)>max(abs(FtChampPart)/2),1));
                
                
                FreqOscPtr = abs(FreqsPart - FreqExpLoc) < 1e-5;
                
                Ecmplx = FtSigPart(FreqOscPtr)./FtEpsPart(FreqOscPtr);
                
                MeanTimeFFT{kl,kc}(k) = mean(Tpart);
                
                StorageE{kl,kc}(k) = real(Ecmplx);
                LossE{kl,kc}(k) = abs(imag(Ecmplx));
                ThickLoc{kl,kc}(k) = mean(BaseDist((k-1)*FitShift+1:(k-1)*FitShift+FitLength));
                
                
                DephasFFT{kl,kc}(k) = mod(angle(FtSigPart(FreqOscPtr))-angle(FtEpsPart(FreqOscPtr))+1,2*pi)-1;
                
            end
        end
        
        
        %% Fit d'un cosinus sur 3 periode
              
        % filtrage sigma et epsilon
        FtSigma = fft(Sigma);
        FtEpsilon  = fft(Epsilon);
        
        
        FreqsFilt = (0:length(FtSigma)-1)*FreqEch/length(FtSigma);
                
        % filt BF
        FtSigma((FreqsFilt<0.5 | FreqsFilt>(FreqEch-0.5))) = 0;
        FtEpsilon((FreqsFilt<0.5 | FreqsFilt>(FreqEch-0.5))) = 0;
        
       
        % filt HFBF
%         FtSigma((FreqsFilt<0.5 | FreqsFilt>(FreqEch-0.5))|(FreqsFilt>2 & FreqsFilt<(FreqEch-2))) = 0;
%         FtEpsilon((FreqsFilt<0.5 | FreqsFilt>(FreqEch-0.5))|(FreqsFilt>2 & FreqsFilt<(FreqEch-2))) = 0;
        
%         


        SigmaFilt = ifft(FtSigma);
        EpsilonFilt = ifft(FtEpsilon);
        
        % no filt
%         SigmaFilt = Sigma;
%         EpsilonFilt = Epsilon;
        
        
        % fit du cosinus
                
        if Ncycle > 5
            
            FitLength = floor(3*Texp*FreqEch);
            FitShift = floor(3*Texp*FreqEch);
            Nfit = floor((Ncycle/3));
            
            MeanTimeCos{kl,kc} = zeros(1,Nfit);
            
            ThickLoc{kl,kc} = zeros(1,Nfit);

            DephasCos{kl,kc} = zeros(1,Nfit);
            ModuleCos{kl,kc} = zeros(1,Nfit);
            
            
            for k = 1:Nfit
                
                Tpart = T((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                
                

                SigFilPartToFit = SigmaFilt((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                EpsFilPartToFit = EpsilonFilt((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                
                
  
                
                MeanTimeCos{kl,kc}(k) = mean(Tpart);

                ThickLoc{kl,kc}(k) = mean(BaseDist((k-1)*FitShift+1:(k-1)*FitShift+FitLength));
                
                

                
               
[AmpEpsFilPartFit,PhaseEpsFilPartFit,EpsFilPartFit,~,OffsetEps] = FitCosFunc(EpsFilPartToFit,Tpart,FreqExp{kl,kc});        

                
           
                
                
[AmpSigFilPartFit,PhaseSigFilPartFit,SigFilPartFit,~,OffsetSig] = FitCosFunc(SigFilPartToFit,Tpart,FreqExp{kl,kc});
                

                
                
                DephasCos{kl,kc}(k) = mod(PhaseSigFilPartFit-PhaseEpsFilPartFit+1,2*pi)-1;
                ModuleCos{kl,kc}(k) = AmpSigFilPartFit/AmpEpsFilPartFit;
  
% %                                 
%                 figure
%                 hold on
%                 yyaxis right
%                 plot(Tpart,EpsFilPartToFit,'-*')
%                 plot(Tpart,EpsFilPartFit(Tpart),'linewidth',2)
%                 ylim([-AmpEpsFilPartFit AmpEpsFilPartFit]*2+OffsetEps)
%                 xlim([Tpart(1) Tpart(end)])
%                 ylabel('Epsilon')
%                 xlabel('Time')
%                 yyaxis left
%                 plot(Tpart,SigFilPartToFit,'-o')
%                 plot(Tpart,SigFilPartFit(Tpart),'linewidth',2)
%                 ylim([-AmpSigFilPartFit AmpSigFilPartFit]*2+OffsetSig)
%                 xlim([Tpart(1) Tpart(end)])
%                 ylabel('Sigma')
%                 
%                 title(['FiltBFHF(0.5-2). Dephas = ' num2str(DephasCos{kl,kc}(k)) ])
                
%                 fig = gcf;
%                 saveas(fig,['C:\Users\PMMH\Desktop\DataValentin\Figures\MécaDictyFigs\PartialCosFits\NoFilt_' Name '_' ExpType '-' num2str(k) '.png'])
%                 close

                
            end
            
        end
        
        
    end
    
    end


end
ThickLoc = {ThickLoc{:}};
StorageE = {StorageE{:}};
LossE = {LossE{:}};
DephasFFT = {DephasFFT{:}};
DephasCos = {DephasCos{:}};
ModuleCos = {ModuleCos{:}};
PlotColor = {PlotColor{:}};


% pour plots

FreqExpR = precisround([FreqExp{:}],0.1);

FreqList = unique(FreqExpR);

FreqList = 1;

nfreq = length(FreqList);


%% Plot fourier

if PLOT 
%     
% figure(1)
% hold on
% title('Storage modulus versus thickness (Fourier)')
% xlabel('Thickness (nm)')
% ylabel('Storage modulus (kPa)')
% ylim([0 18])
% 
% figure(2)
% hold on
% title('Loss modulus versus thickness (Fourier)')
% xlabel('Thickness (nm)')
% ylabel('Loss modulus (kPa)')
% ylim([0 7])

figure(3)
hold on
title('Phase Diff versus time (Fourier)')
xlabel('Time (s)')
ylabel('Phase Diff (rad)')
ylim([-0.5 pi/2])

% 
% figure(4)
% hold on
% title('Storage modulus vs Oscillation frequency (Fourier)')
% ylabel('E'' (kPa)')
% ax = gca;
% ax.XTick = 1:nfreq;
% xlabel('Freq (Hz)')
% xlim([0.5 nfreq+0.5])
% ylim([0 18])

% figure(5)
% hold on
% title('Loss modulus vs Oscillation frequency (Fourier)')
% ylabel('E'''' (kPa?)')
% ax = gca;
% ax.XTick = 1:nfreq;
% xlabel('Freq (Hz)')
% xlim([0.5 nfreq+0.5])
% ylim([0 7])
% 
% figure(6)
% hold on
% title('Phase Diff vs Oscillation frequency (Fourier)')
% ylabel('Phase Diff (rad)')
% ax = gca;
% ax.XTick = 1:nfreq;
% xlabel('Freq (Hz)')
% xlim([0.5 nfreq+0.5])

figure(8000)
hold on
title('DephasFFT vs DephasCos FiltBF')
plot([-0.5 pi/2],[-0.5 pi/2],'k-')
xlabel('DephasCosFiltBF')
ylabel('DephasFFTNoFilt')
xlim([-0.5 pi/2])
ylim(xlim)

for k = 1:nfreq
    ptr = FreqExpR == FreqList(k);
    TimeFFTptr = {MeanTimeFFT{ptr}};
    Storageptr = {StorageE{ptr}};
    Lossptr = {LossE{ptr}};
    ThickLocptr = {ThickLoc{ptr}};
    Dephasptr = {DephasFFT{ptr}};
    DephasCptr = {DephasCos{ptr}};



    PlotColorptr = {PlotColor{ptr}};
    
    for kt = 1:length(ThickLocptr)
        
%         figure(1)
%         hold on
%         plot(ThickLocptr{kt},Storageptr{kt}/1000,'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
%         
%         figure(2)
%         hold on
%         plot(ThickLocptr{kt},Lossptr{kt}/1000,'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
%         
        figure(3)
        hold on
        plot(TimeFFTptr{kt},Dephasptr{kt},'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        
        figure(8000)
        hold on
        plot(DephasCptr{kt},Dephasptr{kt},'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        
    end
%     
%     figure(4)
%     hold on
%     Spread_Julien([Storageptr{:}]/1000,k,0.8,'r',symbollist{k},0.4)
%     plot([k-0.45 k+0.45],[mean([Storageptr{:}]/1000) mean([Storageptr{:}]/1000)],'--g','linewidth',2)
%     ax = gca;
%     ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz'];
%     
%     figure(5)
%     hold on
%     Spread_Julien([Lossptr{:}]/1000,k,0.8,'b',symbollist{k},0.15)
%     plot([k-0.45 k+0.45],[mean([Lossptr{:}]/1000) mean([Lossptr{:}]/1000)],'--g','linewidth',2)
%     ax = gca;
%     ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz'];
%     
%     figure(6)
%     hold on
%     Spread_Julien([Dephasptr{:}],k,0.8,'b',symbollist{k},0.15)
%     plot([k-0.45 k+0.45],[mean([Dephasptr{:}]) mean([Dephasptr{:}])],'--g','linewidth',2)
%     ax = gca;
%     ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz'];

end

end

%% Plot FitCos
 if PLOT
%      
% figure(11)
% hold on
% title('Elastic modulus versus thickness (FitCos)')
% xlabel('Thickness (nm)')
% ylabel('Elastic modulus (kPa)')
% ylim([0 18])

figure(12)
subplot(nfreq,1,1)
hold on
title('Dephas versus time (FitCos)')
subplot(nfreq,1,nfreq)
hold on
xlabel('Time(s)')
ylabel('Phase Diff (rad)')

% 
% figure(13)
% hold on
% title('Elastic modulus vs Oscillation frequency (FitCos)')
% ylabel('E'' (kPa)')
% ax = gca;
% ax.XTick = 1:nfreq;
% xlabel('Freq (Hz)')
% xlim([0.5 nfreq+0.5])
% ylim([0 18])

% figure(14)
% hold on
% title('Phase Diff vs Oscillation frequency (FitCos)')
% ylabel('Phase Diff (rad)')
% ax = gca;
% ax.XTick = 1:nfreq;
% xlabel('Freq (Hz)')
% xlim([0.5 nfreq+0.5])

% 
% figure(15)
% hold on
% title('Elastic modulus versus time (FitCos)')
% xlabel('Average time of sample')
% ylabel('Elastic modulus (kPa)')
% ylim([0 18])


for k = 1:nfreq
    ptr = FreqExpR == FreqList(k);
    TimeCosptr = {MeanTimeCos{ptr}};
    Moduleptr = {ModuleCos{ptr}};
    ThickLocptr = {ThickLoc{ptr}};
    Dephasptr = {DephasCos{ptr}};
    PlotColorptr = {PlotColor{ptr}};
    
    for kt = 1:length(ThickLocptr)
        
%         figure(11)
%         hold on
%         plot(ThickLocptr{kt},Moduleptr{kt}/1000,'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
%         

% ptrdeph = Dephasptr{kt}>pi;
% Dephasptr{kt}(ptrdeph) = Dephasptr{kt}(ptrdeph)-2*pi;

        figure(12)
        subplot(nfreq,1,k)
        hold on
        plot(TimeCosptr{kt},Dephasptr{kt},'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        ylim([-0.5 pi/2])
        
%         figure(15)
%         hold on
%         plot(Timeptr{kt},Moduleptr{kt}/1000,'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})

    end
    
%     figure(13)
%     hold on
%     Spread_Julien([Moduleptr{:}]/1000,k,0.8,'r',symbollist{k},0.4)
%     plot([k-0.45 k+0.45],[mean([Moduleptr{:}]/1000) mean([Moduleptr{:}]/1000)],'--g','linewidth',2)
%     ax = gca;
%     ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz']
%     
%     figure(14)
%     hold on
%     Spread_Julien([Dephasptr{:}],k,0.8,'b',symbollist{k},0.15)
%     plot([k-0.45 k+0.45],[mean([Dephasptr{:}]) mean([Dephasptr{:}])],'--g','linewidth',2)
%     ax = gca;
%     ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz']
end

end





end

function [AmpFit,PhaseFit,Fit,R2,Offset] = FitCosFunc(Data,Tdata,Freq)
 
        fitcos = fittype('A*cos(2*pi*f*x + p) + M','problem',{'f'},'coefficient', {'A', 'p', 'M'});
        
        AmpStart = (max(Data)-min(Data))/2;
        
        [Fit,gof] = fit(Tdata,Data,fitcos,'problem',{Freq}, ...
            'start',[AmpStart 0 mean(Data)],'upper',[Inf pi max(Data)],'lower',[0 -pi min(Data)]);
        
        Fitc = coeffvalues(Fit);
        
        AmpFit = Fitc(1);

        PhaseFit = Fitc(2);
        
        Offset = Fitc(3);

        R2 = gof.rsquare;

end

function out = precisround(n,precis)

out = round(n./precis).*precis;

end