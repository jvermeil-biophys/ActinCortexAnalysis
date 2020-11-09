function Fourier2Meca_Rheology(specif,tags,PLOT,resfolder,figfolder)
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
path = [resfolder filesep 'D2FR'];
sf   = [resfolder filesep 'F2MR'];

datenow = datestr(now,'yy-mm-dd');

mkdir(sf)

savename=['F2MR_' [tags{1:end}] '_' specif];
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
    
    MeanTime = cell(LP,LM);
    StorageE = cell(LP,LM);
    LossE = cell(LP,LM);
    ThickLoc = cell(LP,LM);
    DephasFFT = cell(LP,LM);
    DephasCos = cell(LP,LM);
    DephasFull = cell(LP,LM);
    ModuleCos = cell(LP,LM);
    R2Sig = cell(LP,LM);
    R2Eps = cell(LP,LM);
    Stiffness = cell(LP,LM);
    
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
        
        Freqs = MF{kc}.FreqsFourier;
        FreqExp{kl,kc} = MF{kc}.ExpFreq;
        FreqWindow = MF{kc}.FreqWindow;
        FreqEch = MF{kc}.FreqEch;
        
        Champ = MF{kc}.field;
        
        EChad = MF{kc}.EChad;
        ExpType = MF{kc}.Exp;
        
        FtDist = MF{kc}.FourierDist;
        FtForce  = MF{kc}.FourierForce;
        FtSigma = MF{kc}.FourierSigma;
        FtEpsilon  = MF{kc}.FourierEpsilon;
       
        
        Sigma = MF{kc}.Sigma;
        Epsilon = MF{kc}.Epsilon;


        F = MF{kc}.force;
        RawDist = MF{kc}.RawSignal;
        BaseDist = MF{kc}.BaseSignal;
        OscDist = RawDist - BaseDist;
        
        T = MF{kc}.time;
        
        
        %% Calcul du déphasage sur zone constante en epaisseur
%         
        % Calcul de Sigma et Epsilon filtré large HF et BF
        
        ptrfreq = find((Freqs<FreqEch-FreqExp{kl,kc}/2&Freqs>FreqEch-FreqExp{kl,kc}*2));
        
        
        ptrfreq = [1 ptrfreq length(Freqs)-ptrfreq+2];
        
        FtSigmaFilt = FtSigma;
        FtSigmaFilt(~ismember(FtSigmaFilt,FtSigmaFilt(ptrfreq))) = 0;
        
        FtEpsilonFilt = FtEpsilon;
        FtEpsilonFilt(~ismember(FtSigmaFilt,FtSigmaFilt(ptrfreq))) = 0;
        
        
        
        SigmaFilt = ifft(FtSigmaFilt);
        EpsilonFilt = ifft(FtEpsilonFilt);
    
        
        % calcul sur une épaisseur qui varie peu
        
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
        
if PLOT && 0
        figure(10)
hold on
plot([-pi pi],[-pi pi],'k-')
plot(DephasFull{kl,kc},Dephas1,'ok','markerfacecolor',col)
plot(DephasFull{kl,kc},Dephas2,'sk','markerfacecolor',col)

        figure(20)
hold on
plot(DephasFull{kl,kc},Err1,'ok','markerfacecolor',col)
plot(DephasFull{kl,kc},Err2,'sk','markerfacecolor',col)
end

        %% Calcul dephasage, module et stiffness en fonction de l'épaisseur
        % avec fourier ou cos fit
        
        
        % calcul sur 3 periode glissante
        
        Texp = 1/FreqExp{kl,kc}; % periode en s
        
        Ncycle = T(end)/Texp;
        
        if Ncycle > 5
            
            FitLength = floor(3*Texp*FreqEch);
            FitShift = floor(Texp/4*FreqEch);
            Nfit = round((Ncycle-3)*4+1);
            
            MeanTime{kl,kc} = zeros(1,Nfit);
            StorageE{kl,kc} = zeros(1,Nfit);
            LossE{kl,kc} = zeros(1,Nfit);
            ThickLoc{kl,kc} = zeros(1,Nfit);
            DephasFFT{kl,kc} = zeros(1,Nfit);
            DephasCos{kl,kc} = zeros(1,Nfit);
            ModuleCos{kl,kc} = zeros(1,Nfit);
            R2Sig{kl,kc} = zeros(1,Nfit);
            R2Eps{kl,kc} = zeros(1,Nfit);
            Stiffness{kl,kc} = zeros(1,Nfit);
            
            
            for k = 1:Nfit
                
                ChampPart = Champ((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                Tpart = T((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                
                SigPartToFT = Sigma((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                EpsPartToFT = Epsilon((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                
                ForPartToFit = F((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                DistPartToFit = OscDist((k-1)*FitShift+1:(k-1)*FitShift+FitLength);

                SigFilPartToFit = SigmaFilt((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                EpsFilPartToFit = EpsilonFilt((k-1)*FitShift+1:(k-1)*FitShift+FitLength);
                
                
                FtSigPart = fft(SigPartToFT);
                FtEpsPart = fft(EpsPartToFT);
                FtChampPart = fft(ChampPart-mean(ChampPart));
                
                FreqsPart = (0:length(FtSigPart)-1)*FreqEch/length(FtSigPart);
                
                
                FreqExpLoc = FreqsPart(find(abs(FtChampPart)>max(abs(FtChampPart)/2),1));
                
                
                FreqOscPtr = abs(FreqsPart - FreqExpLoc) < 1e-5;
                
                Ecmplx = FtSigPart(FreqOscPtr)./FtEpsPart(FreqOscPtr);
                
                MeanTime{kl,kc}(k) = mean(Tpart);
                
                StorageE{kl,kc}(k) = real(Ecmplx);
                LossE{kl,kc}(k) = abs(imag(Ecmplx));
                ThickLoc{kl,kc}(k) = mean(BaseDist((k-1)*FitShift+1:(k-1)*FitShift+FitLength));
                
                DephasFFT{kl,kc}(k) = mod(angle(FtSigPart(FreqOscPtr))-angle(FtEpsPart(FreqOscPtr)),2*pi);
                
                
                
                % figure
                % hold on
                
               
[AmpEpsFilPartFit,PhaseEpsFilPartFit,R2E] = FitCosFunc(EpsFilPartToFit,Tpart,FreqExp{kl,kc});        

                
                % yyaxis right
                % plot(Tpart,EpsFilPartToFit)
                % plot(Tpart,EpsFilPartFit(Tpart))
                % ylabel('Epsilon')
                % xlabel('Time')
                
                
[AmpSigFilPartFit,PhaseSigFilPartFit,R2S] = FitCosFunc(SigFilPartToFit,Tpart,FreqExp{kl,kc});
                
                %                 yyaxis left
                % plot(Tpart,SigFilPartToFit)
                % plot(Tpart,SigFilPartFit(Tpart))
                % ylabel('Sigma')
                %
                %
                % fig = gcf;
                % saveas(fig,['C:\Users\PMMH\Desktop\DataValentin\Figures\MécaDictyFigs\PartialCosFits\' Name '_' ExpType '-' num2str(k) '.png'])
                % close
                
                
                DephasCos{kl,kc}(k) = mod(PhaseSigFilPartFit-PhaseEpsFilPartFit,2*pi);
                ModuleCos{kl,kc}(k) = AmpSigFilPartFit/AmpEpsFilPartFit;
                R2Sig{kl,kc}(k) = R2S;
                R2Eps{kl,kc}(k) = R2E;


                % stiffness
                AmpForPartFit = FitCosFunc(ForPartToFit,Tpart,FreqExp{kl,kc});
                AmpDistPartFit = FitCosFunc(DistPartToFit,Tpart,FreqExp{kl,kc});
                
                Stiffness{kl,kc}(k) = AmpForPartFit/AmpDistPartFit;

                
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
R2Sig = {R2Sig{:}};
R2Eps = {R2Eps{:}};
Stiffness = {Stiffness{:}};



if PLOT

%% Plot fourier

figure(1)
hold on
title('Storage modulus versus thickness (Fourier)')
xlabel('Thickness (nm)')
ylabel('Storage modulus (kPa)')
ylim([0 18])

figure(2)
hold on
title('Loss modulus versus thickness (Fourier)')
xlabel('Thickness (nm)')
ylabel('Loss modulus (kPa)')
ylim([0 7])

figure(3)
hold on
title('Phase Diff versus thickness (Fourier)')
xlabel('Thickness (nm)')
ylabel('Phase Diff (rad)')
% ylim([0 7])


FreqExpR = precisround([FreqExp{:}],0.1);

FreqList = unique(FreqExpR);

nfreq = length(FreqList);

figure(4)
hold on
title('Storage modulus vs Oscillation frequency (Fourier)')
ylabel('E'' (kPa)')
ax = gca;
ax.XTick = 1:nfreq;
xlabel('Freq (Hz)')
xlim([0.5 nfreq+0.5])
ylim([0 18])

figure(5)
hold on
title('Loss modulus vs Oscillation frequency (Fourier)')
ylabel('E'''' (kPa?)')
ax = gca;
ax.XTick = 1:nfreq;
xlabel('Freq (Hz)')
xlim([0.5 nfreq+0.5])
ylim([0 7])

figure(6)
hold on
title('Phase Diff vs Oscillation frequency (Fourier)')
ylabel('Phase Diff (rad)')
ax = gca;
ax.XTick = 1:nfreq;
xlabel('Freq (Hz)')
xlim([0.5 nfreq+0.5])



for k = 1:nfreq
    ptr = FreqExpR == FreqList(k);
    Storageptr = {StorageE{ptr}};
    Lossptr = {LossE{ptr}};
    ThickLocptr = {ThickLoc{ptr}};
    Dephasptr = {DephasFFT{ptr}};



    PlotColorptr = {PlotColor{ptr}};
    
    for kt = 1:length(ThickLocptr)
        
        figure(1)
        hold on
        plot(ThickLocptr{kt},Storageptr{kt}/1000,'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        
        figure(2)
        hold on
        plot(ThickLocptr{kt},Lossptr{kt}/1000,'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        
        figure(3)
        hold on
        plot(ThickLocptr{kt},Dephasptr{kt},'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        
        
    end
    
    figure(4)
    hold on
    Spread_Julien([Storageptr{:}]/1000,k,0.8,'r',symbollist{k},0.4)
    plot([k-0.45 k+0.45],[mean([Storageptr{:}]/1000) mean([Storageptr{:}]/1000)],'--g','linewidth',2)
    ax = gca;
    ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz'];
    
    figure(5)
    hold on
    Spread_Julien([Lossptr{:}]/1000,k,0.8,'b',symbollist{k},0.15)
    plot([k-0.45 k+0.45],[mean([Lossptr{:}]/1000) mean([Lossptr{:}]/1000)],'--g','linewidth',2)
    ax = gca;
    ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz'];
    
    figure(6)
    hold on
    Spread_Julien([Dephasptr{:}],k,0.8,'b',symbollist{k},0.15)
    plot([k-0.45 k+0.45],[mean([Dephasptr{:}]) mean([Dephasptr{:}])],'--g','linewidth',2)
    ax = gca;
    ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz'];

end

%% Plot FitCos

figure(11)
hold on
title('Elastic modulus versus thickness (FitCos)')
xlabel('Thickness (nm)')
ylabel('Elastic modulus (kPa)')
ylim([0 18])

figure(12)
subplot(nfreq,1,1)
hold on
title('Dephas versus time (FitCos)')
subplot(nfreq,1,nfreq)
hold on
xlabel('Time(s)')
ylabel('Phase Diff (rad)')


figure(13)
hold on
title('Elastic modulus vs Oscillation frequency (FitCos)')
ylabel('E'' (kPa)')
ax = gca;
ax.XTick = 1:nfreq;
xlabel('Freq (Hz)')
xlim([0.5 nfreq+0.5])
ylim([0 18])

figure(14)
hold on
title('Phase Diff vs Oscillation frequency (FitCos)')
ylabel('Phase Diff (rad)')
ax = gca;
ax.XTick = 1:nfreq;
xlabel('Freq (Hz)')
xlim([0.5 nfreq+0.5])


figure(15)
hold on
title('Elastic modulus versus time (FitCos)')
xlabel('Average time of sample')
ylabel('Elastic modulus (kPa)')
ylim([0 18])

figure(100)
hold on
title('R2 sig vs Elastic modulus (FitCos)')
ylabel('R2')
xlabel('Elastic modulus (kPa)')
xlim([0 18])

figure(101)
hold on
title('R2 eps vs Elastic modulus (FitCos)')
ylabel('R2')
xlabel('Elastic modulus (kPa)')
xlim([0 18])

figure(1000)
hold on
title('Stiffness vs modulus (pN/nm)')
xlabel('Modulus (kPa)')
ylabel('Stiffness (pN/nm)')

for k = 1:nfreq
    ptr = FreqExpR == FreqList(k);
    Timeptr = {MeanTime{ptr}};
    Moduleptr = {ModuleCos{ptr}};
    ThickLocptr = {ThickLoc{ptr}};
    Dephasptr = {DephasCos{ptr}};
    PlotColorptr = {PlotColor{ptr}};
    R2Sigptr = {R2Sig{ptr}};
    R2Epsptr = {R2Eps{ptr}};
    Stiffnessptr = {Stiffness{ptr}};
    
    for kt = 1:length(ThickLocptr)
        
        figure(11)
        hold on
        plot(ThickLocptr{kt},Moduleptr{kt}/1000,'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        

ptrdeph = Dephasptr{kt}>pi;
Dephasptr{kt}(ptrdeph) = Dephasptr{kt}(ptrdeph)-2*pi;

        figure(12)
        subplot(nfreq,1,k)
        hold on
        plot(Timeptr{kt}(1:4:end),Dephasptr{kt}(1:4:end),'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        ylim([-0.5 pi])
        
        figure(15)
        hold on
        plot(Timeptr{kt},Moduleptr{kt}/1000,'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        
        figure(100)
hold on
 plot(Moduleptr{kt}/1000,R2Sigptr{kt},'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})

        figure(101)
hold on
 plot(Moduleptr{kt}/1000,R2Epsptr{kt},'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        

figure(1000)
hold on
 plot(Moduleptr{kt}/1000,Stiffnessptr{kt},'ok','marker',symbollist{k},'markerfacecolor',PlotColorptr{kt})
        
    end
    
    figure(13)
    hold on
    Spread_Julien([Moduleptr{:}]/1000,k,0.8,'r',symbollist{k},0.4)
    plot([k-0.45 k+0.45],[mean([Moduleptr{:}]/1000) mean([Moduleptr{:}]/1000)],'--g','linewidth',2)
    ax = gca;
    ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz']
    
    figure(14)
    hold on
    Spread_Julien([Dephasptr{:}],k,0.8,'b',symbollist{k},0.15)
    plot([k-0.45 k+0.45],[mean([Dephasptr{:}]) mean([Dephasptr{:}])],'--g','linewidth',2)
    ax = gca;
    ax.XTickLabel{k} = [num2str(FreqList(k)) ' Hz']
end

end





end

function [AmpFit,PhaseFit,R2] = FitCosFunc(Data,Tdata,Freq)
 
        fitcos = fittype('A*cos(2*pi*f*x + p) + M','problem',{'f'},'coefficient', {'A', 'p', 'M'});
        
        AmpStart = (max(Data)-min(Data))/2;
        
        [Fit,gof] = fit(Tdata,Data,fitcos,'problem',{Freq}, ...
            'start',[AmpStart 0 mean(Data)],'upper',[Inf pi max(Data)],'lower',[0 -pi min(Data)]);
        
        Fitc = coeffvalues(Fit);
        
        AmpFit = Fitc(1);

        PhaseFit = Fitc(2);

        R2 = gof.rsquare;

end

function out = precisround(n,precis)

out = round(n./precis).*precis;

end


function H0 = H0ChadCalc(E,F0,R,h)


                % A partir de l'équation de Chadwick :
                % F = (pi*E*R*delta^2)/(3*h0)
                %
                % On écrit delta = h0 - h et on resoud une
                % equation du second degré pour h0
                % a partir des valeurs de h (courbe BF) et EChad,
                % et de la force sans l'oscilation F0
                %
                % h0²*pi*EChad*R - h0*(3*F0 + 2*h*pi*E*R) + pi*EChad*R*h² = 0
                
                
   
             
                
                a = pi*E*R; % 'a' de l'eq second degree
                
                b = -(3.*F0 + 2.*a.*h); % 'b' de l'eq second degree
                
                c = a.*h.^2; % 'c' de l'eq second degree
                
                disc = b.^2 - 4.*a.*c; % discriminant eq second degre en h0
                
                
                
                h0m = (-b - sqrt(disc))./(2*a);
                
                if ~all(h0m<h)
                    error('Probleme avec le calcul de h0')
                end
                
                h0p = (-b + sqrt(disc))./(2*a); % racine + de l'eq, la seule plus grande que h
                
                if ~all(h0p>h)
                    error('Probleme avec le calcul de h0')
                end
                
                H0 = h0p*1e9; % épaisseur non déformé en nm
                
end