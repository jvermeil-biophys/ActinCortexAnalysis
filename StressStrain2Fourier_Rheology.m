function StressStrain2Fourier_Rheology(specif,tags,Diam,PLOT,resfolder,figfolder)
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

MeanTimeFFT = cell(LP,1);
StorageE = cell(LP,1);
LossE = cell(LP,1);
ThickLoc = cell(LP,1);
DephasFFT = cell(LP,1);

FreqExp = cell(LP,1);

PlotColor = cell(LP,1);

CatIdxTmp = 0;

for kl = 1:LP
    
    loadname = list{ptrlist(kl)};
    
    fprintf(['Chargement du fichier ' loadname '\n\n'])
    load(loadname);
    
    
    NameList = {};
    
    LM = length(MF);
    
    MeanTimeFFT{kl} = cell(LM,1);
    StorageE{kl} = cell(LM,1);
    LossE{kl} = cell(LM,1);
    ThickLoc{kl} = cell(LM,1);
    DephasFFT{kl} = cell(LM,1);
    CatIdx{kl} = cell(LM,1);
    
    
    FreqExp{kl} = cell(LM,1);
    
    PlotColor{kl} = cell(LM,1);
    
    for kc = 1:LM
        
        %% récuperation des données
        
        Name = MF{kc}.CellName;
        
        if ~ismember(Name,NameList)
            NameList = {NameList{:} Name};
            PlotColor{kl}{kc} = colorlist{1};
            colorlist = {colorlist{2:end}};
            CatIdxTmp = CatIdxTmp +1;
        else
            NameList = {NameList{:} Name};
            ptr = find(strcmp(NameList,Name),1);
            PlotColor{kl}{kc} = PlotColor{kl}{ptr};
        end
        

        FreqEch = MF{kc}.FreqEch;
        
        Champ = MF{kc}.field;
        
        ExpType = MF{kc}.Exp;
        
        Sigma = MF{kc}.Sigma;
        Epsilon = MF{kc}.Epsilon;
        
        
        Dist = MF{kc}.RawSignal;
        Force = MF{kc}.force;
        
        BaseDist = MF{kc}.BaseSignal;
        BaseForce = MF{kc}.BaseForce;
        
        
        
        T = MF{kc}.time;
        
        
        %         if FreqExp{kl}{kc} > 0.5 && FreqExp{kl}{kc} < 2 % travail temporaire sur 1Hz uniquement
        
        %% Calcul module complexe et déphasage sur X periodes avec Fourier
        
        % calcul sur X periode glissante (Y décalage)
        
        
        Texp = 1/MF{kc}.ExpFreq;; % periode en s
        
        Ncycle = T(end)/Texp;
        
        X = 3; % nb de periodes a fitter
        
        Y = X; % nb de periode de décalage de la fenetre entre chaque fit
        
        
        
        if Ncycle > 5 && X <= Ncycle
            
            FitLength = floor(X*Texp*FreqEch); % X periodes
            FitShift = floor(Texp*FreqEch*Y); % shift entre les fit, données indépendantes si >= FitLength
            Nfit = max([floor(Ncycle/Y) + Y - X, 1 ]); % nombre de fit
            
            MeanTimeFFT{kl}{kc} = zeros(1,Nfit);
            
            ThickLoc{kl}{kc} = zeros(1,Nfit);
            
            StorageE{kl}{kc} = zeros(1,Nfit);
            LossE{kl}{kc} = zeros(1,Nfit);
            DephasFFT{kl}{kc} = zeros(1,Nfit);
            CatIdx{kl}{kc} = zeros(1,Nfit);
            
                    FreqExp{kl}{kc} = zeros(1,Nfit);
            
            if PLOT || 0
                figure
                hold on
                ylabel('Module élastique (kPa)')
                xlabel('Nombre d''itération')
                title(['Module élastique avec Fourier pour ' Name '-' ExpType ' sur ' num2str(X) ' périodes ' num2str(Nfit) ' fois (' num2str(Y) ' période(s) de décalage).'])
            end
            
            for k = 1:Nfit
                
                ptrFit = (k-1)*FitShift+1:(k-1)*FitShift+FitLength; % Portion de courbe a fitter
                % Temps moyen sur la zone analysée
                Tpart = T(ptrFit);
                
                % Boucle sur E pour convergence
                
                nloop = 5;
                Etmp(1) =  MF{kc}.EChad;
                
                Start = 1;
                
                while Start == 1 || abs((Etmp(end)-Etmp(end-1))/Etmp(end-1))>0.01
                    
                    
                    DistFit = Dist(ptrFit);
                    ForceFit = Force(ptrFit);
                    
                    BaseDistFit = BaseDist(ptrFit);
                    BaseForceFit = BaseForce(ptrFit);
                    %
                    %
                    [SigPartToFT,EpsPartToFT] = Compute_StressStrain_Rheo(DistFit,ForceFit,Etmp(end),BaseForceFit*1e-12,Diam/2*1e-9,BaseDistFit*1e-9);
                    %
                    
                    % sans la loop
                    %                     % Sigma et epsilon + fft sur la zone
                    %                     SigPartToFT = Sigma(ptrFit);
                    %                     EpsPartToFT = Epsilon(ptrFit);
                    
                    
                    FtSigPart = fft(SigPartToFT);
                    FtEpsPart = fft(EpsPartToFT);
                    
                    FreqsPart = (0:length(FtSigPart)-1)*FreqEch/length(FtSigPart);
                    
                    
                    % Frequence exp sur la portion de courbe a partir du signal
                    % du champ + ptr de cette freq
                    
                    ChampPart = Champ(ptrFit);
                    FtChampPart = fft(ChampPart-mean(ChampPart));
                    FreqExpLoc = FreqsPart(find(abs(FtChampPart)>max(abs(FtChampPart)/2),1));
                    
                    FreqOscPtr = abs(FreqsPart - FreqExpLoc) < 1e-5;
                    
                    
                    % Module complexe a la fréquence selectionée
                    Ecmplx = FtSigPart(FreqOscPtr)./FtEpsPart(FreqOscPtr);
                    
                    Etmp(end+1) = real(Ecmplx);
                    
                    Start = 0;
                    
                    if ~(Etmp(end) > 0)
                        break
                    end
                end
                
                
                
                
                
                if ~(Etmp(end) > 0)
                    
                MeanTimeFFT{kl}{kc}(k) = NaN;
                
                    StorageE{kl}{kc}(k) = NaN;
                    LossE{kl}{kc}(k) = NaN;
                    
                    ThickLoc{kl}{kc}(k) = NaN;
                    
                    
                    DephasFFT{kl}{kc}(k) = NaN;
                    
                    CatIdx{kl}{kc}(k) = NaN;
                    
                else
                    if PLOT || 0
                        plot(Etmp/1000,'-o','linewidth',1.5)
                    end
                    
                MeanTimeFFT{kl}{kc}(k) = mean(Tpart);
                
                    StorageE{kl}{kc}(k) = real(Ecmplx);
                    LossE{kl}{kc}(k) = abs(imag(Ecmplx));
                    
                    ThickLoc{kl}{kc}(k) = mean(BaseDist(ptrFit));
                    
                            FreqExp{kl}{kc}(k) = MF{kc}.ExpFreq;
                    DephasFFT{kl}{kc}(k) = mod(angle(FtSigPart(FreqOscPtr))-angle(FtEpsPart(FreqOscPtr))+1,2*pi)-1;
                    
                    CatIdx{kl}{kc}(k) = CatIdxTmp;
                    
                end
                
                clear Etmp SigPartToFT EpsPartToFT BaseForceFit BaseDistFit ForceFit DistFit FtSigPart FtEpsPart Ecmplx
                
                
            end
        end
        
        
        
    end
    
    %     end
    
    
end

ThickLocPlot = [ThickLoc{:}{:}];
StorageEPlot = [StorageE{:}{:}];
LossEPlot = [LossE{:}{:}];
DephasFFTPlot = [DephasFFT{:}{:}];
CatIdxPlot = [CatIdx{:}{:}];
PlotColorPlot = [PlotColor{:}{:}];
FreqExpR = precisround([FreqExp{:}{:}],0.1);

ptrnan = isnan(CatIdxPlot);

ThickLocPlot(ptrnan) = [];
StorageEPlot(ptrnan) = [];
LossEPlot(ptrnan) = [];
DephasFFTPlot(ptrnan) = [];
CatIdxPlot(ptrnan) = [];
FreqExpR(ptrnan) = [];

StorageEperCell = StorageE{:};
LossEperCell = LossE{:};
MeanTimePerCell = MeanTimeFFT{:};
FreqExpPerCell = FreqExp{:};
CatIdxPerCell = cellfun(@nanmean,CatIdx{:});

FreqList = unique(FreqExpR);

FreqList = FreqList(FreqList~=0);

nfreq = length(FreqList);


%% Plot fourier

if PLOT || 1
    
    %% Modules v Epaisseur
%     figure
%     hold on
%     title('Storage modulus versus thickness (Fourier)')
%     xlabel('Thickness (nm)')
%     ylabel('Storage modulus (kPa)')
%     plot(ThickLocPlot,StorageEPlot/1000,'o')
%     
%     
%     figure
%     hold on
%     title('Loss modulus versus thickness (Fourier)')
%     xlabel('Thickness (nm)')
%     ylabel('Loss modulus (kPa)')
%     plot(ThickLocPlot,LossEPlot/1000,'o')
    
    %% Module elastique v Freq
    figure
    set(0,'DefaultLineMarkerSize',6);
    hold on
    title('Storage modulus as a function of frequency')
    for kfreq = 1:nfreq
        
        ptrfreq = FreqExpR == FreqList(kfreq);
        
        data{kfreq} = StorageEPlot(ptrfreq)/1000;
        datamarkers{kfreq} = 'o';
        datanames{kfreq} = [num2str(FreqList(kfreq)) ' Hz'];
        
                plot([kfreq-0.45 kfreq+0.45],[mean(data{kfreq}) mean(data{kfreq})],'r--','linewidth',2)
                
    end

plotSpread_V(StorageEPlot/1000,'distributionMarkers',datamarkers,'distributionIdx',FreqExpR,'xValues',1:nfreq,'xNames',datanames,'spreadWidth',1)

        %'categoryIdx',...    CatIdxPlot,
        
    xlabel('Frequency')
    ylabel('Storage modulus (kPa) - LogScale')

    ax = gca;
    ax.YScale = 'log';
    
    xlim([0.5 nfreq+0.5])
    
    plotheight = 55;
    
    for ii = 1:length(data)-1
        for jj = ii+1:length(data)
            
            sigtxt = DoParamStats(data{ii},data{jj});
            
            plot([ii jj],[plotheight plotheight],'k-','linewidth',0.5)
            text((ii+jj)/2, plotheight+5, sigtxt,'HorizontalAlignment','center','fontsize',14)
            
            plotheight = plotheight+15;
        end
    end
    
    %% Module visqueux v Freq
    
     figure
    hold on
    title('Loss modulus as a function of frequency')
    for kfreq = 1:nfreq
        
        ptrfreq = FreqExpR == FreqList(kfreq);
        
        data{kfreq} = LossEPlot(ptrfreq)/1000;
        datamarkers{kfreq} = 'o';
        datanames{kfreq} = [num2str(FreqList(kfreq)) ' Hz'];
        
        
                plot([kfreq-0.45 kfreq+0.45],[mean(data{kfreq}) mean(data{kfreq})],'r--','linewidth',2)
                
    end

    xlabel('Frequency')
    ylabel('Storage modulus (kPa) - LogScale')
    xlim([0.5 nfreq+0.5])

plotSpread_V(LossEPlot/1000,'distributionMarkers',datamarkers,'distributionIdx',FreqExpR,'categoryIdx',...
    CatIdxPlot,'xValues',1:nfreq,'xNames',datanames,'spreadWidth',1)


        ax = gca;
    ax.YScale = 'log';
    
        plotheight = 16;
    
    for ii = 1:length(data)-1
        for jj = ii+1:length(data)
            
            sigtxt = DoParamStats(data{ii},data{jj});
            
            plot([ii jj],[plotheight plotheight],'k-','linewidth',0.5)
            text((ii+jj)/2, plotheight+2, sigtxt,'HorizontalAlignment','center','fontsize',14)
            
            plotheight = plotheight+6;
            
        end
    end
    
    %% Module Elastique par cellule en temps
    
    figure
    hold on
       xlabel('Time')
    ylabel('Storage modulus (kPa)')
    
    title('Elastic modulus in time, one color per cell')
    

    CellIdxTmp = 0;
    
    ncells = length(StorageEperCell);
    
    Colors = distinguishable_colors(length(unique(CatIdxPerCell)));
    
    for kcell = 1:ncells
 
        if CellIdxTmp == CatIdxPerCell(kcell)
           
            ptrnonan = find(~isnan(PlotTimeTmp));
            lastTime = PlotTimeTmp(ptrnonan(end));
            PlotTimeTmp = MeanTimePerCell{kcell} + lastTime + 10;
            
        else
            CellIdxTmp = CatIdxPerCell(kcell);
            Color = Colors(CellIdxTmp,:);
            PlotTimeTmp = MeanTimePerCell{kcell} +2*kcell;

        end
        

        plot(PlotTimeTmp,StorageEperCell{kcell}/1000,'ok','markerfacecolor',Color)
        
       
    end
    
end

end



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

function out = precisround(n,precis)

out = round(n./precis).*precis;

end