function [nCompTot, nCompOk] = Data2Meca_J(date,manip,tag,specif,ExperimentalConditions,PLOT,datafolder,resfolder,figurefolder,EXCLUDED)

%% INIT

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked') % docked normal
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',25);

here = pwd;

R2CRITERION = 0.8;

PLOTTMP = 0;

% constantes
DIAMETER = ExperimentalConditions.DiametreBilles; % Diametre de la bille en nm
RADIUS = DIAMETER/2; % Rayon de la bille en nm
SCALE=ExperimentalConditions.Scale;
MAGFIELDCORRECTIONCOEFF = ExperimentalConditions.MagFieldCorrectionCoeff;
RAMPFIELD = ExperimentalConditions.RampField;
FITPERCENT = 0.2;


%% params and settings

datenow = datestr(now,'yy-mm-dd');

% dossier de données et sauvegarde des données  
% path = 'D:\Data\MatFiles\V2D';
% fieldpath =['D:\Data\Raw\' date(7:8) '.' date(4:5) '.' date(1:2)];
% sf   = ['D:\Data\MatFiles' filesep 'D2M'];
fieldpath = [datafolder filesep date(7:8) '.' date(4:5 ) '.' date(1:2)];
path = [resfolder filesep 'V2D'];
sf = [resfolder filesep 'D2M'];


mkdir(sf)

% sauvergarde des figures
subFigureFolder = [figurefolder filesep datenow filesep 'Meca'];
mkdir(subFigureFolder)

if PLOT
    sffDistCurve = [subFigureFolder filesep specif filesep 'DistCurve'];
    mkdir(sffDistCurve)
    
    sffMeca = [subFigureFolder filesep specif filesep 'ForceDist&StressStrain'];
    mkdir(sffMeca)
end



H0abs = {}; % absolu
H0prc = {}; % en percentile de la courbe D3cst
E0 = {}; % elasticité sur le debut des deformation
AlphaNL = {}; % exposant non lineaire
CompNum = {}; % numero de la compression
Hys = {}; % hysteresis
Cort = [];



%% Retrieve datafile
loadname=['V2D_' date '_' manip '_' tag '_' specif '.mat'];

cd(path)

if exist(loadname,'file')
    fprintf(['Chargement du fichier ' loadname '\n\n'])
    
    load(loadname);
    nCells = length(MR); % nombre de cellules
%     cptfit = 0;
    nCompOk = 0;
    nCompTot = 0;
    firstCompressionIndex = 0;
    
    %% Get datas
%     for i=2
    for i=1:nCells
        
        if not(isempty(MR{i}))
            H0abs{i} = []; % absolu
            H0prc{i} = []; % en percentile de la courbe D3cst
            E0{i} = []; % elasticité sur le debut des deformation
            Edim{i} = []; % fit Dimitriadis 
            Echad{i} = []; % fit Chadwick
            Echad_Median{i} = 0;
            Echad_WeightedAvg{i} = 0;
            Echad_Variance{i} = 0;
            EchadConf{i} = []; % Confidence interval width around Echad value
            Delt3{i} = []; % Enfoncement a 5mT
            Delt20{i} = []; % Enfoncement theorique a 20 pN
            ChadsH0{i} = [];
            AlphaNL{i} = []; % exposant non lineaire
            CompNum{i} = [];
            Hys{i} = [];
            

            
            % récupération des données
            currentCellData = MR{i};
            currentCellName{i} = currentCellData.name;
            
            isAllowed = 1;
            for k=1:length(EXCLUDED)
                TEST = contains(currentCellName{i}, EXCLUDED{k});
                if contains(currentCellName{i}, EXCLUDED{k})
                    isAllowed = 0;
                end
            end
            
            indexC = find(currentCellName{i} == 'C');
            index_ = find(currentCellName{i} == '_');
            
            index_ = index_(index_>indexC);
            
            CellN(i) = str2num(currentCellName{i}(indexC+1:index_(1)-1));
            
            if strcmp(ExperimentalConditions.SubstrateCoating, 'Float')
            	distanceIndex = 7;
                D3 = currentCellData.D2*1000 - DIAMETER; % En nm
            else
            	distanceIndex = 3;
                D3 = currentCellData.D3*1000 - DIAMETER; % En nm
            end
            
            
            dz = currentCellData.dz; % En nm
            T = currentCellData.time;
            F = currentCellData.F;
            
            % Donnée force cst [cst = constante]
            CstData = currentCellData.CstData;
            D3cst = CstData(:,distanceIndex)*1000 - DIAMETER; % En nm
            Tcst = CstData(:,2);
            D3smooth = smooth(Tcst,D3cst,7,'sgolay',5);
            
            % calcul du cortex et des centile
            medianCorticalThicknessCst = median(D3cst);
            halfPercentsGrid = [0:0.5:100];
            percentilesCorticalThicknessCst = prctile(D3cst,halfPercentsGrid);
            
            medianCortThick(i) = medianCorticalThicknessCst;
            
            % Données rampe
            FullRampData = currentCellData.RampData{2}; % all data
            D3rmp = FullRampData(:,distanceIndex)*1000 - DIAMETER; % En nm
            Trmp = FullRampData(:,2);
            
            RampData = currentCellData.RampData{1}; % good data
            nComp = length(RampData);
            nCompTot = nCompTot + nComp;
            
            %%%%%FIGURE%%%%%%
            if PLOT
                figure
                hold on
                
                yyaxis left
                axl = gca;
                axl.YColor = [0 0 0];
                ylabel('Cortex thickness (nm)','FontName','Serif')
                plot(T,D3,'k-','linewidth',1)
                plot(Tcst,D3cst,'ok','markerfacecolor','b')
                plot(Trmp,D3rmp,'o','markerfacecolor','r','markersize',4,'markeredgecolor','none')
                
                yyaxis right
                axr = gca; % Get Current Axes
                axr.YColor = [0 0 0];
                plot(T,F,'r-')
                ylabel('Force (pN)','FontName','Serif')
                ylim([0 2500])
                
                title(currentCellName{i})
                
                xlabel('Time(s)','FontName','Serif')
                
                fig = gcf; % Get Current Figure
                saveas(fig,[sffDistCurve filesep currentCellName{i}],'png')
                saveas(fig,[sffDistCurve filesep currentCellName{i}],'fig')
%                 close
            end
            
            if contains(currentCellName{i}, '07-08-20_M1_P1_C7')
                figure
                hold on
                
                yyaxis left
                axl = gca;
                axl.YColor = [0 0 0];
                ylabel('Cortex thickness (nm)','FontName','Serif')
                plot(T,D3,'k-','linewidth',1.5)
                plot(Trmp,D3rmp,'o','markerfacecolor','b','markersize',4,'markeredgecolor','none')
                ylim([0 350])
                
                yyaxis right
                axr = gca; % Get Current Axes
                axr.YColor = [0 0 0];
                plot(T,F,'r-','linewidth',1.5)
                ylabel('Force (pN)','FontName','Serif')
                ylim([0 2700])
                
                xlim([0 60])
                
                title("Compression assay, aSFL 3T3", "fontsize", 30,'FontName','Serif')
                
                xlabel('Time(s)','FontName','Serif')
                
                fig = gcf; % Get Current Figure
                saveas(fig,['C:\Users\JosephVermeil\Desktop\Figures' filesep currentCellName{i}],'png')
                saveas(fig,['C:\Users\JosephVermeil\Desktop\Figures' filesep currentCellName{i}],'fig')
%                 close
            end
            
            %% Mecanique
            
            for ii = 1:nComp
                if not(ischar(RampData{ii}))

                    Btot = RampData{ii}(:,5); % Champ mag sur toute la duréee de la rampe
                    
                    % Les lignes qui suivent servent à detecter le début de la rampe de champ magnetique. 
                    % Pour Bamplitude = Bmax - Bmin, on place le début de la rampe à l'index pour lequel 
                    % B = Bmin + Bamplitude/40 et la fin à l'index où B = Bmax.
                    [maxBtot,indexMaxBtot] = max(Btot);
                    minBtot = min(Btot);
                    thresholdBtot = (maxBtot - minBtot) / 40;
                    [~,indexThresholdBtot] = min(abs(Btot(1:indexMaxBtot) - (minBtot+thresholdBtot)));
                    
                    
                    % l'index bis correspond à la fin de la partie descendante de la rampe.
                    [~,indexThresholdBtotBis] = min(abs(Btot(indexMaxBtot:end) - (minBtot+thresholdBtot)));
                    indexThresholdBtotBis = indexThresholdBtotBis + indexMaxBtot - 1;
                    
                    
%                     % TEST PLOT
%                     figure
%                     pseudoTime = [0:1:length(Btot)-1];
%                     subplot(2, 1, 1)
%                     hold on
%                     plot(pseudoTime, Btot, 'b--')
%                     plot(pseudoTime(indexThresholdBtot:indexMaxBtot), Btot(indexThresholdBtot:indexMaxBtot), 'r-')
%                     plot(pseudoTime(indexMaxBtot:indexThresholdBtotBis), Btot(indexMaxBtot:indexThresholdBtotBis), 'g-')
%                     
%                     subplot(2, 1, 2)
%                     hold on
%                     plot(pseudoTime, RampData{ii}(:,3), 'b--')
%                     plot(pseudoTime(indexThresholdBtot:indexMaxBtot), RampData{ii}(indexThresholdBtot:indexMaxBtot,3), 'r-')
%                     plot(pseudoTime(indexMaxBtot:indexThresholdBtotBis), RampData{ii}(indexMaxBtot:indexThresholdBtotBis,3), 'g-')
                    
                    pointeurBaseField = Btot(1:indexMaxBtot) < MAGFIELDCORRECTIONCOEFF*(RAMPFIELD(1) + 0.05*RAMPFIELD(1)) & Btot(1:indexMaxBtot) > MAGFIELDCORRECTIONCOEFF*(RAMPFIELD(1) - 0.05*RAMPFIELD(1)); % Intersection d'ensembles

                    % Definition of the 2 parts of the ramps
                    Dup = RampData{ii}(indexThresholdBtot:indexMaxBtot,distanceIndex) - DIAMETER/1000; % Distance (en µm) pendant la phase de compression
                    Ddown = RampData{ii}(indexMaxBtot:indexThresholdBtotBis,distanceIndex) - DIAMETER/1000; % Distance (en µm) pendant la phase de décompression

%                     % réhaussement des compressions qui passent en dessous de 30 nm
%                     lowD = min(Dup);
%                     if lowD < 0.03 
%                         Dup = Dup - lowD + 0.03; 
%                         Ddown = Ddown - lowD + 0.03;
%                     end

                    Fup = RampData{ii}(indexThresholdBtot:indexMaxBtot,4); % Force (en pN) pendant la phase de compression
                    Fdown = RampData{ii}(indexMaxBtot:indexThresholdBtotBis,4); % Force (en pN) pendant la phase de décompression

                    %% H0 pour plot

                    InitialThick = mean(Dup(pointeurBaseField)) * 1000; % absolute thickness at the begining of compression

                    %                     [~,ind] = min(abs(percentilesCorticalThicknessCst-AbsoluteInitialThick));
                    %                     percentileOfH0 = halfPercentsGrid(ind); % percentile of data

                    percentileOfH0 = ((InitialThick - min(D3cst)) / (max(D3cst) - min(D3cst))) * 100; % percentage of amplitude


                    %% Analyse mecanique

                    % definition des variables et constantes
                    Hcomp = Dup; % épaisseur du cortex pendant la compression en µm
                    Fcomp = Fup; % force pendant la compression

                    differentialHcomp = diff(Hcomp(5:end));
                    MaxDisplacement = max(differentialHcomp)*1000; % en nm

                    pointeurMonotonie =  differentialHcomp > 0; % 1 quand croissant 0 quand decroissant

                    MaxDisplacementCummulative = MaxDisplacement; % deplacement maximum cumulé sur plus de trois points de suite

                    GrowthPhases = bwconncomp(pointeurMonotonie); % get connected component where pointeurMonotonie == 1 <=> zones de croissance
                    % GrowthPhases.PixelIdxList{k} -> La liste des points de la k ième zone de croissance.

                    for k = 1:GrowthPhases.NumObjects % nombre de zones
                        if length(GrowthPhases.PixelIdxList{k}) > 2 % si plus de trois points croissants d'affilée.
                            MaxDisplacementCummulative =  max([MaxDisplacementCummulative sum(differentialHcomp(GrowthPhases.PixelIdxList{k}))*1000]);
                            % on calcule de combien ca monte sur ces points, si c'est plus grand on note la valeur.
                        end
                    end

                    % On considere que la courbe d'indentation est normalement décroissante 
                    % si il n'y a pas de deplacement de plus de 50nm en 1 pts, ou de remontée 
                    % de plus de 50nm sur 3 pts ou plus.

                    % QUESTION : Quels sont les cas qu'on vire exactement
                    % en faisant ce genre de trucs ???

                    isDecreasing = (MaxDisplacementCummulative < 50) && (MaxDisplacement < 50);



                    %                         figure(1)
                    %                         hold on
                    %                         plot(Hcomp*1000,'g*-')
                    %                         a=1;

                    

                    % smoothing
                    FHsort = sortrows([Hcomp Fcomp]); % B = sortrows(A) sorts the rows of a matrix in ascending order based on the elements in the first column.
                    Fcomp = FHsort(:,2);
                    Hcomp = FHsort(:,1);
                    Fcomp = smooth(Hcomp,Fcomp,7,'sgolay',5);
                    % Goal of 2 lines below : replace the potential negative values with a 0. 
                    % Because otherwise the sqrt() a little below will create ugly complex numbers.
                    Fcomp_isPositive = [Fcomp > 0];
                    Fcomp = Fcomp .* Fcomp_isPositive;
                    %
                    Hcomp = Hcomp (end:-1:1);
                    Fcomp = Fcomp (end:-1:1);


                    % Recuperation des 150 premiers pN d'enfoncement pour estimer H0
                    pointeurFirst150pN = (Fcomp<min(Fcomp)+150);

                    F_H0 = Fcomp(pointeurFirst150pN);
                    H_H0 = Hcomp(pointeurFirst150pN);

                    % estimation de H0 par fit elastique lineaire
                    F_H0el = F_H0.^(1/2);
                    LinFitCoeff_H0 = polyfit(F_H0el, H_H0, 1);
                    H0el = LinFitCoeff_H0(2); % épaisseur de départ

                    isCompressed = all(H0el > Hcomp);


                    %% fit dimitriadis
                    %                              fitdim = fittype(['Ef*16/9*sqrt(R)*(max(x0-x,0).^1.5+0.884*rh*max(x0-x,0).^2+0.781*' ...
                    %                                 'rh^2*max(x0-x,0).^2.5+0.386*rh^3*max(x0-x,0).^3+0.0048*rh^4*max(x0-x,0).^3.5)']...
                    %                                 ,'problem',{'R','rh'}); % fonction a fitter (dimitriadis), R et rh paramètre, Ef et x0 a fitter
                    %                             
                    %                             
                    %                            
                    %                               % on fit de 20 a 500pN
                    %                             
                    %                                FcompFit = Fcomp(Fcomp>20&Fcomp<500);
                    %                                HcompFit = Hcomp(Fcomp>20&Fcomp<500); % en Pa
                    %                            
                    %                             
                    %                             FFF =fit(HcompFit,FcompFit,fitdim,'problem',{Rb, sqrt(Rb)/H0el},'Start',[1000,H0el],...
                    %                                 'Lower',[1 0.8*mean(H_H0)],'Upper',[Inf 1.2*mean(H_H0)]);
                    %                             
                    %                             FFFc=coeffvalues(FFF);
                    %                             
                    %                             EdimFit=FFFc(1); % module
                    %                             H0fitFit=FFFc(2); % épaisseur de départ
                    %                             
                    %                             Fcompfromfit = FFF(HcompFit);
                    %                             
                    %                             R2dim = 1 - sum((FcompFit - Fcompfromfit).^2)/sum((FcompFit - mean(FcompFit)).^2);
                    %                             
                    % 
                    %                                     Edimtmp = EdimFit;
                    %                                     R2dimtmp = R2dim;
                    %                                     
                    %                                     HcompFromFit = HcompFit*1000;
                    %                                     FcompFromFit = FFF(HcompFit);
                    %                                     
                    % 
                    %                             
                    %                             R2dim = R2dimtmp;

                    %% fit chadwick pour delta/2 et h/2

                    fitchad = fittype(['(pi*E*R*(x-h0)^2)/(3*h0)'],'problem',{'R'});

                    % On fitte tant que delta/h < 0.2 [Critère de validité du modèle ?????]

                    InitialThick_um = InitialThick/1000;
                    H3mT = mean(Dup(pointeurBaseField)); % 3 mT car correspond au plat au début de la rampe.
                    pointeurInferieurA20Pourcents = ((InitialThick_um-Hcomp)/InitialThick_um) < FITPERCENT;
                    FcompFit = Fcomp(pointeurInferieurA20Pourcents);
                    HcompFit = Hcomp(pointeurInferieurA20Pourcents); % en Pa

                    try
                        FFF = fit(HcompFit,FcompFit,fitchad,'problem',{RADIUS/1000},'Start',[5000 H0el],... % [1000 H0el] [H0el 1000]
                            'Lower',[-Inf max(HcompFit)],'Upper',[Inf Inf]);
                        % Starting values : 1000 for FcompFit, H0el for H0.
                    catch
                        isCompressed = 0;
                    end

                    FFFcoeff = coeffvalues(FFF);

                    EchadFit = FFFcoeff(1); % module
                    H0chadFit = FFFcoeff(2); % FFFcoeff(2) est une nouvelle, meilleure estimation de H0.

                    FFFci = confint(FFF);
                    EchadFitCI = abs(diff(FFFci(:,1)));
%                     EchadFitCI = FFFci(:,1);
%                     H0chadFitCI = FFFci(:,2);
                    
                    

                    Delta3mT = (H0chadFit - H3mT)*1000; % enfoncement a 3mT en nm

                    Fcompfromfitchad = FFF(HcompFit);

                    StressChad = Fcomp ./ (pi*RADIUS*(H0chadFit - Hcomp));
                    StrainChad = (H0chadFit - Hcomp) / (3*H0chadFit);

                    DistTh = min(HcompFit):0.001:H0chadFit;
                    FcompTh = FFF(DistTh);

                    pointeur20pN = FcompTh > 19 & FcompTh < 21;
                    Delta20pN = (H0chadFit - mean(DistTh(pointeur20pN)))*1000;

                    R2chad = 1 - sum((FcompFit - Fcompfromfitchad).^2)/sum((FcompFit - mean(FcompFit)).^2);

                    Echadtmp = EchadFit;

                    isAGoodFit = (R2chad > R2CRITERION); %&& (Echadtmp<30000);
                    isPositive = (Echadtmp>0);
                    %                             if isnan(Delta5mT)
                    %                                 figure
                    %                                 hold on
                    %                                 plot(HcompFit*1000,FcompFit,'-*b','linewidth',1.1)
                    %                                 plot(HcompFit*1000,Fcompfromfitchad,'-og','linewidth',1.5)
                    %                                 close
                    %                             end                            

                     %% fit chadwick-like
                    % calcul de l'indentation maximum delta
                    Delt = H0el - Hcomp; % indentation en µm

                    % calcul d'epsilon et sigma
                    EpsEq = Delt/(2*H0el); % delta equivalent = delta / 2 -> volume equivalent compresser par un cylindre
                    Sig = Fcomp./(pi*RADIUS*Delt); % stress en Pa

                     %% calcul de l'elasticité a petite deformation
                    %                             
                    %                             E0tmp = [];
                    %                             
                    %                             % on fite sur 100Pa, puis 150, puis 200, etc...
                    %                             
                    %                             step = 50;
                    %                             R2tmp = 0;
                    %                             old_pos = 1;
                    %                             
                    %                             minSig = min(Sig);
                    %                             
                    %                             E0tmp = 0;
                    %                             
                    %                             while max(Sig)> minSig+step+50
                    %                                 
                    %                                 step = step + 50;
                    %                                 
                    %                                 [~,pos] = min(abs(Sig(old_pos:end)-(minSig+step)));
                    %                                 
                    %                                 pos = max(pos+old_pos-1,5);
                    %                                 
                    %                                 EpsFit = EpsEq(1:pos);
                    %                                 SigFit = Sig(1:pos); % en Pa
                    %                                 
                    %                                 old_pos = pos;
                    %                                 
                    %                                 
                    %                                 [Rfit,stats] = robustfit(EpsFit,SigFit);
                    %                                 
                    %                                 E0fit = Rfit(2); % en Pa
                    %                                 
                    %                                 Sigfromfit = Rfit(1) + Rfit(2)*EpsFit;
                    %                                 
                    %                                 R2 = 1 - sum(stats.w.*(SigFit - Sigfromfit).^2)/sum(stats.w.*(SigFit - mean(SigFit)).^2);
                    %                                 
                    %                                 
                    %                                 
                    %                                 if R2 > R2sig && R2 > R2tmp && E0fit > 0
                    %                                     E0tmp = E0fit;
                    %                                     R2tmp = R2;
                    %                                     
                    %                                     EpsFromFit = EpsFit*100;
                    %                                     SigFromFit = Rfit(1) + Rfit(2)*EpsFit;
                    %                                     
                    %                                     
                    %                                 end
                    %                                 
                    %                             end
                    %                             
                    %                             R2 = R2tmp;

                    %% calcul de l'exposant non lineaire
                    Lesp = log(EpsEq);
                    Lsig = log(Sig);
                    Pnl = polyfit(Lesp,Lsig,1);
                    Lsigff = Pnl(2) + Pnl(1)*Lesp;
                    R2al = 1 - sum((Lsig - Lsigff).^2)/sum((Lsig - mean(Lsig)).^2);

                    
                    
                    %% Sauvegarde des valeurs pertinentes
                    dDup = abs(diff(Dup*1000));
                    hasNoLargeJump = (max(dDup) < 1000);
%                     
%                     minD = min(Dup) * 1000;
%                     isAboveZero = (minD > 20);
                    
                    hasSmallConfInt = (abs(EchadFitCI) < 2*abs(Echadtmp));

                    if ~strcmp(ExperimentalConditions.SubstrateCoating, 'Float')
                        currentRampDz = RampData{ii}(:,6);
                        ddz = abs(diff(RampData{ii}(:,6)));
                        hasNormalDz = (max(currentRampDz) < 4) && (max(ddz) < 0.7);
                    else
                        hasNormalDz = 1;
                    end
                    
                    

                    SAVING_CRITERION = isAllowed && isPositive && isCompressed && hasNoLargeJump && hasSmallConfInt && hasNormalDz;
                    
                        % on prend les valeurs pour plot
                    if SAVING_CRITERION
                        H0abs{i} = [H0abs{i}; InitialThick]; % valeur de l'épaisseur du cortex avant compression
                        H0prc{i} = [H0prc{i}; percentileOfH0]; % percentile correspondant
%                        E0{i} = [E0{i}; E0tmp];
%                        Edim{i} = [Edim{i}; Edimtmp];
                        Echad{i} = [Echad{i}; Echadtmp];
                        EchadConf{i} = [EchadConf{i}; EchadFitCI];
                        Delt3{i} = [Delt3{i}; Delta3mT];
                        Delt20{i} = [Delt20{i}; Delta20pN];
                        ChadsH0{i} = [ChadsH0{i}; H0chadFit*1000];
                        CompNum{i} = [CompNum{i}; ii];
                        AlphaNL{i} = [AlphaNL{i}; Pnl(1)]; % Semble inutile ?
                    

                        %% calcul de l'hysteresis

                        AreaComp = abs(trapz(Dup,Fup));
                        AreaRel = abs(trapz(Ddown,Fdown));
                        Hys{i} = [Hys{i}; (AreaComp-AreaRel)/AreaComp*100];
                        nCompOk = nCompOk +1;
                    end
                    
                    %% %%%%%%% FIGURES %%%%%%
                    if PLOT
                        figure
                        if SAVING_CRITERION
                            sgtitle(['Chadwick model for ' currentCellName{i} ' - ' num2str(ii)])
                        else
                            sgtitle(['NOT SAVED DATA - Chadwick model for ' currentCellName{i} ' - ' num2str(ii)])
                        end
                        subplot(121)
                        hold on
                        title(['Force-Thickness, Echad = ' num2str(Echadtmp) 'kPa'])
                        plot(Dup*1000,Fup,'-*b','linewidth',1.1)
                        plot(Ddown*1000,Fdown,'-*r','linewidth',1.1)
                        plot(HcompFit*1000,Fcompfromfitchad,'-m','linewidth',3)
                        coeffPlotHF = ((3*H0chadFit)/(pi*EchadFit*(RADIUS/1000)))^(1/2);
                        plot((H0chadFit - coeffPlotHF*[0:5:Fcomp(1) Fcomp'].^(1/2))*1000,[0:5:Fcomp(1) Fcomp'],'xm','linewidth',1.5,'handlevisibility','off')
                    %                                 plot((LinFitCoeff_H0(2)+LinFitCoeff_H0(1)*[F_H0'].^(1/2))*1000,[F_H0'],'-m','linewidth',3)
                    %                                 plot((LinFitCoeff_H0(2)+LinFitCoeff_H0(1)*[0:5:F_H0(1) F_H0'].^(1/2))*1000,[0:5:F_H0(1) F_H0'],'--m','linewidth',2,'handlevisibility','off')
                    %                                 plot(LinFitCoeff_H0(2)*1000,0,'*m','linewidth',5)
                    %                                 plot(Dup*1000,Fup,'-*b','linewidth',1.1,'handlevisibility','off')
                    %                                 plot(Ddown*1000,Fdown,'-*r','linewidth',1.1,'handlevisibility','off')
                        ylabel('Force (pN)','FontName','Serif')
                        xlabel('Thickness (nm)','FontName','Serif')

                        % xlim([220 500])

                        subplot(122)
                        hold on
                        title(['Stress-Strain, R2 = ' num2str(R2chad)])
                    %                                 plot(EpsEq*100,Sig,'ok','linewidth',1.5)
                    %                                 if exist('EpsFromFit','var')
                    %                                     plot(EpsFromFit,SigFromFit,'-b','linewidth',2)
                    %                                     clear EpsFromFit SigFromFit
                    %                                     legend('Experimantal data',['Linear fit - E = ' num2str(E0tmp)],'location','best')
                    %                                 end
                        plot(StrainChad*100,StressChad,'ok','linewidth',1.5)
                        xlabel('Strain (%)','FontName','Serif')
                        ylabel('Stress (Pa)','FontName','Serif')

                    %                                 xlim([0 25])
                    %                                 ylim([0 700])
                        if SAVING_CRITERION
                            fig = gcf;
                            saveas(fig,[sffMeca filesep currentCellName{i} '-' num2str(ii)],'png')
                            saveas(fig,[sffMeca filesep currentCellName{i} '-' num2str(ii)],'fig')
                            close
                        end
                    end
                    
                    if contains(currentCellName{i}, '07-08-20_M1_P1_C7') && ii == 3
                        figure
                        sgtitle(['Chadwick model fit, aSFL 3T3' newline  'Echad = ' num2str(Echadtmp, '%.2e') 'kPa'], 'fontsize', 30,'FontName','Serif')
                        subplot(121)
                        hold on
                        % title('Force-Thickness', 'fontsize', 18)
                        plot(Dup*1000,Fup,'-ob','linewidth',1.1)
                        plot(Ddown*1000,Fdown,'-or','linewidth',1.1)
                        plot(HcompFit*1000,Fcompfromfitchad,'-m','linewidth',3)
                        coeffPlotHF = ((3*H0chadFit)/(pi*EchadFit*(RADIUS/1000)))^(1/2);
                        plot((H0chadFit - coeffPlotHF*[0:5:Fcomp(1) Fcomp'].^(1/2))*1000,[0:5:Fcomp(1) Fcomp'],'--m','linewidth',1.5,'handlevisibility','off')
                        ylabel('Force (pN)','FontName','Serif')
                        xlabel('Thickness (nm)','FontName','Serif')
                        % xlim([220 500])

                        subplot(122)
                        hold on
                        % title('Stress & Strain from model', 'fontsize', 18)
                        plot(StrainChad*100,StressChad,'ok','linewidth',1.5)
                        xlabel('Strain (%)','FontName','Serif')
                        ylabel('Stress (kPa)','FontName','Serif')
                        % xlim([0 25])
                        % ylim([0 700])

                        fig = gcf; % Get Current Figure
                        saveas(fig,['C:\Users\JosephVermeil\Desktop\Figures' filesep currentCellName{i} '_CompressionNo' num2str(ii)],'png')
                        saveas(fig,['C:\Users\JosephVermeil\Desktop\Figures' filesep currentCellName{i} '_CompressionNo' num2str(ii)],'fig')
                    end

                    CompressionFeatures(firstCompressionIndex + ii).CellName = currentCellName{i};
                    CompressionFeatures(firstCompressionIndex + ii).CellIndex = i;
                    CompressionFeatures(firstCompressionIndex + ii).CompressionIndex = ii;
                    CompressionFeatures(firstCompressionIndex + ii).isAllowed = isAllowed;
                    CompressionFeatures(firstCompressionIndex + ii).isDecreasing = isDecreasing;
                    CompressionFeatures(firstCompressionIndex + ii).isCompressed = isCompressed;
                    CompressionFeatures(firstCompressionIndex + ii).isAGoodFit = isAGoodFit;
                    CompressionFeatures(firstCompressionIndex + ii).hasNoLargeJump = hasNoLargeJump;
                    CompressionFeatures(firstCompressionIndex + ii).hasSmallConfInt = hasSmallConfInt;
                    CompressionFeatures(firstCompressionIndex + ii).hasNormalDz = hasNormalDz;
                    CompressionFeatures(firstCompressionIndex + ii).SAVING_CRITERION = SAVING_CRITERION;
                    CompressionFeatures(firstCompressionIndex + ii).Echad = Echadtmp;
                    CompressionFeatures(firstCompressionIndex + ii).EchadCI = EchadFitCI;
                    CompressionFeatures(firstCompressionIndex + ii).H0 = H0chadFit*1000;
                    CompressionFeatures(firstCompressionIndex + ii).R2chad = R2chad;
                end
            end
            firstCompressionIndex = firstCompressionIndex + nComp;
            %             EchadWeight = (Echadplot./EchadCIvect).^2;
            %             EchadWeightedMean = sum(Echadplot.*EchadWeight)/sum(EchadWeight);
            %             EchadWeightedStd = sqrt(sum(((Echadplot-EchadWeightedMean).^2.*EchadWeight))/sum(EchadWeight));
            if length(Echad{i}) > 0
                Echad_Median{i} = median(Echad{i});
                Echad_Std{i} = std(Echad{i});
                Echad_Weight = (Echad{i}./EchadConf{i}).^2;
                Echad_WeightedMean{i} = sum(Echad{i}.*Echad_Weight)/sum(Echad_Weight);
                Echad_WeightedStd{i} = sqrt(sum(((Echad{i}-Echad_WeightedMean{i}).^2.*Echad_Weight))/sum(Echad_Weight));
                
            end
        end
    end
    
%     % Dernière ligne de CompressionFeatures pour faire le total
%     try
%         CompressionFeatures(firstCompressionIndex + 1).CellName = "Total";
%         CompressionFeatures(firstCompressionIndex + 1).isDecreasing = sum([CompressionFeatures.isDecreasing]);
%         CompressionFeatures(firstCompressionIndex + 1).isCompressed = sum([CompressionFeatures.isCompressed]);
%         CompressionFeatures(firstCompressionIndex + 1).isAGoodFit = sum([CompressionFeatures.isAGoodFit]);
%     catch
%     end
    
end

save([sf filesep 'D2M_' date '_' manip '_' specif '.mat'],'Echad','EchadConf','AlphaNL','Cort','H0prc','H0abs','ChadsH0','H0el','CompNum','Hys','Delt3','Delt20','CompressionFeatures',...
    'Echad_Median','Echad_Std','Echad_WeightedMean','Echad_WeightedStd')
% Deleted variables : 'E0','Edim',

end

