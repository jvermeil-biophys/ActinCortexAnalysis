function outvar = Data2MecaLoop_BigTable(date,manip,specif,tag,DIAMETER,resfolder,figurefolder,PLOT,VERBOSE)

% This version computes h0 first based on the first 15 points, then computes E for
% each data point in the compression using this h0. Then, on the longest
% portion where E is ~constant, the chadwick model is used to fit and get a
% new h0. This is done repeatedly until h0 and E do not change  anymore.

%% INIT

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');

R2CRITERION = 0.9;

SAVING_CRITERIONs = [];

% constante
RADIUS = DIAMETER/2000; % Rayon de la bille en �m

%% params and settings

outvar = [];

% home
% h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';
datenow = datestr(now,'yy-mm-dd');

% dossier de donn�es et sauvegarde des donn�es
% path = 'D:\Data\MatFiles\V2D';
% fieldpath =['D:\Data\Raw\' date(7:8) '.' date(4:5) '.' date(1:2)];
% sf   = ['D:\Data\MatFiles' filesep 'D2M'];
%fieldpath = [datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2)]; % INUTILE, -> datafolder supprim� des variable entr�e

path = [resfolder filesep 'V2D'];
sf = [resfolder filesep 'D2M'];


mkdir(sf)
% Param pour fit
LoopCont = 0;
LoopDef = 1; % conditionn� a loopcont = 0

ContTh = 1000; % en Pa
DefTh = 100; % en % de def

if LoopCont
    params = ['-Stress' num2str(ContTh) '-R2' num2str(R2CRITERION)];
elseif LoopDef
    params = ['-Strain' num2str(DefTh) '-R2' num2str(R2CRITERION)];
end


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
Echad = {}; % fit Chadwick
EchadConf = {}; % Confidence interval width around Echad value
ChadsH0 = {}; % h0 du fit chadwick
CompNum = {}; % numero de compression
StrainRateMed = {}; % stress rate median

%% Loading previous datas

if exist([resfolder filesep 'MecaDataTable.mat'])
    load([resfolder filesep 'MecaDataTable'])
else
    varnames = {'ExpType', 'ExpDay', 'CellID', 'CellName', 'CompNum'   , 'TpsComp', 'RawDataTDFB', 'MaxIndent' , 'MinThickness' , 'Hysteresis', 'CompTime'  , 'InitialThickness', 'SurroundingThickness', 'PreviousThickness', 'H0Chadwick', 'EChadwick', 'CiEChadwick', 'R2Chadwick'    , 'FitParams', 'Validated', 'Comments'};
    vartypes = {'string' , 'string', 'string', 'string'  , 'doublenan' , 'string' , 'cell'       , 'doublenan' , 'doublenan'    , 'doublenan' , 'doublenan' , 'doublenan'       , 'doublenan'           , 'doublenan'        , 'doublenan' , 'doublenan', 'doublenan'  , 'doublenan'     , 'string'   ,'logical'   , 'string'};
    BigTable = table('Size',[0 length(varnames)],'VariableTypes',vartypes,'VariableNames',varnames);
end

TableSize = size(BigTable,1);



%% Analysis
loadname=['V2D_' date '_' manip '_' tag '_' specif '.mat'];

if exist(loadname,'file')

        fprintf(['\nChargement du fichier ' loadname '\n\n'])

    
    load([path filesep loadname]);
    
    % for R248 exp
    if ~exist('MR')
        if exist('MR1')
            MR = MR1;
        elseif exist('MR2')
            MR = MR2;
        elseif exist('MR3')
            MR = MR3;
        end
    end
    
    nCells = length(MR); % nombre de cellules
    
    nCompOk = 0;
    nCompTot = 0;
    
    
    
    for i=1:nCells
        if not(isempty(MR{i}))
            % r�cup�ration des donn�es
            currentCellData = MR{i};
            currentCellName{i} = currentCellData.name;
            
            indexC = find(currentCellName{i} == 'C');
            indexM = find(currentCellName{i} == 'M');
            index_ = find(currentCellName{i} == '_');
            
            index_ = index_(index_>indexC);
            
            CellN(i) = str2num(currentCellName{i}(indexC+1:index_(1)-1));
            
            CellId = currentCellName{i}(indexM:index_(1)-1);
            
            D3 = currentCellData.D3*1000 - DIAMETER;
            dz = currentCellData.dz*1000;
            T = currentCellData.time;
            S = currentCellData.S;
            F = currentCellData.F;
            B = currentCellData.B;
            
            % Donn�e force cst [cst = constante]
            CstData = currentCellData.CstData;
            D3cst = CstData(:,3)*1000-DIAMETER;
            Tcst = CstData(:,2);
            D3smooth = smooth(Tcst,D3cst,7,'sgolay',5);
            
            % calcul du cortex et des centile
            medianCorticalThicknessCst = median(D3cst);
            halfPercentsGrid = [0:0.5:100];
            percentilesCorticalThicknessCst = prctile(D3cst,halfPercentsGrid);
            
            medianCortThick(i) = medianCorticalThicknessCst;
            
            % Donn�es rampe
            FullRampData = currentCellData.RampData{2}; % all data
            D3rmp = FullRampData(:,3)*1000 - DIAMETER;
            Trmp = FullRampData(:,2);
            
            RampData = currentCellData.RampData{1}; % good data
            nComp = length(RampData);
            nCompTot = nCompTot + nComp;
            
            CompDuration = currentCellData.RampData{3};
            
            CstRmp = currentCellData.RampData{4};
            
            
            
            %%%%%FIGURE%%%%%%
            if PLOT
                figcurve = figure;
                hold on
                
                yyaxis left
                axl = gca;
                axl.YColor = [0 0 0];
                ylabel('Cortex thickness (nm)','FontName','Serif')
                plot([T(1) T(end)],[0 0],'--k','handlevisibility','off')
                plot(T,D3,'k-','linewidth',1)
                plot(Tcst,D3cst,'ok','markerfacecolor','b')
                plot(Trmp,D3rmp,'o','markerfacecolor','r','markersize',5,'markeredgecolor','none')
                
                yyaxis right
                axr = gca; % Get Current Axes
                axr.YColor = [0 0 0];
                plot(T,B,'r-o','markersize',3)
                ylabel('Champ (mT)','FontName','Serif')
                
                if contains(specif,'M270')
                    ylim([0 800])
                else
                    ylim([0 2500])
                end
                
                title(currentCellName{i})
                
                xlabel('Time(s)','FontName','Serif')
                
                
                %                 close
            end
            
            %% Mecanique
            
            for ii = 1:nComp
                if not(ischar(RampData{ii}))
                    if VERBOSE
                        cprintf('-black', ['\n\nCompression : ' currentCellName{i} '-' num2str(ii) ' '])
                    end
                    
                    %%% [data saving] %%%
                    
                    % Check is this comp has already been analyzed with
                    % specified parameters
                    
                    CurrentCell = (strcmp(BigTable.ExpType,specif)) & (strcmp(BigTable.ExpDay,date)) & ...
                        (strcmp(BigTable.CellID,CellId)) & (BigTable.CompNum == ii) & ...
                        (strcmp(BigTable.FitParams,params(2:end))) & (strcmp(BigTable.TpsComp,CompDuration{ii})) ;
                    
                    AlreadyIn = sum(CurrentCell);
                    
                    
                    if AlreadyIn > 1
                        
                        error(['DuplicateCompression ! Check lines (' num2str(findAlreadyIn) ') in BigTable'])
                        
                    elseif AlreadyIn == 1
                        
                        CurrentTableLine = find(CurrentCell);
                        
                    else
                        CurrentTableLine = TableSize+1;
                        TableSize = TableSize+1;
                    end
                    
                    
                    % Replace existing or create new data line
                    BigTable(CurrentTableLine,'ExpType') = {specif};
                    BigTable(CurrentTableLine,'ExpDay') = {date};
                    BigTable(CurrentTableLine,'CellID') = {CellId};
                    BigTable(CurrentTableLine,'CellName') = {currentCellName{i}(1:index_(1)-1);};
                    BigTable(CurrentTableLine,'CompNum') = {ii};
                    BigTable(CurrentTableLine,'TpsComp') = {CompDuration{ii}};
                    BigTable(CurrentTableLine,'RawDataTDFB') = {[RampData{ii}(:,2:5)]};
                    BigTable(CurrentTableLine,'FitParams') = {params(2:end)};
                    
                    %%% [/data saving] %%
                    
                    
                    Btot = RampData{ii}(:,5); % Champ mag sur toute la dur�ee de la rampe
                    
                    minBtot = min(Btot);
                    
                    % Les lignes qui suivent servent � detecter le d�but de la rampe de champ magnetique.
                    % Pour Bamplitude = Bmax - Bmin, on place le d�but de la rampe � l'index pour lequel
                    % B = Bmin + Bamplitude/40 et la fin � l'index o� B = Bmax.
                    [maxBtot,indexMaxBtot] = max(Btot);
                    minBtot = min(Btot);
                    thresholdBtot = (maxBtot - minBtot) / 50; % variation de 2%
                    [~,indexThresholdBtot] = min(abs(Btot(1:indexMaxBtot) - (minBtot+thresholdBtot)));
                    
                    
                    % l'index bis correspond � la fin de la partie descendante de la rampe.
                    [~,indexThresholdBtotBis] = min(abs(Btot(indexMaxBtot:end) - (minBtot+thresholdBtot)));
                    indexThresholdBtotBis = indexThresholdBtotBis + indexMaxBtot - 1;
                    
                    
                    pointeur3 = Btot(1:indexMaxBtot) < minBtot+0.3 & Btot(1:indexMaxBtot) > minBtot-0.3; % Intersection d'ensembles
                    
                    % Definition of the 2 parts of the ramps
                    Dup = RampData{ii}(indexThresholdBtot:indexMaxBtot,3) - DIAMETER/1000; % Distance (en �m) pendant la phase de compression
                    Ddown = RampData{ii}(indexMaxBtot:indexThresholdBtotBis,3) - DIAMETER/1000; % Distance (en �m) pendant la phase de d�compression
                    
                    Fup = RampData{ii}(indexThresholdBtot:indexMaxBtot,4); % Force (en pN) pendant la phase de compression
                    Fdown = RampData{ii}(indexMaxBtot:indexThresholdBtotBis,4); % Force (en pN) pendant la phase de d�compression
                    
                    Tup = RampData{ii}(indexThresholdBtot:indexMaxBtot,2);
                    
                    hasEnoughPoints = length(Fup) > 8;
                    
                    %% H0 pour plot
                    
                    InitialThick = mean(RampData{ii}(pointeur3,3)) * 1000 - DIAMETER; % absolute thickness at the begining of compression
                    
                    %% Analyse mecanique
                    
                    % definition des variables et constantes
                    Hcomp = Dup; % �paisseur du cortex pendant la compression en �m
                    Fcomp = Fup; % force pendant la compression
                    Tcomp = Tup;
                    
                    differentialHcomp = diff(Hcomp(5:end));
                    MaxDisplacement = max(differentialHcomp)*1000; % en nm
                    
                    pointeurMonotonie =  differentialHcomp > 0; % 1 quand croissant 0 quand decroissant
                    
                    MaxDisplacementCummulative = MaxDisplacement; % deplacement maximum cumul� sur plus de trois points de suite
                    
                    GrowthPhases = bwconncomp(pointeurMonotonie); % get connected component where pointeurMonotonie == 1 <=> zones de croissance
                    % GrowthPhases.PixelIdxList{k} -> La liste des points de la k i�me zone de croissance.
                    
                    for k = 1:GrowthPhases.NumObjects % nombre de zones
                        if length(GrowthPhases.PixelIdxList{k}) > 2 % si plus de trois points croissants d'affil�e.
                            MaxDisplacementCummulative =  max([MaxDisplacementCummulative sum(differentialHcomp(GrowthPhases.PixelIdxList{k}))*1000]);
                            % on calcule de combien ca monte sur ces points, si c'est plus grand on note la valeur.
                        end
                    end
                    
                    % On considere que la courbe d'indentation est normalement d�croissante
                    % si il n'y a pas de deplacement de plus de 50nm en 1 pts, ou de remont�e
                    % de plus de 50nm sur 3 pts ou plus.
                    
                    % On teste �galement que la tendance globale de la
                    % courbe est d�croissante
                    
                    DecreaseFit = robustfit(Fcomp,Hcomp*1000);
                    
                    % QUESTION : Quels sont les cas qu'on vire exactement
                    % en faisant ce genre de trucs ???
                    
                    isDecreasing = (MaxDisplacementCummulative < 50) && (MaxDisplacement < 50) && DecreaseFit(2) < 0;
                    
                    
                    
                    %                         figure(1)
                    %                         hold on
                    %                         plot(Hcomp*1000,'g*-')
                    %                         a=1;
                    
                    
                    
                    % smoothing
                    FHsort = sortrows([Hcomp Fcomp Tcomp]); % B = sortrows(A) sorts the rows of a matrix in ascending order based on the elements in the first column.
                    Fcomp = FHsort(:,2);
                    Hcomp = FHsort(:,1);
                    Tcomp = FHsort(:,3);
                    Fcomp = smooth(Hcomp,Fcomp,7,'sgolay',5);
                    % Goal of 2 lines below : replace the potential negative values with a 0.
                    % Because otherwise the sqrt() a little below will create ugly complex numbers.
                    Fcomp_isPositive = [Fcomp > 0];
                    Fcomp = Fcomp .* Fcomp_isPositive;
                    %
                    Hcomp = Hcomp(end:-1:1);
                    Fcomp = Fcomp(end:-1:1);
                    Tcomp = Tcomp(end:-1:1);
                    
                    
                    % Recuperation des 10% d'enfoncement relatif au champ minimum pour estimer H0
                    pointeurFirst10pc = (InitialThick-Hcomp*1000)/InitialThick < 0.10;
                    pointeurFirst10pc(find(pointeurFirst10pc == 0, 1):end) = 0;
                    
                    nptsH0el = sum(pointeurFirst10pc);
                    
                    
                    F_H0 = Fcomp(pointeurFirst10pc);
                    H_H0 = Hcomp(pointeurFirst10pc);
                    
                    
                    % estimation de H0 par fit elastique lineaire
                    F_H0el = F_H0.^(1/2);
                    LinFitCoeff_H0 = polyfit(F_H0el, H_H0, 1);
                    H0el = LinFitCoeff_H0(2); % �paisseur de d�part
                    
                    
                    H0prev = 10e15;
                    H0tmp = H0el;
                    
                    cpt = 0;
                    
                    
                    
                    if VERBOSE
                        if LoopCont
                            fprintf(['\nSelection on STRESS (' num2str(ContTh) ' Pa). Npts inital loop : ' num2str(nptsH0el) '\n'])
                        elseif LoopDef
                            fprintf(['\nSelection on STRAIN (' num2str(DefTh) '%%). Npts inital loop : ' num2str(nptsH0el) '\n'])
                        end
                    end
                    
                    while abs((H0prev(end)-H0tmp)/H0prev(end)) > 0.01
                        
                        if LoopCont
                            %% utilisation du h0 pour calculer contrainte
                            Contrainte = Fcomp ./ (2*pi*RADIUS*(H0tmp - Hcomp));
                            
                            MaxCont = max(Contrainte);
                            
                            InferieurAXPourcents = Contrainte < ContTh;
                            
                        elseif LoopDef
                            %% utilisation du h0 pour calculer deformation
                            Deformation = ((H0tmp-Hcomp)/(3*H0tmp));
                            
                            MaxDef = max(Deformation);
                            
                            InferieurAXPourcents = Deformation*100 < DefTh;
                                                        
                        end
                        
                        ptrFit = find(InferieurAXPourcents == 0,1,'first');
                        
                        if isempty(ptrFit)
                            ptrFit = length(Fcomp);
                        end
                        
                        %% fit chadwick pour delta/2 et h/2
                        
                        fitchad = fittype('(pi*E*R*((h0-x)*(x<h0))^2)/(3*h0)','problem',{'R'});
                        
                        % condition on stress or strain
                        FcompFit = Fcomp(1:ptrFit);
                        HcompFit = Hcomp(1:ptrFit);
                        TcompFit = Tcomp(1:ptrFit);
                        
                        % H < H0
                        FcompFit(HcompFit>H0tmp) = [];
                        TcompFit(HcompFit>H0tmp) = [];
                        HcompFit(HcompFit>H0tmp) = [];
                        
%                         % validity of chadwick
%                         CtctRad = sqrt(RADIUS*(H0tmp - HcompFit));
%                         ptrChad = CtctRad>2*(HcompFit/2);
%                         
%                         FcompFit = FcompFit(ptrChad);
%                         HcompFit = HcompFit(ptrChad);
%                         TcompFit = TcompFit(ptrChad);
                        
                        
                        
                        try
                            
                            FFF = fit(HcompFit,FcompFit,fitchad,'problem',{RADIUS},'Start',[1000 H0tmp],... % [1000 H0el] [H0el 1000]
                                'Lower',[-Inf -Inf],'Upper',[Inf Inf]);
                            % Starting values : 1000 for FcompFit, H0tmpfor H0.
                            
                            FFFcoeff = coeffvalues(FFF);
                            
                            EchadFit = FFFcoeff(1); % module
                            H0chadFit = FFFcoeff(2); % FFFcoeff(2) est une nouvelle, meilleure estimation de H0.
                            
                            FFFci = confint(FFF);
                            
                            EchadFitCI = abs(diff(FFFci(:,1)));
                            %                     H0chadFitCI = FFFci(:,2);
                            
                            canBeFitted = 1;
                            
                            H0prev = [H0prev H0tmp];
                            H0tmp = H0chadFit;
                            
                        catch fitError
                            if VERBOSE
                                cprintf('Red', [' --> Cannot be fitted.\nError : ' fitError.message ])
                                
                            end
                            notsavedtext = 'Couldn''t be fitted';
                            canBeFitted = 0;
                            break
                        end
                        
                        cpt = cpt + 1;
                        doesConverge = 1;
                        
                        if cpt >= 16
                            if VERBOSE
                                cprintf('Red', ['-->  Does not converge.\n'])
                            end
                            doesConverge = 0;
                            break
                        end
                        
                    end
                else
                    canBeFitted = 0;
                    notsavedtext = ['Prob from V2D : ' RampData{ii}];
                    
                     %%% [data saving] %%%
                    
                   CurrentCell = (strcmp(BigTable.ExpType,specif)) & (strcmp(BigTable.ExpDay,date)) & ...
                        (strcmp(BigTable.CellID,CellId)) & (BigTable.CompNum == ii) & ...
                        (strcmp(BigTable.FitParams,params(2:end))) & (strcmp(BigTable.TpsComp,CompDuration{ii})) ;
                    
                    AlreadyIn = sum(CurrentCell);
                    
                    
                    if AlreadyIn > 1
                        
                        error(['DuplicateCompression ! Check lines (' num2str(findAlreadyIn) ') in BigTable'])
                        
                    elseif AlreadyIn == 1
                        
                        CurrentTableLine = find(CurrentCell);
                        
                    else
                        CurrentTableLine = TableSize+1;
                        TableSize = TableSize+1;
                    end
                    
                    
                    % Replace existing or create new data line
                    BigTable(CurrentTableLine,'ExpType') = {specif};
                    BigTable(CurrentTableLine,'ExpDay') = {date};
                    BigTable(CurrentTableLine,'CellID') = {CellId};
                    BigTable(CurrentTableLine,'CellName') = {currentCellName{i}(1:index_(1)-1);};
                    BigTable(CurrentTableLine,'CompNum') = {ii};
                    BigTable(CurrentTableLine,'TpsComp') = {CompDuration{ii}};
                    BigTable(CurrentTableLine,'RawDataTDFB') = {[RampData{ii}(:,2:5)]};
                    BigTable(CurrentTableLine,'FitParams') = {params(2:end)};
                    %%% [/data saving] %%
                    
                end
                
                if canBeFitted
                    
                    Fcompfromfitchad = FFF(HcompFit);
                    
                    StressChad = Fcomp(Hcomp<H0chadFit) ./ (pi*RADIUS*(H0chadFit - Hcomp(Hcomp<H0chadFit)));
                    
                    StrainChad = (H0chadFit - Hcomp(Hcomp<H0chadFit)) / (3*H0chadFit);
                    
                    R2chad = 1 - sum((FcompFit - Fcompfromfitchad).^2)/sum((FcompFit - mean(FcompFit)).^2);
                    
                    Echadtmp = EchadFit;
                    
                    isAGoodFit = (R2chad > R2CRITERION);
                    isPositive = (Echadtmp>0) && Echadtmp<150000;
                    
                    % Strain Rate computation
                    Tmid = Tup(1:end-1)+diff(Tup)/2;
                    Tinterp = sort([Tup; Tmid]);
                    Finterp = interp1(Tup,Fup,Tinterp,'spline');
                    Fdot = diff(Fup)./diff(Tup);
                    Fdotinterp = interp1(Tmid,Fdot,Tinterp,'spline');
                    
                    StrainRate = 0.5.*sqrt(1./(6.*pi.*RADIUS.*Finterp.*H0chadFit.*EchadFit)).*Fdotinterp; %
                    
                    
                    Deltadot  = diff(-smooth(Dup,8))./diff(Tup); % pas besoin de H0 pour la d�riv�e
                    
                    
                    ExpStrainRate = Deltadot/(3*H0chadFit);
                    
                    ExpStrainRateInterp = interp1(Tmid,ExpStrainRate,Tinterp,'spline');
                    
                    % hysteresis computation
                    
                    if length(Ddown)>1
                        Hyst = abs(trapz(Dup,Fup))-abs(trapz(Ddown,Fdown));
                    else
                        Hyst = NaN;
                    end
                    
                    %% Sauvegarde des valeurs pertinentes
%                     
                    if length(HcompFit)>3
                        orderedMat = sortrows([TcompFit HcompFit]);
                        dhcomp = abs(diff(orderedMat(:,2)))*1000;
                        d3hcomp = dhcomp(1:end-2)+dhcomp(2:end-1)+dhcomp(3:end);
                        hasNoLargeJump = ((max(dhcomp) < max([0.15*(max(Dup*1000)-min(Dup*1000)) 15])) ...
                            && (max(d3hcomp)< max([0.30*(max(Dup*1000)-min(Dup*1000)) 30]))) || strcmp(specif(end-2:end),'02s') ...
                            || strcmp(specif(end-2:end),'05s');
                    else
                        hasNoLargeJump = 1;
                    end
                    
                    hasSmallConfInt = (abs(EchadFitCI) < 2*abs(Echadtmp));
                    
                    currentRampDz = RampData{ii}(:,6);
                    hasNormalDz = (max(abs(currentRampDz)) < 1);
                    
                    SAVING_CRITERION =  isPositive && hasNormalDz && hasSmallConfInt && hasNoLargeJump && hasEnoughPoints && doesConverge && isAGoodFit;
                    
                    SAVING_CRITERIONs =  [isPositive hasNormalDz hasSmallConfInt hasNoLargeJump hasEnoughPoints doesConverge isAGoodFit];
                    
                    
                    if SAVING_CRITERION
                        if VERBOSE
                            cprintf('*green', ['Successfully saved !'])
                        end
 
                    else
                        %% display of problems
                        
                        nreason = sum(~SAVING_CRITERIONs);
                        
                        
                        notsavedtext = [];
                        
                        if ~isPositive
                            
                            notsavedtext = [notsavedtext 'E is not Positive '];
                            
                            nreason = nreason - 1;
                            
                            
                            if nreason >0
                                
                                notsavedtext = [notsavedtext '& '];
                                
                            end
                            
                        end
                        
                        
                        
                        if ~hasNormalDz
                            
                            notsavedtext = [notsavedtext 'Disaligned vertically '];
                            
                            nreason = nreason - 1;
                            
                            if nreason >0
                                
                                notsavedtext = [notsavedtext '& '];
                                
                            end
                            
                        end
                        
                        
                        
                        
                        
                        if ~hasSmallConfInt
                            
                            notsavedtext = [notsavedtext 'Has too big CI on E '];
                            
                            nreason = nreason - 1;
                            
                            if nreason >0
                                
                                notsavedtext = [notsavedtext '& '];
                                
                            end
                            
                        end
                        
                        
                        if ~hasNoLargeJump
                            
                            notsavedtext = [notsavedtext 'Has large jumps '];
                            
                            nreason = nreason - 1;
                            
                            if nreason >0
                                
                                notsavedtext = [notsavedtext '& '];
                                
                            end
                            
                        end
                        
                        
                        
                        if ~hasEnoughPoints
                            
                            notsavedtext = [notsavedtext 'Not enough points '];
                            
                            nreason = nreason - 1;
                            
                            if nreason >0
                                
                                notsavedtext = [notsavedtext '& '];
                                
                            end
                            
                        end
                        
                        
                        if ~doesConverge
                            
                            notsavedtext = [notsavedtext 'No Convergence on H0 '];
                            
                            nreason = nreason - 1;
                            
                        end
                        
                        
                        if ~isAGoodFit
                            
                            notsavedtext = [notsavedtext 'Bad R2 on fit '];
                            
                            nreason = nreason - 1;
                            
                        end
                        
                        
                        if VERBOSE
                            cprintf('Red', ['Could not be saved because : ' notsavedtext])
                        end
                    end
                    
                    %% %%%%%%% FIGURES %%%%%%
                    if PLOT
                        
                        figure(figcurve);
                        hold on
                        yyaxis left
                        if SAVING_CRITERION
                            
                            plot(RampData{ii}(:,2),RampData{ii}(:,3)*1000-DIAMETER,'o','markerfacecolor','g','markersize',6,'markeredgecolor','none')
                            
                        else
                            
                            text(mean(RampData{ii}(:,2)),(max(RampData{ii}(:,3)*1000-DIAMETER)*1.2),...
                                notsavedtext,...
                                'HorizontalAlignment','center')
                        end
                        
                        saveas(figcurve,[sffDistCurve filesep currentCellName{i}],'png')
                        saveas(figcurve,[sffDistCurve filesep currentCellName{i}],'fig')
                        
                        figure
                        if SAVING_CRITERION
                            sgtitle(['Chadwick model for ' currentCellName{i} ' - ' num2str(ii)])
                        else
                            sgtitle(['NOT SAVED DATA - Chadwick model for ' currentCellName{i} ' - ' num2str(ii)])
                        end
                        subplot(221)
                        hold on
                        title(['Force-Thickness, Echad = ' num2str(Echadtmp/1000,'%.1f') 'kPa'])
                        plot(Dup*1000,Fup,'-*b','linewidth',1.1)
                        plot(Ddown*1000,Fdown,'-*r','linewidth',1.1)
                        plot(HcompFit*1000,Fcompfromfitchad,'-m','linewidth',3)
                        plot(Hcomp*1000,FFF(Hcomp),'xm','linewidth',1.5,'handlevisibility','off')
                        ylabel('Force (pN)','FontName','Serif')
                        xlabel('Thickness (nm)','FontName','Serif')
                        
                        % xlim([220 500])
                        
                        subplot(222)
                        hold on
                        title(['Stress-Strain, ChadR2 = ' num2str(R2chad)])
                        plot(StrainChad*100,StressChad,'ok','linewidth',1.5)
                        xlabel('Chadwick Strain (%)','FontName','Serif')
                        ylabel('Chadwick Stress (Pa)','FontName','Serif')
                        
                        %                                 xlim([0 25])
                        %                                 ylim([0 700])
                        
                        
                        subplot(223)
                        hold on
                        title('Evolution of H0 with loops')
                        plot(0:cpt,[H0prev(2:end) H0tmp]*1000,'ko-','markerfacecolor','c')
                        xlabel('Loop n�','FontName','Serif')
                        ylabel('H0 (nm)','FontName','Serif')
                        
                        
                        subplot(224)
                        hold on
                        title('Evolution of strain rate in time')
                        plot(Tinterp-Tinterp(1),StrainRate*100,'-ko','markerfacecolor','m')
                        plot(Tinterp-Tinterp(1),ExpStrainRateInterp*100,'-ks','markerfacecolor','r')
                        
                        xlabel('Time during compression (s)','FontName','Serif')
                        ylabel('Strain Rate ($s^{-1}$)','FontName','Serif','interpreter','latex')
                        legend('theoretical','experimental')
                        legend('location','best')
                        
                        yl = ylim;
                        ylim([-5 yl(2)]);
                        
                        fig = gcf;
                        saveas(fig,[sffMeca filesep currentCellName{i} '-' num2str(ii)],'png')
                        saveas(fig,[sffMeca filesep currentCellName{i} '-' num2str(ii)],'fig')
                        
                        if ~SAVING_CRITERION
                            close
                        end
                    end
                    
                    
                else
                                    
                    if PLOT
                        figure(figcurve);
                        yyaxis left
                        text(mean(RampData{ii}(:,2)),max(RampData{ii}(:,3))*1000-DIAMETER,...
                            'CANTBEFITTED',...
                            'HorizontalAlignment','center')
                    end
                end
                
                
                
                
                %%% [data saving] %%%                               
                
                if canBeFitted && SAVING_CRITERION
                    
                    BigTable(CurrentTableLine,'SurroundingThickness') = {(median(CstRmp{ii}{1})*1000-DIAMETER)};
                    BigTable(CurrentTableLine,'PreviousThickness') = {(median(CstRmp{ii}{2})*1000-DIAMETER)};
                    BigTable(CurrentTableLine,'InitialThickness') = {InitialThick};
                    BigTable(CurrentTableLine,'CompTime') = {mean(Tup)};
                    BigTable(CurrentTableLine,'MaxIndent') = {InitialThick-min(Dup)*1000};
                    BigTable(CurrentTableLine,'MinThickness') = {min(Dup)*1000};
                    BigTable(CurrentTableLine,'H0Chadwick') = {H0chadFit*1000};
                    BigTable(CurrentTableLine,'EChadwick') = {Echadtmp};
                    BigTable(CurrentTableLine,'CiEChadwick') = {EchadFitCI};
                    BigTable(CurrentTableLine,'R2Chadwick') = {R2chad};      
                    BigTable(CurrentTableLine,'Hysteresis') = {Hyst};               
                    
                    
                    BigTable(CurrentTableLine,'Validated') = {1};
                    BigTable(CurrentTableLine,'Comments') = {'Succes!'};
                    
                    nCompOk = nCompOk +1;
                    
                elseif canBeFitted
                    
                    BigTable(CurrentTableLine,'SurroundingThickness') = {(median(CstRmp{ii}{1})*1000-DIAMETER)};
                    BigTable(CurrentTableLine,'PreviousThickness') = {(median(CstRmp{ii}{2})*1000-DIAMETER)};
                    BigTable(CurrentTableLine,'InitialThickness') = {InitialThick};
                    BigTable(CurrentTableLine,'CompTime') = {mean(Tup)};
                    BigTable(CurrentTableLine,'MaxIndent') = {InitialThick-min(Dup)*1000};
                    BigTable(CurrentTableLine,'MinThickness') = {min(Dup)*1000};
                    BigTable(CurrentTableLine,'H0Chadwick') = {H0chadFit*1000};
                    BigTable(CurrentTableLine,'EChadwick') = {Echadtmp};
                    BigTable(CurrentTableLine,'CiEChadwick') = {EchadFitCI};
                    BigTable(CurrentTableLine,'R2Chadwick') = {R2chad};    
                    BigTable(CurrentTableLine,'Hysteresis') = {Hyst};               
                   
                    BigTable(CurrentTableLine,'Validated') = {0};
                    BigTable(CurrentTableLine,'Comments') = {notsavedtext};
                    
                else
                    
                    BigTable(CurrentTableLine,'SurroundingThickness') = {NaN};
                    BigTable(CurrentTableLine,'PreviousThickness') = {NaN};
                    BigTable(CurrentTableLine,'InitialThickness') = {NaN};
                    BigTable(CurrentTableLine,'CompTime') = {NaN};
                    BigTable(CurrentTableLine,'MaxIndent') = {NaN};
                    BigTable(CurrentTableLine,'MinThickness') = {NaN};
                    BigTable(CurrentTableLine,'H0Chadwick') = {NaN};
                    BigTable(CurrentTableLine,'EChadwick') = {NaN};
                    BigTable(CurrentTableLine,'CiEChadwick') = {NaN};
                    BigTable(CurrentTableLine,'R2Chadwick') = {NaN};    
                    BigTable(CurrentTableLine,'Hysteresis') = {NaN};               
                   
                    
                    BigTable(CurrentTableLine,'Validated') = {0};
                    BigTable(CurrentTableLine,'Comments') = {notsavedtext};
                end
                %%% [/data saving] %%%
                
                
            end
        end
    end   
end

save([resfolder filesep 'MecaDataTable'],'BigTable');

pctgsaved = nCompOk/nCompTot*100;

cprintf('*Blue', ['\n\n' num2str(round(pctgsaved*10)/10) '%% of compressions validated.\n\n\n'])

end

