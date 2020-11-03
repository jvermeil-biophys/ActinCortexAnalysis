function Data2Meca_R248(date,manip,specif,Diam,PLOT)

%% INIT

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');

here = pwd;

R2sig = 0.9;

PLOTTMP = 0;

%% params and settings

% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';


% dossier de données et sauvegarde des données
path = 'D:\Data\MatFile\V2D';
fieldpath =['D:\Data\Raw\' date(7:8) '.' date(4:5) '.' date(1:2)];
sf   = ['D:\Data\MatFile' filesep 'D2M'];

datenow = datestr(now,'yy-mm-dd');

mkdir(sf)

scf = 'C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Données';

mkdir([scf filesep datenow filesep 'D2F'])

% sauvergarde des figures
sfig = 'D:\Data\Figures';
sff = [sfig filesep datenow filesep 'Meca'];
mkdir(sff)

if PLOT
    
    
    sffd = [sff filesep specif filesep 'DistCurve'];
    mkdir(sffd)
    
    sffm = [sff filesep specif filesep 'ForceDist&StressStrain'];
    mkdir(sffm)
end
df = ['C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Figure' filesep 'Meca' filesep specif];

mkdir(df)


Cort = [];

%% Retrieve datafile
loadname=['V2D_' date '_' manip '_R248_' specif '.mat'];

cd(path)

if exist(loadname,'file')
    
    fprintf(['Chargement du fichier ' loadname '\n\n'])
    load(loadname);
    
    
    nCells = length(MR); % nombre de cellules
    
    cptfit = 0;
    
    
    %% Get datas
    for i=1:nCells
        if not(isempty(MR{i}))
            H0abs{1,i} = []; % absolu
            H0prc{1,i} = []; % en percentile de la courbe D3cst
            E0{1,i} = []; % elasticité sur le debut des deformation
            E16mT{1,i} = []; % fit dimitriadis sur le debut des deformation
            AlphaNL{1,i} = []; % exposant non lineaire
            CompNum{1,i} = [];
            Hys{1,i} = [];
            finalSig{1,i}  = [];
            
            
            H0abs{2,i} = []; % absolu
            H0prc{2,i} = []; % en percentile de la courbe D3cst
            E0{2,i} = []; % elasticité sur le debut des deformation
            E16mT{2,i} = []; % fit dimitriadis sur le debut des deformation
            AlphaNL{2,i} = []; % exposant non lineaire
            CompNum{2,i} = [];
            Hys{2,i} = [];
            finalSig{2,i}  = [];
            
            
            H0abs{3,i} = []; % absolu
            H0prc{3,i} = []; % en percentile de la courbe D3cst
            E0{3,i} = []; % elasticité sur le debut des deformation
            E16mT{3,i} = []; % fit dimitriadis sur le debut des deformation
            AlphaNL{3,i} = []; % exposant non lineaire
            CompNum{3,i} = [];
            Hys{3,i} = [];
            finalSig{3,i}  = [];
            
            
            
            % récupération des données
            
            CellDat = MR{i};
            
            name{i} = CellDat.name;
            
            iC = find(name{i} == 'C');
            iUs = find(name{i} == '_');
            
            iUs = iUs(iUs>iC);
            
            CellN(i) = str2num(name{i}(iC+1:iUs(1)-1));
            
            D3 = CellDat.D3*1000 - Diam;
            dz = CellDat.dz;
            T = CellDat.time;
            F = CellDat.F;
            
            % Donnée force cst
            CstData = CellDat.CstData;
            D3cst = CstData(:,3)*1000-Diam;
            Tcst = CstData(:,2);
            D3s = smooth(Tcst,D3cst,7,'sgolay',3);
            
            
            % calcul du cortex et des centile
            CortWCell = median(D3cst);
            grid = [0:0.5:100];
            PRCT = prctile(D3cst,grid);
            
            Cort(i) = CortWCell;
            
            
            % Donnée rampe
            FullRampData = CellDat.RampData{2}; % all data
            D3rmp = FullRampData(:,3)*1000-Diam;
            Trmp = FullRampData(:,2);
            
            
            RampData = CellDat.RampData{1}; % good data
            
            nComp = length(RampData);
            
            
            
            %%%%%FIGURE%%%%%%
            if PLOT && 0
                figure
                hold on
                
                yyaxis left
                axl = gca;
                axl.YColor = [0 0 0];
                ylabel('Cortex thickness (nm)','FontName','Serif')
                plot(T,D3,'k-','linewidth',1)
                plot(Tcst,D3cst,'ok','markerfacecolor','b')
                plot(Trmp,D3rmp,'o','markerfacecolor',[0.7 0.1 0.5],'markersize',4,'markeredgecolor','none')
                
                
                yyaxis right
                axr = gca;
                axr.YColor = [0 0 0];
                plot(T,F)
                ylabel('Force (pN)','FontName','Serif')
                title(name{i})
                
                xlabel('Time(s)','FontName','Serif')
                
                fig = gcf;
                saveas(fig,[sffd filesep name{i}],'png')
                close
            end
            
            %% Mecanique
            
            nCompOk = 0;
            for ii = 1:nComp
                if not(ischar(RampData{ii}))
                    
                    rmpType = mod(ii,3);
                    
                    rmpType(rmpType == 0) = 3;
                    
                    Btot = RampData{ii}(:,5); % Champ mag sur toute la duréee de la rampe
                    [~,Mind] = max(Btot);
                    
                    Bup = Btot(1:Mind);
                    
                    Dup = RampData{ii}(1:Mind,3)-Diam/1000; % Distance (en µm) pendant la phase de compression
                    Ddown = RampData{ii}(Mind:end,3)-Diam/1000; % Distance (en µm) pendant la phase de décompression
                    
                    lowD = min(Dup)-0.03;
                    
                    if lowD < 0
                        Dup = Dup - lowD; % réhaussement des comp qui
                        Ddown = Ddown - lowD;% passent en dessous de 30 nm
                    end
                    
                    Fup = RampData{ii}(1:Mind,4); % Force (en pN) pendant la phase de compression
                    Fdown = RampData{ii}(Mind:end,4); % Force (en pN) pendant la phase de décompression
                    
                    %% H0 pour plot
                    
                    AbsThick = max(Dup)*1000; % absolute thickness at the begining of compression
                    
                    [~,ind] = min(abs(PRCT-AbsThick));
                    prct = grid(ind); % percentile of data
                    
                    
                    prct = (AbsThick-min(D3cst))/(max(D3cst)-min(D3cst))*100; % percentage of amplitude
                    
                    
                    %% Analyse mecanique
                    
                    % definition des variables et constantes
                    Hcomp = Dup; % épaisseur du cortex pendant la compression en µm
                    Fcomp = Fup; % force pendant la compression
                    
                    DH = diff(Hcomp(5:end));
                    MaxDisp = max(DH)*1000; % en nm
                    
                    ptrmono =  DH > 0; % 1 quand croissant 0 quand decroissant
                    
                    MaxDispCum = MaxDisp; % deplacement maximum cumulé sur plus de trois points de suite
                    
                    CC = bwconncomp(ptrmono); % get connected component = zone de croissance
                    
                    for k = 1:CC.NumObjects % nombre de zone
                        if length(CC.PixelIdxList{k})>2 % si plus de trois points croissant d'affilé
                            MaxDispCum =  max([MaxDispCum sum(DH(CC.PixelIdxList{k}))*1000]);
                            % on calcul de combien ca monte sur ces points,
                            % si c'est plus grand on note la valeur
                        end
                    end
                    
                    % on considere "decroissant" si il n'y a pas de
                    % deplacement
                    % de plus de 50nm en 1 pts, ou de remonté de plus de
                    % 50nm sur 3pts ou plus
                    %
                    
                    isDecreasing = (MaxDispCum < 50) & (MaxDisp < 50);
                    
                    if isDecreasing
                        
                        %                         figure(1)
                        %                         hold on
                        %                         plot(Hcomp*1000,'g*-')
                        %                         a=1;
                        
                        % constante
                        Rb = Diam/2000; % Rayon de la bille en µm
                        
                        % smoothing
                        FHsort = sortrows([Hcomp Fcomp]);
                        Fcomp = FHsort(:,2);
                        Hcomp = FHsort(:,1);
                        Fcomp = smooth(Hcomp,Fcomp,7,'sgolay',3);
                        
                        Hcomp = Hcomp(end:-1:1);
                        Fcomp = Fcomp(end:-1:1);
                        
                        
                        % Recuperation des 150 premier pN
                        % d'enfoncement pour estimer H0
                        ptrh0 = (Fcomp<min(Fcomp)+150);
                        
                        Fh0 = Fcomp(ptrh0);
                        Hh0 = Hcomp(ptrh0);
                        
                        % estimation de h0 par fit elastique lineaire
                        Fh0el = Fh0.^(1/2);
                        Pel = polyfit(Fh0el,Hh0,1);
                        H0el = Pel(2); % épaisseur de départ
                        
                        if all(H0el > Hcomp)
                            %% fit maison autour de 16mT
                            
                            fitfunc = fittype(['pi*E*R*(x0-x)^2/(2*x0)']...
                                ,'problem',{'R'}); % fonction a fitter , R et rh paramètre, Ef et x0 a fitter
                            
                            
                            E16mTtmp = [];
                            
                            % on fite sur autour de 16 mT (petite bobine
                            % avec petite aimants)
                            
                           
                            FcompFit =  Fup(Bup>26&Bup<36);
                            HcompFit = Dup(Bup>26&Bup<36); % en Pa
                            
                            
                            FFF =fit(HcompFit,FcompFit,fitfunc,'problem',{Rb},'Start',[1000,H0el]);
                            
                            FFFc=coeffvalues(FFF);
                            
                            E16mTFit=FFFc(1); % module
                            H0fitFit=FFFc(2); % épaisseur de départ
                            
                            HcompFromFit = HcompFit*1000;
                            FcompFromFit = FFF(HcompFit);
                            
                            R216mT = 1 - sum((FcompFit - FcompFromFit).^2)/sum((FcompFit - mean(FcompFit)).^2);
                            

                                E16mTtmp = E16mTFit;
                                R216mTtmp = R216mT;
                                
                                
                                
                                
            
                            
%                             figure
%                             hold on
%                             plot(Dup*1000,Fup,'-*b')
%                             plot(HcompFit*1000,FcompFit,'-og')
%                             plot(HcompFromFit,FcompFromFit,'--m','linewidth',1.5)
%                             close
                            
                            
                            
                            
                            %% fit chadwick-like
                            % calcul de l'indentation maximum delta
                            Delt = H0el - Hcomp; % indentation en µm
                            
                            % calcul d'epsilon et sigma
                            EpsEq = Delt/(2*H0el); % delta equivalent = delta / 2 -> volume equivalent compresser par un cylindre
                            Sig = Fcomp./(pi*Rb*Delt); % stress en Pa
                            
                            
                            % calcul de l'elasticité a petite deformation
                            
                            E0tmp = [];
                            
                            % on fite sur 100Pa, puis 150, puis 200, etc...
                            
                            step = 50;
                            R2tmp = 0;
                            old_pos = 1;
                            
                            minSig = min(Sig);
                            
                            while max(Sig)> minSig+step+50
                                
                                step = step + 50;
                                
                                [~,pos] = min(abs(Sig(old_pos:end)-(minSig+step)));
                                
                                pos = max(pos+old_pos-1,5);
                                
                                EpsFit = EpsEq(1:pos);
                                SigFit = Sig(1:pos); % en Pa
                                
                                old_pos = pos;
                                
                                
                                [Rfit,stats] = robustfit(EpsFit,SigFit);
                                
                                E0fit = Rfit(2); % en Pa
                                
                                Sigfromfit = Rfit(1) + Rfit(2)*EpsFit;
                                
                                R2 = 1 - sum(stats.w.*(SigFit - Sigfromfit).^2)/sum(stats.w.*(SigFit - mean(SigFit)).^2);
                                
                                if R2 > R2sig && R2 > R2tmp && E0fit > 0
                                    E0tmp = E0fit;
                                    R2tmp = R2;
                                    
                                    finalSigtmp = minSig+step;
                                    
                                    EpsFromFit = EpsFit*100;
                                    SigFromFit = Rfit(1) + Rfit(2)*EpsFit;
                                    
                                    
                                end
                                
                            end
                            
                            R2 = R2tmp;
                            
                            
                            %% %%%%%%% FIGURE %%%%%%
                            if PLOT
                                figure
                                hold on
                                title(['Force-Thickness and stress-strain _ ' num2str(2^rmpType) 's Comp'])
                                plot(Dup*1000,Fup,'-*b','linewidth',1.1)
                                plot(Ddown*1000,Fdown,'-*r','linewidth',1.1)
                                plot((Pel(2)+Pel(1)*[Fh0'].^(1/2))*1000,[Fh0'],'--m','linewidth',1.5)
                                
                                
                                plot(Pel(2)*1000,0,'*m','linewidth',3)
                                plot(Dup*1000,Fup,'-*b','linewidth',1.1,'handlevisibility','off')
                                plot(Ddown*1000,Fdown,'-*r','linewidth',1.1,'handlevisibility','off')
                                ylabel('Force (pN)')
                                xlabel('Thickness (nm)','FontName','Serif')
                                
                                
                                legend('Compression','Relaxation','Elastic fit','Fited H0','FontName','Serif')
                                
                                
                                
                                fig = gcf;
                                saveas(fig,[sffm filesep name{i} '-' num2str(ii)],'png')
                                close
                                
                                
                            end
                            
                            
                            %% calcul de l'exposant non lineaire
                            
                            Lesp = log(EpsEq);
                            Lsig = log(Sig);
                            
                            Pnl = polyfit(Lesp,Lsig,1);
                            
                            Lsigff = Pnl(2) + Pnl(1)*Lesp;
                            
                            R2al = 1 - sum((Lsig - Lsigff).^2)/sum((Lsig - mean(Lsig)).^2);
                            
                            
                            if (R2 > R2sig) && (E0tmp>0) && (E0tmp<30000) % && (R2dim > R2sig)
                                
                                                 
                                % on prend les valeurs pour plot
                                
                                H0abs{rmpType,i} = [H0abs{rmpType,i}; AbsThick]; % valeur de l'épaisseur du cortex avant compression
                                H0prc{rmpType,i} = [H0prc{rmpType,i}; prct]; % percentile correspondant
                                
                                E0{rmpType,i} = [E0{rmpType,i}; E0tmp];
                                E16mT{rmpType,i} = [E16mT{i}; E16mTtmp];
                                
                                finalSig{rmpType,i} = [finalSig{rmpType,i}; finalSigtmp];
                                
                                CompNum{rmpType,i} = [CompNum{rmpType,i}; ii];
                                
                                
                                
                                AlphaNL{rmpType,i} = [AlphaNL{rmpType,i}; Pnl(1)];
                                
                                
                                %% calcul de l'hysteresis
                                
                                AreaComp = abs(trapz(Dup,Fup));
                                AreaRel = abs(trapz(Ddown,Fdown));
                                Hys{rmpType,i} = [Hys{rmpType,i}; (AreaComp-AreaRel)/AreaComp*100];
                                
                                nCompOk = nCompOk +1;
                                
                                
                            end
                        end
                    else
                        %                         figure(1)
                        %                         hold on
                        %                         plot(Hcomp*1000,'r*-')
                        %                         a=1;
                    end
                end
            end
        end
    end
end

save([sf filesep 'D2M_' date '_' manip '_' specif '.mat'],'E0','E16mT','finalSig','AlphaNL','Cort','H0prc','H0abs','CompNum','Hys')


end


