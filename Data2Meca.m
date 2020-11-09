function [nCompTot, nCompOk] = Data2Meca(date,manip,specif,tag,Diam,PLOT,resfolder,figfolder)

%% INIT

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',15)
here = pwd;

R2sig = 0;

PLOTTMP = 0;

%% params and settings



% dossier de données et sauvegarde des données
path = [resfolder filesep 'V2D'];
sf   = [resfolder filesep 'D2M'];

datenow = datestr(now,'yy-mm-dd');

mkdir(sf)


% sauvergarde des figures
sfig = figfolder;
sff = [sfig filesep datenow filesep 'Meca'];
mkdir(sff)

if PLOT
    
    
    sffd = [sff filesep specif filesep 'DistCurve'];
    mkdir(sffd)
    
    sffm = [sff filesep specif filesep 'ForceDist&StressStrain'];
    mkdir(sffm)
    
    sfffit = [sff filesep specif filesep 'Fits'];
    mkdir(sfffit)
end


H0abs = {}; % absolu
DimH0 = {}; % en percentile de la courbe D3cst
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
    
    cptfit = 0;
    
    nCompOk = 0;
    nCompTot = 0;
    
    %% Get datas
    for i=1:nCells
        if not(isempty(MR{i}))
            H0abs{i} = []; % absolu
            DimH0{i} = []; % en percentile de la courbe D3cst
            E0{i} = []; % elasticité sur le debut des deformation
            Edim{i} = []; % fit dimitriadis
            Echad{i} = []; % fit dchadwick
            EchadConf{i} = []; % fit dchadwick confidence interval
            MaxDef{i} = []; % deformation maximum
            MinThick{i} = []; % epaisseur minimum pdt comp
            Delt5{i} = []; % Enfoncement a 5mT
            Delt20{i} = []; % Enfoncement theorique a 20 pN
            ChadsH0{i} = [];
            AlphaNL{i} = []; % exposant non lineaire
            CompNum{i} = [];
            Hys{i} = [];
            
            
            
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
            
            nCompTot = nCompTot + nComp;
            
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
            
            
            for ii = 1:nComp
                if not(ischar(RampData{ii}))
                    
                    Btot = RampData{ii}(:,5); % Champ mag sur toute la duréee de la rampe
                    [~,Mind] = max(Btot);
                    
                    ptr45 = Btot(1:Mind)<4.6 & Btot(1:Mind)>4.4; % zone a 4.5 mT de champ
                    
                    
                    Dup = RampData{ii}(1:Mind,3)-Diam/1000; % Distance (en µm) pendant la phase de compression
                    Ddown = RampData{ii}(Mind:end,3)-Diam/1000; % Distance (en µm) pendant la phase de décompression
                    
                    
                    Dup = Dup;
                    
                    
%                     lowD = min(Dup)-0.03;
                    %
                    %                     if lowD < 0
                    %                         Dup = Dup - lowD; % réhaussement des comp qui
                    %                         Ddown = Ddown - lowD;% passent en dessous de 30 nm
                    %                     end
                    
                    Fup = RampData{ii}(1:Mind,4); % Force (en pN) pendant la phase de compression
                    Fdown = RampData{ii}(Mind:end,4); % Force (en pN) pendant la phase de décompression
                    
                    %% H0 pour plot
                    
                    AbsThick = mean(Dup(1:3))*1000; % absolute thickness at the begining of compression

                    
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
                        
                                             
                        Hcomp = Hcomp (end:-1:1);
                        Fcomp = Fcomp (end:-1:1);
                        
                        
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
                            
                            % calcul de la deformation maximum
                            
                            deltamax = H0el - min(Hcomp);
                            maxdef = deltamax/H0el;
                            
                            minthick = min(Hcomp);
                            
                            if minthick > 0
                                %% fit dimitriadis
                                
                                Hdep = AbsThick/1000;
                                
                                fitdim = fittype(['Ef*16/9*sqrt(R)*(max(x0-x,0).^1.5+0.884*rh*max(x0-x,0).^2+0.781*' ...
                                    'rh^2*max(x0-x,0).^2.5+0.386*rh^3*max(x0-x,0).^3+0.0048*rh^4*max(x0-x,0).^3.5)']...
                                    ,'problem',{'R','rh'}); % fonction a fitter (dimitriadis), R et rh paramètre, Ef et x0 a fitter
                                
                                
                                
                                
                                % on fit de tant que delta/h < 0.25
                                FcompFit = Fcomp((Hdep-Hcomp)/Hdep<0.25);
                                HcompFit = Hcomp((Hdep-Hcomp)/Hdep<0.25); % en Pa
                                
                                try
                                FFFdim =fit(HcompFit,FcompFit,fitdim,'problem',{Rb, sqrt(Rb)/H0el},'Start',[1000,H0el],...
                                    'Lower',[1 -Inf],'Upper',[Inf Inf]);
                                catch
                                    themall
                                end 
                                FFFc=coeffvalues(FFFdim);
                                
                                EdimFit=FFFc(1); % module
                                H0dimFit=FFFc(2); % épaisseur de départ
                                
                                Fcompfromfit = FFFdim(HcompFit);
                                
                                R2dim = 1 - sum((FcompFit - Fcompfromfit).^2)/sum((FcompFit - mean(FcompFit)).^2);
                                
                                
                                Edimtmp = EdimFit;
                                R2dimtmp = R2dim;
                                
                                HcompFromFit = HcompFit*1000;
                                FcompFromFit = FFFdim(HcompFit);
                                
                                
                                
                                R2dim = R2dimtmp;
                                
                                %% fit chadwick pour delta/2 et h/2
                                
                                fitchad = fittype(['(pi*E*R*(x-h0)^2)/(3*h0)']...
                                    ,'problem',{'R'}); % fonction a fitter (chadwick), R donné , E et h0 a fitter
                                
                                
                                
                                
                                
                                H45 = mean(Dup(ptr45));
                                
                                % on fit de tant que delta/h < 0.25
                                FcompFit = Fcomp((Hdep-Hcomp)/Hdep<0.25);
                                HcompFit = Hcomp((Hdep-Hcomp)/Hdep<0.25); % en Pa
                                
                                
                                
                                try
                                    FFFchad =fit(HcompFit,FcompFit,fitchad,'problem',{Rb},'Start',[1000 H0el],...
                                        'Lower',[-Inf -Inf],'Upper',[Inf Inf]);
                                FITNOTGOOD = 0;
                                
                                catch
                                    FITNOTGOOD = 1;
                                end
                                
                                FFFc=coeffvalues(FFFchad);
                                
                                EchadFit=FFFc(1); % module
                                Delta5mT = (FFFc(2) - H45)*1000; % enfoncement a 5mT en nm
                                try
                                FFFci = confint(FFFchad);
                                catch
                                    themall
                                end
                                EchadCI = abs(diff(FFFci(:,1)));
                                
                                Fcompfromfitchad = FFFchad(HcompFit);
                                
                                DistTh = min(HcompFit):0.001:FFFc(2);
                                
                                FcompTh = FFFchad(DistTh);
                                
                                ptr20pN = FcompTh>19 & FcompTh<21;
                                
                                Delta20pN = (FFFc(2) - mean(DistTh(ptr20pN)))*1000;
                                
                                H0FitChad = FFFc(2);
                                
                                
                                
                                
%                                 if isnan(Delta5mT)
%                                     figure
%                                     hold on
%                                     plot(HcompFit*1000,FcompFit,'-*b','linewidth',1.1)
%                                     plot(HcompFit*1000,Fcompfromfitchad,'-og','linewidth',1.5)
%                                     
%                                     close
%                                 end
                                
                                if PLOT
                                    figure
                                    hold on
                                    title([name{i} '-' num2str(ii)])
                                    plot(Dup*1000,Fup,'-*b','linewidth',1.1)
                                    plot(Ddown*1000,Fdown,'-*r','linewidth',1.1)
                                    plot(HcompFit*1000,Fcompfromfitchad,'-m','linewidth',1.5)
                                    xlabel('Thickness (nm)')
                                    ylabel('Force (pN)')
                                    
                                    legend({'Compression','Relaxation',['Chadwick fit : E = ' num2str(EchadFit/1000)]})
                                    
                                    
ax = gca;
ax.LineWidth = 1.5;
box on
                                    
                                    fig = gcf;
                                    saveas(fig,[sfffit filesep name{i} '-' num2str(ii) '_Chadwick'],'png')
                                    saveas(fig,[sfffit filesep name{i} '-' num2str(ii) '_Chadwick'],'fig')
                                    close
                                    
                                end
                                
                                
                                R2chad = 1 - sum((FcompFit - Fcompfromfitchad).^2)/sum((FcompFit - mean(FcompFit)).^2);
                                
                                if FITNOTGOOD == 0
                                Echadtmp = EchadFit;
                                else
                                    Echadtmp = -1;
                                end
                                
                                
                                
                                
                                %% calcul stress-strain avec modele chadwick-like
                                % calcul de l'indentation maximum delta
                                Delt = H0el - Hcomp; % indentation en µm
                                
                                % calcul d'epsilon et sigma
                                EpsEq = Delt/(2*H0el); % delta equivalent = delta / 2 -> volume equivalent compresser par un cylindre
                                Sig = Fcomp./(pi*Rb*Delt); % stress en Pa
                                             
                                %% calcul de l'elasticité a petite deformation
                                
                                E0tmp = [];
                                
                                % on fite sur 100Pa, puis 150, puis 200, etc...
                                
                                step = 50;
                                R2tmp = 0;
                                old_pos = 1;
                                
                                minSig = min(Sig);
                                
                                E0tmp = 0;
                                
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
                                        
                                        EpsFromFit = EpsFit*100;
                                        SigFromFit = Rfit(1) + Rfit(2)*EpsFit;
                                        
                                        
                                    end
                                    
                                end
                                
                                R2 = R2tmp;
                                                               
                                %% %%%%%%% FIGURES %%%%%%
                                if PLOT && 0
                                    figure
                                    subplot(121)
                                    hold on
                                    title('Force-Thickness and stress-strain')
                                    plot(Dup*1000,Fup,'-*b','linewidth',1.1)
                                    plot(Ddown*1000,Fdown,'-*r','linewidth',1.1)
                                    plot((Pel(2)+Pel(1)*[Fh0'].^(1/2))*1000,[Fh0'],'-m','linewidth',3)
                                    plot((Pel(2)+Pel(1)*[0:5:Fh0(1) Fh0'].^(1/2))*1000,[0:5:Fh0(1) Fh0'],'--m','linewidth',2,'handlevisibility','off')
                                    
                                    
                                    plot(Pel(2)*1000,0,'*m','linewidth',5)
                                    %                                 plot(Dup*1000,Fup,'-*b','linewidth',1.1,'handlevisibility','off')
                                    %                                 plot(Ddown*1000,Fdown,'-*r','linewidth',1.1,'handlevisibility','off')
                                    ylabel('Force (pN)')
                                    xlabel('Thickness (nm)','FontName','Serif')
                                    
%                                     xlim([220 500])
                                    %                                 if exist('FcompFromFit','var')
                                    %                                 plot(HcompFromFit,FcompFromFit,'-g','linewidth',1.5)
                                    %                                 clear HcompFromFit FcompFromFit
                                    %                                 legend('Compression','Relaxation','Elastic fit','Fited H0','Dimitriadis fit','FontName','Serif')
                                    %                                 else
                                    %                                 legend('Compression','Relaxation','Elastic fit','Fited H0','FontName','Serif')
                                    %                                 end
                                    
                                    % legend('Compression','Relaxation','FontName','Serif')
                                    % legend('Compression','Relaxation','Elastic fit','FontName','Serif')
                                    % legend('Compression','Relaxation','Elastic fit','Fited H0','FontName','Serif')
                                    
                                    
                                    subplot(122)
                                    hold on
                                    title([name{i} ' - ' num2str(ii)])
                                    plot(EpsEq*100,Sig,'ok','linewidth',1.5)
                                    if exist('EpsFromFit','var')
                                        plot(EpsFromFit,SigFromFit,'-b','linewidth',2)
                                        clear EpsFromFit SigFromFit
                                        legend('Experimantal data',['Linear fit - E = ' num2str(E0tmp)],'location','best')
                                    end
                                    xlabel('Strain (%)','FontName','Serif')
                                    ylabel('Stress (Pa)','FontName','Serif')
                                    
%                                     xlim([0 25])
%                                     ylim([0 700])
                                    
                                    fig = gcf;
                                    saveas(fig,[sffm filesep name{i} '-' num2str(ii)],'png')
                                    saveas(fig,[sffm filesep name{i} '-' num2str(ii)],'fig')
                                    close
                                    
                                    
                                    
                                    
                                end
                                                              
                                %% calcul de l'exposant non lineaire
                                
                                Lesp = log(EpsEq);
                                Lsig = log(Sig);
                                
                                Pnl = polyfit(Lesp,Lsig,1);
                                
                                Lsigff = Pnl(2) + Pnl(1)*Lesp;
                                
                                R2al = 1 - sum((Lsig - Lsigff).^2)/sum((Lsig - mean(Lsig)).^2);
                                
                                
                                if (R2chad > R2sig)  && (Echadtmp>0) % && (Echadtmp<50000)
                                    
                                    % on prend les valeurs pour plot
                                    
                                    H0abs{i} = [H0abs{i}; AbsThick]; % valeur de l'épaisseur du cortex avant compression
                                    
                                    
                                    
                                    E0{i} = [E0{i}; E0tmp];
                                    Edim{i} = [Edim{i}; Edimtmp];
                                    Echad{i} = [Echad{i}; Echadtmp];
                                    EchadConf{i} = [EchadConf{i}; EchadCI];
                                    MaxDef{i} = [MaxDef{i}; maxdef];
                                    MinThick{i} = [MinThick{i}; minthick];
                                    Delt5{i} = [Delt5{i}; Delta5mT];
                                    Delt20{i} = [Delt20{i}; Delta20pN];
                                    ChadsH0{i} = [ChadsH0{i}; H0FitChad*1000];
                                    DimH0{i} = [DimH0{i}; H0dimFit]; % percentile correspondant
                                    
                                    
                                    CompNum{i} = [CompNum{i}; ii];
                                    
                                    
                                    
                                    AlphaNL{i} = [AlphaNL{i}; Pnl(1)];
                                    
                                    
                                    %% calcul de l'hysteresis
                                    %
                                    %                                 AreaComp = abs(trapz(Dup,Fup));
                                    %                                 AreaRel = abs(trapz(Ddown,Fdown));
                                    %                                 Hys{i} = [Hys{i}; (AreaComp-AreaRel)/AreaComp*100];
                                    %
                                    %                                 nCompOk = nCompOk +1;
                                    %
                                end
                                
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

save([sf filesep 'D2M_' date '_' manip '_' specif '.mat'],'E0','Edim','Echad','EchadConf','MaxDef','MinThick','AlphaNL','Cort','DimH0','H0abs','ChadsH0','H0el','CompNum','Delt5','Delt20')


end


