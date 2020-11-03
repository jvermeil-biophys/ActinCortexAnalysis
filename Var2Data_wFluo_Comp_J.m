function [dist,force] = Var2Data_wFluo_Comp_J(date,manip,tag,specif,ExperimentalConditions,nimgr,nimgtot,datafolder,resfolder,figurefolder)

warning('off')


set(0,'DefaultFigureWindowStyle','docked')


% save folder
datenow = datestr(now,'yy-mm-dd');
warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',16)


% dossier de données et sauvegarde
fieldpath = [datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2)];
path = [resfolder filesep 'R2V'];
subFigureFolder = [figurefolder filesep datenow filesep 'Meca' filesep specif];
subFigureFolderBF = [subFigureFolder filesep 'B&F'];
if ~exist(subFigureFolder,'dir')
    mkdir(subFigureFolder);
end
if ~exist(subFigureFolderBF,'dir')
    mkdir(subFigureFolderBF);
end


datenow = datestr(now,'yy-mm-dd');

mkdir([resfolder filesep 'V2D'])
mkdir([subFigureFolder])

cd(path)

specifchamp = ['_' tag];

% constantes
DIAMETER = ExperimentalConditions.DiametreBilles; % Diametre de la bille en nm
RADIUS = DIAMETER/2; % Rayon de la bille en nm
SCALE=ExperimentalConditions.Scale;
MAGFIELDCORRECTIONCOEFF = ExperimentalConditions.MagFieldCorrectionCoeff;


specifload=['_' tag 'mT'];
savename=['V2D_' date '_' manip specifchamp '_' specif];
savenamepart=[savename '_Part'] ;
loadname=['R2V_' date '_' manip specifchamp '_' specif '.mat'];


AllRamp = [];
AllForce = [];
Stds = [];

if exist(loadname,'file')
    
    fprintf(['Chargement du fichier ' loadname '\n\n'])
    load(loadname);
    nc = length(MT);
    
    for kc = 1:nc
        
        if not(isempty(MT{kc}))
            VALID = 1;
            name = MT{kc}.exp;
            if strcmp(name(end-1:end),'in')
                fieldname = [name(1:end-2) '_Field.txt'];
            else
                fieldname = [name '_Field.txt'];
            end
            
            fprintf(['Cellule : ' name '\n\n'])
            
            
            
            % donnée partie constante
            fprintf('Récuperation des données et calcul pour Cst... ');
            
            Scst = MT{kc}.Scst;
            
            X1cst = MT{kc}.xcst(:,1);
            X2cst = MT{kc}.xcst(:,2);
            
            Y1cst = MT{kc}.ycst(:,1);
            Y2cst = MT{kc}.ycst(:,2);
            
           
            dzcst = MT{kc}.dzcst/1000;
            
            % correction des points problématique
            
            
            
%             % trop grand saut en DZ
%             ddzcst = diff(dzcst);
%             while max(ddzcst*1000)> 700
%                 ptrdel = find(ddzcst*1000>700)+1;
%                 X1cst(ptrdel) = [];
%                 X2cst(ptrdel) = [];
%                 dzcst(ptrdel) = [];
%                 Y1cst(ptrdel) = [];
%                 Y2cst(ptrdel) = [];
%                 Scst(ptrdel) = [];
%                 ddzcst=diff(dzcst);
%             end
            
            

            
            
            
            % création des variable xp yp
            xpcst = [X1cst,X2cst];
            ypcst = [Y1cst,Y2cst];
            
            % application du facteur d'échelle microscope
            xcst=xpcst*SCALE; ycst=ypcst*SCALE;
            
            % calcul distance en x (D),
            % dans le plan (D2) et l'espace (D3)
            Dcst = abs(xcst(:,2)-xcst(:,1));
            dycst = abs(ycst(:,1)-ycst(:,2));
            D2cst =(Dcst.^2+dycst.^2).^0.5;
            D3cst = (D2cst.^2+dzcst.^2).^0.5;
            
            SrtD3 = sort(D3cst*1000-DIAMETER);
%             Npts = length(SrtD3);
%             P = (Npts:-1:1)/Npts;
%             ptrcort = P>0.9;
%             
%             
%             CortBot = max(SrtD3(ptrcort));
            
            CortBot = prctile(SrtD3,10);
            
            VALID = CortBot < 5000; % 500 auparavant
            
            fprintf(' OK\n');
            
            if VALID
                %% donnée partie compression
                fprintf('Récuperation des données et calcul pour Rmp');
                
                Srmp = MT{kc}.Srmp;

                X1rmp = MT{kc}.xrmp(:,1);
                X2rmp = MT{kc}.xrmp(:,2);
                
                Y1rmp = MT{kc}.yrmp(:,1);
                Y2rmp = MT{kc}.yrmp(:,2);
                
                dzrmp = MT{kc}.dzrmp/1000;
                
                % correction des points problématique
                
                % trop grand saut en DZ
                ddzrmp = diff(dzrmp);
                while max(ddzrmp)> 700
                    ptrdel = find(ddzrmp>700)+1;
                    X1rmp(ptrdel) = [];
                    X2rmp(ptrdel) = [];
                    dzrmp(ptrdel) = [];
                    Y1rmp(ptrdel) = [];
                    Y2rmp(ptrdel) = [];
                    Srmp(ptrdel) = [];
                    ddzrmp=diff(dzrmp);
                end
                

                            
                ncycle = round(Scst(end)/nimgtot);
                
                
                
                 % création des variable xp yp
                xprmp = [X1rmp,X2rmp];
                yprmp = [Y1rmp,Y2rmp];
                
                % application du facteur d'échelle microscope
                xrmp=xprmp*SCALE; yrmp=yprmp*SCALE;
                
                % calcul distance en x (D),
                % dans le plan (D2) et l'espace (D3)
                Drmp = abs(xrmp(:,2)-xrmp(:,1));
                dyrmp = abs(yrmp(:,1)-yrmp(:,2));
                D2rmp =(Drmp.^2+dyrmp.^2).^0.5;
                D3rmp = (D2rmp.^2+dzrmp.^2).^0.5;
                
%                 % CONDITION SUR MAX(dz) POUR FAIRE QUE D2 = D3 !!
%                 if MODE == "2D"
%                     D3rmp = D2rmp;
%                 end
                
                fprintf(' OK\n');
                
                              
                %% Courbe temporelle des distance et force
                fprintf('Reconstitution de la courbe temporelle...');
                
                
                D3all = [D3cst; D3rmp];
                D2all = [D2cst; D2rmp];
                dzall = [dzcst; dzrmp];
                Dall = [Dcst; Drmp];
                Sall = [Scst; Srmp];
                D3ind = [Sall D3all];
                D2ind = [Sall D2all];
                Dind = [Sall Dall];
                dzind = [Sall dzall];
                D3sort = sortrows(D3ind);
                D2sort = sortrows(D2ind);
                Dsort = sortrows(Dind);
                dzsort = sortrows(dzind);
                
                D3 = D3sort(:,2);
                D2 = D2sort(:,2);
                D = Dsort(:,2);
                dz = dzsort(:,2);
                S = D3sort(:,1);
                
                ptrcst = ismember(S,Scst);
                ptrrmp = ismember(S,Srmp);
                
                dzrmp = dz(ptrrmp);
                D2rmp = D2(ptrrmp);
                D3rmp = D3(ptrrmp);
                dzcst = dz(ptrcst);
                D2cst = D2(ptrcst);
                D3cst = D3(ptrcst);
     
                
                fprintf(' OK\n');
                
                
                fprintf('Création des vecteurs champ et temps...');
                
                
                BTMat = dlmread([fieldpath filesep fieldname],'\t',0,0);
                
                B = BTMat(S,1)*MAGFIELDCORRECTIONCOEFF;
                Bcst = BTMat(Scst,1)*MAGFIELDCORRECTIONCOEFF;
                Brmp = BTMat(Srmp,1)*MAGFIELDCORRECTIONCOEFF;
                
                Tfull = (BTMat(:,2)-BTMat(1,2))/1000;
                
                Tfull2 = [Tfull(2:end); Tfull(end)+0.04];
                
                kcycle = [1:ncycle] * nimgtot;                
                
                Tfull2(kcycle) = Tfull2(kcycle-1)+0.04;                    
                
                % Tfull c'est toutes les positions temporelles de toutes
                % les images prises. Mais en fait un certain nombre de ces
                % images appartiennent à un triplet qu'il faut repérer par
                % 1 seule position dans le temps. C'est ce qui est fait
                % dans la liste T, grâce à la liste S qui contient les
                % numéro des frames "centrales" des triplets, et des frames
                % de rampe.
                
                T = Tfull2(S);
                Tcst = (BTMat(Scst,2)-BTMat(1,2))/1000;
                Trmp = (BTMat(Srmp,2)-BTMat(1,2))/1000;
                
                fprintf(' OK\n');
                
  
                
                fprintf('Création des vecteurs force...');
                
                
                
                % calcul de la force dipolaire
                if ExperimentalConditions.SubstrateCoating == 'Float'
                    Dvalid = D2;
                    Dvalid_nm = D2*1000;
                else
                    Dvalid = D3;
                    Dvalid_nm = D3*1000;
                end
                
                M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1
                
                V = 4/3*pi*(DIAMETER/2)^3; % en nm^3
                
                m = M.*10^-9*V; % moment dipolaire en A.nm^2
                
                Bind = 2*10^5*m./(Dvalid_nm.^3); % champ induit par la magnétisation d'une bille voisine
                
                Nvois = 1; % nombre de billes au contact
                
                Btot = B + Nvois*Bind;% B corrigé // B induit de la bille voisine
                
                Mtot = 1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
                
                mtot = Mtot.*10^-9*V;
                
                anglefactor=abs(3*(D./Dvalid).^2-1); % effet angle de la chaine / champ B
                
                F = 3*10^5.*anglefactor.*mtot.^2./Dvalid_nm.^4; % en pN // sans la force des voisins
                
                %% PLOT B&F
                
                figure(kc)
                title(name)
                ax1 = subplot(4,1,1);
                hold on
                plot(T,B,'b-')
%                 xlabel('Time(s)','FontName','Serif')
                xticklabels(ax1,{})
                ylabel('Magnetic Field (mT)','FontName','Serif')
                
                ax2 = subplot(4,1,2);
                hold on
                plot(T,F,'r-')
%                 xlabel('Time(s)','FontName','Serif')
                xticklabels(ax2,{})
                ylabel('Force (pN)','FontName','Serif')
                
                ax3 = subplot(4,1,3);
                hold on
                plot(T,D2,'m-')
                plot(T,D2,'kx')
%                 xlabel('Time(s)','FontName','Serif')
                xticklabels(ax3,{})
                ylabel('D2','FontName','Serif')
                
                ax4 = subplot(4,1,4);
                hold on
                plot(T,dz,'g-')
                plot(T,dz,'kx')
                xlabel('Time(s)','FontName','Serif')
                ylabel('dz','FontName','Serif')
                
                linkaxes([ax1,ax2,ax3,ax4],'x')
                
                fig = gcf;
                
                
%                 FORMER VERSION OF THE GRAPH

%                 figure(kc)
%                 
%                 hold on
% %                 plot(T,B,'b-')
% %                 plot(T,F,'r-')
% %                 % lsline
% %                 % legend(['P-value for correlation : ' num2str(Pmat(1,2),2)])
% %                 xlabel('T (s)')
% %                 ylabel('B / F')
%                 
%                 yyaxis left
%                 axl = gca;  % Get Current Axes
%                 axl.YColor = [0 0 0];
%                 plot(T,B,'b-')
%                 ylabel('Magnetic Field (mT)','FontName','Serif')
%                 ylim([0 300])
%                 
%                 yyaxis right
%                 axr = gca;
%                 axr.YColor = [0 0 0];
%                 plot(T,F,'r-')
%                 ylabel('Force (pN)','FontName','Serif')
%                 ylim([0 1300])
%                 
%                 title(name)
%                 
%                 xlabel('Time(s)','FontName','Serif')
%                 
%                 fig = gcf;
                
                saveas(fig,[subFigureFolderBF filesep 'B&F&dist_' name '.png'],'png')
                saveas(fig,[subFigureFolderBF filesep 'B&F&dist_' name '.fig'],'fig')
%                 close
                
                %%
                
                Fcst = F(ptrcst);
                Frmp = F(ptrrmp);
                
                fprintf(' OK\n');
                
                
                Tcorr = T(end) - T(1) - 2*(ncycle-1);
                Npts = Tcorr/0.65;
                
                VALID = ((length(D3cst)/Npts)>0.5) && (T(end) - T(1))>60 ; 
                
                %% traitement des compressions individuelle et sauvegarde
                if VALID
                fprintf('Découpage...')
                
                FullRmpDat = [Srmp Trmp D3rmp Frmp Brmp dzrmp D2rmp];
                
                RmpDat = cell(1,ncycle);
                
                
                
                for k = 1:ncycle

                    bgn = (k-1)*nimgtot+(nimgtot-nimgr)/2+1;
                    
                    m = ismember(Srmp,bgn:bgn+nimgr-1);
                    Stmp = Srmp(m);
                    D3tmp = D3rmp(m);
                    Ftmp = Frmp(m);
                    Btmp = Brmp(m);
                    Ttmp = Trmp(m);
                    dztmp = dzrmp(m);
                    D2tmp = D2rmp(m);
                    
                    if not(isempty(D3tmp))
                        
                        dD3 = abs(diff(D3tmp));
                        
                        
                        if ~(length(Stmp)<nimgr)&&max(dD3)*1000 < 100
                            
                            kr = size(AllRamp,1);
                            AllRamp(kr+1,:) = D3tmp-D3tmp(1);
                            AllForce(kr+1,:) = Ftmp;
                            Stds(kr+1) = max(diff(D3tmp))*1000;
                            
                        end
                        
                        
                    [~,Mind] = max(Btmp);
                    bordercst = [Stmp(1)-33:Stmp(1)-1 Stmp(end)+1:Stmp(end)+33];
                    ptrborder = ismember(bordercst,Scst);
                    ptrSborder = ismember(Scst,bordercst);
                    Ncst = sum(ptrborder);
                    
                    NcstNeg = sum(D3cst(ptrSborder)*1000-DIAMETER<0);
                    
                   
                        if (Mind < 11)
                            RmpDat{k} = '2fewpts';
    
                        elseif max(dD3)>100
                            
                            RmpDat{k} = 'BadD3Jump';
                            
                        elseif Ncst < 8
                            
                            RmpDat{k} = 'BadCst';
                            
                        elseif NcstNeg > 11
                            
                            RmpDat{k} = 'NegCst';
                            
                        else
                            RmpDat{k} = [Stmp Ttmp D3tmp Ftmp Btmp dztmp D2tmp];
                            
                        end
                    else
                        RmpDat{k} = 'NoPts';
                    end
                    
                end
                
                              
                
                CstDat = [Scst Tcst D3cst Fcst Bcst dzcst D2cst];
                
                fprintf(' OK\n');
                
                
                
                fprintf('Création des variables...');
                
                try
                    RD = [RmpDat{:}];
                    if ischar(RD)
                        VALID = 0;
                    end
                catch
                    RD = [];
                    clear RD
                end
                
                if VALID
                    MR{kc}.name = name;
                    MR{kc}.D2=D2;
                    MR{kc}.D3=D3;
                    MR{kc}.dz=dz;
                    MR{kc}.B=B;
                    MR{kc}.F=F;
                    MR{kc}.time = T;
                    MR{kc}.S = S;
                    MR{kc}.RampData = {RmpDat, FullRmpDat};
                    MR{kc}.CstData  = CstDat;
                    
                    fprintf(' OK\n\n');
                    
                    
                    
                    
                    fprintf('Sauvegarde partielle...');
                    cd([resfolder filesep 'V2D'])
                    save(savenamepart,'MR');
                    cd(path)
                    
                end
                
                
                
                
                
                
                fprintf(' OK\n\n\n');
                else
                    
                cprintf('err','\nPas assez de pts\n')
                end
            else
                
                cprintf('err','\nCortex trop épais\n')
            end
        end
    end
end

% % figure
% % hold on
% % yyaxis left
% figure(3)
% hold on
% plot(Stds,'*')
% 
% figure(1)
% hold on
% plot(AllRamp(Stds<20,:)','-r')
% % plot(mean(AllRamp,1),'linewidth',4)
% % yyaxis right
% figure(2)
% hold on
% plot(AllForce(Stds<20,:)','k-')
% % plot(mean(AllForce,1))

dist = AllRamp(Stds<20,:)';
force = AllForce(Stds<20,:)';

EE = exist('MR');

if EE == 1
    fprintf('Sauvegarde complete...');
    cd([resfolder filesep 'V2D'])
    save(savename,'MR','AllRamp');
    delete([savenamepart '.mat']);
    cd(path)
    fprintf(' OK\n\n');
    
    fprintf([num2str(kc) ' jeux de données traités\n\n\n\n']);
else
    
    fprintf('No existing file in the given list. \nCheck for a path error. \n\n');
    
end