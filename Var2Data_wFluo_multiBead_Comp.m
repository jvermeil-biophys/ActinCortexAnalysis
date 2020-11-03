function [dist,force] = Var2Data_wFluo_multiBead_Comp(date,tag,Bcorrec,manip,specif,nimgr,nimgtot,Db,...
    datafolder,resfolder,resfolderbis)

warning('off')


set(0,'DefaultFigureWindowStyle','docked')

% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';


% % save folder
% sf = 'D:\Data\Figures';
% datenow = datestr(now,'yy-mm-dd');
% sffd = [sf filesep datenow filesep 'Meca' filesep 'ForceDep' filesep specif];

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)


% dossier de données et sauvegarde
path = [resfolder filesep 'R2V'];
fieldpath =[datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2)];
sf   = resfolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'V2D'])

if nargin ==11 && not(isempty(resfolderbis))
scf = resfolderbis;

mkdir([scf filesep datenow filesep 'V2D'])
end

cd(path)

specifchamp = ['_' tag];

pixelconv=1/15.8;

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
            
            
            
            % trop grand saut en DZ
            ddzcst = diff(dzcst);
            while max(ddzcst*1000)> 700
                ptrdel = find(ddzcst*1000>700)+1;
                X1cst(ptrdel) = [];
                X2cst(ptrdel) = [];
                dzcst(ptrdel) = [];
                Y1cst(ptrdel) = [];
                Y2cst(ptrdel) = [];
                Scst(ptrdel) = [];
                ddzcst=diff(dzcst);
            end
            
            

            
            
            
            % création des variable xp yp
            xpcst = [X1cst,X2cst];
            ypcst = [Y1cst,Y2cst];
            
            % application du facteur d'échelle microscope
            xcst=xpcst*pixelconv; ycst=ypcst*pixelconv;
            
            % calcul distance en x (D),
            % dans le plan (D2) et l'espace (D3)
            Dcst = abs(xcst(:,2)-xcst(:,1));
            dycst = abs(ycst(:,1)-ycst(:,2));
            D2cst =(Dcst.^2+dycst.^2).^0.5;
            D3cst = (D2cst.^2+dzcst.^2).^0.5;
            
            SrtD3 = sort(D3cst*1000-Db);
            Npts = length(SrtD3);
            P = (Npts:-1:1)/Npts;
            ptrcort = P>0.9;

            
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
                xrmp=xprmp*pixelconv; yrmp=yprmp*pixelconv;
                
                % calcul distance en x (D),
                % dans le plan (D2) et l'espace (D3)
                Drmp = abs(xrmp(:,2)-xrmp(:,1));
                dyrmp = abs(yrmp(:,1)-yrmp(:,2));
                D2rmp =(Drmp.^2+dyrmp.^2).^0.5;
                D3rmp = (D2rmp.^2+dzrmp.^2).^0.5;
                fprintf(' OK\n');
                
                              
                %% Courbe temporelle des distance et force
                fprintf('Reconstitution de la courbe temporelle...');
                
                
                D3all = [D3cst; D3rmp];
                dzall = [dzcst; dzrmp];
                Dall = [Dcst; Drmp];
                Sall = [Scst; Srmp];
                D3ind = [Sall D3all];
                Dind = [Sall Dall];
                dzind = [Sall dzall];
                D3sort = sortrows(D3ind);
                Dsort = sortrows(Dind);
                dzsort = sortrows(dzind);
                
                D3 = D3sort(:,2);
                D = Dsort(:,2);
                dz = dzsort(:,2);
                S = D3sort(:,1);
                
                ptrcst = ismember(S,Scst);
                ptrrmp = ismember(S,Srmp);
                
                D3rmp = D3(ptrrmp);
                D3cst = D3(ptrcst);
     
                
                fprintf(' OK\n');
                
                
                fprintf('Création des vecteurs champ et temps...');
                
                
                BTMat = dlmread([fieldpath filesep fieldname],'\t',0,0);
                
                B = BTMat(S,1)*Bcorrec;
                Bcst = BTMat(Scst,1)*Bcorrec;
                Brmp = BTMat(Srmp,1)*Bcorrec;
                
                Tfull = (BTMat(:,2)-BTMat(1,2))/1000;
                
                Tfull2 = [Tfull(2:end); Tfull(end)+0.04];
                
                kcycle = [1:ncycle] * nimgtot;                
                
                Tfull2(kcycle) = Tfull2(kcycle-1)+0.04;                    
                
                T = Tfull2(S);
                Tcst = (BTMat(Scst,2)-BTMat(1,2))/1000;
                Trmp = (BTMat(Srmp,2)-BTMat(1,2))/1000;
                
                fprintf(' OK\n');
                
  figure
  hold on
  plot(Tcst,D3cst,'-*')
  plot(Trmp,D3rmp,'*')
  pause
                
                fprintf('Création des vecteurs force...');
                
                
                
                % calcul de la force dipolaire
                
                D3nm = D3*1000; % en nm
                
                M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1
                
                V = 4/3*pi*(Db/2)^3; % en nm^3
                
                m = M.*10^-9*V; % moment dipolaire en A.nm^2
                
                Bind = 2*10^5*m./(D3nm.^3); % champ induit par la magnétisation d'une bille voisine
                
                Nvois = 1; % nombre de billes au contact
                
                Btot = B + Nvois*Bind;% B corrigé // B induit de la bille voisine
                
                Mtot = 1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
                
                mtot = Mtot.*10^-9*V;
                
                anglefactor=abs(3*(D./D3).^2-1); % effet angle de la chaine / champ B
                
                F = 3*10^5.*anglefactor.*mtot.^2./D3nm.^4; % en pN // sans la force des voisins
                
%                 figure
%                 hold on
%                 plot(B)
%                 plot(F)
                
                Fcst = F(ptrcst);
                
                Frmp = F(ptrrmp);
                
                fprintf(' OK\n');
                
                
                Tcorr = T(end) - T(1) - 2*(ncycle-1);
                Npts = Tcorr/0.65;
                
                VALID = ((length(D3cst)/Npts)>0.5) && (T(end) - T(1))>60 ; 
                
                %% traitement des compressions individuelle et sauvegarde
                if VALID
                fprintf('Découpage...')
                
                FullRmpDat = [Srmp Trmp D3rmp Frmp Brmp];
                
                RmpDat = cell(1,ncycle);
                
                
                
                for k = 1:ncycle

                    bgn = (k-1)*nimgtot+(nimgtot-nimgr)/2+1;
                    
                    m = ismember(Srmp,bgn:bgn+nimgr-1);
                    Stmp = Srmp(m);
                    D3tmp = D3rmp(m);
                    Ftmp = Frmp(m);
                    Btmp = Brmp(m);
                    Ttmp = Trmp(m);
                    
                    if not(isempty(D3tmp))
                        
                        dD3 = abs(diff(D3tmp));
                        
                        
                        if ~(length(Stmp)<nimgr)&&max(dD3)*1000<100
                            
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
                    
                    NcstNeg = sum(D3cst(ptrSborder)*1000-Db<0);
                    
                   
                        if (Mind < 11)
                            RmpDat{k} = '2fewpts';
    
                        elseif max(dD3)>100
                            
                            RmpDat{k} = 'BadD3Jump';
                            
                        elseif Ncst < 8
                            
                            RmpDat{k} = 'BadCst';
                            
                        elseif NcstNeg > 11
                            
                            RmpDat{k} = 'NegCst';
                            
                        else
                            RmpDat{k} = [Stmp Ttmp D3tmp Ftmp Btmp];
                            
                        end
                    else
                        RmpDat{k} = 'NoPts';
                    end
                    
                end
                
                              
                
                CstDat = [Scst Tcst D3cst Fcst Bcst];
                
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
                    cd([sf filesep 'V2D'])
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
    cd([sf filesep 'V2D'])
    save(savename,'MR','AllRamp');
    delete([savenamepart '.mat']);
    
    if nargin ==11 && not(isempty(resfolderbis))
    cd([scf filesep datenow filesep 'V2D'])
    save(savename,'MR','AllRamp');
    end
    cd(path)
    fprintf(' OK\n\n');
    
    fprintf([num2str(kc) ' jeux de données traités\n\n\n\n']);
else
    
    fprintf('No existing file in the given list. \nCheck for a path error. \n\n');
    
end
