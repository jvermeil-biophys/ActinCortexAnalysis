function [dist,force] = Var2Data_Comp_R248(date,tag,Bcorrec,manip,specif,nimgbr,nimgar,...
    nimgr1,nimgr2,nimgr3,Db,datafolder,resfolder)

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




datenow = datestr(now,'yy-mm-dd');

mkdir([resfolder filesep 'V2D'])


cd(path)

specifchamp = ['_' tag];

pixelconv=1/15.8; % 100X
% pixelconv=1/9.7; % 63X

specifload=['_' tag 'mT'];

savename=['V2D_' date '_' manip specifchamp '_' specif];

savenamepart=[savename '_Part'] ;

loadname=['R2V_' date '_' manip specifchamp '_' specif '.mat'];


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
                D2all = [D2cst; D2rmp]; % AJOUT
                Sall = [Scst; Srmp];
                D3ind = [Sall D3all];
                Dind = [Sall Dall];
                D2ind = [Sall D2all]; 
                dzind = [Sall dzall];
                D3sort = sortrows(D3ind);
                Dsort = sortrows(Dind);
                D2sort = sortrows(D2ind); 
                dzsort = sortrows(dzind);
                
                D3 = D3sort(:,2);
                D = Dsort(:,2);
                D2 = D2sort(:,2); 
                dz = dzsort(:,2);
                S = D3sort(:,1);
                
                ptrcst = ismember(S,Scst);
                ptrrmp = ismember(S,Srmp);
                
                D3rmp = D3(ptrrmp);
                D3cst = D3(ptrcst);
                
                
                fprintf(' OK\n');
                
                
                fprintf('Création des vecteurs champ et temps...');
                
                BTMat = dlmread([fieldpath filesep fieldname],'\t',0,0);
                
                
                
                NimgtotInit = length(BTMat);
                
                Nbt = NimgtotInit;
                
                n1 = 0;
                n2 = 0;
                n3 = 0;
                
                
                while Nbt>0
                    Nbt = Nbt - nimgbr- nimgar- nimgr1;
                    
                    if Nbt < 0
                        
                        error('Incoherence with cycle sizes or incomplete cycle')
                        
                    elseif Nbt == 0
                        n1 = n1 +1 ;
                    elseif Nbt >0
                        n1 = n1 +1 ;
                        
                        Nbt = Nbt - nimgbr- nimgar- nimgr2;
                        
                        
                        if Nbt < 0
                            
                            error('Incoherence with cycle sizes or incomplete cycle')
                        elseif Nbt == 0
                            n2 = n2 +1 ;
                        elseif Nbt >0
                            Nbt = Nbt - nimgbr- nimgar- nimgr3;
                            n2 = n2 +1 ;
                            
                            if Nbt< 0
                                error('Incoherence with cycle sizes or incomplete cycle')
                                
                            else
                                n3 = n3 +1 ;
                            end
                        end
                    end
                end
                
                
                
                
                
                
                
                Tfull = (BTMat(:,2)-BTMat(1,2))/1000;
                
                Tfull2 = [Tfull(2:end); Tfull(end)+0.04];
                
                
                
                krempl = [nimgbr+nimgar+nimgr1];
                
                nimgcycle = [nimgbr+nimgar+nimgr3 nimgbr+nimgar+nimgr1 nimgbr+nimgar+nimgr2];
                
                for ii = 2:n1+n2+n3
                    krempl(ii) = krempl(ii-1) + nimgcycle(mod(ii,3)+1);
                end
                
                
                
                
                
                Tfull2(krempl) = Tfull2(krempl-1)+0.04;
                
                T = Tfull2(S);
                Tcst = Tfull2(Scst);
                Trmp = Tfull2(Srmp);
                
                
                
                
                Bfull = BTMat(:,1)*Bcorrec;
                
                Bfull2 = [Bfull(2:end); Bfull(end)];
                
                Bfull2(krempl) = Bfull2(krempl-1);
                
                B = Bfull2(S)*Bcorrec;
                Bcst = Bfull2(Scst)*Bcorrec;
                Brmp = Bfull2(Srmp)*Bcorrec;
                
                fprintf(' OK\n');
                
                
                
                fprintf('Création des vecteurs force...');
                
                
                
                % calcul de la force dipolaire
                
                D3nm = D3*1000; % en nm
                
                if contains(specif,'M270')
                    
                    M = 0.74257*1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1 pour M270

                    
                else
                    
                    M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1 pour M450
                    
                end
                
                V = 4/3*pi*(Db/2)^3; % en nm^3
                
                m = M.*10^-9*V; % moment dipolaire en A.nm^2
                
                Bind = 2*10^5*m./(D3nm.^3); % champ induit par la magnétisation d'une bille voisine
                
                Nvois = 1; % nombre de billes au contact
                
                Btot = B + Nvois*Bind;% B corrigé
                
                if contains(specif,'M270')
                    
                    Mtot = 0.74257*1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
                    
                else
                    
                    Mtot = 1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
                    
                end
                
                mtot = Mtot.*10^-9*V;
                
                anglefactor=abs(3*(D./D3).^2-1); % effet angle de la chaine / champ B
                
                F = 3*10^5.*anglefactor.*mtot.^2./D3nm.^4; % en pN
                
                
                
                Fcst = F(ptrcst);
                
                Frmp = F(ptrrmp);
                
                
                fprintf(' OK\n');
                
                
                
                %% traitement des compressions individuelle et sauvegarde
                if VALID
                    
                    figure
                    plot(Tcst,D3cst*1000-Db,'*-')
                    title(name)
                    hold on
                    
                    fprintf('Découpage...')
                    
                    FullRmpDat = [Srmp Trmp D3rmp Frmp Brmp dzrmp];
                    
                    
                    kbeg(1:3:length(krempl)) = krempl(1:3:length(krempl))-121;
                    kbeg(2:3:length(krempl)) = krempl(2:3:length(krempl))-216;
                    kbeg(3:3:length(krempl)) = krempl(3:3:length(krempl))-406;
                    
                    kend(1:3:length(krempl)) = krempl(1:3:length(krempl))-24;
                    kend(2:3:length(krempl)) = krempl(2:3:length(krempl))-24;
                    kend(3:3:length(krempl)) = krempl(3:3:length(krempl))-24;
                    
                    for k = 1:n1+n2+n3
                        
                        
                        m = ismember(Srmp,kbeg(k):kend(k));
                        Stmp = Srmp(m);
                        D3tmp = D3rmp(m);
                        Ftmp = Frmp(m);
                        Btmp = Brmp(m);
                        Ttmp = Trmp(m);
                        dztmp = dzrmp(m);
                        
                        
                        
                        
                        if k == 1
                        mcstbefore = ismember(Scst,1:kbeg(1)-1);  
                         mcstsurround = ismember(Scst,[1:kbeg(1)-1 kend(1)+1:krempl(1)]);
                        else                            
                        mcstbefore = ismember(Scst,kend(k-1)+1:kbeg(k)-1);   
                         mcstsurround = ismember(Scst,[krempl(k-1)+1:kbeg(k)-1 kend(k)+1:krempl(k)]);
                        end
                        
                        
                        plot(Ttmp,D3tmp*1000-Db,'*')
                        
                        if not(isempty(D3tmp))
                            
                            dD3 = abs(diff(D3tmp));
                            
                            
                            
                            
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
                                RmpDat{k} = [Stmp Ttmp D3tmp Ftmp Btmp dztmp];
                                
                            end
                            
                            CstRmp{k} = {D3cst(mcstsurround) D3cst(mcstbefore)};
                            
                        else
                            RmpDat{k} = 'NoPts';
                            
                            CstRmp{k} = 'NoPts';
                        end
                        
                      
                            
                        if length(Stmp) <= nimgr1
                            CompDuration{k} = '1s';
                        elseif (length(Stmp) > nimgr1) && (length(Stmp) <= nimgr2)
                            CompDuration{k} = '2s';
                               else
                            CompDuration{k} = '4s';
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
                        MR{kc}.RampData = {RmpDat, FullRmpDat,CompDuration,CstRmp};
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


EE = exist('MR');

if EE == 1
    fprintf('Sauvegarde complete...');
    cd([resfolder filesep 'V2D'])
    save(savename,'MR');
    delete([savenamepart '.mat']);
    cd(path)
    fprintf(' OK\n\n');
    
    fprintf([num2str(kc) ' jeux de données traités\n\n\n\n']);
else
    
    fprintf('No existing file in the given list. \nCheck for a path error. \n\n');
    
end
% cd(h)