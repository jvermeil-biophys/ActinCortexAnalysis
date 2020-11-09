function Res2Var_wFluo_multiZ_Comp(restype,date,tag,manip,specif,fich,depthoname,nimgbr,nimgr,nimg,...
    AUTO,datafolder,resfolder)

%       Res2Var_wFluo_multiZ_Comp(extension du fichier texte (ex: 'Results' ou 'Inside'), date de la manip au format 'jj-mm-aa', type de manip (ex: '5mT' ou 'R40'),
%       numéro de manip (ex: 'M1' ou 'M2' etc...), nom caractériqtique de la manip (ex: 'Dicty_WT'), puits à analyser (ex: 1 ou [1:3]),
%       nom du depthograph (ex: '19-05-20_DepthographM450'), nombre d'image avant rampe dans une boucle, nombre d'image rampe dans une boucle,
%       nombre d'image total d'une boucle, AUTO, RawdataFolder, MatfileFolder)

set(0,'DefaultFigureWindowStyle','docked')

% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';


% dossier de données et sauvegarde
path =[datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2)];
sf   = resfolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'R2V'])
% 
% scf = 'C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Données';
% 
% mkdir([scf filesep datenow filesep 'R2V'])

cd(path)

pixelconv=1/15.8;

specifchamp = ['_' tag];


specifload=['_' tag];


savename=['R2V_' date '_' manip specifchamp '_' specif];
savenamecst=['R2V_' date '_' manip specifchamp 'cst_' specif ];

savenamepart=[savename '_Part'] ;


if nargin == 10
    button = questdlg('AUTO Mode',...
        'Activate auto mode ?','Yes');
    
    if strcmp(button,'Yes')
        AUTO = 1;
    elseif strcmp(button,'No')
        AUTO = 0;
    else
        Stop
    end
end

% Recuperation des noms de fichiers
[ListDep, ListFor] = get_fich_name(fich);


E = exist([sf filesep 'R2V' filesep savenamepart '.mat'],'file');

if E == 2
    fprintf('\n\nChargement de données partielle précédement traitées...')
    load([sf filesep 'R2V' filesep savenamepart])
    
    c = 1;
    for i = 1:length(ListD)
        if not(isempty(ListD{i}))
            newlistD{c} =  ListD{i};
            c = c +1;
        end
    end
    ListD = newlistD;
    clear c newlistD
    kii = length(MT);
    [ListDep,pDep] = setdiff(ListDep,ListD);
    ListFor = {ListFor{pDep}};
    fprintf([' OK. \n' num2str(kii) ' expériences chargées.\n\n\n']);
else
    kii = 0;
    ListD = {};
    ListF = {};
end

E = exist([sf filesep 'R2V' filesep savename '.mat'],'file');

if E == 2
    
    if not(AUTO)
        %     button = questdlg('Que voulez vous faire du précédent fichier ?','Le fichier existe déjà !',...
        %         'Effacer','Completer','Completer');
        
        button = 'Completer';
        
        if strcmp(button,'Completer')
            fprintf('\n\nChargement de données précédement traitées...')
            load([sf filesep 'R2V' filesep savename])
            
            c = 1;
            for i = 1:length(ListD)
                if not(isempty(ListD{i}))
                    newlistD{c} =  ListD{i};
                    c = c +1;
                end
            end
            ListD = newlistD;
            clear c newlistD
            
            kii = length(MT);
            [ListDep,pDep] = setdiff(ListDep,ListD);
            ListFor = {ListFor{pDep}};
            fprintf([' OK. \n' num2str(kii) ' expériences chargées.\n\n\n']);
        end
    else
        fprintf('\n\n[AUTOMODE] Chargement de données précédement traitées...')
        load([sf filesep 'R2V' filesep savename])
        
        c = 1;
        for i = 1:length(ListD)
            if not(isempty(ListD{i}))
                newlistD{c} =  ListD{i};
                c = c +1;
            end
        end
        ListD = newlistD;
        clear c newlistD
        
        kii = length(MT);
        [ListDep,pDep] = setdiff(ListDep,ListD);
        ListFor = {ListFor{pDep}};
        fprintf([' OK. \n[AUTOMODE] ' num2str(kii) ' expériences chargées.\n\n\n']);
    end
end

nacq=length(ListFor);




for ki=1:nacq
    
    Noim=ListDep{ki}; Nofich=ListFor{ki};
    
    
    stackname = [date '_' manip '_' Nofich specifload '.tif'];
    name=[date '_' manip '_' Noim specifload];
    
    
    resname = [name '_' restype '.txt'];
    
    
    E = exist([path filesep resname],'file');
    
    if E == 2
        
        
        try
            
            kii = kii +1; % compteur de ficheir trouvé
            
            
            
            fprintf([restype ' ' date '_' manip '_' Noim specifload ' : \n\n']);
            
            clear S* X* Y* ptr* D* Area* x* y* p1 p2 N H Mat ...
                B M F dz Txy lxy k1 Ind ;
            
            
            fprintf('Lecture du fichier Results...');
            Mat=dlmread(resname,'\t',1,0);
            fid = fopen(resname);
            HeadAll = fgetl(fid);
            fclose(fid);
            HeadCell = textscan(HeadAll,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
            HeadTable = [HeadCell{:}];
            
            Nxm = find(strcmp(HeadTable,'XM')) +1 ;
            Nym = find(strcmp(HeadTable,'YM')) +1 ;
            Ns = find(strcmp(HeadTable,'Slice')) +1 ;
            
            % On recup STD ou Mean intensity value
            Nstdmean = find(strcmp(HeadTable,'StdDev')) +1 ;
            
            if isempty(Nstdmean)
                
            Nstdmean = find(strcmp(HeadTable,'Mean')) +1 ;
            
            end
            
           
            
            fprintf(' OK\n');
            
            fprintf('Verification de la vidéo...');
            
            S = Mat(:,Ns);
            
            % on sépare les point selon le Z de mesure
            Sfull = 1:max(S);
            
            [Sdown,Smid,Sup,Sramp] = SplitZRampImg(Sfull,nimgbr,nimgr,nimg);
            
            l = min([length(Sdown) length(Smid) length(Sup)]);
            
            Sdown = Sdown(1:l);
            Smid  = Smid(1:l) ;
            Sup   = Sup(1:l)  ;
            
            Sall = [Sdown;Smid;Sup];
            
            stackpath = [path filesep stackname];
            
            cols = [];
            
            for is = 1:max(S)
                Itmp = imread(stackpath,'tiff',is);
                IntTot = sum(sum(Itmp));
                if IntTot == 0
                    [~, col] = find(Sall == is);
                    cols = [cols col];               
                end
                
            end
            
            Sdown(cols) = [];
                    Smid(cols) = [];
                    Sup(cols) = []; 
            

                        Sall = [Sdown;Smid;Sup];
            
            
            
            fprintf(' OK\n');
            
            
            % Séparation des données relative a chaque bille
            fprintf('Tri des données par bille...');
            
            if strcmp(restype,'Results')
                restype2 = '';
            else 
                restype2 = restype;
            end
            
                      distcheck = 16; % max distance (in px) between last know position and next point 
            

            
            [X1tot,Y1tot,S1tot,ptr1,X2tot,Y2tot,S2tot,ptr2] = splitBeads([Mat(:,Nxm) S],...
                [path filesep stackname],[path filesep name],Mat(:,Nym),distcheck,AUTO,restype2);
          
            if sum(S1tot == S2tot) == length(S1tot)
                Stot = S1tot;
                clear S1tot S2tot
            end
            
            
            % pour les points en multiZ on récupère X et STD/Mean selon down/mid/up
            STD1 = nan(3,length(Sdown));
            STD2 = STD1;
            
            % down
            X1d = zeros(size(Sdown));
            X2d = X1d;
            [~,pd1,p1d] = intersect(Sdown,Stot);
            [~,pd2,p2d] = intersect(Sdown,Stot);
            X1d(pd1) = X1tot(p1d);
            STD1(1,pd1) = Mat(ptr1(p1d),Nstdmean);
            X2d(pd2) = X2tot(p2d);
            STD2(1,pd2) = Mat(ptr2(p2d),Nstdmean);
            
            % mid
            X1m = zeros(size(Smid));
            X2m = X1m;
            [~,pm1,p1m] = intersect(Smid,Stot);
            [~,pm2,p2m] = intersect(Smid,Stot);
            X1m(pm1) = X1tot(p1m);
            STD1(2,pm1) = Mat(ptr1(p1m),Nstdmean);
            X2m(pm2) = X2tot(p2m);
            STD2(2,pm2) = Mat(ptr2(p2m),Nstdmean);
            
            % up
            X1u = zeros(size(Sup));
            X2u = X1u;
            [~,pu1,p1u] = intersect(Sup,Stot);
            [~,pu2,p2u] = intersect(Sup,Stot);
            X1u(pu1) = X1tot(p1u);
            STD1(3,pu1) = Mat(ptr1(p1u),Nstdmean);
            X2u(pu2) = X2tot(p2u);
            STD2(3,pu2) = Mat(ptr2(p2u),Nstdmean);
            
            
            m = ismember(STD1+STD2,max(STD1+STD2));
            
            
            X1all = [X1d; X1m; X1u];
            X2all = [X2d; X2m; X2u];
            
            X1cst = X1all(m);
            X2cst = X2all(m);
            Scst = Sall(m);
            
            
            
            % redefinition de p1 et p2 après la récuperation des "bon" x
            % pour les partie cst en force
            p1cst = ismember(Stot,Scst);
            p2cst = ismember(Stot,Scst);
            
            Y1cst = Y1tot(p1cst);
            Y2cst = Y2tot(p2cst);
            
            
            % idem pour les partie avec rampe
            p1rmp = ismember(Stot,Sramp);
            p2rmp = ismember(Stot,Sramp);
            
            X1rmp = X1tot(p1rmp);
            X2rmp = X2tot(p2rmp);
            
            Y1rmp = Y1tot(p1rmp);
            Y2rmp = Y2tot(p2rmp);
            
            STD1rmp = Mat(ptr1(p1rmp),Nstdmean);
            STD2rmp = Mat(ptr2(p2rmp),Nstdmean);
            
            Srmp = Stot(p1rmp);
            
            
            fprintf(' OK\n');
            
            fprintf('Calcul des z et de dz...\n');
            
            
            dzcst = DzCalc_Zcorr_multiZ(X1tot,X1cst,Y1tot,Y1cst,X2tot,...
                X2cst,Y2tot,Y2cst,Stot,Scst,stackname,Sdown,Smid,Sup,depthoname,datafolder);
            
            dzrmp = DzCalc_Zcorr(X1rmp,Y1rmp,X2rmp,Y2rmp,Srmp,stackname,depthoname,datafolder);
            
         
            
           
            
            dzmat = [Scst dzcst';Srmp dzrmp'];
            
            dzmatsort = sortrows(dzmat);
            
            dztot = smooth(dzmatsort(:,2),'rlowess');
            stot = dzmatsort(:,1);
            
            dzrmp = dztot(ismember(stot,Srmp));
            dzcst = dztot(ismember(stot,Scst));
            
            
             dzcst = dzcst*1.33/1.52; % correction pour la difference d'indice optique huile/eau
            
            
            dzrmp = dzrmp*1.33/1.52; % correction pour la difference d'indice optique huile/eau
            
            
            
            fprintf('Z et dz OK\n');
            
            MT{kii}.exp = name;
            MT{kii}.stack = stackname;
            MT{kii}.Sramp = Sramp;
            
            MT{kii}.Scst = Scst;
            MT{kii}.xcst = [X1cst X2cst];
            MT{kii}.ycst = [Y1cst Y2cst];
            MT{kii}.dzcst = dzcst;
            
            MT{kii}.Srmp = Srmp;
            MT{kii}.xrmp = [X1rmp X2rmp];
            MT{kii}.yrmp = [Y1rmp Y2rmp];
            MT{kii}.dzrmp = dzrmp;
            
            ListD{kii} = Noim;
            ListF{kii} = Nofich;
            
            fprintf('\nSauvegarde partielle des positions...');
            cd([sf filesep 'R2V'])
            save(savenamepart,'MT','ListD','ListF');
            cd(path)
            fprintf(' OK\n\n');
            
                    
            
            
        catch ERR
            cprintf('*err',['\n\nFile ' name ' has not been analyzed. \n'])
            cprintf('-err',['Cause : ' ERR.message '\n'])
            cprintf('err',['Error in ' ERR.stack(1).name '\n\n'])
        end
        
        
    end
end


if exist('MT')
    fprintf('Sauvegarde des positions...');
    cd([sf filesep 'R2V'])
    save(savename,'MT','ListD','ListF');
    delete([savenamepart '.mat']);
%     cd([scf filesep datenow filesep 'R2V'])
%     save(savename,'MT','ListD','ListF');
    cd(path)
    fprintf(' OK\n\n\n')
end



%% sauvergade des portions constante pour ananlyse classique
if exist('MT')
    
    for kii = 1:length(MT)
if not(isempty(MT{kii}))
    
            MT2{kii}.exp = MT{kii}.exp;
            MT2{kii}.stack = MT{kii}.stack;
            MT2{kii}.S = MT{kii}.Scst;
            MT2{kii}.x = MT{kii}.xcst;
            MT2{kii}.y = MT{kii}.ycst;
            MT2{kii}.dz = MT{kii}.dzcst;
end
    end
    
    MT = MT2;
    fprintf('Sauvegarde des positions pour la partie cst...');
    cd([sf filesep 'R2V'])
    save(savenamecst,'MT','ListD','ListF');
%     cd([scf filesep datenow filesep 'R2V'])
%     save(savenamecst,'MT','ListD','ListF');
    cd(path)
    fprintf(' OK\n\n\n')
end

% cd(h)
