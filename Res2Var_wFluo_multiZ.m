function Res2Var_wFluo_multiZ(date,tag,manip,specif,fich,resfile,nimg,depthoname,AUTO,...
    datafolder,resfolder)

%       Res2Var is using the imageJ particle analysis data to create
%       individual trajectories for the beads. It is using the txt files
%       created by ImageJ to get the X and Y positions the beads, it then
%       loads the images of the video and computes the DZ between the beads
%       using the depthograph. This particular version (_wFluo_multiZ)
%       is used for analyzing constant field experiments, with triplet 
%       of image in Z for each time points. It can handle the presence of
%       fluo images at the end of each loop.
%       
%       Res2Var_wFluo_multiZ(restype,date,tag,manip,specif,fich,
%       depthoname,nimgbr,nimgr,nimg,...
%       AUTO,datafolder,resfolder)
%
%       date : date of experiment 'dd-mm-yy'
%       tag : tag present on the tif and txt files (ex: '5mT', '5mT_Dicty')
%       manip : 'manip' number (ex: 'M1' ou 'M2' etc...), 
%       specif : experimental condiion (ex: 'Dicty_WT')
%       fich : chambers to analyze for the condition (ex: 1 or [1:3])
%       resfile : extension of text file (ex: 'Results' or  'Inside')
%       nimg : nb img total in a loop, can be put to 3 as default if there
%       is no fluo images in the movies
%       depthoname : depthograph (ex:'19-05-20_DepthographM450')
%       AUTO : automatique mode of analysis (skip file when user input
%       is needed, AUTO variable should be defined in the main file) 
%       datafolder : path to the raw data, should be RawdataFolder in main
%       resfolder : path to where save the matlabdata, should be 
%       MatfileFolder in main

set(0,'DefaultFigureWindowStyle','docked')

warning('off','all')

% dossier de données et sauvegarde
path =[datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2)];
sf   = resfolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'R2V'])

cd(path)

pixelconv=1/15.8;

specifchamp = ['_' tag];


specifload=['_' tag];


savename=['R2V_' date '_' manip specifchamp '_' specif];


savenamepart=[savename '_Part'] ;

if nargin == 8
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

cprintf('Unt',['\n\n' date '_' manip '_' specif '_AUTO=' num2str(AUTO) '\n\n'])

% Recuperation des noms de fichiers
[ListDep, ListFor] = get_fich_name(fich);


E = exist([sf filesep 'R2V' filesep savenamepart '.mat'],'file');

if E == 2
    fprintf('\nChargement de données partielle précédement traitées...')
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
    cprintf('Com', [' OK.']);
    fprintf(['\n' num2str(kii) ' expériences chargées.\n\n\n']);
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
        fprintf('\nChargement de données précédement traitées...')
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
        
        cprintf('Com', [' OK.']);
        fprintf(['\n' num2str(kii) ' expériences chargées.\n\n\n']);
    end
    else
         fprintf('\n[AUTOMODE] Chargement de données précédement traitées...')
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
        cprintf('Com',[' OK.']);
        fprintf(['\n[AUTOMODE] ' num2str(kii) ' expériences chargées.\n\n\n']);
    end
end


nacq=length(ListFor);




for ki=1:nacq
    
    Noim=ListDep{ki}; Nofich=ListFor{ki};
    
    
    stackname = [date '_' manip '_' Nofich specifload '.tif'];
    name=[date '_' manip '_' Noim specifload];
    resname = [name '_' resfile '.txt'];
    
    
    E = exist([path filesep resname],'file');
    
    
    
    if E == 2
        
%         try
            
            
            
            
            fprintf(['Results ' date '_' manip '_' Noim specifload ' : \n\n']);
            
            clear S* X* Y* ptr* D* Area* x* y* p1 p2 N H Mat ...
                B M F dz Txy lxy k1 Ind ;
            
            
            fprintf('Lecture du fichier Results...');
            Mat=dlmread([path filesep resname],'\t',1,0);
            fid = fopen([path filesep resname]);
            HeadAll = fgetl(fid);
            fclose(fid);
            HeadCell = textscan(HeadAll,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
            HeadTable = [HeadCell{:}];
            
            Nxm = find(strcmp(HeadTable,'XM')) +1 ;
            Nym = find(strcmp(HeadTable,'YM')) +1 ;
            Ns = find(strcmp(HeadTable,'Slice')) +1 ;
            Nstd = find(strcmp(HeadTable,'StdDev')) +1 ;
            
            
            Nx = find(strcmp(HeadTable,'X')) +1 ;
            Ny = find(strcmp(HeadTable,'Y')) +1 ;
            Nmaj = find(strcmp(HeadTable,'Major')) +1 ;
            Nmin = find(strcmp(HeadTable,'Minor')) +1 ;
            
            
            cprintf('Com',' OK\n');
            
            
            
              fprintf('Verification de la vidéo...');
            
            S = Mat(:,Ns);
            
            % on sépare les point selon le Z de mesure
            Sfull = 1:max(S);
            
          
            [Sdown,Smid,Sup] = SplitZimage(Sfull,nimg);
            
            [Smd,ptrd,ptrm] = intersect(Sdown+1,Smid);    
            [Smdu,ptrmd,ptru] = intersect(Smd,Sup-1);
            
            
            l = min([length(Sdown) length(Smid) length(Sup)]);
            
            Sdown = Sdown(ptrd(ptrmd));
            Smid  = Smdu ;
            Sup   = Sup(ptru)  ;
            
           Sall = [Sdown;Smid;Sup];
            
            stackpath = [path filesep stackname];
            
            for is = 1:max(S)
                try
                Itmp = imread(stackpath,'tiff',is);
                catch
                    lol
                end
                IntTot = sum(sum(Itmp));
                if IntTot == 0
                    [~, col] = find(Sall == is);
                    Sdown(col) = [];
                    Smid(col) = [];
                    Sup(col) = [];                    
                end
                
            end
            
             Sall = [Sdown;Smid;Sup];
            
            cprintf('Com',' OK\n');
            
            
            
            % Séparation des données relative a chaque bille
            fprintf('Tri des données par bille...');
            
            if strcmp(resfile,'Results')
                tag = '';
            else
                tag= resfile;
            end
            distcheck = 16; % max distance (in px) between last know position and next point 
            
            [X1tot,Y1tot,S1tot,ptr1,X2tot,Y2tot,S2tot,ptr2] = splitBeads([Mat(:,Nxm) S],...
                [path filesep stackname],[path filesep name],Mat(:,Nym),distcheck,AUTO,tag);
            
            if sum(S1tot == S2tot) == length(S1tot)
                Stot = S1tot;
                clear S1tot S2tot
            end
            
            % on sépare les point selon le Z de mesure
           
            
            STD1 = nan(3,length(Sdown));
            STD2 = STD1;
            
            
            X1d = zeros(size(Sdown));
            X2d = X1d;
            [~,pd1,p1d] = intersect(Sdown,Stot);
            [~,pd2,p2d] = intersect(Sdown,Stot);
            X1d(pd1) = X1tot(p1d);
            STD1(1,pd1) = Mat(ptr1(p1d),Nstd);
            X2d(pd2) = X2tot(p2d);
            STD2(1,pd2) = Mat(ptr2(p2d),Nstd);
            
            X1m = zeros(size(Smid));
            X2m = X1m;
            [~,pm1,p1m] = intersect(Smid,Stot);
            [~,pm2,p2m] = intersect(Smid,Stot);
            X1m(pm1) = X1tot(p1m);
            STD1(2,pm1) = Mat(ptr1(p1m),Nstd);
            X2m(pm2) = X2tot(p2m);
            STD2(2,pm2) = Mat(ptr2(p2m),Nstd);
            
            X1u = zeros(size(Sup));
            X2u = X1u;
            [~,pu1,p1u] = intersect(Sup,Stot);
            [~,pu2,p2u] = intersect(Sup,Stot);
            X1u(pu1) = X1tot(p1u);
            STD1(3,pu1) = Mat(ptr1(p1u),Nstd);
            X2u(pu2) = X2tot(p2u);
            STD2(3,pu2) = Mat(ptr2(p2u),Nstd);
            
            
            m = ismember(STD1+STD2,max(STD1+STD2));
            
            
            X1all = [X1d; X1m; X1u];
            X2all = [X2d; X2m; X2u];
            
            X1 = X1all(m);
            X2 = X2all(m);
            S = Sall(m);
            
            % redefinition de p1 et p2 après la récuperation des "bon" x
            p1 = ismember(Stot,S);
            p2 = ismember(Stot,S);
            
            Y1 = Y1tot(p1);
            Y2 = Y2tot(p2);
            
            
            
            cprintf('Com',' OK\n');
            
            fprintf('Calcul de dz...\n');
            
%             Z2 = Zcorr_multiZ(ZStep,X2tot,X2,Y2tot,Y2,Stot,S,stackname,Sdown,Smid,Sup,depthoname,datafolder);
%             Z1 = Zcorr_multiZ(ZStep,X1tot,X1,Y1tot,Y1,Stot,S,stackname,Sdown,Smid,Sup,depthoname,datafolder);
%             
%             dz= abs(Z1 - Z2);
            
            
            dz = DzCalc_Zcorr_multiZ(X1tot,X1,Y1tot,Y1,X2tot,...
                X2,Y2tot,Y2,Stot,S,stackname,Sdown,Smid,Sup,depthoname,datafolder);
            
            dz = smooth(dz,'rlowess');
            

            dz = dz*1.33/1.52; % correction pour la difference d'indice optique huile/eau
            
            cprintf('Com', ' ...OK\n');
            
            kii = kii +1; % compteur de ficheir analysé
            
            MT{kii}.exp = name;
            MT{kii}.stack = stackname;
            MT{kii}.S = S;
            MT{kii}.x = [X1 X2];
            MT{kii}.y = [Y1 Y2];
            MT{kii}.dz = dz;
            MT{kii}.std = [STD1 STD2];
            
            ListD{kii} = Noim;
            ListF{kii} = Nofich;
            
            fprintf('\nSauvegarde partielle des positions...');
            cd([sf filesep 'R2V'])
            save(savenamepart,'MT','ListD','ListF');
            cd(path)
            cprintf('Com', ' OK\n\n')
            
%         catch ERR
%             cprintf('*err',['\n\nFile ' name ' has not been analyzed. \n'])
%             cprintf('-err',['Cause : ' ERR.message '\n'])
%             cprintf('err',['Error in ' ERR.stack(1).name '\n\n'])
%         end
        
    end
end


if exist('MT')
    if AUTO
        fprintf('\nSauvegarde partielle des positions...(END OF AUTOMODE)\n\n');
        cd([sf filesep 'R2V'])
        save(savenamepart,'MT','ListD','ListF');
        cd(path)
    else
        fprintf('Sauvegarde des positions...');
        cd([sf filesep 'R2V'])
        save(savename,'MT','ListD','ListF');
        delete([savenamepart '.mat']);
        cd(path)
        
        cprintf('Com', ' OK\n\n\n')
    end
end