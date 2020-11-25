function Res2Var_Rheology(restype,date,tags,manip,specif,fich,depthoname,...
    AUTO,datafolder,resfolder)

set(0,'DefaultFigureWindowStyle','docked')


% dossier de données et sauvegarde
path =[datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2)];
sf   = resfolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'R2V'])

cd(path)

pixelconv=1/15.8;

specifchamp = ['_' [tags{:}]];

savename=['R2V_' date '_' manip specifchamp '_' specif];

savenamepart=[savename '_Part'] ;


% Recuperation des noms de fichiers
[ListDep, ListFor] = get_fich_name(fich);

for kt = 1:length(tags)
   for kl = 1:length(ListDep)
    ListDep2{kt,kl} = [ListDep{kl} '_' tags{kt}];
    ListFor2{kt,kl} = [ListFor{kl} '_' tags{kt}];    
   end
end

ListDep = {ListDep2{:}};
ListFor = {ListFor2{:}};

E = exist([sf filesep 'R2V' filesep savenamepart '.mat'],'file');

if E == 2
    fprintf('\n\nChargement de données partielle précédement traitées...')
    load([sf filesep 'R2V' filesep savenamepart])
    
    newlistD = {};
    
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

celllist = '';

for ki=1:nacq
    
    Noim=ListDep{ki}; Nofich=ListFor{ki};
                       
        stackname = [date '_' manip '_' Nofich '.tif'];
        name=[date '_' manip '_' Noim];
        
        stackpath = [path filesep stackname];
        
        resname = [name '_' restype '.txt'];
        
        inds = strfind(name,'_');
        
        cellname = name(1:inds(end)-1);
        
        Tag = name(inds(end)+1:end);
        
        E = exist([path filesep resname],'file');
        
        if E == 2
            
            kii = kii + 1;
            
            
            fprintf([restype ' ' date '_' manip '_' Noim ' : \n\n']);
            
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
            Nstd = find(strcmp(HeadTable,'StdDev')) +1 ;
            
            
            Nx = find(strcmp(HeadTable,'X')) +1 ;
            Ny = find(strcmp(HeadTable,'Y')) +1 ;
            Nmaj = find(strcmp(HeadTable,'Major')) +1 ;
            Nmin = find(strcmp(HeadTable,'Minor')) +1 ;
            
            
            fprintf(' OK\n');
            
            
            fprintf('Verification de la vidéo...');
            
            S = Mat(:,Ns);
            
            
            badInd = [];
            for is = 1:max(S)
                
                
                Itmp = imread(stackpath,'tiff',is);
                IntTot = sum(sum(Itmp));
                if IntTot == 0
                    ind = find(S == is);
                    badInd = [badInd ind];
                end
                
            end
            
            S(badInd) = [];
            
            
            
            fprintf(' OK\n');
            
            fprintf('Tri des données par bille...');
            
            [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads([Mat(:,Nxm) S],...
                [path filesep stackname],[path filesep name],Mat(:,Nym),AUTO);
            
            if sum(S1 == S2) == length(S1)
                Stot = S1;
                clear S1tot S2tot
            end
            
            fprintf(' OK\n');
            
            fprintf('Calcul de dz...\n');
            
            dz = DzCalc_Zcorr_Rheo(X1,Y1,X2,Y2,Stot,stackname,depthoname,datafolder);
            
            dz = dz*1.33/1.52;
            
            fprintf(' OK\n');
            
            MT{kii}.cell = cellname;
            
            MT{kii}.tag = Tag;
            MT{kii}.stack = stackname;
            
            MT{kii}.S = Stot;
            MT{kii}.x = [X1 X2];
            MT{kii}.y = [Y1 Y2];
            MT{kii}.dz = dz;       
            
            ListD{kii} = Noim;
            ListF{kii} = Nofich;
            
            
             fprintf('\nSauvegarde partielle des positions...');
            cd([sf filesep 'R2V'])
            save(savenamepart,'MT','ListD','ListF');
            cd(path)
            fprintf(' OK\n\n');
            
        end
    end
    


if exist('MT')
    fprintf('Sauvegarde des positions...');
    cd([sf filesep 'R2V'])
    save(savename,'MT','ListD','ListF');
        delete([savenamepart '.mat']);
fprintf(' OK\n\n\n')
end


end