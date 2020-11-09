function ResInside(argchamp,pc,kimoname,MultiZ)

set(0,'DefaultFigureWindowStyle','docked')

% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';


% dossier de données et sauvegarde
path =['D:\Data\Raw\BilleInternes\' num2str(pc) 'pc'] ;
sf   = 'D:\Data\MatFile';

datenow = datestr(now,'yy-mm-dd');

sfa   = ['D:\Data\Figures' filesep datenow filesep 'InsideBeads'];

mkdir(sfa)

date = get_dates();
manip = {'M1','M2'};
% champ = {'_5mT','_20mT'};
champ = ['_' num2str(argchamp) 'mT'];


pixelconv = 1/15.8;

Cpt = 0;

ListDep = get_fich_name();

E = exist([sf filesep 'ResInside' champ '_' num2str(pc) 'pc_Part.mat'],'file');

if E == 2
    fprintf('\n\nChargement de données précédement traitées...')
    load([sf filesep 'ResInside' champ '_' num2str(pc) 'pc_Part.mat'])
    Cpt = length(MRI);
    fprintf([' OK. \n' num2str(Cpt) ' expériences chargées.\n\n\n']);
else
    Cpt = 0;
    ListN = {};

end


    for km = 1:2
        for kd = 1:length(date)
            for ki = 1:length(ListDep)
                
                Noim=ListDep{ki};
                
                name=[date{kd} '_' manip{km} '_' Noim champ];
                resname = [name '_Inside.txt'];
                stackname = [name '.tif'];
                
                
                E = exist([path filesep resname],'file');
                
                
                
                if E && not(ismember(name,ListN))
                    Cpt = Cpt + 1;
                    
                    
                    fprintf(['Lecture du fichier ' name '...'])
                    cd(path)
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
                    
                    % Séparation des données relative a chaque bille
                    fprintf('Tri des données par bille...');
                    
                    [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads_Fluo([Mat(:,Nxm) Mat(:,Ns)],[path filesep stackname],Mat(:,Nym),0);
                    
                    if MultiZ
                    
                                    % on sépare les point selon le Z de mesure
            Sfull = 1:max(Mat(:,Ns));
            
            
            [Sdown,Smid,Sup] = SplitZimage(Sfull,3);
            
            
            l = min([length(Sdown) length(Smid) length(Sup)]);
            
            Sdown = Sdown(1:l);
            Smid  = Smid(1:l) ;
            Sup   = Sup(1:l)  ;
            
            Sall = [Sdown;Smid;Sup];
            
            STD1 = nan(3,length(Sdown));
            STD2 = STD1;
            
            
            X1d = zeros(size(Sdown));
            X2d = X1d;
            [~,pd1,p1d] = intersect(Sdown,S1);
            [~,pd2,p2d] = intersect(Sdown,S2);
            X1d(pd1) = X1(p1d);
            STD1(1,pd1) = Mat(ptr1(p1d),Nstd);
            X2d(pd2) = X2(p2d);
            STD2(1,pd2) = Mat(ptr2(p2d),Nstd);
            
            X1m = zeros(size(Smid));
            X2m = X1m;
            [~,pm1,p1m] = intersect(Smid,S1);
            [~,pm2,p2m] = intersect(Smid,S2);
            X1m(pm1) = X1(p1m);
            STD1(2,pm1) = Mat(ptr1(p1m),Nstd);
            X2m(pm2) = X2(p2m);
            STD2(2,pm2) = Mat(ptr2(p2m),Nstd);
            
            X1u = zeros(size(Sup));
            X2u = X1u;
            [~,pu1,p1u] = intersect(Sup,S1);
            [~,pu2,p2u] = intersect(Sup,S2);
            X1u(pu1) = X1(p1u);
            STD1(3,pu1) = Mat(ptr1(p1u),Nstd);
            X2u(pu2) = X2(p2u);
            STD2(3,pu2) = Mat(ptr2(p2u),Nstd);
            
            
            m = ismember(STD1+STD2,max(STD1+STD2));
            
            
            X1all = [X1d; X1m; X1u];
            X2all = [X2d; X2m; X2u];
            
            X1 = X1all(m);
            X2 = X2all(m);
            S = Sall(m);
            STD1 = STD1(m);
            STD2 = STD2(m);
            
            dx = abs(X1-X2);
            % redefinition de p1 et p2 après la récuperation des "bon" x
            p1 = ismember(S1,S);
            p2 = ismember(S2,S);
            
            Y1 = Y1(p1);
            Y2 = Y2(p2);
                        dy = abs(Y1-Y2);
                    else
                        
                    % prise de données sur les images a deux billes
                    [S,p1,p2] = intersect(S1,S2);
                    
                    % récuperation des données par billes
                    
                    X1 = X1(p1);
                    X2 = X2(p2);
                    
                    dx = abs(X1-X2);
                    
                    Y1 = Y1(p1);
                    Y2 = Y2(p2);
                    
                    dy = abs(Y1-Y2);
                    
                    STD1 = Mat(ptr1(p1),Nstd);
                    STD2 = Mat(ptr2(p2),Nstd);
                    
                    end
                    
                    
                    fprintf(' OK\n');
                    
                            fprintf('Calcul des z et de dz...');
                            
                    [Z1 Ind] = Zcorr(STD1,X1,Y1,S,stackname,kimoname);
                    Z1 = abs(Z1);
                    Z2 = abs(Zcorr(STD1,X2,Y2,S,stackname,kimoname));
                    
                    dz= abs(Z1 - Z2);
                         fprintf(' OK\n');
                         
                         
                         fprintf('Calcul des distances et temps...')
                         
                         dx = dx*pixelconv*1000;
                         dy = dy*pixelconv*1000;
                         dz = dz*1000;
                        
                         
                         D = dx;
                         D2 = sqrt(dx.^2+dy.^2);
                         D3 = sqrt(dx.^2+dy.^2+dz.^2);
                         
                         
                         T = S*0.2;
                         
                         fprintf(' OK\n')
                         
                         fprintf('Calcul de l''épaisseur et de l''activité')
                         

                         Cort = mean(D3);
                         
%                          Dm = round(D3./20)*20;              
%                          md = mode(Dm);                                    
%                          Cort = md;
                         
                        

                        tau = 1:100; % un tau = 0.5 seconde 
                        
                        t = ceil(T(1)):0.5:floor(T(end));
                        x = interp1(T,D3,t);
                        
                        
                        H = zeros(1,length(tau));
                        
                        for kt = 1:length(tau)
                            it = 1:length(t)-length(tau); % it représente le temps, chaque pas vaut 0.5s
                            sdx = (x(it) - x(it+tau(kt))).^2;
                            H(kt) = sqrt(sum(sdx)/length(it));
                        end
                                           
%                          Itmp = imread(stackname,'tiff',S(Ind));
%                          
%                          imshow(Itmp)
%                          hold on
%                          plot(X1(Ind),Y1(Ind),'*')
%                          plot(X2(Ind),Y2(Ind),'*')
%                          title(['Cort = ' num2str(Cort-4432) 'nm - Acti = ' num2str(std(D3-4432)) 'nm'])
%                          
%                          fig = gcf;
%                          saveas(fig,[sfa filesep 'Inside_' name '.png'],'png');
%                          close(fig)
                        
                         figure
                         hold on
                         plot(T,D3-4432,'-*')
                         title(['Cort = ' num2str(Cort-4432) 'nm - Acti = ' num2str(std(D3-4432)) 'nm'])
                         
                         fig = gcf;
                         saveas(fig,[sfa filesep 'Inside_' name '_Curve.png'],'png');
                         close(fig)
                        
                         fprintf(' OK\n')
                         
                         
                         fprintf('Création des variables')
                         
                         MRI{Cpt}.Dist = D3-4432;
                         MRI{Cpt}.Time = T;
                         MRI{Cpt}.Cort = Cort;
                         MRI{Cpt}.Tau  = tau/2;
                         MRI{Cpt}.Acti = H;
                         MRI{Cpt}.Dz = dz;
                         MRI{Cpt}.stackname = stackname;
                         
                         fprintf(' OK\n\n')
                         
                                     fprintf('\nSauvegarde partielle des positions...');
            ListN{Cpt} = name;

    save([sf filesep 'ResInside' champ '_' num2str(pc) 'pc_Part.mat'],'MRI','ListN')

            fprintf(' OK\n\n');
                         
                    cd(h)
                end
            end
        end
    end
    
    
    save([sf filesep 'ResInside' champ '_' num2str(pc) 'pc.mat'],'MRI','ListN')
    
end
