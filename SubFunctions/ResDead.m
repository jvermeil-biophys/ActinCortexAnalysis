function ResDead(argchamp)

set(0,'DefaultFigureWindowStyle','docked')

% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';


% dossier de données et sauvegarde
path ='D:\Data\Raw\Dead';
sf   = 'D:\Data\MatFile';

datenow = datestr(now,'yy-mm-dd');

sfa   = ['D:\Data\Figures' filesep datenow filesep 'DeadCells'];

mkdir(sfa)

date = get_dates();
manip = {'M1','M2'};
% champ = {'_5mT','_20mT'};
champ = ['_' num2str(argchamp) 'mT'];


pixelconv = 1/15.8;

Cpt = 0;

[ListDep, ListFor] = get_fich_name();
    for km = 1:2
        for kd = 1:length(date)
            for ki = 1:length(ListFor)
                
                Noim=ListDep{ki};
                
                name=[date{kd} '_' manip{km} '_' Noim champ];
                resname = [name '_Dead.txt'];
                stackname = [name '.tif'];
                
                
                E = exist([path filesep resname],'file');
                
                
                if E
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
                    
                    [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads_Fluo([Mat(:,Nxm) Mat(:,Ns)],[path filesep stackname],Mat(:,Nym));
                    
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
                    
                    
                    if not(isempty(Nmin))
                        MajAx1 = Mat(ptr1(p1),Nmaj);
                        MinAx1 = Mat(ptr1(p1),Nmin);
                        
                        ElpsFctr1 = MajAx1./MinAx1;
                        
                        
                        MajAx2 = Mat(ptr2(p2),Nmaj);
                        MinAx2 = Mat(ptr2(p2),Nmin);
                        
                        ElpsFctr2 = MajAx2./MinAx2;
                    end
                    
                    
                    fprintf(' OK\n');
                    
                            fprintf('Calcul des z et de dz...');
                            
                    [Z1 Ind] = Zcorr(STD1,X1,Y1,S,stackname,'Kimo_Old');
                    Z1 = abs(Z1);
                    Z2 = abs(Zcorr(STD1,X2,Y2,S,stackname,'Kimo_Old'));
                    
                    dz= abs(Z1 - Z2);
                         fprintf(' OK\n');
                         
                         
                         fprintf('Calcul des distances et temps...')
                         
                         dx = dx*pixelconv*1000;
                         dy = dy*pixelconv*1000;
                         dz = dz*1000;
                        
                         
                         D = dx;
                         D2 = sqrt(dx.^2+dy.^2);
                         D3 = sqrt(dx.^2+dy.^2+dz.^2);
                         
                         
                         T = S*0.6;
                         
                         fprintf(' OK\n')
                         
                         fprintf('Calcul de l''épaisseur et de l''activité')
                         

                         Cort = mean(D3);
                         
                        

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
%                          saveas(fig,[sfa filesep 'Dead_' name '.png'],'png');
%                          close(fig)
%                         
%                          figure
%                          hold on
%                          plot(T,D3-4432,'-*')
%                          title(['Cort = ' num2str(Cort-4432) 'nm - Acti = ' num2str(std(D3-4432)) 'nm'])
%                          
%                          fig = gcf;
%                          saveas(fig,[sfa filesep 'Dead_' name '_Curve.png'],'png');
%                          close(fig)
%                         
                         fprintf(' OK\n')
                         
                         
                         fprintf('Création des variables')
                         
                         MRI{Cpt}.Dist = D3-4432;
                         MRI{Cpt}.Time = T;
                         MRI{Cpt}.Cort = Cort;
                         MRI{Cpt}.Tau  = tau/2;
                         MRI{Cpt}.Acti = H;
                         MRI{Cpt}.Elps = [ElpsFctr1 ElpsFctr2];
                         MRI{Cpt}.Dz = dz;
                         MRI{Cpt}.stackname = stackname;
                         
                         fprintf(' OK\n\n')
                         
                    cd(h)
                end
            end
        end
    end
    
    
    save([sf filesep 'ResDead' champ],'MRI')
    
end
