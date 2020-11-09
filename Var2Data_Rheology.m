function [dist,force] = Var2Data_Rheology(date,tags,Bcorrec,Bcst,manip,specif,Db,...
    datafolder,resfolder,resfolderbis,PLOT,figfolder)

warning('off')


set(0,'DefaultFigureWindowStyle','docked')

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultAxesFontSize',15)


% dossier de données et sauvegarde
path = [resfolder filesep 'R2V'];
fieldpath =[datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2)];
sf   = resfolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'V2D'])

if not(isempty(resfolderbis))
    scf = resfolderbis;
    
    mkdir([scf filesep datenow filesep 'V2D'])
end

cd(path)

specifchamp = ['_' [tags{:}]];

pixelconv=1/15.8;


savename=['V2D_' date '_' manip specifchamp '_' specif];

savenamepart=[savename '_Part'] ;

loadname=['R2V_' date '_' manip specifchamp '_' specif '.mat'];


if exist(loadname,'file')
    
    fprintf(['Chargement du fichier ' loadname '\n\n'])
    load(loadname);
    nc = length(MT);
    
    for kc = 1:nc
        
        if not(isempty(MT{kc}))
            name = MT{kc}.stack(1:end-4);
            if strcmp(name(end-1:end),'in')
                fieldname = [name(1:end-2) '_Field.txt'];
                timename = [name(1:end-2) '_Time.txt'];
            else
                fieldname = [name '_Field.txt'];
                timename = [name '_Time.txt'];
            end
            
            cellname = MT{kc}.cell;
            exp = MT{kc}.tag;
            
            fprintf(['Cellule : ' name '\nExp : ' exp '\n\n'])
            
            
            
            % donnée partie constante
            fprintf('Récuperation des données et calcul de la distance... ');
            
            S = MT{kc}.S;
            
            X1 = MT{kc}.x(:,1);
            X2 = MT{kc}.x(:,2);
            
            Y1 = MT{kc}.y(:,1);
            Y2 = MT{kc}.y(:,2);
            
            
            dz = MT{kc}.dz'/1000;
            
            %             % correction des points problématique
            %
            %
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
            xp = [X1,X2];
            yp = [Y1,Y2];
            
            % application du facteur d'échelle microscope
            x=xp*pixelconv; y=yp*pixelconv;
            
            % calcul distance en x (D),
            % dans le plan (D2) et l'espace (D3)
            D = abs(x(:,2)-x(:,1));
            dy = abs(y(:,1)-y(:,2));
            D2 =(D.^2+dy.^2).^0.5;
            D3 = (D2.^2+dz.^2).^0.5*1000; % en nm
            
            Dist = D3 - Db; % Distance entre les billes
            
            fprintf(' OK\n');
            
            fprintf('Création des vecteurs champ et temps...');
            
            
            BTMat = dlmread([fieldpath filesep fieldname],'\t',0,0);
            
            B = BTMat(S,1)*Bcorrec+Bcst;
            
            Tcam = (BTMat(S,4)-BTMat(1,4))/1000;
            
            
            TimgMat = dlmread([fieldpath filesep timename],'\t',0,0);
                      
            Timg = (TimgMat(S,1)-TimgMat(1,1));
            
            
            
            
            fprintf(' OK\n');
            
            
            
            fprintf('Création des vecteurs force...');
            
            
            
            % calcul de la force dipolaire
                        
            M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1
            
            V = 4/3*pi*(Db/2)^3; % en nm^3
            
            m = M.*10^-9*V; % moment dipolaire en A.nm^2
            
            Bind = 2*10^5*m./(D3.^3); % champ induit par la magnétisation d'une bille voisine
            
            Nvois = 1; % nombre de billes au contact
            
            Btot = B + Nvois*Bind;% B corrigé // B induit de la bille voisine
            
            Mtot = 1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
            
            mtot = Mtot.*10^-9*V;
            
            anglefactor=abs(3*(D./D3).^2-1); % effet angle de la chaine / champ B
            
            F = 3*10^5.*anglefactor.*mtot.^2./D3.^4; % en pN // sans la force des voisins
            
            %                 figure
            %                 hold on
            %                 plot(B)
            %                 plot(F)
            
            
            
            
            
            fprintf('Création des variables...');
            
            
            MR{kc}.name = name;
            MR{kc}.cellname = cellname;
            MR{kc}.exp = exp;
            MR{kc}.Dist = Dist;
            MR{kc}.D3=D3;
            MR{kc}.dz=dz;
            MR{kc}.B=B;
            MR{kc}.F=F;
            MR{kc}.time = Tcam;
            MR{kc}.timeimg = Timg;
            MR{kc}.S = S;
            
            fprintf(' OK\n\n');
            
            
            %%%%%%%%%%%%%%% Temp Plots %%%%%%%%%%%%%%%%%
            if PLOT
           
            figsavefol = [figfolder filesep datenow];
            mkdir(figsavefol);
            
            figure
            hold on
            
            title(['Distance and Force in time for : ' name])
            
            yyaxis left
            plot(Tcam,Dist)
            ylabel('Distance (nm)')
            
            yyaxis right
            plot(Tcam,F)
            ylabel('Force (pN)')
            
            xlabel('Time (s)')
            
            fig = gcf;
            saveas(fig,[figsavefol filesep name '_Time'],'png')
            saveas(fig,[figsavefol filesep name '_Time'],'fig')
            close
            
            figure
            hold on
            
            title(['Force vs. Distance for : ' name])
            
            plot(Dist,F,'-o')
            xlabel('Distance (nm)')
            ylabel('Force (pN)')
            
            
            fig = gcf;
            saveas(fig,[figsavefol filesep name '_FvD'],'png')
            saveas(fig,[figsavefol filesep name '_FvD'],'fig')
            close
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('Sauvegarde partielle...');
            cd([sf filesep 'V2D'])
            save(savenamepart,'MR');
            cd(path)      
            
            fprintf(' OK\n\n\n');
            
        end
    end
end



EE = exist('MR');

if EE == 1
    fprintf('Sauvegarde complete...');
    cd([sf filesep 'V2D'])
    save(savename,'MR');
    delete([savenamepart '.mat']);
    
    if not(isempty(resfolderbis))
        cd([scf filesep datenow filesep 'V2D'])
        save(savename,'MR','AllRamp');
    end
    cd(path)
    fprintf(' OK\n\n');
    
    fprintf([num2str(kc) ' jeux de données traités\n\n\n\n']);
else
    
    fprintf('No existing file in the given list. \nCheck for a path error. \n\n');
    
end

end
