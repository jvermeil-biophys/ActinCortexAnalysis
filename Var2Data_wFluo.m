function NotSaved = Var2Data_wFluo(date,chmpMag,Bcorrec,manip,specif,Db,...
    datafolder,resfolder,resfolderbis)

warning('off')

NotSaved = {};

set(groot,'DefaultFigureWindowStyle','docked')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot,'DefaultAxesFontSize',30)

load([resfolder filesep 'V2D\Classification.mat'])

% dossier de données et sauvegarde
path = [resfolder filesep 'R2V'];
fieldpath =[datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2)];
sf   = resfolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'V2D'])

if nargin ==9 && not(isempty(resfolderbis))
scf = resfolderbis;

mkdir([scf filesep datenow filesep 'V2D'])
end

specifchamp = ['_' chmpMag];

savename=['V2D_' date '_' manip '_' specif];

savenamepart=[savename '_Part'] ;

loadname=['R2V_' date '_' manip specifchamp '_' specif '.mat'];

% General data
pixelconv=1/15.8;

dxtot = [];
dytot = [];
dztot = [];

if exist([path filesep loadname],'file')
    
    fprintf(['Chargement du fichier ' loadname '\n\n'])
    load([path filesep loadname]);
    nc = length(MT);
    
    for kc = 1:nc
        if not(isempty(MT{kc}))
            Valid = 1;
            name = MT{kc}.exp;
            
            if strcmp(name(end-1:end),'in')
                name = name(1:end-2);
            end
            
            if contains(name,AllMigr)
                class = 'M';
            elseif contains(name,AllEtal)
                class = 'E';
            elseif contains(name,AllGene)
                class = 'G';
            else
                class = 'OK';
            end
            
            if contains(name,AllUne)
                nbBilleInt = 1;
            elseif contains(name,AllDeux)
                nbBilleInt = 2;
            elseif contains(name,AllTrois)
                nbBilleInt = 3;
            elseif contains(name,AllQuatre)
                nbBilleInt = 4;
            elseif contains(name,AllCinq)
                nbBilleInt = 5;
            elseif contains(name,AllSix)
                nbBilleInt = 6;
            else
                nbBilleInt = 0;
%                 error(['No bead number available for ' name])
            end
          
            
            fprintf(['Manip : ' name '\n\n']);
            fprintf(['Classification : ' class '\n\n']);
            fprintf('Récuperation des données sur les billes...');

            fieldname = [name '_Field.txt'];
            
            
            
            S = MT{kc}.S;
            S = unique(S);
            
            FE = exist([fieldpath filesep fieldname],'file');
            
            
           
            X1 = MT{kc}.x(:,1);
            X2 = MT{kc}.x(:,2);
            dxp = X1 - X2; % en pixels
            
       
            
            Y1 = MT{kc}.y(:,1);
            Y2 = MT{kc}.y(:,2);
            dyp = Y1 - Y2; % en pixels
    
    
            
            
            L1 = mean(sqrt((X1(1:end-1)-X1(2:end)).^2+(Y1(1:end-1)-Y1(2:end)).^2));
            
            L2 = mean(sqrt((X2(1:end-1)-X2(2:end)).^2+(Y2(1:end-1)-Y2(2:end)).^2));
            
            L = (L1+L2)/2*pixelconv*1000; % distance moyenne parcourue entree deux images en nm
            
            dz = MT{kc}.dz; % en nanometre
            
            fprintf(' OK\n');
            
            
            fprintf('Création des vecteurs champ et temps...');
            
            if FE
                BTMat = dlmread([fieldpath filesep fieldname],'\t',0,0);
  
                B = Bcorrec*BTMat(S,1);
                T = (BTMat(S,2)-BTMat(1,2))/1000;
            else
                fprintf('\n');
                error(['Missing file : ' fieldname])
                
            end
            
            ImgSep = mean(round(diff(T)/0.1)*0.1);
            
            NptsTot = round((T(end) - T(1))/ImgSep);
%             length(T)
%             plot(T,T,'*')
            
            fprintf(' OK\n');
            
            
            fprintf('Calcul de la distance...');
            
            % conversion en nm des grandeur
            dx = dxp*pixelconv*1000;
            dy = dyp*pixelconv*1000;
            
           
                dxtot = [dxtot dx'];
                dytot = [dytot dy'];
                dztot = [dztot dz'];
                
            % calcul distance
%             [D3, ErrD3] = DCalc(Diam,ErrDiam,dx,ErrX,dy,ErrY,dz,ErrZ);
%             D3tot = D3 + Diam;
            D3tot = sqrt(dx.^2+dy.^2+dz.^2);
            D3 = D3tot - Db;
            Dyz = sqrt(dy.^2+dz.^2);
            D2 = sqrt(dx.^2+dy.^2);
            
            

% figure(200)
% hold on
% plot(mean(dz),mean(dy),'*')

%             figure(65456)
%             hold on
%             minZ = min([mean(Z1) mean(Z2)]);
%             maxZ = max([mean(Z1) mean(Z2)]);            
%             plot3(minZ,maxZ,median(D3),'*')
%             xlabel('Min Z')
%             ylabel('Max Z')
%             zlabel('Med D3')
%             
%             figure(13325)
%             hold on
%             plot(mean([minZ maxZ]),median(D3),'*')
%             xlabel('Mean Z')
%             ylabel('Med D')


            % correction des points problématique
            
            ptrdel2 = find(isnan(dz));
            % mauvais champ mag
            
            if ~contains(specif,'Sep')
            mB = find(B < 0.8*Bcorrec*str2double(chmpMag(1:end-2)));
            ptrdel2 = [ptrdel2 mB'];
            end
             ptrdel2 = unique(ptrdel2);
            

            D3(ptrdel2) = [];
            D3tot(ptrdel2) = [];
            Dyz(ptrdel2) = [];
            T(ptrdel2) = [];
            dz(ptrdel2) = [];
            dy(ptrdel2) = [];
            dx(ptrdel2) = [];
            B(ptrdel2) = [];
            S(ptrdel2) = [];

  

            % trop grand saut en DZ
            ddz = diff(abs(dz)); 
            ddzcount = 0;
            while (max(abs(ddz))> 700)&&(ddzcount<50)
            ptrdel = find(ddz>700)+1;
            ptrdel = [ptrdel; find(ddz<-700)];
            ptrdel = unique(ptrdel);
            D3(ptrdel) = [];
            D3tot(ptrdel) = [];
            Dyz(ptrdel) = [];
            T(ptrdel) = [];
            dz(ptrdel) = [];
            dy(ptrdel) = [];
            dx(ptrdel) = [];
            B(ptrdel) = [];
            S(ptrdel) = [];
            ddz=diff(abs(dz));
            ddzcount = ddzcount +1;
            end
            

            if ddzcount>49
                Valid = 0;           
                cprintf('err', ' Too much ddz problems\n');
                pause(0.5)
            end 
            
            clear ptrdel
            
            % point isolé en 3D (trop haut ou trop bas)
            
          
            dd3 = diff(D3);            
            ptrdeld3 = find(abs(dd3)>300)+1;
            if ~isempty(ptrdeld3)
                ptrdel = [];
               for kp = 1:length(ptrdeld3)
                  if ptrdeld3(kp)>length(dd3) || (abs(dd3(ptrdeld3(kp)))>300 &&...
                          sign(dd3(ptrdeld3(kp)))~=sign(dd3(ptrdeld3(kp)-1)))
                      ptrdel = [ptrdel ptrdeld3(kp)];
                  end
               end
                D3(ptrdel) = [];
            D3tot(ptrdel) = [];
            Dyz(ptrdel) = [];
            T(ptrdel) = [];
            dz(ptrdel) = [];
            dy(ptrdel) = [];
            dx(ptrdel) = [];
            B(ptrdel) = [];
            S(ptrdel) = [];
            end
            
            clear ptrdel
            
            fprintf(' OK\n');
     
        
            fprintf('Création du vecteur force...');
            
            % calcul de la force dipolaire
            
            D3nm = D3+Db; % en nm
            
            M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1

            V = 4/3*pi*(Db/2)^3; % en nm^3
            
            m = M.*10^-9*V; % moment dipolaire en A.nm^2
            try
            Bind = 2*10^5*m./(D3nm.^3); % champ induit par la magnétisation d'une bille voisine
            catch
                continue
            end
            Nvois = 1; % nombre de billes au contact
            
            Btot = B + Nvois*Bind;% B corrigé
            
            Mtot = 1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
            
            mtot = Mtot.*10^-9*V;
            
            anglefactor=abs(3*(dx./D3nm).^2-1); % effet angle de la chaine / champ B
            
            F = 3*10^5.*anglefactor.*mtot.^2./D3nm.^4; % en pN
            
            fprintf(' OK\n');
            
                                            
            fprintf('Validation des données...')
            
            
            
            
            if Valid
            if  (length(T)/NptsTot > 0.5) && ...
                    ((T(end)-T(1) > 120)|| contains(chmpMag,'R80')||...
                    contains(specif,'Sep')||contains(chmpMag,'R40')||...
                    contains(chmpMag,'R90')||...
                    contains(specif,'Inside')||contains(specif,'Outside')) 
                
                Valid = 1;
                fprintf(' OK\n');
            else
                Valid = 0;            
                cprintf('err', ' Not Enough Points\n');
                pause(0.5)
                
            end
            end
            
            if Valid
            fprintf('Création des variables...');
            
            
                MR{kc}.actigood = length(T)/NptsTot > 0.5;
                MR{kc}.name = name;
                MR{kc}.class = class;
                MR{kc}.nbBI = nbBilleInt;
                MR{kc}.stack = MT{kc}.stack;
                MR{kc}.F=F;
                MR{kc}.dx=dx;
                MR{kc}.D2=D2;
                MR{kc}.D3=D3;            
                MR{kc}.dz=dz;
                MR{kc}.dy=dy;
                MR{kc}.Dyz=Dyz;
                MR{kc}.B=B;
                MR{kc}.time = T;
                MR{kc}.S = S;
                MR{kc}.L = L;

                

                
                fprintf(' OK\n\n');
                             
                
                fprintf('Sauvegarde partielle...');
                save([sf filesep 'V2D' filesep savenamepart],'MR');
                fprintf(' OK\n\n\n');
            else
                NotSaved = {NotSaved{:}, name};
            end

            
            clear S X* Y* Z* dx dy dz dxp dyp dzp ptr* D3* B BT* Elps* DistC* ...
                    ErrD T ls S name F Valid ddzcount;
            
        end
        
    end
end

EE = exist('MR');

if EE == 1
    fprintf('Sauvegarde complete...');
    save([sf filesep 'V2D' filesep savename],'MR');
    delete([sf filesep 'V2D' filesep savenamepart '.mat']);
    if nargin ==9 && not(isempty(resfolderbis))
    save([scf filesep datenow filesep 'V2D' filesep savename],'MR');
    end
    fprintf(' OK\n\n');
    
    fprintf([num2str(kc) ' jeux de données traités\n\n\n']);
else
    
    fprintf('No existing file in the given list. \nCheck for a path error. \n\n');
    
end


if isempty(NotSaved)
    NotSaved = {'empty'};
end

NotSaved = NotSaved{:};





