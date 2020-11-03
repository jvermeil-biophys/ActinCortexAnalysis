function Data2StressStrain_Rheology(date,manip,specif,tags,Diam,PLOT,resfolder,figfolder)

%% INIT

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultTextFontSize',20)

%% params and settings
% dossier de données et sauvegarde des données
path = [resfolder filesep 'V2D'];
sf   = [resfolder filesep 'D2SS'];

datenow = datestr(now,'yy-mm-dd');

mkdir(sf)

savename=['D2SS_' date '_' [tags{2:end}] '_' specif];
savenamepart = [savename '_Part'];
kk = 1;

if PLOT
    % sauvergarde des figures
    sfig = figfolder;
    sff = [sfig filesep datenow filesep 'Meca' filesep specif filesep 'Fourier'];
    mkdir(sff)
end


%% Retrieve datafile and analyze
for km = 1:length(manip)
    
    loadname=['V2D_' date '_' manip{km} '_' [tags{:}] '_' specif '.mat'];
    
    cd(path)
    
    if exist(loadname,'file')
        
        fprintf(['Chargement du fichier ' loadname '\n\n'])
        load(loadname);
        
        InfoTable = struct('cellname',{''},'Echad',[]);
        
        
        %% Analyse rampe pour avoir un premier module elastique
        fprintf('Calcul de modules elastique a partir des rampes...');
        for kc = 1:length(MR)
            if strcmp(MR{kc}.exp,'RS') % Rampe
                %% Récupération des données
                Dist = MR{kc}.Dist;
                B = MR{kc}.B;
                F = MR{kc}.F;
                CellName = MR{kc}.cellname;
                
                %% Calcul d'un module elastique a partir de la rampe
                
                PtrChamp0 = find(B<0.5,1); % Premier point de champ plus petit que 0.5mT
                
                [~,PtrChampMax] = max(B); % Champ maximum
                
                DistComp = Dist(PtrChamp0:PtrChampMax)/1000; % Evolution de la distance pendant la compression en µm
                ForceComp = F(PtrChamp0:PtrChampMax); % Evolution de la force pendant la compression
                
                FitStart = find(ForceComp>1,1); % On fit a partir du moment ou la force dépasse 1pN
                
                % fit chadwick pour delta/2 et h/2
                
                fitchad = fittype(['(pi*E*R*(x-h0)^2)/(3*h0)']...
                    ,'problem',{'R'}); % fonction a fitter (chadwick), R param, E et h0 a fitter
                
                FFF =fit(DistComp(FitStart:end),ForceComp(FitStart:end),fitchad,'problem',{Diam/2000},'Start',[1000 DistComp(1)]);
                
                FFFc=coeffvalues(FFF);
                
                Echad=FFFc(1); % module
                
                InfoTable.cellname{end+1} = CellName;
                InfoTable.Echad(end+1)   = Echad;
                
                clear EChad CellName
                
            end
        end
        
        fprintf(' OK\n\n\n');
        
        %% Analyse des oscillations
        for kc = 1:length(MR)
            if ~strcmp(MR{kc}.exp,'RS') % Manip avec oscillation
                
                %% Récupération des données
                
                ExpType = MR{kc}.exp;
                Dist = MR{kc}.Dist; % en nm
                B = MR{kc}.B; % en mT
                F = MR{kc}.F; % en pN
                T = MR{kc}.time; % en s
                Timg = MR{kc}.timeimg; % en s
                S = MR{kc}.S;
                CellName = MR{kc}.cellname;
                
                fprintf(['Manip : ' CellName '_' ExpType '\n\n'])
                
                %% Isolation de la partie oscillante du signal
                % la plus longue sans images qui sautent
                
                fprintf('Isolation des oscillations...');
                
                StartOscInd = find(diff(F)>10,1) + 1;
                % Le saut du champ constant au champ oscillant
                
                SignalS    = S(StartOscInd:end); % Support Image
                SignalT    = T(StartOscInd:end); % Temps
                SignalTimg    = Timg(StartOscInd:end); % Temps
                SignalDist = Dist(StartOscInd:end); % Distance
                SignalF    = F(StartOscInd:end); % Force
                SignalB = B(StartOscInd:end); % Champ
                % Sur la partie oscillante
                
                % Verification de la coehrence des temps entre camera TTL
                % et image Metadonnée
                
                Tdiff = SignalT - SignalTimg;
                if sum(Tdiff) ~= 0
                    error('Incohérence des temps !!!')
                end
                
                % Identification de la partie la plus longue sans saut d'image
                DS = diff(SignalS) == 1; % 1 quand continue, 0 sinon
                zpos = find(~[0 DS' 0]); % position de tout les 0
                [~,idx] = max(diff(zpos)); % plus grande distance entre deux 0
                IndPartC0 = zpos(idx)+50:zpos(idx+1)-2; % Indices de la plus grande partie continue
                % (50pts de moins au début pour eviter la compression initiale)
                
                SignalT_C0    = SignalT(IndPartC0) - SignalT(IndPartC0(1));
                SignalDist_C0 = SignalDist(IndPartC0);
                SignalF_C0    = SignalF(IndPartC0);
                SignalB_C0    = SignalB(IndPartC0);
                
                fprintf(' OK\n\n');
                
                %% Transformée de fourier et filtrage en fréquence des courbes de force et distance
                
                fprintf('Récupération des basses fréquence...');
                
                % Fréquences pertinentes
                
                FreqEch = 1/mean(diff(SignalT_C0)); % Frequence moyenne de prise des images
                
                FreqEch = precisround(FreqEch,1e-4);
                
                
                
                FtChamp = fft(SignalB_C0-mean(SignalB_C0));
                
                Freqs = (0:length(FtChamp)-1)*FreqEch/length(FtChamp);
                
                Freqs = precisround(Freqs,1e-4);
                
                FreqExp = Freqs(find(abs(FtChamp)>max(abs(FtChamp)/2),1));
                
                FreqExp = precisround(FreqExp,1e-4); % frequence de solicitation
                
                
                % Transformée fourier pour filtrage
                
                FtDist = fft(SignalDist_C0);
                
                FtForce = fft(SignalF_C0);
                
                
                
                
                FtDistBF = FtDist;
                FtForceBF = FtForce;
                
                FtDistBF(Freqs>max([0.1 FreqExp/10])&Freqs<FreqEch-max([0.1 FreqExp/10])) = 0;
                FtForceBF(Freqs>max([0.1 FreqExp/10])&Freqs<FreqEch-max([0.1 FreqExp/10])) = 0;
                
                InvFtDistBF = ifft(FtDistBF);
                InvFtForceBF = ifft(FtForceBF);
                
                
                fprintf(' OK\n\n');
                
                fprintf('Calcul de la contrainte et de la déformation...');
                
                %% Calcul de H0 et Delta (indentation) le long de la courbe BF
                
                EChad = InfoTable.Echad(find(strcmp(InfoTable.cellname,CellName))); % E pour cette cellule en Pa
                
                F0 = mean(InvFtForceBF)*1e-12; % Force moyenne sur cette manip (en N) sans oscilation
                
                R = Diam/2*1e-9; % rayon bille en m
                
                h = InvFtDistBF*1e-9; % épaisseur sans oscillation en m
                
                
                H0 = H0ChadCalc(EChad,F0,R,h); % H0 par chadwick inverse
                
                
                
                %% Calcul de la contrainte et de la deformation en temps
                
                DeltaOsc = H0 - SignalDist_C0; % Indentation en nm
                SurfCont = pi*Diam/2*DeltaOsc*1e-6; % Surface de contact en µm²
                
                Sigma = SignalF_C0./SurfCont; % Contrainte en Pa
                Epsilon = DeltaOsc./H0; % Deformation en %
                
                
                
                %% Sauvegarde
                
                fprintf('Sauvegarde partielle...');
                
                MF{kk}.time = SignalT_C0;
                MF{kk}.force = SignalF_C0;
                MF{kk}.field = SignalB_C0;
                MF{kk}.Sigma = Sigma;
                MF{kk}.Epsilon = Epsilon;
                MF{kk}.RawSignal = SignalDist_C0;
                MF{kk}.BaseSignal = InvFtDistBF;
                
                MF{kk}.ExpFreq = FreqExp;
                MF{kk}.FreqsFourier = Freqs;
                MF{kk}.FreqEch = FreqEch;
                
                MF{kk}.CellName = CellName;
                MF{kk}.Exp = ExpType;
                MF{kk}.EChad = EChad;
                
                kk = kk +1;
                
                cd(sf)
                save(savenamepart,'MF');
                cd(path)
                
                fprintf(' OK\n\n\n');
                
                
            end
            
        end
        
    end
end
fprintf('Sauvegarde complete...');
cd(sf)
save(savename,'MF');
delete([savenamepart '.mat'])
cd(path)

fprintf(' OK\n\n\n');

end




function out = precisround(n,precis)

out = round(n./precis).*precis;

end


function H0 = H0ChadCalc(E,F0,R,h)


% A partir de l'équation de Chadwick :
% F = (pi*E*R*delta^2)/(3*h0)
%
% On écrit delta = h0 - h et on resoud une
% equation du second degré pour h0
% a partir des valeurs de h (courbe BF) et EChad,
% et de la force sans l'oscilation F0
%
% h0²*pi*EChad*R - h0*(3*F0 + 2*h*pi*E*R) + pi*EChad*R*h² = 0





a = pi*E*R; % 'a' de l'eq second degree

b = -(3.*F0 + 2.*a.*h); % 'b' de l'eq second degree

c = a.*h.^2; % 'c' de l'eq second degree

disc = b.^2 - 4.*a.*c; % discriminant eq second degre en h0



h0m = (-b - sqrt(disc))./(2*a);

if ~all(h0m<h)
    error('Probleme avec le calcul de h0')
end

h0p = (-b + sqrt(disc))./(2*a); % racine + de l'eq, la seule plus grande que h

if ~all(h0p>h)
    error('Probleme avec le calcul de h0')
end

H0 = h0p*1e9; % épaisseur non déformé en nm

end