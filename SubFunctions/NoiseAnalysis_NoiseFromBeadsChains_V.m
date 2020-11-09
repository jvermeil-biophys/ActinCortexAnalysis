function NoiseAnalysis_NoiseFromBeadsChains_V(datafolder,names,depthoDirPath,depthoname,isImageTriplet,scale,DIAMETER,PLOTVID)

set(0,'DefaultFigureWindowStyle','docked') % 'normal'

% dossier de données et sauvegarde
datenow = datestr(now,'yy-mm-dd');
% DIAMETER = 4503; Argument !

% scale est un argument maintenant !
% scale=1/15.8; % 100X
% scale=1/9.7; % 63X


% Variables globales pour plot final sur toute les vid�os
MedianD2 = [];
StdD2 = [];
MedianD3 = [];
StdD3 = [];

for k = 1:length(names)
    
    filename = names{k};
    
    imageStackName = [filename '.tif'];
    imageStackPath = [datafolder filesep imageStackName];
    doesStackExist = exist(imageStackPath,'file');
    
    
    restype = "Results";
    ijResultsPath = strcat(datafolder,filesep,filename,'_',restype,'.txt');
    doesResultExist = exist(ijResultsPath,'file');
    
    % fieldPath = [datafolder filesep filename '_Field.txt'];
    % doesFieldExist = exist([path filesep fieldPath],'file');
    
    savename = ['NoiseAnalysis_' filename '.mat'];
    
    % on regarde si cette vid�o a d�j� �t� trait�e
    isAnalyzedYet = exist(savename,'file');
    
    nimg = 3; % Images par boucle, lorsque isImageTriplet = 1.
    
    if doesStackExist == 2 && doesResultExist == 2 && isAnalyzedYet == 0 % && doesFieldExist
        
        fprintf('Lecture du fichier Results...');
        ResultsMatrix=dlmread(ijResultsPath,'\t',1,0);
        fid = fopen(ijResultsPath);
        HeadAll = fgetl(fid);
        fclose(fid);
        HeadCell = textscan(HeadAll,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
        HeadTable = [HeadCell{:}];
        % Useful parameters
        Nxm = find(strcmp(HeadTable,'XM')) +1 ;
        Nym = find(strcmp(HeadTable,'YM')) +1 ;
        Ns = find(strcmp(HeadTable,'Slice')) +1 ;
        Nstd = find(strcmp(HeadTable,'StdDev')) +1 ;
        % Deprecated parameters
        Nx = find(strcmp(HeadTable,'X')) +1 ;
        Ny = find(strcmp(HeadTable,'Y')) +1 ;
        Nmaj = find(strcmp(HeadTable,'Major')) +1 ;
        Nmin = find(strcmp(HeadTable,'Minor')) +1 ;
        
        fprintf(' OK\n');
        
        
        fprintf('Verification de la vidéo...');
        
        % on sépare les point selon le Z de mesure
        S = ResultsMatrix(:,Ns);
        Sfull = 1:max(S);
        [Sdown,Smid,Sup] = SplitZimage(Sfull,nimg);
        [Smd,ptrd,ptrm] = intersect(Sdown+1,Smid);
        [Smdu,ptrmd,ptru] = intersect(Smd,Sup-1);
        
        Ntriplets = min([length(Sdown) length(Smid) length(Sup)]);
        Sdown = Sdown(ptrd(ptrmd));
        Smid  = Smdu ;
        Sup   = Sup(ptru)  ;
        Sall = [Sdown;Smid;Sup];
        
        for is = 1:max(S)
            try
                Itmp = imread(imageStackPath,'tiff',is);
            catch
                themall
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
        
        
        % Beads tracking
        fprintf('Tracking des billes...');
        
        % trackingTable contient les coordonnées dans le plan de chaque bille.
        % C'est une structure où chaque bille possède les champs X, Y, ptr et S.
        % Pour le moment seuls les slices ou toutes les billes sont présentes
        % sont analysées ; donc à priori le champ S est identiques pour toutes les billes.
        trackingTable = NoiseAnalysis_TrackMultiBeadsHungarian_V([ResultsMatrix(:,Nxm) S],ResultsMatrix(:,Nym),imageStackPath,PLOTVID);
        NImagesAnalysed = length(trackingTable(1).X);
        
        % On va maintenant "traduire" la trackingTable en une pairsTable qui
        % contient des données pour chaque paire de billes (Npairs = Nbeads-1).
        % pairsTable(1) renvoie aux données de la paire '12' entre la bille 1
        % et la bille 2.
        Nbeads = length(trackingTable);
        Npairs = Nbeads-1;
        for k = 1:Npairs
            pairsTable(k).pair = [num2str(k) num2str(k+1)];
            pairsTable(k).firstBead = k;
            pairsTable(k).secondBead = k+1;
            pairsTable(k).S = trackingTable(k).S';
            % D2raw > distance D2 non améliorée par le choix de la meilleure
            % des 3 images, pour les films constitués de triplets d'images.
            % Pour les films sans triplets on prendra D2 = D2raw.
            pairsTable(k).D2raw = ((trackingTable(k+1).X' - trackingTable(k).X').^2 + (trackingTable(k+1).Y' - trackingTable(k).Y').^2).^(1/2);
            pairsTable(k).D2rawstd = std(pairsTable(k).D2raw);
        end
        
        
        
        % Passage du code abandonné : il n'est pas pertinent de faire des moyennes
        % sur les billes de la même chaine, car bien que leur diamètres soient
        % proches, ils sont suffisemment différents pour que la dispersion des
        % diamètres masque le bruit de mesure.
        
        %     globalDistanceTable.S = distanceTable(1).S;
        %     globalDistanceTable.D2mean = zeros(1, NImagesAnalysed);
        %     globalDistanceTable.D2std = zeros(1, NImagesAnalysed);
        %
        %     for i = 1:NImagesAnalysed
        %         D2_i = zeros(1, Npairs);
        %         for b = 1:Npairs
        %             D2_i(b) = distanceTable(b).D2(i);
        %         end
        %         globalDistanceTable.D2mean(i) = mean(D2_i);
        %         globalDistanceTable.D2std(i) = std(D2_i);
        %     end
        %     globalD2std = std(distanceTable(b).D2);
        %     1000*pixelconv*globalD2std
        
        %     figure
        %     hold on
        %     plot(globalDistanceTable.S, 1000*pixelconv*(globalDistanceTable.D2mean) - DIAMETER, 'b-');
        %     plot(globalDistanceTable.S, 1000*pixelconv*(globalDistanceTable.D2mean + globalDistanceTable.D2std) - DIAMETER, 'r-');
        %     plot(globalDistanceTable.S, 1000*pixelconv*(globalDistanceTable.D2mean - globalDistanceTable.D2std) - DIAMETER, 'r-');
        
        
        % 1er cas : si on a un film avec des triplets d'image, on va pouvoir
        % raffiner la distance D2 et calculer le dz avec ce triplet d'image.
        if isImageTriplet
            for k = 1:Npairs
                % trackingTable(k) > 1e bille (X1, Y1, STD1, etc) ;
                % trackingTable(k+1) > 2e bille (X2, Y2, STD2, etc)
                % pairsTable(k) > paire de ces deux billes.
                
                STD1 = nan(3,length(Sdown));
                STD2 = nan(3,length(Sdown));
                
                X1d = zeros(size(Sdown));
                X2d = zeros(size(Sdown));
                [~,pd1,p1d] = intersect(Sdown,trackingTable(k).S);
                [~,pd2,p2d] = intersect(Sdown,trackingTable(k+1).S);
                X1d(pd1) = trackingTable(k).X(p1d);
                STD1(1,pd1) = ResultsMatrix(trackingTable(k).ptr(p1d),Nstd);
                X2d(pd2) = trackingTable(k+1).X(p2d);
                STD2(1,pd2) = ResultsMatrix(trackingTable(k+1).ptr(p2d),Nstd);
                
                X1m = zeros(size(Smid));
                X2m = zeros(size(Smid));
                [~,pm1,p1m] = intersect(Smid,trackingTable(k).S);
                [~,pm2,p2m] = intersect(Smid,trackingTable(k+1).S);
                X1m(pm1) = trackingTable(k).X(p1m);
                STD1(2,pm1) = ResultsMatrix(trackingTable(k).ptr(p1m),Nstd);
                X2m(pm2) = trackingTable(k+1).X(p2m);
                STD2(2,pm2) = ResultsMatrix(trackingTable(k+1).ptr(p2m),Nstd);
                
                X1u = zeros(size(Sup));
                X2u = zeros(size(Sup));
                [~,pu1,p1u] = intersect(Sup,trackingTable(k).S);
                [~,pu2,p2u] = intersect(Sup,trackingTable(k+1).S);
                X1u(pu1) = trackingTable(k).X(p1u);
                STD1(3,pu1) = ResultsMatrix(trackingTable(k).ptr(p1u),Nstd);
                X2u(pu2) = trackingTable(k+1).X(p2u);
                STD2(3,pu2) = ResultsMatrix(trackingTable(k+1).ptr(p2u),Nstd);
                
                
                m = ismember(STD1+STD2,max(STD1+STD2)); % ici STD et STD2 ont chacun 3 lignes et S colonnes.
                % Dans ce cas max() renvoie 1 ligne avec le max de chaque colonne à chaque fois.
                % et m indique donc la ligne (down/mid/up) sur laquelle il faut lire pour trouver la meilleure valeur de X.
                
                X1all = [X1d; X1m; X1u];
                X2all = [X2d; X2m; X2u];
                
                % Les données importantes pour la suite sont le X1, X2, Y1, Y2, Svalid
                % pour la paire de bille k : on va donc les stocker dans la structure pairsTable.
                pairsTable(k).X1 = X1all(m);
                pairsTable(k).X2 = X2all(m);
                
                pairsTable(k).SValid = Sall(m);
                pointeurValid = ismember(pairsTable(k).S,pairsTable(k).SValid);
                pairsTable(k).Y1 = trackingTable(k).Y(pointeurValid)';
                pairsTable(k).Y2 = trackingTable(k+1).Y(pointeurValid)';
                
                % On recalcule aussi D2 avec ces meilleurs coordonnées.
                pairsTable(k).D2 = ((pairsTable(k).X2 - pairsTable(k).X1).^2 + (pairsTable(k).Y2 - pairsTable(k).Y1).^2).^(1/2);
                pairsTable(k).D2std = std(pairsTable(k).D2);
            end
            
            fprintf('Calcul de dz...\n');
            
            
            
            %         lastFrame = max(distanceTable(1).S);
            for k = 1:Npairs
                pairsTable(k).dz = DzCalc_Zcorr_multiZ(trackingTable(k).X,pairsTable(k).X1',...
                    trackingTable(k).Y,pairsTable(k).Y1',trackingTable(k+1).X,pairsTable(k).X2',...
                    trackingTable(k+1).Y,pairsTable(k).Y2',pairsTable(k).S',pairsTable(k).SValid',...
                    imageStackPath,Sdown,Smid,Sup,depthoname,depthoDirPath);
                % IMPORTANT : En sortie de DzCalc_Zcorr, dz est en nanomètres !
                % Pour avoir dz en pixels ("pixels verticaux") comme D2, on va
                % diviser dz par 1000*scale [nm par pixels].
                pairsTable(k).dz = pairsTable(k).dz'/(1000*scale);
                
            end
            
            
            % 2e cas : On n'a pas de triplets d'images : pas de raffinement des
            % coordonnées, pour D2 on reprend simplement D2raw, et pour le calcul
            % des dz on utilise le DzCalc_Zcorr (et pas DzCalc_Zcorr_multiZ) qui
            % infère la distance dz à partir d'une seule image (contre un triplet
            % sinon).
        else
            for k = 1:Npairs
                pairsTable(k).X1 = trackingTable(k).X';
                pairsTable(k).X2 = trackingTable(k+1).X';
                pairsTable(k).Y1 = trackingTable(k).Y';
                pairsTable(k).Y2 = trackingTable(k+1).Y';
                pairsTable(k).SValid = pairsTable(k).S;
                pairsTable(k).D2 = pairsTable(k).D2raw;
                pairsTable(k).D2std = pairsTable(k).D2rawstd;
            end
            
            fprintf('Calcul de dz...\n');
            
            for k = 1:Npairs
                pairsTable(k).dz = DzCalc_Zcorr(pairsTable(k).X1,pairsTable(k).Y1,pairsTable(k).X2,pairsTable(k).Y2,...
                    pairsTable(k).S,imageStackPath,depthoname,depthoDirPath);
                % IMPORTANT : En sortie de DzCalc_Zcorr, dz est en nanomètres !
                % Pour avoir dz en pixels ("pixels verticaux") comme D2, on va
                % diviser dz par 1000*scale [nm par pixels].
                pairsTable(k).dz = pairsTable(k).dz'/(1000*scale);
            end
        end
        
        % Finalement on calcule la distance 3D (en voxels).
        
        for k = 1:Npairs
            % Correction au cas où trop grand saut en DZ.
            % Je commence par copier les arrays qui vont subir cette
            % "relecture".
            pairsTable(k).dzValid_3D = pairsTable(k).dz;
            pairsTable(k).SValid_3D = pairsTable(k).SValid;
            pairsTable(k).D2Valid_3D = pairsTable(k).D2;
            
            ddz = diff(abs(pairsTable(k).dzValid_3D));
            ddzcount = 0;
            
            MAXIMUM_SLOPE = 5; % Les variations en dz plus brutales que ça seront coupées.
            while (max(abs(ddz)) > MAXIMUM_SLOPE) && (ddzcount<50)
                ptrdel = find(ddz > MAXIMUM_SLOPE) + 1;
                ptrdel = [ptrdel; find(ddz <  -MAXIMUM_SLOPE)];
                ptrdel = unique(ptrdel);
                
                pairsTable(k).SValid_3D(ptrdel) = [];
                pairsTable(k).D2Valid_3D(ptrdel) = [];
                pairsTable(k).dzValid_3D(ptrdel) = [];
                
                ddz=diff(abs(pairsTable(k).dzValid_3D));
                ddzcount = ddzcount +1;
            end
            
            %dz = dz*1.33/1.52; % correction pour la difference d'indice optique huile/eau
            pairsTable(k).D3 = ((pairsTable(k).D2Valid_3D).^2 + (pairsTable(k).dzValid_3D).^2).^(1/2);
            pairsTable(k).D3std = std(pairsTable(k).D3);
        end
        
        
        
        % Pour avoir des résultats plus faciles à interpréter que des pixels,
        % resultTableNanometers sauvegarde tous les résultats importants
        % convertis en nm.
        
        for k = 1:Npairs
            resultTableNanometers(k).isImageTriplet = isImageTriplet;
            resultTableNanometers(k).SValid = pairsTable(k).SValid;
            resultTableNanometers(k).D2raw = pairsTable(k).D2raw * (1000*scale);
            resultTableNanometers(k).D2rawstd = pairsTable(k).D2rawstd * (1000*scale);
            resultTableNanometers(k).D2 = pairsTable(k).D2 * (1000*scale);
            resultTableNanometers(k).D2std = pairsTable(k).D2std * (1000*scale);
            
            resultTableNanometers(k).dz = pairsTable(k).dz * (1000*scale);
            resultTableNanometers(k).SValid_3D = pairsTable(k).SValid_3D;
            resultTableNanometers(k).D2Valid_3D = pairsTable(k).D2Valid_3D * (1000*scale);
            resultTableNanometers(k).dzValid_3D = pairsTable(k).dzValid_3D * (1000*scale);
            resultTableNanometers(k).D3 = pairsTable(k).D3 * (1000*scale);
            resultTableNanometers(k).D3std = pairsTable(k).D3std * (1000*scale);
        end
        
        
        % version boucle
        %         D2median_plot = zeros(1,Npairs);
        %         D2std_plot = zeros(1,Npairs);
        %         D3median_plot = zeros(1,Npairs);
        %         D3std_plot = zeros(1,Npairs);
        %         for k = 1:Npairs
        %             D2median_plot(k) = median(resultTableNanometers(k).D2);
        %             D2std_plot(k) = resultTableNanometers(k).D2std;
        %             D3median_plot(k) = median(resultTableNanometers(k).D3);
        %             D3std_plot(k) = resultTableNanometers(k).D3std;
        %         end
        
        % version sans boucle
        D2median_plot = cellfun(@median,{resultTableNanometers.D2});
        D2std_plot = cellfun(@std,{resultTableNanometers.D2});
        D3median_plot = cellfun(@median,{resultTableNanometers.D3});
        D3std_plot = cellfun(@std,{resultTableNanometers.D3});
        
        % Enregistrement global
        MedianD2 = [MedianD2 D2median_plot];
        StdD2 = [StdD2 D2std_plot];
        MedianD3 = [MedianD3 D3median_plot];
        StdD3 = [StdD3 D3std_plot];
        
        
        %% figure sur les vid�os individuelles
        
        if PLOTVID
            figure
            for k = 1:Npairs
                subplot(ceil(Npairs/3),3,k)
                hold on
                title(['2D raw distance for pair ' pairsTable(k).pair ' ; std = ' num2str(1000*scale*pairsTable(k).D2rawstd)])
                plot(pairsTable(k).S, 1000*scale*(pairsTable(k).D2raw) - DIAMETER, 'b-');
            end
            
            % Comme D2 a été pas mal améliorée depuis D2raw, on la trace à nouveau ici.
            figure
            for k = 1:Npairs
                subplot(ceil(Npairs/3),3,k)
                hold on
                title(['2D distance for pair ' pairsTable(k).pair ' ; std = ' num2str(1000*scale*pairsTable(k).D2std)])
                plot(pairsTable(k).SValid, 1000*scale*(pairsTable(k).D2) - DIAMETER, 'b-');
            end
            
            % Dans tous les cas on finit par plotter la distance 3D.
            figure
            for k = 1:Npairs
                subplot(ceil(Npairs/3),3,k)
                hold on
                title(['3D distance for pair ' pairsTable(k).pair ' ; std = ' num2str(1000*scale*pairsTable(k).D3std)])
                plot(pairsTable(k).SValid_3D, 1000*scale*(pairsTable(k).D3) - DIAMETER, 'b-');
            end
            
            figure
            hold on
            color = 'b';
            subplot(2,2,1)
            [~,edges] = histcounts(D2median_plot,'binwidth',6);
            histogram(D2median_plot,edges,'facecolor',color)
            %     plot([median(D2median_plot) median(D2median_plot)],[0 26],'r--')
            %     ax = gca;
            %     ax.XTick = [0.01 0.05 0.1 0.5 1 3 10 15];
            %     xlim([0 55])
            %     ylim([0 25])
            xlabel('D2 median (nm)')
            ylabel('Count')
            
            subplot(2,2,2)
            [~,edges] = histcounts(D2std_plot,'binwidth',2);
            histogram(D2std_plot,edges,'facecolor',color)
            %     plot([median(D2median_plot) median(D2median_plot)],[0 26],'r--')
            %     ax = gca;
            %     ax.XTick = [0.01 0.05 0.1 0.5 1 3 10 15];
            %     xlim([0 55])
            %     ylim([0 25])
            xlabel('D2 std (nm)')
            ylabel('Count')
            
            subplot(2,2,3)
            [~,edges] = histcounts(D3median_plot,'binwidth',6);
            histogram(D3median_plot,edges,'facecolor',color)
            %     plot([median(D2median_plot) median(D2median_plot)],[0 26],'r--')
            %     ax = gca;
            %     ax.XTick = [0.01 0.05 0.1 0.5 1 3 10 15];
            %     xlim([0 55])
            %     ylim([0 25])
            xlabel('D3 median (nm)')
            ylabel('Count')
            
            subplot(2,2,4)
            [~,edges] = histcounts(D3std_plot,'binwidth',2);
            histogram(D3std_plot,edges,'facecolor',color)
            %     plot([median(D2median_plot) median(D2median_plot)],[0 26],'r--')
            %     ax = gca;
            %     ax.XTick = [0.01 0.05 0.1 0.5 1 3 10 15];
            %     xlim([0 55])
            %     ylim([0 25])
            xlabel('D3 std (nm)')
            ylabel('Count')
            
            fig = gcf;
            saveas(fig,[datafolder filesep 'NoiseAnalysisPlot_' filename '.png'],'png')
            saveas(fig,[datafolder filesep 'NoiseAnalysisPlot_' filename '.fig'],'fig')
            close
        end
        
        %% Sauvergarde des donn�es par vid�o
        save([datafolder filesep savename],'trackingTable','pairsTable','resultTableNanometers');
        
    elseif isAnalyzedYet == 2
        
        load([datafolder filesep savename])
        
        D2median_plot = cellfun(@median,{resultTableNanometers.D2});
        D2std_plot = cellfun(@std,{resultTableNanometers.D2});
        D3median_plot = cellfun(@median,{resultTableNanometers.D3});
        D3std_plot = cellfun(@std,{resultTableNanometers.D3});
        
        % Enregistrement global
        MedianD2 = [MedianD2 cellfun(@median,{resultTableNanometers.D2})- DIAMETER];
        StdD2 = [StdD2 cellfun(@std,{resultTableNanometers.D2})];
        MedianD3 = [MedianD3 cellfun(@median,{resultTableNanometers.D3})- DIAMETER];
        StdD3 = [StdD3 cellfun(@std,{resultTableNanometers.D3})];
    end
    
end

%% Figure r�capitulative des toutes les vid�os

figure

color = 'b';
subplot(2,2,1)
hold on
[~,edges] = histcounts(MedianD2,'binwidth',6);
h = histogram(MedianD2,edges,'facecolor',color);
plot([mean(MedianD2) mean(MedianD2)],[0 max(h.Values)+1],'r--','linewidth',1.5)
xlabel('D2 median (nm)')
ylabel('Count')

subplot(2,2,2)
hold on
[~,edges] = histcounts(StdD2,'binwidth',2);
h = histogram(StdD2,edges,'facecolor',color);
plot([mean(StdD2) mean(StdD2)],[0 max(h.Values)+1],'r--','linewidth',1.5)
xlabel('D2 std (nm)')
ylabel('Count')

subplot(2,2,3)
hold on
[~,edges] = histcounts(MedianD3,'binwidth',6);
h = histogram(MedianD3,edges,'facecolor',color);
plot([mean(MedianD3) mean(MedianD3)],[0 max(h.Values)+1],'r--','linewidth',1.5)
xlabel('D3 median (nm)')
ylabel('Count')

subplot(2,2,4)
hold on
[~,edges] = histcounts(StdD3,'binwidth',2);
h = histogram(StdD3,edges,'facecolor',color);
plot([mean(StdD3) mean(StdD3)],[0 max(h.Values)+1],'r--','linewidth',1.5)
xlabel('D3 std (nm)')
ylabel('Count')


end
