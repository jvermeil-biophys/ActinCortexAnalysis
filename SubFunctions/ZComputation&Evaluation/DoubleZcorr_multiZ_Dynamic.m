function [Z1, Z2] = DoubleZcorr_multiZ_Dynamic(X1tot,X1,Y1tot,Y1,X2tot,...
                X2,Y2tot,Y2,Stot,S,stackname,Sdown,Smid,Sup,depthoname,depthofolder)
% renvoi Z1 et Z2 en micron 
tic

set(0,'DefaultFigureWindowStyle','docked')



Ls = length(S); % nombre de positions
l = 300; % longeur de la ligne a correler
l2 = 350;
Leng = l2/10;

% image de référence et normalization ligne par ligne
load([depthofolder filesep 'EtalonageZ\' depthoname '.mat']);
if exist('Ktot')
    K = Ktot;
end
if exist('fmin')
    f = fmin;
end
Lk = length(K);
I = K(:,ceil(Lk/2)-l:ceil(Lk/2)+l);
Ime = mean(I,2);
RIme = double(repmat(Ime,1,2*l+1));
Id = double(I);
I = Id./RIme;

sk = size(K,1);

% variable finale
Z1 = nan(Ls,1);
Z2 = nan(Ls,1);


Sdone = [];

DICS1t = nan(sk,3);


% on parcours les images
for i = 1:Ls
    if not(ismember(S(i),Sdone))
        
        if ismember(S(i),Sdown)
            Slist = [S(i), S(i)+1, S(i)+2];
        elseif ismember(S(i),Smid)
            Slist = [S(i)-1, S(i), S(i)+1];
        elseif ismember(S(i),Sup)
            Slist = [S(i)-2, S(i)-1, S(i)];
        end
        
        MultiPos1 = zeros(3,2);
        MultiPos2 = zeros(3,2);
        
        for ii = 1:3
                 
            ptrpos = find(Stot == Slist(ii));
            
            if isempty(ptrpos)
                X1tmp = X1(i);
                Y1tmp = Y1(i);
                
                X2tmp = X2(i);
                Y2tmp = Y2(i);
            else
                X1tmp = X1tot(ptrpos);
                Y1tmp = Y1tot(ptrpos);
                
                X2tmp = X2tot(ptrpos);
                Y2tmp = Y2tot(ptrpos);
            end
            
            rX1tmp = round(X1tmp);
            rY1tmp = round(Y1tmp); 
            
            rX2tmp = round(X2tmp);
            rY2tmp = round(Y2tmp);            
            
            

            ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-3 rY1tmp+3],[rX1tmp-Leng-2 rX1tmp+Leng+2]}); % ROI verticale
            
            if length(ROI1) < 75 
                    ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % ROI horizontale
            end
            
            ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-3 rY2tmp+3],[rX2tmp-Leng-2 rX2tmp+Leng+2]}); % ROI verticale
                        
             if length(ROI2) < 75 
                    ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-Leng-2 rY2tmp+Leng+2],[rX2tmp-3 rX2tmp+3]})'; % ROI horizontale
             end

            
             
            ROI1d = double(ROI1);
            
            ROI2d = double(ROI2);
            
            X1inRoi = round((X1tmp - rX1tmp +Leng +2 +1)*10);
            Y1inRoi = round((Y1tmp - rY1tmp +3 +1)*10); 
            
            X2inRoi = round((X2tmp - rX2tmp +Leng +2 +1)*10);
            Y2inRoi = round((Y2tmp - rY2tmp +3 +1)*10); 
            
            
            [X1roi,Y1roi] = meshgrid(rX1tmp-Leng-2:rX1tmp+Leng+2,rY1tmp-3:rY1tmp+3);
            [X1roi2,Y1roi2] = meshgrid(rX1tmp-Leng-2:0.1:rX1tmp+Leng+2,rY1tmp-3:0.1:rY1tmp+3);

            [X2roi,Y2roi] = meshgrid(rX2tmp-Leng-2:rX2tmp+Leng+2,rY2tmp-3:rY2tmp+3);
            [X2roi2,Y2roi2] = meshgrid(rX2tmp-Leng-2:0.1:rX2tmp+Leng+2,rY2tmp-3:0.1:rY2tmp+3);

            ROI1int = interp2(X1roi,Y1roi,ROI1d,X1roi2,Y1roi2,'spline');              
            L1 =  ROI1int(Y1inRoi,X1inRoi-Leng*10:X1inRoi+Leng*10);          

            ROI2int = interp2(X2roi,Y2roi,ROI2d,X2roi2,Y2roi2,'spline');              
            L2 =  ROI2int(Y2inRoi,X2inRoi-Leng*10:X2inRoi+Leng*10);          

            
            
            [Ml1,~] = max(L1(0.5*l2:1.5*l2));
            [ml1,~] = min(L1(0.5*l2:1.5*l2));
            
            [Ml2,~] = max(L2(0.5*l2:1.5*l2));
            [ml2,~] = min(L2(0.5*l2:1.5*l2));
            
            locsL1 = 1:length(L1);
            
            locsL2 = 1:length(L2);
            
            L1mask = (L1 > 0.5*(Ml1-ml1)+ml1)&(locsL1>l2/2)&(locsL1<1.5*l2);
            
            L2mask = (L2 > 0.5*(Ml2-ml2)+ml2)&(locsL2>l2/2)&(locsL2<1.5*l2);
            
            midL1 = round(mean(locsL1(L1mask).*(L1(L1mask)/mean(L1(L1mask)))));
            
            midL2 = round(mean(locsL2(L2mask).*(L2(L2mask)/mean(L2(L2mask)))));
            
            
            %bille 1
            if midL1+l>2*l2+1
                midL1 = 2*l2+1-l;
                cprintf('*err','midL shift !\n\n')                
               
            elseif midL1-l<1
                midL1 = l +1;                
                cprintf('*err','midL shift !\n\n')
            end
            
            try
                L1 = L1(midL1-l:midL1+l);
            catch
                
                cprintf('*err','Error with midL !!!\n\n')
                pause
            end
            
            
            %bille 2
             if midL2+l>2*l2+1
                midL2 = 2*l2+1-l;
                cprintf('*err','midL shift !\n\n')                
               
            elseif midL2-l<1
                midL2 = l +1;                
                cprintf('*err','midL shift !\n\n')
            end
            
            try
                L2 = L2(midL2-l:midL2+l);
            catch
                
                cprintf('*err','Error with midL !!!\n\n')
                pause
            end
 
            
            %bille 1
            L1 = (L1 + L1(end:-1:1))/2;
            
            L1 = L1/mean(L1);
            

            ML1 = repmat(L1,sk,1); % replication pour avoir la meme taille qu'image ref
                        
            DI1 = I - ML1; % difference des deux
                        
            DIC1 = DI1.^2; % au carré
                        
            DICS1 = sum(DIC1,2); % somme des carrée

            DICS1t(:,ii) = DICS1;
            
            x = 1:length(DICS1);
            xq = 1:0.1:length(DICS1);
            
            DICS1i = interp1(x,DICS1,xq);
            
            [~, locs1] = findpeaks(-smooth(DICS1i),xq,'minpeakdistance',5,'sortstr','descend');
            
            MultiPos1(ii,:) = locs1(1:2);
            
            % bille 2
            
            L2 = (L2 + L2(end:-1:1))/2;
            
            L2 = L2/mean(L2);
            

            ML2 = repmat(L2,sk,1); % replication pour avoir la meme taille qu'image ref
                        
            DI2 = I - ML2; % difference des deux
                        
            DIC2 = DI2.^2; % au carré
                        
            DICS2 = sum(DIC2,2); % somme des carrée

            DICS2t(:,ii) = DICS2;
            
            x = 1:length(DICS2);
            xq = 1:0.1:length(DICS2);
            
            DICS2i = interp1(x,smooth(DICS2),xq,'spline');
            
            [~, locs2] = findpeaks(-smooth(DICS2i),xq,'minpeakdistance',5,'sortstr','descend');
            
            
        end
        
        Sdone = [Sdone Slist];

% on veut multipos(:,1/2) décroissant -> l'indice 1 correspond la curvelo,
% la correlation pour l'image du bas de la bille, qui doit conc correler
% vers la fin du depthograph, et donc avoir l'indice le plus grand

% on teste les combinaison des 2 maximum qui verifie les critère d'ordre+
[Bool1, Bool2] = deal(nan(1,8));
Str = cell(1,8);
%bille 1
for k1 = 1:2
        for k2 = 1:2
            for k3 = 1:2
                k = bi2de([k3-1,k2-1,k1-1])+1;
                Str{k} = [num2str(k1) num2str(k2) num2str(k3)];
                Bool1(k) = (MultiPos1(1,k1)-MultiPos1(2,k2) < 25) && (MultiPos1(2,k2)-MultiPos1(3,k3)<25) ... % saut de moins de 500 nm
        && (MultiPos1(1,k1)>MultiPos1(2,k2)) && (MultiPos1(2,k2)>MultiPos1(3,k3)); % saut dans le bon sens
            end
        end
end

%bille 2
for k1 = 1:2
        for k2 = 1:2
            for k3 = 1:2
                k = bi2de([k3-1,k2-1,k1-1])+1;
                Str{k} = [num2str(k1) num2str(k2) num2str(k3)];
                Bool2(k) = (MultiPos2(1,k1)-MultiPos2(2,k2) < 25) && (MultiPos2(2,k2)-MultiPos2(3,k3)<25) ... % saut de moins de 500 nm
        && (MultiPos2(1,k1)>MultiPos2(2,k2)) && (MultiPos2(2,k2)>MultiPos2(3,k3)); % saut dans le bon sens
            end
        end
end

%bille 1
if Bool1(1) == 1
    ShiftLo1 = MultiPos1(1,1)-MultiPos1(2,1);
    ShiftUp1 = MultiPos1(2,1)-MultiPos1(3,1); 
elseif sum(Bool1) == 1
    n = find(Bool1 ==1);
    k1 = str2double(Str{n}(1));
    k2 = str2double(Str{n}(2));
    k3 = str2double(Str{n}(3));
    ShiftLo1 = MultiPos1(1,k1)-MultiPos1(2,k2);
    ShiftUp1 = MultiPos1(2,k2)-MultiPos1(3,k3); 
elseif sum(Bool1) > 1
    %%% que faire ?
    error('Multiple correct combination of correlation 1 !')
elseif sum(Bool1) == 0
    %%% que faire ?
    error('No correct combination of correlation 1 !')
end


%bille 2
if Bool2(1) == 1
    ShiftLo2 = MultiPos2(1,1)-MultiPos2(2,1);
    ShiftUp2 = MultiPos2(2,1)-MultiPos2(3,1); 
elseif sum(Bool2) == 1
    n = find(Bool2 ==1);
    k1 = str2double(Str{n}(1));
    k2 = str2double(Str{n}(2));
    k3 = str2double(Str{n}(3));
    ShiftLo2 = MultiPos2(1,k1)-MultiPos2(2,k2);
    ShiftUp2 = MultiPos2(2,k2)-MultiPos2(3,k3); 
elseif sum(Bool2) > 1
    %%% que faire ?
    error('Multiple correct combination of correlation 2 !')
elseif sum(Bool2) == 0
    %%% que faire ?
    error('No correct combination of correlation 2 !')
end

ShiftLo = round((ShiftLo1+ShiftLo2)/2);
ShiftUp = round((ShiftUp1+ShiftUp2)/2);


% bille1
curvelo1 = DICS1t(:,1); % correlation de l'image du bas de la bille (point blanc qui perd du contraste)
curveloshift1 = [curvelo1; DICS1t(end,1)*ones(ShiftLo+ShiftUp,1)];

curvemid1 = DICS1t(:,2);
curvemidshift1 = [DICS1t(1,2)*ones(ShiftLo,1); curvemid1; DICS1t(end,2)*ones(ShiftUp,1)];

curveup1 = DICS1t(:,3); % correlation de l'image du haut de la bille (anneau)
curveupshift1 = [DICS1t(1,3)*ones(ShiftLo+ShiftUp,1); curveup1];


try
        curve1 = curveloshift1 + curvemidshift1 + curveupshift1;
catch
    lol
end
sup = 1:length(curve1);

supint = 1:0.01:length(curve1);

curve1int = interp1(sup,curve1,supint,'spline');

[~,indmin1] = min(curve1int);
icm1 = supint(indmin1); % indice of curve min

Z1(i) = (((icm1-ShiftLo)-f)*stepnm)/1000; % step en nm entre chaque image
      

% bille 2

curvelo2 = DICS2(:,1); % correlation de l'image du bas de la bille (point blanc qui perd du contraste)
curveloshift2 = [curvelo2; DICS2t(end,1)*ones(ShiftLo+ShiftUp,1)];

curvemid2 = DICS2t(:,2);
curvemidshift2 = [DICS2t(1,2)*ones(ShiftLo,1); curvemid2; DICS2t(end,2)*ones(ShiftUp,1)];

curveup2 = DICS2t(:,3); % correlation de l'image du haut de la bille (anneau)
curveupshift2 = [DICS2t(1,3)*ones(ShiftLo+ShiftUp,1); curveup2];



        curve2 = curveloshift2 + curvemidshift2 + curveupshift2;

sup = 1:length(curve2);

supint = 1:0.01:length(curve2);

curve2int = interp1(sup,curve2,supint,'spline');

[~,indmin2] = min(curve2int);
icm2 = supint(indmin2); % indice of curve min

Z2(i) = ((icm2 - ShiftLo-f)*stepnm)/1000; % step en nm entre chaque image
        
        
 
        if 0
        figure
        plot(curve1)
        hold on
        plot([DICS2t(1,3)*ones(2*Zstppx,1); DICS2t(:,3)])
        plot([DICS2t(1,2)*ones(Zstppx,1); DICS2t(:,2); DICS2t(end,2)*ones(Zstppx,1)])
        plot([DICS2t(:,1);DICS2t(end,1)*ones(2*Zstppx,1)])
        plot(icm,cm,'*r')
        title([stackname '_' num2str(S(i))],'interpreter','none')
        pause(0.05)
        end
        
        clear Slist
    end
end
toc