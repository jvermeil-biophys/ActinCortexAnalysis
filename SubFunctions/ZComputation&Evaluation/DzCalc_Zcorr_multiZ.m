function dzcalc = DzCalc_Zcorr_multiZ(X1tot,X1,Y1tot,Y1,X2tot,...
    X2,Y2tot,Y2,Stot,S,stackname,Sdown,Smid,Sup,depthoname,depthofolder)
% renvoi Z1 et Z2 en micron
% tic

set(0, 'defaultfigureposition', get(0, 'Screensize').*[20 20 0.98 0.9]);


Ls = length(S); % nombre de positions
l = 300; % longeur de la ligne a correler
l2 = 330;
Leng = l2/10;

% image de référence et normalization ligne par ligne
load([depthofolder filesep 'EtalonnageZ\' depthoname '.mat']);
if exist('Ktot')
    K = Ktot;
end
if exist('fmin')
    f = fmin;
end
Lk = length(K);

I = K(:,ceil(Lk/2)-l:ceil(Lk/2)+l);

Id = double(I);
Ime = mean(Id,2);
RIme = double(repmat(Ime,1,2*l+1));
I = Id./RIme;

sk = size(K,1);

firstimg = imread(stackname,'tiff',1);
simg = size(firstimg);

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
        MultiVal1 = zeros(3,2);
        MultiVal2 = zeros(3,2);
        
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

            ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % ROI verticale
            ROI1d = double(ROI1);
            if length(ROI1) < 2*Leng+5
                ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-3 rY1tmp+3],[rX1tmp-Leng-2 rX1tmp+Leng+2]}); % ROI horizontale
                ROI1d = double(ROI1);
            end
            
            if length(ROI1) < 2*Leng+5
                ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % ROI verticale
                ROI1d = double(ROI1);
                
                [Xmeshtmp,Ymeshtmp] =  meshgrid(rY1tmp-Leng-2:rY1tmp+Leng+2,rX1tmp-3:rX1tmp+3);
                
                [x,y] = find(Xmeshtmp>0&Xmeshtmp<simg(2)&Ymeshtmp>0&Ymeshtmp<simg(1));
                
                ROItmp = interp2(Xmeshtmp(unique(x),unique(y)),Ymeshtmp(unique(x),unique(y)),ROI1d,Xmeshtmp,Ymeshtmp,'spline');
                
                ROI1d = ROItmp;
                
            end
            
            
            
            ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-Leng-2 rY2tmp+Leng+2],[rX2tmp-3 rX2tmp+3]})';% ROI verticale
            ROI2d = double(ROI2);
            if length(ROI2) < 2*Leng+5
                ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-3 rY2tmp+3],[rX2tmp-Leng-2 rX2tmp+Leng+2]});% ROI horizontale
                ROI2d = double(ROI2);
                
            end
            
            if length(ROI2) < 2*Leng+5
                ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-Leng-2 rY2tmp+Leng+2],[rX2tmp-3 rX2tmp+3]})';% ROI verticale
                
                ROI2d = double(ROI2);
                
                [Xmeshtmp,Ymeshtmp] =  meshgrid(rY2tmp-Leng-2:rY2tmp+Leng+2,rX2tmp-3:rX2tmp+3);
                
                [x,y] = find(Xmeshtmp>0&Xmeshtmp<simg(2)&Ymeshtmp>0&Ymeshtmp<simg(1));
                

                ROItmp = interp2(Xmeshtmp(unique(x),unique(y)),Ymeshtmp(unique(x),unique(y)),ROI2d,Xmeshtmp,Ymeshtmp,'spline');

                ROI2d = ROItmp;
                
            end
            
            
            
            
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
                %                 cprintf('*err','midL shift !\n\n')
                
            elseif midL1-l<1
                midL1 = l +1;
                %                 cprintf('*err','midL shift !\n\n')
            end
            
            try
                L1 = L1(midL1-l:midL1+l);
            catch
                
                cprintf('*err','Error with midL !!!\n\n')
                
            end
            
            
            %bille 2
            if midL2+l>2*l2+1
                midL2 = 2*l2+1-l;
                %                 cprintf('*err','midL shift !\n\n')
                
            elseif midL2-l<1
                midL2 = l +1;
                %                 cprintf('*err','midL shift !\n\n')
            end
            
            try
                L2 = L2(midL2-l:midL2+l);
            catch
                
                cprintf('*err','Error with midL !!!\n\n')
                pause
            end
            
            L1 = (L1 + L1(end:-1:1))/2;
            L2 = (L2 + L2(end:-1:1))/2;
            
            %bille 1
            [~,MPtmp,MVtmp] = ComputeZCorrel(L1,I);
            
            MultiPos1(ii,:) = MPtmp;
            MultiVal1(ii,:) = MVtmp;
            
            clear MVtmp MPtmp
            % bille 2
            [~,MPtmp,MVtmp] = ComputeZCorrel(L2,I);
            
            MultiPos2(ii,:) = MPtmp;
            MultiVal2(ii,:) = MVtmp;
            
            clear MVtmp MPtmp
        end
        
        Sdone = [Sdone Slist];
        
        % on veut multipos(:,1/2) décroissant -> l'indice 1 correspond la curvelo,
        % la correlation pour l'image du bas de la bille, qui doit conc correler
        % vers la fin du depthograph, et donc avoir l'indice le plus grand
        
        % sorting correlations with alignement constraints (+700 nm max)
        [Bool1,Bool2,ValSum1,ValSum2] = deal(NaN(1,8));
        Str = cell(1,8);
        
        for k1 = 1:2
            for k2 = 1:2
                for k3 = 1:2
                    k = bi2de([k3-1,k2-1,k1-1])+1;
                    Str{k} = [num2str(k1) num2str(k2) num2str(k3)];
                    Bool1(k) = (MultiPos1(1,k1)>MultiPos1(2,k2)) &&...
                        (MultiPos1(2,k2)>MultiPos1(3,k3))&&(MultiPos1(1,k1)<35 + MultiPos1(2,k2)) &&...
                        (MultiPos1(2,k2)< 35 +MultiPos1(3,k3)); % saut dans le bon sens
                    ValSum1(k) = MultiVal1(1,k1)+MultiVal1(2,k2)+MultiVal1(3,k3);
                    Bool2(k) = (MultiPos2(1,k1)>MultiPos2(2,k2)) &&...
                        (MultiPos2(2,k2)>MultiPos2(3,k3))&&(MultiPos2(1,k1)<35+MultiPos2(2,k2)) &&...
                        (MultiPos2(2,k2)<35+MultiPos2(3,k3)); % saut dans le bon sens
                    ValSum2(k) = MultiVal2(1,k1)+MultiVal2(2,k2)+MultiVal2(3,k3);
                end
            end
        end
        clear k1 k2 k3
        
        
        % getting good Zs with aligned correlation
        
        
        [Zup1,Wup1,Zmid1,Wmid1,Zdown1,Wdown1] = Corr2Zs(Bool1,Str,MultiPos1,MultiVal1,ValSum1,f);
        [Zup2,Wup2,Zmid2,Wmid2,Zdown2,Wdown2] = Corr2Zs(Bool2,Str,MultiPos2,MultiVal2,ValSum2,f);
        
        dzup(i) = Zup2-Zup1;
        Wdzup(i) = Wup1+Wup2;
        dzmid(i) = Zmid2-Zmid1;
        Wdzmid(i) = Wmid1+Wmid2;
        dzdown(i) = Zdown2-Zdown1;
        Wdzdown(i) = Wdown1+Wdown2;
        
    end
end

dzcalc = sum([dzup.*Wdzup; dzmid.*Wdzmid; dzdown.*Wdzdown],1)./(sum([Wdzup; Wdzmid; Wdzdown],1))*stepnm;


% toc