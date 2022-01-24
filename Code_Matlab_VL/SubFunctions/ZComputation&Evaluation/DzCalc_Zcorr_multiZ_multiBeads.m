function dzcalc = DzCalc_Zcorr_multiZ_multiBeads(Xbigtot,Xbig,Ybigtot,Ybig,Xsmalltot,...
                Xsmall,Ysmalltot,Ysmall,Stot,S,stackname,Sdown,Smid,Sup,depthonameBig,depthonameSmall,depthofolder)
% renvoi dz 
% tic

set(0, 'defaultfigureposition', get(0, 'Screensize').*[20 20 0.98 0.9]);


lsmall = 80; % longeur de la ligne a correler
l2small = 90;
Lengsmall = l2small/10;


% image de référence et normalization ligne par ligne
load([depthofolder filesep 'EtalonnageZ\' depthonameSmall '.mat']);
if exist('Ktot')
    KSmall = Ktot;
elseif exist('K')    
    KSmall = K;
end
if exist('fmin')
    fSmall = fmin;
elseif exist('f')
    fSmall = f;
end
Lksmall = size(KSmall,2);
Ismall = KSmall(:,ceil(Lksmall/2)-lsmall:ceil(Lksmall/2)+lsmall);
Idsmall = double(Ismall);
Imesmall = mean(Idsmall,2);
RImesmall = double(repmat(Imesmall,1,2*lsmall+1));
Ismall = Idsmall./RImesmall;



lbig = 450; % longeur de la ligne a correler
l2big = 500;
Lengbig = l2big/10;

% image de référence et normalization ligne par ligne
load([depthofolder filesep 'EtalonnageZ\' depthonameBig '.mat']);
if exist('Ktot')
    KBig = Ktot;
elseif exist('K')    
    KBig = K;
end
if exist('fmin')
    fBig = fmin;
elseif exist('f')
    fBig = f;
end
Lkbig = size(KBig,2);
Ibig = KBig(:,ceil(Lkbig/2)-lbig:ceil(Lkbig/2)+lbig);
Idbig = double(Ibig);
Imebig = mean(Idbig,2);
RImebig = double(repmat(Imebig,1,2*lbig+1));
Ibig = Idbig./RImebig;



Sdone = [];

Ls = length(S); % nombre de positions

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

%% Grosse bille

       

        
        MultiPos1 = zeros(3,2);
        MultiVal1 = zeros(3,2);
        
        for ii = 1:3
                 
            ptrpos = find(Stot == Slist(ii));
            
            if isempty(ptrpos)

                X1tmp = Xbig(i);

                Y1tmp = Ybig(i);

            else
                X1tmp = Xbigtot(ptrpos);
                Y1tmp = Ybigtot(ptrpos);
                
            end
            
            rX1tmp = round(X1tmp);
            rY1tmp = round(Y1tmp); 
                     
%             
%             imshow(imread(stackname,'tiff',Slist(ii)))
%             hold on
%             plot([X2tmp-3 X2tmp-3 X2tmp+3 X2tmp+3 X2tmp-3],[Y2tmp-Leng-2 Y2tmp+Leng+2 Y2tmp+Leng+2 Y2tmp-Leng-2 Y2tmp-Leng-2],'y-*')

            ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-Lengbig-2 rY1tmp+Lengbig+2],[rX1tmp-3 rX1tmp+3]})'; % ROI verticale
            
            if length(ROI1) < 2*Lengbig+5
                    ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-3 rY1tmp+3],[rX1tmp-Lengbig-2 rX1tmp+Lengbig+2]}); % ROI horizontale
            end
            
            
            ROI1d = double(ROI1);
            
           
            X1inRoi = round((X1tmp - rX1tmp +Lengbig +2 +1)*10);
            Y1inRoi = round((Y1tmp - rY1tmp +3 +1)*10); 
            
        
            [X1roi,Y1roi] = meshgrid(rX1tmp-Lengbig-2:rX1tmp+Lengbig+2,rY1tmp-3:rY1tmp+3);
            [X1roi2,Y1roi2] = meshgrid(rX1tmp-Lengbig-2:0.1:rX1tmp+Lengbig+2,rY1tmp-3:0.1:rY1tmp+3);


            ROI1int = interp2(X1roi,Y1roi,ROI1d,X1roi2,Y1roi2,'spline');              
            L1 =  ROI1int(Y1inRoi,X1inRoi-Lengbig*10:X1inRoi+Lengbig*10);          

          
            [Ml1,~] = max(L1(0.5*l2big:1.5*l2big));
            [ml1,~] = min(L1(0.5*l2big:1.5*l2big));
            
            
            locsL1 = 1:length(L1);

            
            L1mask = (L1 > 0.5*(Ml1-ml1)+ml1)&(locsL1>l2big/2)&(locsL1<1.5*l2big);
            
            midL1 = round(mean(locsL1(L1mask).*(L1(L1mask)/mean(L1(L1mask)))));
                        

            if midL1+lbig>2*l2big+1
                midL1 = 2*l2big+1-lbig;
%                 cprintf('*err','midL shift !\n\n')                
               
            elseif midL1-lbig<1
                midL1 = lbig +1;                
%                 cprintf('*err','midL shift !\n\n')
            end
            
            try
                L1 = L1(midL1-lbig:midL1+lbig);
            catch
                
                cprintf('*err','Error with midL !!!\n\n')
                pause
            end
            
            
            L1 = (L1 + L1(end:-1:1))/2;

            
            [~,MPtmp,MVtmp] = ComputeZCorrel(L1,Ibig);        
              
            MultiPos1(ii,:) = MPtmp;
            MultiVal1(ii,:) = MVtmp;

            clear MVtmp MPtmp
           
        end
        


%% Petite bille

        
        MultiPos2 = zeros(3,2);
        MultiVal2 = zeros(3,2);
        
        for ii = 1:3
                 
            ptrpos = find(Stot == Slist(ii));
            
            if isempty(ptrpos)
                X2tmp = Xsmall(i);
                Y2tmp = Ysmall(i);
            else
                X2tmp = Xsmalltot(ptrpos);
                Y2tmp = Ysmalltot(ptrpos);
            end
            
            rX2tmp = round(X2tmp);
            rY2tmp = round(Y2tmp);            
%             
%             imshow(imread(stackname,'tiff',Slist(ii)))
%             hold on
%             plot([X2tmp-3 X2tmp-3 X2tmp+3 X2tmp+3 X2tmp-3],[Y2tmp-Leng-2 Y2tmp+Leng+2 Y2tmp+Leng+2 Y2tmp-Leng-2 Y2tmp-Leng-2],'y-*')

   
            
             ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-Lengsmall-2 rY2tmp+Lengsmall+2],[rX2tmp-3 rX2tmp+3]})';% ROI verticale
                        
             if length(ROI2) < 2*Lengsmall+5
                     ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-3 rY2tmp+3],[rX2tmp-Lengsmall-2 rX2tmp+Lengsmall+2]});% ROI horizontale
             end

            
             
            
            ROI2d = double(ROI2);
            
            
            X2inRoi = round((X2tmp - rX2tmp +Lengsmall +2 +1)*10);
            Y2inRoi = round((Y2tmp - rY2tmp +3 +1)*10); 
            
            

            [X2roi,Y2roi] = meshgrid(rX2tmp-Lengsmall-2:rX2tmp+Lengsmall+2,rY2tmp-3:rY2tmp+3);
            [X2roi2,Y2roi2] = meshgrid(rX2tmp-Lengsmall-2:0.1:rX2tmp+Lengsmall+2,rY2tmp-3:0.1:rY2tmp+3);

try
   
            ROI2int = interp2(X2roi,Y2roi,ROI2d,X2roi2,Y2roi2,'spline');      

catch
    themall
end
            L2 =  ROI2int(Y2inRoi,X2inRoi-Lengsmall*10:X2inRoi+Lengsmall*10);          

            
            
            
            [Ml2,~] = max(L2(0.5*l2small:1.5*l2small));
            [ml2,~] = min(L2(0.5*l2small:1.5*l2small));
            
  
            
            locsL2 = 1:length(L2);
            
 
            L2mask = (L2 > 0.5*(Ml2-ml2)+ml2)&(locsL2>l2small/2)&(locsL2<1.5*l2small);
            

            midL2 = round(mean(locsL2(L2mask).*(L2(L2mask)/mean(L2(L2mask)))));
            
            

            
            %bille 2
             if midL2+lsmall>2*l2small+1
                midL2 = 2*l2small+1-lsmall;
%                 cprintf('*err','midL shift !\n\n')                
               
            elseif midL2-lsmall<1
                midL2 = lsmall +1;                
%                 cprintf('*err','midL shift !\n\n')
            end
            
            try
                L2 = L2(midL2-lsmall:midL2+lsmall);
            catch
                
                cprintf('*err','Error with midL !!!\n\n')
                pause
            end
 
            L2 = (L2 + L2(end:-1:1))/2;
            
            % bille 2
            [~,MPtmp,MVtmp] = ComputeZCorrel(L2,Ismall);        
              
            MultiPos2(ii,:) = MPtmp;
            MultiVal2(ii,:) = MVtmp;
            
             clear MVtmp MPtmp
        end
        
        Sdone = [Sdone Slist];

        
        %% calcul dz
        
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
    
    
    [Zup1,Wup1,Zmid1,Wmid1,Zdown1,Wdown1] = Corr2Zs(Bool1,Str,MultiPos1,MultiVal1,ValSum1,fBig);
    [Zup2,Wup2,Zmid2,Wmid2,Zdown2,Wdown2] = Corr2Zs(Bool2,Str,MultiPos2,MultiVal2,ValSum2,fSmall);

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