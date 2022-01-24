function dzcalc = DzCalc_Zcorr_multiBeads(Xbig,Ybig,...
                Xsmall,Ysmall,S,stackname,depthonameBig,depthonameSmall,depthofolder);
            
      
            
            
% renvoi dz en nm
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
 

%% Grosse bille

       

        
        MultiPos1 = zeros(3,2);
        MultiVal1 = zeros(3,2);
        


                X1tmp = Xbig(i);

                Y1tmp = Ybig(i);


            
            rX1tmp = round(X1tmp);
            rY1tmp = round(Y1tmp); 
                     
%             
%             imshow(imread(stackname,'tiff',Slist(ii)))
%             hold on
%             plot([X2tmp-3 X2tmp-3 X2tmp+3 X2tmp+3 X2tmp-3],[Y2tmp-Leng-2 Y2tmp+Leng+2 Y2tmp+Leng+2 Y2tmp-Leng-2 Y2tmp-Leng-2],'y-*')

            ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-Lengbig-2 rY1tmp+Lengbig+2],[rX1tmp-3 rX1tmp+3]})'; % ROI verticale
            
            if length(ROI1) < 2*Lengbig+5
                    ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-3 rY1tmp+3],[rX1tmp-Lengbig-2 rX1tmp+Lengbig+2]}); % ROI horizontale
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

            
            [~,MP1,~] = ComputeZCorrel(L1,Ibig);        
              
  
            zbig = fBig - MP1(1);
        

        


%% Petite bille


                X2tmp = Xsmall(i);
                Y2tmp = Ysmall(i);

            
            rX2tmp = round(X2tmp);
            rY2tmp = round(Y2tmp);            
%             
%             imshow(imread(stackname,'tiff',S(i)))
%             hold on
%             plot([X2tmp-3 X2tmp-3 X2tmp+3 X2tmp+3 X2tmp-3],[Y2tmp-Leng-2 Y2tmp+Leng+2 Y2tmp+Leng+2 Y2tmp-Leng-2 Y2tmp-Leng-2],'y-*')

   
            
             ROI2 = imread(stackname,'tiff',S(i),'pixelregion',{[rY2tmp-Lengsmall-2 rY2tmp+Lengsmall+2],[rX2tmp-3 rX2tmp+3]})';% ROI verticale
                        
             if length(ROI2) < 2*Lengsmall+5
                     ROI2 = imread(stackname,'tiff',S(i),'pixelregion',{[rY2tmp-3 rY2tmp+3],[rX2tmp-Lengsmall-2 rX2tmp+Lengsmall+2]});% ROI horizontale
             end

            
             
            
            ROI2d = double(ROI2);
            
            
            X2inRoi = round((X2tmp - rX2tmp +Lengsmall +2 +1)*10);
            Y2inRoi = round((Y2tmp - rY2tmp +3 +1)*10); 
            
            

            [X2roi,Y2roi] = meshgrid(rX2tmp-Lengsmall-2:rX2tmp+Lengsmall+2,rY2tmp-3:rY2tmp+3);
            [X2roi2,Y2roi2] = meshgrid(rX2tmp-Lengsmall-2:0.1:rX2tmp+Lengsmall+2,rY2tmp-3:0.1:rY2tmp+3);


   
            ROI2int = interp2(X2roi,Y2roi,ROI2d,X2roi2,Y2roi2,'spline');      


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
            [~,MP2,~] = ComputeZCorrel(L2,Ismall);        
              
   
            zsmall = fSmall - MP2(1);


dzcalc(i) = (zsmall-zbig)*stepnm;

end



% toc