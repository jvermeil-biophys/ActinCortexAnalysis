function dzcalc = DzCalc_Zcorr(X1,Y1,X2,Y2,S,stackname,depthoname,depthofolder)
% renvoi dz en micron 
tic

set(0,'DefaultFigureWindowStyle','docked')



Ls = length(S); % nombre de positions
% l = 450; % longeur de la ligne a correler
l = 300; % longeur de la ligne a correler
l2 = 500;
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





% on parcours les images
for i = 1:Ls
   
        
        MultiPos1 = zeros(3,2);
        MultiPos2 = zeros(3,2);
        MultiVal1 = zeros(3,2);
        MultiVal2 = zeros(3,2);
        
  
                 
  
                X1tmp = X1(i);
                Y1tmp = Y1(i);
                
                X2tmp = X2(i);
                Y2tmp = Y2(i);

            
            rX1tmp = round(X1tmp);
            rY1tmp = round(Y1tmp); 
            
            rX2tmp = round(X2tmp);
            rY2tmp = round(Y2tmp);            
            
            

            ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-3 rY1tmp+3],[rX1tmp-Leng-2 rX1tmp+Leng+2]}); % ROI verticale
            
            if length(ROI1) < 2*Leng+5
                    ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % ROI horizontale
            end
            
            ROI2 = imread(stackname,'tiff',S(i),'pixelregion',{[rY2tmp-3 rY2tmp+3],[rX2tmp-Leng-2 rX2tmp+Leng+2]}); % ROI verticale
                        
             if length(ROI2) < 2*Leng+5
                    ROI2 = imread(stackname,'tiff',S(i),'pixelregion',{[rY2tmp-Leng-2 rY2tmp+Leng+2],[rX2tmp-3 rX2tmp+3]})'; % ROI horizontale
             end
             if length(ROI2) < 2*Leng+5
                ROI2 = double(ROI2);
                bonusLength = 2*Leng+5 - size(ROI2,2);
                bonus = ROI2(:,end) * ones(1,bonusLength);
                ROI2 = horzcat(ROI2,bonus);
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

            try
                ROI1int = interp2(X1roi,Y1roi,ROI1d,X1roi2,Y1roi2,'spline');              
                L1 =  ROI1int(Y1inRoi,X1inRoi-Leng*10:X1inRoi+Leng*10);          
            catch
                themall
            end
            try
                ROI2int = interp2(X2roi,Y2roi,ROI2d,X2roi2,Y2roi2,'spline');              
                L2 =  ROI2int(Y2inRoi,X2inRoi-Leng*10:X2inRoi+Leng*10);          
            catch
                themall
            end
            
            
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
                pause
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
 
            %bille 1
            [~,MP1,MV1] = ComputeZCorrel(L1,I);        
              

            % bille 2
            [~,MP2,MV2] = ComputeZCorrel(L2,I);        
              
        
            Z1 = f-MP1(1);
            Z2 = f-MP2(1);
  dzcalc(i) = Z1-Z2;
    
end

dzcalc = -dzcalc*stepnm;


toc