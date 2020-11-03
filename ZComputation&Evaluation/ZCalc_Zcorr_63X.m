function Z1 = ZCalc_Zcorr_63X(X1,Y1,S,stackname,depthoname,depthofolder)

% calcul de Z test pour 63X

set(0,'DefaultFigureWindowStyle','docked')



Ls = length(S); % nombre de positions
l = 300; % longeur de la ligne a correler
l2 = 320;
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
   
  
  
                X1tmp = X1(i);
                Y1tmp = Y1(i);

            
            rX1tmp = round(X1tmp);
            rY1tmp = round(Y1tmp);           
            
            

            ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-3 rY1tmp+3],[rX1tmp-Leng-2 rX1tmp+Leng+2]}); % ROI verticale
            
            if length(ROI1) < 2*Leng+4
                    ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % ROI horizontale
            end


            
             
            ROI1d = double(ROI1);
           
            
            X1inRoi = round((X1tmp - rX1tmp +Leng +2 +1)*10); % pour le positionnement du centre 
            Y1inRoi = round((Y1tmp - rY1tmp +3 +1)*10); % dans l'image interpolé x10

            
            
        [X1roi,Y1roi] = meshgrid(rX1tmp-Leng-2:rX1tmp+Leng+2,rY1tmp-3:rY1tmp+3);
            [X1roi2,Y1roi2] = meshgrid(rX1tmp-Leng-2:0.1:rX1tmp+Leng+2,rY1tmp-3:0.1:rY1tmp+3);


            ROI1int = interp2(X1roi,Y1roi,ROI1d,X1roi2,Y1roi2,'spline');              
            L1 =  ROI1int(Y1inRoi,X1inRoi-Leng*10:X1inRoi+Leng*10);          

            
            
            [Ml1,~] = max(L1(0.5*l2:1.5*l2));
            [ml1,~] = min(L1(0.5*l2:1.5*l2));
            
            locsL1 = 1:length(L1);
            
           
            
            L1mask = (L1 > 0.5*(Ml1-ml1)+ml1)&(locsL1>l2/2)&(locsL1<1.5*l2);
            
            midL1 = round(mean(locsL1(L1mask).*(L1(L1mask)/mean(L1(L1mask)))));
            
            
            
            
            
            %bille 1
            if midL1+l>2*l2+1
                midL1 = 2*l2+1-l;
                cprintf('*err','midL shift !\n\n')  
                figtitle = ['Image ' num2str(S(i)) ' - MidLShift'];              
               
            elseif midL1-l<1
                midL1 = l +1;                
                cprintf('*err','midL shift !\n\n')
                figtitle = ['Image ' num2str(S(i)) ' - MidLShift'];
            else
                figtitle = ['Image ' num2str(S(i))];
            end
            
            try
                L1cut = L1(midL1-l:midL1+l);
            catch
                
                cprintf('*err','Error with midL !!!\n\n')
                pause
            end
            
            
%             figure
%             hold on
%             subplot(311)
%             hold on
%             imshow(ROI1)
%             plot(X1inRoi/10,Y1inRoi/10,'*')
%             title(figtitle)
%             subplot(312)
%             hold on
%             imshow(uint16(ROI1int))
%             plot(X1inRoi,Y1inRoi,'*')
%             plot(X1inRoi-Leng*10:X1inRoi+Leng*10,ones(20*Leng+1)*Y1inRoi,'r-.')
%             subplot(313)
%             hold on
%             plot(L1)
%             plot(midL1,L1(midL1),'*')
%             

            %bille 1

            [~,MP1,MV1] = ComputeZCorrel(L1cut,I);        
              

           
              
            
        
            Z1(i) = (f-MP1(1))*stepnm;
            


    
end

end

