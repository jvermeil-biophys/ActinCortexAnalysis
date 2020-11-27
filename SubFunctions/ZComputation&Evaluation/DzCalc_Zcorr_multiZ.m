function dzcalc = DzCalc_Zcorr_multiZ(X1tot,X1,Y1tot,Y1,X2tot,...
    X2,Y2tot,Y2,Stot,S,stackname,Sdown,Smid,Sup,depthoname,depthofolder)
%
% Computes dz in micron by correlating a line of pixel from the beads image
% to a reference depthograph. This version uses Z triplet of images taken
% with a piezo objective in order to be more precise. On each image of the
% triplet, a pixel line is correlated to the depthograph, and then dz is
% averaged on the three image of a triplet.
%
% dzcalc = DzCalc_Zcorr_multiZ(X1tot,X1,Y1tot,Y1,X2tot,...
%   X2,Y2tot,Y2,Stot,S,stackname,Sdown,Smid,Sup,depthoname,depthofolder)
%
% OUTPUTS :
% dzcalc : distane dz in microns between the two beads
%
% INPUTS :
% X1, Y1 : corrdinates of bead 1
% X2, Y2 : corrdinates of bead 2
% S : list of images to analyze
% stackname : name of the tiff stack to analyze
% depthoname : name of the depthograph to use
% depthofolder : folder where to find the depthograph
%
%   /!\ /!\ Works for M450 and M270, not for MyOne. You need to change line 26 in the code  /!\ /!\


Ls = length(S); % number of images
l = 300; % target length (after interpolation) of pixel line to correlate
l2 = 1.1*l; % longer line to get extra points for a better interpolation
Leng = l2/10; % actual length of the pixel line that will be retrieved from the image

% loading depthograph
load([depthofolder filesep 'EtalonnageZ\' depthoname '.mat']);

% some old depthographs are saved with 'Ktot' and 'fmin' instaed of 'K' and 'f'
if exist('Ktot')
    K = Ktot;
end
if exist('fmin')
    f = fmin;
end

% normalization of depthograph
Lk = length(K);
I = K(:,ceil(Lk/2)-l:ceil(Lk/2)+l);
Id = double(I);
Ime = mean(Id,2);
RIme = double(repmat(Ime,1,2*l+1));
I = Id./RIme; % normalized depthograph


firstimg = imread(stackname,'tiff',1); % first image of the stack
simg = size(firstimg); % size in pixel of images in this stack

Sdone = []; % list of analyzed images


for i = 1:Ls % going through images
    if not(ismember(S(i),Sdone)) % checking that image has not been analyzed yet
        
        % loading of the 3 images that constitute the Z triplet
        if ismember(S(i),Sdown)
            Slist = [S(i), S(i)+1, S(i)+2];
        elseif ismember(S(i),Smid)
            Slist = [S(i)-1, S(i), S(i)+1];
        elseif ismember(S(i),Sup)
            Slist = [S(i)-2, S(i)-1, S(i)];
        end
        
        % initiailzing variable for Z position and correlation values on
        % each of the three images
        MultiPos1 = zeros(3,2);
        MultiPos2 = zeros(3,2);
        MultiVal1 = zeros(3,2);
        MultiVal2 = zeros(3,2);
        
        for ii = 1:3 % going through the three images in the Z triplet
            
            % getting beads position
            ptrpos = find(Stot == Slist(ii)); % if the beads have been tracked on this image, we take their tracked positions
            
            if ~isempty(ptrpos)
                X1tmp = X1tot(ptrpos);
                Y1tmp = Y1tot(ptrpos);
                
                X2tmp = X2tot(ptrpos);
                Y2tmp = Y2tot(ptrpos);
                
            else % is the beads have not been tracked, we take the position on the 'best focus' image of the triplet
                X1tmp = X1(i);
                Y1tmp = Y1(i);
                
                X2tmp = X2(i);
                Y2tmp = Y2(i);
            end
            
            % rounding position because image is saved as a matrix without subpixel precision
            rX1tmp = round(X1tmp);
            rY1tmp = round(Y1tmp);
            
            rX2tmp = round(X2tmp);
            rY2tmp = round(Y2tmp);
            
            
            
            % In order to get a pixel line to do the correlation on the depthograph,
            % we first select a rectangular ROI centered on the beads position. This
            % rectangle is  longer than the target line length, and has a width of
            % 5 pixels. There is then a 2D interpolation of this rectangle on a mesh
            % with 10x smaller mesh size. The final line taken in the interpolated
            % ROI to get the subpixel precision.
            
            
            % rectangular ROI around the first bead
            ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % vertical ROI; better because no interference from neighbouring beads
            
            ROI1d = double(ROI1); % from UINT16 image format to double
            
            if length(ROI1) < 2*Leng+5 % if the bead is to close to an edge an the ROI cannot be complete
                ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-3 rY1tmp+3],[rX1tmp-Leng-2 rX1tmp+Leng+2]}); % horizontal ROI
                
                ROI1d = double(ROI1); % from UINT16 image format to double
            end
            
            
            if length(ROI1) < 2*Leng+5 % if the beads happen to be in a corner of the image and the horizontal ROI is not complete either
                ROI1 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % we get back to the vertical ROI
                ROI1d = double(ROI1); % from UINT16 image format to double
                
                % Interpolation to complete the missing pixel
                
                [Xmeshtmp,Ymeshtmp] =  meshgrid(rY1tmp-Leng-2:rY1tmp+Leng+2,rX1tmp-3:rX1tmp+3); % theoretical mesh of a complete ROI
                
                [x,y] = find(Xmeshtmp>0&Xmeshtmp<=simg(1)&Ymeshtmp>0&Ymeshtmp<=simg(2)); % Actual fraction of the ROI
                
 try
                ROItmp = interp2(Xmeshtmp(unique(x),unique(y)),Ymeshtmp(unique(x),unique(y)),ROI1d,Xmeshtmp,Ymeshtmp,'spline'); % interpolation
 catch
     themall
 end
                ROI1d = ROItmp;
                
            end
            
            
            
            % rectangular ROI around the second bead
            ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-Leng-2 rY2tmp+Leng+2],[rX2tmp-3 rX2tmp+3]})'; % vertical ROI; better because no interference from neighbouring beads
            
            ROI2d = double(ROI2); % from UINT16 image format to double
            
            if length(ROI2) < 2*Leng+5 % if the bead is to close to an edge an the ROI cannot be complete
                ROI2 = imread(stackname,'tiff',Slist(ii),'pixelregion',{[rY2tmp-3 rY2tmp+3],[rX2tmp-Leng-2 rX2tmp+Leng+2]}); % horizontal ROI
                
                ROI2d = double(ROI2); % from UINT16 image format to double
            end
            
            
            if length(ROI2) < 2*Leng+5 % if the beads happen to be in a corner of the image and the horizontal ROI is not complete either
                ROI2 = imread(stackname,'tiff',S(i),'pixelregion',{[rY2tmp-Leng-2 rY2tmp+Leng+2],[rX2tmp-3 rX2tmp+3]})'; % we get back to the vertical ROI
                ROI2d = double(ROI2); % from UINT16 image format to double
                
                % Interpolation to complete the missing pixel
                
                [Xmeshtmp,Ymeshtmp] =  meshgrid(rY2tmp-Leng-2:rY2tmp+Leng+2,rX2tmp-3:rX2tmp+3); % theoretical mesh of a complete ROI
                
                [x,y] = find(Xmeshtmp>0&Xmeshtmp<=simg(1)&Ymeshtmp>0&Ymeshtmp<=simg(2)); % Actual fraction of the ROI
                
                ROItmp = interp2(Xmeshtmp(unique(x),unique(y)),Ymeshtmp(unique(x),unique(y)),ROI2d,Xmeshtmp,Ymeshtmp,'spline'); % interpolation
                
                ROI2d = ROItmp;
                
            end
            
            % interpolation of the ROIs on a finer grid
            [X1roi,Y1roi] = meshgrid(rX1tmp-Leng-2:rX1tmp+Leng+2,rY1tmp-3:rY1tmp+3); % starting grid
            [X2roi,Y2roi] = meshgrid(rX2tmp-Leng-2:rX2tmp+Leng+2,rY2tmp-3:rY2tmp+3);
            
            [X1roi2,Y1roi2] = meshgrid(rX1tmp-Leng-2:0.1:rX1tmp+Leng+2,rY1tmp-3:0.1:rY1tmp+3); % finer grid
            [X2roi2,Y2roi2] = meshgrid(rX2tmp-Leng-2:0.1:rX2tmp+Leng+2,rY2tmp-3:0.1:rY2tmp+3);
            
            ROI1int = interp2(X1roi,Y1roi,ROI1d,X1roi2,Y1roi2,'spline'); % interpolation
            ROI2int = interp2(X2roi,Y2roi,ROI2d,X2roi2,Y2roi2,'spline');
            
            
            % getting back the subpixel position for X and Y in the interpolated rectangle with more resolution
            X1inRoi = round((X1tmp - rX1tmp +Leng +2 +1)*10);
            Y1inRoi = round((Y1tmp - rY1tmp +3 +1)*10);
            
            X2inRoi = round((X2tmp - rX2tmp +Leng +2 +1)*10);
            Y2inRoi = round((Y2tmp - rY2tmp +3 +1)*10);
            
            % getting pixel line for correlation
            L1 =  ROI1int(Y1inRoi,X1inRoi-Leng*10:X1inRoi+Leng*10);   % bead 1
            L2 =  ROI2int(Y2inRoi,X2inRoi-Leng*10:X2inRoi+Leng*10);     % bead 2
            
            % 'symetrization' of the pixel line, to better correlate with
            % the kimograph, also symetric.
            L1 = (L1(l/10+1:end-l/10) + L1(end-l/10:-1:l/10+1))/2;
            L2 = (L2(l/10+1:end-l/10) + L2(end-l/10:-1:l/10+1))/2;
            
            
            % bead 1
            [~,MPtmp,MVtmp] = ComputeZCorrel(L1,I); % computing line correlation on kimograph
            
            MultiPos1(ii,:) = MPtmp; % storing position results for image ii of the Z triplet
            MultiVal1(ii,:) = MVtmp; % storing correlation values for image ii of the Z triplet
            
            clear MVtmp MPtmp
            % bille 2
            [~,MPtmp,MVtmp] = ComputeZCorrel(L2,I); % computing line correlation on kimograph
            
            MultiPos2(ii,:) = MPtmp; % storing position results for image ii of the Z triplet
            MultiVal2(ii,:) = MVtmp; % storing correlation values for image ii of the Z triplet
            
            clear MVtmp MPtmp
        end
        
        Sdone = [Sdone Slist]; % Adding the three images of the triplet to the list of analyzed images
        

        % sorting correlations with alignement constraints 
        [Bool1,Bool2,ValSum1,ValSum2] = deal(NaN(1,8));
        Str = cell(1,8);
        
        for k1 = 1:2
            for k2 = 1:2
                for k3 = 1:2
                    k = bi2de([k3-1,k2-1,k1-1])+1;
                    Str{k} = [num2str(k1) num2str(k2) num2str(k3)]; % name of combination. '111' is first correlation value on each image of triplet. 
                    % '122' is first correlation of first image and second correlation on second and third image
                    
                    % for this specific combination of correlation, check if : 
                    Bool1(k) = (MultiPos1(1,k1)>MultiPos1(2,k2)) && ... % bottom Z is lower than middle Z
                        (MultiPos1(2,k2)>MultiPos1(3,k3))&& ... % middle Z is lower than top Z
                        (MultiPos1(1,k1)<35 + MultiPos1(2,k2)) &&... % not more than 700nm (35 voxel) between bottom and middle Z
                        (MultiPos1(2,k2)< 35 +MultiPos1(3,k3));  % not more than 700 nm (35 voxel) between middle and top Z
                    
                    ValSum1(k) = MultiVal1(1,k1)+MultiVal1(2,k2)+MultiVal1(3,k3); % sum of correlation values for this combination
                    
                    % for this specific combination of correlation, check if : 
                    Bool2(k) = (MultiPos2(1,k1)>MultiPos2(2,k2)) &&... % bottom Z is lower than middle Z
                        (MultiPos2(2,k2)>MultiPos2(3,k3))&& ... % middle Z is lower than top Z
                        (MultiPos2(1,k1)<35+MultiPos2(2,k2)) &&... % not more than 700nm (35 voxel) between bottom and middle Z
                        (MultiPos2(2,k2)<35+MultiPos2(3,k3)); % not more than 700 nm (35 voxel) between middle and top Z
                    
                    ValSum2(k) = MultiVal2(1,k1)+MultiVal2(2,k2)+MultiVal2(3,k3); % sum of correlation values for this combination
                end
            end
        end
        clear k1 k2 k3
        
        
        % getting best Zs with aligned correlation 
        
        [Zup1,Wup1,Zmid1,Wmid1,Zdown1,Wdown1] = Corr2Zs(Bool1,Str,MultiPos1,MultiVal1,ValSum1,f);
        [Zup2,Wup2,Zmid2,Wmid2,Zdown2,Wdown2] = Corr2Zs(Bool2,Str,MultiPos2,MultiVal2,ValSum2,f);
        
        % saving dz and corresponding sum of correlation values for each
        % image.
        dzup(i) = Zup2-Zup1;
        Wdzup(i) = Wup1+Wup2;
        dzmid(i) = Zmid2-Zmid1;
        Wdzmid(i) = Wmid1+Wmid2;
        dzdown(i) = Zdown2-Zdown1;
        Wdzdown(i) = Wdown1+Wdown2;
        
    end
end

dzcalc = sum([dzup.*Wdzup; dzmid.*Wdzmid; dzdown.*Wdzdown],1)./(sum([Wdzup; Wdzmid; Wdzdown],1))*stepnm; % computing dz on as the weighted average of dz on each Z triplet, with the correlation values as weight.

end