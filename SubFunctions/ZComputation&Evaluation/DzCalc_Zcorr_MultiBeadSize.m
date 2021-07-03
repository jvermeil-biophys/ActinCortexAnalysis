function dzcalc = DzCalc_Zcorr_MultiBeadSize(X1,Y1,X2,Y2,S,stackname,depthonameSB,depthonameLB,DIAMETERL,DIAMETERS,DELTA,depthofolder)

%
% Computes dz in micron by correlating a line of pixel from the beads image
% to a reference depthograph.
%
% dzcalc = DzCalc_Zcorr(X1,Y1,X2,Y2,S,stackname,depthoname,depthofolder)
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
l = 300; % target length (after interpolation) of pixel line to correlate /!\ Change to 100 for MyOne beads /!\
l2 = 1.1*l; % longer line to get extra points for a better interpolation
Leng = l2/10; % actual length of the pixel line that will be retrieved from the image

R = DIAMETERL/2;
r = DIAMETERS/2;

%% PART I

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


for i = 1:Ls % going through images
    
    % getting position of beads
    X1tmp = X1(i);
    Y1tmp = Y1(i);
    
    % rounding position because image is saved as a matrix without subpixel precision
    rX1tmp = round(X1tmp);
    rY1tmp = round(Y1tmp);

    
    % In order to get a pixel line to do the correlation on the depthograph,
    % we first select a rectangular ROI centered on the beads position. This
    % rectangle is  longer than the target line length, and has a width of
    % 5 pixels. There is then a 2D interpolation of this rectangle on a mesh
    % with 10x smaller mesh size. The final line taken in the interpolated
    % ROI to get the subpixel precision.
    
    
    % rectangular ROI around the first bead
    ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % vertical ROI; better because no interference from neighbouring beads
    
    ROI1d = double(ROI1); % from UINT16 image format to double
    
    if length(ROI1) < 2*Leng+5 % if the bead is to close to an edge an the ROI cannot be complete
        ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-3 rY1tmp+3],[rX1tmp-Leng-2 rX1tmp+Leng+2]}); % horizontal ROI
        
        ROI1d = double(ROI1); % from UINT16 image format to double
    end
    
    
    if length(ROI1) < 2*Leng+5 % if the beads happen to be in a corner of the image and the horizontal ROI is not complete either
        ROI1 = imread(stackname,'tiff',S(i),'pixelregion',{[rY1tmp-Leng-2 rY1tmp+Leng+2],[rX1tmp-3 rX1tmp+3]})'; % we get back to the vertical ROI
        ROI1d = double(ROI1); % from UINT16 image format to double
        
        % Interpolation to complete the missing pixel
        
        [Xmeshtmp,Ymeshtmp] =  meshgrid(rY1tmp-Leng-2:rY1tmp+Leng+2,rX1tmp-3:rX1tmp+3); % theoretical mesh of a complete ROI
        
        [x,y] = find(Xmeshtmp>0&Xmeshtmp<simg(2)&Ymeshtmp>0&Ymeshtmp<simg(1)); % Actual fraction of the ROI
        
        ROItmp = interp2(Xmeshtmp(unique(x),unique(y)),Ymeshtmp(unique(x),unique(y)),ROI1d,Xmeshtmp,Ymeshtmp,'spline'); % interpolation
        
        ROI1d = ROItmp;
        
    end

    
    % interpolation of the ROIs on a finer grid
    [X1roi,Y1roi] = meshgrid(rX1tmp-Leng-2:rX1tmp+Leng+2,rY1tmp-3:rY1tmp+3); % starting grid
    
    [X1roi2,Y1roi2] = meshgrid(rX1tmp-Leng-2:0.1:rX1tmp+Leng+2,rY1tmp-3:0.1:rY1tmp+3); % finer grid
    
    ROI1int = interp2(X1roi,Y1roi,ROI1d,X1roi2,Y1roi2,'spline'); % interpolation
    
    
    % getting back the subpixel position for X and Y in the interpolated rectangle with more resolution
    X1inRoi = round((X1tmp - rX1tmp +Leng +2 +1)*10);
    Y1inRoi = round((Y1tmp - rY1tmp +3 +1)*10);
    
    % getting pixel line for correlation
    L1 =  ROI1int(Y1inRoi,X1inRoi-Leng*10:X1inRoi+Leng*10);   % bead 1
    
    % 'symetrization' of the pixel line, to better correlate with
    % the kimograph, also symetric
    L1 = (L1(l/10+1:end-l/10) + L1(end-l/10:-1:l/10+1))/2;    
    
    % bead 1
    [~,MP1,~] = ComputeZCorrel(L1,I); % computing line correlation on kimograph

  
    
    
%% PART II


% loading depthograph
load([depthofolder filesep 'EtalonnageZ\' depthoname '.mat']);

% some old depthographs are saved with 'Ktot' and 'fmin' instaed of 'K' and 'f'
if exist('Ktot')
    K = Ktot;
end
if exist('fmin')
    f = fmin;
end

focalZ2 = f;

% normalization of depthograph
Lk = length(K);
I = K(:,ceil(Lk/2)-l:ceil(Lk/2)+l);
Id = double(I);
Ime = mean(Id,2);
RIme = double(repmat(Ime,1,2*l+1));
I = Id./RIme; % normalized depthograph


firstimg = imread(stackname,'tiff',1); % first image of the stack
simg = size(firstimg); % size in pixel of images in this stack


for i = 1:Ls % going through images
    
    % getting position of beads

    X2tmp = X2(i);
    Y2tmp = Y2(i);
    
    % rounding position because image is saved as a matrix without subpixel precision

    
    rX2tmp = round(X2tmp);
    rY2tmp = round(Y2tmp);
    
    % In order to get a pixel line to do the correlation on the depthograph,
    % we first select a rectangular ROI centered on the beads position. This
    % rectangle is  longer than the target line length, and has a width of
    % 5 pixels. There is then a 2D interpolation of this rectangle on a mesh
    % with 10x smaller mesh size. The final line taken in the interpolated
    % ROI to get the subpixel precision.

    
    
    
    % rectangular ROI around the second bead
    ROI2 = imread(stackname,'tiff',S(i),'pixelregion',{[rY2tmp-Leng-2 rY2tmp+Leng+2],[rX2tmp-3 rX2tmp+3]})'; % vertical ROI; better because no interference from neighbouring beads
    
    ROI2d = double(ROI2); % from UINT16 image format to double
    
    if length(ROI2) < 2*Leng+5 % if the bead is to close to an edge an the ROI cannot be complete
        ROI2 = imread(stackname,'tiff',S(i),'pixelregion',{[rY2tmp-3 rY2tmp+3],[rX2tmp-Leng-2 rX2tmp+Leng+2]}); % horizontal ROI
        
        ROI2d = double(ROI2); % from UINT16 image format to double
    end
    
    
    if length(ROI2) < 2*Leng+5 % if the beads happen to be in a corner of the image and the horizontal ROI is not complete either
        ROI2 = imread(stackname,'tiff',S(i),'pixelregion',{[rY2tmp-Leng-2 rY2tmp+Leng+2],[rX2tmp-3 rX2tmp+3]})'; % we get back to the vertical ROI
        ROI2d = double(ROI2); % from UINT16 image format to double
        
        % Interpolation to complete the missing pixel
        
        [Xmeshtmp,Ymeshtmp] =  meshgrid(rY2tmp-Leng-2:rY2tmp+Leng+2,rX2tmp-3:rX2tmp+3); % theoretical mesh of a complete ROI
        
        [x,y] = find(Xmeshtmp>0&Xmeshtmp<simg(2)&Ymeshtmp>0&Ymeshtmp<simg(1)); % Actual fraction of the ROI
        
        ROItmp = interp2(Xmeshtmp(unique(x),unique(y)),Ymeshtmp(unique(x),unique(y)),ROI2d,Xmeshtmp,Ymeshtmp,'spline'); % interpolation
        
        ROI2d = ROItmp;
        
    end
    
    % interpolation of the ROIs on a finer grid
    [X2roi,Y2roi] = meshgrid(rX2tmp-Leng-2:rX2tmp+Leng+2,rY2tmp-3:rY2tmp+3);

    [X2roi2,Y2roi2] = meshgrid(rX2tmp-Leng-2:0.1:rX2tmp+Leng+2,rY2tmp-3:0.1:rY2tmp+3);

    ROI2int = interp2(X2roi,Y2roi,ROI2d,X2roi2,Y2roi2,'spline');
    
    % getting back the subpixel position for X and Y in the interpolated rectangle with more resolution
   
    X2inRoi = round((X2tmp - rX2tmp +Leng +2 +1)*10);
    Y2inRoi = round((Y2tmp - rY2tmp +3 +1)*10);
    
    % getting pixel line for correlation
    L2 =  ROI2int(Y2inRoi,X2inRoi-Leng*10:X2inRoi+Leng*10);     % bead 2
    
    % 'symetrization' of the pixel line, to better correlate with
    % the kimograph, also symetric
    L2 = (L2(l/10+1:end-l/10) + L2(end-l/10:-1:l/10+1))/2;
    
    
    
    % bead 2
    [~,MP2,~] = ComputeZCorrel(L2,I); % computing line correlation on kimograph



%% PART III
    
    beadSizeMismatchOffset = (DELTA - (R-r))/stepnm;
    Z1 = f-MP1(1); % Z for first bead, computed as the position on the kimograph relative to the focus
    Z2 = f-MP2(1); % Z for second bead, computed as the position on the kimograph relative to the focus
    dzcalc(i) = Z2-Z1 - beadSizeMismatchOffset; % Z distance between the two beads in voxel
    
end

dzcalc = -dzcalc*stepnm; % cponversion in nanometers

end


