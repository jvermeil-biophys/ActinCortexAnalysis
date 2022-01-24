function [K,f,f_HD,stepnm] = MakeDeptho(STD,X,Y,S,stackname,savename,stepnm)

% step en nm entre les images

% on allonge pour avoir les imager non traquées


Iinfo = imfinfo(stackname,'tiff');

Lstack = length(Iinfo);

Mstd = (STD>0.5*max(STD));

Xgood = mean(X(Mstd));
Ygood = mean(Y(Mstd));


lengthInput = 55; % longeur de la ligne a correler 55 pour 100X 33 pour 63X
l2 = lengthInput*10;
halfWidth = 500; % 500 pour 100X 300 pour 63X
DepthoMatrix = zeros(Lstack,2*halfWidth+1);

for i = 1:Lstack
    
    currentImage = imread(stackname,'tiff',i);
    ptr = find(S==i);
    
    if length(ptr) == 1
        Xtmp = X(ptr);
        Ytmp = Y(ptr);
    else
        Xtmp = Xgood;
        Ytmp = Ygood;
    end
    
    rXtmp = round(Xtmp);
    rYtmp = round(Ytmp);
    
    XinRoi = round((Xtmp - rXtmp +lengthInput +2 +1)*10);
    YinRoi = round((Ytmp - rYtmp +3 +1)*10);
    
    try
        ROI = currentImage(rYtmp-3:rYtmp+3,rXtmp-lengthInput-2:rXtmp+lengthInput+2);
    catch
        a = 1
    end
    
    %%%%%% affichage de la selection pour les images selectionnée
    if ismember(i,[])
    figure
    imshow(currentImage)
    hold on
    % plot(X,Y,'*-')
    plot(Xtmp,Ytmp,'*')
    rectangle('Position',[rXtmp-lengthInput-2, rYtmp-3, 2*lengthInput+4, 6],'edgecolor','r')
    end
    
    
    ROId = double(ROI);
    
    [Xroi,Yroi] = meshgrid(rXtmp-lengthInput-2:rXtmp+lengthInput+2,rYtmp-3:rYtmp+3);
    [Xroi2,Yroi2] = meshgrid(rXtmp-lengthInput-2:0.1:rXtmp+lengthInput+2,rYtmp-3:0.1:rYtmp+3);
    
    ROIint = interp2(Xroi,Yroi,ROId,Xroi2,Yroi2,'spline');
    
    u16ROIint = uint16(ROIint);
    
    DepthoLineOutput =  u16ROIint(YinRoi,XinRoi-lengthInput*10:XinRoi+lengthInput*10);
    
    DepthoLineOutput = double(DepthoLineOutput);
    
    [Ml,iMl] = max(DepthoLineOutput(0.5*l2:1.5*l2));
    [ml,iml] = min(DepthoLineOutput(0.5*l2:1.5*l2));
    
    
    locsL = 1:length(DepthoLineOutput);
    Lmask = (DepthoLineOutput > 0.5*(Ml-ml)+ml)&(locsL>0.5*l2)&(locsL<1.5*l2);
    
    
    midL = round(mean(locsL(Lmask).*(DepthoLineOutput(Lmask)/mean(DepthoLineOutput(Lmask)))));
    
    
    if midL+halfWidth>2*l2+1
        midL = 2*l2+1-halfWidth;
    elseif midL-halfWidth<1
        midL = halfWidth +1;
    elseif isnan(midL)
        midL = 701;
    end
    

        DepthoLineOutput = DepthoLineOutput(midL-halfWidth:midL+halfWidth);

    
    DepthoLineOutput = (DepthoLineOutput + DepthoLineOutput(end:-1:1))/2;
    
    DepthoMatrix(i,:) = DepthoLineOutput;
    
    
end



Lk = DepthoMatrix(:,halfWidth);
Lks = smooth(Lk,'rlowess');
Lks_HD = interp(Lks,4);

% Experimental new code !
for i = 100:301
    if(Lk(i) < 0.9 * Lks(i))
        DepthoMatrix(i,:) = DepthoMatrix(i-1,:);
    end
end

[~,f] = max(Lks);
[~,f_HD] = max(Lks_HD);

figure
imshow(uint16(DepthoMatrix));
hold on
plot(halfWidth,f,'r.')

figure
hold on
plot(Lks,'k')
plot(Lk,'m')
plot(f,Lks(f),'*r')

K = DepthoMatrix;
save([savename '.mat'],'K','f','f_HD','stepnm');

imwrite(uint16(DepthoMatrix),[savename '.tif'],'tiff');
