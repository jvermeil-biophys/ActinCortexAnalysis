function [K,f,stepnm] = MakeDepthoMy1(STD,X,Y,S,stackname,savename,stepnm)

% step en nm entre les images

% on allonge pour avoir les imager non traquées


Iinfo = imfinfo(stackname,'tiff');

Lstack = length(Iinfo);

Mstd = (STD>0.5*max(STD));

Xgood = mean(X(Mstd));
Ygood = mean(Y(Mstd));


Leng = 11; % longeur de la ligne a correler 55 pour 100X 33 pour 63X
l2 = Leng*10;
l = 100; % 5
K = zeros(Lstack,2*l+1);

for i = 1:Lstack
    
    I = imread(stackname,'tiff',i);
    
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
    
    XinRoi = round((Xtmp - rXtmp +Leng +2 +1)*10);
    YinRoi = round((Ytmp - rYtmp +3 +1)*10);
    
    try
        ROI = I(rYtmp-3:rYtmp+3,rXtmp-Leng-2:rXtmp+Leng+2);
    catch
        a= 1
    end
    
    %%%%%% affichage de la selection pour les images selectionnée
    if ismember(i,[])
    figure
    imshow(I)
    hold on
%     plot(X,Y,'*-')
    plot(Xtmp,Ytmp,'*')
    rectangle('Position',[rXtmp-Leng-2 rYtmp-3  2*Leng+4 6],'edgecolor','r')
    end
    
    
    ROId = double(ROI);
    
    [Xroi,Yroi] = meshgrid(rXtmp-Leng-2:rXtmp+Leng+2,rYtmp-3:rYtmp+3);
    [Xroi2,Yroi2] = meshgrid(rXtmp-Leng-2:0.1:rXtmp+Leng+2,rYtmp-3:0.1:rYtmp+3);
    
    ROIint = interp2(Xroi,Yroi,ROId,Xroi2,Yroi2,'spline');
    
    u16ROIint = uint16(ROIint);
    
    L =  u16ROIint(YinRoi,XinRoi-Leng*10:XinRoi+Leng*10);
    
    L = double(L);
    
    [Ml,iMl] = max(L(0.5*l2:1.5*l2));
    [ml,iml] = min(L(0.5*l2:1.5*l2));
    
    
    locsL = 1:length(L);
    Lmask = (L > 0.5*(Ml-ml)+ml)&(locsL>0.5*l2)&(locsL<1.5*l2);
    
    
    midL = round(mean(locsL(Lmask).*(L(Lmask)/mean(L(Lmask)))));
    
    
    if midL+l>2*l2+1
        midL = 2*l2+1-l;
    elseif midL-l<1
        midL = l +1;
    elseif isnan(midL)
        midL = 701;
    end
    

        L = L(midL-l:midL+l);

    
    L = (L + L(end:-1:1))/2;
    
    K(i,:) = L;
    
    
end


Lk = K(:,l);
Lks = smooth(Lk,'rlowess');

[~,f] = max(Lks);

figure
imshow(imadjust(uint16(K)));
hold on
plot(l,f,'r.')

figure
hold on
plot(Lks,'k')
plot(Lk,'m')
plot(f,Lks(f),'*r')

save([savename '.mat'],'K','f','stepnm');

imwrite(uint16(K),[savename '.tif'],'tiff');
