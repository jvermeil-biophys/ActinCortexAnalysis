function whosWho = AsymmetricBeadPair_findWhosWho(X1,Y1,X2,Y2,S,stackpath,diameterL,diameterS)

H = ceil(diameterL/2);
W = ceil(diameterL/2);

ROI1 = imread(stackpath,'tiff',S,'pixelregion',{[ceil(Y1)-H ceil(Y1)+H],[ceil(X1)-W ceil(X1)+W]});
% figure(1)
% colormap(gray(65535));
% I1 = imadjust(ROI1);
% image(I1)

[counts1,binLocations1] = imhist(ROI1);
counts1bis = smooth(counts1,3);
max1 = max(counts1bis);
% figure(2)
% plot(binLocations1, counts1bis)

ROI2 = imread(stackpath,'tiff',S,'pixelregion',{[ceil(Y2)-H ceil(Y2)+H],[ceil(X2)-W ceil(X2)+W]});
% figure(3)
% colormap(gray(65535));
% I2 = imadjust(ROI2);
% image(I2)

[counts2,binLocations2] = imhist(ROI2);
counts2bis = smooth(counts2,3);
max2 = max(counts2bis);
% figure(4)
% plot(binLocations2, counts2bis)

whosWho = (max1 > max2);
% The height of the highest peak correspond to the amount of black pixels in the ROI, which will be higher for the
% largest bead.