function [Xc,Yc,ImgStack] = findcenter(path,stackname,stacksize,th)

ImgStack = cell(1,stacksize);

for ks = 1:stacksize
    
    % récupération image
    Itmp = imread([path filesep stackname '.tif'],'tiff',ks);
    Itmp = imgaussfilt(Itmp);
    
    clear Img imgBin Spots CellMask
    
    % Binarization
    Img = double(Itmp);
    ImgBin = Img > th;
    ImgBin = imfill(ImgBin,'holes');
    ImgBin = imerode(ImgBin,strel('square',3));
    
    % detections des objets sur l'image binaire
    Spots = bwlabel(ImgBin);    
    nspots = max(max(Spots));
    SpotSizes = zeros(nspots,1);
    for ksp = 1:nspots
        SpotSizes(ksp) = sum(sum(Spots == ksp));
    end
    
    [~,CellInd] = max(SpotSizes); % Label du plus gros objet, la cellule
    
    try
    Spots(Spots ~= CellInd) = 0;
    catch
        lol
    end
    
    CellMask = logical(imfill(imdilate(Spots,strel('square',3)),'holes'));
    
    ImgStack{ks} = Itmp.*uint8(CellMask);
    
%     imshow(ImgStack{ks})
    
    % utilisation de la distance map pour trouver le centre
    D = bwdist(~CellMask);
    [Y,X] = find(D == max(max(D)));
    
    Xc(ks) = round(mean(X));
    Yc(ks) = round(mean(Y));
    
        figure(1)
        imshow(ImgStack{ks})
        pause
        hold off
    
    figure(1)
    imshow(CellMask)
    pause
    hold off
    
    figure(1)
    imshow(uint8(D))
    hold on
    plot(Xc,Yc,'*')    
    pause
    hold off
    
end



Xc = smooth(Xc);
Yc = smooth(Yc);


end