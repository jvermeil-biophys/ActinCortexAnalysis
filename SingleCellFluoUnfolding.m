function [Xc,Yc,ImgUf] = SingleCellFluoUnfolding(Img)
% chemin vers la stack, nom de la stack (sans .tif),threshold pour
% binarisation,

warning('off','all')
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');


%% repérage du centre

% Binarization
    ImgD = double(Img);
    ImgBin = ImgD > mode(ImgD((ImgD>0)))/2;
    ImgBin = imfill(ImgBin,'holes');
    ImgBin = imerode(ImgBin,strel('square',3.));
    
    CellBin = bwpropfilt(ImgBin,'perimeter',1,4);
    
    % utilisation de la distance map pour trouver le centre
    D = bwdist(~CellBin);
    [Y,X] = find(D == max(max(D)));
    
    Xc = round(mean(X));
    Yc = round(mean(Y));
    
%% dépliage
    ImgUf = RadialUnfolding(Img,Xc,Yc);

end
