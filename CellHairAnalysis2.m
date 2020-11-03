function CellHairAnalysis2(path,savepath,figurepath,stackname,savename,PLOT)
% chemin vers la stack, nom de la stack (sans .tif),threshold pour
% binarisation,

warning('off','all')
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');

% dossier de sauvegarde
sf = savepath;
ff = [figurepath filesep datestr(now,'yy-mm-dd') filesep 'EdgesPlot'];
ff2 = [figurepath filesep datestr(now,'yy-mm-dd') filesep 'EdgesPlotEroded'];
df = ['C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Figure' filesep 'FluoHairs'];
mkdir(sf)
mkdir(df)
mkdir(ff)
mkdir(ff2)

n = length(stackname);

    ElipsePerim = [];
    Solidity = [];
    PerimRatio = [];

for kn = 1:n


Img= imread([path filesep stackname{kn} '.tif'],'tiff',1);

%% fitting elipse

    % Binarization
    ImgGauss = imgaussfilt(Img);
    
    Th = adaptthresh(ImgGauss,0.6,'NeighborhoodSize',21,'Statistic','gaussian');
    ImgBin = imbinarize(ImgGauss,Th);   
    ImgBin = imfill(ImgBin,'holes');
    
    CellBin = bwpropfilt(ImgBin,'Area',1,4); 
    
%     figure
%     imshow(CellBin)
    
    CellEdge = edge(CellBin);
    
    [x,y] = find(CellEdge == 1); 
    
    figure
    imshow(imadjust(ImgGauss))
    hold on
    plot(y,x,'.');
    
    
    ImgBin = ImgDouble > mode(Img((Img>0)));
    ImgBin = imfill(ImgBin,'holes');

    CellBin = bwpropfilt(ImgBin,'perimeter',1,4);

     stats = regionprops(CellBin,'Perimeter','Area','MajorAxisLength','MinorAxisLength','Solidity');

    % fited elipse perimeter, area and circularity
    a = stats.MajorAxisLength/2;
    b = stats.MinorAxisLength/2;    
    h = (a-b)^2/(a+b)^2;
    perimE = pi*(a+b)*(1+ 3*h/(10+sqrt(4-3*h)));  
    
%% getting real perimeter through unfolding

ImgCropped = uint8(CellBin.*double(Img));

[Xc,Yc,UfImg] = SingleCellFluoUnfolding(ImgCropped);

MaxLine = max(UfImg)/2;

MaxMat = repmat(MaxLine,size(UfImg,1),1);

ImgDiff = UfImg-MaxMat;

ImgDiffBin = ImgDiff >0;

CellCortBin = bwpropfilt(ImgDiffBin,'area',[200 Inf],4);
CellCortBin = imfill(CellCortBin,'holes');

imshow(CellCortBin)

for k = 1:361
    CellEdge(k) = find(ImgDiffBin(:,k),1);
end

% 
% figure
% imshow(imadjust(UfImg))
% hold on 
% plot(1:361,CellEdge,'r')

theta = deg2rad(1:361);
r = size(UfImg,1)-CellEdge+24;

X = Yc- r.*cos(theta);
Y = Xc -r.*sin(theta);

% figure
% imshow(imadjust(Img))
% hold on
% plot(Y,X,'r.-','markersize',7)
% plot(Xc,Yc,'*')

    ElipsePerim = [ElipsePerim stats.Perimeter/perimE];
    Solidity = [Solidity stats.Solidity];
    
    
    
    
   
     if PLOT &0
    figure
%     histogram(Img1)
    hold on
    imshow(imadjust(imgaussfilt(Img)))
    plot(y,x,'.')
    
    legend('Perimeter','ErodedPerimeter')
    fig = gcf;
    saveas(fig,[ff filesep savename '_' stackname{kn} '_Img' num2str(ki) '.png'])
    close
 end
    
    
end
    
  
    MR.Names{kn} = [savename '_' stackname{kn}];

   
    


    MR.Solidity = Solidity;
    MR.PerimRatio = PerimRatio;
    MR.ElipsePerim = ElipsePerim;
    
save([sf filesep savename '.mat'],'MR')

end
