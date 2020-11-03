function CellHairAnalysis(path,savepath,figurepath,stackname,savename,PLOT)
% chemin vers la stack, nom de la stack (sans .tif),threshold pour
% binarisation,

warning('off','all')
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');


% dossier de sauvegarde
sf = savepath;
ff = [figurepath filesep datestr(now,'yy-mm-dd') filesep 'EdgesPlot'];
ff2 = [figurepath filesep datestr(now,'yy-mm-dd') filesep 'MaxLine'];
df = ['C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Figure' filesep 'FluoHairs'];
mkdir(sf)
mkdir(df)
mkdir(ff)
mkdir(ff2)

n = length(stackname);

    ElipsePerim = [];
    IntVar = [];
    IntPeaksSize = [];
   

for kn = 1:n

    


ptr = Inf;

% récupération du nombre d'images
Iinfo = imfinfo([path filesep stackname{kn} '.tif']);
stacksize = length(Iinfo)-1;

Res = Iinfo(1).XResolution;



Img= imread([path filesep stackname{kn} '.tif'],'tiff',1);
% Img1 = imgaussfilt(Img1);

 % Binarization
    ImgGauss = imgaussfilt(Img);
    
    Th = adaptthresh(ImgGauss,0.6,'NeighborhoodSize',21,'Statistic','gaussian');
    ImgBin = imbinarize(ImgGauss,Th);  
    ImgBin = imfill(ImgBin,'holes');
    
    ImgEro = imerode(ImgBin,strel('square',4));
    ImgDil = imdilate(ImgEro,strel('square',4));
    ImgDil = imdilate(ImgDil,strel('square',4));
    
    ImgBack = uint8(double(Img).*~ImgDil);
   
%     imwrite(ImgBack,[ff2 filesep savename '_' stackname{kn} '_Backgroung.tif'],'TIFF')
% 
%     
    
    CellBin = bwpropfilt(ImgBin,'Area',1,4); 
    CellBin = imdilate(CellBin,strel('sphere',2));
    
    ImgCropped = uint8(CellBin.*double(Img));

    CellEdge = edge(CellBin);
    
    [x,y] = find(CellEdge == 1); 
    
 if PLOT 
    figure
%     histogram(Img1)
    hold on
    imshow(imadjust(ImgGauss))
    plot(y,x,'.')
    
    legend('Perimeter')
    fig = gcf;
    saveas(fig,[ff filesep savename '_' stackname{kn} '.png'])
    close
 end
    
     stats = regionprops(CellBin,'Perimeter','Area','MajorAxisLength','MinorAxisLength','Solidity');

     

    % fited elipse perimeter, area and circularity
    a = stats.MajorAxisLength/2;
    b = stats.MinorAxisLength/2;    
    h = (a-b)^2/(a+b)^2;
    perimE = pi*(a+b)*(1+ 3*h/(10+sqrt(4-3*h)));  
    

    ElipsePerim = [ElipsePerim stats.Perimeter/perimE];
    
    
    
 % variation d'intensité
[Xc,Yc,UfImg] = SingleCellFluoUnfolding(ImgCropped);

[~,rMax] = max(double(UfImg));

Lim = size(ImgCropped,1);


xMax = Yc-(size(UfImg,1)-rMax+49).*cos(deg2rad(1:361));
yMax = Xc-(size(UfImg,1)-rMax+49).*sin(deg2rad(1:361));

rdxMax = round(xMax);
rdyMax = round(yMax);


for k = 1:361

MaxLine(k) = mean(mean(ImgCropped(rdxMax(k):rdxMax(k)+1,rdyMax(k):rdyMax(k)+1)));

end
% 
MaxLine = MaxLine-mean(mean(double(ImgBack)));
MaxLine = MaxLine/mean(MaxLine);

MaxLineS = smooth(MaxLine,7,'sgolay',3);

[Pks,PosPks,wPks,pPks] = findpeaks(MaxLineS,'MinPeakProminence',0.4*mean(MaxLine));


 if PLOT 
    figure
    hold on
    subplot(121)
    hold on
    imshow(imadjust(ImgCropped))
    plot(yMax(PosPks),xMax(PosPks),'ro','markersize',8)
    plot(yMax,xMax,'b.','markersize',3)
    plot(yMax(1),xMax(1),'m<','markerfacecolor','m')
    title([savename ':' stackname{kn}])
    subplot(222)
    plot(MaxLine)
    title('Original Data')
    set(gca,'FontSize',12);
    subplot(224)
    findpeaks(MaxLineS,'MinPeakProminence',0.4*mean(MaxLine),'annotate','extent')
    title('Smoothed data for peaks detection')
    set(gca,'FontSize',12);
    
    fig = gcf;
    saveas(fig,[ff2 filesep savename '_' stackname{kn} '_MaxLine''.png'])
    close
 end


% 
% for k = 2:length(MaxLine)-1
%     try
%     submat = double(UfImg(yM(k)-1:yM(k)+1,k-1:k+1));
%     catch
%         jdhegfgge
%     end
%     AvgMaxLine(k-1) = mean(submat(:));
% end
% 
% 
% AvgMaxLine = AvgMaxLine/mean(AvgMaxLine);


% IntVar = [IntVar std(AvgMaxLine)];


IntVar = [IntVar mean(mean(double(ImgBack)))];
IntPeaksSize = [IntPeaksSize pPks'];
clear *Pks
    
  
MR.Names{kn} = [savename '_' stackname{kn}];

   
    
end

    MR.ElipsePerim = ElipsePerim;
    MR.IntVar = IntVar;
    MR.IntPeaksSize = IntPeaksSize;
    
save([sf filesep savename '.mat'],'MR')

end
