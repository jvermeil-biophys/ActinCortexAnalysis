% fuse depthograph together
clear all

path = 'D:\MagneticPincherData\Raw\EtalonnageZ\MultiZCorrection';
savepath = 'D:\MagneticPincherData\Raw\EtalonnageZ\';

date = '21-07-08';
%genericName = [date '_Deptho_M2'];
genericName = [date '_Deptho270_M2'];

folderContent = dir('D:\MagneticPincherData\Raw\EtalonnageZ\MultiZCorrection');
nFiles = length(folderContent);

ndeptho = 0;

for i=1:nFiles
    currentFileName = folderContent(i).name;
    ext = '.tif';
    A = contains(currentFileName, genericName);
    B = contains(currentFileName, ext);
    if contains(currentFileName, genericName) && contains(currentFileName, ext)
        ndeptho = ndeptho+1;
    end
end

num = cell(1, ndeptho);
for i=1:ndeptho
    num{i} = num2str(i);
end

for kdepth = 1:ndeptho
    test = [path filesep genericName '_' num{kdepth} '.mat'];
    load([path filesep genericName '_' num{kdepth} '.mat'])
    Kim{kdepth} = K;
    Foc(kdepth) = f;
    LengK(kdepth) = size(K,1);

end

LK = unique(LengK);
if length(LK) ~= 1 
    error('Different kimo sizes')
end

fmin = min(Foc);
fmax = max(Foc);

for kdepth = 1:ndeptho
    fcur = Foc(kdepth);
    Kim{kdepth} = Kim{kdepth}(fcur-fmin+1:fcur+LK-fmax,:);


end

Ktot = zeros(size(Kim {1}));

for kdepth = 1:ndeptho

    Ktot = Ktot + double(Kim{kdepth});

end

Ktot = Ktot/ndeptho;

clear K f

K = uint16(Ktot);
f = fmin;

figure
imshow(imadjust(uint16(Ktot)))

save([savepath filesep genericName '.mat'],'K','f','stepnm');

imwrite(uint16(K),[savepath filesep genericName '.tif'])

clear K* f*
    

