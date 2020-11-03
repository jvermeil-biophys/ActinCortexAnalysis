function CellFluoUnfolding(path,savepath,stackname,th,savename)
% chemin vers la stack, nom de la stack (sans .tif),threshold pour
% binarisation,

warning('off','all')
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');

% dossier de sauvegarde
ff = savepath;
mkdir(ff)




%% depliage de l'image


if exist([path filesep stackname '_Fluo_Unfolded.tif'])
    delete([path filesep stackname '_Fluo_Unfolded.tif'])
end

ptr = Inf;

% récupération du nombre d'images
Iinfo = imfinfo([path filesep stackname '.tif']);
stacksize = length(Iinfo)-1;

Res = Iinfo(1).XResolution;

% réupération de la stack, selection de la cellules seulement et calcul des centres
[Xc,Yc,ImgStack] = findcenter(path,stackname,stacksize,th);

savenamepic = [ff filesep savename '_Pic.tif'];
savename1   = [ff filesep savename '_CleanedUp.tif'];
savename2   = [df filesep savename '_CleanedUp.tif'];

if exist(savename1)
    delete(savename1)
    delete(savename2)
end

for ks = 1:stacksize
    imwrite(ImgStack{ks}, savename1, 'WriteMode', 'append', 'Compression','none','resolution',Res);
end

slice = round(rand(1)*stacksize);
imwrite(ImgStack{slice}, savenamepic,'Compression','none','resolution',Res);

copyfile(savename1, savename2)

% Création de l'image dépliée a la bonne taille
[ly,lx] = getoutputsize(ImgStack{1},Xc(1),Yc(1));
UfI = uint8(zeros([lx,ly,stacksize]));

% Dépliage image par image
for ks = 1:stacksize
    Itmp = ImgStack{ks};
    
    
    %             figure(1)
    %             imshow(Itmp)
    %             hold on
    %             plot(Xc,Yc,'*')
    
    UfI(:,:,ks) = RadialUnfolding(Itmp,Xc(ks),Yc(ks));
    
    %     figure(2)
    %     imshow(UfI(:,:,ks))
    %     pause
    
    UfImask = UfI(:,:,ks) <2;
    A = all(UfImask,2);
    Aptr = find(diff(A) ~= 0);
    ptr = min([ptr,min(Aptr)]);
    
    
end

UfI = UfI(ptr+1:end,:,:);

savename3 = [ff filesep savename '.tif'];
savename4 = [df filesep savename '.tif'];

if exist(savename3)
    delete(savename3)
    delete(savename4)
end

for ks = 1:stacksize
    imwrite(UfI(:,:,ks), savename3, 'WriteMode', 'append', 'Compression','none','resolution',Res);

end
    copyfile(savename3, savename4)

end


function [Lx,Ly] = getoutputsize(pic,Xc,Yc)

pic = double(pic);

lx = size(pic,1);
ly = size(pic,2);

Lx = 361; % 360 deg

% rayon max = du centre au coin le plus loin
Ly = ceil(sqrt(2)*max(lx,ly))+1;

end
