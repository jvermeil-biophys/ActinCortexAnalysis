function [Sdown,Smid,Sup,Sramp,Sfluo, idxRamp] = SplitZRampImg(Sfull,nimgbr,nimgr,nimg,wFluo,nBlackImagesOfLoopEnd)

%
% Splits the image list Sfull between list for down, middle and up images
% corresponding to the triplets taken in Z when not in a ramp.  
% And a list of image that are during a ramp. Improvement on SplitZimage
% for the case where there is non-triplet images in the middle of the loop.
%
%
% OUTPUTS : 
% Sdown : list of images that are the first of a Z triplet 
% Smid : list of images that are the second of a Z triplet 
% Sup : list of images that are the last of a Z triplet 
% Sramp : list of images in the ramps
%
% INPUTS :
% Sfull : vector of image numbers to split
% nimgbr : number of images before a ramp in a loop
% nimgr : number of images in a ramp in a loop
% nimg : number of images per loop
%

% initializing masks
N = length(Sfull);
N_loop = ceil(N/nimg);
Mdown = false(1,N_loop*nimg);
Mmid  = Mdown;
Mup   = Mdown;
Mramp = Mdown;
Mfluo = Mdown;
countRamp = zeros(1,N_loop*nimg);

nbImgFluo = 0;
if wFluo
    nbImgFluo = 1;
end
    
for i_loop = 1:N_loop
    for i = 1:3:nimgbr-2 % going by three for the triplets of images before ramp
        position = i + (i_loop-1)*nimg;
        Mdown(position) = true;
        Mmid(position+1) = true;
        Mup(position+2) = true;

    end

    for i = nimgbr+1 : nimgbr+nimgr - nBlackImagesOfLoopEnd(i_loop) % ramp images
        position = i + (i_loop-1)*nimg;
        Mramp(position) = true;
        countRamp(position) = i_loop;
    end

    for i = nimgbr+nimgr+1 - nBlackImagesOfLoopEnd(i_loop):3:nimg-nbImgFluo - nBlackImagesOfLoopEnd(i_loop) % going by three for the triplets of images after ramp
        position = i + (i_loop-1)*nimg;
        Mdown(position) = true;
        Mmid(position+1) = true;
        Mup(position+2) = true;
    end

    %    Mup = Mup|(mod(Sfull,nimg) == 0); % last image of a loop is always an up image

    if wFluo
        position = nimg - nBlackImagesOfLoopEnd(i_loop) + (i_loop-1)*nimg;
        Mfluo(position) = true;
    end
end

% for i = 1:3:nimgbr-2 % going by three for the triplets of images before ramp
% 
%     Mdown = Mdown|(mod(Sfull,nimg) == i);
%     Mmid = Mmid|(mod(Sfull,nimg) == i+1);
%     Mup = Mup|(mod(Sfull,nimg) == i+2);
% 
% end
% 
% for i = nimgbr+1:nimgbr+nimgr - nBlackImagesEachLoop % ramp images
%     Mramp =  Mramp|(mod(Sfull,nimg) == i);
% end
% 
% for i = nimgbr+nimgr+1 - nBlackImagesEachLoop:3: nimg - nbImgFluo - nBlackImagesEachLoop % going by three for the triplets of images after ramp
% 
%     Mdown = Mdown|(mod(Sfull,nimg) == i);
%     Mmid = Mmid|(mod(Sfull,nimg) == i+1);
%     Mup = Mup|(mod(Sfull,nimg) == i+2);
% 
% end
% 
% %    Mup = Mup|(mod(Sfull,nimg) == 0); % last image of a loop is always an up image
% 
% if wFluo
%     Mfluo = Mfluo|(mod(Sfull+nBlackImagesEachLoop*ones(1,max(Sfull)),nimg) == 0);
% end

Sdown = Sfull(Mdown);
Smid  = Sfull(Mmid) ;
Sup   = Sfull(Mup)  ;
Sramp = Sfull(Mramp);
Sfluo = Sfull(Mfluo);
idxRamp = countRamp;


end