function [Sdown,Smid,Sup] = SplitZimage(Sfull,nimg)

%
% Splits the image list Sfull between list for down, middle and up images
% corresponding to the triplets taken in Z. nimg is the number of image in
% a loops and is usefull when fluo images are not present and nimg is not a
% multiple of 3
%
% [Sdown,Smid,Sup] = SplitZimage(Sfull,nimg)
%
% OUTPUTS : 
% Sdown : list of images that are the first of a Z triplet 
% Smid : list of images that are the second of a Z triplet 
% Sup : list of images that are the last of a Z triplet 
%
% INPUTS :
% Sfull : vector of image numbers to split
% nimg : number of images per loop
%

% initializing masks
Mdown = false(1,length(Sfull));
Mmid  = Mdown;
Mup   = Mdown;



for i = 1:3:nimg % going by three for the triplets of images
 
    % using modulo function to get one image out of three in each list
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2); % does not work if nimg is multiple of 3

end

if (nimg/3) == fix(nimg/3) % if there is no fluo images (nimg is multiple of 3)
    Mup = Mup|(mod(Sfull,nimg) == 0); % godd formula for this case
end

% creating vectors
Sdown = Sfull(Mdown);
Smid  = Sfull(Mmid) ;
Sup   = Sfull(Mup)  ;

end