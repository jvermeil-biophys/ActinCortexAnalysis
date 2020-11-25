function [Sdown,Smid,Sup] = SplitZimage(Sfull,nimg)

% [Sdown,Smid,Sup] = SplitZimage(Sfull,nimg)
%
% Splits the image list Sfull between list for down, middle and up images
% corresponding to the triplets taken in Z. nimg is the number of image in
% a loops and is usefull when fluo images are not present and nimg is not a
% multiple of 3

Mdown = logical(zeros(1,length(Sfull)));
Mmid  = Mdown;
Mup   = Mdown;



for i = 1:3:nimg
 
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end

if (nimg/3) == fix(nimg/3)
    Mup = Mup|(mod(Sfull,nimg) == 0);
end

Sdown = Sfull(Mdown);
Smid  = Sfull(Mmid) ;
Sup   = Sfull(Mup)  ;




end