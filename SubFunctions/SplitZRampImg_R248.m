function [Sdown,Smid,Sup,Sramp] = SplitZRampImg_R248(Sfull,nimgbr,nimgar,nimgr1,nimgr2,nimgr3)

%
% Splits the image list Sfull between list for down, middle and up images
% corresponding to the triplets taken in Z when not in a ramp.  
% And a list of image that are during a ramp. Improvement on SplitZRampImg
% for the case where there is three different type of ramps.
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
% nimgar : number of images after a ramp in a loop
% nimgr1: number of images in first ramp
% nimgr2: number of images in second ramp
% nimgr3: number of images in third ramp
%

nimg = 3*nimgbr+nimgr1+nimgr2+nimgr3+3*nimgar; % total number of image


% initializing masks
Mdown = false(1,length(Sfull));
Mmid  = Mdown;
Mup   = Mdown;
Mramp   = Mdown;



for i = 1:3:nimgbr-2 % before ramp 1
 
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end

for i = nimgbr+1:nimgbr+nimgr1 % ramp 1
    Mramp =  Mramp|(mod(Sfull,nimg) == i);
end

for i = nimgbr+nimgr1+1:3:2*nimgbr+nimgr1+nimgar % between ramp 1 and 2
    
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end


for i = 2*nimgbr+nimgr1+nimgar+1:2*nimgbr+nimgr1+nimgar+nimgr2 % ramp 2
    Mramp =  Mramp|(mod(Sfull,nimg) == i);
end

for i = 2*nimgbr+nimgr1+nimgar+nimgr2+1:3:3*nimgbr+nimgr1+2*nimgar+nimgr2 % between ramp 2 and 3
    
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end


for i = 3*nimgbr+nimgr1+2*nimgar+nimgr2+1:3*nimgbr+nimgr1+2*nimgar+nimgr2+nimgr3 % ramp 3
    Mramp =  Mramp|(mod(Sfull,nimg) == i);
end

for i = 3*nimgbr+nimgr1+2*nimgar+nimgr2+nimgr3+1:3:3*nimgbr+nimgr1+3*nimgar+nimgr2+nimgr3-1 % after ramp 3
    
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end


    Mup = Mup|(mod(Sfull,nimg) == 0); % last image of a loop is always an up image


Sdown = Sfull(Mdown);
Smid  = Sfull(Mmid) ;
Sup   = Sfull(Mup)  ;
Sramp = Sfull(Mramp);



end