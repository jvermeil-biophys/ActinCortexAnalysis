function [Sdown,Smid,Sup] = SplitZimage(Sfull,nimg)

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