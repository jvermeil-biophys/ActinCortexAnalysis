function [Sdown,Smid,Sup,Sramp] = SplitZRampImg(Sfull,nimgbr,nimgr,nimg)

Mdown = logical(zeros(1,length(Sfull)));
Mmid  = Mdown;
Mup   = Mdown;
Mramp   = Mdown;



for i = 1:3:nimgbr-2
 
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end

for i = nimgbr+1:nimgbr+nimgr
    Mramp =  Mramp|(mod(Sfull,nimg) == i);
end

for i = nimgbr+nimgr+1:3:nimg-1
    
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end


    Mup = Mup|(mod(Sfull,nimg) == 0);


Sdown = Sfull(Mdown);
Smid  = Sfull(Mmid) ;
Sup   = Sfull(Mup)  ;
Sramp = Sfull(Mramp);



end