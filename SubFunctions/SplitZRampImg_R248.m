function [Sdown,Smid,Sup,Sramp] = SplitZRampImg_R248(Sfull,nimgbr,nimgar,nimgr1,nimgr2,nimgr3)

nimg = 3*nimgbr+nimgr1+3*nimgar+nimgr2+nimgr3;

Mdown = logical(zeros(1,length(Sfull)));
Mmid  = Mdown;
Mup   = Mdown;
Mramp   = Mdown;



for i = 1:3:nimgbr-2
 
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end

for i = nimgbr+1:nimgbr+nimgr1
    Mramp =  Mramp|(mod(Sfull,nimg) == i);
end

for i = nimgbr+nimgr1+1:3:2*nimgbr+nimgr1+nimgar
    
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end


for i = 2*nimgbr+nimgr1+nimgar+1:2*nimgbr+nimgr1+nimgar+nimgr2
    Mramp =  Mramp|(mod(Sfull,nimg) == i);
end

for i = 2*nimgbr+nimgr1+nimgar+nimgr2+1:3:3*nimgbr+nimgr1+2*nimgar+nimgr2
    
    Mdown = Mdown|(mod(Sfull,nimg) == i);
    Mmid = Mmid|(mod(Sfull,nimg) == i+1);
    Mup = Mup|(mod(Sfull,nimg) == i+2);

end


for i = 3*nimgbr+nimgr1+2*nimgar+nimgr2+1:3*nimgbr+nimgr1+2*nimgar+nimgr2+nimgr3
    Mramp =  Mramp|(mod(Sfull,nimg) == i);
end

for i = 3*nimgbr+nimgr1+2*nimgar+nimgr2+nimgr3+1:3:3*nimgbr+nimgr1+3*nimgar+nimgr2+nimgr3-1
    
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