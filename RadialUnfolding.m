function Newpic = RadialUnfolding(Oldpic,Xc,Yc)

Oldpic = double(Oldpic);

lx = size(Oldpic,1);
ly = size(Oldpic,2);

Lx = 361; % 360 deg

% rayon max = du centre au coin le plus loin
Ly = ceil(sqrt(2)*max(lx,ly))+1;

Newpic =  zeros(Ly,Lx);
NewpicMask = Newpic;

for kx = 1:lx
    for ky = 1:ly
        [R,Theta] = CircCoord(kx,ky,Xc,Yc);
        if kx > Xc
            Theta = Theta + 91;
        else
            Theta = Theta + 271;
        end
        if R > 50
            R = R-49;
            Rm = round(R);
            Thetam = round(Theta);  
            sgR = sign(R-Rm);
            sgT = sign(Theta-Thetam);
            rR = abs(R-Rm);
            rT = abs(Theta-Thetam);

            Newpic(Rm,Thetam) = Newpic(Rm,Thetam) + abs(1-rR)*abs(1-rT)*Oldpic(ky,kx);
            NewpicMask(Rm,Thetam) = NewpicMask(Rm,Thetam) + abs(1-rR)*abs(1-rT);

            Newpic(Rm+sgR,Thetam) = Newpic(Rm+sgR,Thetam) + abs(rR)*abs(1-rT)*Oldpic(ky,kx);
            NewpicMask(Rm+sgR,Thetam) = NewpicMask(Rm+sgR,Thetam) + abs(rR)*abs(1-rT);
            
            Newpic(Rm,Thetam+sgT) = Newpic(Rm,Thetam+sgT) + abs(1-rR)*abs(rT)*Oldpic(ky,kx);
            NewpicMask(Rm,Thetam+sgT) = NewpicMask(Rm,Thetam+sgT) + abs(1-rR)*abs(rT);
            
            Newpic(Rm+sgR,Thetam+sgT) = Newpic(Rm+sgR,Thetam+sgT) + abs(rR)*abs(rT)*Oldpic(ky,kx);
            NewpicMask(Rm+sgR,Thetam+sgT) = NewpicMask(Rm+sgR,Thetam+sgT) + abs(rR)*abs(rT);
            
        end  
    end
end

Newpic = Newpic./NewpicMask;
Newpic = imgaussfilt(uint8(rot90(Newpic,2)));

end

