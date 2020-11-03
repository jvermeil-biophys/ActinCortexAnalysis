

function Theta = Angle(x,y,Xc,Yc)
Theta = radtodeg(atan((y-Yc)./(x-Xc)));
end