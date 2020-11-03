
function [R,Theta] = CircCoord(x,y,Xc,Yc)
R = Rayon(x,y,Xc,Yc);
Theta = Angle(x,y,Xc,Yc);
end