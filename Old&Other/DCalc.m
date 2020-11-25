Dx = 4554; 
Dy = 709;
Dz = 482;

ErrDx = sqrt(2)*2;
ErrDy = sqrt(2)*2;
ErrZ = 45;
ErrDz = ErrZ*sqrt(2);
ErrR = 30;

ErrD3d = sqrt((Dx^2*ErrDx^2 + Dy^2*ErrDy^2 + Dz^2*ErrDz^2)/(Dx^2+Dy^2+Dz^2))

ErrD = sqrt(ErrD3d^2 + ErrR^2)
