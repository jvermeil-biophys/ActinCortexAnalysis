path = "C:\Users\JosephVermeil\Desktop\testSizeDetector_im1";
X1 = 804.883;
Y1 = 277.774;
X2 = 863.673;
Y2 = 281.649;
S = 1;
scale = 15.8;
diameterL = 4.5*scale;
diameterS = 2.7*scale;
whosWho = AsymmetricBeadPair_findWhosWho(X1,Y1,X2,Y2,S,path,diameterL,diameterS);

DIAMETER = "4500_2700";
if (contains(DIAMETER, '_'))
    diameters = split(DIAMETER,"_");
    diameters = str2double(diameters);
    iL = 1 + (max(diameters) == diameters(2));
    iS = 3 - iL;
    diameterL = diameters(iL);
    diameterS = diameters(iS);
end