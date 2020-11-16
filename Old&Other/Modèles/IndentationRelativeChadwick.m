% Chadwick : F = (2*pi*E*R*Delta^2)/(3*H0)
%
% dans notre cas,on a une bille de chaque 
% coté et donc on considere que le 
% Delta de chadwick est delta/2 et que
% le H0 de chadwick est h0/2 avec
% delta et h0 les valeurs qui correspondent
% a notre situation experimentale
% l'equation est donc :
% 
% F = (pi*E*R*delta^2)/(3*h0)
%
% On écrit delta = h0 - h et on resoud une
% equation du second degré pour h0
% a partir des valeurs mediane de h et E, 
% et de la force moyenne de 70pN
%
% pi*E*R*h0² - h0*(3*F + 2*h*pi*E*R) + pi*E*R*h² = 0
%



E = 7e3; % module moyen en Pa
R = 4438e-9/2; % rayon bille en m
F = 70e-12; % Force a 5mT en N
h = 220e-9; % epaisseur mediane à 5mT en m


a = pi*E*R; % 'a' de l'eq second degree

b = -(3*F + 2*a*h); % 'b' de l'eq second degree

c = a*h^2; % 'c' de l'eq second degree

disc = b^2 - 4*a*c; % discriminant eq second degre en h0

h0m = (-b - sqrt(disc))/(2*a);
h0p = (-b + sqrt(disc))/(2*a) % racine + de l'eq, la seule plus grande que h


delta = h0p - h % indentation (m)

relind = delta/h0p % indentation relative 

rad = sqrt(R*delta)*1e9  % radius of contact area (nm)

surf = pi*rad^2*1e-6 % surface de contact en µm²

