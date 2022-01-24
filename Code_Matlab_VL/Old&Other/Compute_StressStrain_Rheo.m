function [Sigma,Epsilon] = Compute_StressStrain_Rheo(H,F,Eapprox,Fno,R,Hno)

% Compute stress and strain for rheological analysis of pinching beads experiments.
% Computes stress and strain on an oscilating curves using thickness
% without iscillation to compute the uncompressed thickness.
% Uses chadwick eq + module approx to find a the uncompressed thickness.

% Arguments : (CompressedThickness, Force, Module approx, 
% Force NoOscillation, Bead radius, Compressedthickness NoOscillation)
%
% Unités : (nm, pN, Pa, N, m, m)


%% Calcul de H0 en temps

 H0 = H0ChadCalc(Eapprox,Fno,R,Hno); % H0 en m par chadwick inverse
 
 
 %% Calcul de la contrainte et de la deformation en temps
                
                DeltaOsc = H0 - H; % Indentation en nm
                SurfCont = 2*pi*R*1e9*DeltaOsc*1e-6; % Surface de contact en µm²
                
                Sigma = F./SurfCont; % Contrainte en Pa
                Epsilon = DeltaOsc./(3*H0); % Deformation en %
 
end

 
 
 
function H0 = H0ChadCalc(E,F0,R,h)


% A partir de l'équation de Chadwick :
% F = (pi*E*R*delta^2)/(3*h0)
%
% On écrit delta = h0 - h et on resoud une
% equation du second degré pour h0
% a partir des valeurs de h (courbe BF) et EChad,
% et de la force sans l'oscilation F0
%
% h0²*pi*EChad*R - h0*(3*F0 + 2*h*pi*E*R) + pi*EChad*R*h² = 0



a = pi*E*R; % 'a' de l'eq second degree

b = -(3.*F0 + 2.*a.*h); % 'b' de l'eq second degree

c = a.*h.^2; % 'c' de l'eq second degree

disc = b.^2 - 4.*a.*c; % discriminant eq second degre en h0



h0m = (-b - sqrt(disc))./(2*a);

if ~all(h0m<h)
    error('Probleme avec le calcul de h0')
end

h0p = (-b + sqrt(disc))./(2*a); % racine + de l'eq, la seule plus grande que h

if ~all(h0p>h)
    error('Probleme avec le calcul de h0')
end

H0 = h0p*1e9; % épaisseur non déformé en nm

end