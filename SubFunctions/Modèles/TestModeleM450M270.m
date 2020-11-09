
%% Paramètres

% Rayon de bille (en nm)

R1 = 2250; % M450

R2 = 1350; % M270

% Parametre manip venant du papier

H0 = 220; % épaisseur mediane du cortex a 5mT

H0cor = H0/87*100; % correction pôur l'enfoncement a 5mT

Delta50 = 120; % indentation max moyennée sur toute les comp a 50mT

% Variable modèle

a1approx = sqrt(R1*Delta50) % a est le rayon de la surface de contact entre la bille et le cortex

a1correct = sqrt(R1*Delta50 - (Delta50/2)^2)


a2approx = sqrt(R2*Delta50)

a2correct = sqrt(R2*Delta50 - (Delta50/2)^2)

% --> Meme avec un très grand delta il y a moins d'1% de diff pour a, on
% prend la version approximée


%% Condition de compression complète du cortex : a > h0/2 
% avec a = sqrt(R*delta), ou delta est l'indentation des deux billes
% on cherche le premier delta pour lequel la couche est
% completement comprimée. On prend aussi H0corrigé pour h0


% cas M450, on veut a1>h0/2, soit delta>h0²/4R1

deltamin1 = H0cor^2/(4*R1)


% cas M450/M270, on veut a2+a1>h0, soit delta>h0²/(sqrt(R1)+sqrt(R2))²

deltamin12 = H0cor^2/(sqrt(R1)+sqrt(R2))^2


% cas M270, on veut a2>h0/2, soit delta>h0²/4R2

deltamin2 = H0cor^2/(4*R2)

% --> dans tout les cas la couche est déja en compression complète a 5mT


%% Force entre les deux billes

B = 0:0.5:50;

% M450

D3 = H0 + 2*R1;

M1 = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1

V1 = 4/3*pi*R1^3; % en nm^3

m1 = M1.*10^-9*V1; % moment dipolaire en A.nm^2

Bind1 = 2*10^5*m1./(D3.^3); % champ induit par la magnétisation d'une bille voisine

Nvois = 1; % nombre de billes au contact

Btot1 = B + Nvois*Bind1;% B corrigé // B induit de la bille voisine

Mtot1 = 1.05*1600*(0.001991*Btot1.^3+17.54*Btot1.^2+153.4*Btot1)./(Btot1.^2+35.53*Btot1+158.1);

mtot1 = Mtot1.*10^-9*V1;



F1 = 3*10^5.*mtot1.^2./D3.^4; % en pN // sans la force des voisins




% M270 // on prend l'aimantation 2* plus faible plus le moment

D32 = H0 + 2*R2;

M2 = 0.5*(1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1)); % aimantation A.m^-1

V2 = 4/3*pi*R2^3; % en nm^3

m2 = M2.*10^-9*V2; % moment dipolaire en A.nm^2

Bind1 = 2*10^5*m2./(D32.^3); % champ induit par la magnétisation d'une bille voisine

Nvois = 1; % nombre de billes au contact

Btot2 = B + Nvois*Bind1;% B corrigé // B induit de la bille voisine

Mtot2 = 1.05*1600*(0.001991*Btot2.^3+17.54*Btot2.^2+153.4*Btot2)./(Btot2.^2+35.53*Btot2+158.1);

mtot2 = Mtot2.*10^-9*V2;



F2 = 3*10^5.*mtot2.^2./D32.^4; % en pN // sans la force des voisins





                
                
                
                
                
