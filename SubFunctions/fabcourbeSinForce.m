function fabcourbeSinForce(meanF,ampF,FreqExp,lengthManip)

% fabcourbeSinForce(meanF,ampF,FreqExp,lengthManip)
% Genere un fichier texte pour labview avec une force sinusoidale 
% d'amplitude ampF (pN) et de moyenne meanF (pN)
% avec une frequence FreqExp (Hz)
% sur un temps de lengthManip (minute)

f=10000; % frequence d'echantillonage (10,000 Hz)

%% calcul de la force dipolaire; correspondance entre B et F de 5 a 40 mT

                Db = 4503; % diam bille en nm

                Dnm = Db + 200; % dist entre centre des billes en nm
                
                B = 5:40; % champ en mT
                
                M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1
                
                V = 4/3*pi*(Db/2)^3; % en nm^3
                
                m = M.*10^-9*V; % moment dipolaire en A.nm^2
                
                Bind = 2*10^5*m./(Dnm.^3); % champ induit par la magnétisation d'une bille voisine
                
                Nvois = 1; % nombre de billes au contact
                
                Btot = B + Nvois*Bind;% B corrigé // B induit de la bille voisine
                
                Mtot = 1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
                
                mtot = Mtot.*10^-9*V;
                
%                 anglefactor=abs(3*(D./D3).^2-1); % effet angle de la chaine / champ B
                
                F = 3*10^5.*mtot.^2./Dnm.^4; % en pN // sans correction d'angle



%% calcul de Bexp pour la force voulue Fexp

T = (0:1/f:lengthManip*60);


Fexp = meanF + ampF*sin(2*pi*T*FreqExp);


Bexp = interp1(F,B,Fexp,'spline')*1000;

% /!\ le champ doit etre donné en µT /!\


figure
hold on
yyaxis left
plot(T,Fexp)
ylabel('Force (pN)')
yyaxis right
plot(T,Bexp)
ylabel('Champ mag (µT)')
xlabel('T (s)') 


dlmwrite(['C:\Users\Valentin\Dropbox\TheseValentin\TxtLabview\SinusForce_mean' num2str(meanF) 'pN_amp' num2str(ampF) 'pN_freq' num2str(FreqExp) 'Hz_length' num2str(lengthManip) 'min.txt'],Bexp','newline','pc','precision','%1.0f');

end