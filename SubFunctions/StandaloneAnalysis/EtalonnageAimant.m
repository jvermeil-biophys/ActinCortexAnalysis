%% data exp 
Distance = [0 5 10 15 20 30 50]; % en mm

Champ = [500 520 566; 114 100 115; 30 27 28; 10 11 12; 5 5 5.5; 1.5 1.8 1.8; 0.3 0.4 0.3]; % en mT

Db = 4503; % diametre M450

%% Calcul gradient
Gradient = abs(diff(Champ)'./diff(Distance)); % en mT/mm = en T/m

DistG = diff(Distance)/2 + Distance(1:end-1); % en mm

%% Calcul du champ aux points de gradient
ChampInt = interp1(Distance,Champ,DistG,'spline'); % en mT


%% plot champ et gradient
figure
plot(Distance,Champ','-o')
xlabel('Distance à l''aimant (mm)')
ylabel('Champ magnétique mésuré (mT)')

figure
plot(DistG,Gradient,'-o')
xlabel('Distance à l''aimant (mm)')
ylabel('Gradient de champ (mT/mm)')


%% Calcul de la magnétisation 

B = mean(ChampInt,2); % Champ Mag

M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1

V = 4/3*pi*(Db/2)^3; % en nm^3
                
m = M.*10^-9*V; % moment dipolaire en A.nm^2


%% Calcul de la force subie

GradB = mean(Gradient)'; % en T/m

F = m.*GradB*1e-6 % en pN

