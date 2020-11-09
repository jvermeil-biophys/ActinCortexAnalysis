%% data exp 
Distance = [0 5 10 15 20 30 50]; % en mm
XL = [Distance(1) Distance(end)]; % pour plot

Champ = [500 520 566; 114 100 115; 30 27 28; 10 11 12; 5 5 5.5; 1.5 1.8 1.8; 0.3 0.4 0.3]; % en mT
% a la main sur feuille quadrillée pour distance


Db = 4503; % diametre M450

%% Calcul gradient
Gradient = abs(diff(Champ)'./diff(Distance)); % en mT/mm = en T/m

DistG = diff(Distance)/2 + Distance(1:end-1); % en mm

%% Calcul du champ aux points de gradient
ChampInt = interp1(Distance,Champ,DistG,'spline'); % en mT


%% Calcul de la magnétisation 

B = mean(ChampInt,2); % Champ Mag

M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % aimantation A.m^-1

V = 4/3*pi*(Db/2)^3; % en nm^3
                
m = M.*10^-9*V; % moment dipolaire en A.nm^2


%% Calcul de la force subie

GradB = mean(Gradient)'; % en T/m

F = m.*GradB*1e-6; % en pN

%% plots
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Serif')
set(0,'DefaultAxesFontWeight','normal')


figure

subplot(2,2,1)
plot(Distance,Champ','-o','linewidth',2)
% xlabel('Distance à l''aimant (mm)')
ylabel('Champ magnétique (mT)')
xlim(XL)

subplot(2,2,2)
plot(DistG,Gradient,'-o','linewidth',2)
% xlabel('Distance à l''aimant (mm)')
ylabel('Gradient (mT/mm)')
xlim(XL)

subplot(2,2,3)
plot(DistG,m,'-o','linewidth',2)
xlabel('Distance à l''aimant (mm)')
ylabel('Aimantation moyenne (A.m$^{-1}$)')
xlim(XL)

subplot(2,2,4)
plot(DistG,F,'-o','linewidth',2)
xlabel('Distance à l''aimant (mm)')
ylabel('Force moyenne (pN)')
xlim(XL)
ax = gca;
% ax.YScale = 'log';
% ax.XTick = [1 10 100 1000];





