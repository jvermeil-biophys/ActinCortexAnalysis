clear all
close all
clc

set(0,'DefaultFigureWindowStyle','docked')

%% Constantes 

% Valeur liées aux aimants
W =  [5 5]; % width of magnets in [mm]
R =  [17.5 30]; % radius of magnets in [mm]
Br = [1.32 1.32 1.37 1.32 1.32 1.32 1.32]; % Remanence field in [T]

% Champ mag voulus
Bexp = [5 10 20 40 80 100]; % en [mT]


%% Définitions des variables

dz = 0.05; % en [mm]
z = 10 :dz:50; % distance aimant - mesure en [mm]

Bz = zeros(length(z),length(R)); % chmp mag
Gz = zeros(length(z),length(R)); % grad chmp mag


% Champ mag au milieu de l'axe entre deux aimants séparés de 2*z
for i = 1:length(R)
Bz(:,i) = 2*Bcalc(z,R(i),W(i),Br(i)); % en [mT]
Gz(:,i) = Gcalc(z,R(i),W(i),Br(i),4); % en [mT/mm]
end


figure(1)
hold on
xlabel('Distance between the magnets (mm)')
ylabel('Manetic field intensity on the axis between the magnets (mT)')
title('Magnetic field produced between two cylindrical magnets')
lgd1 = legend('on');

for i = 1:length(R)
    plot(2*z,Bz(:,i),'-','linewidth', 1.5)
    lgd1.String{end} = ['Model for D = ' num2str(2*R(i)) 'mm and W = ' num2str(W(i)) '.'];
end

hold off




figure(2)
hold on
xlabel('Magnetic field between the magnets (mT)')
ylabel('Manetic gradient between the magnets (mT/mm)')
title('Magnetic gradient vs. magnetic field')
lgd2 = legend('on');

for i = 1:length(R)
    plot(Bz(:,i),Gz(:,i),'-','linewidth', 1.5)
    lgd2.String{end} = ['Model for D = ' num2str(2*R(i)) 'mm and W = ' num2str(W(i)) '.'];
end
















function B = Bcalc(z,R,W,Br) % calcul le champ sur l'axe d'un aimant défini
% par R,W,Br à une distance z
cosa2 = (z + W)./sqrt((z + W).^2 + R^2);
cosa1 = z./sqrt(z.^2 + R^2);
B = Br/2*(cosa2-cosa1)*10^3; % en [mT]
end

function G = Gcalc(z,R,W,Br,L) % calcul le gradient moyen sur une 
% longueur 2*L centré entre deux aimants définis par R,W,Br et séparé de z
G = (Bcalc(z+L,R,W,Br) + Bcalc(z-L,R,W,Br) - 2*Bcalc(z,R,W,Br))/L;
end
