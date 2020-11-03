clear all
close all
clc

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',10)
%% CHAMP MAXIMUM VOULU = 0.3 mT/mm

%% Experimental constants
W = [5 5]; % width of magnets in [mm]
R = [17.5 30]; % radius of magnets in [mm]

zexp = [4.9 6.0 7.1 8.2 9.3 10.4 12.6 15.9 20.4]; % measure points on the axis in [mm]

dz = 0.1;
z = 4.5:dz:60;

Bz = zeros(length(z),length(R));

Br = [1.32 1.32 1.37 1.32 1.32 1.32 1.32 1.32 1.37]; % Remanence field in T

col = [[1 0 0]' rand(3,length(R)-1)];

shift = (0:length(R)-1)*0;


figure(1)
hold on
xlabel('Distance from the magnet (mm)')
ylabel('Manetic field intensity (mT)')
title('Magnetic field produced by cylindrical magnet')
lgd1 = legend('on');

figure(2)
hold on
xlabel('Distance from the magnet (mm)')
ylabel('Manetic gradient (mT/mm)')
title('Magnetic gradient produced by cylindrical magnet')
lgd2 = legend('on');

%% Magnet model 

for i = 1:length(R)
cosa2 = (z + W(i))./sqrt((z + W(i)).^2 + R(i)^2);
cosa1 = z./sqrt(z.^2 + R(i)^2);
Bz(:,i) = Br(i)/2*(cosa2-cosa1);

figure(1)
plot(z,Bz(:,i)*10^3+shift(i),'-','color',col(:,i),'linewidth', 1.5)
lgd1.String{end} = ['Model for D = ' num2str(2*R(i)) 'mm and W = ' num2str(W(i)) '. Shifted of +' num2str(shift(i)) 'mT'];

DBzDz = (Bz(3:end,i)-Bz(1:end-2,i))/(2*dz);

figure(2)
plot(z(2:end-1),abs(DBzDz),'-','color',col(:,i),'linewidth', 1.5)
lgd2.String{end} = ['Model for D = ' num2str(2*R(i)) 'mm and W = ' num2str(W(i))];

end

% %% Experimental data
% 
% Bze = [195.4 177   159.8 144  129   115   92.3  66.3 43.4;...
%        247.2 222.7 200.7 181  161.8 145.3 109.7 80.5 52.9;...
%        156   133.5 113.2 95.8 80.8  69.2  50    32   19  ;...
%        196   166.8 142   120  101.7 81.3  62.7  40.4 24.3;...
%        244.5 207.7 175   150  127   108.5 80    52.4 31.8];
% 
% for i = 1:size(Bze,1)
%     figure(1)
%     plot(zexp,Bze(i,:)+shift(i+1),'-o','color',col(:,i+1),'markerfacecolor',col(:,i+1))
%     lgd1.String{end} = ['Data for D = ' num2str(2*R(i+1)) 'mm and W = ' num2str(W(i+1)) '. Shifted of +' num2str(shift(i+1)) 'mT'];
%     
%     
% end
% 
% hold off



figure(3)
hold on
xlabel('Distance from the magnet (mm)')
ylabel('Manetic field intensity (mT)')
title('Magnetic field produced by two cylindrical magnets')
lgd1 = legend('on');

figure(4)
hold on
xlabel('Distance from the center of the system (mm)')
ylabel('Manetic gradient (mT/mm)')
title('Magnetic gradient produced by two cylindrical magnets')
lgd2 = legend('on');

%% Double magnet model 
    Zmax = 10;

dz = 0.1;
zz = 4.5:dz:Zmax-4.5;

z2 = Zmax - zz;

Bz2 = zeros(length(zz),length(R));

for i = 1:length(R)
 
    
cosa2 = (zz + W(i))./sqrt((zz + W(i)).^2 + R(i)^2);
cosa1 = zz./sqrt(zz.^2 + R(i)^2);

cosa22 = (z2 + W(i))./sqrt((z2 + W(i)).^2 + R(i)^2);
cosa11 = z2./sqrt(z2.^2 + R(i)^2);

Bz2(:,i) = Br(i)/2*(cosa2-cosa1+cosa22-cosa11);
 
figure(3)
plot(zz-Zmax/2,Bz2(:,i)*10^3+shift(i),'-','linewidth', 1.5)
lgd1.String{end} = ['Model for D = ' num2str(2*R(i)) 'mm and W = ' num2str(W(i)) '.'];

DBz2Dz = (Bz2(3:end,i)-Bz2(1:end-2,i))*10^3/(2*dz);

figure(4)
plot(zz(2:end-1)-Zmax/2,abs(DBz2Dz),'-','linewidth', 1.5)
lgd2.String{end} = ['Model for D = ' num2str(2*R(i)) 'mm and W = ' num2str(W(i))];

end

