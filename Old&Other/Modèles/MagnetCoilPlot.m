close all
clear all

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)

% Valeur du champ pour le montage bobine + aimants (Diam. 60mm Thick. 5mm)

% Commande labview (en mT)

ChmpLab = [-50 -40 -30 -20 -15 -10 -5 0 5 10 15 20 30 40 50];



% Champ mesuré (en mT) pour différentes distance entre les bobines (en mm)

Dist = [35 55 80];

ChmpExp = [4.4 15.8 28.1 38.8 44.6 50.4 56.1 61.9 67.4 73.5 78.9 84.7 96.6 108.8 120.8;...
           3.3 9.7  16.1 22.5 25.7 28.9 32.1 35.3 38.5 41.7 44.9 48.1 54.5 60.9  67.3 ;...
           2.5 6.2  9.4  12.8 14.5 16.2 17.9 19.6 21.3 23   24.8 26.5 29.9 33.3  36.7 ];
       
       
 plot(ChmpLab,ChmpExp,'*-')
 xlabel('Commande labview (mT)')
 ylabel('Champ mesuré (mT)')
 title('Champ pour deux Bobines-Aimants a différentes distances')
  legend('D = 35mm','D = 55mm','D = 80mm')
 
 
 Dce = diff(ChmpExp,1,2);
 
 Dcl = diff(ChmpLab);
 
 DchmpVal = Dce./Dcl;
 
 Dchmp = mean(DchmpVal,2);
 
 figure
 hold on
 yyaxis left
 plot(Dist,Dchmp,'*-')
 xlabel('Distance entre les bobines')
 ylabel('Coeff labview')
 yyaxis right
 plot(Dist,ChmpExp(:,8),'*-')
 ylabel('Champ des aimants (mT)')
