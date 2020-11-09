% Donnée sur les manips de compression ou on regarde la relation
% hauteur/stiffness en ajoutant ou retirant artificielment de la hauteur

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextInterpreter','Latex')



height25 =    [130    110   80     50   30    0    -10   -30  -50  -70    -90   -110    -130   -150  -170  -200 -220  -240  -260 -280   -300];
powexp25 =    [-2.3  -2.2  -1.9  -1.7 -1.5  -1.3  -1.2  -1.1  -1  -0.85  -0.8   -0.79   -0.74 -0.69 -0.76  -0.7 -0.54 -0.6 -0.68 -0.93 -0.95];
corrcoef25 =  -0.001*[463 461 442 427 421 441 440 433 395 365 361 334 315 337 274 381 384 423 403 497 597];
ncomp =     [246   245   246    246  248   246   244   237  232  228   220     196    179    156    118   86   73    46    36    30    26 ];
weightedE25 = [9.2   8.9   8.5    8.0  7.4   7.0   6.8   6.3  5.9  5.4   5.2    4.9     4.5    4.2   4.0    3.6  3.3   3.1   3.0   2.8   2.5];

SEpowexp25 = [0.18 0.17 0.16 0.14 0.13 0.11 0.11 0.11 0.10 0.10 0.09 0.1 0.1 0.1 0.11 0.12 0.12 0.17 0.23 0.25 0.23];
CIpowexp25 = tinv(0.95,ncomp-2).*SEpowexp25;

height100    = [-300  -260  -220  -170  -130  -90   -50  -30  -10  0];
powexp100    = [-0.84 -0.77 -0.82 -0.92 -0.91 -0.92 -1.2 -1.2 -1.2 -1.3];
weightedE100 = [2.3   2.5   2.6   3.4   4.1   4.7   5.5  5.8  6.2  6.4];


selectheight = [0 -30 -50 -90 -110 -150 -200 -260 -300];
powexpselect = [-1.3 -1.4 -1.5 -1.7 -1.9 -2.1 -2.5 -2.4 -2.7];
Eselect = [7 7 7 7 7 7.1 6.3 5.4 5.3];


figure
hold on
plot(height25,powexp25,'ok','markerfacecolor','b')
plot(height100,powexp100,'ok-','markerfacecolor','g')
plot(selectheight,powexpselect,'k-o','markerfacecolor','r')
title('Power law vs height added, min $>0$')
xlabel('Added height (nm)')
ylabel('Power law exponent')

legend('Exponent vs height 25%','Exponent vs height 100%','filtered height 25%')

ptr0 = height25==0;

figure
hold on
plot([150 -320],[powexp25(ptr0) powexp25(ptr0)],'r-','linewidth',1.3)
plot([150 -320],powexp25(ptr0)+[CIpowexp25(ptr0) CIpowexp25(ptr0)],'r--','linewidth',0.5)
plot([150 -320],powexp25(ptr0)-[CIpowexp25(ptr0) CIpowexp25(ptr0)],'r--','linewidth',0.5)
errorbar(height25,powexp25,CIpowexp25,'k-')
plot(height25,powexp25,'ok','markerfacecolor','b')
title('Power law vs height added, min $>0$')
xlabel('Added height (nm)')
ylabel('Power law exponent')
xlim([-320 150])



figure
hold on
plot(height25,corrcoef25,'ok-','markerfacecolor','b')
title('Corrcoef vs height added, min $>0$')
xlabel('Added height (nm)')
ylabel('Correlation coefficient')

figure
hold on
plot(height25,ncomp,'ok-','markerfacecolor','b')
title('Number of compression with a $>0$ minimum')
xlabel('Added height (nm)')
ylabel('Number of compression')

figure
hold on
plot(height25,weightedE25,'ok-','markerfacecolor','b')
plot(height100,weightedE100,'ok-','markerfacecolor','g')
plot(selectheight,Eselect,'k-o','markerfacecolor','r')
title('Weighted modulus, min $>0$')
xlabel('Added height (nm)')
ylabel('Weighted average E (kPa)')

legend('E vs height 25%','E vs height 100%','filtered height 25%')

