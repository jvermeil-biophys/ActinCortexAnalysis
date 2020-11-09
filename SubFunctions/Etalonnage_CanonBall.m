%% Données

Pos = [-4 -3 -2 -1 0 1 2 3 4]';

Chmp5 = [1.4 2.4 3.3 4.3 5.3 6.3 7.3 8.3 9.3]';

Chmp10 = [6.4 7.3 8.2 9.2 10.2 11.2 12.2 13.3 14.3]';

Chmp20 = [16.3 17.2 18.1 19.1 20 21.1 22.1 23.1 24.2]';

Chmp30 = [27.5 28 28.6 29.1 29.8 30.4 31.2 31.9 32.8]';

Chmp50 = [47.3 47.7 48.2 48.8 49.4 50.1 50.9 51.7 52.6]';


%% Fits

[LM5_1,gof] = fit(Pos(1:5),Chmp5(1:5),'poly1');
Rsq5_1 = gof.rsquare;
[LM5_2,gof] = fit(Pos(5:end),Chmp5(5:end),'poly1');
Rsq5_2 = gof.rsquare;

MidField5 = mean([LM5_1.p2 LM5_2.p2])
Grad5 = mean([LM5_1.p1 LM5_2.p1])


[LM10_1,gof] = fit(Pos(1:5),Chmp10(1:5),'poly1');
Rsq10_1 = gof.rsquare;
[LM10_2,gof] = fit(Pos(5:end),Chmp10(5:end),'poly1');
Rsq10_2 = gof.rsquare;

MidField10 = mean([LM10_1.p2 LM10_2.p2])
Grad10 = mean([LM10_1.p1 LM10_2.p1])


[LM20_1,gof] = fit(Pos(1:5),Chmp20(1:5),'poly1');
Rsq20_1 = gof.rsquare;
[LM20_2,gof] = fit(Pos(5:end),Chmp20(5:end),'poly1');
Rsq20_2 = gof.rsquare;

MidField20 = mean([LM20_1.p2 LM20_2.p2])
Grad20 = mean([LM20_1.p1 LM20_2.p1])


[LM30_1,gof] = fit(Pos(1:5),Chmp30(1:5),'poly1');
Rsq30_1 = gof.rsquare;
[LM30_2,gof] = fit(Pos(5:end),Chmp30(5:end),'poly1');
Rsq30_2 = gof.rsquare;

MidField30 = mean([LM30_1.p2 LM30_2.p2])
Grad30 = mean([LM30_1.p1 LM30_2.p1])


[LM50_1,gof] = fit(Pos(1:5),Chmp50(1:5),'poly1');
Rsq50_1 = gof.rsquare;
[LM50_2,gof] = fit(Pos(5:end),Chmp50(5:end),'poly1');
Rsq50_2 = gof.rsquare;

MidField50 = mean([LM50_1.p2 LM50_2.p2])
Grad50 = mean([LM50_1.p1 LM50_2.p1])



%% Plots

figure
hold on
title('Etalonnage 5mT')
plot(Pos, Chmp5,'b-*')
plot(Pos(1:5),LM5_1(Pos(1:5)),'r--','linewidth',1.5)
plot(Pos(5:end),LM5_2(Pos(5:end)),'m--','linewidth',1.5)
legend('Data','Fit1','Fit2')

figure
hold on
title('Etalonnage 10mT')
plot(Pos, Chmp10,'b-*')
plot(Pos(1:5),LM10_1(Pos(1:5)),'r--','linewidth',1.5)
plot(Pos(5:end),LM10_2(Pos(5:end)),'m--','linewidth',1.5)

figure
hold on
title('Etalonnage 20mT')
plot(Pos, Chmp20,'b-*')
plot(Pos(1:5),LM20_1(Pos(1:5)),'r--','linewidth',1.5)
plot(Pos(5:end),LM20_2(Pos(5:end)),'m--','linewidth',1.5)

figure
hold on
title('Etalonnage 30mT')
plot(Pos, Chmp30,'b-*')
plot(Pos(1:5),LM30_1(Pos(1:5)),'r--','linewidth',1.5)
plot(Pos(5:end),LM30_2(Pos(5:end)),'m--','linewidth',1.5)

figure
hold on
title('Etalonnage 50mT')
plot(Pos, Chmp50,'b-*')
plot(Pos(1:5),LM50_1(Pos(1:5)),'r--','linewidth',1.5)
plot(Pos(5:end),LM50_2(Pos(5:end)),'m--','linewidth',1.5)
