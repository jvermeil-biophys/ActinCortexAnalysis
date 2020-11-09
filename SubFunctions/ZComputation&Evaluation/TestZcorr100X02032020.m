% test du Z au 100X du 02-03-2020


figure(3)
hold on
title('300nm')

figure(5)
hold on
title('500nm')

figure(8)
hold on
title('800nm')


names = {'03_1-1' '03_1-2' '03_1-3' '05_1-1' '05_1-2' '05_1-3' '05_2-1' '05_2-2' '05_2-3' '08_1-1' '08_1-2' '08_1-3'};



for k = 1:length(names)

[VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.03.02\ManipTest100X_' names{k} '_Results']);

Z1 = ZCalc_Zcorr_100X(X1,Y1,S,['D:\Data\Raw\20.03.02\ManipTest100X_' names{k}],'02-03-20_Depthograph100X','D:\Data\Raw');

Sfull = 1:S(end);
AntiS = Sfull(not(ismember(Sfull,S)));

figure
subplot(221)
hold on
plot(S,Z1,'-*')
plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')
title(names{k})

dz = diff(Z1);

subplot(222)
hold on
plot(Z1(1:end-1),dz,'*')

subplot(2,2,[3 4])
hold on
histogram(dz,'binwidth',20)
xlim([0 max(dz)+10])



figure(str2num(names{k}(2)))
plot(Z1(1:end-1),dz,'*')

end