% test du Z au 63X du 02-03-2020


figure(3)
hold on
title('300nm')

figure(5)
hold on
title('500nm')

figure(8)
hold on
title('800nm')


names = {'03_1-1' '03_1-2' '05_1-1' '05_1-2' '05_2-1' '05_2-2' '08_1-1' '08_1-2'};



for k = 1:length(names)

[VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.03.02\ManipTest63X_' names{k} '_Results']);

Z1 = ZCalc_Zcorr_63X(X1,Y1,S,['D:\Data\Raw\20.03.02\ManipTest63X_' names{k}],'02-03-20_Depthograph63X','D:\Data\Raw');

Sfull = 1:S(end);
AntiS = Sfull(not(ismember(Sfull,S)));

figure
hold on
plot(S,Z1,'-*')
plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')

title(names{k})

figure(str2num(names{k}(2)))
plot(diff(Z1),'*')
end
