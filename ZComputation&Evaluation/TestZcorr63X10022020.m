
% test du Z au 63X du 10-02-2020


figure(1)
hold on
title('500nm')

figure(2)
hold on
title('800nm')

figure(3)
hold on
title('1300nm')


names = {'05_1-1' '05_1-2'};

% 
% for k = 1:length(names)
% 
% [VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.02.10\TestManipPifoc63X_' names{k} '_Results']);
% 
% Z1 = ZCalc_Zcorr_63X(X1,Y1,S,['D:\Data\Raw\20.02.10\TestManipPifoc63X_' names{k}],'10-02-20_Depthograph','D:\Data\Raw');
% 
% Sfull = 1:S(end);
% AntiS = Sfull(not(ismember(Sfull,S)));
% 
% figure
% hold on
% plot(S,Z1,'-*')
% plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')
% 
% title(names{k})
% 
% dZ1 = Z1(3:end)-Z1(1:end-2);
% 
% ptr = find(ismember(S,[5 10 20 25 35 40 50 55 65 70 80 85 95 100 110 115 125 130 140 145 155 160 170 175]-1))
% 
%  figure(1)
% plot(Z1(ptr(1:end-1)),dZ1(ptr(1:end-1)),'*')
% 
% ylim([0 550])
% 
% 
% 
% end
% 
% 
% names = {'08_2-1' '08_2-2'};
% 
% 
% for k = 1:length(names)
% 
% [VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.02.10\TestManipPifoc63X_' names{k} '_Results']);
% 
% Z1 = ZCalc_Zcorr_63X(X1,Y1,S,['D:\Data\Raw\20.02.10\TestManipPifoc63X_' names{k}],'10-02-20_Depthograph','D:\Data\Raw');
% 
% Sfull = 1:S(end);
% AntiS = Sfull(not(ismember(Sfull,S)));
% 
% figure
% hold on
% plot(S,Z1,'-*')
% plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')
% 
% title(names{k})
% 
% dZ1 = Z1(3:end)-Z1(1:end-2);
% 
% ptr = find(ismember(S,[5 10 20 25 35 40 50 55 65 70 80 85 95 100 110 115 125 130]-1))
% 
%  figure(2)
% plot(Z1(ptr(1:end-1)),dZ1(ptr(1:end-1)),'*')
% 
% 
% 
% ylim([0 850])
% 
% 
% end


names = {'13_2-1' '13_2-2' '13_2-3'};


for k = 1:length(names)

[VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.02.10\TestManipPifoc63X_' names{k} '_Results']);

Z1 = ZCalc_Zcorr_63X(X1,Y1,S,['D:\Data\Raw\20.02.10\TestManipPifoc63X_' names{k}],'10-02-20_Depthograph','D:\Data\Raw');

Sfull = 1:S(end);
AntiS = Sfull(not(ismember(Sfull,S)));

figure
hold on
plot(S,Z1,'-*')
plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')

title(names{k})

dZ1 = Z1(3:end)-Z1(1:end-2);

ptr = find(ismember(S,[5 10 20 25 35 40 50 55 65 70 80 85 95 100 110 115 125 130]-1))

 figure(3)
plot(Z1(ptr(1:end-1)),dZ1(ptr(1:end-1)),'*')


ylim([0 1500])


end
