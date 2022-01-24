%test du Z au 63X du 04-02-2020

figure(1)
hold on
title('300nm')

figure(2)
hold on
title('500nm')

figure(3)
hold on
title('800nm')

figure(4)
hold on
title('1200nm')

names = {'03um_1-1' '03um_1-2' '03um_1-3'};


for k = 1:length(names)

[VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Manip_' names{k} '_Results']);

Z1 = ZCalc_Zcorr_63X(X1,Y1,S,['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Manip_' names{k}],'04-02-20_Depthograph','D:\Data\Raw');


Sfull = 1:S(end);
AntiS = Sfull(not(ismember(Sfull,S)));

figure
hold on
plot(S,Z1,'-*')
plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')


title(names{k})

dz = diff(Z1);

figure(1)
plot(Z1(dz>0),dz(dz>0),'*')

end


names = {'05um_1-1' '05um_1-2' '05um_1-3' '05um_1-4'};
for k = 1:length(names)

[VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Manip_' names{k} '_Results']);

Z1 = ZCalc_Zcorr_63X(X1,Y1,S,['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Manip_' names{k}],'04-02-20_Depthograph','D:\Data\Raw');


Sfull = 1:S(end);
AntiS = Sfull(not(ismember(Sfull,S)));

figure
hold on
plot(S,Z1,'-*')
plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')


title(names{k})

dz = diff(Z1);

figure(2)
plot(Z1(dz>0),dz(dz>0),'*')

end




names = {'08um_1-1' '08um_1-2' '08um_1-3' '08um_1-4' };
for k = 1:length(names)

[VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Manip_' names{k} '_Results']);

Z1 = ZCalc_Zcorr_63X(X1,Y1,S,['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Manip_' names{k}],'04-02-20_Depthograph','D:\Data\Raw');


Sfull = 1:S(end);
AntiS = Sfull(not(ismember(Sfull,S)));

figure
hold on
plot(S,Z1,'-*')
plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')


title(names{k})

dz = diff(Z1);

figure(3)
plot(Z1(dz>0),dz(dz>0),'*')

end





names = {'12um_1-1' '12um_1-2' '12um_1-3' '12um_1-4' };
for k = 1:length(names)

[VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Manip_' names{k} '_Results']);

Z1 = ZCalc_Zcorr_63X(X1,Y1,S,['D:\Data\Raw\20.02.04\TestNewPifoc_63X_Manip_' names{k}],'04-02-20_Depthograph','D:\Data\Raw');


Sfull = 1:S(end);
AntiS = Sfull(not(ismember(Sfull,S)));

figure
hold on
plot(S,Z1,'-*')
plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')


title(names{k})

dz = diff(Z1);

figure(4)
plot(Z1(dz>0),dz(dz>0),'*')

end
