clear all

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',23)

name = {'1mT_1-1','1mT_1-2','1mT_2-1','1mT_2-2','1mT_2-3','1mT_3-1','1mT_3-2','1mT_3-3'...
    ,'5mT_2-1','5mT_2-2','5mT_2-3','5mT_2-4','5mT_2-5','5mT_2-6','5mT_2-7','5mT_2-8','5mT_2-9'...
    ,'5mT_2-10','5mT_2-11','5mT_2-12'};

for kn = 1:length(name)

[VarName1,Area,StdDev,XM,YM,Slice] = importfilegeneric(['D:\Data\Raw\19.06.11_BrownianTest\' name{kn} '.txt']);

for k = 1:max(Slice)
    ptrslice = find(Slice==k);
    if length(ptrslice) == 2
        dx = XM(ptrslice(1))-XM(ptrslice(2));
        dy = YM(ptrslice(1))-YM(ptrslice(2));
        d(k) = sqrt(dx^2+dy^2);
    end
end

d = d(100:400)/15.8*1000;

figure
plot(d)
title([name{kn} ' - Mean : ' num2str(mean(d)) ' - Std : ' num2str(std(d))])

end


clear all

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',23)

name = {'1mT_1-1','1mT_1-2','1mT_1-3','1mT_1-4','1mT_2-1','1mT_2-2','1mT_2-3','1mT_2-4',...
    '1mT_2-5','1mT_2-6','1mT_2-7','3mT_1-1','3mT_1-2','3mT_1-3','3mT_1-4','3mT_1-5'...
    ,'3mT_1-6','3mT_1-7'};

for kn = 1:length(name)

[VarName1,Area,StdDev,XM,YM,Slice] = importfilegeneric(['D:\Data\Raw\19.06.17_BrownianTest\' name{kn} '.txt']);

for k = 1:max(Slice)
    ptrslice = find(Slice==k);
    if length(ptrslice) == 2
        dx = XM(ptrslice(1))-XM(ptrslice(2));
        dy = YM(ptrslice(1))-YM(ptrslice(2));
        d(k) = sqrt(dx^2+dy^2);
    end
end

d = d(100:400)/15.8*1000;

figure
plot(d)
title([name{kn} ' - Mean : ' num2str(mean(d)) ' - Std : ' num2str(std(d))])

end