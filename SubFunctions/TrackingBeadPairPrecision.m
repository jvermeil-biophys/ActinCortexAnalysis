%% INIT

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)

here = pwd;

T = T(Slice);
T = T - T(1);
T = T/1000;

pixelconv=1/15.8;

XM = XM*pixelconv*1000;
YM = YM*pixelconv*1000;

XM1 = XM1*pixelconv*1000;
YM1 = YM1*pixelconv*1000;


D = sqrt((XM-XM1).^2 + (YM-YM1).^2);




DT = diff(T);

DTptr = find(DT > 1);

DTptr = [0; DTptr; length(T)];
ndt =  length(DTptr);

for k = 1:ndt-1
    
    Tcut(:,k) = T(DTptr(k)+1:DTptr(k+1));
   
    Dcut(:,k) = D(DTptr(k)+1:DTptr(k+1));

      
end

StdD = std(Dcut);

DevD = mean(StdD);



figure
plot(Tcut,Dcut);
title('Dist vs. T')

Tcom = [Tcut(:,1), Tcut(:,1), Tcut(:,1), Tcut(:,1)];

f1 = figure
plot(Tcom,Dcut)
title(['Dist @ 3mT. Std = ' num2str(DevD) 'nm'])
xlabel('T (s)')
ylabel('Distance (nm)')


f1b = figure
plot(StdD,'o','markerfacecolor','k')
ax = gca;
ax.FontSize = 10;

inset(f1,f1b,0.2)