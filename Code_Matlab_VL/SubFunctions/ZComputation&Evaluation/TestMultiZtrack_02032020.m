clear all
grid = -1500:10:1500;
PLOTperDepth = 1;
PLOTperDay = 1;


ndepth = 0;


%% 02 03 20 100X entre eux
Zerr = [];
allnum = {'1-1' '1-2' '1-3' '2-1' '2-2' '3-1' '3-2' '3-3' '3-4'}

for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = MultiZevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\02-03-20_Depthograph100X_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\02-03-20_Depthograph100X_Anti-' totest '.mat'],grid);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 02 03 20','fontsize',30)
hold on
plot([-3000 4000],[0 0],'k--')
plot([0 0],[-3000 4000],'k--')
plot(grid,grid,'r')
plot(grid,FakeCalc,'*')
xlabel('ZposReal','FontSize',20)
ylabel('ZposCalc','FontSize',20)
axes('Position',[.6 .18 .3 .3])
box on
hold on
plot(grid,Zerrme,'*')
plot(grid,ZerrFit,'r')
xlabel('Fake ZposCalc (nm)')
ylabel('Zerr (Real-Calc, nm)')
end

ndepth = ndepth+1;
AllZerrme(:,ndepth) = Zerrme;

clear Ze* P