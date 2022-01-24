grid = -1500:10:1500;
PLOTperDepth = 1;
PLOTperDay = 1;


ndepth = 0;

%% 25 02 20 entre eux
Zerr = [];
allnum = {'100nmUp63XSat-1' '100nmUp63XSat-2' '50nmUp63XSat-1' '50nmUp63XSat-2' '20nmUp63XSat-1' '20nmUp63XSat-2'}

for k= 1:length(allnum)
    totest = allnum{k};
    for kk = 1:length(allnum)
        if kk ~=k
            against = allnum{kk};
            
            Zerr(:,k,kk) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\25-02-20_Deptho_' totest '.mat'],...
                ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\25-02-20_Deptho_' against '.mat'],grid,PLOTperDepth);
            
        end
    end
end

Zerrme = nanmean(Zerr,2)';


P = polyfit(grid,Zerrme,15);
ZerrFit = polyval(P,grid);

FakeCalc = grid-Zerrme;
figure
title('Mean err 10 02 20','fontsize',30)
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


ndepth = ndepth+1;
AllZerrme(:,ndepth) = Zerrme;

clear Ze* P