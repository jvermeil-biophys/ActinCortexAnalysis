grid = -1500:10:1500;
PLOTperDepth = 1;
PLOTperDay = 1;


ndepth = 0;

%% 17-05-18 entre eux
Zerr = [];
allnum = {'1','2','3','4','5','6','7','8'};
for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\17-05-18_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\17-05-18_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 17-05-18','fontsize',30)
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

%% 12-07-18 entre eux
Zerr = [];
allnum = {'1','2-1','2-2','3'};
for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\12-07-18_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\12-07-18_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 12-07-18','fontsize',30)
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

%% 18-07-18 entre eux
Zerr = [];
allnum = {'1-1','1-2','2-1','2-2'};
for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\18-07-18_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\18-07-18_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 18-07-18','fontsize',30)
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

%% 31-07-18 entre eux
Zerr = [];
allnum = {'1','2','3','4'};
for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\31-07-18_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\31-07-18_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 31-07-18','fontsize',30)
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

%% 28 08 18 entre eux
Zerr = [];
allnum = {'1','2','3'};
for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\28-08-18_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\28-08-18_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 28 08 18','fontsize',30)
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

%% 30 10 18 entre eux
Zerr = [];
allnum = {'1','2','3','4','5'};
for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\30-10-18_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\30-10-18_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 30 10 18','fontsize',30)
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

%% 12 12 18 entre eux
Zerr = [];
allnum = {'1','2-1','2-2','3','4'};
% allnum = {'7'};
for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\12-12-18_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\12-12-18_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 12 12 18','fontsize',30)
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

%% 15 03 19 entre eux
Zerr = [];
allnum = {'1','2','3','4','6','7','8','9-1','9-2'};
% allnum = {'7'};
for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\15-03-19_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\15-03-19_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 15 03 19','fontsize',30)
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

%% 04 02 20 entre eux
Zerr = [];
allnum = {'1-1' '1-2' '1-3' '2-1' '2-2' '3-1' '3-2' '3-3' '3-4'}

for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\04-02-20_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\04-02-20_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
    P = polyfit(grid,Zerrme,15);
    ZerrFit = polyval(P,grid);
    
FakeCalc = grid-Zerrme;
figure
title('Mean err 04 02 20','fontsize',30)
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

%% 10 02 20 entre eux
Zerr = [];
allnum = {'1' '2' '3' '4' '5' '6-1' '6-2'}

for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\10-02-20_Depthograph_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\10-02-20_Depthograph_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
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
end

ndepth = ndepth+1;
AllZerrme(:,ndepth) = Zerrme;

clear Ze* P

%% 02 03 20 100X entre eux
Zerr = [];
allnum = {'1-1' '1-2' '1-3' '2-1' '2-2' '3-1' '3-2' '3-3' '3-4'}

for k= 1:length(allnum)
totest = allnum{k};

Zerr(:,k) = Zevaluation(['D:\Data\Raw\EtalonnageZ\MultiZCorrection\02-03-20_Depthograph100X_' totest '.mat'],...
    ['D:\Data\Raw\EtalonnageZ\MultiZCorrection\02-03-20_Depthograph100X_Anti-' totest '.mat'],grid,PLOTperDepth);

end
Zerrme = nanmean(Zerr,2)';

if PLOTperDay
    
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
end

ndepth = ndepth+1;
AllZerrme(:,ndepth) = Zerrme;

clear Ze* P


%% mean over all deptho days

ZerrFinal = nanmean(AllZerrme,2)';

P = polyfit(grid,ZerrFinal,20);
Zfit = polyval(P,grid);

figure
hold on
plot(grid,AllZerrme,'*','handlevisibility','off')
plot(grid,ZerrFinal,'r','linewidth',3)
plot(grid,Zfit,'b','linewidth',1.5)
xlabel('Computed Z (nm)','fontsize',18)
ylabel('Z error (nm)','fontsize',18)
ylim([-400 300])


figure(16654)
hold on
plot(grid,Zfit,'b','linewidth',1.5)