function DeadPlot

close all

warning('off')
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)



% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';

Champ = [5 20];

% dossier de données et sauvegarde
path = 'D:\Data\MatFile';

datenow = datestr(now,'yy-mm-dd');

sfa   = ['D:\Data\Figures' filesep datenow];
sfs   = ['D:\Data\Figures' filesep datenow];

mkdir(sfa)
mkdir(sfs)


% variable
Actot = [];

Cortot = [];
Grouping = [];

cpt = [0 0];



cd(path)
for i = 1:2
    champ = Champ(i);
    Actot{i} = [];
    
load(['ResDead_' num2str(champ) 'mT'])

for k = 1:length(MRI)
    
    D3   = MRI{k}.Dist;
    Dz   = MRI{k}.Dz;
    T    = MRI{k}.Time;
    EF1  = MRI{k}.Elps(:,1);
    EF2  = MRI{k}.Elps(:,2);
    Cort = MRI{k}.Cort - 4432;
    Acti = MRI{k}.Acti;
    Tau  = MRI{k}.Tau;
    
    Cortot = [Cortot Cort];
    Actot{i} = [Actot{i}; Acti];
    Grouping = [Grouping i];
    
    cpt(i) = cpt(i) +1;
    
%     EF = max(EF1,EF2);
%     
%     m1 = EF<1.3;
%     m2 = ((EF>1.3)+(EF<1.6))==2;
%     m3 = (EF>1.6);
    
    
%     figure
%     plot(T,D3,'k')
%     hold on
%     plot(T(m1),D3(m1),'*g')
%     plot(T(m2),D3(m2),'*y')
%     plot(T(m3),D3(m3),'*r')
    
    
end
err(:,i) = std(Actot{i},1);
randX{i} = (rand(cpt(i),1)-0.5)/5;
figure(1)
hold on
plot(randX{i}+i,Cortot(Grouping == i),'o','linewidth',5,'markersize',3)
end





figure(1)
hold on
boxplot(Cortot,Grouping,'plotstyle','traditional','colors','k','widths',0.2,'symbol','');
ylabel('Beads Distance (nm)')
title(['Two beads pinching a dead cell'])
ax = gca;
ax.FontSize = 30;
ax.FontWeight = 'bold';
ax.XTickLabel = [{'5mT'},{'20mT'}] ;


fig = gcf;
saveas(fig,[sfs filesep 'BeadsDistance_dead'],'png')

% figure
% hold on
% title(['Two beads inside a dead cell at ' num2str(champ) 'mT'])
% boxplot(Dztot)
% ylabel('Beads Dz (nm)')
% 
% fig = gcf;
% saveas(fig,[sfs filesep 'BeadsDz_dead_' num2str(champ) 'mT'],'png')


figure(2)
hold on
title(['Two beads pinching a dead cell'])
for i = 1:2
pe = errorbar(Tau,mean(Actot{i},1),err(:,i)/sqrt(cpt(i)),'k');
pe.HandleVisibility = 'off';
plot(Tau,mean(Actot{i},1),'*-','linewidth',2)
end
ax = gca;
ax.FontSize = 30;
ax.FontWeight = 'bold';
xlabel('Tau (s)')
ylabel('Activity (nm)')
legend('5mT','20mT','location','northwest')


fig = gcf;
saveas(fig,[sfa filesep 'Activity_dead'],'png')


close all