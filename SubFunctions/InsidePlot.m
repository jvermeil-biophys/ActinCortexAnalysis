function InsidePlot(Champ,pc)

% close all

warning('off')
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)



% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';

% Champ = [5 20];

% dossier de données et sauvegarde
path = 'D:\Data\MatFile';

datenow = datestr(now,'yy-mm-dd');

sfa   = ['D:\Data\Figures' filesep datenow];
sfs   = ['D:\Data\Figures' filesep datenow];

mkdir(sfa)
mkdir(sfs)


% variable
Actot = [];
ActStd = [];
Cortot = [];
Grouping = [];

cpt = [0 0];

col = [0 0 1; 1 0 0];

cd(path)
for i = 1:length(Champ)
    for ii = 1:length(pc)
        champ = Champ(i);
        Actot{i,ii} = [];
        
        E = exist(['ResInside_' num2str(champ) 'mT_' num2str(pc(ii)) 'pc.mat'],'file')
        
        if E
        
        load(['ResInside_' num2str(champ) 'mT_' num2str(pc(ii)) 'pc.mat'])
        
        for k = 1:length(MRI)
            
            D3   = MRI{k}.Dist;
            Dz   = MRI{k}.Dz;
            T    = MRI{k}.Time;
            Cort = MRI{k}.Cort - 4432;
            Acti = MRI{k}.Acti;
            Tau  = MRI{k}.Tau;
            Stackname = MRI{k}.stackname;
            
            if ((ii == 1)&&(std(D3) < 700)) || ((std(D3) < 700) && not((strcmp(Stackname(1:11),'29-05-18_M1')||...
                    strcmp(Stackname(1:11),'04-06-18_M1')||strcmp(Stackname(1:11),'07-06-18_M2'))))
                
                Dm = round(D3./10)*10;              
                    md = mode(Dm);                                    
                    Cort = md;
            Cortot = [Cortot Cort];
            Actot{i,ii} = [Actot{i,ii}; Acti];
            ActStd = [ActStd std(D3)];
            Grouping = [Grouping ii];
            
            cpt(ii) = cpt(ii) +1;
            end
            
            
        end
        if size(Acti,1)>1
            err(i,:) = std(Acti);
        else
            err(i,:) = zeros(1,100);
        end
        end
    end
    
end





figure(1)
hold on
for ii = 1:length(pc)
    randX{ii} = (rand(cpt(ii),1)-0.5)/5;
    plot(randX{ii}+ii,Cortot(Grouping == ii),'o','linewidth',5,'markersize',3,'color',col(ii,:))
end
boxplot(Cortot,Grouping,'plotstyle','traditional','colors','k','widths',0.2,'symbol','');
ylabel('Beads Distance (nm)')
title(['Two beads inside a cell'])
ax = gca;
ax.FontSize = 30;
ax.FontWeight = 'bold';
ax.XTickLabel = [{[num2str(pc(1)) '% Optiprep']},{[num2str(pc(2)) '% Optiprep']}] ;


fig = gcf;
saveas(fig,[sfs filesep 'BeadsDistance_Inside_Boxplot'],'png')



figure(2)
hold on
bar([0 7.5], [mean(Cortot(Grouping == 1)) mean(Cortot(Grouping == 2))])
errorbar([0 7.5], [mean(Cortot(Grouping == 1)) mean(Cortot(Grouping == 2))],[std(Cortot(Grouping == 1))/sqrt(sum(Grouping ==1)) std(Cortot(Grouping == 2))/sqrt(sum(Grouping ==1))],'.k')
ylabel('Beads Distance (nm)')
title(['Two beads inside a cell'])
ax = gca;
ax.FontSize = 30;
ax.FontWeight = 'bold';
ax.XTick = [0 7.5];
ax.XTickLabel = [{[num2str(pc(1)) '% Optiprep']},{[num2str(pc(2)) '% Optiprep']}] ;


fig = gcf;
saveas(fig,[sfs filesep 'BeadsDistance_Inside_Barplot'],'png')
% close all

end