function Deadplot()
% close all

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

H = [];

cd(path)
for i = 1:length(Champ)
    champ = Champ(i);
    Actot{i} = [];
    Cortot{i} = [];
    H{i} = [];
    
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
        Stackname = MRI{k}.stackname;
        
        Cortot{i} = [Cortot{i} Cort];
        Actot{i} = [Actot{i}; Acti];
        Grouping = [Grouping i];
        
        cpt(i) = cpt(i) +1;
        
        
        % peak detection
        D3s = smooth(D3,7,'sgolay',3);
        
        [Pks,TpsPks,wPks,pPks] =  findpeaks(D3s,T,'MinPeakProminence',0);
%         findpeaks(D3s,T,'annotate','extents','MinPeakProminence',0)
        
        bPks = Pks - 0.95*pPks; % Niveau bas des pics
        
        hPks = Pks -0.2*pPks; % niveau de mesure de la largeur
        
        IndPks = find(ismember(T,TpsPks) == 1); % Position des maximums dans les vecteurs
        
        %                     figure
        %                     plot(wPks,'*')
        %                     hold on
        
        for ii = 1:length(wPks)
            
            
            ind1 = IndPks(ii);
            
            while (ind1>1) && (D3s(ind1) > hPks(ii))
                ind1 = ind1 -1;
            end
            T1 = T(ind1);
            
            ind2 = IndPks(ii);
            while (ind2<length(D3s)) && (D3s(ind2) > hPks(ii))
                ind2 = ind2 + 1;
            end
            T2 = T(ind2);
            
            wPks(ii) = T2-T1;
        end
        
        %                     plot(wPks,'*')
        
        NpD = length(IndPks); % nombre de pics detecté
        
        %                     close
        
        % variable contenant les infos individuelles sur les pics
        peakI = struct('height',pPks,'width',wPks,'top',Pks,'bot',bPks,'loc',IndPks);
        
        peakI.locs = {}; % pour les indices de tout les pics
        
        peakI.good = zeros(1,NpD); % validité du pics (pour le tri)
        peakI.info = cell(1,NpD); % info en cas de non validité
        
        % on récupère les positions complètes de chaque pics
        
        for ip = 1:NpD
            
            i1 = IndPks(ip); % avant
            i2 = i1; % après
            
            while (i1>1) && (D3s(i1) > bPks(ip))
                i1 = i1 - 1;
            end
            
            while (i2<length(D3s)) && (D3s(i2) > bPks(ip))
                i2 = i2 + 1;
            end
            
            peakI.locs{ip} = i1:i2;
            
        end
        
        
        %% tri des pics
        
        % conditions d'exclusion d'un pics : grande différence
        % de hauteur entre deux point,
        % manque de donnée sur 5+ secondes,  trop grand
        % dz
        
        
        for ip = 1:NpD
            
            if (max(abs(diff(D3s(peakI.locs{ip}))))>0.8*pPks(ip))
                peakI.info{ip} = 'gapH';
                
            elseif (max(diff(T(peakI.locs{ip}))) > 5)
                peakI.info{ip} = 'GapT';
                
            elseif (max(Dz(peakI.locs{ip})) > 900)
                peakI.info{ip} = 'Dz';
                
            else
                peakI.good(ip) = 1;
                
            end
            
            if not(isempty(peakI.good))
                NpG = sum(peakI.good); % nombre de pics good/gardé
            else
                NpG = 0;
            end
            
            
        end
        
        if NpD == 0
            NpG = 0;
        end
        
        peakI.good = logical(peakI.good);
        
        
        
  
          
          
        H{i} = [H{i} peakI.height(peakI.good)'];
        
    
            clear peakI 
             
    end
    err(i,:) = std(Acti);
    

end





figure
hold on
for i = 1:length(Champ)
    randX = (rand(cpt(i),1)-0.5)/5;
    plot(randX+i,Cortot{i},'o','linewidth',5,'markersize',3)
end
boxplot([Cortot{:}],Grouping,'plotstyle','traditional','colors','k','widths',0.2,'symbol','');
ylabel('Beads Distance (nm)')
title(['Experiment on a dead cell - Cortex'])
ax = gca;
ax.FontSize = 30;
ax.FontWeight = 'bold';

fig = gcf;
saveas(fig,[sfs filesep 'BeadsDistance_Dead'],'png')
close


for i = 1:length(Champ)
figure
hold on
title(['Dead cell peaks'])
histogram(H{i})
ylabel('Count')
xlabel('Height (nm)')

fig = gcf;
saveas(fig,[sfs filesep 'DeadCellsPeaks_' num2str(Champ(i)) 'mT'],'png')
close

% 
% figure
% hold on
% title(['Experiment on a dead cell - Activity'])
% 
%     pe = errorbar(Tau,mean(Actot{i},1),err(i,:)/sqrt(cpt(i)),'k');
%     pe.HandleVisibility = 'off';
%     plot(Tau,mean(Actot{i},1),'*-','linewidth',2)
% 
% ax = gca;
% ax.FontSize = 30;
% ax.FontWeight = 'bold';
% xlabel('Tau (s)')
% ylabel('Activity (nm)')
% legend('5mT','location','northwest')
% 
% 
% fig = gcf;
% saveas(fig,[sfs filesep 'Activity_Dead'],'png')
% close
% end

% close all

end