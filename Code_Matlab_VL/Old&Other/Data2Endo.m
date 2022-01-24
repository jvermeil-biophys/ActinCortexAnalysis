function Data2Endo(Inhib,PLOT,...
    MatfileFolder,figurefolder,DropboxDataFolder)



%% INIT
set(groot,'DefaultFigureWindowStyle','normal')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot,'DefaultAxesFontSize',30)
warning('off','all')

% data folder
df = [MatfileFolder filesep 'V2D'];

% save folder
date=get_dates();

sf   = MatfileFolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'D2P'])

scf = DropboxDataFolder;

mkdir([scf filesep datenow filesep 'D2P'])



% pour les figures
ff = [figurefolder filesep datenow];
mkdir([ff filesep 'CurvePlots' filesep Inhib])

%% Data retrieval and plot
nd = size(date,2);

ic = 0;


AllDistances = {};

for kd=1:nd
    for km=1:5
        
        %% Data retrieval and peaks detection
        datafile=['V2D_' date{kd} '_M' num2str(km) '_' Inhib '.mat'];
        
        
        % verification de l'existence du fichier
        E = exist([df filesep datafile],'file');
        
        if E == 2
            
            
            load([df filesep datafile],'MR');
            fprintf(['File loaded : \n' datafile '. \n\n']);
            % nombre de courbe
            s=size(MR); nc=s(2);
            
            
            
            for kc=1:nc
                if not(isempty(MR{kc}))
                    
                    % Name of the experiment and image
                    
                    AcqName = MR{kc}.name;
                    Class = MR{kc}.class;
                    
                    expdate = AcqName(1:8); % date of experiment
                    
                    CultureTag = findcultureday(expdate);
                    
                    % recuperation des donnes de positions (en nm)
                    T = MR{kc}.time;
                    S = MR{kc}.S;
                    D3 = double(MR{kc}.D3);
                       
                    ic = ic+1;
                    
                    figure
                    title('Where does the endocytosis happen ? Click on it !')
                    plot(S,D3,'-o')
                    
                    [x y] = ginput(1);

                    close
                    Dbefore{ic} = D3(S<x);
                    Dafter{ic} = D3(S>x);
                    
                    MedBefore(ic) = median(Dbefore{ic});
                    MedAfter(ic) = median(Dafter{ic});
                    
                    StdBefore(ic) = std(Dbefore{ic});
                    StdAfter(ic) = std(Dafter{ic});
                    
                    name{ic} = AcqName;
                    
%                     
%                     figure
%                     hold on                        
%                         plot(S,D3,'-k','linewidth',2,'handlevisibility','off')
%                         plot(S,D3,'ok','linewidth',0.1,'markerfacecolor','b')
%                         xlabel('Image n°') 
%                         ylabel('Bead separation (nm)')
%                         
%                     
%                     fig = gcf;
%                         saveas(fig,[ff filesep 'CurvePlots' filesep Inhib filesep AcqName '_Image.png'],'png')
%                         close
%                         
%                         figure
%                     hold on                        
%                         plot(T,D3,'-k','linewidth',2,'handlevisibility','off')
%                         plot(T,D3,'ok','linewidth',0.1,'markerfacecolor','b')
%                         xlabel('Time (s)') 
%                         ylabel('Bead separation (nm)')
%                         
%                     
%                     fig = gcf;
%                         saveas(fig,[ff filesep 'CurvePlots' filesep Inhib filesep AcqName '_Time'],'png')
%                         close
                        
                    
                    
                end
            end
        end
    end
end

figure
subplot(3,4,[1,2,5,6,9,10])
hold on
title('Median bead separation before and after endocytosis')
Spread_Julien(MedBefore,1,0.7,'b','o')
Spread_Julien(MedAfter,2,0.7,'r','o')
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Before','After'};
ylabel('Median Bead Separation (nm)')
xlim([0 3])
listsub = [3 4 7 8 11 12];
ax = gca;
ax.FontSize = 18;



for kc = 1:6
subplot(3,4,listsub(kc))
hold on
title(name{kc})
    sortD = sort(Dbefore{kc}-MedBefore(kc));
Ptmp = (length(sortD):-1:1)/length(sortD);
              
                plot(sortD,Ptmp,'b')
                
                    sortD2 = sort(Dafter{kc}-MedAfter(kc));
Ptmp2 = (length(sortD2):-1:1)/length(sortD2);
              
                plot(sortD2,Ptmp2,'r')
                
                ax = gca;
                ax.FontSize = 10;
                                                
legend('Before','After')
ylabel('P(h)')
xlabel('Height above median bead separation')
end


figure
hold on
title('Relative variablity of bead separation')
Spread_Julien(StdBefore./StdAfter,1,0.7,'k','o')
ylabel('StdBefore/StdAfter')



