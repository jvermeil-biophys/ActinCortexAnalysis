function Data2Peaks_Endo(Inhib,PLOT,alignlvl,...        
    MatfileFolder,figurefolder,DropboxDataFolder)


%% INIT
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)
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


%% Data retrieval and plot
nd = size(date,2);

ic = 1;


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
                %                 try
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
                        Dx = double(MR{kc}.dx);
                        Dy = double(MR{kc}.dy);
                        Dz = double(MR{kc}.dz);
                        F = double(MR{kc}.F);

                        AllDistances{end+1} = D3;
                        
                        D3s = smooth(T,D3,7,'sgolay',3);
                    
                    if MR{kc}.actigood
                        
                     
                        %% detection des pics
                        
                        clear Pks* TpsPks* wPks* pPks*
                        
                        
                        [hPks,TpsPks,wPks,pPks] =  findpeaks(D3s,T,'MinpEakdistance',30);
                        
                        
                        bPks = hPks - 0.80*pPks; % Niveau bas des pics
                        
 
                        IndPks = find(ismember(T,TpsPks) == 1); % Position des maximums dans les vecteurs
                        
                        
                        NpD = length(IndPks); % nombre de pics detecté
                        
                        % variable contenant les infos individuelles sur les pics
                        peakI = struct('prom',pPks,'width',wPks,'height',hPks,'loc',IndPks,'promnorm',pPks./(hPks-pPks));

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
                            peakI.Dy(ip)= mean(Dy(peakI.locs{ip}));
                            peakI.Dz(ip)= mean(Dz(peakI.locs{ip}));
                            
%                             if pPks(ip) > 500
%                                 figure
%                                 subplot(121)
%                                 hold on
% %                                 plot(T(i1-20:i2+20),D3(i1-20:i2+20),'k')
%                                 plot(T(i1:IndPks(ip)),D3(i1:IndPks(ip)),'-*b')
%                                 plot(T(IndPks(ip):i2),D3(IndPks(ip):i2),'-or')
%                                 xlabel('Temps')
%                                 ylabel('D3')
%                                 
%                                 subplot(122)
%                                 hold on
% %                                 plot(T(i1-20:i2+20),D3(i1-20:i2+20),'k')
%                                 plot(Dz(i1:IndPks(ip)),Dy(i1:IndPks(ip)),'-*b')
%                                 plot(Dz(IndPks(ip):i2),Dy(IndPks(ip):i2),'-or')
%                                 xlabel('Dz')
%                                 ylabel('Dy')
%                                 
%                             end
                            
                        end
                        
                        
                        
                        %% tri des pics
                        
                        % conditions d'exclusion d'un pics : grande différence
                        % de hauteur entre deux point (60% de la prominence),
                        % manque de donnée sur 5+ secondes,
                        
                        %                     figure
                        %                     yyaxis left
                        %                                         findpeaks(D3s,T,'annotate','extents','MinPeakProminence',30)
                        %                     fig = gcf;
                        %                     hold on
                        %                     plot(T,D3,'*k','handlevisibility','off')
                        
                        
                        for ip = 1:NpD
                            
                            if (max(abs(diff(D3s(peakI.locs{ip}))))>0.6*pPks(ip))
                                peakI.info{ip} = 'gapH';
                                %                            plot(T(peakI.locs{ip}),D3s(peakI.locs{ip}),'b--','linewidth',2,'handlevisibility','off')
                                
                            elseif (max(diff(T(peakI.locs{ip}))) > 5)
                                peakI.info{ip} = 'GapT';
                                %                            plot(T(peakI.locs{ip}),D3s(peakI.locs{ip}),'r--','linewidth',2,'handlevisibility','off')
                                
                            else
                                peakI.good(ip) = 1;
                                
                            end
                            
                            if not(isempty(peakI.good))
                                NpG = sum(peakI.good); % nombre de pics good/gardé
                            else
                                NpG = 0;
                            end
                            
                            
                        end
                        
                        
                        %                     saveas(fig,[ffp filesep Cell filesep 'PeaksDet_' AcqName],'png')
                        %                     close
                        
                                           peakI.good = logical(peakI.good);
                    else
                        peakI = [];
                        NpD = 0;
                    end
                    
                    
                    
                    
                    
                    SrtD3 = sort(D3);
                    Npts = length(SrtD3);
                    P = (Npts:-1:1)/Npts;

                    ptrcort = P>alignlvl;
                    
                    AlignW = max(SrtD3(ptrcort));
                    
                    CortBot = prctile(D3,10);
                    CortTop = prctile(D3,90);
                                      
                    MedCortW = median(D3);
                    
                    ProbaD3 = SrtD3 - AlignW; 
                    
                    NormD3 = D3./mean(D3);
                    
                    Ptmp = (length(ProbaD3):-1:1)/length(ProbaD3);            
                   
                    
                    if strcmp(Class,'OK')
                            colclass = 'g';
                        elseif strcmp(Class,'G')
                            colclass = 'b';
                        elseif strcmp(Class,'E')
                            colclass = 'm';
                        elseif strcmp(Class,'M')
                            colclass = 'r';
                        end
                        
%                          figure(2325)
%                         hold on
%                         plot(T,D3-mean(D3),colclass)
                    
                    if PLOT && CortBot < 500
                        
                        
                        %% plot avec cortex
                        
                        mkdir([ff filesep 'CurvePlots' filesep Inhib])
                        figure
                        hold on
                        if contains(Inhib,'Endo')
                        plot(S,D3,'-k','linewidth',2,'handlevisibility','off')
                        plot(S,D3,'ok','linewidth',0.1,'markerfacecolor',colclass)
                        xlabel('Image n°') 
                        else
                        plot(T,D3,'-k','linewidth',2,'handlevisibility','off')
                        plot(T,D3,'ok','linewidth',0.1,'markerfacecolor',colclass)
                        xlabel('T (s)')
                        end
                        ylabel('D3 (nm)')
                        if contains(Inhib,'Sep')
                            yyaxis right
                            plot(T,F,'linewidth',2.5)
                            ylabel('Force (pN)')
                        end
                        fig = gcf;
                        saveas(fig,[ff filesep 'CurvePlots' filesep Inhib filesep AcqName '_Zoom'],'png')
                        close
                        
                       
                        
                    end
                    
                    %
                    if NpD == 0
                        NpG = 0;
                    end
                    
 
                    %
                    
                    %% Store peaks data
                    
                    if CortBot < 500
                   
                   
                   [~,i] = min(abs(Ptmp-0.1));            
                   h10pc(ic) = ProbaD3(i);                   

                        
                    MP{ic}.npeaksD = NpD;
                    MP{ic}.npeaksG = NpG;
                    MP{ic}.CortBot = CortBot; 
                    MP{ic}.CortTop = CortTop; 
                    MP{ic}.MedCortW = MedCortW;  
                    MP{ic}.h10pc = h10pc(ic);                     
                    MP{ic}.time = T;
                    MP{ic}.dx = Dx;
                    MP{ic}.dz = Dz;
                    MP{ic}.dy = Dy;
                    MP{ic}.dist = D3;
                    MP{ic}.normdist = NormD3;
                    MP{ic}.dtwD3 = (D3-CortBot)./(mean(D3-CortBot));
                    MP{ic}.peaks = peakI;
                    MP{ic}.name = AcqName;
                    MP{ic}.class = Class;
                    MP{ic}.CultureTag = CultureTag;
                    
                    
                    
                    ic = ic+1;
                    end
                    
                    clear NpD NpG T Dz Dy D3 peakI AcqName *Pks
                end
                %                 catch
                %                         pause
                %                 end
            end
        end
    end
    
    %% Saving data
    Sname = ['D2P_' date{kd} '_' Inhib];
    Sname2 = ['Hxxpc_' date{kd} '_' Inhib];
    
    if exist('MP','var')     
        savename2 = [sf filesep 'D2P' filesep Sname2];
        save(savename2,'h10pc')
        
        
        savename = [sf filesep 'D2P' filesep Sname];
        save(savename,'MP')
        savename = [scf filesep datenow filesep 'D2P' filesep Sname];
        save(savename,'MP')
        fprintf(['File saved : \n' Sname '.mat \n\n']);
        clear MP h10pc
        ic = 1;
    end
    
    
    
    
    
end

save('D:\Data\MatFile\Blebbi_AllD.mat','AllDistances')

end

