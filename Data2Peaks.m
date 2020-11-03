function Data2Peaks(Inhib,PLOT,...
    MatfileFolder,figurefolder,resfolderbis)


%% INIT
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)
warning('off','all')

% data folder
df = [MatfileFolder filesep 'V2D'];

% save folder
date=get_dates();

sf   = MatfileFolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'D2P'])

if nargin ==5 && not(isempty(resfolderbis))
scf = resfolderbis;

mkdir([scf filesep datenow filesep 'D2P'])

end

% pour les figures
ff = [figurefolder filesep datenow];

CommonTimes = {};
CurveOnes = {};
CurveTwos = {};
NameOnes = {};
NameTwos = {};
Lags = {};
Corrs = {};

%% Data retrieval and plot
nd = size(date,2);

ic = 1;

AllMeanF = [];
AllL = [];
AllCortActi = [];

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
            
            
            
            % trouvé les cellules avec deux interfaces
            codelist = {};
            
            for k = 1:nc
                if not(isempty(MR{k}))
                    codelist = {codelist{:} MR{k}.name(10:end-4)};
                    
                end
            end
            
            if contains(Inhib,'Outside')
                ptrDI = [];
            else
                ptrDI = find(contains(codelist,'-')); % nom dans la liste avec une possible double interface
                
            end
            
            
            
            for kc=1:nc
                %                 try
                if not(isempty(MR{kc}))
                    
                    % Name of the experiment and image
                    
                    AcqName = MR{kc}.name;
                    code = AcqName(10:end-4);
                    Class = MR{kc}.class;
                    
                    if contains(Inhib,'Outside')
                        nbBI = 0;
                    else
                        nbBI = MR{kc}.nbBI;
                    end
                    
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
                    
                    meanF = mean(F);
                    
                    % Si cette courbe est dans ptrDI on cherche la
                    % courbe associéé
                    
%                     if ismember(kc,ptrDI)
%                         
%                         ptrDI = ptrDI(ptrDI~=kc);
%                         
%                         pos = strfind(code,'-')-1;
%                         
%                         kassoc = find(contains(codelist,code(1:pos)));
%                         
%                         kassoc = kassoc(kassoc~=kc);
%                         
%                         ptrDI = ptrDI(ptrDI ~= kassoc);
%                         
%                         % on fait la correlation entre les deux courbe
%                         
%                         D3other = MR{kassoc}.D3;
%                         Tother = MR{kassoc}.time;
%                         AcqNameOther = MR{kassoc}.name;
%                         
%                         CommonTime = max([T(1),Tother(1)]):0.5:min([T(end),Tother(end)]);
%                         
%                         D3int = interp1(T,D3,CommonTime);
%                         D3otherint = interp1(Tother,D3other,CommonTime);
%                         
%                         [c,lags] = xcorr(D3int-mean(D3int),D3otherint-mean(D3otherint)/mean(D3otherint),'normalized');
%                         
%                         
%                         figure
%                         hold on
%                         subplot(211)
%                         hold on
%                         title(Inhib)
%                         plot(CommonTime,D3int,'b-*')
%                         plot(CommonTime,D3otherint,'r-*')
%                         legend(AcqName,AcqNameOther,'interpreter','none')
%                         xlabel('Time (s)')
%                         ylabel('Thickness (nm)')
%                         subplot(212)
%                         plot(lags/2,c,'-o')
%                         grid on
%                         ylim([-0.5 0.5])
%                         xlabel('Lag (s)')
%                         ylabel('Cross correlation')
%                         
%                         
%                         
%                         CommonTimes{end+1} = CommonTime;
%                         CurveOnes{end+1} = D3int;
%                         CurveTwos{end+1} = D3otherint;
%                         NameOnes{end+1} = AcqName;
%                         NameTwos{end+1} = AcqNameOther;
%                         Lags{end+1} = lags/2;
%                         Corrs{end+1} = c;
%                         
%                         
%                     end
%                     
                    
                    
                                       
                    
                    AllDistances{end+1} = D3;
                    
                    D3s = smooth(T,D3,7,'sgolay',3);
                    
                    if MR{kc}.actigood
                        
                        
                        %% detection des pics
                        
                        clear Pks* TpsPks* wPks* pPks*
                        
                        
                        [hPks,TpsPks,wPks,pPks] =  findpeaks(D3s,T,'MinPeakProminence',8);
                        
                        
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
                    
                    ptrcort = P>0.5;
                    
                    AlignW = max(SrtD3(ptrcort));
                    
                    CortBot = prctile(D3,10);
                    CortTop = prctile(D3,90);
                    
                    CortActi = CortTop-CortBot;
                    AllCortActi = [AllCortActi CortActi];
                    
                    %% fluctuation par moitié de courbe quand + de 10minutes
                    
                    if T(end) - T(1) > 600
                       
                        Tmid = T(1) + (T(end)-T(1))/2;
                        
                        CortActiBeg = prctile(D3(T<Tmid),90)-prctile(D3(T<Tmid),10);
                        CortActiEnd = prctile(D3(T>Tmid),90)-prctile(D3(T>Tmid),10);
                        
                        ActiBegEnd = 1;
                        
                    else
                        
                        ActiBegEnd = 0;
                        
                    end
                    
                    MedCortW = median(D3);
                    
                    CortAssym = ((CortTop-MedCortW)-(MedCortW-CortBot))/CortActi;
                                        
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
                    
                    if PLOT && (CortBot < 500 || contains(Inhib,'RPE1'))
                        
                        
                        %% plot avec cortex
                        
                        mkdir([ff filesep 'CurvePlots' filesep Inhib])
                        figure
                        hold on
                        if contains(Inhib,'Sep')
                            ptr = find(F < 1);
                            patch([T(ptr(1))-0.7 T(ptr(end))+10 T(ptr(end))+10 T(ptr(1))-0.7],[0 0 3000 3000],[0.8 0.8 0.8],'facealpha',0.4,'edgecolor','none')
                            
                            
                        end
                        if contains(Inhib,'Endo')||contains(Inhib,'Div')
                            plot(S,D3,'-k','linewidth',2,'handlevisibility','off')
                            plot(S,D3,'ok','linewidth',0.1,'markerfacecolor',colclass)
                            xlabel('Image n°')
                            namename = 'Images';
                        else
                            plot(T,D3,'-k','linewidth',2,'handlevisibility','off')
                            plot(T,D3,'ok','linewidth',0.1,'markerfacecolor',colclass)
                            xlabel('T (s)')
                            namename= 'Time';
                        end
                        ylabel('D3 (nm)')
                        
                        %                         plot([0 T(end)],[MedCortW MedCortW],'--k')
                        %
                        fig = gcf;
                        saveas(fig,[ff filesep 'CurvePlots' filesep Inhib filesep AcqName '_' namename],'png')
                        close
                        
                        
                        
                    end
                    
                    %
                    if NpD == 0
                        NpG = 0;
                    end
                    
                    
                    %
                    
                    %% Store peaks data
                    
                    if (CortBot < 500 || contains(Inhib,'RPE1'))
                        
                        
                        [~,i] = min(abs(Ptmp-0.1));
                        h10pc(ic) = ProbaD3(i);
                        
                        
                        MP{ic}.npeaksD = NpD;
                        MP{ic}.npeaksG = NpG;
                        MP{ic}.CortBot = CortBot;
                        MP{ic}.CortTop = CortTop;
                        MP{ic}.CortActi = CortActi;
                        MP{ic}.ActiBegEnd = ActiBegEnd;
                        if ActiBegEnd
                            MP{ic}.CortActiBeg = CortActiBeg;
                            MP{ic}.CortActiEnd = CortActiEnd;
                        else
                            MP{ic}.CortActiBeg = NaN;
                            MP{ic}.CortActiEnd = NaN;                            
                        end
                        MP{ic}.CortAssym = CortAssym;
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
                        MP{ic}.nbBI = nbBI;
                        MP{ic}.CultureTag = CultureTag;
                        MP{ic}.meanF = meanF;
                        
                            AllMeanF = [AllMeanF meanF];
                        
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
%
% save('D:\Data\MatFile\Blebbi_AllD.mat','AllDistances')




% sauvegarde cellule double interfaces
save([sf filesep 'D2P' filesep 'DoubleInt_' Inhib],'CommonTimes','CurveOnes', 'CurveTwos',...
    'NameOnes', 'NameTwos', 'Lags', 'Corrs')





end

