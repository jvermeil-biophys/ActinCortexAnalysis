function Peaks2Stats_J(date,Cell,PLOT,alignlvl,matfilefolder,figurefolder,resfolderbis)

%% INIT

set(0,'DefaultFigureWindowStyle','docked')

warning('off','all')


% data folder
df = [matfilefolder filesep 'D2P'];

% save folder
sf   = matfilefolder;

datenow = datestr(now,'yy-mm-dd');

mkdir([sf filesep 'P2S'])

if nargin == 6 && not(isempty(resfolderbis))
    scf = resfolderbis;
    
    mkdir([scf filesep datenow filesep 'P2S'])
end

% curve folder
ff = [figurefolder filesep datenow];

ffp = [ff filesep 'PointsDist' filesep Cell];

ffp2 = [ff filesep 'TemporalAnalysis'];

mkdir(ffp)

mkdir(ffp2)




nd = size(date,2);

Ctot = [];
Mins = [];
Maxs = [];



%% Data retrieval and peaks sorting

% variable globale
Ttot = 0; % temps total d'exp
NpTot = 0; % nb tot de pics


P   = []; % proeminence des pics
PN   = []; % proeminence des pics normalisé par le cortex local
H   = []; % hauteur des pics
PeakCort   = []; % épaisseur cortex par pics
W   = []; % largeur des pics

D3DataPts = []; % tt les points de données
CBotTot = [];
CTopTot = [];
CActTot = [];
CActBegTot = [];
CActEndTot = [];
CAssTot = [];
Timeinds = [];
CstdTot = [];
CmeanTot = [];
MedCwTot = [];
nbBilleInt = [];

Distances = {};
YDistances = {};
XDistances = {};
Times = {};
Names = {};

CTagList = {};

Data = struct('type',Cell);
NC = 0;
NC2 = 0;

H10pc= [];

AllD3 = [];
AllD3N = [];

dtwD3s = {};

for kd=1:nd
    
    datafile0=['Hxxpc_' date{kd} '_' Cell '.mat'];
    E0 = exist([df filesep datafile0]);
    
    if E0 == 2
        
        load([df filesep datafile0])
        fprintf(['File Loaded : \n' datafile0 '\n'])
        
        H10pc = [H10pc h10pc];
        
        clear h10pc
        
    end
end


nrem = max([1 round(0.05*length(H10pc))]);

Sh = sort(H10pc);

listMedianThickness = [];
listFluctuations = [];

if length(H10pc)<600
    Hmin = -Inf;
    Hmax= Inf;
else
    Hmax = Sh(end-nrem+1)-0.0001;
    Hmin = Sh(nrem) +0.0001;
end

DataTime = {};


for kd=1:nd
    datafile=['D2P_' date{kd} '_' Cell '.mat'];
    
    
    % verification de l'existence du fichier
    E = exist([df filesep datafile]);
    
    if E == 2
        
        
        load([df filesep datafile])
        fprintf(['File loaded : \n' datafile '. \n']);
        
        % nombre de courbe
        s=size(MP); nc=s(2);
        
        
        for ni = 1:nc
            
            NC = NC + 1;
            % names
            try
                AcqName = MP{ni}.name;
            catch
                themall
            end
            CTag = MP{ni}.CultureTag;
            
            Class = MP{ni}.class;
            if strcmp(Class,'OK')
                colclass = 'g';
            elseif strcmp(Class,'G')
                colclass = 'b';
            elseif strcmp(Class,'E')
                colclass = 'm';
            elseif strcmp(Class,'M')
                colclass = 'r';
            end
            
            CortBot = MP{ni}.CortBot;
            CBotTot = [CBotTot CortBot];
            
            CortTop = MP{ni}.CortTop;
            CTopTot = [CTopTot CortTop];
            
            CortActi = MP{ni}.CortActi;
            CActTot = [CActTot CortActi];
            
            CortActiBeg = MP{ni}.CortActiBeg;
            CActBegTot = [CActBegTot CortActiBeg];
            
            CortActiEnd = MP{ni}.CortActiEnd;
            CActEndTot = [CActEndTot CortActiEnd];
            
            CortAssym = MP{ni}.CortAssym;
            CAssTot = [CAssTot CortAssym];
            
            posptr= strfind(AcqName,'_C')+2;
            
            TimeInd = str2num(AcqName(posptr:posptr+1));
            
            if isempty(TimeInd)
                TimeInd = str2num(AcqName(posptr));
            end
            
            Timeinds = [Timeinds TimeInd];
            
            CTagList = {CTagList{:} CTag};
            
            MedCortW = MP{ni}.MedCortW;
            MedCwTot = [MedCwTot MedCortW];
            
            nbBI = MP{ni}.nbBI;
            nbBilleInt = [nbBilleInt nbBI];
            
            
            
            
            Data.cell{NC}.name = AcqName;
            Data.cell{NC}.CTag = CTag;
            D3 = MP{ni}.dist;
            DY = MP{ni}.dy;
            DX = abs(MP{ni}.dx);
            NormD3 = MP{ni}.normdist;
            
            
            
            
            Cmean = mean(D3);
            CmeanTot = [CmeanTot Cmean];
            Cstd = std(D3);
            CstdTot = [CstdTot Cstd];
            
            dtwD3s{NC} = MP{ni}.dtwD3;
            
            % temps et dist data
            T = MP{ni}.time;
            
            AllD3 =  [AllD3 D3'];
            AllD3N =  [AllD3N NormD3'];
            
            
            dT = diff(T);
            ptrdT = dT<3;
            
            Ttot = Ttot + sum(dT(ptrdT));
            Data.cell{NC}.time = T;
            Data.cell{NC}.Ttot = T(end) - T(1);
            
            NpTot = NpTot + MP{ni}.npeaksG;
            
            Data.cell{NC}.NpC = MP{ni}.npeaksG;
            
            peakI = MP{ni}.peaks;
            
            
            
            if ~isempty(peakI)
                Data.cell{NC}.NpC = sum(peakI.height(peakI.good)>60);
                
                P  = [P  peakI.prom(peakI.good)'];
                H  = [H  peakI.height(peakI.good)'];
                W  = [W  peakI.width(peakI.good)'];
                PN = [PN peakI.promnorm(peakI.good)'];
                
                PeakCort = [PeakCort repmat(CortBot,1,sum(peakI.good))];
                
                
                Data.cell{NC}.prom = peakI.prom(peakI.good);
                Data.cell{NC}.height = peakI.height(peakI.good);
                Data.cell{NC}.width = peakI.width(peakI.good);
            end
            
            
            
            
            
            if (MP{ni}.h10pc < Hmax) &&  (MP{ni}.h10pc > Hmin)
                
                            
                NC2 = NC2 +1;
                
                % % % Autocorrelation
                
                %
                if (T(end) - T(1) > 300) &&  ((T(end)-T(1))/(length(T)-1) <= 1)
                    
                   
                    [C,L] = AUTOCORR(T,D3,2);
                    
                    ptr = L<301;
                    
                    C = C(ptr);
                    L = L(ptr);
                    
                    Ctot = [Ctot C];
                    
                    
                    [valmn, locmn] = findpeaks(-C,'NPeaks',1,'minpeakprominence',0.01);
                    
                    [valmx, locmx] = findpeaks(C,'NPeaks',1,'minpeakprominence',0.01);
                    
                    if PLOT
                        figure
                        hold on
                        title([Cell ' : Autocorrelation of ' AcqName])
                        plot(L,C,'b-*')
                        xlabel('Lag (sec)')
                        ylabel('AutoCorrelation')
                        
                        plot(L(locmn),C(locmn),'*r','linewidth',5)
                        plot(L(locmx),C(locmx),'*g','linewidth',5)
                        
                        xlim([0 350])
                        ylim([-0.8 1.05])
                        
                        fig = gcf;
                        saveas(fig,[ffp2 filesep 'Autocorr_' Cell ' - ' AcqName '.png'],'png')
                        saveas(fig,[ffp2 filesep 'Autocorr_' Cell ' - ' AcqName '.fig'],'fig')
                        close
                    end
                    Mins = [Mins L(locmn)];
                    Maxs = [Maxs L(locmx)];
                    
                    
                end
                
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                
                
                
                Distances{end+1} = D3;
                YDistances{end+1} = DY;
                XDistances{end+1} = DX;
                Times{end+1} = T;
                Names{end+1} = AcqName;
                
                DataTime{NC2} = [T D3];
                
                D3Centered = D3' - prctile(D3,100*(1-alignlvl));
                
                
                
                D3DataPts = [D3DataPts D3Centered];
                D3Data{NC2} = sort(D3Centered);
                if CortBot>0
                    D3DataN{NC2} = sort(NormD3);
                else
                    D3DataN{NC2} = NaN;
                end
                
                Data.cell{NC2}.dist = D3;
                Data.Npts(NC2) = length(D3);
                
                
                
                
                if PLOT
                    
                    plotable = sort(D3Centered);
                    Ptmp = (length(plotable):-1:1)/length(plotable);
                    figure(1)
                    hold on
                    plot(plotable,Ptmp*100,'color',colclass)
                    ylabel('P(h)')
                    xlabel('Height above cortex')
                    
                                       
                    
                    title(AcqName)
                    
                end
                
                
                
                
            else
                
                AcqName = MP{ni}.name;
                T = MP{ni}.time;
                D3 = MP{ni}.dist;
                CortBot = MP{ni}.CortBot;
                
                D3Centered = D3' - CortBot;
                
                if PLOT
                    
                    plotable = sort(D3Centered);
                    Ptmp = (length(plotable):-1:1)/length(plotable);
                    figure(1)
                    hold on
                    plot(plotable,Ptmp*100,colclass)
                    ylabel('P(h)')
                    xlabel('Height above cortex')
                    
                    title(AcqName)
                end
                
                
                
            end
            
            if PLOT
                
                plotable = sort(D3Centered);
                
                Ptmp = (length(plotable):-1:1)/length(plotable);
                
                figure(1)
                hold on
                plot(plotable,Ptmp*100,'color',colclass)
                ylabel('P(h)')
                xlabel('Height above cortex')
                
                figure(2)
                hold on
                xlim([-1 2])
                listMedianThickness = [listMedianThickness MedCortW];
                plot(rand(1),MedCortW,'ok','markerfacecolor',colclass,'markersize',10)
                ylabel('Cortex Thickness (nm)')
                
                figure(3)
                hold on
                xlim([-1 2])
                listFluctuations = [listFluctuations CortActi];
                plot(rand(1),CortActi,'ok','markerfacecolor',colclass,'markersize',10)
                ylabel('Thickness variation amplitude (nm)')             
            end
        end
    end
end

figure(2)
hold on
plot([0 1],[median(listMedianThickness) median(listMedianThickness)],'r--','linewidth',2)

figure(3)
hold on
plot([0 1],[median(listFluctuations) median(listFluctuations)],'r--','linewidth',2)
%% Storing Data
MS.freq = NpTot/Ttot*60; % pics / min
MS.timetot = Ttot;
MS.proms = P;
MS.promsN = PN;
MS.heights = H;
MS.peaksCortW = PeakCort;
MS.D3Data = D3Data;
MS.D3DataN = D3DataN;
MS.AllD3 = AllD3;
MS.AllD3N = AllD3N;
MS.D3DataPts = sort(D3DataPts);
MS.widths = W;
MS.data = Data;
MS.CortBots = CBotTot;
MS.CortBot = mean(CBotTot);
MS.CortTops = CTopTot;
MS.CortTop = mean(CTopTot);
MS.CortActisBeg = CActBegTot;
MS.CortActisEnd = CActEndTot;
MS.CortActis = CActTot;
MS.CortAssyms = CAssTot;
MS.CTags = CTagList;
MS.MedCortWs = MedCwTot;
MS.nbBI = nbBilleInt;
MS.Cortmeans = CmeanTot;
MS.Cortsdts = CstdTot;
MS.MedCortW = mean(MedCwTot);
MS.H10pc = H10pc(H10pc>Hmin&H10pc<Hmax);

%% Saving data

save([sf filesep 'P2S' filesep 'TimeDatas_' Cell],'Distances','YDistances','Times','Names')
save([sf filesep 'P2S' filesep 'P2S_' Cell],'MS')

if nargin == 6 && not(isempty(resfolderbis))
save([scf filesep datenow filesep 'P2S' filesep 'P2S_'  Cell],'MS')
end
fprintf(['File saved : \nP2S_' Cell '.mat\n\n'])
