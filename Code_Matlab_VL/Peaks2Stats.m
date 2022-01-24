function Peaks2Stats(specif,PLOT,alignlvl,resfolder,figurefolder)

%
% Peaks2Stats is used for analyzing constant field curves. It load the data
% from Data2Peaks and pulls together data from all the cells. It also
% performs curves autocorrelation and activity curves.
%
% Peaks2Stats(specif,PLOT,alignlvl,matfilefolder,figurefolder)
%
% specif : experimental condiion (ex: 'Dicty_WT')
% PLOT : 1 = plot autocorr and activity, 0 = no plots
% alignlvl : alignement level for activty curve, usually 0.5 (50%) or 0.9 (90%)
% resfolder : path to the matlab data, should be MatFileFolder in main
% figurefolder : path to figure saving folder, should be FigureFolder
% in main
%


% figure options
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)
warning('off','all')

% data folder
df = [resfolder filesep 'D2P'];

% data save folder
sf   = resfolder;
mkdir([sf filesep 'P2S'])

% figure folder
datenow = datestr(now,'yy-mm-dd'); % today
ff = [figurefolder filesep datenow];
mkdir(ff)


%% Data retrieval and peaks sorting

% global variables
Ttot = 0; % total time of experiment
NpTot = 0; % total number of points


P   = []; % peaks prominence
H   = []; % peaks height
PeakCort   = []; % cortex thickness par cell for each peak
W   = []; % peak width

AllD3 = []; % all datapoints pulled together
D3DataPts = []; % all datapoints centered on alignlvl
CBotTot = []; % cortex bottom level (10%)
CTopTot = []; % cortex top level (90%)
CActTot = []; % cortex fluctuations amplitude
CActBegTot = []; % cortex fluctuations amplitude for first half of long curves
CActEndTot = []; % cortex fluctuations amplitude for second half of long curves
CAssTot = []; % cortex fluctations assymetry
MedCTot = []; % median cortex thickness
ActiCurves = {}; % activity curves

PeaksDet = 0; % number of peaks detected
PeaksKept = 0; % number of peaks kept

% autocorr data
Mins = []; % First minimum of autocorrelation
Maxs = []; % first maximum of autocorrelation
Lags = {}; % list of lags
Corrs = {}; % list of correlation values

% for time datas saving
Distances = {}; % all 3D distances
YDistances = {}; % all Dy
XDistances = {}; % all Dx
ZDistances = {}; % all Dz
Times = {}; % all time vectors
Names = {}; % call names


datafile=['D2P_' specif '.mat']; % file to load


% check if file exists
E = exist([df filesep datafile],'file');

if E == 2
    
    load([df filesep datafile])
    fprintf(['\nFile loaded : \n' datafile '. \n']);
    
    % number of curves
    nc=size(MP,2);
    
    
    for ni = 1:nc
        
        AcqName = MP{ni}.name;  % cell name       
        
        CortBot = MP{ni}.CortBot;
        CBotTot = [CBotTot CortBot]; % cortex bottom level (10%)    
        
        CortTop = MP{ni}.CortTop;
        CTopTot = [CTopTot CortTop]; % cortex top level (90%)
        
        CortActi = MP{ni}.CortActi;
        CActTot = [CActTot CortActi]; % cortex fluctuations amplitude
        
        CortActiBeg = MP{ni}.CortActiBeg; 
        CActBegTot = [CActBegTot CortActiBeg]; % cortex fluctuations amplitude for first half of long curves
        
        CortActiEnd = MP{ni}.CortActiEnd;
        CActEndTot = [CActEndTot CortActiEnd]; % cortex fluctuations amplitude for second half of long curves
        
        CortAssym = MP{ni}.CortAssym;
        CAssTot = [CAssTot CortAssym]; % cortex fluctations assymetry
                   
        PeaksDet = PeaksDet + MP{ni}.npeaksD; % number of peaks detected        
        PeaksKept = PeaksKept + MP{ni}.npeaksG; % number of peaks kept
           
        MedCortW = MP{ni}.MedCortW;
        MedCTot = [MedCTot MedCortW]; % median cortex thickness
        
        D3 = MP{ni}.dist; % 3D distances
        DY = MP{ni}.dy;
        DX = abs(MP{ni}.dx);
        DZ = MP{ni}.dz;
        
       T = MP{ni}.time; % time            
        
       % total time (jumps of more than 3 seconds are not considered)
        dT = diff(T);
        ptrdT = dT<3;        
        Ttot = Ttot + sum(dT(ptrdT));
        
        NpTot = NpTot + MP{ni}.npeaksG; % total number of peaks
        
        peakI = MP{ni}.peaks; % peaks data from data2peaks
        
        if ~isempty(peakI)
            
            P  = [P  peakI.prom(peakI.good)']; % peaks prominence
            H  = [H  peakI.height(peakI.good)']; % peaks height
            W  = [W  peakI.width(peakI.good)']; % peak width
            PeakCort = [PeakCort repmat(CortBot,1,sum(peakI.good))]; % cortex thickness par cell for each peak
            
        end
        
        
       % Autocorrelation
       
       MinTime = 300; % minimum time to do autocorrelation
       Fs = 2; % autocorrelation sampling frequency
       
        if (T(end) - T(1) > MinTime+1) && ((T(end)-T(1))/(length(T)-1) <= 1) % done only for curves longer than MinTime and with consistency in time
            
            
            [C,L] = FFTandAUTOCORR(T,D3,Fs); % performs curve autocorrelation. 
            % This function also does FFT but here it is desactivated as it is unsused
            
            ptr = L<MinTime+1; % only taking the first MinTime seconds
            
            C = C(ptr); % Correlation value
            L = L(ptr); % lags            
            
            [~, locmn] = findpeaks(-C,'NPeaks',1,'minpeakprominence',0.01); % first minimum
            
            [~, locmx] = findpeaks(C,'NPeaks',1,'minpeakprominence',0.01); % first maximum
            
            if PLOT % plot autocorr curve and makrs first minimum and maximum
                figure
                hold on
                title([specif ' : Autocorrelation of ' AcqName])
                plot(L,C,'b-*')
                xlabel('Lag (sec)')
                ylabel('AutoCorrelation')
                
                plot(L(locmn),C(locmn),'*r','linewidth',5)
                plot(L(locmx),C(locmx),'*g','linewidth',5)
                
                xlim([0 350])
                ylim([-0.8 1.05])
                
                fig = gcf; % save fig
                saveas(fig,[ff filesep 'Autocorr_' specif ' - ' AcqName '.png'],'png')
                saveas(fig,[ff filesep 'Autocorr_' specif ' - ' AcqName '.fig'],'fig')
                close
                
            end
            
            Lags = [Lags(:)' L']; % lag from
            Corrs = [Corrs(:)' C']; % correlation
            Mins = [Mins L(locmn)]; % First minimum of autocorrelation
            Maxs = [Maxs L(locmx)]; % First maximum of autocorrelation
            
           
        end
        
        Distances{end+1} = D3; % all 3D distances
        YDistances{end+1} = DY; % all Dy
        XDistances{end+1} = DX; % all Dx
        ZDistances{end+1} = DZ; % all Dz
        Times{end+1} = T; % time
        Names{end+1} = AcqName; % cell name
                
        D3Centered = D3' - prctile(D3,100*(1-alignlvl));
        
        AllD3 =  [AllD3 D3']; % all datapoints pulled together
        D3DataPts = [D3DataPts D3Centered]; % all datapoints centered on alignlvl
        
        % activity curve
        ActiCurveX = sort(D3Centered); 
        ActiCurveY = (length(ActiCurveX):-1:1)/length(ActiCurveX)*100;

        ActiCurves = {ActiCurves{:} {ActiCurveX; ActiCurveY}}; % activity curves
        
        if PLOT % ploting activity curves, cortex thickness and fluctuations (no saving)

            
            figure
            hold on
            title(specif)
            plot(ActiCurveX,ActiCurveY,'b')
            ylabel('P(h)')
            xlabel('Height above cortex')
            
            figure
            hold on
            title(specif)
            xlim([-1 2])
            plot(rand(1),MedCortW,'ok','markerfacecolor','b','markersize',10)
            ylabel('Cortex Thickness (nm)')
            
            figure
            hold on
            title(specif)
            xlim([-1 2])
            plot(rand(1),CortActi,'ok','markerfacecolor','b','markersize',10)
            ylabel('Thickness variation amplitude (nm)')
            
        end
        
    end
    
    
% Removing peaks smaller than 15nm (noise)
peakslim = 15;

Psort = P(P>peakslim); 
Hsort = H(P>peakslim);
Wsort = W(P>peakslim);
PeakCortsort = PeakCort(P>peakslim);


%% Storing Data
MS.freq = NpTot/Ttot*60; % global peaks frequency
MS.timetot = Ttot;
MS.proms = Psort;
MS.heights = Hsort;
MS.peaksCortW = PeakCortsort;
MS.AllD3 = AllD3;
MS.D3DataPts = sort(D3DataPts);
MS.widths = Wsort;
MS.CortBots = CBotTot;
MS.CortBot = mean(CBotTot);
MS.CortTops = CTopTot;
MS.CortTop = mean(CTopTot);
MS.CortActisBeg = CActBegTot;
MS.CortActisEnd = CActEndTot;
MS.CortActis = CActTot;
MS.CortAssyms = CAssTot;
MS.MedCortWs = MedCTot;
MS.Acticurves = ActiCurves;
MS.lags = Lags;
MS.corrs = Corrs;
MS.corrmins = Mins;
MS.corrmaxs = Maxs;
MS.mintimecorr = MinTime;
MS.freqsamp = Fs; 


fprintf(['\nTotal number of peaks kepts for plotting : ' ...
    num2str(PeaksKept) '/' num2str(PeaksDet) ' (' num2str(PeaksKept/PeaksDet*100) '%%)\n']);


%% Saving data
save([sf filesep 'P2S' filesep 'TimeDatas_' specif],'Distances','YDistances','XDistances','Times','Names')
save([sf filesep 'P2S' filesep 'P2S_' specif],'MS')
fprintf(['\nFile saved : \nP2S_' specif '.mat\n\n\n'])
    
    
else
    
    cprintf('err',['\nFile : ' datafile ' not found !\n'])
    
end


