function Data2Peaks(specif,PLOT,...
    resfolder,figurefolder)

%
% Data2Peaks is used for analyzing constant field curves. It detects peaks
% and validated them based on several criterions. it also computes the
% median thickness of the cortex and the fluctuations amplitude. For curves
% that are longer than 10 minutes it also computes fluctuations for the
% first and second half of the curves. If the cell has two pinching they
% are computed and saved together for later comparison.
%
% Data2Peaks(specif,PLOT,...
%    datafolder,figurefolder)
%
% specif : experimental condiion (ex: 'Dicty_WT')
% PLOT : 1 = plot some graphs, 0 = thou shall not plot
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
df = [resfolder filesep 'V2D'];

% listing file to analyze
dfContents = dir(df);
filelist = {dfContents.name};
ListOfInterest = filelist(contains(filelist,specif));

% save folder
datenow = datestr(now,'yy-mm-dd');
ff = [figurefolder filesep datenow]; % figure folder
mkdir([resfolder filesep 'D2P']) % matfile save folder

% initialisation of variables for double pinched cells
CommonTimes = {};
CurveOnes = {};
CurveTwos = {};
NameOnes = {};
NameTwos = {};
Lags = {};
Corrs = {};

%% Data retrieval and plot
nf = length(ListOfInterest); % number of files to analyze

ic = 1;

for kf=1:nf
    
    %% Data retrieval and peaks detection
    
    % loading file
    datafile=ListOfInterest{kf};
    
    load([df filesep datafile],'MR');
    fprintf(['File loaded : \n' datafile '. \n\n']);
    
    % number of curves
    nc=size(MR,2);
    
    % fidn cell with double pinching
    codelist = {};
    
    for k = 1:nc
        if not(isempty(MR{k}))
            ptr_ = strfind(MR{k}.name,'_');
            codelist = [codelist(:)' {MR{k}.name(ptr_(1)+1:ptr_(4)-1)}]; % get the cell ID Mx_Py_Cz
        else
            codelist = [codelist(:)' {''}];
        end
    end
    
    if contains(specif,'Outside') || contains(specif,'Inside') % no double interface for inside or outside experiments
        
        ptrDI = [];
        
    else
        
        ptrDI = find(contains(codelist,'-')); % find cell IDs with double pinching (identified by Cz having the form Cz-1 or Cz-2)
        
    end
    
    
    
    for kc=1:nc
        if not(isempty(MR{kc}))
            
            % Name of the experiment and image
            AcqName = MR{kc}.name;
            
            ptr_ = strfind(MR{kc}.name,'_');
            code = MR{kc}.name(ptr_(1)+1:ptr_(4)-1); % get the cell ID Mx_Py_Cz
            
            % get data
            T = MR{kc}.time;
            S = MR{kc}.S;
            D3 = double(MR{kc}.D3);
            Dx = double(MR{kc}.dx);
            Dy = double(MR{kc}.dy);
            Dz = double(MR{kc}.dz);
            
            %% if cell has double pinching we look for the second
            % pinching of the same cell to do cross corelation
            
            if ismember(kc,ptrDI)
                
                ptrDI = ptrDI(ptrDI~=kc); % remove current cell from list
                
                pos = strfind(code,'-')-1; % find cell code (from curve code)
                
                kassoc = find(contains(codelist,code(1:pos))); % find the two curves
                
                if length(kassoc) > 1 % sometimes the second curve was not validated in V2D
                                      % here we check that there is indedd two curves from the same cell
                
                kassoc = kassoc(kassoc~=kc); % get index of second curve
                
                ptrDI = ptrDI(ptrDI ~= kassoc); % remove second curve from list
                
                % load data from second curve and do autocorrelation
                
                D3other = MR{kassoc}.D3;
                Tother = MR{kassoc}.time;
                AcqNameOther = MR{kassoc}.name;
                
                CommonTime = max([T(1),Tother(1)]):0.5:min([T(end),Tother(end)]); % % time vector common to the two curves
                
                % interpolation on common time
                D3int = interp1(T,D3,CommonTime);
                D3otherint = interp1(Tother,D3other,CommonTime);
                
                [c,lags] = xcorr(D3int-mean(D3int),D3otherint-mean(D3otherint)/mean(D3otherint),'normalized');
                
                % saving data
                CommonTimes{end+1} = CommonTime;
                CurveOnes{end+1} = D3int;
                CurveTwos{end+1} = D3otherint;
                NameOnes{end+1} = AcqName;
                NameTwos{end+1} = AcqNameOther;
                Lags{end+1} = lags/2;
                Corrs{end+1} = c;
                
                end
            end
            
            %% peaks detection
            
            % smoothing data for peak detection
            D3s = smooth(T,D3,7,'sgolay',3);
            
            
            [hPks,TpsPks,wPks,pPks] =  findpeaks(D3s,T,'MinPeakProminence',15); % peaks detection
            
            
            bPks = hPks - 0.80*pPks; % base level of peaks
            
            
            IndPks = find(ismember(T,TpsPks) == 1); % Position of peaks maxima
            
            
            NpD = length(IndPks); % number of detected peaks
            
            % Structure with peaks info
            peakI = struct('prom',pPks,'width',wPks,'height',hPks,'loc',IndPks,'promnorm',pPks./(hPks-pPks));
            
            %% peaks sorting/valiation
            
            
            peakI.good = zeros(1,NpD); % Good peaks list
            peakI.info = cell(1,NpD); % Reason for bad peaks
            
            % Getting indices of all points in the peaks
            
            peakI.locs = {};
            
            for ip = 1:NpD
                
                i1 = IndPks(ip); % first index of peak
                i2 = i1; % last index of peaks
                
                while (i1>1) && (D3s(i1) > bPks(ip))
                    i1 = i1 - 1;
                end
                
                while (i2<length(D3s)) && (D3s(i2) > bPks(ip))
                    i2 = i2 + 1;
                end
                
                peakI.locs{ip} = i1:i2; % peak indices
                peakI.Dy(ip)= mean(Dy(peakI.locs{ip})); % Dy during peak
                peakI.Dz(ip)= mean(Dz(peakI.locs{ip})); % Dz during peak
                
                
            end
            
            
            
            % exclusion conditions for peaks : 60% of the peak
            % prominence is between two consecutive points (gapH) or
            % there is more thzn 5 second between two consecutive
            % points (gapT)
            
            
            for ip = 1:NpD
                
                if (max(abs(diff(D3s(peakI.locs{ip}))))>0.6*pPks(ip))
                    peakI.info{ip} = 'gapH';
                    
                elseif (max(diff(T(peakI.locs{ip}))) > 5)
                    peakI.info{ip} = 'gapT';
                    
                else
                    peakI.good(ip) = 1;
                    
                end
            end
            
            if not(isempty(peakI.good))
                NpG = sum(peakI.good); % number of good peaks
            else
                NpG = 0;
            end
            
            peakI.good = logical(peakI.good);
            
            
            %% Activity related computations
            
            CortBot = prctile(D3,10);
            CortTop = prctile(D3,90);
            
            CortActi = CortTop-CortBot; % fluctuations amplitude
            
            % if curve longer than 10 minute comparison betwen first half
            % and second half fluctuations
            if T(end) - T(1) > 600
                
                Tmid = T(1) + (T(end)-T(1))/2;
                
                CortActiBeg = prctile(D3(T<Tmid),90)-prctile(D3(T<Tmid),10);
                CortActiEnd = prctile(D3(T>Tmid),90)-prctile(D3(T>Tmid),10);
                
                ActiBegEnd = 1;
                
            else
                
                ActiBegEnd = 0;
                
            end
            
            % median cortex thickness
            MedCortW = median(D3);
            
            % cortex fluctuations assymetry
            CortAssym = ((CortTop-MedCortW)-(MedCortW-CortBot))/CortActi;
            
            
            if PLOT % curve plotting
                
                % plot the curves with both image n° and time x coordinate
                
                mkdir([ff filesep 'CurvePlots' filesep specif])
                
                figure
                hold on
                plot(S,D3,'-ko','linewidth',2,'markerfacecolor','g','handlevisibility','off')
                ylabel('D3 (nm)')
                xlabel('Image n°')
                namename = 'Images';
                
                fig = gcf;
                saveas(fig,[ff filesep 'CurvePlots' filesep specif filesep AcqName '_' namename],'png')
                close
                
                figure
                hold on
                plot(T,D3,'-ko','linewidth',2,'markerfacecolor','b','handlevisibility','off')
                ylabel('D3 (nm)')
                xlabel('Time (s)')
                namename = 'Time';
                
                fig = gcf;
                saveas(fig,[ff filesep 'CurvePlots' filesep specif filesep AcqName '_' namename],'png')
                close
                
            end
            
            
            if NpD == 0 % if no peaks detected
                NpG = 0; % no good peaks kept
            end
            
            
            
            
            %% Store data
            
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
            MP{ic}.time = T;
            MP{ic}.dx = Dx;
            MP{ic}.dz = Dz;
            MP{ic}.dy = Dy;
            MP{ic}.dist = D3;
            MP{ic}.peaks = peakI;
            MP{ic}.name = AcqName;
            
            
            ic = ic+1;
        end
    end
end

%% Saving data
Sname = ['D2P_' specif];

if exist('MP','var')
    
    
    savename = [resfolder filesep 'D2P' filesep Sname];
    save(savename,'MP')
    
    fprintf(['File saved : \n' Sname '.mat \n\n']);
end

% separate saving for double pinching
if ~isempty(CommonTimes)
    save([resfolder filesep 'D2P' filesep 'DoubleInt_' specif],'CommonTimes','CurveOnes', 'CurveTwos',...
        'NameOnes', 'NameTwos', 'Lags', 'Corrs')
end

end

