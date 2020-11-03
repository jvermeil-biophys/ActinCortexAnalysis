clear all

%% Paths & options
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';
RawdataFolder = 'D:\Data\Raw\';
MatfileFolder = 'D:\Data\MatFile';
DropboxDataFolder = 'C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Données';
FigureFolder = 'D:\Data\Figures';

matdf = [MatfileFolder filesep 'FCA'];
mkdir(matdf);
matsf = [MatfileFolder filesep 'FCP'];
mkdir(matsf);


%% get data
load([matdf filesep 'FluoChannelData.mat']);


%% find and process peaks
for k= 1:length(Data)
   clear N profile xprofile
    N = size(Data(k).areas,1);
    
    for n = 1:N
    profile{n} = Data(k).lines{n};
    xprofile{n} = Data(k).xlines{n};
    
    avgint = mean(profile{n});
    
    [pksVal pksLoc pksWidth pksProm] = findpeaks(profile{n},xprofile{n},'MinPeakProminence',0);
    
    pksValNorm = pksVal/avgint;
    
    pksPromNorm = pksProm/avgint;
    
    arealength = Data(k).areas(n,3)/Data(k).Res;
    
    Data(k).peaksdata{n} = [pksVal' pksValNorm' pksLoc' pksWidth' pksProm' pksPromNorm'];
    Data(k).avgint(n) = avgint;
    Data(k).cond = Data(k).name(1);
    Data(k).arealength(n) = arealength;
    
    end
    
end

    save([matsf filesep 'FluoChannelDataProcessed.mat'],'Data')
