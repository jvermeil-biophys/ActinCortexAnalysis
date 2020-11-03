clear all
close all

%% INIT

warning('off','all')

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextInterpreter','Latex')

here = pwd;

%% Paths & options
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';
RawdataFolder = 'D:\Data\Raw\';
MatfileFolder = 'D:\Data\MatFile';
DropboxDataFolder = 'C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Données';
FigureFolder = 'D:\Data\Figures';

matdf = [MatfileFolder filesep 'FCP'];
mkdir(matdf);
matsf = [MatfileFolder filesep 'FCFig'];
mkdir(matsf);

% Couleurs par drogues
Cct = [0 109 219]./255; % control

Cl  = [73 0 146]./255; % Latranculin

Cb  = [219 209 0]./255; % blebbi
Cy  = [255 255 109 ]./255; % y27
Ccl = [146 0 0]./255; % calyculinA 10nM
Ccl10 = [146 0 0]./255./2; % calyculinA 1nM

Cs  = [36 255 36]./255; % smifh2

Cck = [182 109 255]./255;

Cin = [0 73 73]./255; % inside
Co  = [255 182 119]./255; % outside

%% Get datas
load([matdf filesep 'FluoChannelDataProcessed.mat'],'Data')

conds = [Data(:).cond];

ptrctrl = conds == 'C';
ptrblebbi = conds == 'B';

CtrlPdata = [Data(ptrctrl).peaksdata];
BlebbiPdata = [Data(ptrblebbi).peaksdata];

CtrlInt = [];
CtrlIntNorm = [];
CtrlPromInt = [];
CtrlPromIntNorm = [];
CtrlSize = [];

for k = 1:length(CtrlPdata)
    CtrlInt = [CtrlInt; CtrlPdata{k}(:,1)];
    CtrlIntNorm = [CtrlIntNorm; CtrlPdata{k}(:,2)];
    CtrlPromInt = [CtrlPromInt; CtrlPdata{k}(:,5)];
    CtrlPromIntNorm = [CtrlPromIntNorm; CtrlPdata{k}(:,6)];
    CtrlSize = [CtrlSize; CtrlPdata{k}(:,4)];
    
end

BlebbiInt = [];
BlebbiIntNorm = [];
BlebbiPromInt = [];
BlebbiPromIntNorm = [];
BlebbiSize = [];

for k = 1:length(BlebbiPdata)
    BlebbiInt = [BlebbiInt; BlebbiPdata{k}(:,1)];
    BlebbiIntNorm = [BlebbiIntNorm; BlebbiPdata{k}(:,2)];
    BlebbiPromInt = [BlebbiPromInt; BlebbiPdata{k}(:,5)];
    BlebbiPromIntNorm = [BlebbiPromIntNorm; BlebbiPdata{k}(:,6)];
    BlebbiSize = [BlebbiSize; BlebbiPdata{k}(:,4)];
    
end

CtrlLength = sum([Data(ptrctrl).arealength]);


BlebbiLength = sum([Data(ptrblebbi).arealength]);

CtrlClength  = [Data(ptrctrl).Clength];
BlebbiClength  = [Data(ptrblebbi).Clength];


proms = [0.01:0.01:0.8];

for kp = 1:length(proms)
    
Cptrprom = CtrlPromIntNorm>proms(kp);
CtrlNumb(kp) = sum(Cptrprom);

Bptrprom = BlebbiPromIntNorm>proms(kp);
BlebbiNumb(kp) = sum(Bptrprom);

end
%% Plots
% distribution of peaks data
figure(1)
hold on
histogram(CtrlInt,'binwidth',150,'facecolor',Cct,'facealpha',0.3)
histogram(BlebbiInt,'binwidth',150,'facecolor',Cb,'facealpha',0.3)
ylabel('Count')
xlabel('Peaks Intensity (a.u.)')


figure(2)
hold on
histogram(CtrlIntNorm,'binwidth',0.1,'facecolor',Cct,'facealpha',0.3)
histogram(BlebbiIntNorm,'binwidth',0.1,'facecolor',Cb,'facealpha',0.3)
ylabel('Count')
xlabel('Peaks Normalized Intensity')


figure(3)
hold on
histogram(CtrlPromInt,'binwidth',30,'facecolor',Cct,'facealpha',0.3)
histogram(BlebbiPromInt,'binwidth',30,'facecolor',Cb,'facealpha',0.3)
ylabel('Count')
xlabel('Peaks Prominence Intensity (a.u.)')


figure(4)
hold on
histogram(CtrlPromIntNorm,'binwidth',0.01,'facecolor',Cct,'facealpha',0.3)
histogram(BlebbiPromIntNorm,'binwidth',0.01,'facecolor',Cb,'facealpha',0.3)
ylabel('Count')
xlabel('Peaks Normalized Prominence Intensity (a.u.)')


figure(5)
hold on
histogram(CtrlSize,'binwidth',0.1,'facecolor',Cct,'facealpha',0.3)
histogram(BlebbiSize,'binwidth',0.1,'facecolor',Cb,'facealpha',0.3)
ylabel('Count')
xlabel('half prom width (µm)')

% cell sizes 
figure(6)
hold on
Spread_Julien(CtrlClength,1,0.7,Cct,'o',5)
Spread_Julien(BlebbiClength,2,0.7,Cb,'o',5)

% peaks frequencies
figure(7)
hold on
errorbar(proms*100,BlebbiNumb/BlebbiLength,sqrt(BlebbiNumb)/BlebbiLength,'k.')
plot(proms*100,BlebbiNumb/BlebbiLength,'ko','markerfacecolor',Cb,'markersize',8)
errorbar(proms*100,CtrlNumb/CtrlLength,sqrt(CtrlNumb)/CtrlLength,'k.')
plot(proms*100,CtrlNumb/CtrlLength,'ko','markerfacecolor',Cct,'markersize',8)

ylabel('Peaks per micron of cell cortex')
xlabel('Threshold for peak selection (\% of average intensity)')
yl = ylim;
ylim([0 yl(2)])

