clear all
close all
clc

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)

% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';

pixelconv=1/15.8;

% cd('G:\Raw\17.11.03\TailleBille3') % not found ?
% cd('G:\Raw\BeadsSizes\BeadSize_18.10.25') % 36 results
% cd('G:\Raw\BeadsSizes\BeadSize_17.10.12') % 14-39 results et mauvais std
% cd('G:\Raw\BeadsSizes\BeadSize_17.10.13') % 40 results
% 
% cd('G:\Raw\BeadsSizes\BeadSize_19.05.07\BeadSizeNoCoat\5mT') % 1:19
% titlestr = 'M-450 beads size new batch nocoat 5mT in water';
% 
% cd('G:\Raw\BeadsSizes\BeadSize_19.05.07\BeadSizeNoCoat\50mT') % 1:19
% titlestr = 'M-450 beads size new batch nocoat 50mT in water';
% 
% cd('G:\Raw\BeadsSizes\BeadSize_19.05.07\BeadSizeSerum\5mT') % 1:32
% titlestr = 'M-450 beads size new batch serum 5mT in water';

% cd('G:\Raw\BeadsSizes\BeadSize_19.05.07\BeadSizeSerum\50mT') % 1:26
% titlestr = 'M-450 beads size new batch serum 50mT in water';

% cd('G:\Raw\BeadsSizes\BeadSize_fromTrueOutside020419') % 1:5
% titlestr = 'M-450 beads from trueoutside exp in water';

% cd('G:\Raw\BeadsSizes\BeadSize_19.05.09\OldBatchNocoating5mT') % 1:22
% titlestr = 'M-450 beads from old batch nocoat 5mT';

cd('G:\Raw\BeadsSizes\BeadSize_19.05.20\BeadSize\BeadSizeNewBeads_Medium_5mT') % 1:42
titlestr = 'M-450 beads from new batch new serum 5mT in med';

% cd('G:\Raw\BeadsSizes\BeadSize_19.05.20\BeadSize\BeadSizeNewBeads_Medium_50mT') % 1:38
% titlestr = 'M-450 beads from new batch new serum 50mT in med';
% 
% cd('G:\Raw\BeadsSizes\BeadSize_19.05.20\BeadSize\BeadSizeNewBeads_Water_5mT') % 1:46
% titlestr = 'M-450 beads from new batch new serum 5mT in water';
% % 
% cd('G:\Raw\BeadsSizes\BeadSize_19.05.20\BeadSize\BeadSizeNewBeads_Water_50mT') % 1:23
% titlestr = 'M-450 beads from new batch new serum 50mT in water';



AllD = [];
AllStd = [];
AllAlpha = [];

for i = 1:42
    resname = ['Results' num2str(i) '.txt'];
    Mat=dlmread(resname,'\t',1,0);
     fid = fopen(resname);
    HeadAll = fgetl(fid);
    fclose(fid);
    HeadCell = textscan(HeadAll,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
    HeadTable = [HeadCell{:}];
    
    Nxm = find(strcmp(HeadTable,'XM')) +1 ;
    Nym = find(strcmp(HeadTable,'YM')) +1 ;
    
    Nstd = find(strcmp(HeadTable,'StdDev')) +1 ;
    Nmax = find(strcmp(HeadTable,'Max')) +1 ;
    Nmin = find(strcmp(HeadTable,'Min')) +1 ;
    
    Std = Mat(:,Nstd);
    
    Max = Mat(:,Nmax);
    Min = Mat(:,Nmin);
    
    m = (Std > 8000)&(Max<65534)&(Min<32768);
   
    
    
    X = Mat(m,Nxm)*pixelconv;
    Y = Mat(m,Nym)*pixelconv;
    
    A = sortrows([X,Y]);
    
    X = A(:,1);
    Y = A(:,2);
    

       
    [Dx,Dy] = deal([]);
    for kk = 1:length(X)
    Dx = [Dx; X-circshift(X,kk)];
    Dy = [Dy; Y-circshift(Y,kk)];
    end
    
    alpha = rad2deg(atan(Dy./Dx));
    
    D = sqrt(Dx.^2+Dy.^2);
    
    mD = D < 4.6 & D> 4.3;
    mAl = alpha > -20 & alpha < 20;
    
    D = unique(D(mD&mAl));
    
    
    AllAlpha = ([AllAlpha unique(alpha(mD))']);
    
    AllStd = [AllStd Std'];
    
    AllD = [AllD D'];
    
    clear m* Std alpha D
end

figure
histogram(AllStd)
legend(['n = ' num2str(length(AllStd)) ' Mean = ' num2str(mean(AllStd)) ' Std = ' num2str(std(AllStd))])
ylabel('Count')
xlabel('Bead Size (µm)')
title('All intensity std','interpreter','none')

figure
histogram(AllAlpha)
legend(['n = ' num2str(length(AllAlpha)) ' Mean = ' num2str(mean(AllAlpha)) ' Std = ' num2str(std(AllAlpha))])
ylabel('Count')
xlabel('Bead Size (µm)')
title('All angles','interpreter','none')

figure
h = histogram(AllD,30,'normalization','probability')
xlim([4.3 4.6])
title(['n = ' num2str(length(AllD)) ' Mean = ' num2str(mean(AllD)) ' Std = ' num2str(std(AllD))])
ylabel('Frequency (%)')
xlabel('Bead Size (µm)')

ylim([0 0.18])

ax = gca;
ax.LineWidth = 1.5;
ax.YTickLabels = [0 2 4 6 8 10 12 14 16 18]
box on



cd(h)
