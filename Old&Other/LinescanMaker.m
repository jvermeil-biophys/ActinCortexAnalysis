clear all
close all
clc

Cct = [0 109 219]./255; % control

Cl  = [73 0 146]./255; % Latranculin

namesL = {'1-1' '1-237' '1-473' '2-1' '2-257' '2-497' '3-1' '3-336' '3-513' ...
    '4-1' '4-224' '4-401' '5-1' '5-184' '5-508' '6-1' '6-178' '6-479'};

for  kn = 1:length(namesL)
    
    [X,Y] = importfilelinescan(['C:\Users\Valentin\Dropbox\TheseValentin\Papiers\'...
        'GraphePourFigure\Data\Fig2\LatA-' namesL{kn} '_LineScan.txt']);
    
    Y = Y/mean(Y);
    
% detection cell
ind = findchangepts(Y,'MaxNumChanges',4);

pt1 = min(ind);
pt2 = max(ind);

while Y(pt1)>0.1
    pt1 = pt1-1;
end
Xbegin = X(pt1);


while Y(pt2)>0.1
    pt2 = pt2+1;
end
Xend = X(pt2);

% normalization in X
XX = (X - Xbegin)/(Xend-Xbegin);

figure(1)
hold on
plot(XX,Y)

XXX = -0.2:0.01:1.2;
YY = interp1(XX,Y,XXX);

Yall(kn,:) = YY;

clear X XX Y YY pt* Xend Xbegin

end



YmeanL = mean(Yall,1);

plot(XXX,YmeanL,'r','linewidth',2)


namesC = {'1-1' '1-200' '1-300' '2-1' '2-75' '2-200' '3-1' '3-147' '3-460' ...
    '4-1' '4-95' '4-469' '5-1' '5-233' '5-469' '6-1' '6-178' '6-355'};
namesC = {'1-1' '1-200' '1-300' '2-1' '2-75' '2-200' '3-1' '3-147' '3-460' ...
    '4-1' '4-95' '4-339' '5-1' '5-233' '5-469' '6-1' '6-178' '6-355'};

for  kn = 1:length(namesC)
    
    [X,Y] = importfilelinescan(['C:\Users\Valentin\Dropbox\TheseValentin\Papiers\'...
        'GraphePourFigure\Data\Fig2\Ctrl-' namesC{kn} '_LineScan.txt']);

    Y = Y/mean(Y);
% detection cell
ind = findchangepts(Y,'MaxNumChanges',4);

pt1 = min(ind);
pt2 = max(ind);

while Y(pt1)>0.1
    pt1 = pt1-1;
end
Xbegin = X(pt1);


while Y(pt2)>0.1
    pt2 = pt2+1;
end
Xend = X(pt2);

% normalization in X
XX = (X - Xbegin)/(Xend-Xbegin);

figure(2)
hold on
plot(XX,Y)

XXX = -0.2:0.01:1.2;
YY = interp1(XX,Y,XXX);

Yall(kn,:) = YY;

clear X XX Y YY begining ending

end



YmeanC = mean(Yall,1);
figure(2)
plot(XXX,YmeanC,'r','linewidth',2)

figure
hold on
plot(XXX,YmeanC,'linewidth',3,'color',Cct)
plot(XXX,YmeanL,'linewidth',3,'color',Cl)
ylabel('Linescan average')
xlabel('Distance normalized to cell diameter')
legend({'Ctrl n = 6*3','LatA n = 6*3'})