clear all
close all

pixelconv=1/15.8;

[VarName, Area, StdDev, XM, YM, Slice] = importfile100620(['D:\Data\Raw\20.06.09\TailleResults.txt']);

nS = max(Slice);

D1 = [];

for k = 1:nS
    ptrslice = Slice == k;
    
A = sortrows([XM(ptrslice),YM(ptrslice)]);
    
    X = A(:,1)*pixelconv;
    Y = A(:,2)*pixelconv;
    

       
    [Dx,Dy] = deal([]);
    for kk = 1:length(X)
    Dx = [Dx; X-circshift(X,kk)];
    Dy = [Dy; Y-circshift(Y,kk)];
    end
    
    
    D1 = [D1 sqrt(Dx.^2+Dy.^2)'];
        
    
end

 mD = D1 < 1.3 & D1> 0.9;

D1 = unique(D1(mD));

figure
histogram(D1)
hold on
title(['Mean = ' num2str(mean(D1)) ' - Std = ' num2str(std(D1))])



[VarName, Area, StdDev, XM, YM, Slice] = importfile100620(['D:\Data\Raw\20.06.11\TailleBille11-06-20.txt']);

nS = max(Slice);

D2 = [];

for k = 1:nS
    ptrslice = Slice == k;
    
A = sortrows([XM(ptrslice),YM(ptrslice)]);
    
    X = A(:,1)*pixelconv;
    Y = A(:,2)*pixelconv;
    

       
    [Dx,Dy] = deal([]);
    for kk = 1:length(X)
    Dx = [Dx; X-circshift(X,kk)];
    Dy = [Dy; Y-circshift(Y,kk)];
    end
    
    
    D2 = [D2 sqrt(Dx.^2+Dy.^2)'];
        
    
end

 mD = D2 < 1.3 & D2> 0.9;

D2 = unique(D2(mD));

figure
histogram(D2)
hold on
title(['Mean = ' num2str(mean(D2)) ' - Std = ' num2str(std(D2))])



figure
histogram([D1 D2])
hold on
title(['Mean = ' num2str(mean([D1 D2])) ' - Std = ' num2str(std([D1 D2]))])
