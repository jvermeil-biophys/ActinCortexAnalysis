clear all
close all

pixelconv = 1/15.8;
stepnm = 20*1.34/1.52;

RealDist = (4504+1093)/2;

folder = 'D:\Data\Raw\20.06.17';




%% Bille 1
names1 = {'1' '2-1' '2-2' '3' '4' '5' '6' '7' '8' '9' '10-1' '10-2' '11' '12'  '13' '14' '15' '16' '17' '18' '19' '20'};



for  k =1:length(names1)
    
    filename = [folder filesep 'DistM450My1_' names1{k}];
    stackname = [folder filesep 'DistM450My1_' names1{k} '.tif'];
    
    [VarName1, Area, StdDev, XM, YM, Slice] = importfilemultibeads([filename '_Results.txt']);
    
    [X1,Y1,A1,ptr1,X2,Y2,A2,ptr2,S] = splitBeads_M450My1([XM Slice],...
        stackname,filename,YM,Area,4,0,'Bead1');
    
    
    D1{k} = sqrt((X1-X2).^2+(Y1-Y2).^2)*pixelconv*1000 - RealDist; % en nm


    D1{k} = D1{k}';
    

    RelMvtSmall1{k} = sqrt((X2-X2(1)).^2+(Y2-Y2(1)).^2)*pixelconv*1000;
    
end



%% Bille 2


names2 = {'1' '2-1' '2-2' '3' '4' '5' '7' '8' '10-1' '10-2' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'};



for  k =1:length(names2)
    
    filename = [folder filesep 'DistM450My1_' names2{k}];
    stackname = [folder filesep 'DistM450My1_' names2{k} '.tif'];
    
    [VarName1, Area, StdDev, XM, YM, Slice] = importfilemultibeads([filename '_Results.txt']);
    
    [X1,Y1,A1,ptr1,X2,Y2,A2,ptr2,S] = splitBeads_M450My1([XM Slice],...
        stackname,filename,YM,Area,4,0,'Bead2');
    
    
    D2{k} = sqrt((X1-X2).^2+(Y1-Y2).^2)*pixelconv*1000 - RealDist - 1093; % en nm

D2{k} = D2{k}';
    
    RelMvtBig{k} = sqrt((X1-X1(1)).^2+(Y1-Y1(1)).^2)*pixelconv*1000;
    RelMvtSmall2{k} = sqrt((X2-X2(1)).^2+(Y2-Y2(1)).^2)*pixelconv*1000;
    
end

%% Distance entre deux petite billes


names2 = {'1' '2-1' '2-2' '3' '4' '5' '7' '8' '10-1' '10-2' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'};



for  k =1:length(names2)
    
    filename = [folder filesep 'DistM450My1_' names2{k}];
    stackname = [folder filesep 'DistM450My1_' names2{k} '.tif'];
    
    [VarName1, Area, StdDev, XM, YM, Slice] = importfilemultibeads([filename '_Results.txt']);
    
    [X1,Y1,A1,ptr1,X2,Y2,A2,ptr2,S] = splitBeads_M450My1([XM Slice],...
        stackname,filename,YM,Area,4,0,'2Small');
    
    
    D3{k} = sqrt((X1-X2).^2+(Y1-Y2).^2)*pixelconv*1000 - 1093; % en nm

D3{k} = D3{k}';
    
    
end

%% plot



figure(1)
hold on
title('Bille 1 bleu  Bille2 rouge  2petite vert')
xlabel('Slice')
ylabel('MeasuredDistance-RealDistance (nm)')

for k =1:length(names1)
    Mean1(k) = mean(D1{k});
    Std1(k) = std(D1{k});
    plot(D1{k},'-*b')
end


for k =1:length(names2)
    Mean2(k) = mean(D2{k});
    Std2(k) = std(D2{k});
    plot(D2{k},'-*r')
end


for k =1:length(names2)
    Mean3(k) = mean(D3{k});
    Std3(k) = std(D3{k});
    
    plot(D3{k},'-*g')
end


figure
hold on
title('Dist au point de départ individuel des billes, petite bleu(1) vert(2) grosse rouge')
for k =1:length(names2)
    plot(RelMvtBig{k},'-*r')
    plot(RelMvtSmall1{k},'-*b')
    plot(RelMvtSmall2{k},'-*g')
end



figure
hold on
histogram([D1{:}],'binwidth',20,'normalization','pdf')
histogram([D2{:}],'binwidth',20,'normalization','pdf')


figure
title('Mean dist')
hold on
Spread_Julien(Mean1,1,0.8,'b','o')
plot([0.55 1.45],[mean(Mean1) mean(Mean1)],'k--')
Spread_Julien(Mean2,2,0.8,'r','o')
plot([1.55 2.45],[mean(Mean2) mean(Mean2)],'k--')
Spread_Julien(Mean3,3,0.8,'g','o')
plot([2.55 3.45],[mean(Mean3) mean(Mean3)],'k--')
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'Bille1' 'Bille2' '2small'}

figure
title('Std dist')
hold on
Spread_Julien(Std1,1,0.8,'b','o')
Spread_Julien(Std2,2,0.8,'r','o')
Spread_Julien(Std3,3,0.8,'g','o')
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'Bille1' 'Bille2' '2small'}