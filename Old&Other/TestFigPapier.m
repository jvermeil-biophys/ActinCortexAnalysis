Cdm = [0 109 219]./255; % DMSO
Cl  = [73 0 146]./255; % Latranculin
Co  = [180 180 180]./255; % outside

ThickCt = [];
ForCt = [];

[ThickCt,ForCt] = Var2Data_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,'1s',4415,MatfileFolder,FigureFolder,ThickCt,ForCt);
[ThickCt,ForCt] = Var2Data_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,'1s',4415,MatfileFolder,FigureFolder,ThickCt,ForCt);
[ThickCt,ForCt] = Var2Data_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickCt,ForCt);
[ThickCt,ForCt] = Var2Data_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickCt,ForCt);
[ThickCt,ForCt] = Var2Data_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickCt,ForCt);
[ThickCt,ForCt] = Var2Data_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickCt,ForCt);

ThickLat = [];
ForLat = [];

[ThickLat,ForLat] = Var2Data_Comp('25-09-18','R80',0.56,'M1','DC-Comp_LatA500',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickLat,ForLat);

ThickOut = [];
ForOut = [];

[ThickOut,ForOut] = Var2Data_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Outside',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickOut,ForOut);
[ThickOut,ForOut] = Var2Data_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Outside',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickOut,ForOut);
[ThickOut,ForOut] = Var2Data_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Outside',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickOut,ForOut);
[ThickOut,ForOut] = Var2Data_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Outside',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickOut,ForOut);

ThickBare= [];
ForBare = [];

[ThickBare,ForBare] = Var2Data_Comp('23-01-19','R80',0.6,'M1','DC-Comp_Naked',95,179,'1s',4438,MatfileFolder,FigureFolder,ThickBare,ForBare);


MeanThickCt = mean(ThickCt,2);
StdThickCt = std(ThickCt-ThickCt(1,:),0,2);
MeanForCt = mean(ForCt,2);

MeanThickLat = mean(ThickLat,2);
StdThickLat = std(ThickLat-ThickLat(1,:),0,2);
MeanForLat = mean(ForLat,2);

MeanThickOut = mean(ThickOut,2);
StdThickOut = std(ThickOut-ThickOut(1,:),0,2);
MeanForOut = mean(ForOut,2);

MeanThickBare = mean(ThickBare,2);
StdThickBare = std(ThickBare-ThickBare(1,:),0,2);
MeanForBare = mean(ForBare,2);




% comp only align
figure
hold on
plot(MeanThickCt(1:47)-MeanThickCt(1),MeanForCt(1:47),'linewidth',7,'color',Cdm)
plot(MeanThickCt(1:47)-MeanThickCt(1),MeanForCt(1:47)+StdThickCt(1:47),'--','linewidth',1,'color',Cdm)
plot(MeanThickCt(1:47)-MeanThickCt(1),MeanForCt(1:47)-StdThickCt(1:47),'--','linewidth',1,'color',Cdm)
fill([MeanThickCt(1:47)-MeanThickCt(1); MeanThickCt(47:-1:1)-MeanThickCt(1); MeanThickCt(1)-MeanThickCt(1)],...
    [MeanForCt(1:47)+StdThickCt(1:47); MeanForCt(47:-1:1)-StdThickCt(47:-1:1); MeanForCt(1)+StdThickCt(1)],...
    Cdm,'facealpha',0.5,'edgecolor','none')

plot(MeanThickLat(1:47)-MeanThickLat(1),MeanForLat(1:47),'linewidth',7,'color',Cl)
plot(MeanThickLat(1:47)-MeanThickLat(1),MeanForLat(1:47)+StdThickLat(1:47),'--','linewidth',1,'color',Cl)
plot(MeanThickLat(1:47)-MeanThickLat(1),MeanForLat(1:47)-StdThickLat(1:47),'--','linewidth',1,'color',Cl)
fill([MeanThickLat(1:47)-MeanThickLat(1); MeanThickLat(47:-1:1)-MeanThickLat(1); MeanThickLat(1)-MeanThickLat(1)],...
    [MeanForLat(1:47)+StdThickLat(1:47); MeanForLat(47:-1:1)-StdThickLat(47:-1:1); MeanForLat(1)+StdThickLat(1)],...
    Cl,'facealpha',0.5,'edgecolor','none')

plot(MeanThickOut(1:47)-MeanThickOut(1),MeanForOut(1:47),'linewidth',7,'color',Co)
plot(MeanThickOut(1:47)-MeanThickOut(1),MeanForOut(1:47)+StdThickOut(1:47),'--','linewidth',1,'color',Co)
plot(MeanThickOut(1:47)-MeanThickOut(1),MeanForOut(1:47)-StdThickOut(1:47),'--','linewidth',1,'color',Co)
fill([MeanThickOut(1:47)-MeanThickOut(1); MeanThickOut(47:-1:1)-MeanThickOut(1); MeanThickOut(1)-MeanThickOut(1)],...
    [MeanForOut(1:47)+StdThickOut(1:47); MeanForOut(47:-1:1)-StdThickOut(47:-1:1); MeanForOut(1)+StdThickOut(1)],...
    Co,'facealpha',0.5,'edgecolor','none')

plot(MeanThickBare(1:47)-MeanThickBare(1),MeanForBare(1:47),'linewidth',7,'color','k')
plot(MeanThickBare(1:47)-MeanThickBare(1),MeanForBare(1:47)+StdForBare(1:47),'--','linewidth',1,'color','k')
plot(MeanThickBare(1:47)-MeanThickBare(1),MeanForBare(1:47)-StdForBare(1:47),'--','linewidth',1,'color','k')
fill([MeanThickBare(1:47)-MeanThickBare(1); MeanThickBare(47:-1:1)-MeanThickBare(1); MeanThickBare(1)-MeanThickBare(1)],...
    [MeanForBare(1:47)+StdForBare(1:47); MeanForBare(47:-1:1)-StdForBare(47:-1:1); MeanForBare(1)+StdForBare(1)],...
    'k','facealpha',0.5,'edgecolor','none')

xlabel('Indentation (nm)')
ylabel('Force (pN)')


% comp and rel
figure
hold on
plot(MeanThickCt,MeanForCt,'linewidth',3,'color',Cdm)

plot(MeanThickLat,MeanForLat,'linewidth',3,'color',Cl)

plot(MeanThickOut,MeanForOut,'linewidth',3,'color',Co)

plot(MeanThickBare,MeanForBare,'linewidth',3,'color','k')

xlabel('Thickness (nm)')
ylabel('Force (pN)')



% comp and rel align
figure
hold on
plot(MeanThickCt-MeanThickCt(1),MeanForCt,'linewidth',3,'color',Cdm)

plot(MeanThickLat-MeanThickLat(1),MeanForLat,'linewidth',3,'color',Cl)

plot(MeanThickOut-MeanThickOut(1),MeanForOut,'linewidth',3,'color',Co)

plot(MeanThickBare-MeanThickBare(1),MeanForBare,'linewidth',3,'color','k')

xlabel('Indentation (nm)')
ylabel('Force (pN)')
