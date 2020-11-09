%% load RelevéFluoInCellOrdonnéCtrl/Nase

figure(1)
hold on
title('Max Intensity')
Spread_Julien(MaxCtrl,1,3000,0.7,'b')
Spread_Julien(MaxNase,2,3000,0.7,'r')
xlim([0 3])
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Ctrl','Neuraminidase'};

yCmax = max([MaxCtrl; MaxNase])*1.1;

pC = ranksum(MaxCtrl,MaxNase);
if 5/100 > pC && pC > 1/100
    txtC = ['*']  ;
elseif 1/100 > pC && pC > 1/1000
    txtC = ['**'];
elseif 1/1000 > pC
    txtC = ['***'];
else
    txtC = ['NS'];
end

line([1 2], [yCmax yCmax],'color','k', 'linewidth', 1.3)
text(3/2,yCmax*1.05,txtC,'HorizontalAlignment','center','fontsize',25)

ylim([0 yCmax*1.3])


figure(2)
hold on
title('Mean Intensity in cell')
Spread_Julien(MeanCtrl,1,1000,0.7,'b')
Spread_Julien(MeanNase,2,1000,0.7,'r')
xlim([0 3])
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Ctrl','Neuraminidase'};


yCmax = max([MeanCtrl; MeanNase])*1.1;

pC = ranksum(MeanCtrl,MeanNase);
if 5/100 > pC && pC > 1/100
    txtC = ['*']  ;
elseif 1/100 > pC && pC > 1/1000
    txtC = ['**'];
elseif 1/1000 > pC
    txtC = ['***'];
else
    txtC = ['NS'];
end

line([1 2], [yCmax yCmax],'color','k', 'linewidth', 1.3)
text(3/2,yCmax*1.05,txtC,'HorizontalAlignment','center','fontsize',25)
ylim([0 yCmax*1.3])

figure(3)
hold on
title('Std Intensity in cell')
Spread_Julien(StdDevCtrl,1,400,0.7,'b')
Spread_Julien(StdDevNase,2,400,0.7,'r')
xlim([0 3])
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Ctrl','Neuraminidase'};

yCmax = max([StdDevCtrl; StdDevNase])*1.1;

pC = ranksum(StdDevCtrl,StdDevNase);
if 5/100 > pC && pC > 1/100
    txtC = ['*']  ;
elseif 1/100 > pC && pC > 1/1000
    txtC = ['**'];
elseif 1/1000 > pC
    txtC = ['***'];
else
    txtC = ['NS'];
end

line([1 2], [yCmax yCmax],'color','k', 'linewidth', 1.3)
text(3/2,yCmax*1.05,txtC,'HorizontalAlignment','center','fontsize',25)
ylim([0 yCmax*1.3])



%% load RelevéFluoOrdonnéCtrl/Nase
figure(4)
hold on
title('Mode Intensity')
Spread_Julien(ModeCtrl,1,100,0.7,'b')
Spread_Julien(ModeNase,2,100,0.7,'r')
xlim([0 3])

ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Ctrl','Neuraminidase'};

yCmax = max([ModeCtrl; ModeNase])*1.1;

pC = ranksum(ModeCtrl,ModeNase);
if 5/100 > pC && pC > 1/100
    txtC = ['*']  ;
elseif 1/100 > pC && pC > 1/1000
    txtC = ['**'];
elseif 1/1000 > pC
    txtC = ['***'];
else
    txtC = ['NS'];
end

line([1 2], [yCmax yCmax],'color','k', 'linewidth', 1.3)
text(3/2,yCmax*1.05,txtC,'HorizontalAlignment','center','fontsize',25)

ylim([0 yCmax*1.3])



figure(5)
hold on
title('Mean Intensity')
Spread_Julien(MeanCtrl,1,500,0.7,'b')
Spread_Julien(MeanNase,2,500,0.7,'r')
xlim([0 3])
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Ctrl','Neuraminidase'};

yCmax = max([MeanCtrl; MeanNase])*1.1;

pC = ranksum(MeanCtrl,MeanNase);
if 5/100 > pC && pC > 1/100
    txtC = ['*']  ;
elseif 1/100 > pC && pC > 1/1000
    txtC = ['**'];
elseif 1/1000 > pC
    txtC = ['***'];
else
    txtC = ['NS'];
end

line([1 2], [yCmax yCmax],'color','k', 'linewidth', 1.3)
text(3/2,yCmax*1.05,txtC,'HorizontalAlignment','center','fontsize',25)
ylim([0 yCmax*1.3])




figure(6)
hold on
title('Std Intensity')
Spread_Julien(StdDevCtrl,1,200,0.7,'b')
Spread_Julien(StdDevNase,2,200,0.7,'r')
xlim([0 3])
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Ctrl','Neuraminidase'};


yCmax = max([StdDevCtrl; StdDevNase])*1.1;

pC = ranksum(StdDevCtrl,StdDevNase);
if 5/100 > pC && pC > 1/100
    txtC = ['*']  ;
elseif 1/100 > pC && pC > 1/1000
    txtC = ['**'];
elseif 1/1000 > pC
    txtC = ['***'];
else
    txtC = ['NS'];
end

line([1 2], [yCmax yCmax],'color','k', 'linewidth', 1.3)
text(3/2,yCmax*1.05,txtC,'HorizontalAlignment','center','fontsize',25)
ylim([0 yCmax*1.3])