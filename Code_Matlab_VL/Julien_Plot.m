%figure beebox

% you need to import the ExportSummaryTable you sent me

E=[];H=[];X=[];
k=[4 2];  %that's the vector to choose whch conditions you want to print
%and in which order. The increment of the loop will go through this vector
%k=[4 3 2 1];

for j=1:length(k)
    i=k(j);

    A=ExportSummaryStructure(i).meanEfitChadwick;
    X=[X; j*ones(size(A))];
    E=[E; A];
    B=ExportSummaryStructure(i).meanH0fitChadwick;
    H=[H; B];

end
% a loop to create append the values of E, H and X in two vectors. X will
% have the show a digit with the condition number.


figure (1);
hold off;
subplot(1,2,2)
a1=gca; % this is get current axes, useful to manipulate a graph.
beeswarm(X,log10(E));
hold on;
h=boxplot(log10(E),X,'Colors','k','Width',0.5,'symbol','');
% for E both beeswarm and boxplot are drawn with log10(E) to have a regular
% spacing of the beeswarm.

set(h,'LineWidth',0.7)
a1.YTick=3:5;ylim([2.7 4.8]);
yyt=[2+log10(1:9) 3+log10(1:9) 4+log10(1:9) 5];
a1.YAxis.MinorTick='on';
a1.YAxis.MinorTickValues=yyt;
% this trick add manualy some log ticks to make it look like a regular log
% axis you can add more points to yyt if your data go below 100 or above 10^5
a1.YAxis.TickLength=[.03 .03];
a1.YTickLabel={'1' '10' '100'};
ylabel('Cortex elastic modulus (kPa)');
a1.XTick=1:2; a1.XTickLabel={'Floating' 'Adhered'};
%a1.XTick=1:4; a1.XTickLabel={'Non-adhered' 'Non-adhered \n + Linker 'Adhered' 'Adhered \n +Linker'};
set(a1,'FontSize',12)
box on

subplot(1,2,1);

beeswarm(X(H<1000),H(H<1000));
% to remove "outoutliers"

hold on;
h=boxplot(H(H<1000),X(H<1000),'Colors','k','Width',0.5,'symbol','');
set(h,'LineWidth',0.7)
ylim([0 1000]);
a2=gca;
a2.YTick=0:200:1000;
ylabel('Cortex thickness (nm)');
a2.XTick=1:2; a2.XTickLabel={'Floating' 'Adhered'};
a2.YAxis.TickLength=[.02 .02];
%a2.XTick=1:4; a2.XTickLabel={'Non-adhered' 'Non-adhered \n + Linker 'Adhered' 'Adhered \n +Linker'};
set(a2,'FontSize',12)
box on
