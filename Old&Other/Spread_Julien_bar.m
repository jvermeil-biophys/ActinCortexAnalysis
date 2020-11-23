function Spread_Julien_bar(YData,s,color,ybin)
%% Input: donnée à plotter, abscisse, ,couleur,'Marker'

Ysort=sort(YData);
[N,edges,bin]=histcounts(Ysort,min(Ysort):ybin:max(Ysort)); 



hold on
for i=1:length(N)
    plot([0 N(i)]/sum(N)*4/2+s,[edges(i)+ybin/2 edges(i)+ybin/2],'LineWidth',3.5,'color',color)
end
