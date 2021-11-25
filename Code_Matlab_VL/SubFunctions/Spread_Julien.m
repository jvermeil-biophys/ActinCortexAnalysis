function Spread_Julien(YData,XCenter,spread,color,symbol,ybin)
%% Input: donnée à plotter, abscisse, taille des bin, plotcolor

if nargin == 4
    symbol = 'o';
end

if nargin < 6
    
    
    
try
ybin =median(diff(sort(YData)))*12;
Ysort=sort(YData);
[N,~,bin]=histcounts(Ysort,min(Ysort)-1:ybin:max(Ysort));   % ici
catch
ybin =median(diff(sort(YData)))*2;
Ysort=sort(YData);
[N,~,bin]=histcounts(Ysort,min(Ysort)-1:ybin:max(Ysort));
end
else
   Ysort=sort(YData);
   
   binvect = [0 ybin];
   
   i = 3;
   while min(Ysort)-1+sum(binvect) < max(Ysort)
   
   binvect(i) = ybin*(1.2)^(i-1);
   i = i+1;
   end
   
   
[N,~,bin]=histcounts(Ysort,min(Ysort)-1+cumsum(binvect));
end


nx=zeros(size(Ysort));
shift=zeros(size(Ysort));
for i=1:length(N)
    ptr=find(bin==i);
    if ~isempty(ptr)
        if iseven(length(ptr))
            nx(bin==i)= (1:length(ptr));
            shift(bin==i) = -0.5; 
        else
            nx(bin==i)= 0:length(ptr)-1;
        end
    end
end



x = XCenter + ((-1).^nx.*(round(nx/2)+shift)+(rand(1,length(nx))*0-0))/max([max(N),8])*spread ;

plot(x,Ysort,symbol,'linewidth',0.2,'markeredgecolor','k','markerfacecolor',color,'color',color,'markersize',6,'handlevisibility','off')

