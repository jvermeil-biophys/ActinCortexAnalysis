ZStep =[400 500 800 1000 1500 2000];

date = datestr(now,'yy-mm-dd');

load('D:\Data\Raw\EtalonageZ\New\Kimo.mat')
mkdir(['D:\Data\Figures\' date '\KimoCrossCorr'])

for k = 1:length(ZStep)

os = ZStep(k)/stepnm;

Lk = size(K,1);

x = 1:(Lk-2*os);
x = x + os - f;

x = x*stepnm;

clear dk*

for i = os+1:Lk-os
    
    dk12(i-os) = sqrt(sum(K(i-os,:)-K(i,:)).^2);
    dk23(i-os) = sqrt(sum(K(i,:)-K(i+os,:)).^2);
    dk13(i-os) = sqrt(sum(K(i-os,:)-K(i+os,:)).^2);
    
    dk = dk12 + dk23 + dk13;
    
end


figure
plot(x,dk12)
hold on
plot(x,dk23)
plot(x,dk13)
plot(x,dk)
plot([((110-f)*stepnm) ((250-f)*stepnm)],[0 0],'--c','linewidth',2)
plot(0,0,'*r')
plot([((250-f)*stepnm) ((250-f)*stepnm)],[0 1.3*max(dk)],'--c','linewidth',1.2)
plot([((110-f)*stepnm) ((110-f)*stepnm)],[0 1.3*max(dk)],'--c','linewidth',1.2)

title(['CrossCorrelation for ZStep = ' num2str(ZStep(k))])
legend('B+F','F+H','B+H','Sum','TrackRange','Focus','location','northwest')
hold off

fig = gcf;

saveas(fig,['D:\Data\Figures\' date '\KimoCrossCorr\CrossCorrZS' num2str(ZStep(k))],'png')


close(fig)
end
