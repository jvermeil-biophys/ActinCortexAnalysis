%clear
close all

%load resuniq.mat


N=length(All);
%field=zeros(N,1);
%gradB=field;

B=All(:,5);  %B en miliTesla 
% a modifier pour importer le champ en mT

gradB=All(:,6); %gradB en tesla /m

mu0=pi*4e-7;

R=1.3455; % rayon des billes en micron à changer selon la taille

%H0=field*1e-3/mu0;

x=All(:,2); %x en µm
% cette conversion est celle pixel/ µm pour une orca flash et un objectif
% 100x

t=All(:,4);  %t en secondes
% a changer pour importer le temps en secondes

U=x./t*1e-6;   %U vitesse en m/s

V=4/3*pi*R^3*1e-18;   %V volume en m^3

Fmes=6*pi*R.*1e-3.*U*1e6;  %Force mesurée à partir de la viscosité en pN

figure(1); hold on;
plot(B,Fmes,'o')
xlabel('B (mT)'); ylabel('Force (pN)')

M=Fmes./(V*gradB*1e12);

figure(2); hold on;
plot(B,M,'o')
xlabel('B (mT)'); ylabel('M (A.m-1)')

uB = unique(B);

for jk = 1:length(uB)
    
    Mf(jk) = mean(M(B == uB(jk)));
    
end

B0=0:0.1:100;
M0 = (0.001991*B0.^3+17.54*B0.^2+153.4*B0)./(B0.^2+35.53*B0+158.1)*1600*1.05;  
%Valeurs tabulées de la magnétisation pour les M450

AdjustM = fittype('A*(0.001991*x.^3+17.54*x.^2+153.4*x)./(x.^2+35.53*x+158.1)*1600*1.05');

fittedM = fit(B,M,AdjustM);
Coeff = fittedM.A;

CI = confint(fittedM);

PolyM = fittype('A*(B*x.^3+C*x.^2+D*x+H)./(E*x.^2+F*x+G)');
fittedPolyM = fit([0; B],[0; M],PolyM,'lower',[0 0 0 0 0 0 0 0]); % ,'StartPoint',[1600 0.001991 17.54 153.4 1 35.53 158.1 0]

TanHfit = fittype('A*tanh(B*x)');
fittedTanH = fit([0; B],[0; M],TanHfit,'StartPoint',[10000 0.01]);



plot(B0,M0,'-*')

plot(B0,fittedM(B0),'-g','linewidth',2)
plot(B0,fittedPolyM(B0),'--r','linewidth',1.5)


plot(B0,fittedTanH(B0),'-m','linewidth',2)

plot(B0,CI(1)*M0,'--g','handlevisibility','off')
plot(B0,CI(2)*M0,'--g','handlevisibility','off')

% % mean
% MfInt = interp1([0 uB(1:3)'],[0 Mf(1:3)],1:0.1:20,'spline');
% plot(1:0.1:20,MfInt,'-m','linewidth',2)



legend('Data M270','M450 Tabulées',['Fit pour M270 - Mag = ' num2str(Coeff) '*M450'],'Fit poly3/poly2')
legend('location','best')



%% calcul d'erreur 
% Delta=((M(2:3:end)+M(3:3:end))/2-interp1(Btab,Muniq,B(1:3:end),'spline'))./((M(2:3:end)+M(3:3:end))/2);
% figure(5); plot(B(1:3:end),Delta,'o')
% figure(2); hold on;
% plot(B,M,'o')
% plot([0 B0],[0 M0]/2,'r-')
% plot([0 B0],[0 M0],'b-')
% legend('Mesurée','Polystyrène/2','Polystyrène')
% xlabel('B (mT)','Fontsize',15); ylabel('Aimantation','Fontsize',15)
% hold off
% 
% figure(3); hold on;
% plot(B0,M0,'r-')
% plot(Btab,Muniq,'rd')
% %legend('Mesurée','Polystyrène/2','Polystyrène')
% xlabel('B (mT)','Fontsize',15); ylabel('Aimantation','Fontsize',15)
% 
% 
% for i=1:3:length(All)
% 
% figure(3);
% %plot(B(i),M(i),'bo-')
% errorbar(B(i),(M(i+1)+M(i+2))/2,M(i)-M(i+1),M(i+2)-M(i),'bo-')
% 
% figure(4); hold on;
% plot(B(i),Fmes(i),'bo')
% xlabel('B (mT)'); ylabel('Force pN')
% end

