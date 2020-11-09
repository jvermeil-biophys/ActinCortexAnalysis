%% INIT & param
clear all
close all
warning('off','all')
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',17)
set(0,'DefaultAxesFontName','Serif')
set(0,'DefaultAxesFontWeight','normal')

% home
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';

Cct = [0 109 219]./255; % control

Cl  = [73 0 146]./255; % Latranculin

%% data

E = {};

H0= {};

% dossier de données et sauvegarde des données
path = 'D:\Data\MatFile\D2M';

List = ls(path);

for i = 3:size(List,1)
    if contains(List(i,:),'Comp_Ctrl')
        
        fprintf(['Chargement du fichier ' List(i,:) '\n\n'])
        load([path filesep List(i,:)]);
        

        E = {E{:} E0{:}};

        H0= {H0{:} H0abs{:}};

        
    end
end

%% fit data and param
Lres = 150; % nanometer, thickness of residual layer
Eres = 40; % kPa

% data to fit
Efit = vertcat(E{:})/1000; % Etot

H0fit = vertcat(H0{:}); % total thickness
LactFit = H0fit - Lres; % actin thickness


ptrpos =~(H0fit>500&Efit>6);


%% fit global
% 
% fitfunction = fittype(['(x+Lr)*(Ea*(Er))/(x*(Er)+Lr*Ea)'],'problem',{'Lr','Er'}); % x = Lact
% 
% [FitParam Goodness] = fit(LactFit(ptrpos),Efit(ptrpos),fitfunction,'problem',{Lres,Eres},'Start',[1]);
% 
% Fitco=coeffvalues(FitParam);
%                             
% Eact=Fitco(1); % module de l'actine
% 
% 
% % pour plot
% LactPlot = min(LactFit(ptrpos)):1:max(LactFit(ptrpos));
% 
% EtotPlot = FitParam(LactPlot);
% 
% % pour R2
% EtotR2 = FitParam(LactFit(ptrpos));
% 
% R2 = Goodness.rsquare;
% 
% % figure
% 
% figure(1)
% hold on
% % plot([0 10],[Eact Eact],'g--','linewidth',1.5)
% % plot([0 10],[Eres Eres],'--','color',Cl,'linewidth',1.5)
% plot(H0fit(ptrpos),Efit(ptrpos),'ko','markerfacecolor',Cct)
% plot(LactFit(ptrpos)+Lres,EtotR2,'ro','linewidth',2)
% xlabel('Cortex Thickness (nm)')
% ylabel('E')
% grid on
% 
% legend('Data',['Fit - R2 : ' num2str(R2)])


%% fit par cell

fitfunctionCell = fittype(['x*(Ea*(Er))/((x-Lr)*(Er)+Lr*Ea)']); % x = Ltot

EactAll = [];
EresAll = [];
LresAll = [];

for ke = 1:length(E)
    if length(E{ke})>5
        
        EFit = E{ke}/1000;
        
        LtotFit = H0{ke};
        
        [FitParamCell Goodness] = fit(LtotFit,EFit,fitfunctionCell,'Start',[1 10 100]);
        
        FitCoeffsCell = coeffvalues(FitParamCell);
        
        EactFit = FitCoeffsCell(1);
        EresFit = FitCoeffsCell(2); 
        LresFit = FitCoeffsCell(3); 
        
        EactAll = [EactAll EactFit];
        EresAll = [EresAll EresFit];
        LresAll = [LresAll LresFit];
        
        R2 = Goodness.rsquare;
        
        EFitter = FitParamCell([LresFit:1:1000]);
        
        figure
        grid on
        hold on
        plot(LtotFit,EFit,'ko','markerfacecolor',Cct)
        plot([LresFit:1:1000],EFitter,'r--','linewidth',1.5)
        
        title(['Eact = ' num2str(EactFit) 'kPa - Eres = ' num2str(EresFit) 'kPa - Lres = ' num2str(LresFit) 'nm'])
        
        legend({'Data',['Fit - R2 = ' num2str(R2)]})
        
        xlim([100 1000])
        ylim([0 35])
        
    end
    
end