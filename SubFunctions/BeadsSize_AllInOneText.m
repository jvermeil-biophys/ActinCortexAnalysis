% cleanup
clear all
close all
clc

% definition de parametre basique pour figure
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',20)

% localisation des codes matlab
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab'; % a changer une fois pour correspondre a ton pc


%% %%%%%%%%%%%%%%%%%% Parametre a changer

pixelconv=1/15.8; % consersion pixel micron pour le 100X

TailleT
heo = 2.7;

cd('D:\Data\Raw\20.07.21_TailleM270')
% chemin vers le fichier texte

% nom du fichier texte
filename = 'Size_Results_Valentin.txt'

% Titre du graphe
titlestr = 'M-270 BSA coated, 35mT, in BSA';

%%%%%%%%%%%%%%%%%%%% Fin des parametre a changer



%% Lecture du fichier texte

Mat=dlmread(filename,'\t',1,1); % lit le fichier texte dans un matrice
% a partir de la deuxieme colone et de la deuxieme ligne

fid = fopen(filename); % file identifiant
HeadAll = fgetl(fid); % lecture de la premiere ligne (celle ou sont
% nommé les variables)
fclose(fid); % fermeture de l'id du fichier

HeadCell = textscan(HeadAll,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
% Lecture des différent nom de colone

HeadTable = [HeadCell{:}];


%% Identification des numéros de colone d'interet
Nxm = find(strcmp(HeadTable,'XM'));
Nym = find(strcmp(HeadTable,'YM'));

Nstd = find(strcmp(HeadTable,'StdDev'));
Nmax = find(strcmp(HeadTable,'Max'));
Nmin = find(strcmp(HeadTable,'Min'));

Nslice = find(strcmp(HeadTable,'Slice'));



%% declaration des variable a remplir
AllD = [];
AllStd = [];
AllAlpha = [];


%% boucle sur chacune des images
Slice = Mat(:,Nslice);

for i = 1:max(Slice)
    
    ptrslice = find(Slice == i); % travailler sur l'image courante
    
    Std = Mat(ptrslice,Nstd); % deviation standard d'intensité dans la tache lumineuse
    
    Max = Mat(ptrslice,Nmax); % maximum d'intensité dans la tache
    Min = Mat(ptrslice,Nmin); % minimum d'intensité dans la tache
    
    m = (Std > 5000)&(Max<65534); %&(Min<32768); % Condition a verifier sur std, Max et Min
    
    
    % recuperation des X et Y qui respectent les condition definient dans m
    X = Mat(ptrslice(m),Nxm)*pixelconv;
    Y = Mat(ptrslice(m),Nym)*pixelconv;
    
    % clasement selon les X croissants
    A = sortrows([X,Y]);
    
    X = A(:,1);
    Y = A(:,2);
     
    

    Xmat = repmat(X',length(X),1);
    Dx = X-Xmat; % Distance croisées entre tout les points de X
 


    Ymat = repmat(Y',length(Y),1);
    Dy = Y-Ymat; % Distance croisées entre tout les points de Y



    alpha = rad2deg(atan(Dy./Dx));
    
    D = sqrt(Dx.^2+Dy.^2);
    
    mD = D < TailleTheo+0.5 & D> TailleTheo - 0.5; % condition sur la taille
    mAl = alpha > -20 & alpha < 20; % condition sur l'angle
    
    D = unique(D(mD&mAl)); % Distances triée
    
    
    AllAlpha = ([AllAlpha unique(alpha(mD))']);
    
    AllStd = [AllStd Std'];
    
    AllD = [AllD D'];
    
    clear m* Std alpha D
end

figure
histogram(AllStd)
legend(['n = ' num2str(length(AllStd)) ' Mean = ' num2str(mean(AllStd)) ' Std = ' num2str(std(AllStd))])
ylabel('Count')
xlabel('Intensity Std')
title('All intensity std','interpreter','none')

figure
histogram(AllAlpha)
legend(['n = ' num2str(length(AllAlpha)) ' Mean = ' num2str(mean(AllAlpha)) ' Std = ' num2str(std(AllAlpha))])
ylabel('Count')
xlabel('Angle between beads')
title('All angles','interpreter','none')

figure
hold on
h = histogram(AllD,30,'normalization','probability')
xlim([TailleTheo-0.3 TailleTheo+0.3])
title(['n = ' num2str(length(AllD)) ' - Mean = ' num2str(mean(AllD)) ' µm - Std = ' num2str(std(AllD)*1000) ' nm'])
ylabel('Frequency (%)')
xlabel('Bead Size (µm)')



ax = gca;
ax.LineWidth = 1.5;
ax.YTickLabels = [0 2 4 6 8 10 12 14 16 18]
box on



% cd(h)
