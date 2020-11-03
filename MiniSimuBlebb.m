%% Mini simulation de blebbing quand le cortex est fin

SurfCont = 0.4; % surface de la zone sondée en µm²

Rcell = 10; % rayon de la cellule en micron

SurfCell = 4*pi*Rcell;

Nunit = round(SurfCell/SurfCont); % nombre de zone ou un blebb peut apparaitre


Tsim = 10; % temps de simulation en minute

BlebbSeuil = 0.75; % seuil en fraction de temps ou le cortex est fin

Pblebb = 0.5; % probabilité de blebber quand le cortex est fin

Nblebb = 0; % nombre de blebb sur le temps de simu

for kt = 1:Tsim*6 % un pts temporel toute les 10s
    
    Nblebb = Nblebb + round(sum(rand(1,Nunit)<BlebbSeuil/100)*Pblebb);
    
end


Nblebb