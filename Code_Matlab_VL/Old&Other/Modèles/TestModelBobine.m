mu0 = 4*pi*10^(-4); % [mT.m/A] perméabilité magnétique du vide
I   = 2;            % [A] courant
N   = 750;          % nb spire


L   = 35*10^(-3) ;  % [m] longeur bobine
E   = 20*10^(-3) ;  % [m] distance du bord de la bobine au point de mesure

R1  = 22*10^(-3);   % [m] rayon interne
R2  = 42*10^(-3);   % [m] rayon externe

R   = (R1+R2)/2 ;   % [m] rayon moyen

% Une bobine actuel, sans coeur, valeur du champ en mT
fun1 = @(n,l,i,r,e) mu0.*n./l.*i./2.*(1./sqrt(r.^2./(l+e).^2+1)-1./sqrt(r.^2./e.^2+1));

fun1(N,L,I,R,E)

Es  = [35:80]*10^(-3);

% figure
hold on
plot(Es*1000,fun1(N,L,I,R,Es))




% 
% D   = 20*10^(-3); % [m] Distance des bobine à x = 0
% 
% % Deux bobines séparée de 2*D centrée sur x = 0, on mesure en Xm
% fun2 = @(n,l,i,r,Xm,D) fun1(n,l,i,r,(D-Xm))+fun1(n,l,i,r,(D+Xm));
% 
% X = [-D*10^3:D*10^3]*10^(-3);
% 
% 
% Ds = [10:50]*10^(-3);
% 
% Ds = D;
% 
% figure
% hold on
% for i = 1:length(Ds)
% Xs = [-Ds(i)*10^3:Ds(i)*10^3]*10^(-3);
% plot(Xs*1000,fun2(N,L,I,R,Xs,Ds(i)))
% end




 

