function ZerrGrid = Zevaluation(nameExp,nameRef,grid,PLOT)

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',10)

try
load(nameExp)
Kexp = double(K); % Kimo a tester
fexp = f; % focus de test
stepnmexp = stepnm; % pas du deptho a tester
clear K f stepnm
catch
load(nameExp)
Kexp = double(Ktot); % Kimo a tester
fexp = fmin; % focus de test
stepnmexp = stepnm; % pas du deptho a tester
clear Ktot fmin stepnm
end

try
load(nameRef)
Kref = double(K); % kimo de ref
fref = f; % focus de ref
stepnmref = stepnm; % pas du deptho de ref
clear K f stepnm
catch
load(nameRef)
Kref = double(Ktot); % kimo de ref
fref = fmin; % focus de ref
stepnmref = stepnm; % pas du deptho de ref
clear Ktot fmin stepnm
end

% K(1,:) = positiv Z; K(end,:) = negativ Z

Krme = mean(Kref,2);
SizesRef = size(Kref);
Lref = SizesRef(1);

NormRef = repmat(Krme,1,SizesRef(2));

KrefNorm = Kref./NormRef;

SizesExp = size(Kexp);
Lexp = SizesExp(1);

xbeg = max([1 fexp-75]); % pos de départ
xend = min([fexp+75 Lexp]); % pos de fin

ind = 0;

for x = xbeg:xend
    
    ind = ind+1;
    ZposReal(ind) = (fexp-x)*stepnmexp;
    
    if x == fexp
        indzero = ind;
    end
    
    clear Mat
    
   

L = Kexp(x,:);
LNorm = L./mean(L);

MLNorm = repmat(LNorm,Lref,1);

DiffK = KrefNorm - MLNorm;

DiffKsq = DiffK.^2;

DiffKsqsum = sum(DiffKsq,2);

% DiffKsqsum(1) = correlation avec le haut de la bille (anneau)


curve = smooth(DiffKsqsum);
sup = 1:length(curve);

supint = 1:0.01:length(curve);

curveint = interp1(sup,curve,supint,'spline');

curveintExtended = [max(curveint) curveint max(curveint)];
supintExtended = [0 supint supint(end)+1];

[~, locs] = findpeaks(-smooth(curveintExtended),supintExtended,'minpeakdistance',5,'sortstr','descend');

% figure
% findpeaks(-smooth(curveintExtended),supintExtended,'minpeakdistance',5,'sortstr','descend','annotate','extent')


if length(locs)>1

locs = locs(1:2)-1;

Ztest = (fref-locs)*stepnmref;

Zerrtmp = abs(ZposReal(ind) - Ztest);

[~,itmp] = min(Zerrtmp);

Zcalc = Ztest(itmp);

else
    
  Zcalc = (fref-locs)*stepnmref;
  
end

% Zcalc = (fref-locs(1))*stepnmref;

try
ZposCalc(ind) = Zcalc;
catch
    themall
end

clear Zcalc

% figure(12584)
% plot(curve,'g')
% hold on
% plot(icm,cm,'*r')
% plot(curvelo)
% plot(curvemid)
% plot(curveup)
% pause(0.2)
% hold off
    
clear curve*

end

% ZposCalc = ZposCalc - ZposCalc(indzero);


Zerr = ZposReal-ZposCalc;

ZposCalcUnique = unique(ZposCalc);

for i=1:length(ZposCalcUnique)
    ZerrUnique(i) = mean(Zerr(ZposCalc == ZposCalcUnique(i)));
end

if PLOT
figure
hold on
plot([-3000 4000],[0 0],'k--')
plot([0 0],[-3000 4000],'k--')
plot(ZposReal,ZposCalc,'*')
plot(ZposReal,ZposReal,'r')
xlabel('ZposReal','FontSize',20)
ylabel('ZposCalc','FontSize',20)
title(nameRef)


axes('Position',[.6 .18 .3 .3])
box on
hold on
plot(ZposCalc,Zerr,'-*')
plot(ZposCalcUnique,ZerrUnique,'-*')
xlabel('ZposCalc (nm)')
ylabel('Zerr (Real-Calc, nm)')

axes('Position',[.15 .72 .3 .1])
hold on
imshow(imadjust(uint16(Kref)))

axes('Position',[.15 .62 .3 .1])
hold on
imshow(imadjust(uint16(Kexp)))

hold off
end


ZerrGrid = interp1(ZposCalcUnique,ZerrUnique,grid);

