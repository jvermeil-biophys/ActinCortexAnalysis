function ZerrGrid = MultiZevaluation(nameExp,nameRef,grid)

% close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',10)

try
load(nameExp)
Kexp = double(K); % Kimo a tester
fexp = f; % focus de test
clear K f
catch
    load(nameExp)
Kexp = double(Ktot); % Kimo a tester
fexp = fmin; % focus de test
clear Ktot fmin
end

try
load(nameRef)
Kref = double(K); % kimo de ref
fref = f; % focus de ref
clear K f
catch
load(nameRef)
Kref = double(Ktot); % kimo de ref
fref = fmin; % focus de ref
clear Ktot fmin
end

% K(1,:) = positiv Z; K(end,:) = negativ Z

Krme = mean(Kref,2);
SizesRef = size(Kref);
Lref = SizesRef(1);

NormRef = repmat(Krme,1,SizesRef(2));

KrefNorm = Kref./NormRef;

SizesExp = size(Kexp);
Lexp = SizesExp(1);

Zstep = 500;

nstep = Zstep/stepnm;

refmax = Lref - fref;
refmin = fref -1;


xbeg = max([nstep+1,fexp-refmin]); % pos de départ
xend = min([Lexp-nstep,fexp+refmax]); % pos de fin

ind = 0;

for x = xbeg:xend
    
    ind = ind+1;
    ZposReal(ind) = (fexp-x)*stepnm;
    
    if x == fexp
        indzero = ind;
    end
    
    clear Mat
    
   
    for i = -1:1
Li = Kexp(x+i*nstep,:);
LiNorm = Li./mean(Li);

MLiNorm = repmat(LiNorm,Lref,1);

DiffK = KrefNorm - MLiNorm;

DiffKsq = DiffK.^2;

DiffKsqsum = sum(DiffKsq,2);

% DiffKsqsum(1) = correlation avec le haut de la bille (anneau)

Mat(:,-i+2) = DiffKsqsum;


clear  Diff* Li* MLi*
    end
   
curvelo = [Mat(1,3)*ones(2*nstep,1); Mat(:,3)];
curvemid = [Mat(1,2)*ones(nstep,1); Mat(:,2); Mat(end,2)*ones(nstep,1)];
curveup = [Mat(:,1);Mat(end,1)*ones(2*nstep,1)];

    
curve = smooth(curvelo+curvemid+curveup);
sup = 1:length(curve);

supint = 1:0.01:length(curve);

curveint = interp1(sup,curve,supint,'spline');

[~,indmin] = min(curveint);
icmi = supint(indmin);


ZposCalc(ind) = (fref-(icmi-nstep))*stepnm;

% 
% figure
% plot(curve,'g')
% hold on
% plot(icmi,cm,'*r')
% plot(curvelo)
% plot(curvemid)
% plot(curveup)
% hold off
    
clear curve*

end

ZposCalc = ZposCalc - ZposCalc(indzero);




figure
hold on
plot([-3000 4000],[0 0],'k--')
plot([0 0],[-3000 4000],'k--')
plot(ZposReal,ZposCalc,'*')
plot(ZposReal,ZposReal,'r')
xlabel('ZposReal','FontSize',20)
ylabel('ZposCalc','FontSize',20)
title(nameRef)

Zerr = ZposReal-ZposCalc;

ZposCalcUnique = unique(ZposCalc);

for i=1:length(ZposCalcUnique)
    ZerrUnique(i) = mean(Zerr(ZposCalc == ZposCalcUnique(i)));
end

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

ZerrGrid = interp1(ZposCalcUnique,ZerrUnique,grid);

