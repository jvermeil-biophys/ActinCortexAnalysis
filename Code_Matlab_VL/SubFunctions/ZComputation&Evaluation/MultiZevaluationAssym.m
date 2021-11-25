function ZerrGrid = MultiZevaluationAssym(nameExp,nameRef,grid,PLOT,assymup,assymdown)

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

nstepup = nstep+assymup;
nstepdown = nstep+assymdown;

nsteps = [nstepup 0 -nstepdown];

refmax = Lref - fref;
refmin = fref -1;

xbeg = max([nstepdown+1 fexp-75]); % pos de départ
xend = min([fexp+75 Lexp-nstepup]); % pos de fin


ind = 0;

for x = xbeg:xend
    
    ind = ind+1;
    ZposReal(ind) = (fexp-x)*stepnm;
    
    if x == fexp
        indgood = ind;
    end
    
    clear Mat
    
   
    for i = 1:3
Li = Kexp(x+nsteps(i),:);
LiNorm = Li./mean(Li);

MLiNorm = repmat(LiNorm,Lref,1);

DiffK = KrefNorm - MLiNorm;

DiffKsq = DiffK.^2;

DiffKsqsum = sum(DiffKsq,2);

% DiffKsqsum(1) = correlation avec le haut de la bille (anneau)

Mat(:,i) = DiffKsqsum;


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

ZposCalc = ZposCalc - ZposCalc(indgood);





Zerr = ZposReal-ZposCalc;

ZposCalcUnique = unique(ZposCalc);

for i=1:length(ZposCalcUnique)
    ZerrUnique(i) = mean(Zerr(ZposCalc == ZposCalcUnique(i)));
end

if PLOT
figure
hold on
plot([-1500 1500],[0 0],'k--')
plot([0 0],[-1500 1500],'k--')
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

