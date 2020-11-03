function ZerrGrid = MultiZevaluationAssymGuess(nameExp,nameRef,grid,PLOT,assymup,assymdown)

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

Krme = median(Kref,2);
SizesRef = size(Kref);
Lref = SizesRef(1);

NormRef = repmat(Krme,1,SizesRef(2));

KrefNorm = Kref./NormRef;

SizesExp = size(Kexp);
Lexp = SizesExp(1);

nstepup = assymup;
nstepdown = assymdown;

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
        LiNorm = Li./median(Li);
        
        MLiNorm = repmat(LiNorm,Lref,1);
        
        DiffK = KrefNorm - MLiNorm;
        
        DiffKsq = DiffK.^2;
        
        DiffKsqsum = sum(DiffKsq,2);
        
        % DiffKsqsum(1) = correlation avec le haut de la bille (anneau)
        
        Mat(:,i) = DiffKsqsum;
        
        ix = 1:length(DiffKsqsum);
        ixq = 1:0.1:length(DiffKsqsum);
        
        iDiffKsqsum = interp1(ix,smooth(DiffKsqsum),ixq,'spline');
        
        iDiffKsqsumExtended = [max(iDiffKsqsum) iDiffKsqsum max(iDiffKsqsum)];
        ixqExtended = [0 ixq ixq(end)+1];
        
        [val,locs] = findpeaks(-iDiffKsqsumExtended,ixqExtended,'minpeakdistance',5,'sortstr','descend');
        
        if length(locs) < 2
            locs(2) = NaN;
            val(2) = NaN;
        end
        
        MultiPos(i,:) = locs(1:2)-1;
        MultiVal(i,:) = -val(1:2);
        
        clear  Diff* Li* MLi*
    end
    
    
    Bool = deal(nan(1,8));
    Str = cell(1,8);
    
    for k1 = 1:2
        for k2 = 1:2
            for k3 = 1:2
                k = bi2de([k3-1,k2-1,k1-1])+1;
                Str{k} = [num2str(k1) num2str(k2) num2str(k3)];
                Bool(k) = (MultiPos(1,k1)>=MultiPos(2,k2)) &&...
                    (MultiPos(2,k2)>=MultiPos(3,k3)); % saut dans le bon sens
                ValSum(k) = MultiVal(1,k1)+MultiVal(2,k2)+MultiVal(3,k3);
            end
        end
    end
    clear k1 k2 k3
    
    if Bool(1) == 1 % la meilleure correlation est bien alignée, on la garde
        ShiftLo = round(MultiPos(1,1)-MultiPos(2,1));
        ShiftUp = round(MultiPos(2,1)-MultiPos(3,1));
    elseif sum(Bool) == 1 % il y a une seule autre possibilité, on la prend
        n = find(Bool ==1);
        k1 = str2double(Str{n}(1));
        k2 = str2double(Str{n}(2));
        k3 = str2double(Str{n}(3));
        ShiftLo = round(MultiPos(1,k1)-MultiPos(2,k2));
        ShiftUp = round(MultiPos(2,k2)-MultiPos(3,k3));
        clear k1 k2 k3
    elseif sum(Bool) > 1 % il y a plusieur autre possibilité, on regarde
        ptr = find(Bool == 1);
        vect1chg = [2 3 5]; % un seul changement
        mask1chg = ismember(vect1chg,ptr);
        ptr1chg = find(mask1chg == 1);
        if sum(mask1chg) == 1 % il n'y a qu'une seule autre possibilité incluant 1 seule changement, on la prend
            k1 = str2double(Str{vect1chg(ptr1chg)}(1));
            k2 = str2double(Str{vect1chg(ptr1chg)}(2));
            k3 = str2double(Str{vect1chg(ptr1chg)}(3));
            ShiftLo = round(MultiPos(1,k1)-MultiPos(2,k2));
            ShiftUp = round(MultiPos(2,k2)-MultiPos(3,k3));
            clear k1 k2 k3
        elseif sum(mask1chg) == 0 % pas de possibilité a 1 changement
            vect2chg = [4 6 7]; % deux changements
            mask2chg = ismember(vect2chg,ptr);
            ptr2chg = find(mask2chg == 1);
            if sum(mask2chg) == 1 % il n'y a qu'une seule autre possibilité incluant 1 seule changement, on la prend
                k1 = str2double(Str{vect2chg(ptr2chg)}(1));
                k2 = str2double(Str{vect2chg(ptr2chg)}(2));
                k3 = str2double(Str{vect2chg(ptr2chg)}(3));
                ShiftLo = round(MultiPos(1,k1)-MultiPos(2,k2));
                ShiftUp = round(MultiPos(2,k2)-MultiPos(3,k3));
                clear k1 k2 k3
            else % il y a deux possibilité a deux changement, on prend la meilleur correl
            [~,im] = min(ValSum(vect2chg(ptr2chg)));
            k1 = str2double(Str{vect2chg(ptr2chg(im))}(1));
            k2 = str2double(Str{vect2chg(ptr2chg(im))}(2));
            k3 = str2double(Str{vect2chg(ptr2chg(im))}(3));
            ShiftLo = round(MultiPos(1,k1)-MultiPos(2,k2));
            ShiftUp = round(MultiPos(2,k2)-MultiPos(3,k3));
            clear k1 k2 k3 im
            end
        else % plusieur combianaison a un changement possible, on prend celle avec la meilleur correlation totale
            [~,im] = min(ValSum(vect1chg(ptr1chg)));            
            k1 = str2double(Str{vect1chg(ptr1chg(im))}(1));
            k2 = str2double(Str{vect1chg(ptr1chg(im))}(2));
            k3 = str2double(Str{vect1chg(ptr1chg(im))}(3));
            ShiftLo = round(MultiPos(1,k1)-MultiPos(2,k2));
            ShiftUp = round(MultiPos(2,k2)-MultiPos(3,k3));
            clear k1 k2 k3 im
        end
    elseif sum(Bool) == 0
%         %%% que faire ?
%         error('No correct combination of correlation !')
% rien ne marche, on fait comme si de rien n'était ?
        ShiftLo = max([round(MultiPos(1,1)-MultiPos(2,1)) 0]);
        ShiftUp = max([round(MultiPos(2,1)-MultiPos(3,1)) 0]);
    end
    
    
    curvelo = Mat(:,1);  % correlation de l'image du bas de la bille (point blanc qui perd du contraste)
    curveloshift = [curvelo; Mat(end,1)*ones(ShiftLo+ShiftUp,1)];
    
    curvemid = Mat(:,2);
    curvemidshift = [Mat(1,2)*ones(ShiftLo,1); curvemid; Mat(end,2)*ones(ShiftUp,1)];
    
    curveup = Mat(:,3);% correlation de l'image du haut de la bille (anneau)
    curveupshift = [Mat(1,3)*ones(ShiftLo+ShiftUp,1); curveup];
   
    
%     figure
%      hold on
%     plot(curvelo)
%     plot(curvemid)
%     plot(curveup)
%     
    
    curve = curveloshift + curvemidshift + curveupshift;
%     curve = curvelo;
    sup = 1:length(curve);
    
    supint = 1:0.01:length(curve);
    
    curveint = interp1(sup,curvemidshift,supint,'spline');
    
    [~,indmin] = min(curveint);
    icmi = supint(indmin);
    
    
    ZposCalc(ind) = (fref-(icmi-ShiftLo))*stepnm;
%     ZposCalc(ind) = (fref-(icmi-assymdown))*stepnm;
    
    
%     figure(12584)
%     plot(curve,'g')
%     hold on
%     plot(curveloshift)
%     plot(curvemidshift)
%     plot(curveupshift)
%     pause
%     hold off
    
    clear curve*
    
end

% ZposCalc = ZposCalc - ZposCalc(indgood);





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

