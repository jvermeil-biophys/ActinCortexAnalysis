function [Zdist,dzerr] = MultiZdzEvaluationAssymGuess(nameExp,nameRef,PLOT,assymup,assymdown,dz)

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
    fref = f;
    clear K f
catch
    load(nameRef)
    Kref = double(Ktot); % kimo de ref
    fref = fmin;
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



nsteps = [assymup 0 -assymdown];


xbeg = max([fexp-fref+assymdown+1 assymdown+1]); % pos de départ
xend = min([fexp-fref+Lref-1-assymup-dz Lexp-assymup-dz]); % pos de fin


ind = 0;

for x = xbeg:xend
    
    ind = ind+1;
    
    x1 = x;
    x2 = x+dz;

   Mat = NaN(Lref,3); 
    
    % correlating line with depthograph
    for i = 1:3
        Li = Kexp(x1+nsteps(i),:);        
        [Curve,MPtmp,MVtmp] = ComputeZCorrel(Li,KrefNorm);
        
        Mat1(:,i) = Curve;  
        % Curve(1) = correlation avec le haut de la bille (anneau)
               
        MultiPos1(i,:) = MPtmp;
        MultiVal1(i,:) = MVtmp;
        
        clear  Li Curve MPtmp MVtmp
        
        Li = Kexp(x2+nsteps(i),:);        
        [Curve,MPtmp,MVtmp] = ComputeZCorrel(Li,KrefNorm);
        
        if 0
            figure(124)
            hold on
            plot(Curve)
        end
        
        Mat2(:,i) = Curve;  
               
        MultiPos2(i,:) = MPtmp;
        MultiVal2(i,:) = MVtmp;
        
        clear  Li Curve MPtmp MVtmp
    end
    
    % sorting correlations with alignement constraints
    [Bool1,Bool2,ValSum1,ValSum2] = deal(NaN(1,8));
    Str = cell(1,8);
    
    for k1 = 1:2
        for k2 = 1:2
            for k3 = 1:2
                k = bi2de([k3-1,k2-1,k1-1])+1;
                Str{k} = [num2str(k1) num2str(k2) num2str(k3)];
                Bool1(k) = (MultiPos1(1,k1)>MultiPos1(2,k2)) &&...
                    (MultiPos1(2,k2)>MultiPos1(3,k3))&&(MultiPos1(1,k1)<35 + MultiPos1(2,k2)) &&...
                    (MultiPos1(2,k2)< 35 +MultiPos1(3,k3)); % saut dans le bon sens
                ValSum1(k) = MultiVal1(1,k1)+MultiVal1(2,k2)+MultiVal1(3,k3);
                Bool2(k) = (MultiPos2(1,k1)>MultiPos2(2,k2)) &&...
                    (MultiPos2(2,k2)>MultiPos2(3,k3))&&(MultiPos2(1,k1)<35+MultiPos2(2,k2)) &&...
                    (MultiPos2(2,k2)<35+MultiPos2(3,k3)); % saut dans le bon sens
                ValSum2(k) = MultiVal2(1,k1)+MultiVal2(2,k2)+MultiVal2(3,k3);
            end
        end
    end
    clear k1 k2 k3
    
    % getting good Zs with aligned correlation
    
    
    [Zup1,Wup1,Zmid1,Wmid1,Zdown1,Wdown1] = Corr2Zs(Bool1,Str,MultiPos1,MultiVal1,ValSum1,fref);
    [Zup2,Wup2,Zmid2,Wmid2,Zdown2,Wdown2] = Corr2Zs(Bool2,Str,MultiPos2,MultiVal2,ValSum2,fref);

    dzup = abs(Zup2-Zup1);
    Wdzup = Wup1+Wup2;
    Zdistup = (Zup1+Zup2)/2;
    dzmid = abs(Zmid2-Zmid1);
    Wdzmid = Wmid1+Wmid2;
    Zdistmid = (Zmid1+Zmid2)/2;
    dzdown = abs(Zdown2-Zdown1);
    Wdzdown = Wdown1+Wdown2;
    Zdistdown = (Zdown1+Zdown2)/2;
    
    dzcalc(ind) = (dzup.*Wdzup+ dzmid.*Wdzmid + dzdown.*Wdzdown)/(Wdzup+Wdzmid+Wdzdown);
    Zdist(ind) = (fref-(x+dz/2))*stepnm;

%     if (dzcalc(ind)-dz)*stepnm >200
%         pause(0.1)
%     end
end


dzerr = (dz-dzcalc)*stepnm;



% P = polyfit(Zdist,dzerr,10);

if PLOT
figure
hold on
plot(Zdist,dzerr,'*')
% plot(sort(Zdist),polyval(P,sort(Zdist)),'b','linewidth',2)
end

end


