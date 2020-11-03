% test du Z au 100X du 03-03-2020

% 500 nm de saut a differente hauteur absolue


%% lame et lamelle
names = {'1-1' '1-2' '2-1' '2-2' '3-1' '3-2'};

Alldz = [];
AllZ = [];
Alldzm = [];

for k = 1:length(names)
    
    [VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.03.03\LameLamelle_' names{k} '_Results']);
    
    Z1 = ZCalc_Zcorr_100X(X1,Y1,S,['D:\Data\Raw\20.03.03\LameLamelle_' names{k}],'02-03-20_Depthograph100X','D:\Data\Raw');
    
    Sfull = 1:S(end);
    AntiS = Sfull(not(ismember(Sfull,S)));
    
    figure
    hold on
    plot(S,Z1,'-*')
    plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')
    title(names{k})
    
    
    Sfull = 1:max(S);
    
    [Sdown,Smid,Sup] = SplitZimage(Sfull,3);
    
    Sdown = intersect(Sdown,S);
    Smid  = intersect(Smid,S);
    Sup   = intersect(Sup,S);
    
    [Smd,ptrmid] = intersect(Sdown+1,Smid);
    [Smdu,ptruimd,ptrmdiu] = intersect(Smd,Sup-1);
    
    Sdown = Sdown(ptrmid(ptruimd));
    Smid  = Smdu ;
    Sup   = Sup(ptrmdiu)  ;
    
    [~,~,ptrd] = intersect(Sdown,S);
    [~,~,ptrm]  = intersect(Smid,S);
    [~,~,ptru]   = intersect(Sup,S);
    
    Zd = Z1(ptrd);
    Zm = Z1(ptrm);
    Zu = Z1(ptru);
    
    dzd = Zm - Zd;
    dzu = Zu - Zm;
    
   dzm = diff(Zm);
    
    Alldz = [Alldz dzd dzu];
    Alldzm = [Alldzm dzm];
    AllZ = [AllZ Zm Zm];
    
    figure
    subplot(1,4,[1 3])
    hold on
    plot([min(Zm)-200 max(Zm)+200],[500 500],'k--','linewidth',2)
    plot(Zm,dzd,'b*')
    plot(Zm,dzu,'r*')
    xlim([min(Zm)-200 max(Zm)+200])
    legend({'Expected jump size' 'Bottom-mid' 'Mid-up'})
    ylabel('Jump size (nm)')
    xlabel('Position of midlle image from focus (nm)')
    title(['Lame et lamelle ' names{k}])
   
    subplot(1,4,4)
    hold on
    title('Precision sur la position médiane')
    histogram(dzm,'binwidth',5)
    xlim([-200 200])
    
end



figure
hold on
title('All jumps for Lame et Lamelle')
plot(AllZ,Alldz,'k*')
plot([min(AllZ)-200 max(AllZ)+200],[500 500],'r--','linewidth',1)
xlim([min(AllZ)-200 max(AllZ)+200])

figure
hold on
title(['Distribution of jumps for Lame et Lamelle.      Mean : ' num2str(mean(Alldz)) '        Std err : ' num2str(std(Alldz)/sqrt(length(Alldz)))])
histogram(Alldz)

figure
hold on
title(['Distribution of median pos differential for Lame et Lamelle. '])
histogram(Alldzm,'binwidth',2)
xlim([-200 200])


%% lamelle + PDMS
names = {'2' '3' '4' '5' '6' '7' '8'};

Alldz = [];
AllZ = [];
Alldzm = [];

for k = 1:length(names)
    
    [VarName1, Area, StdDev, X1, Y1, S]  = importfile04022020(['D:\Data\Raw\20.03.03\LamellePDMS_' names{k} '_Results']);
    
    Z1 = ZCalc_Zcorr_100X(X1,Y1,S,['D:\Data\Raw\20.03.03\LamellePDMS_' names{k}],'02-03-20_Depthograph100X','D:\Data\Raw');
    
    Sfull = 1:S(end);
    AntiS = Sfull(not(ismember(Sfull,S)));
    
    figure
    hold on
    plot(S,Z1,'-*')
    plot(AntiS,ones(1,length(AntiS))*mean(Z1),'or')
    title(names{k})
    
    
    Sfull = 1:max(S);
    
    [Sdown,Smid,Sup] = SplitZimage(Sfull,3);
    
    Sdown = intersect(Sdown,S);
    Smid  = intersect(Smid,S);
    Sup   = intersect(Sup,S);
    
    [Smd,ptrmid] = intersect(Sdown+1,Smid);
    [Smdu,ptruimd,ptrmdiu] = intersect(Smd,Sup-1);
    
    Sdown = Sdown(ptrmid(ptruimd));
    Smid  = Smdu ;
    Sup   = Sup(ptrmdiu)  ;
    
    [~,~,ptrd] = intersect(Sdown,S);
    [~,~,ptrm]  = intersect(Smid,S);
    [~,~,ptru]   = intersect(Sup,S);
    
    Zd = Z1(ptrd);
    Zm = Z1(ptrm);
    Zu = Z1(ptru);
    
    dzd = Zm - Zd;
    dzu = Zu - Zm;
    
   dzm = diff(Zm);
    
    Alldz = [Alldz dzd dzu];
    Alldzm = [Alldzm dzm];
    AllZ = [AllZ Zm Zm];
    
    figure
    subplot(1,4,[1 3])
    hold on
    plot([min(Zm)-200 max(Zm)+200],[500 500],'k--','linewidth',2)
    plot(Zm,dzd,'b*')
    plot(Zm,dzu,'r*')
    xlim([min(Zm)-200 max(Zm)+200])
    legend({'Expected jump size' 'Bottom-mid' 'Mid-up'})
    ylabel('Jump size (nm)')
    xlabel('Position of midlle image from focus (nm)')
    title(['Lamelle PDMS ' names{k}])
   
    subplot(1,4,4)
    hold on
    title('Precision sur la position médiane')
    histogram(dzm,'binwidth',5)
    xlim([-200 200])
    
end



figure
hold on
title('All jumps for Lamelle et PDMS')
plot(AllZ,Alldz,'k*')
plot([min(AllZ)-200 max(AllZ)+200],[500 500],'r--','linewidth',1)
xlim([min(AllZ)-200 max(AllZ)+200])

figure
hold on
title(['Distribution of jumps for Lamelle et PDMS.      Mean : ' num2str(mean(Alldz)) '        Std err : ' num2str(std(Alldz)/sqrt(length(Alldz)))])
histogram(Alldz)

figure
hold on
title(['Distribution of median pos differential for Lamelle et PDMS'])
histogram(Alldzm,'binwidth',2)
xlim([-200 200])
