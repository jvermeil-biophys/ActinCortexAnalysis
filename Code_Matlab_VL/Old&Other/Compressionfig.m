%% Choix des conditions
BARE = 1;
DMSO = 1;
LatA = 1;
BLEB = 0;
CK66 = 0;
SMIF = 0;
SERU = 0;
INSI = 1;

% couleurs par conditions
Cct = [0 109 219]./255; % control

Cl  = [73 0 146]./255; % Latranculin

Cb  = [219 209 0]./255; % blebbi
Cy  = [255 255 109 ]./255; % y27
Ccl = [146 0 0]./255; % calyculinA 10nM
Ccl10 = [146 0 0]./255./2; % calyculinA 1nM

Cs  = [36 255 36]./255; % smifh2

Cck = [182 109 255]./255; % ck666

Cin = [0 73 73]./255; % inside
Co  = [255 182 119]./255; % outside

close all

AllForces = [];

figure(1)
figure(2)
T = (1:95)*21; % en ms

%% DC LifeAct compression experiements
if DMSO
    Dists = [];
    Forces = [];
    
    [dist,force] = Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4416);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    [dist,force] =Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4416);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    [dist,force] = Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    [dist,force] = Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    [dist,force] = Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M2','DC-Comp_Ctrl',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    [dist,force] = Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Ctrl',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    figure(1)
    hold on
    Dists = Dists*1000;
    MeanDist = mean(Dists,2)';
    plot(T,MeanDist,'color',Cct,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T,ErrTop,'--','color',Cct,'handlevisibility','off')
    plot(T,ErrBot,'--','color',Cct,'handlevisibility','off')
    fill([T T(end:-1:1) T(1)],...
        [ErrTop ErrBot(end:-1:1) ErrTop(1)],Cct,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'DMSO';
    
    
    
    figure(3)
    hold on
    MeanDist = mean(Dists,2)';
    plot(T(1:48),MeanDist(1:48),'color',Cct,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T(1:48),ErrTop(1:48),'--','color',Cct,'handlevisibility','off')
    plot(T(1:48),ErrBot(1:48),'--','color',Cct,'handlevisibility','off')
    fill([T(1:48) T(48:-1:1) T(1)],...
        [ErrTop(1:48) ErrBot(48:-1:1) ErrTop(1)],Cct,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'DMSO';
    
    figure
    plot(T,Dists,'.','color',Cct,'handlevisibility','off')
    hold on
    plot(T,MeanDist,'color',Cct,'linewidth',5)
    ylabel('Indentation (nm)')
    xlabel('T (ms)')
    leg10 = legend('location','best');
    leg10.String{end} = 'DMSO';
    
    figure(2)
    hold on
    MeanForce = mean(Forces,2)';
    MeanDist = -MeanDist;
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(MeanDist(1:48),MeanForce(1:48),'color',Cct,'linewidth',5)
    plot([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'--','color',Cct,'handlevisibility','off')
    fill([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],Cct,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    plot(MeanDist(48:end),MeanForce(48:end),'r','handlevisibility','off')
    leg2 = legend('location','best');
    leg2.String{end} = 'DMSO';
    
    
    
    
    AllForces = [AllForces Forces];
    
    save('C:\Users\Valentin\Dropbox\TheseValentin\Papiers\MatlabFig\Data\CompData_DMSO.mat','MeanDist','ErrTop','ErrBot');
    
    clear Forces Dists dist force Mean*
end

%% Bare beads compression experiment
if BARE
    Dists = [];
    Forces = [];
    
    [dist,force] = Var2Data_wFluo_Comp('23-01-19','R80',0.56,'M1','DC-Comp_Naked',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    figure(1)
    hold on
    Dists = Dists*1000;
    MeanDist = mean(Dists,2)';
    plot(T,MeanDist,'-k','linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T,ErrTop,'--k','handlevisibility','off')
    plot(T,ErrBot,'--k','handlevisibility','off')
    fill([T T(end:-1:1) T(1)],...
        [ErrTop ErrBot(end:-1:1) ErrTop(1)],'k','facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'Bare beads';
    
    figure(3)
    hold on
    MeanDist = mean(Dists,2)';
    plot(T(1:48),MeanDist(1:48),'color','k','linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T(1:48),ErrTop(1:48),'--','color','k','handlevisibility','off')
    plot(T(1:48),ErrBot(1:48),'--','color','k','handlevisibility','off')
    fill([T(1:48) T(48:-1:1) T(1)],...
        [ErrTop(1:48) ErrBot(48:-1:1) ErrTop(1)],'k','facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'Bare beads';
    
    figure
    plot(T,Dists,'.k','handlevisibility','off')
    hold on
    plot(T,MeanDist,'k','linewidth',5)
    ylabel('Indentation (nm)')
    xlabel('T (ms)')
    leg11 = legend('location','best');
    leg11.String{end} = 'Bare beads';
    
    figure(2)
    hold on
    MeanForce = mean(Forces,2)';
    MeanDist = -MeanDist;
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(MeanDist(1:48),MeanForce(1:48),'k','linewidth',5)
    plot([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'--k','handlevisibility','off')
    fill([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'k','facealpha',0.5,'edgecolor','none','handlevisibility','off');
    plot(MeanDist(48:end),MeanForce(48:end),'r','handlevisibility','off')
    leg2 = legend('location','best');
    leg2.String{end} = 'Bare beads';
    
    
    
    
    
    AllForces = [AllForces Forces];
    
    save('C:\Users\Valentin\Dropbox\TheseValentin\Papiers\MatlabFig\Data\CompData_BareBeads.mat','MeanDist','ErrTop','ErrBot');
    
    clear Forces Dists dist force Mean*
end

%% Latranculin 500nM compression experiements
if LatA
    Dists = [];
    Forces = [];
    
    [dist,force] = Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_LatA500',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    figure(1)
    hold on
    Dists = Dists*1000;
    MeanDist = mean(Dists,2)';
    plot(T,MeanDist,'-','color',Cl,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T,ErrTop,'--','color',Cl,'handlevisibility','off')
    plot(T,ErrBot,'--','color',Cl,'handlevisibility','off')
    fill([T T(end:-1:1) T(1)],...
        [ErrTop ErrBot(end:-1:1) ErrTop(1)],Cl,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'Latranculin A';
    
    figure(3)
    hold on
    MeanDist = mean(Dists,2)';
    plot(T(1:48),MeanDist(1:48),'color',Cl,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T(1:48),ErrTop(1:48),'--','color',Cl,'handlevisibility','off')
    plot(T(1:48),ErrBot(1:48),'--','color',Cl,'handlevisibility','off')
    fill([T(1:48) T(48:-1:1) T(1)],...
        [ErrTop(1:48) ErrBot(48:-1:1) ErrTop(1)],Cl,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'latranculin A';
    
    figure
    plot(T,Dists,'.','color',Cl,'handlevisibility','off')
    hold on
    plot(T,MeanDist,'color',Cl,'linewidth',5)
    ylabel('Indentation (nm)')
    xlabel('T (ms)')
    leg12 = legend('location','best');
    leg12.String{end} = 'Latranculin A';
    
    figure(2)
    hold on
    MeanForce = mean(Forces,2)';
    MeanDist = -MeanDist;
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(MeanDist(1:48),MeanForce(1:48),'color',Cl,'linewidth',5)
    plot([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'--','color',Cl,'handlevisibility','off')
    fill([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],Cl,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    plot(MeanDist(48:end),MeanForce(48:end),'r','handlevisibility','off')
    leg2 = legend('location','best');
    leg2.String{end} = 'Latranculin A';
    
    
    
    
    AllForces = [AllForces Forces];
    
    save('C:\Users\Valentin\Dropbox\TheseValentin\Papiers\MatlabFig\Data\CompData_LatA.mat','MeanDist','ErrTop','ErrBot');
    
    clear Forces Dists dist force Mean*
end

%% Inside beads compression experiements
if INSI
    Dists = [];
    Forces = [];
    
    [dist,force] = Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4416);
    Dists = [Dists dist];
    Forces = [Forces force];
    [dist,force] = Var2Data_wFluo_Comp('28-08-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4416);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    [dist,force] = Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    [dist,force] = Var2Data_wFluo_Comp('17-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    [dist,force] =  Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    [dist,force] = Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    [dist,force] =  Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    [dist,force] =  Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_Inside',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    [dist,force] = Var2Data_wFluo_Comp('12-12-18','R80',0.56,'M1','DC-Comp_Inside',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    
    figure(1)
    hold on
    Dists = Dists*1000;
    MeanDist = mean(Dists,2)';
    plot(T,MeanDist,'-','color',Cin,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T,ErrTop,'--','color',Cin,'handlevisibility','off')
    plot(T,ErrBot,'--','color',Cin,'handlevisibility','off')
    fill([T T(end:-1:1) T(1)],...
        [ErrTop ErrBot(end:-1:1) ErrTop(1)],Cin,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'Inside';
    
    figure
    plot(T,Dists,'.','color',Cin,'handlevisibility','off')
    hold on
    plot(T,MeanDist,'color',Cin,'linewidth',5)
    ylabel('Indentation (nm)')
    xlabel('T (ms)')
    leg12 = legend('location','best');
    leg12.String{end} = 'Inside';
    
    figure(2)
    hold on
    MeanForce = mean(Forces,2)';
    MeanDist = -MeanDist;
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(MeanDist(1:48),MeanForce(1:48),'color',Cin,'linewidth',5)
    plot([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'--','color',Cin,'handlevisibility','off')
    fill([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],Cin,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    plot(MeanDist(48:end),MeanForce(48:end),'r','handlevisibility','off')
    leg2 = legend('location','best');
    leg2.String{end} = 'Inside';
    
    
    
    
    AllForces = [AllForces Forces];
    
    save('C:\Users\Valentin\Dropbox\TheseValentin\Papiers\MatlabFig\Data\CompData_Inside.mat','MeanDist','ErrTop','ErrBot');
    
end

%% Serum beads compression experiment
if SERU
    Dists = [];
    Forces = [];
    
    [dist,force] = Var2Data_wFluo_Comp('23-01-19','R80',0.56,'M1','DC-Comp_SerumJ1',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    figure(1)
    hold on
    Dists = Dists*1000;
    MeanDist = mean(Dists,2)';
    plot(T,MeanDist,'-','color',Co,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T,ErrTop,'--','color',Co,'handlevisibility','off')
    plot(T,ErrBot,'--','color',Co,'handlevisibility','off')
    fill([T T(end:-1:1) T(1)],...
        [ErrTop ErrBot(end:-1:1) ErrTop(1)],Co,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'Serum beads';
    
    figure
    plot(T,Dists,'.g','handlevisibility','off')
    hold on
    plot(T,MeanDist,'g','linewidth',5)
    ylabel('Indentation (nm)')
    xlabel('T (ms)')
    leg11 = legend('location','best');
    leg11.String{end} = 'Serum beads';
    
    figure(2)
    hold on
    MeanForce = mean(Forces,2)';
    MeanDist = -MeanDist;
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(MeanDist(1:48),MeanForce(1:48),'color',Co,'linewidth',5)
    plot([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'--','color',Co,'handlevisibility','off')
    fill([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],Co,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    plot(MeanDist(48:end),MeanForce(48:end),'r','handlevisibility','off')
    leg2 = legend('location','best');
    leg2.String{end} = 'Serum beads';
    
    
    
    
    
    AllForces = [AllForces Forces];
    
    clear Forces Dists dist force Mean*
end

%% Blebbistatin compression experiment
if BLEB
    Dists = [];
    Forces = [];
    
    [dist,force] = Var2Data_wFluo_Comp('30-10-18','R80',0.56,'M1','DC-Comp_Blebbi',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    figure(1)
    hold on
    Dists = Dists*1000;
    MeanDist = mean(Dists,2)';
    plot(T,MeanDist,'-','color',Cb,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T,ErrTop,'--','color',Cb,'handlevisibility','off')
    plot(T,ErrBot,'--','color',Cb,'handlevisibility','off')
    fill([T T(end:-1:1) T(1)],...
        [ErrTop ErrBot(end:-1:1) ErrTop(1)],Cb,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'Blebbi';
    
    figure
    plot(T,Dists,'.','color',Cb,'handlevisibility','off')
    hold on
    plot(T,MeanDist,'color',Cb,'linewidth',5)
    ylabel('Indentation (nm)')
    xlabel('T (ms)')
    leg11 = legend('location','best');
    leg11.String{end} = 'Blebbi';
    
    figure(2)
    hold on
    MeanForce = mean(Forces,2)';
    MeanDist = -MeanDist;
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(MeanDist(1:48),MeanForce(1:48),'color',Cb,'linewidth',5)
    plot([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'--','color',Cb,'handlevisibility','off')
    fill([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],Cb,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    plot(MeanDist(48:end),MeanForce(48:end),'r','handlevisibility','off')
    leg2 = legend('location','best');
    leg2.String{end} = 'Blebbi';
    
    
    
    
    
    AllForces = [AllForces Forces];
    
    clear Forces Dists dist force Mean*
end

%% CK666 compression experiment
if CK66
    Dists = [];
    Forces = [];
    
    [dist,force] = Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M1','DC-Comp_CK666',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    [dist,force] = Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M2','DC-Comp_CK666',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    figure(1)
    hold on
    Dists = Dists*1000;
    MeanDist = mean(Dists,2)';
    plot(T,MeanDist,'-','color',Cck,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T,ErrTop,'--','color',Cck,'handlevisibility','off')
    plot(T,ErrBot,'--','color',Cck,'handlevisibility','off')
    fill([T T(end:-1:1) T(1)],...
        [ErrTop ErrBot(end:-1:1) ErrTop(1)],Cck,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'CK666';
    
    figure
    plot(T,Dists,'.','color',Cck,'handlevisibility','off')
    hold on
    plot(T,MeanDist,'color',Cck,'linewidth',5)
    ylabel('Indentation (nm)')
    xlabel('T (ms)')
    leg11 = legend('location','best');
    leg11.String{end} = 'CK666';
    
    figure(2)
    hold on
    MeanForce = mean(Forces,2)';
    MeanDist = -MeanDist;
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(MeanDist(1:48),MeanForce(1:48),'color',Cck,'linewidth',5)
    plot([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'--','color',Cck,'handlevisibility','off')
    fill([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],Cck,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    plot(MeanDist(48:end),MeanForce(48:end),'r','handlevisibility','off')
    leg2 = legend('location','best');
    leg2.String{end} = 'CK666';   
    
    
    
    
    AllForces = [AllForces Forces];
    
    clear Forces Dists dist force Mean*
end

%% SMIFH2 compression experiment
if SMIF
    Dists = [];
    Forces = [];
    
    [dist,force] = Var2Data_wFluo_Comp('24-09-18','R80',0.56,'M2','DC-Comp_SMIFH2',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    [dist,force] = Var2Data_wFluo_Comp('25-09-18','R80',0.56,'M1','DC-Comp_SMIFH2',95,179,4441);
    Dists = [Dists dist];
    Forces = [Forces force];
    
    
    figure(1)
    hold on
    Dists = Dists*1000;
    MeanDist = mean(Dists,2)';
    plot(T,MeanDist,'-','color',Cs,'linewidth',5)
    Err = std(Dists,[],2)'/sqrt(size(Dists,2));
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(T,ErrTop,'--','color',Cs,'handlevisibility','off')
    plot(T,ErrBot,'--','color',Cs,'handlevisibility','off')
    fill([T T(end:-1:1) T(1)],...
        [ErrTop ErrBot(end:-1:1) ErrTop(1)],Cs,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    leg1 = legend('location','best');
    leg1.String{end} = 'SMIFH2';
    
    figure
    plot(T,Dists,'.','color',Cs,'handlevisibility','off')
    hold on
    plot(MeanDist,'color',Cs,'linewidth',5)
    ylabel('Indentation (nm)')
    xlabel('T (ms)')
    leg11 = legend('location','best');
    leg11.String{end} = 'SMIFH2';
    
    figure(2)
    hold on
    MeanForce = mean(Forces,2)';
    MeanDist = -MeanDist;
    ErrTop = MeanDist+Err;
    ErrBot = MeanDist-Err;
    plot(MeanDist(1:48),MeanForce(1:48),'color',Cs,'linewidth',5)
    plot([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],'--','color',Cs,'handlevisibility','off')
    fill([ErrBot(1:48) ErrTop(48:-1:1) ErrBot(1)],...
        [MeanForce(1:48) MeanForce(48:-1:1) MeanForce(1)],Cs,'facealpha',0.5,'edgecolor','none','handlevisibility','off');
    plot(MeanDist(48:end),MeanForce(48:end),'r','handlevisibility','off')
    leg2 = legend('location','best');
    leg2.String{end} = 'SMIFH2';   
    
    
    
    
    AllForces = [AllForces Forces];
    
    clear Forces Dists dist force Mean*
end

%% Finishing figures

figure(1)
title('Average indentation curve')
ylabel('Indentation (nm)')
yyaxis right
ylabel('Dipolar Force (pN)')
    xlabel('T (ms)')
MeanForce = mean(AllForces,2);
plot(T,MeanForce,'r')
ax = gca;

ax.YColor = [0 0 0];
yyaxis left
leg1 = legend('location','best');
leg1.String{end} = 'Average Force';

figure(3)
title('Average indentation curve')
ylabel('Indentation (nm)')
yyaxis right
ylabel('Dipolar Force (pN)')
    xlabel('T (ms)')
MeanForce = mean(AllForces,2);
plot(T(1:48),MeanForce(1:48),'r')
ax = gca;

ax.YColor = [0 0 0];
yyaxis left
leg3 = legend('location','best');
leg3.String{end} = 'Average Force';


figure(2)
hold on
xlabel('Indentation (nm)')
ylabel('Dipolar Force (pN)')
legend('location','best')
title('Force-Indentation curve')




