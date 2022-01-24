names = {'4s_800nm' '8s_800nm' '4s_1600nm' '8s_1600nm'};


depthoname = 'Deptho_DCmed-7.5pc_301018';

Zdifftot = [];
Zreftot = [];
Zcalctot = [];

nfiles = length(names);

for kn = 1:nfiles
    resfile =  ['D:\Data\Raw\18.10.30\Tracking_' names{kn} '_Results.txt'];
    piffile =  ['D:\Data\Raw\18.10.30\Tracking_' names{kn} '_Pifoc.txt'];
    stackname = ['D:\Data\Raw\18.10.30\Tracking_' names{kn} '.tif'];
    
    [VarName1,Area,Mean,StdDev,Min,Max,X,Y,XM,YM,Major,...
        Minor,Angle1,IntDen,RawIntDen,Slice] = importfileDeptho2(resfile);
    
    [Sdown,Smid,Sup] = SplitZimage(Slice',3);
    
    [Smd,ptrd,ptrm] = intersect(Sdown+1,Smid);
    [Smdu,ptrmd,ptru] = intersect(Smd,Sup-1);
    
    
    l = min([length(Sdown) length(Smid) length(Sup)]);
    
    Sdown = Sdown(ptrd(ptrmd));
    Smid  = Smdu ;
    Sup   = Sup(ptru)  ;
    
    % bottom images
    ptrSb = ismember(Slice,Sdown);
    
    Sb = Slice(ptrSb);
    Xb = XM(ptrSb);
    Yb = YM(ptrSb);
    Stdb = StdDev(ptrSb);
    
    % middle images
    ptrSm = ismember(Slice,Smid);
    
    Sm = Slice(ptrSm);
    Xm = XM(ptrSm);
    Ym = YM(ptrSm);
    Stdm = StdDev(ptrSm);
    
    % upper images
    ptrSu = ismember(Slice,Sup);
    
    Su = Slice(ptrSu);
    Xu = XM(ptrSu);
    Yu = YM(ptrSu);
    Stdu = StdDev(ptrSu);
    
    
    Zcalcb= Zcorr(Stdb,Xb,Yb,Sb,stackname,depthoname);
    
    Zcalcb = -Zcalcb*1000;
    [~,posb] = max(Zcalcb);
    
    Zcalcm= Zcorr(Stdm,Xm,Ym,Sm,stackname,depthoname);
    
    Zcalcm = -Zcalcm*1000;
    [~,posm] = max(Zcalcm);
    
    Zcalcu= Zcorr(Stdu,Xu,Yu,Su,stackname,depthoname);
    
    Zcalcu = -Zcalcu*1000;
    [~,posu] = max(Zcalcu);
    
   [Zreffull,Tfull] = importpifoc(piffile);
    
   T = Tfull(Sm);
   
    Zrefb = Zreffull(Sb)*1000;
    Zrefm = Zreffull(Sm)*1000;
    Zrefu = Zreffull(Su)*1000;
    
    Zrefb = Zrefb - Zrefb(posb) + Zcalcb(posb);    
    Zrefm = Zrefm - Zrefm(posm) + Zcalcm(posm);    
    Zrefu = Zrefu - Zrefu(posu) + Zcalcu(posu);
    
    Zdiffb = abs(Zrefb - Zcalcb);
    Zdiffm = abs(Zrefm - Zcalcm);
    Zdiffu = abs(Zrefu - Zcalcu);
    
    figure
    plot(T,Zrefb,'*-')
    hold on
    plot(T,Zcalcb,'-*')
    title(['Mean of diff = ' num2str(mean(Zdiffb)) 'nm'])
    xlabel('T (s)')
    ylabel('Z (nm)')
    legend('Zrefb','Ztrackb')
    
        figure
    plot(T,Zrefm,'*-')
    hold on
    plot(T,Zcalcm,'-*')
    title(['Mean of diff = ' num2str(mean(Zdiffm)) 'nm'])
    xlabel('T (s)')
    ylabel('Z (nm)')
    legend('Zrefm','Ztracm')
    
        figure
    plot(T,Zrefu,'*-')
    hold on
    plot(T,Zcalcu,'-*')
    title(['Mean of diff = ' num2str(mean(Zdiffu)) 'nm'])
    xlabel('T (s)')
    ylabel('Z (nm)')
    legend('Zrefu','Ztracku')
    
%     
%     Zdifftot = [Zdifftot Zdiff'];
%     Zreftot = [Zreftot Zref'];
%     Zcalctot = [Zcalctot Zcalc'];
    
    
    %     figure
    %     plot(Tfull,Zreffull,'*-')
    %     hold on
    %     plot(Tfull,Zcalcfull,'-*')
    %     title(['Mean of diff = ' num2str(mean(Zdifffull*1000)) 'nm'])
    %     xlabel('T (s)')
    %     ylabel('Z (nm)')
    %     fig = gcf;
    %     saveas(fig,['D:\Data\Figures\18-10-29\TrackingFullFig_' names{kn} '.fig'],'fig')
    %     saveas(fig,['D:\Data\Figures\18-10-29\TrackingFullFig_' names{kn} '.png'],'png')
    %     close
    
end
% 
% ptrzerr = (Zdifftot>10)&(Zdifftot<1000)&(Zcalctot>-600);
% 
% Zerr = mean(Zdifftot(ptrzerr));
% 
% figure(1000)
% hold on
% plot(Zcalctot,Zdifftot,'r*')
% plot(Zcalctot(ptrzerr),Zdifftot(ptrzerr),'bo')
% title(['Tracking Precision 181030 - Mean error : ' num2str(Zerr) 'nm'])
% ylabel('Tracking error (nm)')
% xlabel('Absolute Z from focus (nm)')


