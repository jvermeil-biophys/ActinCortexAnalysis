% names = {'1.3AN-2-1' '1.3AN-2-2' '1.3AN-3' '1.3AN-4-1' '1.3AN-4-2'};
names = {'1.4AN-1' '1.4AN-2' '1.4AN-3' '1.4AN-4'};

% depthoname = 'Deptho_DCmed-7.5pc_1.3AN_251018';
depthoname = 'Deptho_DCmed-7.5pc_1.4AN_251018'

nfiles = length(names);

for kn = 1:nfiles
    resfile =  ['D:\Data\Raw\18.10.25\Tracking_100X_' names{kn} '_Results.txt'];
    piffile =  ['D:\Data\Raw\18.10.25\Tracking_100X_' names{kn} '_Pifoc.txt'];
    stackname = ['D:\Data\Raw\18.10.25\Tracking_100X_' names{kn} '.tif'];
    
    [VarName1,Area,Mean,StdDev,Min,Max,X,Y,XM,YM,Major,...
        Minor,Angle1,IntDen,RawIntDen,Slice] = importfileDeptho2(resfile);
    
    [Sdown,Smid,Sup] = SplitZimage(Slice',3);
    
    [Smd,ptrd,ptrm] = intersect(Sdown+1,Smid);
    [Smdu,ptrmd,ptru] = intersect(Smd,Sup-1);
    
    
    l = min([length(Sdown) length(Smid) length(Sup)]);
    
    Sdown = Sdown(ptrd(ptrmd));
    Smid  = Smdu ;
    Sup   = Sup(ptru)  ;
    
    ptrS = ismember(Slice,Smid);
    
    S = Slice(ptrS);
    X1 = XM(ptrS);
    Y1 = YM(ptrS);
    
    
    Zcalc = Zcorr_multiZ(500,XM,X1,YM,Y1,Slice,S,stackname,Sdown,Smid,Sup,depthoname);
    
    Zcalc =  -Zcalc;
    
    [Zref,T] = importpifoc(piffile);
    
    Zref = Zref(S) - Zref(S(1)) + Zcalc(1);
    
    T = T(S) - T(S(1));
    
    Zdiff = abs(Zref-Zcalc);
    
    
    plot(T,Zref,'*-')
    hold on
    plot(T,Zcalc,'-*')
    title(['Mean of diff = ' num2str(mean(Zdiff*1000)) 'nm'])
    xlabel('T (s)')
    ylabel('Z (nm)')
    legend('Zref','Ztrack')
    fig = gcf;
    saveas(fig,['D:\Data\Figures\TrackingZ\TrackingFig_' names{kn} '.fig'],'fig')
    saveas(fig,['D:\Data\Figures\TrackingZ\TrackingFig_' names{kn} '.png'],'png')
    close
    
        figure(1001)
    hold on
    plot(Zref,Zdiff,'*')
    title('Tracking Precision 181025')
    ylabel('Tracking error (nm)')
xlabel('Absolute Z from focus (nm)')
    
end