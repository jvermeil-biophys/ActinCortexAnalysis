names = {'2' '3' '4'};


depthoname = 'Deptho_DCmed-7.5pc_301018';


nfiles = length(names);

for kn = 1:nfiles
    resfile =  ['D:\Data\Raw\18.10.31\Tracking_50s_8µm-' names{kn} '_Results.txt'];
    piffile =  ['D:\Data\Raw\18.10.31\Tracking_50s_8µm-' names{kn} '_Pifoc.txt'];
    stackname = ['D:\Data\Raw\18.10.31\Tracking_50s_8µm-' names{kn} '.tif'];
    
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
    
    
    Zcalc = -Zcalc;
    
    
    [Zreffull,Tfull] = importpifoc(piffile);
    
    Zref = Zreffull(S) - Zreffull(S(1)) + Zcalc(1);
   Zreffull = Zreffull(Smid) - Zreffull(Smid(1)) + Zcalc(1);
    
    T = Tfull(S) - Tfull(S(1));
    Tfull = Tfull(Smid) - Tfull(Smid(1));
    Zdiff = abs(Zref-Zcalc);
    
    
    figure
    plot(Tfull,Zreffull,'*-')
    hold on
    plot(T,Zcalc,'-*')
    title(['Mean of diff = ' num2str(mean(Zdiff*1000)) 'nm'])
    xlabel('T (s)')
    ylabel('Z (nm)')
    legend('Zref','Ztrack')
    fig = gcf;
    saveas(fig,['D:\Data\Figures\TrackingZ\TrackingFig_50s_8µm-' names{kn} '.fig'],'fig')
    saveas(fig,['D:\Data\Figures\TrackingZ\TrackingFig_50s_8µm-' names{kn} '.png'],'png')
    close
    
    figure(1000)
    hold on
    plot(Zref,Zdiff,'*')
    title('Tracking Precision 181031')
    ylabel('Tracking error (nm)')
    xlabel('Absolute Z from focus (nm)')
    

    
end