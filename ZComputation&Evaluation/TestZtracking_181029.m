names = {'0.3um_2s_1' '0.3um_2s_2' '0.3um_5s_1' '0.3um_5s_2' '0.3um_10s_1' '0.3um_10s_2'...
    '0.8um_2s_1' '0.8um_2s_2' '0.8um_5s_1' '0.8um_5s_2' '0.8um_10s_1' '0.8um_10s_2'...
    '1.5um_2s_1' '1.5um_2s_2' '1.5um_5s_1' '1.5um_5s_2' '1.5um_10s_1' '1.5um_10s_2'...
    };

% names = {'0.3um_5s_2'}

depthoname = 'Deptho_DCmed-7.5pc_1.4AN_291018';


nfiles = length(names);

for kn = 1:nfiles
    resfile =  ['D:\Data\Raw\18.10.29\Tracking_' names{kn} '_Results.txt'];
    piffile =  ['D:\Data\Raw\18.10.29\Tracking_' names{kn} '_Pifoc.txt'];
    stackname = ['D:\Data\Raw\18.10.29\Tracking_' names{kn} '.tif'];
    
    [VarName1,Area,StdDev,Min,Max,X,Y,XM,YM,Slice] = importfileDeptho(resfile);
    
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
    ylabel('Tracking error (nm)')
    xlabel('Absolute Z from focus (nm)')
end