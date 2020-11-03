Sfull = 1:max(Slice);
[Sdown,Smid,Sup] = SplitZimage(Sfull,3);

l = min([length(Sdown) length(Smid)  length(Sup)]);

        Sdown = Sdown(1:l);
        Smid = Smid(1:l);
        Sup = Sup(1:l);
        
        Sall = [Sdown;Smid;Sup];
        
        STD = nan(3,length(Sdown));
        
        
        Xd = zeros(size(Sdown));
        [~,pd1,p1d] = intersect(Sdown,Slice);
        Xd(pd1) = XM(p1d);
        STD(1,pd1) = StdDev(p1d);

        
        Xm = zeros(size(Smid));
        [~,pm1,p1m] = intersect(Smid,Slice);
        Xm(pm1) = XM(p1m);
        STD(2,pm1) = StdDev(p1m);

        
        Xu = zeros(size(Sup));
        [~,pu1,p1u] = intersect(Sup,Slice);
        Xu(pu1) = XM(p1u);
        STD(3,pu1) = StdDev(p1u);

        
        
        m = ismember(STD,max(STD));
        
        
        Xall = [Xd; Xm; Xu];

        
        X = Xall(m);
        S = Sall(m);
        
        % redefinition de p après la récuperation des "bon" x
        p = ismember(Slice,S);
               
        Y = Y(p);


Sm = logical(sum(m));

        

Zref = Zref(Smid(Sm)) - Zref(Smid(1));

T = T(Smid(Sm)) - T(Smid(1))/1000;






Z = Zcorr_multiZ(1000,X,Y,S,'D:\Data\Raw\ZCaracterisation\Tracks\FixCellTrack1-1_Z1_250ms',...
    Sdown,Smid,Sup,'Kimo_FixCellB_120318');

% MultZ = Zcorr(StdDev,XM,YM,Slice,'D:\Data\Raw\ZCaracterisation\Tracks\FixCellTrack1-1_Z1_250ms','Kimo_FixCellB_120318');

Z = Z - Z(1);

plot(T,Zref,'*-')
hold on
plot(T,Z,'-*')

Err = mean(abs(Z-Zref));

set(0,'DefaultTextInterpreter','none');

title(['FixCellTrack1-1_Z1_250ms Vs. Kimo_FixCellB_120318. MeanErr = ' num2str(Err) 'µm'])
xlabel('T (s)')
ylabel('Z (µm)')

legend('Zref','Ztrack')

fig = gcf;

saveas(fig,'D:\Data\Figures\FixCellTrack1-1_Z1_250ms_Vs_Kimo_FixCellB_120318.png','png')

close(fig)






Z = Zcorr_multiZ(1000,X,Y,S,'D:\Data\Raw\ZCaracterisation\Tracks\FixCellTrack1-1_Z1_250ms',...
    Sdown,Smid,Sup,'Kimo_DCmed_030218');

Z = Z - Z(1);

plot(T,Zref,'*-')
hold on
plot(T,Z,'-*')

Err = mean(abs(Z-Zref));

set(0,'DefaultTextInterpreter','none');

title(['FixCellTrack1-1_Z1_250ms Vs. Kimo_DCmed_030218. MeanErr = ' num2str(Err) 'µm'])
xlabel('T (s)')
ylabel('Z (µm)')

legend('Zref','Ztrack')

fig = gcf;

saveas(fig,'D:\Data\Figures\FixCellTrack1-1_Z1_250ms_Vs_Kimo_DCmed_030218.png','png')

close(fig)




Z = Zcorr_multiZ(1000,X,Y,S,'D:\Data\Raw\ZCaracterisation\Tracks\FixCellTrack1-1_Z1_250ms',...
    Sdown,Smid,Sup,'Kimo_FixCellC_120318');

Z = Z - Z(1);

plot(T,Zref,'*-')
hold on
plot(T,Z,'-*')

Err = mean(abs(Z-Zref));

set(0,'DefaultTextInterpreter','none');

title(['FixCellTrack1-1_Z1_250ms Vs. Kimo_FixCellC_120318. MeanErr = ' num2str(Err) 'µm'])
xlabel('T (s)')
ylabel('Z (µm)')

legend('Zref','Ztrack')

fig = gcf;

saveas(fig,'D:\Data\Figures\FixCellTrack1-1_Z1_250ms_Vs_Kimo_FixCellC_120318.png','png')

close(fig)
