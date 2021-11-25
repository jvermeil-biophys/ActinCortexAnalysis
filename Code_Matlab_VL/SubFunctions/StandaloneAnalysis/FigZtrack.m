Zref = Zref(Slice) - Zref(Slice(1));

T = (T(Slice) - T(Slice(1)))/1000;

Z = Zcorr(StdDev,XM,YM,Slice,'D:\Data\Raw\ZCaracterisation\Tracks\FixCellTrack1-C_Z0','Kimo_FixCellB_120318');

Z = Z - Z(1);

plot(T,Zref,'*-')
hold on
plot(T,Z,'-*')

Err = mean(abs(Z-Zref));

set(0,'DefaultTextInterpreter','none');

title(['FixCellTrack1-C_Z0 Vs. Kimo_FixCellB_120318. MeanErr = ' num2str(Err) 'µm'])
xlabel('T (s)')
ylabel('Z (µm)')

legend('Zref','Ztrack')

fig = gcf;

saveas(fig,'D:\Data\Figures\FixCellTrack1-C_Z0_Vs_Kimo_FixCellB_120318.png','png')

close(fig)


Z = Zcorr(StdDev,XM,YM,Slice,'D:\Data\Raw\ZCaracterisation\Tracks\FixCellTrack1-C_Z0','Kimo_DCmed_030218');

Z = Z - Z(1);

plot(T,Zref,'*-')
hold on
plot(T,Z,'-*')

Err = mean(abs(Z-Zref));

set(0,'DefaultTextInterpreter','none');

title(['FixCellTrack1-C_Z0 Vs. Kimo_DCmed_030218. MeanErr = ' num2str(Err) 'µm'])
xlabel('T (s)')
ylabel('Z (µm)')

legend('Zref','Ztrack')

fig = gcf;

saveas(fig,'D:\Data\Figures\FixCellTrack1-C_Z0_Vs_Kimo_DCmed_030218.png','png')

close(fig)




Z = Zcorr(StdDev,XM,YM,Slice,'D:\Data\Raw\ZCaracterisation\Tracks\FixCellTrack1-C_Z0','Kimo_FixCellC_120318');

Z = Z - Z(1);

plot(T,Zref,'*-')
hold on
plot(T,Z,'-*')

Err = mean(abs(Z-Zref));

set(0,'DefaultTextInterpreter','none');

title(['FixCellTrack1-C_Z0 Vs. Kimo_FixCellC_120318. MeanErr = ' num2str(Err) 'µm'])
xlabel('T (s)')
ylabel('Z (µm)')

legend('Zref','Ztrack')

fig = gcf;

saveas(fig,'D:\Data\Figures\FixCellTrack1-C_Z0_Vs_Kimo_FixCellC_120318.png','png')

close(fig)
