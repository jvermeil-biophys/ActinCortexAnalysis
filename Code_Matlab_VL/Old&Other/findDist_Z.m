function dz = findDist_Z_corr(dt,Sxy,H);

dth = dt/2;

Sf1 = Sxy + dth;
Sf2 = Sxy - dth;

xH = linspace(1,length(H),length(H))';
dz = abs(interp1(xH,H,Sf1) - interp1(xH,H,Sf2));

% figure
% hold on
% plot(xH,H)
% plot(Sf2,interp1(xH,H,Sf2),'-*')
% plot(Sf1,interp1(xH,H,Sf1),'-*')
% plot(Sxy,interp1(xH,H,Sxy),'-*')
% hold off

end
