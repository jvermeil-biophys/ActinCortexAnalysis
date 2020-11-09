% optiprep/iodixanol concentration (w/v)
ConcIo = [0 3 5      10      15 20      30 40      50      60];
ConcOp = [0 5 8.3333 16.6666 25 33.3333 50 66.6666 83.3333 100];

% refractive indexes
RefrDC  = [1.335 1.339 1.343 1.350 1.357 1.364 1.379 1.395 1.411 1.426]; % for DC medium sample @37°C
RefrDIC = [1.338 1.343 1.346 1.353 1.360 1.373 1.383 1.398 1.413 1.429]; % for HL5 medium sample @20°C


% Osmolarity measurments
OsmDC  = [309 311 320 324 318 309 318 331 336]; % for DC medium samples

OsmDIC = [221 234 244 255 260 267 272 276 285 336]; % for HL5 medium samples

       
figure
plot(ConcOp,RefrDC,'o','markersize',3,'linewidth',8)
hold on
plot(ConcOp,RefrDIC,'s','markersize',3,'linewidth',8)
title('Refractive index of DC medium and HL5 medium with various Optiprep concentration')
xlabel('Optiprep concentration (% v/v)')
ylabel('Refractive index (nD)')
legend('Dicty medium', 'DC medium','location','northwest')


figure
hold on
% plot(ConcOp,OsmDC,'o','markersize',5,'linewidth',8)
plot(ConcOp,OsmDIC,'s','markersize',5,'linewidth',8)
% title('Osmolarity of DC medium and HL5 medium with various Optiprep concentration')
title('Osmolarity of HL5 medium with various Optiprep concentration')
xlabel('Optiprep concentration (% v/v)')
ylabel('Osmolarity (mosm)')
% legend('Dicty medium', 'DC medium','location','northwest')
legend('Dicty medium','location','northwest')

figure
hold on
% plot(RefrDC,OsmDC,'o','markersize',5,'linewidth',8)
plot(RefrDIC,OsmDIC,'s','markersize',5,'linewidth',8)
% title('Osmolarity vs. Refractive index for DC medium and HL5 medium')
title('Osmolarity vs. Refractive index for HL5 medium')
ylabel('Osmolarity (mosm)')
xlabel('Refractive index (nD)')
% legend('Dicty medium', 'DC medium','location','northwest')
legend('Dicty medium','location','northwest')



