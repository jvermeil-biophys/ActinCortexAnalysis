clear MR h0

load('D:\Data\MatFile\GPU_N=300_M=2048_dz=0.01_dx=0.004_a=-8_b=4.9_eps=1.25_rho0=1.18_D=0.001_xi=0.mat')
MR{1}.name = 'NicoActiUp';
MR{1}.class = 'OK';
MR{1}.stack = 'NicoActiUp';
MR{1}.F= h0;
MR{1}.dx= h0;
MR{1}.dy= h0;
MR{1}.dz= h0;
MR{1}.Dyz= h0;
MR{1}.time= 1:length(h0);
MR{1}.S= 1:length(h0);
MR{1}.B= h0;
MR{1}.D2= h0;
MR{1}.D3= h0;
MR{1}.actigood= 1;
save('D:\Data\MatFile\V2D\V2D_16-05-19_M1_NicoActiUp.mat','MR')

clear MR h0

load('D:\Data\MatFile\GPU_N=300_M=2048_dz=0.01_dx=0.004_a=-8_b=4.9_eps=1.25_rho0=1.18_D=0.001_xi=0_n11.mat')
MR{1}.name = 'NicoActiUpVar';
MR{1}.class = 'OK';
MR{1}.stack = 'NicoActiUpVar';
MR{1}.F= h0;
MR{1}.dx= h0;
MR{1}.dy= h0;
MR{1}.dz= h0;
MR{1}.Dyz= h0;
MR{1}.time= 1:length(h0);
MR{1}.S= 1:length(h0);
MR{1}.B= h0;
MR{1}.D2= h0;
MR{1}.D3= h0;
MR{1}.actigood= 1;
save('D:\Data\MatFile\V2D\V2D_16-05-19_M1_NicoActiUpVar.mat','MR')

clear MR h0

load('D:\Data\MatFile\GPU_N=300_M=2048_dz=0.01_dx=0.004_a=-7.5_b=4.9_eps=1.25_rho0=1.18_D=0.001_xi=0.mat')
MR{1}.name = 'NicoActiDown';
MR{1}.class = 'OK';
MR{1}.stack = 'NicoActiDown';
MR{1}.F= h0;
MR{1}.dx= h0;
MR{1}.dy= h0;
MR{1}.dz= h0;
MR{1}.Dyz= h0;
MR{1}.time= 1:length(h0);
MR{1}.S= 1:length(h0);
MR{1}.B= h0;
MR{1}.D2= h0;
MR{1}.D3= h0;
MR{1}.actigood= 1;
save('D:\Data\MatFile\V2D\V2D_16-05-19_M1_NicoActiDown.mat','MR')

clear MR h0

load('D:\Data\MatFile\GPU_N=300_M=2048_dz=0.01_dx=0.004_a=-7.5_b=4.9_eps=1.25_rho0=1.18_D=0.001_xi=0_n11.mat')
MR{1}.name = 'NicoActiDownVar';
MR{1}.class = 'OK';
MR{1}.stack = 'NicoActiDownVar';
MR{1}.F= h0;
MR{1}.dx= h0;
MR{1}.dy= h0;
MR{1}.dz= h0;
MR{1}.Dyz= h0;
MR{1}.time= 1:length(h0);
MR{1}.S= 1:length(h0);
MR{1}.B= h0;
MR{1}.D2= h0;
MR{1}.D3= h0;
MR{1}.actigood= 1;
save('D:\Data\MatFile\V2D\V2D_16-05-19_M1_NicoActiDownVar.mat','MR')

clear MR h0