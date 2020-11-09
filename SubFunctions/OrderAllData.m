 ML141
ListLA = {'1' '2' '5'};
ListLA2 = 'ml_';
type = {'DMSO','ML141'};
chmpstr = {'5' '20'};
OrderData(ListLA,ListLA2,type,chmpstr)

% CK666
ListLA = {'LA2' 'LA5' 'LA8'};
ListLA2 = 'LAck_';
type = {'DMSO','CK666'};
chmpstr = {'5' '20'};
OrderData(ListLA,ListLA2,type,chmpstr)

% SMIFH2
ListLA = {'LA2' 'LA3' 'LA8'};
ListLA2 = 'LAsm_';
type = {'DMSO','SMIFH2'};
chmpstr = {'5' '20'};
OrderData(ListLA,ListLA2,type,chmpstr)

% Blebbi
ListLA = {'LA3' 'LA4' 'LA5' 'LA7'};
ListLA2 = 'LAb_';
type = {'DMSO','Blebbi-50','Blebbi-100'};
type2 = {'DMSO','Blebbi','Blebbi'};
chmpstr = {'5' '20'};
OrderData(ListLA,ListLA2,type,chmpstr,type2)

ListLA = {'LA3' 'LA7'};
ListLA2 = 'LAb50_';
type = {'DMSO','Blebbi-50'};
type2 = {'DMSO','Blebbi50'};
chmpstr = {'5' '20'};
OrderData(ListLA,ListLA2,type,chmpstr,type2)

ListLA = {'LA4' 'LA5'};
ListLA2 = 'LAb100_';
type = {'DMSO','Blebbi-100'};
type2 = {'DMSO','Blebbi100'};
chmpstr = {'5' '20'};
OrderData(ListLA,ListLA2,type,chmpstr,type2)

% DMSO
ListLA = {'LA1' 'LA2' 'LA3' 'LA4' 'LA5' 'LA7' 'LA8'};
ListLA2 = 'LA';
type = {'DMSO'};
type2 = {'DMSO'};
chmpstr = {'5' '20'};
OrderData(ListLA,ListLA2,type,chmpstr,type2)


%% Nouvelle manip

% % Ctrl et Ctrl fluo
% ListLA = {'DCLA' 'DCLA'};
% ListLA2 = 'DCLA';
% type = {'CtrlF'};
% type2 = {'_CTRL'};
% chmpstr = {'5'};
% OrderData(ListLA,ListLA2,type,chmpstr,type2)
% 
% % Ctrl et Ctrl fluo
% ListLA = {'DCLA' 'DCLA'};
% ListLA2 = 'DCLA';
% type = {'CtrlNF'};
% type2 = {'_CTRL'};
% chmpstr = {'5'};
% OrderData(ListLA,ListLA2,type,chmpstr,type2)


% % Removal in comp2 
% load('D:\Data\MatFile\V2D\V2D-Flu_15-01-18_M1_R80_MYO1_Comp2.mat')
% NewMR{1} = MR{1};
% NewMR{2} = MR{3};
% NewMR{3} = MR{4};
% NewMR{4} = MR{5};
% NewMR{5} = MR{7};
% NewMR{6} = MR{8};
% NewMR{7} = MR{9};
% NewMR{8} = MR{10};
% MR = NewMR;
% save('D:\Data\MatFile\V2D\V2D-Flu_15-01-18_M1_R80_MYO1_Comp2.mat','MR')