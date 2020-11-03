function FigureMaker(Param,datafolder,figurefolder)
close all

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',27)

% figures saving location
ff = 'D:\Data\Figures';

date = datestr(now,'yy-mm-dd');
ff = [ff filesep date];
df = ['C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Figure'];
mkdir(ff)


% Couleurs par drogues
Cb5  = [255 128 0  ]./255;
Cy5 = [200 200  0  ]./255;
Cp5  = [210 160 0  ]./255;

Cca5 = [0 128 255]./255;

Cs5  = [0   255  90]./255;

Cm5  = [200 100 80  ]./255;
Cc5  = [255 20   0  ]./255;

Cd5  = [160 160 160]./255;

Caw5  = [140 140 140]./255;
Cak5 = [0 255 255]./255;


Cl5 = [110 50 120]./255;
CL5 = [90 40 70]./255;

Stats = 1;
% Param = 0; % dans les input de la fonction




if Param == 1
Titlestr = 'CortexWidth-Param';
elseif Param == 0    
Titlestr = 'CortexWidth-NonParam';
end
listctrl = [1 6];
ncond = {6 1};
Cw{1}     = Data2CortexStats('DCLA_Ctrl'   ,'Ctrl'     ,1 ,Cd5);
Cw{2}     = Data2CortexStats('DCLA_LatA500','LatA500'  ,2 ,Cl5);
Cw{3}     = Data2CortexStats('DCLA_CK666'  ,'CK666'    ,3 ,Cc5);
Cw{4}     = Data2CortexStats('DCLA_ML141'  ,'ML141'    ,4 ,Cm5);
Cw{5}     = Data2CortexStats('DCLA_SMIFH2' ,'SMIFH2'   ,5 ,Cs5);
Cw{6}     = Data2CortexStats('DCLA_Blebbi' ,'Blebbi'   ,6 ,Cb5);
Cw{7}     = Data2CortexStats('DCLA_Y27'    ,'Y27'      ,7 ,Cp5);




% if Param == 1
% Titlestr = 'CortexWidthNase-Param';
% elseif Param == 0    
% Titlestr = 'CortexWidthNase-NonParam';
% end
% listctrl = [1 3];
% ncond = {1 1};
% linectrl = 1;
% Cw{1} = Data2CortexStats('DCLA_Ctrl'    ,'Ctrl'     ,1 ,Cd5);
% Cw{2} = Data2CortexStats('DCLA_Nase'    ,'Nase'     ,2 ,Cp5);
% Cw{3} = Data2CortexStats('DCLA_LatA500' ,'LatA500'  ,3 ,Cl5);
% Cw{4} = Data2CortexStats('DCLA_LatNase' ,'LatNase'   ,4 ,Cp5);
% 

% 
% if Param == 1
% Titlestr = 'CortexWidth-Inside-Param';
% elseif Param == 0    
% Titlestr = 'CortexWidth-Inside-NonParam';
% end
% listctrl = [1 2];
% ncond = {2 1};
% linectrl = 1;
% Cw{1}   = Data2CortexStats('DCLA_Ctrl'   ,'Ctrl'     ,1 ,Cd5 );
% Cw{2}    = Data2CortexStats('DC_Inside_7-5pc','Inside75'  ,2 ,Caw5);
% Cw{3}    = Data2CortexStats('DC_Inside_0pc','Inside75'  ,3 ,Cak5);
% % % 

% 
% 
% if Param == 1
% Titlestr = 'CortexWidth-BeadsSizes-Param';
% elseif Param == 0    
% Titlestr = 'CortexWidth-BeadsSizes-NonParam';
% end
% listctrl = [1 2];
% ncond = {2 1};
% linectrl = 1;
% Cw{1}    = Data2CortexStats('DC-Comp_Naked','Naked'  ,1 ,'k');
% Cw{2}    = Data2CortexStats('DC-Comp_SerumJ6','SerumJ6'  ,2 ,'b');
% % 
% 
% 
% if Param == 1
% Titlestr = 'CortexWidth-FibroLat-Param';
% elseif Param == 0    
% Titlestr = 'CortexWidth-FibroLat-NonParam';
% end
% listctrl = [1 3];
% ncond = {3 1};
% linectrl = 1;
% Cw{1}    = Data2CortexStats('DCLA_Ctrl'   ,'Ctrl'     ,1 ,Cd5 );
% Cw{2}    = Data2CortexStats('DCLA_LatA500' ,'Lat'  ,2 ,Cl5);
% Cw{3}    = Data2CortexStats('MFT16_Vim+Lat' ,'Vim+Lat'  ,3 ,Cak5);
% Cw{4}    = Data2CortexStats('MFT16_NoVim+Lat' ,'NoVim+Lat'  ,4 ,Caw5);
% 


% Param interne de plot
coord = Cw;
for i = 1:length(Cw)
    
    coord{i}(:) = i;
    
end

% Pos ref

ll = length(listctrl);

for kl = 1:ll
    yCmax(kl) = max([Cw{listctrl(kl):listctrl(kl)+ncond{kl}}])*1.05;
yCstep(kl) = 0.02/1.1*yCmax(kl);
end

%% test statistique
if Stats
    for ki = 1:ll
        for k=1:ncond{ki}
            i = listctrl(ki);
            ii = i+k;
            
            
            if Param
            pC = ttest2(Cw{i},Cw{ii});
            pC = 1-pC;
            else
                pC = ranksum(Cw{i},Cw{ii});
            end
            
            figure(1000)
            
            TxtC = 1;
            
            
            if 5/100 > pC && pC > 1/100
                txtC = ['*']  ;
            elseif 1/100 > pC && pC > 1/1000
                txtC = ['**'];
            elseif 1/1000 > pC
                txtC = ['***'];
            else
                txtC = ['NS'];
                TxtC = 1;
            end
            
            if TxtC
                line([i+0.1 ii-0.1], [yCmax(ki)+(2*k-2)*yCstep(ki) yCmax(ki)+(2*k-2)*yCstep(ki)],'color','k', 'linewidth', 1.3)
                text((i+ii)/2,yCmax(ki) + (2*k-1)*yCstep(ki),txtC,'HorizontalAlignment','center','fontsize',11)

                xlim([0 length(Cw)+1])
            end
        end
    end
end

xlim([0 length(Cw)+1])
plot([0 length(Cw)+1],[0 0],'k--')

ymin = min([Cw{:}])-20;

ylim([ymin max(yCmax+(2*[ncond{:}].*yCstep))])

ax = gca;
xtic = ax.XTick;
xlab = ax.XTickLabel;

ax.XTick = [0 xtic];
ax.XTickLabel = ['Condition :\newline     n ='; xlab]

savenamepng = [Titlestr '.png'];
savenamefig = [Titlestr '.fig'];

figure(1000)
title([Titlestr])
    fig = gcf;
    saveas(fig,[ff filesep savenamepng],'png')
    saveas(fig,[df filesep savenamepng],'png')    
    saveas(fig,[ff filesep savenamefig],'fig')
    saveas(fig,[df filesep savenamefig],'fig')
% close all
end
