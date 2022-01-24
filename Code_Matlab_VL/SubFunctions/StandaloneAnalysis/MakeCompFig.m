PLOTFC = 0;
DISTFC = 1;
INDENT = 0;
close all

% Couleurs par drogues
Cct = [0 109 219]./255; % control

Cl  = [73 0 146]./255; % Latranculin

Cb  = [219 209 0]./255; % blebbi
Cy  = [255 255 109 ]./255; % y27
Ccl = [146 0 0]./255; % calyculinA 10nM
Ccl10 = [146 0 0]./255./2; % calyculinA 1nM

Cs  = [36 255 36]./255; % smifh2

Cck = [182 109 255]./255;

Cin = [0 73 73]./255; % inside
Co  = [255 182 119]./255; % outside


% listctrl = [1];
% ncond = {4};
% E{1} = FitComp2Plots('DC-Comp_Naked',PLOTFC,DISTFC,INDENT,'k');
% E{2} = FitComp2Plots('DC-Comp_SerumJ6',PLOTFC,DISTFC,INDENT,'m');
% E{3} = FitComp2Plots('DC-Comp_Ctrl',PLOTFC,DISTFC,INDENT,Cd5);
% E{4} = FitComp2Plots('DC-Comp_LatA500',PLOTFC,DISTFC,INDENT,Cl5);
% E{5} = FitComp2Plots('DC-Comp_Inside',PLOTFC,DISTFC,INDENT,Cin);


listctrl = [1];
ncond = {4};
E{1} = FitComp2Plots('DC-Comp_Ctrl',PLOTFC,DISTFC,INDENT,Cct);
E{2} = FitComp2Plots('DC-Comp_LatA500',PLOTFC,DISTFC,INDENT,Cl);
E{3} = FitComp2Plots('DC-Comp_SiCtrl+LatA',PLOTFC,DISTFC,INDENT,Cin);
E{4} = FitComp2Plots('DC-Comp_SiVim+LatA',PLOTFC,DISTFC,INDENT,Co);



coord = E;
for i = 1:length(E)
    
    coord{i}(:) = i;
    
end

% Pos ref

ll = length(listctrl);

for kl = 1:ll
    
    yCmax(kl) = max([E{listctrl(kl):listctrl(kl)+ncond{kl}}])*1.05;
    
    yCstep(kl) = 0.02/1.1*yCmax(kl);
    
    
end

%% test statistique

for ki = 1:ll
    for k=1:ncond{ki}
        i = listctrl(ki);
        ii = i+k;
        
        
        
        pC = ranksum(E{i},E{ii});
        
        figure(3000)
        
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
            line([i+0.1 ii-0.1], [yCmax(ki)+(2*k-2)*yCstep(ki) yCmax(ki)+(2*k-2)*yCstep(ki)],'color','k', 'linewidth', 1.3,'handlevisibility','off')
            text((i+ii)/2,yCmax(ki) + (2*k-1)*yCstep(ki),txtC,'HorizontalAlignment','center','fontsize',11)
            
            xlim([0 length(E)+1])
        end
    end
end



% sauvergarde des figures
sf = 'D:\Data\Figures';

datenow = datestr(now,'yy-mm-dd');

figfolder = [sf filesep datenow filesep 'Meca' filesep 'GlobalRes'];

mkdir(figfolder)
if DISTFC
    saveas(figure(3000),[figfolder filesep 'ElinDistributions.png'])
    close(figure(3000))
end

if PLOTFC
    saveas(figure(250),[figfolder filesep 'EvH_AllFit.png'])
    close(figure(250))
    saveas(figure(260),[figfolder filesep 'ENvH_AllFit.png'])
    close(figure(260))
end