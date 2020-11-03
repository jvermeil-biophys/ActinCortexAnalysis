function CellHairPlot(datapath,conditions,col)

%% INIT

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',23)


%% load files and get datas

nCond = length(conditions);

list = ls(datapath);
list = list(3:end,:);

nFiles = size(list,1);


ElipsePerim = cell(1,nCond);
IntVar = cell(1,nCond);
IntPeaksSize = cell(1,nCond);

% getting data and making figures
for kC = 1:nCond
    

    ElipsePerim{kC} = [];
    IntVar{kC} = [];
    IntPeaksSize{kC} = [];

    
    for kF = 1:nFiles
        
        if contains(list(kF,:),conditions{kC})
            
            load([datapath filesep list(kF,:)])
            
            ElipsePerim{kC} = [ElipsePerim{kC} MR.ElipsePerim];

            IntVar{kC} = [IntVar{kC} MR.IntVar];
            
             IntPeaksSize{kC} = [IntPeaksSize{kC} MR.IntPeaksSize];
            

        end
    end
    

    
    figure(1)
    hold on
    title('ElipsePerim')
    ylabel('ElipsePerim')
    Spread_Julien(ElipsePerim{kC},kC,0.7,col{kC}) 
    plot([kC-0.35 kC+0.35],[mean(ElipsePerim{kC}) mean(ElipsePerim{kC})],'k','linewidth',2)
    
    if kC ~= 1
       [~,p] = ttest2(ElipsePerim{1},ElipsePerim{kC});
        
       pC(kC-1) = p;
       
        if 5/100 > pC(kC-1) && pC(kC-1) > 1/100
            txtC{kC-1} = ['*']  ;
        elseif 1/100 > pC(kC-1) && pC(kC-1) > 1/1000
            txtC{kC-1} = ['**'];
        elseif 1/1000 > pC(kC-1)
            txtC{kC-1} = ['***'];
        else
            txtC{kC-1} = ['NS'];
        end
    end
    

    
        
    
    figure(2)
    hold on
    title('Intensity Variations')
    ylabel('Intensity Variations')
    Spread_Julien(IntVar{kC},kC,0.7,col{kC}) 
    plot([kC-0.35 kC+0.35],[mean(IntVar{kC}) mean(IntVar{kC})],'k','linewidth',2)
    
    if kC ~= 1
        pS(kC-1) = 1-ttest2(IntVar{1},IntVar{kC});
        
        if 5/100 > pS(kC-1) && pS(kC-1) > 1/100
            txtS{kC-1} = ['*']  ;
        elseif 1/100 > pS(kC-1) && pS(kC-1) > 1/1000
            txtS{kC-1} = ['**'];
        elseif 1/1000 > pS(kC-1)
            txtS{kC-1} = ['***'];
        else
            txtS{kC-1} = ['NS'];
        end
    end
    
    
end


    
[~,Edges] = histcounts(IntPeaksSize{1});




for kC = 1:nCond
    
    
[N{kC},~] = histcounts(IntPeaksSize{kC},'BinEdges',Edges);

end

for kC = 2:nCond
figure
hold on

bar(Edges(1:end-1)+mean(diff(Edges)),N{1}/sum(N{1}),'FaceColor',col{1},'EdgeColor','w',...
            'FaceAlpha',0.5,'barwidth',0.99)
errorbar(Edges(1:end-1)+mean(diff(Edges)),N{1}/sum(N{1}),sqrt(N{1})/sum(N{1}),'k.')

bar(Edges(1:end-1)+mean(diff(Edges)),N{kC}/sum(N{kC}),'FaceColor',col{kC},'EdgeColor','w',...
            'FaceAlpha',0.5,'barwidth',0.99)
errorbar(Edges(1:end-1)+mean(diff(Edges)),N{kC}/sum(N{kC}),sqrt(N{kC})/sum(N{kC}),'k.')

ylim([0 0.5])
end
% stats
yCmax = max([ElipsePerim{:}])*1.05;
ySmax = max([IntVar{:}])*1.05;

yCmin = min([ElipsePerim{:}])*0.95;
ySmin = min([IntVar{:}])*0.95;

yCstep = yCmax*0.03;
ySstep = ySmax*0.03;
   
 figure(1)    
 hold on
 for ii = 2:nCond          
            
 line([1.1 ii-0.1], [yCmax+yCstep yCmax+yCstep],'color','k', 'linewidth', 1.3)
 text((1+ii)/2,yCmax + 2*yCstep,txtC{ii-1},'HorizontalAlignment','center','fontsize',11)

 yCmax = yCmax + 2*yCstep;


 end
ylim([yCmin yCmax*1.05])

 
figure(2)    
 hold on
 for ii = 2:nCond          
            
 line([1.1 ii-0.1], [ySmax+ySstep ySmax+ySstep],'color','k', 'linewidth', 1.3)
 text((1+ii)/2,ySmax + 2*ySstep,txtS{ii-1},'HorizontalAlignment','center','fontsize',11)

 ySmax = ySmax + 2*ySstep;


 end
ylim([ySmin ySmax*1.05])



for k= 1:2
figure(k)
xlim([0 nCond+1])
ax = gca;
ax.XTick = 1:nCond;
ax.XTickLabel = conditions;
end

end