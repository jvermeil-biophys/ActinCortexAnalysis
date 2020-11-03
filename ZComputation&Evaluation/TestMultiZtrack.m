clear all

assym = [3 7 10 18 25 35];

for assymup = assym
assymdown = assymup;
dzs = [3 8 15 20 25 30 45 60];

ranges = [2000 Inf];


PLOTperDepth = 0;
PLOTperDay = 0;

ndepth = 0;

for dz = dzs
    AllMeanZPos = [];
Alldzerr =[];
    %% 17-05-18 entre eux
    
    allnum = {'1','2','3','4','5','6','7','8'};
    for k= 1:length(allnum)
        totest = allnum{k};
        
        [MeanZPos{k},dzerr{k}] = MultiZdzEvaluationAssymGuess(['D:\Data\Raw\EtalonageZ\MultiZCorrection\17-05-18_Depthograph_' totest '.mat'],...
            ['D:\Data\Raw\EtalonageZ\MultiZCorrection\17-05-18_Depthograph_Anti-' totest '.mat'],PLOTperDepth,assymup,assymdown,dz);
        
    end
    
    MeanZPos = [MeanZPos{:}];
    dzerr = [dzerr{:}];
    
    AllMeanZPos = [AllMeanZPos MeanZPos];
    Alldzerr = [Alldzerr dzerr];
    
    
    
    if PLOTperDay
        for range = ranges
            ptrrange = abs(MeanZPos)<range;
            
            outlim = round(2*std(dzerr(ptrrange)));
            
            figure
            hold on
            histogram(dzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(dzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(dzerr(ptrrange))<outlim)-sum(abs(dzerr(ptrrange))<outlim))/length(abs(dzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('12-07-18','fontsize',30)
            set(gca,'YTick',[])
        end
        
        figure
        plot(MeanZPos,dzerr,'*')
    end
    
    clear MeanZPos dzerr ptrrange
    
    %% 12-07-18 entre eux
    
    allnum = {'1','2-1','2-2','3'};
    for k= 1:length(allnum)
        totest = allnum{k};
        
        [MeanZPos{k},dzerr{k}] = MultiZdzEvaluationAssymGuess(['D:\Data\Raw\EtalonageZ\MultiZCorrection\12-07-18_Depthograph_' totest '.mat'],...
            ['D:\Data\Raw\EtalonageZ\MultiZCorrection\12-07-18_Depthograph_Anti-' totest '.mat'],PLOTperDepth,assymup,assymdown,dz);
        
    end
    
    MeanZPos = [MeanZPos{:}];
    dzerr = [dzerr{:}];
    
    AllMeanZPos = [AllMeanZPos MeanZPos];
    Alldzerr = [Alldzerr dzerr];
    
    
    
    if PLOTperDay
        for range = ranges
            ptrrange = abs(MeanZPos)<range;
            
            outlim = round(2*std(dzerr(ptrrange)));
            
            figure
            hold on
            histogram(dzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(dzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(dzerr(ptrrange))<outlim)-sum(abs(dzerr(ptrrange))<outlim))/length(abs(dzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('12-07-18','fontsize',30)
            set(gca,'YTick',[])
        end
        
                figure
        plot(MeanZPos,dzerr,'*')
    end
    
    clear MeanZPos dzerr ptrrange
    
    %% 18-07-18 entre eux
    
    allnum = {'1-1','1-2','2-1','2-2'};
    for k= 1:length(allnum)
        totest = allnum{k};
        
        [MeanZPos{k},dzerr{k}] = MultiZdzEvaluationAssymGuess(['D:\Data\Raw\EtalonageZ\MultiZCorrection\18-07-18_Depthograph_' totest '.mat'],...
            ['D:\Data\Raw\EtalonageZ\MultiZCorrection\18-07-18_Depthograph_Anti-' totest '.mat'],PLOTperDepth,assymup,assymdown,dz);
        
    end
    
    MeanZPos = [MeanZPos{:}];
    dzerr = [dzerr{:}];
    
    AllMeanZPos = [AllMeanZPos MeanZPos];
    Alldzerr = [Alldzerr dzerr];
    
    
    
    if PLOTperDay
        for range = ranges
           ptrrange = abs(MeanZPos)<range;
            
            outlim = round(2*std(dzerr(ptrrange)));
            
            figure
            hold on
            histogram(dzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(dzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(dzerr(ptrrange))<outlim)-sum(abs(dzerr(ptrrange))<outlim))/length(abs(dzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('12-07-18','fontsize',30)
            set(gca,'YTick',[])
        end
    end
    clear MeanZPos dzerr ptrrange
    
    %% 31-07-18 entre eux
    
    allnum = {'1','2','3','4'};
    for k= 1:length(allnum)
        totest = allnum{k};
        
        [MeanZPos{k},dzerr{k}] = MultiZdzEvaluationAssymGuess(['D:\Data\Raw\EtalonageZ\MultiZCorrection\31-07-18_Depthograph_' totest '.mat'],...
            ['D:\Data\Raw\EtalonageZ\MultiZCorrection\31-07-18_Depthograph_Anti-' totest '.mat'],PLOTperDepth,assymup,assymdown,dz);
        
    end
    
    MeanZPos = [MeanZPos{:}];
    dzerr = [dzerr{:}];
    
    AllMeanZPos = [AllMeanZPos MeanZPos];
    Alldzerr = [Alldzerr dzerr];
    
    
    
    if PLOTperDay
        for range = ranges
            ptrrange = abs(MeanZPos)<range;
            
            outlim = round(2*std(dzerr(ptrrange)));
            
            figure
            hold on
            histogram(dzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(dzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(dzerr(ptrrange))<outlim)-sum(abs(dzerr(ptrrange))<outlim))/length(abs(dzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('12-07-18','fontsize',30)
            set(gca,'YTick',[])
        end
    end
    clear MeanZPos dzerr ptrrange
    
    %% 28 08 18 entre eux
    
    allnum = {'1','2','3'};
    for k= 1:length(allnum)
        totest = allnum{k};
        
        [MeanZPos{k},dzerr{k}] = MultiZdzEvaluationAssymGuess(['D:\Data\Raw\EtalonageZ\MultiZCorrection\28-08-18_Depthograph_' totest '.mat'],...
            ['D:\Data\Raw\EtalonageZ\MultiZCorrection\28-08-18_Depthograph_Anti-' totest '.mat'],PLOTperDepth,assymup,assymdown,dz);
        
    end
    
    MeanZPos = [MeanZPos{:}];
    dzerr = [dzerr{:}];
    
    AllMeanZPos = [AllMeanZPos MeanZPos];
    Alldzerr = [Alldzerr dzerr];
    
    
    
    if PLOTperDay
        for range = ranges
            ptrrange = abs(MeanZPos)<range;
            
            outlim = round(2*std(dzerr(ptrrange)));
            
            figure
            hold on
            histogram(dzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(dzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(dzerr(ptrrange))<outlim)-sum(abs(dzerr(ptrrange))<outlim))/length(abs(dzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('12-07-18','fontsize',30)
            set(gca,'YTick',[])
        end
    end
    
    clear MeanZPos dzerr ptrrange
    
    %% 30 10 18 entre eux
    
    allnum = {'1','2','3','4','5'};
    for k= 1:length(allnum)
        totest = allnum{k};
        
        [MeanZPos{k},dzerr{k}] = MultiZdzEvaluationAssymGuess(['D:\Data\Raw\EtalonageZ\MultiZCorrection\30-10-18_Depthograph_' totest '.mat'],...
            ['D:\Data\Raw\EtalonageZ\MultiZCorrection\30-10-18_Depthograph_Anti-' totest '.mat'],PLOTperDepth,assymup,assymdown,dz);
        
    end
    
    MeanZPos = [MeanZPos{:}];
    dzerr = [dzerr{:}];
    
    AllMeanZPos = [AllMeanZPos MeanZPos];
    Alldzerr = [Alldzerr dzerr];
    
    
    
    if PLOTperDay
        for range = ranges
            ptrrange = abs(MeanZPos)<range;
            
            outlim = round(2*std(dzerr(ptrrange)));
            
            figure
            hold on
            histogram(dzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(dzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(dzerr(ptrrange))<outlim)-sum(abs(dzerr(ptrrange))<outlim))/length(abs(dzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('12-07-18','fontsize',30)
            set(gca,'YTick',[])
        end
    end
    
    clear MeanZPos dzerr ptrrange
    
    %% 12 12 18 entre eux
    
    allnum = {'1','2-1','2-2','3','4'};
    % allnum = {'7'};
    for k= 1:length(allnum)
        totest = allnum{k};
        
        [MeanZPos{k},dzerr{k}] = MultiZdzEvaluationAssymGuess(['D:\Data\Raw\EtalonageZ\MultiZCorrection\12-12-18_Depthograph_' totest '.mat'],...
            ['D:\Data\Raw\EtalonageZ\MultiZCorrection\12-12-18_Depthograph_Anti-' totest '.mat'],PLOTperDepth,assymup,assymdown,dz);
        
    end
    
    MeanZPos = [MeanZPos{:}];
    dzerr = [dzerr{:}];
    
    AllMeanZPos = [AllMeanZPos MeanZPos];
    Alldzerr = [Alldzerr dzerr];
    
    
    
    if PLOTperDay
        for range = ranges
            ptrrange = abs(MeanZPos)<range;
            
            outlim = round(2*std(dzerr(ptrrange)));
            
            figure
            hold on
            histogram(dzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(dzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(dzerr(ptrrange))<outlim)-sum(abs(dzerr(ptrrange))<outlim))/length(abs(dzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('12-07-18','fontsize',30)
            set(gca,'YTick',[])
        end
    end
    
    clear MeanZPos dzerr ptrrange
    
    %% 15 03 19 entre eux
    
    allnum = {'1','2','3','4','6','7','8','9-1','9-2'};
    % allnum = {'7'};
    for k= 1:length(allnum)
        totest = allnum{k};
        
        [MeanZPos{k},dzerr{k}] = MultiZdzEvaluationAssymGuess(['D:\Data\Raw\EtalonageZ\MultiZCorrection\15-03-19_Depthograph_' totest '.mat'],...
            ['D:\Data\Raw\EtalonageZ\MultiZCorrection\15-03-19_Depthograph_Anti-' totest '.mat'],PLOTperDepth,assymup,assymdown,dz);
        
    end
    
    MeanZPos = [MeanZPos{:}];
    dzerr = [dzerr{:}];
    
    AllMeanZPos = [AllMeanZPos MeanZPos];
    Alldzerr = [Alldzerr dzerr];
    
    
    
    if PLOTperDay
        for range = ranges
            ptrrange = abs(MeanZPos)<range;
            
            outlim = round(2*std(dzerr(ptrrange)));
            
            figure
            hold on
            histogram(dzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(dzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(dzerr(ptrrange))<outlim)-sum(abs(dzerr(ptrrange))<outlim))/length(abs(dzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('12-07-18','fontsize',30)
            set(gca,'YTick',[])
        end
    end
    
    clear MeanZPos dzerr ptrrange
    
    %% mean over all deptho days
    
    for range = ranges
        ptrrange = abs(AllMeanZPos)<range;
            
            outlim = round(2*std(Alldzerr(ptrrange)));
            outlim = 100;
            
            figure
            hold on
            histogram(Alldzerr(ptrrange),'normalization','probability')
            ylabel('% of points')
            xlabel('Error on dz (nm)')
            set(gca,'fontsize',20)
            title(['Mean>0 : ' num2str(mean(abs(Alldzerr(ptrrange)))) 'nm   Outliers (' num2str(outlim) 'nm) : ' num2str((length(abs(Alldzerr(ptrrange))<outlim)-sum(abs(Alldzerr(ptrrange))<outlim))/length(abs(Alldzerr(ptrrange))<outlim)*100) '%   ' ...
                'Range : ' num2str(range) '   dz : ' num2str(dz*20) 'nm   AssymUD : ' num2str(assymup*20) '-' num2str(assymdown*20) 'nm.' ])
            yyaxis right
            ylabel('AllDays','fontsize',30)
            set(gca,'YTick',[])
    end
    clear AllMeanZPos Alldzerr ptrrange
    
end
end
