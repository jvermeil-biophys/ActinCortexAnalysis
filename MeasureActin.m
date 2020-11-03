function AnalyzeActin(path,dates,specif,tags,PLOTPEAKS)

%% INIT

warning('off','all')

close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',15)

sf   = 'D:\Data\MatFile\MA';

datenow = datestr(now,'yy-mm-dd');

mkdir(sf)

scf = 'C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Données';

mkdir([scf filesep datenow filesep 'MA'])

% figures saving location
ff = 'D:\Data\Figures';

ff = [ff filesep datenow filesep 'FluoPeaks'];
df = ['C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Figure' filesep 'FluoPeaks'];
mkdir(ff)
mkdir(df)



color{1} = [0 0.7 0.3];
color{2} = [1 0.2 0.5];

AllInts = [];
AllDens = [];
AllSizs = [];

%% Main
for kd = 1:length(dates)
    for kt = 1:length(tags)
        %% Chargement de l'image déplié
        imgname = [path filesep dates{kd} '_' specif '_S' tags{kt} '.tif'];
        if exist(imgname)
            Iinfo = imfinfo(imgname);
            stacksize = length(Iinfo);
            
            Res = Iinfo(1).XResolution;
            Scale = 1000/Res; % in nanometer per pixel
            Scale = 75; % in nanometer per pixel
            
            Itmp = imread([path filesep dates{kd} '_' specif '_S' tags{kt} '.tif'],'tiff',1);
            
            
            [lx,ly] = size(Itmp);
            
            UfI = uint8(zeros([lx,ly,stacksize]));
            
            for ks = 1:stacksize
                UfI(:,:,ks) = imread([path filesep dates{kd} '_' specif '_S' tags{kt} '.tif'],'tiff',ks);
            end
            
            name = [dates{kd} '_' specif '_S' tags{kt}];
            MA{kt}.name = name;
            
            %% Mesure de la quantitée d'actin dans une zone du cortex            
            AllPix = double(UfI(1:end));
            CelPix = AllPix(AllPix>0);
            PixTh = 2*mode(CelPix);
            
            pos = [40 80 120 160 200 240 280 320];
            MA{kt}.pos = pos;
            

            
            
            for kpos = 1:length(pos)
                for ks = 1:stacksize
                    box = double(UfI(:,pos(kpos)-3:pos(kpos)+3,ks));
                    
                    bw = box>PixTh;
                    bw2 = bwpropfilt(bw,'Area',1);
                    
                    boxprofile = mean(box,2);
                    
                    nbpx = round(Res);
                    rem = Res-nbpx;
                    
                    level = (sum(boxprofile(end-nbpx+1:end)) + rem*boxprofile(end-nbpx))/Res;


                    sup = find(boxprofile>level);

                    ActinSizPx{kpos}(ks) = (length(boxprofile) - sup(1));
                    
                    
%                     figure
%                     hold on
%                     plot(boxprofile,'-*')
%                     plot(length(boxprofile)-nbpx:length(boxprofile),...
%                         boxprofile(end-nbpx:end),'ro','linewidth',2)
%                     plot([0 length(boxprofile)],[level level])
%                     plot(sup(1),boxprofile(sup(1)),'rs','linewidth',2)
%                     hold off
                    
                    
             
                    
                    ActinInt(ks)  = sum(box(bw2));
                    ActinDen(ks) = mean(box(bw2));
                    
                end
                
                ActinSiz = (ActinSizPx{kpos}-min(ActinSizPx{kpos}))*Scale;
                
                ActinDen(isnan(ActinDen)) = 0;
                
                ActinIntS = smooth(ActinInt,5,'sgolay',3);
                ActinDenS = smooth(ActinDen,5,'sgolay',3);
                ActinSizS = smooth(ActinSiz,5,'sgolay',3);          
                
                ActinIntN = ActinIntS/(mean(ActinIntS))*100;
                
                IntCentered = (ActinInt - mean(ActinInt))/(mean(ActinInt))*100;
                DenCentered = ActinDen - PixTh;
                
                AllInts = [AllInts IntCentered];
                AllDens = [AllDens DenCentered];
                AllSizs = [AllSizs ActinSiz];
                
                T = 0:length(ActinIntN)-1;
                MA{kt}.time{kpos} = T;
                
                
                [Pks,TpsPks,wPks,pPks] =  findpeaks(ActinDenS,T);

                MA{kt}.curveDen{kpos} = ActinDenS;

                MA{kt}.peaksDen{kpos} = [Pks TpsPks' wPks' pPks];
                
                clear Pks TpsPks wPks pPks
                

                [Pks,TpsPks,wPks,pPks] =  findpeaks(ActinSizS,T,'MinPeakProminence',2*Scale);  

                MA{kt}.curveSiz{kpos} = ActinSizS;

                MA{kt}.peaksSiz{kpos} = [Pks TpsPks' wPks' pPks];
                
                clear Pks TpsPks wPks pPks
                
                
                [Pks,TpsPks,wPks,pPks] =  findpeaks(ActinIntN,T);
                %             findpeaks(ActinIntN,T,'annotate','extents')
                
                
                MA{kt}.curveInt{kpos} = ActinIntN;

                MA{kt}.peaksInt{kpos} = [Pks TpsPks' wPks' pPks];
                
                
                


                
               %                                     findpeaks(ActinIntN,T,'annotate','extents','MinPeakProminence',20)
                %                                     xlabel('T (s)')
                %                                     ylabel('Intensity (norm)')
                %                                     title([name '-' num2str(pos(kpos)) 'deg'])
                %                                     fig = gcf;
                %                                     saveas(fig,[ff filesep name '-FluoCurve'  num2str(pos(kpos)) 'deg_' specif '.png'],'png')
                %
                %                                     saveas(fig,[df filesep name '-FluoCurve'  num2str(pos(kpos)) 'deg_' specif '.fig'],'fig')
                %                                     close
                %             %
                clear ActinInt* T Pks TpsPks pPks wPks ActinDen* ActinSiz ActinSizS
            end
            
            if PLOTPEAKS
                
                
                mkdir([ff filesep name '_Peaks'])
                mkdir([df filesep name '_Peaks'])
                for ks = 1:stacksize
                    fig = figure;
                    imshow(imadjust(UfI(:,:,ks)))
                    hold on
                    title([name '_T = ' num2str(ks)])
                    
                    for kpos = 1:length(pos)
                        plot([pos(kpos)-4 pos(kpos)-4],[0 lx],'--g')
                        plot([pos(kpos)+4 pos(kpos)+4],[0 lx],'--g')
                        
                        plot([pos(kpos) pos(kpos)],[lx lx-ActinSizPx{kpos}(ks)],'r','linewidth',2)
                        
                    end
                    
                    

                    saveas(fig, [ff filesep name '_Peaks' filesep num2str(ks) '.tif'],'tiff');
                    saveas(fig, [df filesep name '_Peaks' filesep num2str(ks) '.tif'],'tiff');

                    close(fig)
                end
            end
            
            clear ActinSizPx
        end        
    end
end

Prob.allints= AllInts;
Prob.alldens= AllDens;
Prob.allsizs= AllSizs;

save([sf filesep 'MA_' specif],'MA','Prob')
save([scf filesep datenow filesep 'MA' filesep 'MA_' specif],'MA','Prob')
end