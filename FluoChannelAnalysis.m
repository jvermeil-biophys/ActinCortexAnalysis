clear all

%% Paths & options
h='C:\Users\Valentin\Dropbox\TheseValentin\Matlab';
RawdataFolder = 'D:\Data\Raw\';
MatfileFolder = 'D:\Data\MatFile';
DropboxDataFolder = 'C:\Users\Valentin\Dropbox\TheseValentin\Données_et_images\Données';
FigureFolder = 'D:\Data\Figures';

matsf = [MatfileFolder filesep 'FCA'];
mkdir(matsf);


%% Get data

% datafolder = [RawdataFolder '19.07.04_ChannelDC\'];
% 
% namelist = {'B_44-Z8' 'B_48-Z7' 'B_52-Z8' 'B_54-Z9' 'B_55-Z7' ...
%     'C_33-Z8' 'C_36-Z9' 'C_39-Z4' 'C_42-Z6' 'C_43-Z8'};
% 
datafolder = [RawdataFolder '19.10.01_ChannelDC\'];

namelist = {'B_22-Z9','B_23-Z7','B_24-Z9','B_28-Z7','B_30-Z10','B_31-Z8','B_34-Z8','B_37-Z7',...
    'C_1-Z8','C_2-Z8','C_4-Z7','C_5-Z7','C_6-Z8','C_8-Z8','C_10-Z7','C_11-Z7','C_12-Z7',...
    'C_14-Z8','C_15-Z8','C_17-Z7','C_19-Z8','C_20-Z8','C_25-Z7','C_29-Z8'};


if exist([matsf filesep 'FluoChannelData.mat'])
    load([matsf filesep 'FluoChannelData.mat']);
end

for k = 1:length(namelist)
    
    if exist('Data','var')
        
        if sum(contains({Data(:).name},namelist{k}))
            ind = find(contains({Data(:).name},namelist{k})==1);
            f = figure('units','normalized','outerposition',[0 0 1 1]);
            hold on
            s1 = subplot(2,2,[3 4]);
            imshow(imadjust(Data(ind).image))
            hold on
            for N = 1:size(Data(ind).areas,1)
                rectangle('position',Data(ind).areas(N,:),'edgecolor','r')
            end
            checkbutton = questdlg('This image has already been analyzed !',['Image :' namelist{k}],...
                'Replace','Skip','Skip');
            close(f)
        else
            ind = k;
            checkbutton = 'Replace';
        end
        
    else
        checkbutton = 'Replace';
        ind = k;
    end

    
    if strcmp(checkbutton,'Replace')
    % récupération du nombre d'images
    Iinfo = imfinfo([datafolder namelist{k} '.tif']);
    Res = Iinfo(1).XResolution; % en pixel/µm
    
    I = imread([datafolder namelist{k} '.tif']);
    
    
    N = 0; % nombre de selection
    
    
    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    s1 = subplot(2,2,[1 2]);
    imshow(imadjust(I));
    button2 = questdlg('How many areas do you want to select ?','Area selection','0','1','2','0');
    narea = str2num(button2);
    
    while narea > 0
        N= N +1;
        s1 = subplot(2,2,[1 2]);
        imshow(imadjust(I));
        
        rect(N,:) = getrect(s1);
        
        rectangle('Position',rect(N,:),'edgecolor','r')
        
        Irect = I(rect(N,2):rect(N,2)+rect(N,4),rect(N,1):rect(N,1)+rect(N,3));
        
        figure(f1);
        s2 = subplot(2,2,3);
        imshow(imadjust(Irect))
        
        figure(f1);
        s3 = subplot(2,2,4);
        Iline{N} = mean(double(Irect),1);
        xpx = 1:length(Iline{N});
        xum{N} = xpx/Res;
        plot(xum{N},Iline{N})
        xlabel('Along the cell (µm)')
        ylabel('Intensity (A.U.)')
        
        clear Irect
        
        saveas(f1,[matsf filesep 'Selection_' namelist{k} '-' num2str(N) '.png'],'png');
        saveas(f1,[matsf filesep 'Selection_' namelist{k} '-' num2str(N) '.fig'],'fig');
        
        cla(s2)
        cla(s3)
        
        narea = narea-1;
    end
    
     f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    s4 = subplot(2,2,[1 2]);
    imshow(imadjust(I));
    button2 = questdlg('Select the cell length','Length selection','OK');
    
        
        rect2 = getrect(s4);
        
        Clength = rect2(3)/Res;
    
    
    close(f1)
    
    
    Data(ind).name = namelist{k};
    Data(ind).areas = rect;
    Data(ind).lines = Iline;
    Data(ind).xlines = xum;
    Data(ind).image = I;
    Data(ind).Res = Res;
    Data(ind).Clength = Clength;
    save([matsf filesep 'FluoChannelData.mat'],'Data')
    
    clear rect rect2 I Res
    end
end

