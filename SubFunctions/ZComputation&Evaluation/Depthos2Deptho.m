% fuse depthograph together
clear all

path = 'D:\Data\Raw\EtalonnageZ\MultiZCorrection';

savepath = 'D:\Data\Raw\EtalonnageZ\';

date = '17-11-20';

num = {'1' '2' '3' '4'}
ndeptho = length(num);

    for kdepth = 1:ndeptho
        
        load([path '\' date '_DepthoM450_' num{kdepth} '.mat'])
        Kim{kdepth} = K;
        Foc(kdepth) = f;
        LengK(kdepth) = size(K,1);
        
    end
    
    LK = unique(LengK);
    if length(LK) ~= 1 
        error('Different kimo sizes')
    end
    
    fmin = min(Foc);
    fmax = max(Foc);
    
    for kdepth = 1:ndeptho
        fcur = Foc(kdepth);
        Kim{kdepth} = Kim{kdepth}(fcur-fmin+1:fcur+LK-fmax,:);

        
    end
    
    Ktot = zeros(size(Kim {1}));
    
    for kdepth = 1:ndeptho
        
        Ktot = Ktot + double(Kim{kdepth});
        
    end
    
    Ktot = Ktot/ndeptho;
    
    clear K f
    
    K = uint16(Ktot);
    f = fmin;
    
    figure
    imshow(imadjust(uint16(Ktot)))
    
    save([savepath '\' date '_DepthographM450.mat'],'K','f','stepnm');
    
    imwrite(uint16(K),[savepath '\' date '_DepthographM450.tif'])
    
    clear K* f*
    

