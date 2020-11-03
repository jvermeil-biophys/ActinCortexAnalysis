%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Assemble les depthographes faits par DepthoMaker.m pour une journée de
% manip en un depthograph final
% 
% A faire pour chaque jour de manip ou des depthographs intermediaires ont été
% faits. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% GeneralFolder = 'C:\Users\BioMecaCell\Desktop\ConstantFieldAnalysis_Joseph\'; % a changer suivant le stockage des données
GeneralFolder = 'D:\PMMH\Matlab Code\ConstantFieldAnalysis_Joseph\';
% GeneralFolder = 'C:\Users\BioMecaCell\Desktop\Deptho_test_pifoc\';

% le reste de la structure de dossier est a respecter

path = [GeneralFolder 'Data_Joseph\Raw\EtalonnageZ\Intermediate'];

savepath = [GeneralFolder  'Data_Joseph\Raw\EtalonnageZ'];

date = '08-07-20'; % date de la manip

num = {'1' '2' '3' '4'}; % liste des Depthographs intermediaires
ndeptho = length(num);

    for kdepth = 1:ndeptho
        
        load([path '\' date '_Deptho63X_P1_FallenBead_' num{kdepth} '.mat'])
        % load([path '\' date '_Depthograph_' num{kdepth} '.mat'])
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
        try
        Kim{kdepth} = Kim{kdepth}(fcur-fmin+1:fcur+LK-fmax,:);
        catch
            aaaa=1
        end
        
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
    
    save([savepath '\' date '_Depthograph63x.mat'],'K','f','stepnm');
    
    imwrite(uint16(K),[savepath '\' date '_Depthograph63x.tif'])
    
    clear K* f*
    

