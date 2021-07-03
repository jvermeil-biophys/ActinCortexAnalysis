function Res2Var_wFluo_multiZ__MultiBeadSize_Comp(restype,date,tag,manip,specif,fich,depthonameSB,depthonameLB,tableExperimentalConditions,nimgbr,nimgr,nimg,wFluo,...
    AUTO,datafolder,resfolder)

%       Res2Var is using the imageJ particle analysis data to create
%       individual trajectories for the beads. It is using the txt files
%       created by ImageJ to get the X and Y positions the beads, it then
%       loads the images of the video and computes the DZ between the beads
%       using the depthograph. This particular version (_wFluo_multiZ_Comp)
%       is used for analyzing compression experiments, where the constant
%       field part between compression have triplet of images in Z. It can
%       handle the presence of fluo images at the end of each loop.
%       This code is specific for
%       experiments made with 100X oil objectiv (see end of code where
%       optical indices correction is applied)
%
%       Res2Var_wFluo_multiZ_Comp(restype,date,tag,manip,specif,fich,
%       depthoname,nimgbr,nimgr,nimg,...
%       AUTO,datafolder,resfolder)
%
%       restype : extension of text file (ex: 'Results' or  'Inside')
%       date : date of experiment 'yy-mm-dd'
%       tag : tag present on the tif and txt files (ex: 'R40', 'R90_Multi')
%       manip : 'manip' number (ex: 'M1' ou 'M2' etc...),
%       specif : experimental condiion (ex: 'Dicty_WT')
%       fich : chambers to analyze for the condition (ex: 1 or [1:3])
%       depthoname : depthograph (ex:'20-05-27_DepthographM450')
%       nimgbr : nb img in a loop before a ramp
%       nimgr : nb img in a ramp
%       nimg : nb img total in a loop
%       AUTO : automatique mode of analysis (skip file when user input
%       is needed, AUTO variable should be defined in the main file)
%       datafolder : path to the raw data, should be RawdataFolder in main
%       resfolder : path to where save the matlabdata, should be
%       MatfileFolder in main


set(0,'DefaultFigureWindowStyle','docked') % defining default figure window for plots

warning('off','all') % removing wrning (orange text) display in consol

%% Extracting parameters from tableExperimentalConditions
ExperimentalConditions = getExperimentalConditions(tableExperimentalConditions, date, manip);
if strcmp(ExperimentalConditions.experimentType, "DEFAULT")
    cprintf('red',['\nWARNING : Experimental conditions not found in ' 'D:\Matlab Analysis\ActinCortexAnalysis\ExperimentalData' filesep 'ExperimentalConditions.csv' ' \n' 'Working now with default parameters, which are probably not accurate ! \n'])
end

OPTICALINDEXCORRECTION = ExperimentalConditions.opticalIndexCorrection;
MultiBeadSizes = false;
BEADTYPE = ExperimentalConditions.beadType;
DIAMETER = ExperimentalConditions.beadDiameter;
if (contains(BEADTYPE, '_'))
    MultiBeadSizes = true
    diameters = split(DIAMETER,"_");
    diameters = str2double(diameters);
    iL = 1 + (max(diameters) == diameters(2));
    iS = 3 - iL;
    DIAMETERL = diameters(iL);
    DIAMETERS = diameters(iS);
    types = split(BEADTYPE,"_");
    BEADTYPEL = types(iL);
    BEADTYPES = types(iS);

end

DELTA = ExperimentalConditions.beadsBrightSpotDelta;
% WITHFLUO = boolean(ExperimentalConditions.withFluoImages);

%%
% Displaying automode status
if AUTO
    cprintf('red','\n\n[AUTOMODE : ON]\n\n')
else
    cprintf('red','\n\n[AUTOMODE : OFF]\n\n')
end

% data and saving folders
path =[datafolder filesep date(1:2) '.' date(4:5) '.' date(7:8)]; % path to folder with raw data
sf   = resfolder; % save folder

mkdir([sf filesep 'R2V']) % create saving folder

if wFluo
    saveFluoFolder = [datafolder filesep date(1:2) '.' date(4:5) '.' date(7:8) '_Fluo'];
    if ~exist(saveFluoFolder, 'dir')
        mkdir(saveFluoFolder)
    end
end

savename=['R2V_' date '_' manip '_' tag '_' specif]; % saving name for matlab data file
savenamecst=['R2V_' date '_' manip '_' tag 'cst_' specif ];
savenamepart=[savename '_Part'] ; % saving name for saving of partial data file


cprintf('Unt',['\n\n' date '_' manip '_' specif '_AUTO=' num2str(AUTO) '\n\n'])

% retrieving list of all possible file names
[ListDep, ListFor] = get_fich_name(fich);



% Checking for previously saved data
E1 = exist([sf filesep 'R2V' filesep savenamepart '.mat'],'file');
E2 = exist([sf filesep 'R2V' filesep savename '.mat'],'file');

if E1 == 2 || E2 == 2
    % Loading of previously saved data, and marking of said data as already
    % analyzed to prevent reanalysis
    
    fprintf('\nLoading of previously saved data...')
    
    if E1 == 2
        load([sf filesep 'R2V' filesep savenamepart])
    elseif E2 == 2
        load([sf filesep 'R2V' filesep savename])
    end
    
    %%% identifying already analyzed data and making a list
    count = 1; % counter
    for i = 1:length(ListD)
        if not(isempty(ListD{i}))
            newlistD{count} =  ListD{i};
            count = count +1;
        end
    end
    ListD = newlistD;
    clear c newlistD
    kii = length(MT);
    [ListDep,pDep] = setdiff(ListDep,ListD);
    ListFor = {ListFor{pDep}};
    %%% end of identification
    
    cprintf('Com', [' OK.']);
    fprintf(['\n' num2str(kii) ' analyzed videos loaded.\n\n\n']);
    
else
    % if no previously analyzed data, list starts empty
    kii = 0;
    ListD = {};
    ListF = {};
end


nacq=length(ListFor); % number of possible names for acquired videos


for ki=1:nacq
    
    
    % defining specific tags for videos and results file names
    Noim=ListDep{ki}; Nofich=ListFor{ki};
    
    % creating names for loading
    stackname = [date '_' manip '_' Nofich '_' tag '.tif'];
    name=[date '_' manip '_' Noim '_' tag];
    resname = [name '_' restype '.txt'];
    
    % checking existence
    E = exist([path filesep resname],'file');
    
    if E == 2
        if AUTO
            % automode with try and catch skip videos when user input is
            % needed
            
            try
                
                [Srmp, X1rmp, X2rmp, Y1rmp, Y2rmp, dzrmp, idxRamp, Scst, X1cst, X2cst, Y1cst, Y2cst, dzcst, Sfluo, nBlackImagesOfLoopEnd] = doAnalysis(date, manip, Noim, path, resname, tag, restype, stackname,...
                                                                                                                                                        name, nimg, nimgbr, nimgr, wFluo,depthonameSB,depthonameLB, datafolder, ...
                                                                                                                                                        OPTICALINDEXCORRECTION, DIAMETERS, DIAMETERL, DELTA, AUTO);
                
                kii = kii +1; % compteur de ficheir trouv�
                
                fieldname = [name '_Field.txt'];
                
                fprintf(['\nLoading of ' fieldname '...']);
                BTMat = dlmread([path filesep fieldname],'\t',0,0);
                cprintf('Com', ' OK\n\n')
                
                BTMat = correctionFieldData(BTMat, nimg, nimgbr, nimgr, nBlackImagesOfLoopEnd); 
                % If N black images detected in given loop, this removes N data points at the end of the ramp and add N
                % "mock data points" (zeros) at the black image location in the stack (they won't be analysed).
                
                % registering analyzed data in global matrix
                MT{kii}.exp = name;
                MT{kii}.stack = stackname;
                
                MT{kii}.Scst = Scst;
                MT{kii}.xcst = [X1cst X2cst];
                MT{kii}.ycst = [Y1cst Y2cst];
                MT{kii}.dzcst = dzcst;
                
                MT{kii}.Srmp = Srmp;
                MT{kii}.xrmp = [X1rmp X2rmp];
                MT{kii}.yrmp = [Y1rmp Y2rmp];
                MT{kii}.dzrmp = dzrmp;
                
                MT{kii}.isInRamp = isInRamp;
                MT{kii}.idxRamp = idxRamp;
                MT{kii}.BTMat = BTMat;
                
                % marking file as analyzed
                ListD{kii} = Noim;
                ListF{kii} = Nofich;
                
                % partial saving
                fprintf('\nPartial saving...');
                save([sf filesep 'R2V' filesep savenamepart],'MT','ListD','ListF');
                cprintf('Com', ' OK\n\n')
                
                
            catch ERR
                cprintf('*err',['\n\nFile ' name ' has not been analyzed. \n'])
                cprintf('-err',['Cause : ' ERR.message '\n'])
                cprintf('err',['Error in ' ERR.stack(1).name '\n\n'])
            end
            
            
        else
            [Srmp, X1rmp, X2rmp, Y1rmp, Y2rmp, dzrmp, idxRamp, Scst, X1cst, X2cst, Y1cst, Y2cst, dzcst, Sfluo, nBlackImagesOfLoopEnd] = doAnalysis(date, manip, Noim, path, resname, tag, restype, stackname,...
                                                                                                                                                        name, nimg, nimgbr, nimgr, wFluo, depthonameSB,depthonameLB, datafolder, ...
                                                                                                                                                        OPTICALINDEXCORRECTION, DIAMETERS, DIAMETERL, DELTA, AUTO);
                                                                                                                                                    
            kii = kii +1; % compteur de ficheir trouv�
            
            fieldname = [name '_Field.txt'];
            
            fprintf(['\nLoading of ' fieldname '...']);
            BTMat = dlmread([path filesep fieldname],'\t',0,0);
            cprintf('Com', ' OK\n\n')
            
            BTMat = correctionFieldData(BTMat, nimg, nimgbr, nimgr, nBlackImagesOfLoopEnd);
            
            % registering analyzed data in global matrix
            MT{kii}.exp = name;
            MT{kii}.stack = stackname;
            
            MT{kii}.Scst = Scst;
            MT{kii}.xcst = [X1cst X2cst];
            MT{kii}.ycst = [Y1cst Y2cst];
            MT{kii}.dzcst = dzcst;
            
            MT{kii}.Srmp = Srmp;
            MT{kii}.xrmp = [X1rmp X2rmp];
            MT{kii}.yrmp = [Y1rmp Y2rmp];
            MT{kii}.dzrmp = dzrmp;
            
            MT{kii}.idxRamp = idxRamp;
            MT{kii}.BTMat = BTMat;
            
            % marking file as analyzed
            ListD{kii} = Noim;
            ListF{kii} = Nofich;
            
            % partial saving
            fprintf('Partial saving...');
            save([sf filesep 'R2V' filesep savenamepart],'MT','ListD','ListF');
            cprintf('Com', ' OK\n\n\n')
            
            % Saving fluo images in a separated folder if necessary
            if wFluo
                fprintf('Saving fluo images...');
                %saveFluoSubFolder = [datafolder filesep date(7:8) '.' date(4:5) '.' date(1:2) '_Fluo'];
                saveFluoSubFolder = [saveFluoFolder filesep name];
                stackpath = [path filesep stackname];
                if ~exist(saveFluoSubFolder, 'dir')
                   mkdir(saveFluoSubFolder)
                end
                N = length(Sfluo);
                for i = 1:N
                    currentImageFluo = imread(stackpath,'tiff',Sfluo(i));
                    imwrite(currentImageFluo, [saveFluoSubFolder filesep name '_Fluo_' num2str(Sfluo(i)) '.tif'])
                end
                fprintf(' OK\n\n\n')
            end
            
        end
        
    end
end


% saving of all analyzed data and removing partial matlab date file
% (if not in automode)
if exist('MT')
    if AUTO
        fprintf('\nPartial saving...(END OF AUTOMODE)\n\n');
        save([sf filesep 'R2V' filesep savenamepart],'MT','ListD','ListF');
    else
        fprintf('Complete saving...');
        save([sf filesep 'R2V' filesep savename],'MT','ListD','ListF');
        delete([sf filesep 'R2V' filesep savenamepart '.mat']);
        
    end
    cprintf('Com', ' OK\n\n\n')
end


% saving of constant parts of curve
if exist('MT')
    
    for kii = 1:length(MT)
        if not(isempty(MT{kii}))
            
            MT2{kii}.exp = MT{kii}.exp;
            MT2{kii}.stack = MT{kii}.stack;
            MT2{kii}.S = MT{kii}.Scst;
            MT2{kii}.x = MT{kii}.xcst;
            MT2{kii}.y = MT{kii}.ycst;
            MT2{kii}.dz = MT{kii}.dzcst;
            MT2{kii}.BTMat = MT{kii}.BTMat;

        end
    end
    
    MT = MT2;
    fprintf('Saving of constant parts...');
    save([sf filesep 'R2V' filesep savenamecst],'MT','ListD','ListF');
    fprintf(' OK\n\n\n')
    
else
    cprintf('err','No file found ! Check parameters.')
end



end






%%%%% Analysis function

function [Srmp, X1rmp, X2rmp, Y1rmp, Y2rmp, dzrmp, idxRamp, Scst, X1cst, X2cst, Y1cst, Y2cst, dzcst, Sfluo, nBlackImagesOfLoopEnd] = doAnalysis(date, manip, Noim, path, resname, tag, restype, stackname,...
                                                                                                                                                        name, nimg, nimgbr, nimgr, wFluo, depthonameSB,depthonameLB, datafolder, ...
                                                                                                                                                        OPTICALINDEXCORRECTION, DIAMETERS, DIAMETERL, DELTA, AUTO)

fprintf([restype ' ' date '_' manip '_' Noim '_' tag ' : \n\n']);


% reading imageJ file
fprintf('Reading ImageJ results file...');
Mat=dlmread([path filesep resname],'\t',1,0);
fid = fopen([path filesep resname]);

%%% identification of columns of interest
HeadAll = fgetl(fid);
fclose(fid);
HeadCell = textscan(HeadAll,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
HeadTable = [HeadCell{:}];

Nxm = find(strcmp(HeadTable,'XM')) +1 ;
Nym = find(strcmp(HeadTable,'YM')) +1 ;
Ns = find(strcmp(HeadTable,'Slice')) +1 ;
Nstdmean = find(strcmp(HeadTable,'StdDev')) +1 ;

if isempty(Nstdmean)
    
    Nstdmean = find(strcmp(HeadTable,'Mean')) +1 ;
    
end
%%% end of identification


cprintf('Com', ' OK\n');


% checking video for black images
fprintf('Checking video...');

S = Mat(:,Ns);
N_loop = ceil(max(S)/nimg);
Sfull = 1:N_loop*nimg;

% Are there black images at the end of each loop ? And if yes, how many ?
stackpath = [path filesep stackname];
nBlackImagesOfLoopEnd = zeros(N_loop, 1);
for i_loop = 1:N_loop
    i_image = i_loop*nimg;
    countBlackImages = 0;
    currentImage = imread(stackpath,'tiff',i_image);
    pixelSum = sum(sum(currentImage));
    while pixelSum == 0
        countBlackImages = countBlackImages + 1;
        i_image = i_image - 1;
        currentImage = imread(stackpath,'tiff',i_image);
        pixelSum = sum(sum(currentImage));
    end
    nBlackImagesOfLoopEnd(i_loop) = countBlackImages;
end

% % Is there always the same number of black images at the end of every loop ?
% if max(nBlackImagesOfLoopEnd) == min(nBlackImagesOfLoopEnd)
%     nBlackImagesEachLoop = nBlackImagesOfLoopEnd(1);
% else
%     "Cas merdique";
% end
    
%%% Creating list of images to separate image triplets for
% constant part, and identifying ramp part without image
% triplets
[Sdown, Smid, Sup, Sramp, Sfluo, idxRamp] = SplitZRampImg(Sfull,nimgbr,nimgr,nimg,wFluo,nBlackImagesOfLoopEnd);

L = min([length(Sdown) length(Smid) length(Sup)]);

Sdown = Sdown(1:L);
Smid  = Smid(1:L) ;
Sup   = Sup(1:L)  ;

Sall = [Sdown;Smid;Sup];
%%% end of creating list


cprintf('Com', ' OK\n');


% Spliting data from imageJ file per bead
fprintf('Creating beads 2D trajectories...');

% creating tag for analysis log file name
if strcmp(restype,'Results')
    logtag = '';
else
    logtag = restype;
end

distcheck = 16; % max distance (in px) between last know position and next point



% creating beads trajectories in 2D
[X1tot,Y1tot,S1tot,ptr1,X2tot,Y2tot,S2tot,ptr2] = splitBeads([Mat(:,Nxm) S], stackpath, [path filesep name], Mat(:,Nym), distcheck, AUTO, logtag);

% checking that each image has both beads tracked on it
if sum(S1tot == S2tot) == length(S1tot)
    Stot = S1tot;
    clear S1tot S2tot
else
    error('Problem in splitbeads, some images don''t have both beads')
end


%%% for each triplet, getting the image where both beads are at their best
% focus to keep 2D data from this image. This is done by getting images
% with the higher STD (aka most dynamics gray level, corresponding to best
% focus)
STD1 = nan(3,length(Sdown));
STD2 = STD1;

% down
X1d = zeros(size(Sdown));
X2d = X1d;
[~,pd1,p1d] = intersect(Sdown,Stot);
[~,pd2,p2d] = intersect(Sdown,Stot);
X1d(pd1) = X1tot(p1d);
STD1(1,pd1) = Mat(ptr1(p1d),Nstdmean);
X2d(pd2) = X2tot(p2d);
STD2(1,pd2) = Mat(ptr2(p2d),Nstdmean);

% mid
X1m = zeros(size(Smid));
X2m = X1m;
[~,pm1,p1m] = intersect(Smid,Stot);
[~,pm2,p2m] = intersect(Smid,Stot);
X1m(pm1) = X1tot(p1m);
STD1(2,pm1) = Mat(ptr1(p1m),Nstdmean);
X2m(pm2) = X2tot(p2m);
STD2(2,pm2) = Mat(ptr2(p2m),Nstdmean);

% up
X1u = zeros(size(Sup));
X2u = X1u;
[~,pu1,p1u] = intersect(Sup,Stot);
[~,pu2,p2u] = intersect(Sup,Stot);
X1u(pu1) = X1tot(p1u);
STD1(3,pu1) = Mat(ptr1(p1u),Nstdmean);
X2u(pu2) = X2tot(p2u);
STD2(3,pu2) = Mat(ptr2(p2u),Nstdmean);


m = ismember(STD1+STD2,max(STD1+STD2));


X1all = [X1d; X1m; X1u];
X2all = [X2d; X2m; X2u];

X1cst = X1all(m);
X2cst = X2all(m);
Scst = Sall(m);



% redefinition de p1 et p2 apres la recuperation des "bon" x
% pour les partie cst en force
p1cst = ismember(Stot,Scst);
p2cst = ismember(Stot,Scst);

Y1cst = Y1tot(p1cst);
Y2cst = Y2tot(p2cst);


% idem pour les partie avec rampe
p1rmp = ismember(Stot,Sramp);
p2rmp = ismember(Stot,Sramp);

X1rmp = X1tot(p1rmp);
X2rmp = X2tot(p2rmp);

Y1rmp = Y1tot(p1rmp);
Y2rmp = Y2tot(p2rmp);



STD1rmp = Mat(ptr1(p1rmp),Nstdmean);
STD2rmp = Mat(ptr2(p2rmp),Nstdmean);

Srmp = Stot(p1rmp);

ptrTot = ismember(Sfull,Stot);
idxRamp = idxRamp(ptrTot);
idxRamp = idxRamp(p1rmp);


% m2cst = ismember(Sfull,Scst);
% m2rmp = ismember(Sfull,Srmp);
% m2 = m2cst | m2rmp;
% isInRamp = isInRamp(m2);
% idxRamp = idxRamp(m2);

%%% end of getting best focus

cprintf('Com', ' OK\n');

% whosWho = true if the big bead corresponds to pos X1, Y1 and false otherwise. Thus, to ensure we have big bead in 1
% and small in 2, we'll switch between 1 & 2 if whosWho is false.
whosWho = AsymmetricBeadPair_findWhosWho(X1cst(0),Y1cst(0),X2cst(0),Y2cst(0),Scst(0),stackpath, DIAMETERL, DIAMETERS);
if whosWho == false
    X1tot_tmp = X1tot;
    Y1tot_tmp = Y1tot;
    X1cst_tmp = X1cst;
    Y1cst_tmp = Y1cst;
    X1rmp_tmp = X1rmp;
    Y1rmp_tmp = Y1rmp;
    
    X1tot = X2tot;
    Y1tot = Y2tot;
    X1cst = X2cst;
    Y1cst = Y2cst;
    X1rmp = X2rmp;
    Y1rmp = Y2rmp;
    
    X2tot = X1tot_tmp;
    Y2tot = Y1tot_tmp;
    X2cst = X1cst_tmp;
    Y2cst = Y1cst_tmp;
    X2rmp = X1rmp_tmp;
    Y2rmp = Y1rmp_tmp;
end

%%% Computing dz
fprintf('Computing z and dz...');

% dz for triplet in constant part
dzcst = DzCalc_Zcorr_multiZ_MultiBeadSize(X1tot,X1cst,Y1tot,Y1cst,X2tot,...
    X2cst,Y2tot,Y2cst,Stot,Scst,stackpath,Sdown,Smid,Sup,depthonameSB,depthonameLB,DIAMETERL,DIAMETERS,DELTA,datafolder);

% dz during ramps
dzrmp = DzCalc_Zcorr_MultiBeadSize(X1rmp,Y1rmp,X2rmp,Y2rmp,Srmp,stackpath,depthonameSB,depthonameLB,DIAMETERL,DIAMETERS,DELTA,datafolder);

% sorting dzs
dzmat = [Scst dzcst';Srmp dzrmp'];

dzmatsort = sortrows(dzmat);

dztot = smooth(dzmatsort(:,2),'rlowess');
stot = dzmatsort(:,1);

dzrmp = dztot(ismember(stot,Srmp));
dzcst = dztot(ismember(stot,Scst));


dzcst = dzcst*OPTICALINDEXCORRECTION; % correction for optical index changes between cell medium and objective immersion fluid (example : OPTICALINDEXCORRECTION = 1.33/1.52 for cells in water and objective in oil)
dzrmp = dzrmp*OPTICALINDEXCORRECTION; 

cprintf('Com', ' OK\n');
%%% end of computing dz

end


% Correct the file data to get good correspondance between images and data in case there are black images at the end of
% each loop
function BTMat = correctionFieldData(BTMat, nimg, nimgbr, nimgr, nBlackImagesOfLoopEnd)
    if max(nBlackImagesOfLoopEnd) > 0
        [nRows, nCols] = size(BTMat);
        N = floor(nRows/nimg);
        for i_loop = 1:N
            BTMat = [BTMat(1:i_loop*nimg,:) ; zeros(nBlackImagesOfLoopEnd(i_loop), nCols) ; BTMat(i_loop*nimg+1:end,:)];
            index_endOfRamp = (i_loop-1)*nimg+nimgbr+nimgr;
            BTMat(index_endOfRamp-nBlackImagesOfLoopEnd(i_loop)+1:index_endOfRamp,:) = [];
        end    
    end   
end


function whosWho = AsymmetricBeadPair_findWhosWho(X1,Y1,X2,Y2,S,stackpath,diameterL,diameterS)

H = ceil(diameterL/2);
W = ceil(diameterL/2);

ROI1 = imread(stackpath,'tiff',S,'pixelregion',{[ceil(Y1)-H ceil(Y1)+H],[ceil(X1)-W ceil(X1)+W]});
% figure(1)
% colormap(gray(65535));
% I1 = imadjust(ROI1);
% image(I1)

[counts1,binLocations1] = imhist(ROI1);
counts1bis = smooth(counts1,3);
max1 = max(counts1bis);
% figure(2)
% plot(binLocations1, counts1bis)

ROI2 = imread(stackpath,'tiff',S,'pixelregion',{[ceil(Y2)-H ceil(Y2)+H],[ceil(X2)-W ceil(X2)+W]});
% figure(3)
% colormap(gray(65535));
% I2 = imadjust(ROI2);
% image(I2)

[counts2,binLocations2] = imhist(ROI2);
counts2bis = smooth(counts2,3);
max2 = max(counts2bis);
% figure(4)
% plot(binLocations2, counts2bis)

whosWho = (max1 > max2);
