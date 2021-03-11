function Res2Var_wFluo_multiZ(date,tag,manip,specif,fich,resfile,nimg,depthoname,wFluo,AUTO,tableExperimentalConditions,...
    datafolder,resfolder)

%       Res2Var is using the imageJ particle analysis data to create
%       individual trajectories for the beads. It is using the txt files
%       created by ImageJ to get the X and Y positions the beads, it then
%       loads the images of the video and computes the DZ between the beads
%       using the depthograph. This particular version (_wFluo_multiZ)
%       is used for analyzing constant field experiments, with triplet
%       of image in Z for each time points. It can handle the presence of
%       fluo images at the end of each loop. This code is specific for
%       experiments made with 100X oil objectiv (see end of code where
%       optical indices correction is applied)
%
%       Res2Var_wFluo_multiZ(restype,date,tag,manip,specif,fich,
%       depthoname,nimgbr,nimgr,nimg,...
%       AUTO,datafolder,resfolder)
%
%       date : date of experiment 'dd-mm-yy'
%       tag : tag present on the tif and txt files (ex: '5mT', '5mT_Dicty')
%       manip : 'manip' number (ex: 'M1' ou 'M2' etc...),
%       specif : experimental condiion (ex: 'Dicty_WT')
%       fich : chambers to analyze for the condition (ex: 1 or [1:3])
%       resfile : extension of text file (ex: 'Results' or  'Inside')
%       nimg : nb img total in a loop, can be put to 3 as default if there
%       is no fluo images in the movies
%       depthoname : depthograph (ex:'19-05-20_DepthographM450')
%       AUTO : automatique mode of analysis (skip file when user input
%       is needed, AUTO variable should be defined in the main file)
%       datafolder : path to the raw data, should be RawdataFolder in main
%       resfolder : path to where save the matlabdata, should be
%       MatfileFolder in main

set(0,'DefaultFigureWindowStyle','docked') % defining default figure window for plots

warning('off','all') % removing warning (orange text) display in consol

%% Extracting parameters from tableExperimentalConditions
ExperimentalConditions = getExperimentalConditions(tableExperimentalConditions, date, manip);
if strcmp(ExperimentalConditions.experimentType, "DEFAULT")
    cprintf('red',['\nWARNING : Experimental conditions not found in ' 'D:\Matlab Analysis\ActinCortexAnalysis\ExperimentalData' filesep 'ExperimentalConditions.csv' ' \n' 'Working now with default parameters, which are probably not accurate ! \n'])
end

OPTICALINDEXCORRECTION = ExperimentalConditions.opticalIndexCorrection;
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
    c = 1; % counter
    for i = 1:length(ListD)
        if not(isempty(ListD{i}))
            newlistD{c} =  ListD{i};
            c = c +1;
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
    resname = [name '_' resfile '.txt'];
    
    % checking existence
    E = exist([path filesep resname],'file');
    
    
    if E == 2
        if AUTO
            % automode with try and catch skip videos when user input is
            % needed
            try
                
                % main analysis code (present at the end of present code file)
                [S, X1, X2, Y1, Y2, dz, STD1, STD2] = doAnalysis(date, manip, Noim, tag, path, resname, resfile, stackname, name, nimg, depthoname, datafolder, OPTICALINDEXCORRECTION, AUTO);
                
                kii = kii +1; % number of analyzed videos
                
                fieldname = [name '_Field.txt'];
                
                fprintf(['\nLoading of ' fieldname '...']);
                BTMat = dlmread([path filesep fieldname],'\t',0,0);
                cprintf('Com', ' OK\n\n')
                
                % registering analyzed data in global matrix
                MT{kii}.exp = name;
                MT{kii}.stack = stackname;
                MT{kii}.S = S;
                MT{kii}.x = [X1 X2];
                MT{kii}.y = [Y1 Y2];
                MT{kii}.dz = dz;
                MT{kii}.std = [STD1 STD2];
                MT{kii}.BTMat = BTMat;
                
                % marking file as analyzed
                ListD{kii} = Noim;
                ListF{kii} = Nofich;
                
                % partial saving
                fprintf('\nPartial saving...');
                save([sf filesep 'R2V' filesep savenamepart],'MT','ListD','ListF');
                cprintf('Com', ' OK\n\n')
                
            catch ERR % display of reason for skipping analysis
                cprintf('*err',['\n\nFile ' name ' has not been analyzed. \n'])
                cprintf('-err',['Cause : ' ERR.message '\n'])
                cprintf('err',['Error in ' ERR.stack(1).name '\n\n'])
            end
        else
            
            % main analysis code (present at the end of present code file)
            [S, X1, X2, Y1, Y2, dz, STD1, STD2] = doAnalysis(date, manip, Noim, tag, path, resname, resfile, stackname, name, nimg, depthoname, datafolder, OPTICALINDEXCORRECTION, AUTO);
            
            kii = kii +1; % number of analyzed videos
            
            fieldname = [name '_Field.txt'];
            
            fprintf(['\nLoading of ' fieldname '...']);
            BTMat = dlmread([path filesep fieldname],'\t',0,0);
            cprintf('Com', ' OK\n\n')
            
            % registering analyzed data in global matrix
            MT{kii}.exp = name;
            MT{kii}.stack = stackname;
            MT{kii}.S = S;
            MT{kii}.x = [X1 X2];
            MT{kii}.y = [Y1 Y2];
            MT{kii}.dz = dz;
            MT{kii}.std = [STD1 STD2];
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
                slice = nimg + 1;
                Sfluo = [];
                while slice <= max(S) + 1
                    Sfluo = [Sfluo slice];
                    slice = slice + nimg + 1;
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

end


%%%%% Analysis function

function [S, X1, X2, Y1, Y2, dz, STD1, STD2] = doAnalysis(date, manip, Noim, tag, path, resname, resfile, stackname, name, nimg, depthoname, datafolder, OPTICALINDEXCORRECTION, AUTO)
fprintf(['Results ' date '_' manip '_' Noim '_' tag ' : \n\n']);

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
Nstd = find(strcmp(HeadTable,'StdDev')) +1 ;
%%% end of identification

cprintf('Com',' OK\n');


% checking video for black images
fprintf('Checking video...');

S = Mat(:,Ns);

%%% Creating list of images to separate image triplets
Sfull = 1:max(S);

[Sdown,Smid,Sup] = SplitZimage(Sfull,nimg);

[Smd,ptrd,ptrm] = intersect(Sdown+1,Smid);
[Smdu,ptrmd,ptru] = intersect(Smd,Sup-1);

l = min([length(Sdown) length(Smid) length(Sup)]);

Sdown = Sdown(ptrd(ptrmd)); % list of images 1, for lowest image in Z
Smid  = Smdu ; % list of images 2, for middle image in Z
Sup   = Sup(ptru); % list of images 3, for upper image in Z

Sall = [Sdown;Smid;Sup]; % table of all images in lists
%%% end of creating list

stackpath = [path filesep stackname]; % complete path for loading tif stack

cols = [];

for is = 1:max(S)
    
    Itmp = imread(stackpath,'tiff',is); % loading image
    
    % checking if image is black
    IntTot = sum(sum(Itmp));    
    if IntTot == 0
        [~, col] = find(Sall == is);
        cols = [cols col];
    end
    
end

Sdown(cols) = [];
Smid(cols) = [];
Sup(cols) = [];


Sall = [Sdown;Smid;Sup];

cprintf('Com',' OK\n');



% Spliting data from imageJ file per bead
fprintf('Creating beads 2D trajectories...');

% creating tag for analysis log file name
if strcmp(resfile,'Results')
    logtag = '';
else
    logtag= resfile;
end

distcheck = 16; % max distance (in px) between last know position and next point

% creating beads trajectories in 2D
[X1tot,Y1tot,S1tot,ptr1,X2tot,Y2tot,S2tot,ptr2] = splitBeads([Mat(:,Nxm) S],...
    stackpath,[path filesep name],Mat(:,Nym),distcheck,AUTO,logtag);

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


X1d = zeros(size(Sdown));
X2d = X1d;
[~,pd1,p1d] = intersect(Sdown,Stot);
[~,pd2,p2d] = intersect(Sdown,Stot);
X1d(pd1) = X1tot(p1d);
STD1(1,pd1) = Mat(ptr1(p1d),Nstd);
X2d(pd2) = X2tot(p2d);
STD2(1,pd2) = Mat(ptr2(p2d),Nstd);

X1m = zeros(size(Smid));
X2m = X1m;
[~,pm1,p1m] = intersect(Smid,Stot);
[~,pm2,p2m] = intersect(Smid,Stot);
X1m(pm1) = X1tot(p1m);
STD1(2,pm1) = Mat(ptr1(p1m),Nstd);
X2m(pm2) = X2tot(p2m);
STD2(2,pm2) = Mat(ptr2(p2m),Nstd);

X1u = zeros(size(Sup));
X2u = X1u;
[~,pu1,p1u] = intersect(Sup,Stot);
[~,pu2,p2u] = intersect(Sup,Stot);
X1u(pu1) = X1tot(p1u);
STD1(3,pu1) = Mat(ptr1(p1u),Nstd);
X2u(pu2) = X2tot(p2u);
STD2(3,pu2) = Mat(ptr2(p2u),Nstd);

m = ismember(STD1+STD2,max(STD1+STD2));

X1all = [X1d; X1m; X1u];
X2all = [X2d; X2m; X2u];

X1 = X1all(m);
X2 = X2all(m);
S = Sall(m);

p1 = ismember(Stot,S);
p2 = ismember(Stot,S);

Y1 = Y1tot(p1);
Y2 = Y2tot(p2);
%%% end of getting best focus

cprintf('Com',' OK\n');

%%% Computing dz
fprintf('Computing z and dz...');

dz = DzCalc_Zcorr_multiZ(X1tot,X1,Y1tot,Y1,X2tot,...
    X2,Y2tot,Y2,Stot,S,stackpath,Sdown,Smid,Sup,depthoname,datafolder);

dz = smooth(dz,'rlowess');

dz = dz*OPTICALINDEXCORRECTION; % correction for optical index changes between medium and objective immersion liquid

cprintf('Com', ' OK\n');
%%% end of computing dz

end