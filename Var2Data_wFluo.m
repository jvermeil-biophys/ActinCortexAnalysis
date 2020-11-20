function NotSaved = Var2Data_wFluo(date,tag,Bcorrec,manip,specif,Db,resfolder)

%
%       Var2Data is taking the .mat file (from Res2Var) containing the trajectories of each
%       beads, and computes the distance between those beads centers, and the
%       thickness of the pinched object by substractiong the diameter of the
%       beads. It also gets the info in the _Field.txt file (previously saved by
%       Res2Var with the rest of the data) to create time and magnetic field
%       variables. It then computes the attraction force between the two beads
%       based on the magnetic field strength and the beads distance. This
%       version (_wFluo) is used for the analysis of constant field
%       experiments, with possible fluo images taken in addition to BF.
%       This Var2Data also removes some curves from further analysis if they don't
%       meet goodness criterions (see code for details).
%
%       NotSaved = Var2Data_wFluo(date,chmpMag,Bcorrec,manip,specif,Db,...
%         datafolder,resfolder)
% 
%       OUTPUT
%       NotSaved : Variable with the name of curves that were discarded
%
%       INPUTS
%       date : date of experiment 'dd-mm-yy'
%       tag : tag present on the tif and txt files (ex: '5mT', '5mT_Dicty')
%       Bcorrec : correction factor for the magnetic field (based on
%       calibration with teslameter)
%       manip : 'manip' number (ex: 'M1' ou 'M2' etc...),
%       specif : experimental condiion (ex: 'Dicty_WT')
%       Db : beads diameter in nm
%       datafolder : path to the raw data, should be RawdataFolder in main
%       resfolder : path to where save the matlabdata, should be
%       MatfileFolder in main
%

warning('off') % removing warning (orange text) display in consol

NotSaved = {}; % Initializing a variable to store names of discarded data

% default options for figures
set(groot,'DefaultFigureWindowStyle','docked')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot,'DefaultAxesFontSize',30)

% folder with initial data and were to save analysed data
path = [resfolder filesep 'R2V'];

datenow = datestr(now,'yy-mm-dd'); % today

mkdir([resfolder filesep 'V2D']) % create directory

savename=['V2D_' date '_' manip '_' specif]; % saving name

savenamepart=[savename '_Part'] ; % partial saving name

loadname=['R2V_' date '_' manip '_' tag '_' specif '.mat']; % inital data

% conversion factor from pixel to microns
pixelconv=1/15.8; % 100X
% pixelconv=1/9.7; % 63X

if exist([path filesep loadname],'file')
    
    fprintf(['\n\nLoading file ' loadname '\n'])
    load([path filesep loadname]);
    nc = length(MT); % numbers of curves
    
    for kc = 1:nc
        if not(isempty(MT{kc}))
            
            Valid = 1; % default status for saving is OK
            
            name = MT{kc}.exp;
            
            if strcmp(name(end-1:end),'in')
                name = name(1:end-2);
            end         
          
            
            fprintf(['\n\nExperiment : ' name '\n\n']);
            fprintf('Retrieving beads data...');
            
            S = MT{kc}.S; % list of image with beads tracked on it
            S = unique(S); % shouldn't be necessary, didn't had time to check
   
            % X data for each beads
            X1 = MT{kc}.x(:,1);
            X2 = MT{kc}.x(:,2);
            dxp = X1 - X2; % in pixels          
       
            % Y data for each beads
            Y1 = MT{kc}.y(:,1);
            Y2 = MT{kc}.y(:,2);
            dyp = Y1 - Y2; % in pixels   
    
            % Z data for both beads
            dz = MT{kc}.dz; % in nm
            
            cprintf('Com', [' OK.\n']);
            
            
            fprintf('Creation of vectors for time and magnetic field...');
            

            BTMat = MT{kc}.BTMat; % retrieving data from _Field.txt file
  
            B = Bcorrec*BTMat(S,1); % adjusted magnetic field
            T = (BTMat(S,2)-BTMat(1,2))/1000; % time starting at 0 in seconds

            
            ImgSep = median(round(diff(T)/0.1)*0.1); % approximation of time between two images
            
            NptsTot = round((T(end) - T(1))/ImgSep); % estimate of total number of images given duration of curve

            cprintf('Com', [' OK.\n']);
            
            
            fprintf('Beads distance computation...');
            
            % Conversion from pixels to nm
            dx = dxp*pixelconv*1000; % along the field
            dy = dyp*pixelconv*1000; % perpendicular to the field in the observation plan
                
            % distance computation
            D3tot = sqrt(dx.^2+dy.^2+dz.^2); % center to center distance
            D3 = D3tot - Db; % surface to surface distance (object thickness)
            Dyz = sqrt(dy.^2+dz.^2); % distance in the plan perpendicular to the field
            D2 = sqrt(dx.^2+dy.^2); % distance in the observation plan
            
            
            % Removal of problematic individual points
            % ptrdel = marked for deletion
            
            ptrdel = find(isnan(dz)); % problem with computed Z
            
            % problem with magnetic field
            if ~contains(specif,'Sep') && ~contains(tag,'R') % not used if experiment is observation of beads separation when field is turned off, or if constant data from compression experiment
                
                mB = find(B < 0.8*Bcorrec*str2double(tag(1:end-2))); % if magnetic field is lower than 80% of expected value (indication of coils or BOP problems during experiments)
                ptrdel = [ptrdel mB'];
                
            end
           
            
            % Isolated points         
          
            dd3 = diff(D3);            
            ptrdeld3 = find(abs(dd3)>300)+1; % points that are more than 300 nm from previous point
            if ~isempty(ptrdeld3)
               for kp = 1:length(ptrdeld3)
                  if ptrdeld3(kp)>length(dd3) || (abs(dd3(ptrdeld3(kp)))>300 &&... % if last point of curve or point also more than 300 nm from the next point 
                          sign(dd3(ptrdeld3(kp)))~=sign(dd3(ptrdeld3(kp)-1))) % in the other direction (meaning it's isolated)
                      ptrdel = [ptrdel ptrdeld3(kp)];
                  end
               end
            end
            
           
                        
            
            % Big jumps in dz
          ddz = diff(dz); 
            ddzcount = 0; 
            
            mDz = []; % mask for bad dz
            
            dztmp = dz;
            
            while (max(abs(ddz))> 700)&&(ddzcount<50) % there are some big jumps in dz
            mDz = [mDz; find(ddz>700)+1]; % jumps up position
            mDz = [mDz; find(ddz<-700)+1]; % jumps down position
            
            dztmp(mDz) = dztmp(mDz-1);
            
            ddz=diff(dztmp); % recalculating ddz to see if more jumping points
            ddzcount = ddzcount +1;
            end
            
            ptrdel = [ptrdel; mDz];
            
            
            % deletion
            ptrdel = unique(ptrdel);
            D3(ptrdel) = [];
            D3tot(ptrdel) = [];
            Dyz(ptrdel) = [];
            T(ptrdel) = [];
            dz(ptrdel) = [];
            dy(ptrdel) = [];
            dx(ptrdel) = [];
            B(ptrdel) = [];
            S(ptrdel) = [];

            if ddzcount>49 % if there are too many successive points in a jump curve is considered invalid
                Valid = 0;           
                cprintf('err', ' Too much ddz problems\n');
                pause(0.5) % for visualition
            else

            cprintf('Com', [' OK.\n']);
            
            end
     
        
            fprintf('Computation of force...');
                  
            if contains(specif,'M270') % for M270 beads
                
                M = 0.74257*1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % magnetization [A.m^-1]

                
            else % for M450 beads
                
                M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % magnetization [A.m^-1]
                
            end
            
            V = 4/3*pi*(Db/2)^3; % bead volume [nm^3]
            
            m = M.*10^-9*V; % dipolar moment [A.nm^2]
           
            Bind = 2*10^5*m./(D3tot.^3); % induction field from one bead on the other

            Nvois = 1; % number of beads in contact
            
            Btot = B + Nvois*Bind;% corrected B for neighbours
            
            if contains(specif,'M270') % for M270 beads
                
                Mtot = 0.74257*1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1); % corrected magnetization

            else % for M450 beads
                
                Mtot = 1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1); % corrected magnetization
                
            end
            
            mtot = Mtot.*10^-9*V; % corrected dipolar moment
            
            anglefactor=abs(3*(dx./D3tot).^2-1); % angle between chain and field direction
            
            F = 3*10^5.*anglefactor.*mtot.^2./D3tot.^4; % force [pN]
            
            cprintf('Com', [' OK.\n']);
            
                                            
            fprintf('Data validation...')   
            
            if Valid % if still valid
            if  (length(T)/NptsTot > 0.5) && ... % check if too many points have been removed
                    ((T(end)-T(1) > 120)|| contains(tag,'R')||... % check if duration of curve is long enough, overlook duration check if inside/outside/ramp/separation experiment
                    contains(specif,'Sep')||...
                    contains(specif,'Inside')||contains(specif,'Outside')) 
                                cprintf('Com', [' OK.\n']);
                Valid = 1; % still valid
            else
                Valid = 0; % not valid    
                cprintf('err', ' Not Enough Points\n');
                pause(0.5)
                
            end
            
            else
                cprintf('err', ' Not saved\n');
            end
            

            
            if Valid
            fprintf('Variable creation...');
            
                MR{kc}.name = name;
                MR{kc}.stack = MT{kc}.stack;
                MR{kc}.F=F;
                MR{kc}.dx=dx;
                MR{kc}.D2=D2;
                MR{kc}.D3=D3;            
                MR{kc}.dz=dz;
                MR{kc}.dy=dy;
                MR{kc}.Dyz=Dyz;
                MR{kc}.B=B;
                MR{kc}.time = T;
                MR{kc}.S = S;
              

                cprintf('Com', [' OK.\n\n']);                             
                
                fprintf('Partial saving...');
                save([resfolder filesep 'V2D' filesep savenamepart],'MR');
               cprintf('Com', [' OK.\n']);
            else
                NotSaved = {NotSaved{:}, name}; % if not saved, store name
            end
            
        end
        
    end
    
else
    error(['Couldn''t find the file : ' [path filesep loadname]]) % problem with file name
end

EE = exist('MR');

if EE == 1
    fprintf('\n\nComplete saving...');
    save([resfolder filesep 'V2D' filesep savename],'MR');
    delete([resfolder filesep 'V2D' filesep savenamepart '.mat']); % removing partial save file
    cprintf('Com', [' OK.\n\n']);
    
    fprintf([num2str(kc) ' data set analyzed.\n\n\n\n']);
else
    
    fprintf('No existing file in the given list. \nCheck for a path error. \n\n');
    
end


if isempty(NotSaved)
    NotSaved = {'empty'};
end

NotSaved = NotSaved{:};





