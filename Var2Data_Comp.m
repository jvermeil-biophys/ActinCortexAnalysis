function Var2Data_Comp(date,tag,Bcorrec,manip,specif,nimgr...
    ,nimgtot,CompDur,Db,resfolder,figurefolder)

%
%       Var2Data is taking the .mat file (from Res2Var) containing the trajectories of each
%       beads, and computes the distance between those beads centers, and the
%       thickness of the pinched object by substractiong the diameter of the
%       beads. It also gets the info in the _Field.txt file (previously saved by
%       Res2Var with the rest of the data) to create time and magnetic field
%       variables. It then computes the attraction force between the two beads
%       based on the magnetic field strength and the beads distance. This
%       version (_Comp) is used to analyze classical compression
%       experiments where all compressions have the same length.
%       This Var2Data also removes some points from the curves but doesn't
%       remove whole curves like Var2Data_wFluo.
%
%       Var2Data_Comp(date,tag,Bcorrec,manip,specif,nimgr...
%          ,nimgtot,CompDur,Db,datafolder,resfolder,figurefolder)
%
%       date : date of experiment 'dd-mm-yy'
%       tag : tag present on the tif and txt files (ex: '5mT', '5mT_Dicty')
%       Bcorrec : correction factor for the magnetic field (based on
%       calibration with teslameter)
%       manip : 'manip' number (ex: 'M1' ou 'M2' etc...),
%       specif : experimental condiion (ex: 'Dicty_WT')
%       nimgr : number of images in a ramp
%       nimgtot : number of images in a loop
%       CompDur : Duration of compression as a string and in seconds (ex : '2s', '12s')
%       Db : beads diameter in nm
%       resfolder : path to where are the matlabdata, should be
%       MatfileFolder in main
%       figurefolder : path to figure saving folder, should be FigureFolder
%       in main
%

warning('off') % removing warning (orange text) display in consol

% default options for figures
set(groot,'DefaultFigureWindowStyle','docked')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot,'DefaultAxesFontSize',30)


% folder with initial data and were to save analysed data
path = [resfolder filesep 'R2V']; % data from previous program

datenow = datestr(now,'yy-mm-dd'); % today

subFigureFolder = [figurefolder filesep datenow filesep 'Meca' filesep specif]; % where to save figs
mkdir(subFigureFolder);

mkdir([resfolder filesep 'V2D']); % where to save analyzed data

savename=['V2D_' date '_' manip '_' tag '_' specif]; % name for saving

loadname=['R2V_' date '_' manip '_' tag '_' specif '.mat']; % name to load previous data


% conversion factor from pixel to microns
pixelconv=1/15.8; % 100X
% pixelconv=1/9.7; % 63X

if exist([path filesep loadname],'file')
    
    fprintf(['\n\nLoading file ' loadname '\n'])
    load([path filesep loadname]);
    nc = length(MT); % numbers of curves
    
    for kc = 1:nc
        
        if not(isempty(MT{kc}))
            VALID = 1;
            name = MT{kc}.exp;
            
            fprintf(['Cell : ' name '\n\n'])
            
            %% Constant part
            fprintf('Retrieving data for constant part... ');
            
            Scst = MT{kc}.Scst; % list of image with beads tracked on it
            
            % X data for each beads
            X1cst = MT{kc}.xcst(:,1);
            X2cst = MT{kc}.xcst(:,2); % in pixels
            
            % Y data for each beads
            Y1cst = MT{kc}.ycst(:,1);
            Y2cst = MT{kc}.ycst(:,2); % in pixels
            
            % Z data for both beads
            dzcst = MT{kc}.dzcst/1000; % in nm
            
            % conversion in nm
            xpcst = [X1cst,X2cst];
            ypcst = [Y1cst,Y2cst];
            xcst=xpcst*pixelconv;
            ycst=ypcst*pixelconv;
            
            % Computation of distance between beads
            Dcst = abs(xcst(:,2)-xcst(:,1));
            D2cst = ((xcst(:,2)-xcst(:,1))^2 + (ycst(:,2)-ycst(:,1))^2)^0.5;
            dycst = abs(ycst(:,1)-ycst(:,2));
            D2cst =(Dcst.^2+dycst.^2).^0.5;
            D3cst = (D2cst.^2+dzcst.^2).^0.5;
            
            cprintf('Com', [' OK.\n']);
            
            
            %% Compression part
            fprintf('Retrieving data for compressions');
            
            Srmp = MT{kc}.Srmp; % list of image with beads tracked on it
            
            % X data for each beads
            X1rmp = MT{kc}.xrmp(:,1);
            X2rmp = MT{kc}.xrmp(:,2); % in pixels
            
            % Y data for each beads
            Y1rmp = MT{kc}.yrmp(:,1);
            Y2rmp = MT{kc}.yrmp(:,2); % in pixels
            
            % Z data for both beads
            dzrmp = MT{kc}.dzrmp/1000; % in nm
            
            % conversion in nm
            xprmp = [X1rmp,X2rmp];
            yprmp = [Y1rmp,Y2rmp];
            xrmp=xprmp*pixelconv;
            yrmp=yprmp*pixelconv;
            
            % Computation of distance between beads
            Drmp = abs(xrmp(:,2)-xrmp(:,1));
            dyrmp = abs(yrmp(:,1)-yrmp(:,2));
            D2rmp =(Drmp.^2+dyrmp.^2).^0.5;
            D3rmp = (D2rmp.^2+dzrmp.^2).^0.5;
            
            cprintf('Com', [' OK.\n']);
            
            
            %% Time curve and force computation
            fprintf('Reconstruction of data in time...');
            
            % putting cs and ramp data together
            D3all = [D3cst; D3rmp];
            D2all = [D2cst; D2rmp];
            dzall = [dzcst; dzrmp];
            Dall = [Dcst; Drmp];
            Sall = [Scst; Srmp];
            
            % referencing by image number
            D3ind = [Sall D3all];
            D2ind = [Sall D2all];
            Dind = [Sall Dall];
            dzind = [Sall dzall];
            
            % sorting by image number
            D3sort = sortrows(D3ind);
            D2sort = sortrows(D2ind);
            Dsort = sortrows(Dind);
            dzsort = sortrows(dzind);
            
            % retrieving data sorted in time
            D3 = D3sort(:,2);
            D2 = D2sort(:,2);
            Dx = Dsort(:,2);
            dz = dzsort(:,2);
            S = D3sort(:,1);            
            
            cprintf('Com', [' OK.\n']);
            
            
            fprintf('Creation of time and magnetic field vectors...');
            
            BTMat = MT{kc}.BTMat; % loading data gathered from _Field.txt file
            
             Bfull = BTMat(:,1)*Bcorrec; % magnetic field for all images
            
            Tfull = (BTMat(:,2)-BTMat(1,2))/1000;% timepoints of all images
            
            % Here is a correction for an error in labview :
            % at the begining of each loop there is an extra signal from
            % the camera registered as an image which does not exist. Here we
            % remove it and compensate by adding the right time point at
            % the end of the loop
            Tfull2 = [Tfull(2:end); Tfull(end)+0.05];
            ncycle = round(Scst(end)/nimgtot);
            kcycle = [1:ncycle] * nimgtot;
            Tfull2(kcycle) = Tfull2(kcycle-1)+0.05;            
            Bfull2 = [Bfull(2:end); Bfull(end)];
            Bfull2(kcycle) = Bfull2(kcycle-1);
            
            T = Tfull2(S); % timepoints for images kept (1 point per triplet)
            B = Bfull2(S); % field for for images kept (1 point per triplet)
            
            ptrcst = ismember(S,Scst);
            ptrrmp = ismember(S,Srmp);
            
            Tcst = T(ptrcst); % time for cst part
            Trmp = T(ptrrmp); % time for ramp part

            
            Bcst = B(ptrcst); % for cst part
            Brmp = B(ptrrmp); % for ramp

            
            cprintf('Com', [' OK.\n']);
            
            fprintf('Creation of force vector...');
            
            
            
            % Dipolar force
            
            D3nm = D3*1000; % in nm
            Dxnm = Dx*1000; % in nm
            
            if contains(specif,'M270') % for M270 beads
                
                M = 0.74257*1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % magnetization [A.m^-1]
                
            else % for M450 beads
                
                M = 1.05*1600*(0.001991*B.^3+17.54*B.^2+153.4*B)./(B.^2+35.53*B+158.1); % magnetization [A.m^-1]
                
            end
            
            V = 4/3*pi*(Db/2)^3; % bead volume [nm^3]
            
            m = M.*10^-9*V; % dipolar moment [A.nm^2]
            
            Bind = 2*10^5*m./(D3nm.^3); % induction field from one bead on the other
            
            Nvois = 1; % number of beads in contact
            
            Btot = B + Nvois*Bind;% corrected B for neighbours
            
            if contains(specif,'M270')
                
                Mtot = 0.74257*1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
                
            else
                
                Mtot = 1.05*1600*(0.001991*Btot.^3+17.54*Btot.^2+153.4*Btot)./(Btot.^2+35.53*Btot+158.1);
                
            end
            
            mtot = Mtot.*10^-9*V; % corrected dipolar moment
            
            anglefactor=abs(3*(Dxnm./D3nm).^2-1); % angle between chain and field direction
            
            F = 3*10^5.*anglefactor.*mtot.^2./D3nm.^4; % force [pN]
            
            
            Fcst = F(ptrcst); % force during constant part
            Frmp = F(ptrrmp); % force during ramps
            
            
            cprintf('Com', [' OK.\n']);

            
            %% Ramp splitting
            fprintf('Splitting ramps...')
            
            FullRmpDat = [Srmp Trmp D3rmp Frmp Brmp dzrmp]; % all data relative to ramps
            
            RmpDat = cell(1,ncycle); % variable for storing individual ramps
            CstRmp = cell(1,ncycle); % variable for constant data relative to ramp
            
            cyclevect = 1:ncycle; % list of loops
            
            bgnrmp = (cyclevect-1)*nimgtot+(nimgtot-nimgr)/2+1; % begining of ramps
            bgncst = (cyclevect-1)*nimgtot+1; % begining of constant parts
            
            
            for k = cyclevect
                
                m = ismember(Srmp,bgnrmp(k):bgnrmp(k)+nimgr-1); % image indicices for compression n°k
                Stmp = Srmp(m);
                D3tmp = D3rmp(m);
                Ftmp = Frmp(m);
                Btmp = Brmp(m);
                Ttmp = Trmp(m);
                dztmp = dzrmp(m);
                
                % removal of points with a big dz jump
                ddz = diff(dztmp);
                while max(abs(ddz))> 700
                    ptrdel = find(abs(ddz)>700)+1;
                    Stmp(ptrdel) = [];
                    D3tmp(ptrdel) = [];
                    Ftmp(ptrdel) = [];
                    Btmp(ptrdel) = [];
                    Ttmp(ptrdel) = [];
                    dztmp(ptrdel) = [];
                    ddz=diff(dztmp);
                    length(ptrdel)
                end
                
                % getting positions of constant data points before and after
                % ramp n°k
                mcstsurround = ismember(Scst,[bgncst(k):bgnrmp(k)-1 bgnrmp(k)+nimgr:k*nimgtot]);
                
                % getting positions of constant data points after ramp n°k-1
                % and before ramp n°k
                if k == 1
                    mcstbefore = ismember(Scst,bgncst(k):bgnrmp(k)-1);
                else
                    mcstbefore = ismember(Scst,[bgnrmp(k-1)+nimgr:(k-1)*nimgtot bgncst(k):bgnrmp(k)-1]);
                end
                
                
                % saving data for ramp n°k
                if not(isempty(D3tmp))
                    
                    Cstmed = median(D3cst(mcstsurround));
                    
                    
                    if Cstmed < 0
                        
                        RmpDat{k} = 'NegCst'; % curve if in the negative value around the ramp
                        
                    else
                        
                        RmpDat{k} = [Stmp Ttmp D3tmp Ftmp Btmp dztmp];
                        
                    end
                    
                    CstRmp{k} = {D3cst(mcstsurround) D3cst(mcstbefore)}; % saving constant data relativ to ramp n°k
                    
                else
                    RmpDat{k} = 'NoPts'; % no points for this ramp
                    
                    CstRmp{k} = 'NoPts';
                end
                
                CompDuration{k} = CompDur; % duration of ramp
                
            end
            
            
            
            CstDat = [Scst Tcst D3cst Fcst Bcst]; % saving all constant data
            
            cprintf('Com', [' OK.\n']);
            
            % if all ramps have problem, mark the curve as invalid
            try
                RD = [RmpDat{:}];
                if ischar(RD)
                    VALID = 0;
                end
            catch
                RD = [];
                clear RD
            end
            
            
            
            if VALID % storing data in global matrix MR
                
                fprintf('Variable creation...');
                
                MR{kc}.name = name;
                MR{kc}.D3=D3;
                MR{kc}.D2=D2;
                MR{kc}.dz=dz;
                MR{kc}.B=B;
                MR{kc}.F=F;
                MR{kc}.time = T;
                MR{kc}.S = S;
                MR{kc}.RampData = {RmpDat, FullRmpDat,CompDuration, CstRmp};
                MR{kc}.CstData  = CstDat;
                
                cprintf('Com', [' OK.\n\n']);
                
        
                
            end
            
        else
            
            cprintf('err',['\nNo data for : '  name '\n'])
        end
        
        
    end
end


EE = exist('MR','var');

if EE == 1
    fprintf('Saving data...');
    save([resfolder filesep 'V2D' filesep savename],'MR');
    cprintf('Com', [' OK.\n\n']);
    
    fprintf([num2str(kc) ' datasets analyzed\n\n\n\n']);
else
    
    fprintf('No existing file in the given list. \nCheck for a path error. \n\n');
    
end