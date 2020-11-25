function [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads(X_ind,stackname,filename,Y,distcheck,AUTO,varargin)

% [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads(X_ind,stackname,filename,Y,distcheck,AUTO,varargin)
%
% X_ind is a matrix corresponding to [X Slices] 
% X : list of X coordinates of detected objects
% Slices : list of images where the position objects were detected
% stackname : name of the .tif file
% Y : list of Y coordinates of detected objects
% distcheck : maximum distance that a bead can move from on image to
% the next without asking user input for verification
% AUTO = 1 : automode, AUTO = 0 no automode. Automode skips user input by
% returning an error. In R2V if automode is ON the error will skip the
% analysis of the current video through a try/catch
% varargin : possibility of defining a tag to add to the name of the log
% file
%
% Separates the list of 2D position (X, Y) into tracks for two beads. Takes
% user input when beads have disapeard, aretacts are detected on the image,
% or beads have moved to far from one image to the next. User input is
% taken graphically by opening the problematic image. When in automode
% (AUTO = 1) an error is return when user input is needed so that the R2V
% code skips this video.
%

set(0, 'defaultfigureposition', get(0, 'Screensize')); % figure option

nimg = max(X_ind(:,2)); % number of images to analyze

if length(varargin) == 1 % defining saving tag
    tag = varargin{1};
else
    tag = '';
end

if exist([filename '_AnalysisLog' tag '.mat']) % checking if user input has previously been saved
    load([filename '_AnalysisLog' tag '.mat']);
else
    Slog = []; % if not initializing user input log
    ActionLog = [];
end

%% Initialisation of X and Y variables
img = 1; % starting at image 1
while not(exist('X1','var')) % as long as variables have not been initialized,
    ptr_img = find(X_ind(:,2) == img);  % look for an image with at least two beads on it
    if length(ptr_img) == 2 % Two beads on the image
        
        % initialisation of X and Y variable with position of those two
        % beads
        
        % bead 1
        X1 = X_ind(ptr_img(1),1); % X pos
        S1 = img; % Image number in stack
        ptr1 = ptr_img(1); % Indice of current image in the data matrix from Results.txt
        Y1 = Y(ptr_img(1)); % Y pos
        
        % bead 2
        X2 = X_ind(ptr_img(2),1);
        S2 = img;
        ptr2 = ptr_img(2);
        Y2 = Y(ptr_img(2));
        
    elseif length(ptr_img) > 2 % more than two beads on the image => user input necessary
        
        if ismember(img,Slog) % Check if user input for this image is saved
            % saved position clicked by user
            Xim = ActionLog(img).pos(:,1);
            Yim = ActionLog(img).pos(:,2);
        else % if not ask for user input
            
            if AUTO % Skip if AUTOMODE is on
                error('User input skipped in auto mode')
            end
            
            % Display current image
            I = imread(stackname,'tiff',img);
            I = imadjust(I);
            figure(1000)
            imshow(I)
            set(gcf, 'Position', get(0, 'Screensize'));
            hold on
            plot(X_ind(ptr_img),Y(ptr_img),'o','color',[0 0.3 0],'markersize',5)
            title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(img) '/' num2str(nimg)])
            
            % Collect user input (clicking on the image)
            [Xim,Yim] = ginput(2);
            
            % Save user input for this image
            Slog = [Slog img];
            ActionLog(img).pos = [Xim,Yim];
            save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
            close
            
        end
       
          % Finding the two beads closest to the user clicked
                % positions
                
                Xmat = repmat(X_ind(ptr_img,1),1,2); % X of the objects on current image
                Ymat = repmat(Y(ptr_img),1,2); % Y of the objects on current image
                
                Xclickmat = repmat([Xim(1) Xim(2)],length(ptr_img),1); % X of the beads clicked by user
                Yclickmat = repmat([Yim(1) Yim(2)],length(ptr_img),1); % Y of the beads clicked by user
                
                Dist2mat = (Xmat-Xclickmat).^2+(Ymat-Yclickmat).^2; % Square distance matrix (cost matrix for the matching algorithm)
                
                A = HungarianMatch(Dist2mat); % Best match for each line
                
                % Assigning new beads positions
                X1 = X_ind(ptr_img(A == 1),1);
                Y1 = Y(ptr_img(A == 1),1);
                S1 = img;
                ptr1 = ptr_img(A == 1);
                
                X2 = X_ind(ptr_img(A == 2),1);
                Y2 = Y(ptr_img(A == 2),1);
                S2 = img;
                ptr2 = ptr_img(A == 2);  
    end
    img = img + 1;
end


%% Going through the rest of images
for i=img:nimg
    
    ptr_i = find(X_ind(:,2) == i); % number of object detected on current image
    
    if length(ptr_i) >= 2 % two or more objects
        
        % Finding the two objects closest to previous beads positions
        
        Xmat = repmat(X_ind(ptr_i,1),1,2); % X of the objects on current image
        Ymat = repmat(Y(ptr_i),1,2); % Y of the objects on current image
        
        Xposmat = repmat([X1(end) X2(end)],length(ptr_i),1); % X of the beads on previous image
        Yposmat = repmat([Y1(end) Y2(end)],length(ptr_i),1); % Y of the beads on previous image
        
        Dist2mat = (Xmat-Xposmat).^2+(Ymat-Yposmat).^2; % Square distance matrix (cost matrix for the matching algorithm)
        
        A = HungarianMatch(Dist2mat); % Best match for each line
        
        Dist1 = sqrt(Dist2mat(A == 1,1)); % distance of the closest position to previous bead 1
        Dist2 = sqrt(Dist2mat(A == 2,2)); % distance of the closest position to previous bead 2
        
        if Dist1 < distcheck && Dist2 < distcheck % if the two found objects are
            % close enough to the previous positions, they are validated as beads
            % of interest on this image
            
            X1(end+1) = X_ind(ptr_i(A == 1),1);
            Y1(end+1) = Y(ptr_i(A == 1),1);
            S1(end+1) = i;
            ptr1(end+1) = ptr_i(A == 1);
            
            X2(end+1) = X_ind(ptr_i(A == 2),1);
            Y2(end+1) = Y(ptr_i(A == 2),1);
            S2(end+1) = i;
            ptr2(end+1) = ptr_i(A == 2);
            
        else % if not user is asked to click near the correct objects
            
            if ismember(i,Slog) % check if user input for this image already exist
                button = ActionLog(i).button;
            else
                
                if AUTO % skip if AUTOMODE is ON
                    error('User input skipped in auto mode')
                end
                
                
                % Display image
                I = imread(stackname,'tiff',i);
                I = imadjust(I);
                
                figure(1000)
                imshow(I)
                set(gcf, 'Position', get(0, 'Screensize'));
                hold on
                plot(X1,Y1,'c--') % Bead 1 trajectory
                plot(X2,Y2,'m--') % Bead 2 trajectory
                plot(X1(end),Y1(end),'b*','markersize',5) % Bead 1 last position
                plot(X2(end),Y2(end),'r*','markersize',5) % Bead 2 last position
                plot(X_ind(ptr_i),Y(ptr_i),'o','color',[0 0.3 0],'markersize',5) % Objects detected on the image
                title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(i) '/' num2str(nimg)])
                
                % ask if the two beads of interest are detected on this
                % image
                button = MFquestdlg([0.75 0.1],'Can you select the beads ?',...
                    'Selection or next','No');
                % save answer
                ActionLog(i).button = button;
                
            end
            
            if strcmp(button,'Yes') % If beads of interest are present on the image
                
                if ismember(i,Slog) % check if user input for this image already exist
                    Xim = ActionLog(i).pos(:,1);
                    Yim = ActionLog(i).pos(:,2);
                else
                    
                    figure(1000)
                    % User clicks on the beads
                    [Xim,Yim] = ginput(2);
                    close
                    % Saving clicking positions
                    Slog = [Slog i];
                    ActionLog(i).pos = [Xim,Yim];
                    save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
                end
                
                % Finding the two beads closest to the user clicked
                % positions
                
                Xmat = repmat(X_ind(ptr_i,1),1,2); % X of the objects on current image
                Ymat = repmat(Y(ptr_i),1,2); % Y of the objects on current image
                
                Xclickmat = repmat([Xim(1) Xim(2)],length(ptr_i),1); % X of the beads clicked by user
                Yclickmat = repmat([Yim(1) Yim(2)],length(ptr_i),1); % Y of the beads clicked by user
                
                Dist2mat = (Xmat-Xclickmat).^2+(Ymat-Yclickmat).^2; % Square distance matrix (cost matrix for the matching algorithm)
                
                A = HungarianMatch(Dist2mat); % Best match for each line
                
                
                % These are the beads locations on the current image
                newX = [X_ind(ptr_i(A == 1),1); X_ind(ptr_i(A == 2),1)];
                newY = [Y(ptr_i(A == 1),1);Y(ptr_i(A == 2),1)];
                newptr = [ptr_i(A == 1);ptr_i(A == 2)];            
                
                
                % Assigning the new beads to the closest bead position on
                % the previous image
                
                newXmat = repmat(newX,1,2); % X of the objects on current image
                newYmat = repmat(newY,1,2); % Y of the objects on current image
                
                Xposmat = repmat([X1(end) X2(end)],2,1); % X of the beads on previous image
                Yposmat = repmat([Y1(end) Y2(end)],2,1); % Y of the beads on previous image
                
                newDist2mat = (newXmat-Xposmat).^2+(newYmat-Yposmat).^2; % Square distance matrix (cost matrix for the matching algorithm)
        
                newA = HungarianMatch(newDist2mat); % Best match for each line
                        
                
                % Assigning new beads to end of trajectories
                X1(end+1) = newX(newA == 1);
                Y1(end+1) = newY(newA == 1);
                S1(end+1) = i;
                ptr1(end+1) = newptr(newA == 1);
                
                X2(end+1) = newX(newA == 2);
                Y2(end+1) = newY(newA == 2);
                S2(end+1) = i;
                ptr2(end+1) = newptr(newA == 2);
                
            elseif ismember(i,Slog) % beads not present on the image and info already saved
                
                close
                
            else % beads not present on the image and info needs to be saved saved
                
                close
                
                Slog = [Slog i];
                save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
            end
        end
    end
end

% save final version of log
save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');

% transpose variable (because needed in R2V)
X1 = X1';
X2 = X2';
Y1 = Y1';
Y2 = Y2';
S1 = S1';
S2 = S2';
end

