function [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads(X_ind,stackname,filename,Y,distcheck,AUTO,varargin)

% [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads(X_ind,stackname,filename,Y,distcheck,AUTO,varargin)
%
% X_ind is a matrix corresponding to [X Slices]
% 
% Separates the list of 2D position (X, Y) into tracks for two beads. Takes
% user input when beads have disapeard, aretacts are detected on the image,
% or beads have moved to far from one image to the next. User input is
% taken graphically by opening the problematic image. When in automode
% (AUTO = 1) an error is return when user input is needed so that the R2V
% code skips this video.
% 



set(0, 'defaultfigureposition', get(0, 'Screensize'));

nimg = max(X_ind(:,2));

if length(varargin) == 1
    tag = varargin{1};
else
    tag = '';
end

if exist([filename '_AnalysisLog' tag '.mat'])
    load([filename '_AnalysisLog' tag '.mat']);
else
    Slog = [];
    ActionLog = [];
end

%% on commence par initialiser XY1 et XY2
img = 1;
while not(exist('X1','var'))
    ptr_img = find(X_ind(:,2) == img);
    if length(ptr_img) == 2 % deux billes sur l'image
        
        % on crée les variables pour les premières billes
        
        X1 = X_ind(ptr_img(1),1);
        S1 = img;
        ptr1 = ptr_img(1);
        Y1 = Y(ptr_img(1));
        
        X2 = X_ind(ptr_img(2),1);
        S2 = img;
        ptr2 = ptr_img(2);
        Y2 = Y(ptr_img(2));
        
    elseif length(ptr_img) > 2 % plus de deux billes sur l'images
        
        if ismember(img,Slog)
            Xim = ActionLog(img).pos(:,1);
            Yim = ActionLog(img).pos(:,2);
        else
            
            if AUTO
                error('User input skipped in auto mode')
            end
            
            % afficher l'image de la stack et demander de cliquer sur les
            % deux billes d'interet
            I = imread(stackname,'tiff',img);
            I = imadjust(I);
            figure(1000)
            imshow(I)
            set(gcf, 'Position', get(0, 'Screensize'));
            hold on            
            plot(X_ind(ptr_img),Y(ptr_img),'o','color',[0 0.3 0],'markersize',5)
            title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(img) '/' num2str(nimg)])
            btn = uicontrol('Style', 'pushbutton', 'String', 'No beads',...
                'Position', [20 20 50 20],...
                'Callback', 'close');
            [Xim,Yim] = ginput(2);
            
            
            Slog = [Slog img];
            ActionLog(img).pos = [Xim,Yim];
            
            save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
        end
        
        if exist('Xim','var')
            
            for ii = 1:length(ptr_img)
                
                if abs(sqrt((X_ind(ptr_img(ii),1)-Xim(1))^2 + (Y(ptr_img(ii))-Yim(1))^2)) > 16 ...
                        && abs(sqrt((X_ind(ptr_img(ii),1)-Xim(2))^2 + (Y(ptr_img(ii))-Yim(2))^2)) > 16
                    % c'est une tache ou une autre bille
                    
                elseif abs(sqrt((X_ind(ptr_img(ii),1)-Xim(1))^2 + (Y(ptr_img(ii))-Yim(1))^2)) ...
                        < abs(sqrt((X_ind(ptr_img(ii),1)-Xim(2))^2 + (Y(ptr_img(ii))-Yim(2))^2))
                    % c'est la bille 1
                    
                    X1 = X_ind(ptr_img(ii),1);
                    S1 = img;
                    ptr1 = ptr_img(ii);
                    Y1 = Y(ptr_img(ii));
                    
                else
                    %  c'est la bille 2
                    X2 = X_ind(ptr_img(ii),1);
                    S2 = img;
                    ptr2 = ptr_img(ii);
                    Y2 = Y(ptr_img(ii));
                end
            end
            close
        end
    end
    img = img + 1;
end


%% on parcours le reste des images
for i=img:nimg
    
    ptr_i = find(X_ind(:,2) == i); % on trouve 0, 1, 2 billes ou plus
    % sur chaque images
    
    if length(ptr_i) == 2 % deux bille sur l'image
                  
        
        
        if ((abs(sqrt((X_ind(ptr_i(1),1)-X1(end))^2 + (Y(ptr_i(1))-Y1(end))^2)) > distcheck ...  % distance entre premier point et ancienne bille 1
                && abs(sqrt((X_ind(ptr_i(2),1)-X1(end))^2 + (Y(ptr_i(2))-Y1(end))^2)) > distcheck)... % distance entre deuxieme point et ancienne bille 1. % distance entre les deux  nouveau point
                ) ...
                || ((abs(sqrt((X_ind(ptr_i(1),1)-X2(end))^2 + (Y(ptr_i(1))-Y2(end))^2)) > distcheck ...
                && abs(sqrt((X_ind(ptr_i(2),1)-X2(end))^2 + (Y(ptr_i(2))-Y2(end))^2)) > distcheck)...
                )
            
            
            
            % Si une des deux billes est loin des deux positions précédente on
            % demande l'information à l'utilisateur
            
            if ismember(i,Slog)
                button = ActionLog(i).button;
            else
                
                if AUTO
                    error('User input skipped in auto mode')
                end
                
                
                % afficher l'image de la stack et demander de cliquer sur les
                % deux billes d'interet
                I = imread(stackname,'tiff',i);
                I = imadjust(I);
                
                figure(1000)
                imshow(I)
                set(gcf, 'Position', get(0, 'Screensize'));
                hold on
                plot(X1,Y1,'c--')
                plot(X2,Y2,'m--')
                plot(X1(end),Y1(end),'b*','markersize',5)
                plot(X2(end),Y2(end),'r*','markersize',5)
                plot(X_ind(ptr_i),Y(ptr_i),'o','color',[0 0.3 0],'markersize',5)
                title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(i) '/' num2str(nimg)])
                
              
                button = MFquestdlg([0.75 0.1],'Can you select the beads ?',...
                    'Selection or next','No');
                
                ActionLog(i).button = button;
                
            end
            
            if strcmp(button,'Yes')
                
                if ismember(i,Slog)
                    Xim = ActionLog(i).pos(:,1);
                    Yim = ActionLog(i).pos(:,2);
                else
                    
                    figure(1000)
                    [Xim,Yim] = ginput(2);
                    close
                    
                    Slog = [Slog i];
                    ActionLog(i).pos = [Xim,Yim];
                    save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
                end
                
                a = repmat(X_ind(ptr_i,1),1,2);
                b = repmat(Xim',length(ptr_i),1);
                
                c = repmat(Y(ptr_i),1,2);
                d = repmat(Yim',length(ptr_i),1);
                
                dist = sqrt((a-b).^2 + (c-d).^2);
                
                [~,Ind1] = min(dist(:,1));
                [~,Ind2] = min(dist(:,2));
                
                
                if abs(sqrt((X_ind(ptr_i(Ind1),1)-X1(end))^2 + (Y(ptr_i(Ind1))-Y1(end))^2)) ...
                        < abs(sqrt((X_ind(ptr_i(Ind1),1)-X2(end))^2 + (Y(ptr_i(Ind1))-Y2(end))^2))
                    % c'est la bille 1 la plus proche de X1
                    
                    X1 = [X1 X_ind(ptr_i(Ind1),1)];
                    S1 = [S1 i];
                    ptr1 = [ptr1 ptr_i(Ind1)];
                    Y1 = [Y1 Y(ptr_i(Ind1))];
                    
                    X2 = [X2 X_ind(ptr_i(Ind2),1)];
                    S2 = [S2 i];
                    ptr2 = [ptr2 ptr_i(Ind2)];
                    Y2 = [Y2 Y(ptr_i(Ind2))];
                    
                else
                    % c'est la bille 2
                    X1 = [X1 X_ind(ptr_i(Ind2),1)];
                    S1 = [S1 i];
                    ptr1 = [ptr1 ptr_i(Ind2)];
                    Y1 = [Y1 Y(ptr_i(Ind2))];
                    
                    X2 = [X2 X_ind(ptr_i(Ind1),1)];
                    S2 = [S2 i];
                    ptr2 = [ptr2 ptr_i(Ind1)];
                    Y2 = [Y2 Y(ptr_i(Ind1))];
                end
                
            elseif ismember(i,Slog)
                close
            else
                close
                Slog = [Slog i];
                save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
            end
            
            

            
        elseif abs(sqrt((X_ind(ptr_i(1),1)-X1(end))^2 + (Y(ptr_i(1))-Y1(end))^2)) < ...
                abs(sqrt((X_ind(ptr_i(1),1)-X2(end))^2 + (Y(ptr_i(1))-Y2(end))^2))
            % la bille 1 est plus proche de X1
            
            X1 = [X1 X_ind(ptr_i(1),1)];
            S1 = [S1 i];
            ptr1 = [ptr1 ptr_i(1)];
            Y1 = [Y1 Y(ptr_i(1))];
            
            X2 = [X2 X_ind(ptr_i(2),1)];
            S2 = [S2 i];
            ptr2 = [ptr2 ptr_i(2)];
            Y2 = [Y2 Y(ptr_i(2))];
            
        else
            % la bille 1 est plus proche de X2
            
            X2 = [X2 X_ind(ptr_i(1),1)];
            S2 = [S2 i];
            ptr2 = [ptr2 ptr_i(1)];
            Y2 = [Y2 Y(ptr_i(1))];
            
            X1 = [X1 X_ind(ptr_i(2),1)];
            S1 = [S1 i];
            ptr1 = [ptr1 ptr_i(2)];
            Y1 = [Y1 Y(ptr_i(2))];
            
        end
        
        
    elseif length(ptr_i) > 2 %  plus de 2 billes sur l'image
        
        X1_tmp = [];
        X2_tmp = [];
        
        for ii = 1:length(ptr_i) % pour chacune des billes
            
            
            if abs(sqrt((X_ind(ptr_i(ii),1)-X1(end))^2 + (Y(ptr_i(ii))-Y1(end))^2)) > distcheck ... % distance entre la bille et ancienne pos 1
                    && abs(sqrt((X_ind(ptr_i(ii),1)-X2(end))^2 + (Y(ptr_i(ii))-Y2(end))^2)) > distcheck % distance entre la bille et ancienne pos 2
                % c'est trop loin d'une ancienne position
                
            elseif abs(sqrt((X_ind(ptr_i(ii),1)-X1(end))^2 + (Y(ptr_i(ii))-Y1(end))^2)) < distcheck 
                   
                % la bille est a moins de ~280 nm de l'ancienne position 1
                
                X1_tmp = [X1_tmp X_ind(ptr_i(ii),1)];
                S1_tmp = i;
                ptr1_tmp = ptr_i(ii);
                Y1_tmp = Y(ptr_i(ii));
                
            elseif  abs(sqrt((X_ind(ptr_i(ii),1)-X2(end))^2 + (Y(ptr_i(ii))-Y2(end))^2)) < distcheck
                
                
                % la bille est a moins de ~280 nm de l'ancienne position 2
                
                X2_tmp = [X2_tmp X_ind(ptr_i(ii),1)];
                S2_tmp = i;
                ptr2_tmp = ptr_i(ii);
                Y2_tmp = Y(ptr_i(ii));
            end
        end
        
        
        if (length(X1_tmp) == 1) && (length(X2_tmp) == 1) && not(X1_tmp == X2_tmp)
            X1 = [X1 X1_tmp];
            Y1 = [Y1 Y1_tmp];
            S1 = [S1 S1_tmp];
            ptr1 = [ptr1 ptr1_tmp];
            
            X2 = [X2 X2_tmp];
            Y2 = [Y2 Y2_tmp];
            S2 = [S2 S2_tmp];
            ptr2 = [ptr2 ptr2_tmp];
            
        else
            
            
            % afficher l'image de la stack et demander de cliquer sur les
            % deux billes d'interet
            
            
            if ismember(i,Slog)
                button = ActionLog(i).button;
            else
                
                if AUTO
                    error('User input skipped in auto mode')
                end
                
                
                % afficher l'image de la stack et demander de cliquer sur les
                % deux billes d'interet
                I = imread(stackname,'tiff',i);
                I = imadjust(I);
                
                figure(1000)
                imshow(I)
                set(gcf, 'Position', get(0, 'Screensize'));  
                hold on              
                plot(X1,Y1,'c--')
                plot(X2,Y2,'m--')
                plot(X1(end),Y1(end),'b*','markersize',5)
                plot(X2(end),Y2(end),'r*','markersize',5)
                plot(X_ind(ptr_i),Y(ptr_i),'o','color',[0 0.3 0],'markersize',5)
                title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(i) '/' num2str(nimg)])
                
                button = MFquestdlg([0.75 0.1],'Can you select the beads ?',...
                    'Selection or next','No');                
                
                ActionLog(i).button = button;
                
            end
            
            if strcmp(button,'Yes')
                
                if ismember(i,Slog)
                    Xim = ActionLog(i).pos(:,1);
                    Yim = ActionLog(i).pos(:,2);
                else
                    
                    figure(1000)
                    [Xim,Yim] = ginput(2);
                    close
                    
                    Slog = [Slog i];
                    ActionLog(i).pos = [Xim,Yim];
                    save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
                end
                
                a = repmat(X_ind(ptr_i,1),1,2);
                b = repmat(Xim',length(ptr_i),1);
                
                c = repmat(Y(ptr_i),1,2);
                d = repmat(Yim',length(ptr_i),1);
                
                dist = sqrt((a-b).^2 + (c-d).^2);
                
                [~,Ind1] = min(dist(:,1));
                [~,Ind2] = min(dist(:,2));
                
                
                if abs(sqrt((X_ind(ptr_i(Ind1),1)-X1(end))^2 + (Y(ptr_i(Ind1))-Y1(end))^2)) ...
                        < abs(sqrt((X_ind(ptr_i(Ind1),1)-X2(end))^2 + (Y(ptr_i(Ind1))-Y2(end))^2))
                    % c'est la bille 1 la plus proche de X1
                    
                    X1 = [X1 X_ind(ptr_i(Ind1),1)];
                    S1 = [S1 i];
                    ptr1 = [ptr1 ptr_i(Ind1)];
                    Y1 = [Y1 Y(ptr_i(Ind1))];
                    
                    X2 = [X2 X_ind(ptr_i(Ind2),1)];
                    S2 = [S2 i];
                    ptr2 = [ptr2 ptr_i(Ind2)];
                    Y2 = [Y2 Y(ptr_i(Ind2))];
                    
                else
                    % c'est la bille 2
                    X1 = [X1 X_ind(ptr_i(Ind2),1)];
                    S1 = [S1 i];
                    ptr1 = [ptr1 ptr_i(Ind2)];
                    Y1 = [Y1 Y(ptr_i(Ind2))];
                    
                    X2 = [X2 X_ind(ptr_i(Ind1),1)];
                    S2 = [S2 i];
                    ptr2 = [ptr2 ptr_i(Ind1)];
                    Y2 = [Y2 Y(ptr_i(Ind1))];
                end
                
            elseif ismember(i,Slog)
                close
            else
                close
                Slog = [Slog i];
                save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
            end
            
            
            
            
            
        end
        
    end
    
end

save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');


X1 = X1';
X2 = X2';
Y1 = Y1';
Y2 = Y2';
S1 = S1';
S2 = S2';
end


