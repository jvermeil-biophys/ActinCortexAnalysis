function [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads_Fluo(X_ind,stackname,filename,Y,AUTO,varargin)
% X_ind doit contenir les positions X en pixels

nimg = max(X_ind(:,2));

if length(varargin) == 1
    tag = varargin{1};
else
    tag = '';
end

if exist([filename '_AnalysisLog' tag '.mat'], 'file')
    load([filename '_AnalysisLog' tag '.mat']);
else
    Slog = [];
    ActionLog = [];
end

% On va commencer par initialiser XY1 et XY2

img = 1; % Compteur pour la position dans le stack d'images

while not(exist('X1','var'))
    
    ptr_img = find(X_ind(:,2) == img); % Renvoie le nombre positions X associ�es � l'image 'img'
    
    if length(ptr_img) == 2 % Exactement deux billes sur l'image      
        % on cr�e les variables pour les premi�res billes
        X1 = X_ind(ptr_img(1),1);
        S1 = img;
        ptr1 = ptr_img(1);
        Y1 = Y(ptr_img(1));
        
        X2 = X_ind(ptr_img(2),1);
        S2 = img;
        ptr2 = ptr_img(2);
        Y2 = Y(ptr_img(2));
        
    elseif length(ptr_img) > 2 % Plus de deux billes sur l'images
        
        if ismember(img,Slog) % Est-ce que l'utilisateur a d�j� point� les billes dans un run pr�c�dent ?
            Xim = ActionLog(img).pos(:,1); % Si oui, on charge leur position
            Yim = ActionLog(img).pos(:,2);
        else
            if AUTO % Si plus de 2 billes, le mode AUTO ne convient pas, on renvoie une erreur.
                error('User input skipped in auto mode')
            end
            
            % Afficher l'image de la stack et demander de cliquer sur les
            % deux billes d'interet
            I = imread(stackname,'tiff',img);
            I = imadjust(I); % Ajuster luminosit� & contraste
            figure
            imshow(I)
            hold on
            plot(X_ind(ptr_img),Y(ptr_img),'o','color',[1 0.3 0],'markersize',5) % Afficher les centres des billes d�tect�es.
            title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(img) '/' num2str(nimg)])
            % Etape o� l'utilisateur s�lectionne les billes sur l'image
            btn = uicontrol('Style', 'pushbutton', 'String', 'No beads',...
                'Position', [20 20 50 20],'Callback', 'close');
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


% on parcours le reste des images
for i=img:nimg
    
    ptr_i = find(X_ind(:,2) == i); % on trouve 0, 1, 2 billes ou plus
    % sur chaque images
    
    if length(ptr_i) == 2 % deux bille sur l'image
        
        
        if ((abs(sqrt((X_ind(ptr_i(1),1)-X1(end))^2 + (Y(ptr_i(1))-Y1(end))^2)) > 16 ...
                && abs(sqrt((X_ind(ptr_i(2),1)-X1(end))^2 + (Y(ptr_i(2))-Y1(end))^2)) > 16)...
                && (abs(sqrt((X_ind(ptr_i(1),1)-X_ind(ptr_i(2),1))^2 + (Y(ptr_i(1))-Y(ptr_i(2)))^2)) < 100)...
                ) ...
                || ((abs(sqrt((X_ind(ptr_i(1),1)-X2(end))^2 + (Y(ptr_i(1))-Y2(end))^2)) > 16 ...
                && abs(sqrt((X_ind(ptr_i(2),1)-X2(end))^2 + (Y(ptr_i(2))-Y2(end))^2)) > 16)...
                && (abs(sqrt((X_ind(ptr_i(1),1)-X_ind(ptr_i(2),1))^2 + (Y(ptr_i(1))-Y(ptr_i(2)))^2)) < 100)...
                )
            
            
            
            % Si une des deux billes est loin des deux positions pr�c�dente on
            % demande l'information � l'utilisateur
            
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
                
                figure(1001)
                imshow(I)
                hold on
                plot(X1,Y1,'c--')
                plot(X2,Y2,'m--')
                plot(X1(end),Y1(end),'b*','markersize',5)
                plot(X2(end),Y2(end),'r*','markersize',5)
                plot(X_ind(ptr_i),Y(ptr_i),'o','color',[1 0.3 0],'markersize',5)
                title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(i) '/' num2str(nimg)])
                
                button = questdlg('Can you select the beads ?',...
                    'Selection or next','No');
                
                
                ActionLog(i).button = button;
                
            end
            
            if strcmp(button,'Yes')
                
                if ismember(i,Slog)
                    Xim = ActionLog(i).pos(:,1);
                    Yim = ActionLog(i).pos(:,2);
                else
                    
                    figure(1001)
                    % A nouveau un pointage des billes par l'utilisateur
                    % dans ce cas.
                    [Xim,Yim] = ginput(2);
                    close
                    
                    Slog = [Slog i];
                    ActionLog(i).pos = [Xim,Yim];
                    save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
                end
                
                % Attribution des deux points choisis par l'user � la
                % cat�gorie 'bille 1' ou 'billes 2'
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
            
% Q: A QUOI SERT CE TEST ??  A: C'est un cas ou on va rien faire et passer � l'image suivante     
        elseif (abs(sqrt((X_ind(ptr_i(1),1)-X1(end))^2 + (Y(ptr_i(1))-Y1(end))^2)) > 16 ...
                && abs(sqrt((X_ind(ptr_i(2),1)-X1(end))^2 + (Y(ptr_i(2))-Y1(end))^2)) > 16) ...
                || (abs(sqrt((X_ind(ptr_i(1),1)-X2(end))^2 + (Y(ptr_i(1))-Y2(end))^2)) > 16 ...
                && abs(sqrt((X_ind(ptr_i(2),1)-X2(end))^2 + (Y(ptr_i(2))-Y2(end))^2)) > 16)
            

% Les 2 cas simples ou les 2 billes ont peu boug�
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
        
        
    elseif length(ptr_i) >= 2 % deux ou plus bille sur l'image
        
        X1_tmp = [];
        X2_tmp = [];
        
        for ii = 1:length(ptr_i)
            
            
            if abs(sqrt((X_ind(ptr_i(ii),1)-X1(end))^2 + (Y(ptr_i(ii))-Y1(end))^2)) > 16 ...
                    && abs(sqrt((X_ind(ptr_i(ii),1)-X2(end))^2 + (Y(ptr_i(ii))-Y2(end))^2)) > 16
                % c'est une tache ou une autre bille
                
            elseif abs(sqrt((X_ind(ptr_i(ii),1)-X1(end))^2 + (Y(ptr_i(ii))-Y1(end))^2)) < ...
                    abs(sqrt((X_ind(ptr_i(ii),1)-X2(end))^2 + (Y(ptr_i(ii))-Y2(end))^2))
                
                
                X1_tmp = [X1_tmp X_ind(ptr_i(ii),1)];
                S1_tmp = i;
                ptr1_tmp = ptr_i(ii);
                Y1_tmp = Y(ptr_i(ii));
                
            else
                
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
                hold on
                plot(X1,Y1,'c--')
                plot(X2,Y2,'m--')
                plot(X1(end),Y1(end),'b*','markersize',5)
                plot(X2(end),Y2(end),'r*','markersize',5)
                plot(X_ind(ptr_i),Y(ptr_i),'o','color',[1 0.3 0],'markersize',5)
                title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(i) '/' num2str(nimg)])
                
                button = questdlg('Can you select the beads ?',...
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


