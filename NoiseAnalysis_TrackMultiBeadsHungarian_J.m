function trackingTable = NoiseAnalysis_TrackMultiBeadsHungarian_J(X_ind,Y,stackname)
% X_ind doit contenir les positions X en pixels

nImages = max(X_ind(:,2));

% On va commencer par initialiser
initialImage = 1;
initialPointeurBeads = find(X_ind(:,2) == initialImage);
Nbeads = length(initialPointeurBeads);

while Nbeads == 0 || initialImage == 1
    initialPointeurBeads = find(X_ind(:,2) == initialImage); % Renvoie le nombre positions X associées à l'image 'currentImage'
    Nbeads = length(initialPointeurBeads);
    if Nbeads > 0
        XYPointeur_init = [X_ind(initialPointeurBeads, 1) Y(initialPointeurBeads) initialPointeurBeads];
        XYPointeur_init = sortrows(XYPointeur_init);
        for b = 1:Nbeads
            trackingTable(b).X = [XYPointeur_init(b,1)];
            trackingTable(b).Y = [XYPointeur_init(b,2)];
            trackingTable(b).ptr = [XYPointeur_init(b,3)];
            trackingTable(b).S = [initialImage];
        end
    initialImage = initialImage + 1;
    end
end


% on parcours le reste des images
for i=initialImage:nImages
    
    currentPointeurBeads = find(X_ind(:,2) == i);
    NbeadsCurrentImage = length(currentPointeurBeads);
    
    if NbeadsCurrentImage == Nbeads
        lastCoordinates = zeros(Nbeads, 2);
        for b = 1:Nbeads
            lastCoordinates(b, 1) = trackingTable(b).X(end);
            lastCoordinates(b, 2) = trackingTable(b).Y(end);
        end
        currentCoordinates = [X_ind(currentPointeurBeads, 1) Y(currentPointeurBeads)];
        matchingIndices = leastSquaresMatching(lastCoordinates', currentCoordinates');
        for b = 1:Nbeads
            try
                trackingTable(b).X = [trackingTable(b).X currentCoordinates(matchingIndices(b),1)];
                trackingTable(b).Y = [trackingTable(b).Y currentCoordinates(matchingIndices(b),2)];
                trackingTable(b).ptr = [trackingTable(b).ptr currentPointeurBeads(matchingIndices(b))];
                trackingTable(b).S = [trackingTable(b).S i];
            catch
                S;
            end
        end
    end
end

% Ici on plotte les trajectoires trackées sur la dernière image de la stack.
lastImage = imread(stackname,'tiff',nImages);
lastImage = imadjust(lastImage);
figure
imshow(lastImage)
hold on
colorTable = ['r', 'b', 'g', 'y', 'c', 'm'];
Ncolor = length(colorTable);
for b = 1:Nbeads
    plot(trackingTable(b).X,trackingTable(b).Y, [colorTable(mod(b,Ncolor)+1) '--'])
    plot(trackingTable(b).X(end),trackingTable(b).Y(end), [colorTable(mod(b,Ncolor)+1) 'o'],'markersize',5)
end
end

function matchingIndices = leastSquaresMatching(coordinates1, coordinates2)
% coordinates1 = [1 1; 1 2; 3 1; 8 9]';
% coordinates2 = [8 9; 3 1; 1 2]';

dist2 = @(x,y) dot(y-x,y-x);

% Construction de la matrice de coûts
% Le coût est ici simplement la distance d'une bille à l'autre, au carré.
costMat = zeros(size(coordinates1, 2), size(coordinates2, 2));
for i = 1:size(coordinates1, 2)
    for j = 1:size(coordinates2, 2)
        costMat(i,j) = dist2(coordinates1(:,i), coordinates2(:,j));
    end
end

% Algorithme hongrois - A = assign; C = total cost
% A contient donc un tableau qui permet de "matcher" les points de
% coordinates1 avec leur "successeur" de coordinates2.
% C contient le coût total du matching qui minimise le problème.
% HungarianMatch est une implémentation du hungarian matching trouvée sur
% le forum Mathworks (sous le nom original de munkres()).
[A,C] = HungarianMatch(costMat); 

matchingIndices = A;
end



% BOUTS DE CODE DE splitBeads() SUR LAQUELLE CETTE FONCTION EST BASEE.


%         if ((abs(sqrt((X_ind(pointeurCurrentImage(1),1)-X1(end))^2 + (Y(pointeurCurrentImage(1))-Y1(end))^2)) > 16 ...
%                 && abs(sqrt((X_ind(pointeurCurrentImage(2),1)-X1(end))^2 + (Y(pointeurCurrentImage(2))-Y1(end))^2)) > 16)...
%                 && (abs(sqrt((X_ind(pointeurCurrentImage(1),1)-X_ind(pointeurCurrentImage(2),1))^2 + (Y(pointeurCurrentImage(1))-Y(pointeurCurrentImage(2)))^2)) < 100)...
%                 ) ...
%                 || ((abs(sqrt((X_ind(pointeurCurrentImage(1),1)-X2(end))^2 + (Y(pointeurCurrentImage(1))-Y2(end))^2)) > 16 ...
%                 && abs(sqrt((X_ind(pointeurCurrentImage(2),1)-X2(end))^2 + (Y(pointeurCurrentImage(2))-Y2(end))^2)) > 16)...
%                 && (abs(sqrt((X_ind(pointeurCurrentImage(1),1)-X_ind(pointeurCurrentImage(2),1))^2 + (Y(pointeurCurrentImage(1))-Y(pointeurCurrentImage(2)))^2)) < 100)...
%                 )
%             
%             
%             
%             % Si une des deux billes est loin des deux positions précédente on
%             % demande l'information à l'utilisateur
%             
%             if ismember(i,Slog)
%                 button = ActionLog(i).button;
%             else
%                 
%                 if AUTO
%                     error('User input skipped in auto mode')
%                 end
%                 
%                 
%                 % afficher l'image de la stack et demander de cliquer sur les
%                 % deux billes d'interet
%                 I = imread(stackname,'tiff',i);
%                 I = imadjust(I);
%                 
%                 figure(1001)
%                 imshow(I)
%                 hold on
%                 plot(X1,Y1,'c--')
%                 plot(X2,Y2,'m--')
%                 plot(X1(end),Y1(end),'b*','markersize',5)
%                 plot(X2(end),Y2(end),'r*','markersize',5)
%                 plot(X_ind(pointeurCurrentImage),Y(pointeurCurrentImage),'o','color',[0 0.3 0],'markersize',5)
%                 title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(i) '/' num2str(nImages)])
%                 
%                 button = questdlg('Can you select the beads ?',...
%                     'Selection or next','No');
%                 
%                 
%                 ActionLog(i).button = button;
%                 
%             end
%             
%             if strcmp(button,'Yes')
%                 
%                 if ismember(i,Slog)
%                     Xim = ActionLog(i).pos(:,1);
%                     Yim = ActionLog(i).pos(:,2);
%                 else
%                     
%                     figure(1001)
%                     % A nouveau un pointage des billes par l'utilisateur
%                     % dans ce cas.
%                     [Xim,Yim] = ginput(2);
%                     close
%                     
%                     Slog = [Slog i];
%                     ActionLog(i).pos = [Xim,Yim];
%                     save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
%                 end
%                 
%                 % Attribution des deux points choisis par l'user à la
%                 % catégorie 'bille 1' ou 'billes 2'
%                 a = repmat(X_ind(pointeurCurrentImage,1),1,2);
%                 b = repmat(Xim',length(pointeurCurrentImage),1);
%                 
%                 c = repmat(Y(pointeurCurrentImage),1,2);
%                 d = repmat(Yim',length(pointeurCurrentImage),1);
%                 
%                 dist = sqrt((a-b).^2 + (c-d).^2);
%                 
%                 [~,Ind1] = min(dist(:,1));
%                 [~,Ind2] = min(dist(:,2));
%                 
%                 
%                 if abs(sqrt((X_ind(pointeurCurrentImage(Ind1),1)-X1(end))^2 + (Y(pointeurCurrentImage(Ind1))-Y1(end))^2)) ...
%                         < abs(sqrt((X_ind(pointeurCurrentImage(Ind1),1)-X2(end))^2 + (Y(pointeurCurrentImage(Ind1))-Y2(end))^2))
%                     % c'est la bille 1 la plus proche de X1
%                     
%                     X1 = [X1 X_ind(pointeurCurrentImage(Ind1),1)];
%                     S1 = [S1 i];
%                     ptr1 = [ptr1 pointeurCurrentImage(Ind1)];
%                     Y1 = [Y1 Y(pointeurCurrentImage(Ind1))];
%                     
%                     X2 = [X2 X_ind(pointeurCurrentImage(Ind2),1)];
%                     S2 = [S2 i];
%                     ptr2 = [ptr2 pointeurCurrentImage(Ind2)];
%                     Y2 = [Y2 Y(pointeurCurrentImage(Ind2))];
%                     
%                 else
%                     % c'est la bille 2
%                     X1 = [X1 X_ind(pointeurCurrentImage(Ind2),1)];
%                     S1 = [S1 i];
%                     ptr1 = [ptr1 pointeurCurrentImage(Ind2)];
%                     Y1 = [Y1 Y(pointeurCurrentImage(Ind2))];
%                     
%                     X2 = [X2 X_ind(pointeurCurrentImage(Ind1),1)];
%                     S2 = [S2 i];
%                     ptr2 = [ptr2 pointeurCurrentImage(Ind1)];
%                     Y2 = [Y2 Y(pointeurCurrentImage(Ind1))];
%                 end
%                 
%             elseif ismember(i,Slog)
%                 close
%             else
%                 close
%                 Slog = [Slog i];
%                 save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
%             end
%             
% % A QUOI SERT CE TEST ??            
%         elseif (abs(sqrt((X_ind(pointeurCurrentImage(1),1)-X1(end))^2 + (Y(pointeurCurrentImage(1))-Y1(end))^2)) > 16 ...
%                 && abs(sqrt((X_ind(pointeurCurrentImage(2),1)-X1(end))^2 + (Y(pointeurCurrentImage(2))-Y1(end))^2)) > 16) ...
%                 || (abs(sqrt((X_ind(pointeurCurrentImage(1),1)-X2(end))^2 + (Y(pointeurCurrentImage(1))-Y2(end))^2)) > 16 ...
%                 && abs(sqrt((X_ind(pointeurCurrentImage(2),1)-X2(end))^2 + (Y(pointeurCurrentImage(2))-Y2(end))^2)) > 16)
%             
% 
% % Les 2 cas simples ou les 2 billes ont peu bougé
%         elseif abs(sqrt((X_ind(pointeurCurrentImage(1),1)-X1(end))^2 + (Y(pointeurCurrentImage(1))-Y1(end))^2)) < ...
%                 abs(sqrt((X_ind(pointeurCurrentImage(1),1)-X2(end))^2 + (Y(pointeurCurrentImage(1))-Y2(end))^2))
%             % la bille 1 est plus proche de X1
%             
%             X1 = [X1 X_ind(pointeurCurrentImage(1),1)];
%             S1 = [S1 i];
%             ptr1 = [ptr1 pointeurCurrentImage(1)];
%             Y1 = [Y1 Y(pointeurCurrentImage(1))];
%             
%             X2 = [X2 X_ind(pointeurCurrentImage(2),1)];
%             S2 = [S2 i];
%             ptr2 = [ptr2 pointeurCurrentImage(2)];
%             Y2 = [Y2 Y(pointeurCurrentImage(2))];
%             
%         else
%             % la bille 1 est plus proche de X2
%             
%             X2 = [X2 X_ind(pointeurCurrentImage(1),1)];
%             S2 = [S2 i];
%             ptr2 = [ptr2 pointeurCurrentImage(1)];
%             Y2 = [Y2 Y(pointeurCurrentImage(1))];
%             
%             X1 = [X1 X_ind(pointeurCurrentImage(2),1)];
%             S1 = [S1 i];
%             ptr1 = [ptr1 pointeurCurrentImage(2)];
%             Y1 = [Y1 Y(pointeurCurrentImage(2))];
%             
%         end
%         
%         
%     elseif length(pointeurCurrentImage) >= 2 % deux ou plus bille sur l'image
%         
%         X1_tmp = [];
%         X2_tmp = [];
%         
%         for ii = 1:length(pointeurCurrentImage)
%             
%             
%             if abs(sqrt((X_ind(pointeurCurrentImage(ii),1)-X1(end))^2 + (Y(pointeurCurrentImage(ii))-Y1(end))^2)) > 16 ...
%                     && abs(sqrt((X_ind(pointeurCurrentImage(ii),1)-X2(end))^2 + (Y(pointeurCurrentImage(ii))-Y2(end))^2)) > 16
%                 % c'est une tache ou une autre bille
%                 
%             elseif abs(sqrt((X_ind(pointeurCurrentImage(ii),1)-X1(end))^2 + (Y(pointeurCurrentImage(ii))-Y1(end))^2)) < ...
%                     abs(sqrt((X_ind(pointeurCurrentImage(ii),1)-X2(end))^2 + (Y(pointeurCurrentImage(ii))-Y2(end))^2))
%                 
%                 
%                 X1_tmp = [X1_tmp X_ind(pointeurCurrentImage(ii),1)];
%                 S1_tmp = i;
%                 ptr1_tmp = pointeurCurrentImage(ii);
%                 Y1_tmp = Y(pointeurCurrentImage(ii));
%                 
%             else
%                 
%                 X2_tmp = [X2_tmp X_ind(pointeurCurrentImage(ii),1)];
%                 S2_tmp = i;
%                 ptr2_tmp = pointeurCurrentImage(ii);
%                 Y2_tmp = Y(pointeurCurrentImage(ii));
%             end
%         end
%         
%         
%         if (length(X1_tmp) == 1) && (length(X2_tmp) == 1) && not(X1_tmp == X2_tmp)
%             X1 = [X1 X1_tmp];
%             Y1 = [Y1 Y1_tmp];
%             S1 = [S1 S1_tmp];
%             ptr1 = [ptr1 ptr1_tmp];
%             
%             X2 = [X2 X2_tmp];
%             Y2 = [Y2 Y2_tmp];
%             S2 = [S2 S2_tmp];
%             ptr2 = [ptr2 ptr2_tmp];
%             
%         else
%             
%             
%             % afficher l'image de la stack et demander de cliquer sur les
%             % deux billes d'interet
%             
%             
%             if ismember(i,Slog)
%                 button = ActionLog(i).button;
%             else
%                 
%                 if AUTO
%                     error('User input skipped in auto mode')
%                 end
%                 
%                 
%                 % afficher l'image de la stack et demander de cliquer sur les
%                 % deux billes d'interet
%                 I = imread(stackname,'tiff',i);
%                 I = imadjust(I);
%                 figure(1000)
%                 imshow(I)
%                 hold on
%                 plot(X1,Y1,'c--')
%                 plot(X2,Y2,'m--')
%                 plot(X1(end),Y1(end),'b*','markersize',5)
%                 plot(X2(end),Y2(end),'r*','markersize',5)
%                 plot(X_ind(pointeurCurrentImage),Y(pointeurCurrentImage),'o','color',[0 0.3 0],'markersize',5)
%                 title(['Click the 2 beads to track. ' stackname ' Slice ' num2str(i) '/' num2str(nImages)])
%                 
%                 button = questdlg('Can you select the beads ?',...
%                     'Selection or next','No');
%                 
%                 
%                 ActionLog(i).button = button;
%                 
%             end
%             
%             if strcmp(button,'Yes')
%                 
%                 if ismember(i,Slog)
%                     Xim = ActionLog(i).pos(:,1);
%                     Yim = ActionLog(i).pos(:,2);
%                 else
%                     
%                     figure(1000)
%                     [Xim,Yim] = ginput(2);
%                     close
%                     
%                     Slog = [Slog i];
%                     ActionLog(i).pos = [Xim,Yim];
%                     save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
%                 end
%                 
%                 a = repmat(X_ind(pointeurCurrentImage,1),1,2);
%                 b = repmat(Xim',length(pointeurCurrentImage),1);
%                 
%                 c = repmat(Y(pointeurCurrentImage),1,2);
%                 d = repmat(Yim',length(pointeurCurrentImage),1);
%                 
%                 dist = sqrt((a-b).^2 + (c-d).^2);
%                 
%                 [~,Ind1] = min(dist(:,1));
%                 [~,Ind2] = min(dist(:,2));
%                 
%                 
%                 if abs(sqrt((X_ind(pointeurCurrentImage(Ind1),1)-X1(end))^2 + (Y(pointeurCurrentImage(Ind1))-Y1(end))^2)) ...
%                         < abs(sqrt((X_ind(pointeurCurrentImage(Ind1),1)-X2(end))^2 + (Y(pointeurCurrentImage(Ind1))-Y2(end))^2))
%                     % c'est la bille 1 la plus proche de X1
%                     
%                     X1 = [X1 X_ind(pointeurCurrentImage(Ind1),1)];
%                     S1 = [S1 i];
%                     ptr1 = [ptr1 pointeurCurrentImage(Ind1)];
%                     Y1 = [Y1 Y(pointeurCurrentImage(Ind1))];
%                     
%                     X2 = [X2 X_ind(pointeurCurrentImage(Ind2),1)];
%                     S2 = [S2 i];
%                     ptr2 = [ptr2 pointeurCurrentImage(Ind2)];
%                     Y2 = [Y2 Y(pointeurCurrentImage(Ind2))];
%                     
%                 else
%                     % c'est la bille 2
%                     X1 = [X1 X_ind(pointeurCurrentImage(Ind2),1)];
%                     S1 = [S1 i];
%                     ptr1 = [ptr1 pointeurCurrentImage(Ind2)];
%                     Y1 = [Y1 Y(pointeurCurrentImage(Ind2))];
%                     
%                     X2 = [X2 X_ind(pointeurCurrentImage(Ind1),1)];
%                     S2 = [S2 i];
%                     ptr2 = [ptr2 pointeurCurrentImage(Ind1)];
%                     Y2 = [Y2 Y(pointeurCurrentImage(Ind1))];
%                 end
%                 
%             elseif ismember(i,Slog)
%                 close
%             else
%                 close
%                 Slog = [Slog i];
%                 save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
%             end
%             
%             
%             
%             
%             
%         end
%         
%     end
%     
% end
% 
% save([filename '_AnalysisLog' tag '.mat'],'Slog','ActionLog');
% 
% 
% X1 = X1';
% X2 = X2';
% Y1 = Y1';
% Y2 = Y2';
% S1 = S1';
% S2 = S2';
% end




