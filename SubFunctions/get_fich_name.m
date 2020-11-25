function [ListDep, ListFor] = get_fich_name(listP)

% 
% [ListDep, ListFor] = get_fich_name(listP)
%
% creates a list of possible identifiers for cells based on the numbers of
% selected chambers (listP)

nc = 30; % vid�o par puits
ni = 15; % interface par vid�o


listC = 1:nc; 
listI = 0:ni; 

ListDep = {};
ListFor = {};

for kp = 1:length(listP)
    for kc = listC
        for ki = listI
            num = (ki+1) + (kc-1)*(ni+1) + (kp-1)*((ni+1)*nc);
            if ki == 0
                ListDep{num} = ['P' num2str(listP(kp)) '_C' num2str(kc)];
            else
                ListDep{num} = ['P' num2str(listP(kp)) '_C' num2str(kc) '-' num2str(ki)];             
            end            
                ListFor{num} = ['P' num2str(listP(kp)) '_C' num2str(kc)];
        end
    end
end



end
