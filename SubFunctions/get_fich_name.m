function [ListDep, ListFor] = get_fich_name(listP)

% 
% [ListDep, ListFor] = get_fich_name(listP)
%
% creates a list of possible identifiers (Px_Cy or Px_Cy-z) for cells based 
% list of chambers to analyze. 
%
% ListP = list of chambers (ex : 1, 1:3, [1 3])
%

nc = 100; % max number of vidéos 
ni = 5; % max number of interface


listC = 1:nc; 
listI = 0:ni; 

% initialize variables
ListDep = {};
ListFor = {};

% loop over all possible names
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
