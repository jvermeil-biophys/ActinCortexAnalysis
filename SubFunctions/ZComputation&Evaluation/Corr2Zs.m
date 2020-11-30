function [Zup,Wup,Zmid,Wmid,Zdown,Wdown] = Corr2Zs(Bool,Str,MultiPos,MultiVal,ValSum,fref)

% 
% Computes the best Z values for each beads based on correlation results 
% 
% [Zup,Wup,Zmid,Wmid,Zdown,Wdown] = Corr2Zs(Bool,Str,MultiPos,MultiVal,ValSum,fref)
%
% OUTPUTS : 
% Zup, Zmid, Zdown : Z position for up, mid, down image
% Wup, Wmid, Wdown : corresponding correlation value 
%
% INPUTS : 
% Bool : Boolean marker for the respect of alignement position between the
% image triplets
% Str : name of each combination
% MultiPos : list of positions for best two correlation
% Multival : list of correlation values
% ValSum : sum of correlation values for each combination
% fref : reference focus of depthograph
%

% Find combination that respect order
ptr = find(Bool ==1);
if isempty(ptr)
    ptr=1:8;
end

% amongst combination that respect order, take best correlation
[~,i] = max(ValSum(ptr));
n = ptr(i);

% assign result
k1 = str2double(Str{n}(1));
k2 = str2double(Str{n}(2));
k3 = str2double(Str{n}(3));

Zup   = fref-MultiPos(1,k1);
Zmid  = fref-MultiPos(2,k2);
Zdown = fref-MultiPos(3,k3);

Wup   = MultiVal(1,k1);
Wmid  = MultiVal(2,k2);
Wdown = MultiVal(3,k3);

end