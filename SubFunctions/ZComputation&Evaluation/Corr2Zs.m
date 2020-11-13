function [Zup,Wup,Zmid,Wmid,Zdown,Wdown] = Corr2Zs(Bool,Str,MultiPos,MultiVal,ValSum,fref)

% [Zup,Wup,Zmid,Wmid,Zdown,Wdown] = Corr2Zs(Bool,Str,MultiPos,MultiVal,ValSum,fref)
% 
% Computes the best Z values for each beads based on correlation results 
% 

ptr = find(Bool ==1);
if isempty(ptr)
    ptr=1:8;
end

[~,i] = max(ValSum(ptr));
n = ptr(i);

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