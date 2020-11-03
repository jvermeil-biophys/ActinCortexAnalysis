function [xp,yp,Txy] = findSlice_xy(S1,Sf1,X1,Y1,S2,Sf2,X2,Y2)

[S,p1,p2] = intersect(S1,S2);
X = [X1(p1) X2(p2)];
Y = [Y1(p1) Y2(p2)];

Sf = [Sf1, Sf2];
Sxy = mean(Sf,2);
lxy = length(Sxy);

MS = repmat(S,1,lxy);
MSxy = repmat(Sxy',length(S),1);

DMS = abs(MSxy - MS);

[mDMS,Ind] = min(DMS);

xp = X(Ind,:)';
yp = Y(Ind,:)';

Txy = S(Ind);

% 
% MS1 = repmat(S1,1,length(Sf1));
% MSf1 = repmat(Sf1',length(S1),1);
% 
% DMS1 = abs(MSf1 - MS1);
% 
% [mDMS1,Ind1] = min(DMS1);
% 
% MS2 = repmat(S2,1,length(Sf2));
% MSf2 = repmat(Sf2',length(S2),1);
% 
% DMS2 = abs(MSf2 - MS2);
% 
% [mDMS2,Ind2] = min(DMS2);
% 
% xp(1,:) = X1(Ind1);
% xp(2,:) = X2(Ind2);
% 
% yp(1,:) = Y1(Ind1);
% yp(2,:) = Y2(Ind2);

end
