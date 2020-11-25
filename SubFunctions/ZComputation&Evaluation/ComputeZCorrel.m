function [Curve,MultiPos,MultiVal] = ComputeZCorrel(Li,KrefNorm)

% 
% [Curve,MultiPos,MultiVal] = ComputeZCorrel(Li,KrefNorm)
% 
% computes the correlation between one line of pixel from one bead and the
% depthograph. Returns correlation curve, and best two correlation with
% positions and value.
% 
% 

Lref = size(KrefNorm,1);

LiNorm = Li./mean(Li);

MLiNorm = repmat(LiNorm,Lref,1);

DiffK = KrefNorm - MLiNorm;

DiffKsq = DiffK.^2;

DiffKsqmean = mean(DiffKsq,2);

Curve = DiffKsqmean;

% figure
% plot(Curve)

ix = 1:length(Curve);
ixq = 1:0.1:length(Curve);

iDiffKsqsum = interp1(ix,smooth(Curve),ixq,'spline');

iDiffKsqsumExtended = [max(iDiffKsqsum) iDiffKsqsum max(iDiffKsqsum)];
ixqExtended = [0 ixq ixq(end)+1];

[val,locs] = findpeaks(-iDiffKsqsumExtended,ixqExtended,'minpeakdistance',20,'sortstr','descend');

if length(locs) < 2
    locs(2) = NaN;
    val(2) = 10;
end

MultiPos = locs(1:2)-1;
MultiVal = 1./abs(val(1:2));


end