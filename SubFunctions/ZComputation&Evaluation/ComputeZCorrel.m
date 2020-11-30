function [Curve,MultiPos,MultiVal] = ComputeZCorrel(Li,KrefNorm)

% 
% Computes the correlation between one line of pixel from one bead and the
% depthograph. Returns correlation curve, and best two correlation with
% positions and value. The idea between getting the two best correlation
% position is that for images of beads near the focus in Z, the depthograph
% is very symetrical, thus we take the info of the best correlation on each
% side, and then decide which to take based on the other images in the Z
% triplet when there are some. Otherwise best correlation is used (this is
% done in DzCalcCorr).
% 
% [Curve,MultiPos,MultiVal] = ComputeZCorrel(Li,KrefNorm)
%
% OUTPUTS : 
% Curve : Correlation curve
% MultiPos : position on two best correlations
% MultiVal : value of two best correlations
%
% INPUTS : 
% Li : line of pixel to correlate
% KrefNorm : normalized depthograph
% 

Lref = size(KrefNorm,1); % length of kimograph

LiNorm = Li./mean(Li); % normalized line of pixel

MLiNorm = repmat(LiNorm,Lref,1); % replication of the line of pixels as many times as there are line in the kimograph

DiffK = KrefNorm - MLiNorm; % difference between kimograph and replicated line of pixel

DiffKsq = DiffK.^2; % squared diff

DiffKsqmean = mean(DiffKsq,2); % average of squared distance for each line of the kimograph

Curve = DiffKsqmean; % better name ;)

ix = 1:length(Curve); % X coord of curve
ixq = 1:0.1:length(Curve); % finer X coord of curve

iDiffKsqsum = interp1(ix,smooth(Curve),ixq,'spline'); % fine interpolation of curve (for sub-voxel Z resoltuion)

iDiffKsqsumExtended = [max(iDiffKsqsum) iDiffKsqsum max(iDiffKsqsum)]; % Adding a point on both side of the curve to avoid peak detection at the edges
ixqExtended = [0 ixq ixq(end)+1]; % extended X coord of curve

[val,locs] = findpeaks(-iDiffKsqsumExtended,ixqExtended,'minpeakdistance',20,'sortstr','descend'); % minimum detection (peaks of inverted curve)

if length(locs) < 2 % if only one maximum the second position in NaN and the second correlation value is fixed to very much a lot in order to never be selected
    locs(2) = NaN;
    val(2) = Inf;
end

MultiPos = locs(1:2)-1; % storing
MultiVal = 1./abs(val(1:2)); % storing


end