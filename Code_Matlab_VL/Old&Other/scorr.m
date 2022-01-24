function [Sup scn] = scorr(Sup1,Sup2,Dat1,Dat2,sig1,sig2,w,methode,P)

% create a support for all data
S = union(Sup1,Sup2);
ls = length(S);

% resample data on the new support
D1r = interp1(Sup1,Dat1,S,'pchip');
D2r = interp1(Sup2,Dat2,S,'pchip');

% normalize signals
D1rn = (D1r-min(D1r))/(max(D1r-min(D1r)));
D2rn = (D2r-min(D2r))/(max(D2r-min(D2r)));

%% computes the correlation on an interval
% of size 2*w+1

scn = zeros(ls,1);
for i = 1:ls
    % Affect variables in function of position
    % in the initial vector
    if i < w+1
        D1 = D1rn(1:i+w);
        D2 = D2rn(1:i+w);        
    elseif i+w > ls
        D1 = D1rn(i-w:end);
        D2 = D2rn(i-w:end);
    else
        D1 = D1rn(i-w:i+w);
        D2 = D2rn(i-w:i+w);
    end
    % compute the correlation
    if strcmp(methode,'static')
    scn(i) = xcorr(D1,D2,0);
    elseif strcmp(methode,'statistic')
        scn(i) = (mean((D1-mean(D1)).*(D2-mean(D2))))/...
            (sqrt(mean((D1-mean(D1)).^2))*sqrt(mean((D2-mean(D2)).^2)));
    else
        error('Wrong methode, it must be ''static'' or ''statistic''.');
    end
    
end

Sup = S;

if P
% display results
figure

subplot(312)
plot(S,D2r,'color',[0.6 0.2 0.8])
ylabel(sig2)

subplot(313)
plot(S,scn,'r')
ylabel('normalized correlation coeff')
xlabel('Time points')

subplot(311)
plot(S,D1r,'b')
ylabel(sig1)
if strcmp(methode,'static')
title(['Static correlation on ' num2str(w) ' time points'])
elseif strcmp(methode,'statistic')
    title(['Statistic correlation on ' num2str(w) ' time points'])
end
end

end
