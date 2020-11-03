function [Corr,LagsF] = AUTOCORR(T,D,Fs)

%% parameters definition
% Fs = 2; % sampling frequency in Hertz
Ts = 1/Fs; % sampling period in second

%% interpolation of signal on a regular grid
Tgrid = 1:Ts:T(end); % time


Dgrid = interp1(T,D,Tgrid);
Dgrid = fillmissing(Dgrid,'linear');

% plot interpolation

% figure
% hold on
% plot(T,D,'-*','linewidth',5)
% plot(Tgrid,Dgrid,'-','linewidth',2)
% title('Original Signal and Interpolation on a Grid for autocorr')


%% Autocorrelation

Dgrid = Dgrid - mean(Dgrid);

[Corr,Lags] = xcorr(Dgrid,'coeff');

Corr = Corr(Lags>=0);
Lags = Lags(Lags>=0);

LagsF =Lags/Fs;

end