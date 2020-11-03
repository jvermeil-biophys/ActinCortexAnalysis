clear all
close all
load('D:\Data\MatFile\PourNico.mat')


% parameters definition
Fs = 2; % sampling frequency in Hertz
Ts = 1/Fs; % sampling period in second

% interpolation of signal on a regular grid
T = 1:Ts:T1(end); % time

if isodd(T)
    T = T(1:end-1);
end

D = interp1(T1,D1,T);
D = fillmissing(D,'linear');

% D = sin(T);

D = D - mean(D);


[Corr,Lags] = xcorr(D,round(T(end)/2),'coeff');

figure(1)
hold on
plot(Lags(Lags>0)/Fs,Corr(Lags>0))
xlabel('Lag (sec)')
ylabel('AutoCorrelation')


% parameters definition
Fs = 2; % sampling frequency in Hertz
Ts = 1/Fs; % sampling period in second

% interpolation of signal on a regular grid
T = 1:Ts:T2(end); % time

if isodd(T)
    T = T(1:end-1);
end

D = interp1(T2,D2,T);
D = fillmissing(D,'linear');


D = D - mean(D);

[Corr,Lags] = xcorr(D,round(T(end)/2),'coeff');

figure(1)
hold on
plot(Lags(Lags>0)/Fs,Corr(Lags>0))
xlabel('Lag (sec)')
ylabel('AutoCorrelation')
title('Autocorrelation for signal 1 and 2')