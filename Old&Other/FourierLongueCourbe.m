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


figure
hold on
plot(T1,D1,'-*','linewidth',5)
plot(T,D,'-','linewidth',2)
title('Original Signal and Interpolation on a Grid for FFT')

% 
D = rand(1,length(T));

L = length(D); % length of signal


FTD = fft(D); % fourrier transform of D

P2 = abs(FTD/L); % two sided spectrum

P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); % one sided spectrum

Power = P1.^2/L;

f = Fs*(0:(L/2))/L; % frequency domain for P1

figure
hold on
plot(f(1:end),Power(1:end),'-')
title('Single sided power spectrum without 0 for Data1')
xlabel('Frequency (Hz)')
ylabel('Power')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear T D FTD f L

% interpolation of signal on a regular grid
T = 1:Ts:T2(end); % time

if isodd(T)
    T = T(1:end-1);
end

D = interp1(T2,D2,T);
D = fillmissing(D,'linear');

% find(isnan(D))
% find(isinf(D))

figure
hold on
plot(T2,D2,'-*','linewidth',5)
plot(T,D,'-','linewidth',2)
title('Original Signal and Interpolation on a Grid for FFT')


% parameters definition
L = length(D); % length of signal


FTD = fft(D); % fourrier transform of D

P2 = abs(FTD/L); % two sided spectrum

P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); % one sided spectrum

Power = P1.^2/L;

f = Fs*(0:(L/2))/L; % frequency domain for P1

figure
hold on
plot(f(1:end),Power(1:end),'-')
title('Single sided power spectrum without 0 for Data2')
xlabel('Frequency (Hz)')
ylabel('Power')