function [f,Power] = FFTbasic(T,D,Fs)

Ts = 1/Fs; % sampling period in second

%% interpolation of signal on a regular grid
Tgrid = 1:Ts:T(end); % time

if isodd(Tgrid)
    Tgrid = Tgrid(1:end-1);
end

Dgrid = interp1(T,D,Tgrid);
Dgrid = fillmissing(Dgrid,'linear');

%% FFT
L = length(Dgrid); % length of signal

FTD = fft(Dgrid); % fourrier transform of D

P2 = abs(FTD/L); % two sided spectrum

P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); % one sided spectrum

Power = P1.^2/L;

f = Fs*(0:(L/2))/L; % frequency domain for P1



end