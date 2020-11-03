function [Corr,LagsF] = FFTandAUTOCORR(T,D,Fs)

%% parameters definition
% Fs = 2; % sampling frequency in Hertz
Ts = 1/Fs; % sampling period in second

%% interpolation of signal on a regular grid
Tgrid = 1:Ts:T(end); % time

if isodd(Tgrid)
    Tgrid = Tgrid(1:end-1);
end

Dgrid = interp1(T,D,Tgrid);
Dgrid = fillmissing(Dgrid,'linear');

% plot interpolation

% figure
% hold on
% plot(T,D,'-*','linewidth',5)
% plot(Tgrid,Dgrid,'-','linewidth',2)
% title('Original Signal and Interpolation on a Grid for FFT')


%% FFT
% L = length(Dgrid); % length of signal
% 
% FTD = fft(Dgrid); % fourrier transform of D
% 
% P2 = abs(FTD/L); % two sided spectrum
% 
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1); % one sided spectrum
% 
% Power = P1.^2/L;
% 
% f = Fs*(0:(L/2))/L; % frequency domain for P1
% 
% if POOL
% figure(1)
% else
%     figure
% end
% hold on
% plot(f,Power,'-*r')
% title(['Single sided power spectrum for ' Name])
% xlabel('Frequency (Hz)')
% ylabel('Power')

%% Autocorrelation

Dgrid = Dgrid - mean(Dgrid);

[Corr,Lags] = xcorr(Dgrid,'coeff');

Corr = Corr(Lags>0);
Lags = Lags(Lags>0);

LagsF =Lags/Fs;

end