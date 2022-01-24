function CasTestActi()

%% INIT
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultAxesFontSize',30)

%% Params
RandPos = 0;
RandDur = 0;
RandHeight = 0;
base = 200;
timelength = 600;
sampling = 0.6;
heightdurfreq = [30 8 1]; %; 80  8 0.7; 200 12 0.1];

T = (1:round(timelength/sampling))*sampling;

% Baseline
Curve = ones(1,round(timelength/sampling))*base;

% PeaksAddition
ncat = size(heightdurfreq,1);

bigT = [-T(end:-1:1) T T(2:end)+T(end)];


for kc = 1:ncat
    freq = heightdurfreq(kc,3)/60;
    npeaks = round(timelength*3*freq);
    startpos = (1:npeaks)/freq-600;
    peakspos = startpos + (rand(1,npeaks)-0.5)*2*RandPos/(2*freq);
    
    peaks = zeros(1,length(bigT));
    
    for kp = 1:npeaks
        dur = heightdurfreq(kc,2)*((rand(1)-0.5)*RandDur+1);
        height = heightdurfreq(kc,1)*abs(randn(1)*RandHeight+1.2);
        peaks = peaks + gaussmf(bigT,[dur peakspos(kp)])*height;
        start = round(rand(1)*1900)+1;
        
    end
    Curve = Curve + peaks(start:start+999);
    clear peaks
end


figure(1)
subplot(211)
hold on
plot(T,Curve,'-*')
hold on
xlabel('T (s)')
ylabel('H (nm)')

% courbe d'activité

tau = 1:timelength; % un tau = 0.5 seconde

t = ceil(T(1)):0.5:floor(T(end));
x = interp1(T,Curve,t);
Acti = tau*0;

for kt = 1:length(tau)
    it = 1:length(t)-tau(kt); % it représente le temps, chaque pas vaut 0.5s
    sdx = (x(it) - x(it+tau(kt))).^2;
    Acti(kt) = sqrt(sum(sdx)/length(it));
end

figure(1)
subplot(212)
hold on
plot(tau/2,Acti,'*-')
ylabel('Activity (nm)')
xlabel('Tau (s)')

end