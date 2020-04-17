%% Advanced Signal Processing
%% ECG/RRI exercises both in Sections 2 and 3
%% 2.5 - Real world signals: ECG from iAmp experiment

load('RAW.mat');
figure;
plot(time,data);

% split into sections
ecg_uncon = data(138400:285400);
figure;
plot(time(138400:285400),ecg_uncon);
xlabel('time (s)');

ecg_con15 = data(285401:416900);
figure;
plot(time(285401:416900),ecg_con15);
xlabel('time (s)');

ecg_con50 = data(495800:640900);
figure;
plot(time(495800:640900),ecg_con50);
xlabel('time (s)');

% convert sections to rri signals
[rri_uncon,fs_rri_uncon] = ECG_to_RRI(ecg_uncon,fs);
[rri_con15,fs_rri_con15] = ECG_to_RRI(ecg_con15,fs);
[rri_con50,fs_rri_con50] = ECG_to_RRI(ecg_con50,fs);

%% PDE for unconstrained breathing

alpha = 1;
N = 10;

hr_uncon = 60./rri_uncon;
hr_uncon_ave = (alpha/N).*hr_uncon(1:N:length(hr_uncon));
hr_uncon_ave = filter(alpha*ones(1,10)/10,1,hr_uncon);
plot(hr_uncon_ave);hold on
plot(hr_uncon);

figure('Position',[250 250 800 200]);histogram(hr_uncon,60,'Normalization','pdf');
title('PDE for HR in unconstrained breathing trial (original data)');
xlim([0,100]);
xlabel('Heart rate (bpm)');
ylabel('Probability');
grid on;

figure('Position',[250 250 800 200]);histogram(hr_uncon_ave,60,'Normalization','pdf');
title("PDE for HR in unconstrained breathing trial (averaged data with \alpha = "+alpha+")");
xlim([0,100]);
xlabel('Heart rate (bpm)');
ylabel('Probability');
grid on;


%% PACFs of all trials

% normalise data
rri_uncon_norm = detrend(rri_uncon);
rri_con15_norm = detrend(rri_con15);
rri_con50_norm = detrend(rri_con50);


%% ACF sequences of all trials

acf_uncon = xcorr(rri_uncon,'unbiased');
acf_con15 = xcorr(rri_con15,'unbiased');
acf_con50 = xcorr(rri_con50,'unbiased');

acf_uncon_norm = xcorr(rri_uncon_norm,'unbiased');
acf_con15_norm = xcorr(rri_con15_norm,'unbiased');
acf_con50_norm = xcorr(rri_con50_norm,'unbiased');

% plot figures

figure('Position',[250,250,800,400]);
plot(-length(acf_uncon)/2:length(acf_uncon)/2-1,520.*acf_uncon_norm,'DisplayName','original');
title('ACF of RRI signal for unconstrained breathing')
xlim([-500,500])
xlabel('\tau');
ylabel('R(\tau)')
grid on;

figure('Position',[250,250,800,400]);
plot(-length(acf_con15_norm)/2:length(acf_con15_norm)/2-1,430.*acf_con15_norm,'DisplayName','original');
title('ACF of RRI signal for constrained breathing at 15 bpm')
xlabel('\tau');
ylabel('R(\tau)')
xlim([-450,450])
grid on;

figure('Position',[250,250,800,400]);
plot(-length(acf_con50_norm)/2:length(acf_con50_norm)/2-1,175.*acf_con50_norm,'DisplayName','original');
title('ACF of RRI signal for constrained breathing at 50 bpm')
xlabel('\tau');
ylabel('R(\tau)')
xlim([-525,525])
grid on;

%% PACFs of all trials

pacf_length = 20;

% create pacf vectors
pacf_uncon = zeros(pacf_length,1);
pacf_con15 = zeros(pacf_length,1);
pacf_con50 = zeros(pacf_length,1);
pacf_uncon_norm = zeros(pacf_length,1);
pacf_con15_norm = zeros(pacf_length,1);
pacf_con50_norm = zeros(pacf_length,1);

for i = 1:pacf_length
    uncon_coeffs = aryule(rri_uncon, i);
    pacf_uncon(i) = -uncon_coeffs(i+1);
    con15_coeffs = aryule(rri_con15, i);
    pacf_con15(i) = -con15_coeffs(i+1);
    con50_coeffs = aryule(rri_con50, i);
    pacf_con50(i) = -con50_coeffs(i+1);
end

for i = 1:pacf_length
    uncon_coeffs_norm = aryule(rri_uncon_norm, i);
    pacf_uncon_norm(i) = -uncon_coeffs_norm(i+1);
    con15_coeffs_norm = aryule(rri_con15_norm, i);
    pacf_con15_norm(i) = -con15_coeffs_norm(i+1);
    con50_coeffs_norm = aryule(rri_con50_norm, i);
    pacf_con50_norm(i) = -con50_coeffs_norm(i+1);
end

k = linspace(1, pacf_length, pacf_length);
figure('Position',[250 250 800 200]);
stem(k,pacf_uncon,'DisplayName','Original');
hold on;
stem(k,pacf_uncon_norm,'DisplayName','Zero-mean');
title("PACF of RRI for unconstrained breathing");
xlabel('k');
ylabel('\alpha (k)');
grid on;
legend show;

figure('Position',[250 250 800 200]);
stem(k,pacf_con15,'DisplayName','Original');
hold on;
stem(k,pacf_con15_norm,'DisplayName','Zero-mean');
title("PACF of RRI for constrained breathing at 15 bpm");
xlabel('k');
ylabel('\alpha (k)');
grid on;
legend show;

figure('Position',[250 250 800 200]);stem(k,pacf_con50,'DisplayName','Original');hold on;stem(k,pacf_con50_norm,'DisplayName','Zero-mean');
title("PACF of RRI for constrained breathing at 50 bpm");
xlabel('k');
ylabel('\alpha (k)');
grid on;
legend show;

%% AIC, MDL and AICc

signal = rri_uncon_norm;
p = 10;
err = zeros(N,p);
order = ones(1,p);
N = length(signal);
x = zeros(N,10);

for i = 1:p
    ar = aryule(signal(1:N), i);
    x(:,i) = filter(-1*ar(1:end),1,signal);
end

err=zeros(N,p);
order = ones(1,p);
err(1,:) = (x(1,:)-order*signal(1)).^2;

for i=2:N 
    err(i,:) = (x(i,:)-order*signal(i)).^2 + err(i-1,:);
end

mdl = log(err(N,(1:p))) + (1:p)*log(N)/N;
aic = log(err(N,(1:p))) + (1:p)*2/N;
aicc = zeros(1,p);

for i=1:10
    aicc(i) = aic(i) + (2*(i))*((i)+1)/(N-(i)-1);
end

figure('Position',[10,10,800,200]);
hold on;
plot(aic,'DisplayName','AIC');
plot(mdl,'DisplayName','MDL');
plot(aicc,'DisplayName','AIC_{c}');
grid on;
xlabel('Model order p');
ylabel('Magnitude');
title(['MDL, AIC and corrected AIC of RRI for unconstrained breathing model orders p = 0 to p = ' num2str(p)]);
legend show;

%% 3.5 - Respiratory Sinus Arrhythmia

pgm_uncon = pgmest(rri_uncon_norm);
pgm_con15 = pgmest(rri_con15_norm);
pgm_con50 = pgmest(rri_con50_norm);

figure('Position',[250 250 800 400]);
plot(-0.5:1/length(pgm_uncon):0.5-1/length(pgm_uncon),pgm_uncon,'LineWidth',0.8);
title('Periodogram of RRI for unconstrained breathing');
xlabel('Normalised frequency f');
ylabel('Periodogram P(f)');
grid on;
xlim([-0.3,0.3]);

figure('Position',[250 250 800 400]);
plot(-0.5:1/length(pgm_con15):0.5-1/length(pgm_con15),pgm_con15,'LineWidth',0.8);
title('Periodogram of RRI for breathing constrained at 15 bpm');
xlabel('Normalised frequency f');
ylabel('Periodogram P(f)');
grid on;
xlim([-0.3,0.3]);

figure('Position',[250 250 800 400]);
plot(-0.5:1/length(pgm_con50):0.5-1/length(pgm_con50),pgm_con50,'LineWidth',0.8);
title('Periodogram of RRI for breathing constrained at 50 bpm');
xlabel('Normalised frequency f');
ylabel('Periodogram P(f)');
grid on;
xlim([-0.3,0.3]);

%% averaged periodograms

sig = rri_con50_norm;
sig_seg = zeros(8,50);
sig_seg_pgm = zeros(8,50);

for i = 1:8
    sig_seg(i,:) = sig(50*(i-1)+1:50*i);
    sig_seg_pgm(i,:) = pgmest(sig_seg(i,:));
end

avg_pgm = mean(sig_seg_pgm);
figure('Position',[250 250 800 400]);
plot(-0.5:1/length(avg_pgm):0.5-1/length(avg_pgm),avg_pgm,'LineWidth',0.8);
title('Average periodogram of RRI for breathing constrained at 50bpm with window size 50');
xlabel('Normalised frequency f');
ylabel('Periodogram P(f)');
grid on;
xlim([-0.3,0.3])


% periodogram function
function [p_out] = pgmest(x)
    N = length(x);
    p_out = (1/N)*abs(fft(x)).^2; 
end




