%% Advanced Signal Processing
%% Coursework Assignment 3 - Spectral Estimation and Modelling
%% 3.1 - Periodogram estimate

N = 256;
x = randn(1,N);
x_pgm = pgmest(x);

figure('Position',[100,100,400,300]);
plot(1:N,x);title("WGN for data length N = "+ N);
xlabel('n');ylabel('x[n]');grid on;
figure('Position',[100,100,400,300]);plot(-1:2/N:1-2/N,x_pgm);title("Periodogram for N = "+ N);xlabel('Normalised frequency (\times2\pirad)');ylabel('P_{X}(f)');grid on;

%% Averaged periodogram estimates

% filter coefficients
a = 1;
b = 0.2*[1,1,1,1,1];

x_pgm_filt = fftshift(filter(b,a,x_pgm));
figure('Position',[100,100,400,300]);
plot(-1:2/N:1-2/N,x_pgm_filt);
title("Filtered periodogram for N = "+ N);
xlabel('Normalised frequency (\times2\pirad)');
ylabel('P_{X}(f)');
grid on;
ylim([0,12]);


y = randn(1,1024);
y_seg = ones(8,128);
y_seg_pgm = ones(8,128);

%% split signal into 8 segments

figure;
for i = 1:8
    y_seg(i,:) = y(128*(i-1)+1:128*i);
    y_seg_pgm(i,:) = pgmest(y_seg(i,:));
    plots = subplot(4,2,i);
    plot(-1:2/128:1-2/128,y_seg_pgm(i,:));
    xlabel('Normalised frequency (\times2\pirad)');
    ylabel('P_{X}(f)');
    grid on;
    sgtitle('Periodograms for 8 segments of 1024-sample WGN process');
end

avg_pgm = mean(y_seg_pgm);

figure('Position',[100,100,400,300]);
plot(-1:2/128:1-2/128,avg_pgm);
title("Average of 8 periodograms for 128-sample WGN segments");
xlabel('Normalised frequency (\times2\pirad)');
ylabel('P_{X}(f)');
grid on;

%% 3.2 - Spectrum of AR processes

z = randn(1,1024);

figure('Position',[100,100,800,200]);
plot(1:1024,z);
title("1024-sample WGN process");
xlabel('n');
ylabel('x[n]');
grid on;
xlim([0,1024]);

%% filter coefficients

a = [1,0.9];
b = 1;
z_filt = filter(b,a,z);
z_filt_pgm = pgmest(z_filt);
[h,w] = freqz(b,a,512);

figure('Position',[100,100,600,400]);
plot(w/(2*pi),abs(h).^2,'LineWidth',0.8);
title('PSD and periodogram of filtered 1024-sample WGN signal');xlabel('Normalised frequency f (\times2\piradians)');ylabel('P_{Y}(f)');grid on;
hold on;
plot(0:1/1024:0.5,z_filt_pgm(1:513),'LineWidth',0.8);
legend('PSD','Periodogram');

%% model based PSD estimation

z_acf = xcorr(z_filt, 'unbiased');
a_1 = -z_acf(2)/z_acf(1);
var = z_acf(1) + a_1 * z_acf(2);
[h2,w2] = freqz(var^(1/2), [1,a_1], 512);
z_pgm = pgmest(z);

figure('Position',[100,100,400,300]);
plot(0:1/1024:0.5,z_filt_pgm(1:513),'LineWidth',0.7);
grid on;
title('Model-based PSD and periodogram')
hold on;
plot(w/(2*pi),abs(h2).^2,'LineWidth',0.7);
xlabel('Normalized Frequency (\times2\pirad)')
ylabel('Amplitude')
grid on;
legend('Periodogram','Model-based PSD');

%% PSD estimates for sunspot data

load sunspot.dat
sun = sunspot(:,2);
sun = zscore(sun);
sun_acf = xcorr(sun,'unbiased');

% AR(1)
a_sun_o1 = -sun_acf(2)/sun_acf(1);
var_sun_o1 = sun_acf(1) + a_sun_o1 * sun_acf(2);
[h_sun_o1,w_sun_o1] = freqz(var_sun_o1^(1/2), [1,a_sun_o1], 288);

% AR(2)
a_sun_o2 = -sun_acf(3)/sun_acf(2);
var_sun_o2 = sun_acf(2) + a_sun_o2 * sun_acf(3);
[h_sun_o2,w_sun_o2] = freqz(var_sun_o2^(1/2), [1,a_sun_o2], 288);

% AR(3)
a_sun_o3 = -sun_acf(4)/sun_acf(3);
var_sun_o3 = sun_acf(3) + a_sun_o3 * sun_acf(4);
[h_sun_o3,w_sun_o3] = freqz(var_sun_o3^(1/2), [1,a_sun_o3], 288);

sun_pgm = pgmest(sun);

figure('Position',[100,100,480,360]);plot(0:1/288:0.5,sun_pgm(1:145),'LineWidth',0.7);grid on;title('Model-based PSD and periodogram for standardised sunspot data')
hold on;plot(w_sun_o1/(2*pi),abs(h_sun_o1).^2,'LineWidth',0.7);
plot(w_sun_o2/(2*pi),abs(h_sun_o2).^2,'LineWidth',0.7);
plot(w_sun_o3/(2*pi),abs(h_sun_o3).^2,'LineWidth',0.7);
xlim([0,0.5]);
xlabel('Normalized Frequency (\times2\pirad)');
ylabel('Amplitude')
legend('show');
grid on;
legend('Periodogram','AR(1) model-based PSD','AR(2) model-based PSD','AR(3) model-based PSD');

%% 3.3 - LSE of AR coefficients

load sunspot.dat
sun = sunspot(:,2);
sun = zscore(sun);

%% 3.4 - Dial-tone pad

num = randperm(9,8);
num = [0,2,0,num];

tones = [1336, 1209, 1336, 1447, 1209, 1336, 1447, 1209, 1336, 1447;
         941,  697,  697,  697,  770,  770,  770,  852,  852,  852];
     
x = linspace(0, 0.25, 8192); % gives sampling frequency of f 32768 Hz
dial_seq = []; % output sequence

% generate sequence of 10 digits and 10 gaps
for i=1:10
    dig = sin(2*pi*tones(1, num(i)+1) * x) + sin(2 * pi * tones(2, num(i)+1) * x); 
    gap = zeros(1, 8192);
    dial_seq = [dial_seq,dig,gap];
end

% final digit separate since no gap after
dig = sin(2*pi*tones(1, num(11)+1) * x) + sin(2 * pi * tones(2, num(11)+1) * x); 
dial_seq = [dial_seq,dig];

% noise
sd_noise = 10;
wgn = randn(1,length(dial_seq))*sd_noise;
dial_seq = dial_seq + wgn;

figure('Position',[100,100,800,600]);
plot(5.25/length(dial_seq):5.25/length(dial_seq):5.25,dial_seq-5.25/length(dial_seq));
xlim([0,5.25]);
grid on;
title("Dial tone sequence for London landline number "+num2str(num)+" (with noise of SD \sigma_{n} = "+sd_noise+")");
xlabel("time (s)");
ylabel("y[n]");

%% FFT 

% calculate spectrogram data for time and frequency
[s,f,t] = spectrogram(dial_seq, hann(8192), 0, 8192, 32768);

% plot in dB for magnitude spectrum against frequency
figure('Position',[100,100,800,600]);
hold on;
plot(f, mag2db(abs(s(:,1)))); 
plot(f, mag2db(abs(s(:,3))));
title("Magnitude spectrum of selected segments of signal (with noise of SD \sigma_{n} = "+sd_noise+")")
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Number 0','Number 2');
grid on;

%% spectrogram

figure('Position',[100,100,800,600]);hold on;
spectrogram(dial_seq,hann(8192), 0, 8192, 32768, 'yaxis');
ylim([0.25,1.75]);xlim([0,5.25])
title("Spectrogram of dial tone signal (with noise of SD \sigma_{n} = "+sd_noise+")");
xlabel('time (seconds)');

%% pgm function definition

function [pgm_out] = pgm(x)
    N = length(x);
    pgm_out = (1/N)*abs(fft(x)).^2; 
end

