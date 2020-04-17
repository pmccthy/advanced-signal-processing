% Advanced Signal Processing 
% Coursework Assignment 2 - Linear Stochastic Modelling
%% 2.1 - ACF of uncorrelated and correlated sequences

x = randn(1,1000);
tau = -999:999;
acf = xcorr(x,'unbiased');

figure('Position',[10,10,800,300])
plot(tau,acf);
xlabel('\tau');
ylabel('R_{X}[\tau]');
title('Unbiased ACF estimate for WGN');
grid on;
xlim([-999,999]);

y = filter(ones(9,1),[1],x);
acf_y = xcorr(y,'unbiased');
figure('Position',[10,10,800,300]);
stem(tau(980:1020),acf_y(980:1020));
xlabel('\tau');ylabel('R_{Y}[\tau]');
title('ACF of WGN filtered by 5th order MA filter');
grid on;
xlim([-20,20]);
ylim([-2,10]);

%% 2.2 - CCF

ccf = xcorr(x,y,'unbiased');

figure('Position',[10,10,800,300]);
stem(tau(980:1020),ccf(980:1020));
xlabel('\tau');
ylabel('R_{XY}[\tau]');
title('CCF of unfiltered WGN and WGN filtered by 9th order MA filter');
grid on;
xlim([-20,20]);

%% 2.3 - AR modelling
N = 1000;

a_1 = -2.5 + 5.*rand(N,100);
a_2 = -1.5 + 3.*rand(N,100);
x = rand(1,N);

stable = [];

% check if signal converges and plot asterisk * if it does
for i = 1:N
    if (a_1(i)+a_2(i)<1 && (a_2(i)-a_1(i)<1) && (a_2(i)>-1) && (a_2(i)<1))
        stable = [stable,[a_1(i);a_2(i)]];
    end
end

hold on;
plot(stable(1,:),stable(2,:),'*','DisplayName','Stable');
xlim([-2.5,2.5]);
ylim([-1.5,1.5]);
grid on;
xlabel('a_{1}');
ylabel('a_{2}');
title('Stable filter coefficients for N=1000')

%% Sunspot time series
load sunspot.dat;

figure;
plot(sunspot(:,1),sunspot(:,2));
title('Sunspot Time Series');
xlabel('Year');
ylabel('Wolf Number');
grid on;

N = 5;
sun = sunspot(:,2)';
[a,b] = xcorr(sun(1:N),'unbiased');

figure('Position',[10,10,800,300]);
stem(b,a);
xlabel('\tau');
ylabel('R_{S}[\tau]');
title('ACF for sunspot data length N=5');
grid on;
ylim([0,500]);

sun_norm = zscore(sun(1:N));
[c,d] = xcorr(sun_norm(1:N),'unbiased');

figure('Position',[10,10,800,300]);
stem(d,c);
xlabel('\tau');
ylabel('R_{S}[\tau]');
title('ACF for standardised sunspot data length N=5');
grid on;

%% Yule-Walker equations

% calculate AR coefficients for orders up to 10 using YW eqns
ar1 = aryule(sun,1);
ar2 = aryule(sun,2);
ar3 = aryule(sun,3);
ar4 = aryule(sun,4);
ar5 = aryule(sun,5);
ar6 = aryule(sun,6);
ar7 = aryule(sun,7);
ar8 = aryule(sun,8);
ar9 = aryule(sun,9);
ar10 = aryule(sun,10);

% plot coeffs
figure('Position',[10,10,800,300]);
hold on
grid on
plot(ar1,'DisplayName','1');
plot(ar2,'DisplayName','2');
plot(ar3,'DisplayName','3');
plot(ar4,'DisplayName','4');
plot(ar5,'DisplayName','5');
plot(ar6,'DisplayName','6');
plot(ar7,'DisplayName','7');
plot(ar8,'DisplayName','8');
plot(ar9,'DisplayName','9');
plot(ar10,'DisplayName','10');
xlabel('Model order')
ylabel('Model coefficients')
title('Yule-Walker coefficients for sunspot data for models up to order p=10');
legend show;

%% PACF

load sunspot.dat

p = 10; % model order 
N = 100;

sun = sunspot(:,2);
sun_norm = zscore(sun(1:N));
pcf_lengt = 10;
pcf = zeros(pcf_length,1);
pcf_norm = zeros(pcf_length,1);

for i = 1:pcf_length
    coefs = aryule(sun, i);
    pcf(i) = -coefs(i+1);
end

for i = 1:pcf_length
    coefs_norm = aryule(sun_norm, i);
    pcf_norm(i) = -coefs_norm(i+1);
end

x = linspace(1, pcf_length, pcf_length);
figure('Position',[10,10,800,300]);
stem(x,pcf,'DisplayName','PACF');
hold on;
stem(x,pcf_norm,'DisplayName','Standardised PACF');
grid on;
xlabel('k');
ylabel('\alpha (k)');
title('PACF of sunspot time series');
legend('show');

%% determining model order

p = 10;
sun = sunspot(:,2);
sun_norm = zscore(sun);
N = length(sun_norm(:,1));
x = zeros(N,10);

for i = 1:10
    ar = aryule(sun_norm(1:N), i);
    x(:,i) = filter(-1*ar(1:end), 1 ,sun_norm);
end

err = zeros(N,p);
o = ones(1,p);
err(1,:) = (x(1,:)-o*sun_norm(1)).^2;

for i=2:N 
    err(i,:) = (x(i,:)-o*sun_norm(i)).^2 + err(i-1,:);
end

mdl = log(err(N,(1:p))) + (1:p)*log(N)/N;
aic = log(err(N,(1:p))) + (1:p)*2/N;
aicc = zeros(1,p);

for i=1:10
    aicc(i) = aic(i) + (2*(i))*((i)+1)/(N-(i)-1);
end

figure('Position',[10,10,800,300]);
hold on;
plot(aic,'DisplayName','AIC');
plot(mdl,'DisplayName','MDL');
plot(aicc,'DisplayName','AIC_{c}')
%axis([0 10 -0.5 7]); % Overview
axis([1 10 6 7]); % Zoom
grid on;
xlabel('Model order p');
ylabel('Magnitude');
title(['MDL, AIC and corrected AIC for model orders p = 0 to p = ' num2str(p)]);
legend show

%% AR model

clear all; close all; clc;

load sunspot.dat
sun = sunspot(:,2);

N = 50; % data sample length
n = 1:N;

% model order 1, prediction horizon 1
mo1_ph1 = zeros(N,1); % vector to hold prediciton
mod_1_1 = ar(sun(1:N), 1, 'yw'); % create coefficients for AR model of order 1 using Yule-Walker equations
mo1_ph1 = predict(mod_1_1, sun(1:N),1); % predict next time step (prediction horizon 1)

% model order 1, prediction horizon 2
mo1_ph2 = zeros(N,1); 
mod_1_2 = ar(sun(1:N), 1, 'yw'); 
mo1_ph2 = predict(mod_1_2, sun(1:N),2); 

% model order 1, prediction horizon 5
mo1_ph5 = zeros(N,1); 
mod_1_5 = ar(sun(1:N), 1, 'yw'); 
mo1_ph5 = predict(mod_1_5, sun(1:N),5); 

% model order 1, prediction horizon 10
mo1_ph10 = zeros(N,1); 
mod_1_10 = ar(sun(1:N), 1, 'yw'); 
mo1_ph10 = predict(mod_1_10, sun(1:N),10); 

% model order 2, prediction horizon 1
mo2_ph1 = zeros(N,1); 
mod_2_1 = ar(sun(1:N), 2, 'yw'); 
mo2_ph1 = predict(mod_2_1, sun(1:N),1); 

% model order 2, prediction horizon 2
mo2_ph2 = zeros(N,1); 
mod_2_2 = ar(sun(1:N), 2, 'yw'); 
mo2_ph2 = predict(mod_2_1, sun(1:N),2); 

% model order 2, prediction horizon 5
mo2_ph5 = zeros(N,1); 
mod_2_5 = ar(sun(1:N), 2, 'yw'); 
mo2_ph5 = predict(mod_2_5, sun(1:N),5); 

% model order 2, prediction horizon 10
mo2_ph10 = zeros(N,1); 
mod_2_10 = ar(sun(1:N), 2, 'yw'); 
mo2_ph10 = predict(mod_2_10, sun(1:N),10); 

% model order 10, prediction horizon 1
mo10_ph1 = zeros(N,1); 
mod_10_1 = ar(sun(1:N), 10, 'yw'); 
mo10_ph1 = predict(mod_10_1, sun(1:N),1); 

% model order 10, prediction horizon 2
mo10_ph2 = zeros(N,1); 
mod_10_2 = ar(sun(1:N), 10, 'yw'); 
mo10_ph2 = predict(mod_10_2, sun(1:N),2);

% model order 10, prediction horizon 5
mo10_ph5 = zeros(N,1); 
mod_10_5 = ar(sun(1:N), 10, 'yw'); 
mo10_ph5 = predict(mod_10_5, sun(1:N),5);

% model order 10, prediction horizon 10
mo10_ph10 = zeros(N,1); 
mod_10_10 = ar(sun(1:N), 10, 'yw'); 
mo10_ph10 = predict(mod_10_10, sun(1:N),10);

figure('Position',[10,10,800,400]);
hold on;
plot(n+1700,sun(1:50),'DisplayName','Original data');
plot(n+1700,mo10_ph10,'DisplayName','Predicted data m=10,p=10');
grid on;
legend;
xlabel('year n');
ylabel('Wolf number');
title('True and predicted sunspot data');

%% 2.4 - CRLB

load('data/NASDAQ.mat')
price = NASDAQ.Close;
date = NASDAQ.Date;
plot(date,price);

p = 10;
price_norm = zscore(price);
N = length(price_norm(:,1));
x = zeros(N,10);

for i = 1:p
    ar = aryule(price_norm(1:N), i);
    x(:,i) = filter(-1*ar(1:end), 1, price_norm);
end

err=zeros(N,p);
o = ones(1,p);
err(1,:) = (x(1,:)-o*price_norm(1)).^2;

for i=2:N 
    err(i,:) = (x(i,:)-o*price_norm(i)).^2 + err(i-1,:);
end

mdl = log(err(N,(1:p))) + (1:p)*log(N)/N;
aic = log(err(N,(1:p))) + (1:p)*2/N;
aicc = zeros(1,p);

for i=1:p
    aicc(i) = aic(i) + (2*(i))*((i)+1)/(N-(i)-1);
end

figure('Position',[10,10,800,300]);
hold on;
plot(aic,'DisplayName','AIC');
plot(mdl,'DisplayName','MDL');
plot(aicc,'DisplayName','AIC_{c}')
grid on;
xlabel('Model order p');
ylabel('Magnitude');
title(['MDL and AIC and corrected AIC for model orders p = 0 to p = ' num2str(p)]);
legend show;

pcf_length = 10;
pcf = zeros(pcf_length,1);
pcf_norm = zeros(pcf_length,1);

for i = 1:pcf_length
    coefs = aryule(price, i);
    pcf(i) = -coefs(i+1);
end

for i = 1:pcf_length
    coefs_norm = aryule(price_norm, i);
    pcf_norm(i) = -coefs_norm(i+1);
end

x = linspace(1, pcf_length, pcf_length);
figure('Position',[10,10,800,300]);
stem(x,pcf,'DisplayName','PACF');
hold on;
stem(x,pcf_norm,'DisplayName','Standardised PACF');
grid on;
xlabel('k');
ylabel('\alpha (k)');
title('PACF of NASDAQ time series');
legend('show');

%% Heatmaps

n = 51:50:1001;
var = 51:50:1001;
acf_price = xcorr(price_norm,'unbiased');

[N,variance] = meshgrid(n,var);
crlb_var = 2*(variance.^2)./N;
crlb_a_1 = (variance)./(N*acf_price(924));

figure(1)
h = heatmap(n, var, crlb_var);
h.Colormap = parula;
h.ColorScaling = 'log';
title('Heatmap of noise variance')
xlabel('Number of data points N')
ylabel('True variance \sigma^{2}')

figure(2)
h = heatmap(n, var, crlb_a_1);
h.Colormap = parula;
h.ColorScaling = 'log';
title('Heatmap of a_{1}');
xlabel('Number of data points N')
ylabel('True variance \sigma^{2}')
