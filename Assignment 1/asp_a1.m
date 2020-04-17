% Advanced Signal Processing 
% Coursework Assignment 1 - Random Signals and Stochastic Processes
%% 1.1 - Statistical Estimation
%% a. UNIFORM DISTRIBUTIONS 
%% Ex 1a 
clear all; close all; clc;

x = rand(1,1000);

x = x';
figure('Position', [10 10 800 200]);
plot(x);xlabel('n');
ylabel('x[n]');
title('Realisation of uniform RV');
grid on;

% theoretical mean of uniform distribution between 0 and 1 would be 0.5
theo_mean = 0.5;
samp_mean = mean(x);
err_mean = theo_mean-samp_mean;
fprintf('Error: %1.2f\n',err_mean)

%% Ex 2a
% theoretical SD of uniform distribution between 0 and 1 would be
% 1/sqrt(12)
theo_sd = 1/sqrt(12);
samp_sd = std(x);
err_sd = theo_sd-samp_sd
fprintf('Error: %1.2f\n',err_sd)

%% Ex 3a
clear all; close all; clc;
theo_mean = 0.5;
theo_sd = 1/sqrt(12);
for i = 1:10
    x = rand(1,1000); 
    x = x';
    samp_mean(i) = mean(x);
    samp_sd(i) = std(x);
end

figure('Position',[100,100,600,400]);
plot(samp_mean,'Marker','X','LineStyle','none','LineWidth',2);
ylim([0.45,0.55]);yline(theo_mean,'r','LineWidth',2);
xlabel('Realisation number');
ylabel('Sample mean');legend('Sample mean','Theoretical mean');
title('Theoretical and sample means');
grid on;

figure('Position',[100,100,600,400]);
plot(samp_sd,'Marker','X','LineStyle','none','LineWidth',2);
ylim([0.26,0.3]);
yline(theo_sd,'r','LineWidth',2);
xlabel('Realisation number');
ylabel('Sample SD');legend('Sample SD','Theoretical SD');
title('Theoretical and sample SD');
grid on;

%% Ex 4a
clear all; close all; clc;
theo_mean = 0.5;
theo_sd = 1/sqrt(12);

x = rand(1,1000);

% still need to normalise histograms
figure;histogram(x,10,'Normalization','probability');
histfit(x);yline(0.1,'r','LineWidth',2);
ylabel('Pro');xlabel('Sample value');
legend('Sample histogram (estimated pdf)','Theoretical pdf');
title('Normalised histogram for uniform RV');
grid on;

%% b. GAUSSIAN DISTRIBUTIONS 
%% Ex 1b 
clear all; close all; clc;

x = randn(1,1000);

x = x';
figure('Position', [10 10 800 200]);
plot(x);
xlabel('n');
ylabel('x[n]');

% theoretical mean of unit SD zero-mean Gaussian distribution is 0
theo_mean = 0;
samp_mean = mean(x);
err_mean = theo_mean-samp_mean;
fprintf('Error: %1.2f',err_mean)

% Ex 2b
% theoretical SD of unit SD zero-mean Gaussian distribution is 1
theo_sd = 1
samp_sd = std(x);
err_sd = theo_sd-samp_sd;
fprintf('Error: %1.2f',err_sd)

%% Ex 3b
clear all; close all; clc;

theo_mean = 0;
theo_sd = 1;

for i = 1:10
    x = randn(1,1000); 
    x = x';
    samp_mean(i) = mean(x);
    samp_sd(i) = std(x);
end

figure('Position',[100,100,600,400]);
plot(samp_mean,'Marker','X','LineStyle','none','LineWidth',2);
ylim([-0.1,0.1]);yline(theo_mean,'r','LineWidth',2);
xlabel('Realisation number');
ylabel('Sample mean');
legend('Sample mean','Theoretical mean');
title('Theoretical and sample means');
grid on;

figure('Position',[100,100,600,400]);
plot(samp_sd,'Marker','X','LineStyle','none','LineWidth',2);
ylim([0.9,1.1]);yline(theo_sd,'r','LineWidth',2);
xlabel('Realisation number');
ylabel('Sample SD');
legend('Sample SD','Theoretical SD');
title('Theoretical and sample SD');
grid on;

%% Ex 4b
clear all; close all; clc;
theo_mean = 0;
theo_sd = 1;

x = randn(1,1000);

figure;
h = histfit(x,10);
h(1).FaceColor = [0, 0.4470, 0.7410];
grid on;
y_tick = get(gca, 'YTick');
set(gca, 'YTick', y_tick, 'YTickLabel', y_tick/numel(x));
ylabel('Probability');
xlabel('Sample value');
legend('Sample histogram (estimated pdf)','Theoretical pdf');
title('Normalised histogram for Gaussian RV');

%% 1.2 - Stochastic Processes 
%% Ex 1

clear all; close all; clc;

% create and plot random processed

% RP1
randproc1 = rp1(100,100);
mean1 = mean(randproc1);
sd1 = std(randproc1);

figure('Position', [10 10 500 200]);
plot(mean1);y
label('x_1[n]');
xlabel('Sample number n');
title('Ensemble mean for RP1');
grid on;

figure('Position', [10 10 500 200]);
plot(sd1);
ylabel('x_1[n]');xlabel('Sample number n');
title('Ensemble SD for RP1');
grid on;

% RP2
randproc2 = rp2(100,100);
mean2 = mean(randproc2);
sd2 = std(randproc2);

figure('Position', [10 10 500 200]);
plot(mean2);
ylabel('x_2[n]');
xlabel('Sample number n');
title('Ensemble mean for RP2');
grid on;
figure('Position', [10 10 500 200]);
plot(sd2);ylabel('x_2[n]');
xlabel('Sample number n');
title('Ensemble SD for RP2');
grid on;

% RP3
randproc3 = rp3(100,100);
mean3 = mean(randproc3);
sd3 = std(randproc3);

figure('Position', [10 10 500 200]);
plot(mean3);ylabel('x_3[n]');
xlabel('Sample number n');
title('Ensemble mean for RP3');
grid on;

figure('Position', [10 10 500 200]);
plot(sd3);ylabel('x_3[n]');
xlabel('Sample number n');
title('Ensemble SD for RP3');
grid on;



%% ensemble means and SDs
clc; close all;

for i = 1 : 100
	mean1(i) = mean(randproc1(:,i));
    sd1(i) = std(randproc1(:,i));
    mean2(i) = mean(randproc2(:,i));
    sd2(i) = std(randproc2(:,i));
    mean3(i) = mean(randproc3(:,i));
    sd3(i) = std(randproc3(:,i));
end

figure('Position', [10 10 400 200]);
plot(mean1(1,:));
ylabel('m_{x1}[n]');
('Sample number n');

figure('Position', [10 10 400 200]);
plot(sd1(1,:));
ylabel('\sigma_{x1}[n]');
('Sample number n');

figure('Position', [10 10 400 200]);
plot(mean2(1,:));
ylabel('m_{x2}[n]');
('Sample number n');

figure('Position', [10 10 400 200]);
plot(sd2(1,:));
ylabel('\sigma_{x1}[n]');
('Sample number n');

figure('Position', [10 10 400 200]);
plot(mean2(1,:));
ylabel('m_{x3}[n]');
('Sample number n');

figure('Position', [10 10 400 200]);
plot(sd3(1,:));
ylabel('\sigma_{x1}[n]');
('Sample number n');

%% Ex 2

randproc4 = rp1(4,1000);
randproc5 = rp2(4,1000);
randproc6 = rp3(4,1000);

% means and SDs for each realisation
for i = 1 : 4
	mean4(i) = mean(randproc4(i,:));
    sd4(i) = std(randproc4(i,:));
    mean5(i) = mean(randproc5(i,:));
    sd5(i) = std(randproc5(i,:));
    mean6(i) = mean(randproc6(i,:));
    sd6(i) = std(randproc6(i,:));
end


%% 1.3 - Estimation of probability distributions
%% Ex 1
    
%randproc1 = rp1(1,10000);
%randproc2 = rp2(1,10000);
randproc3 = rp3(1,10000);
set(gcf,'renderer','Painters')
figure;
pdf_est(randproc3);
title('PDF estimate for RP3, N=10,000')

%% function definitions

% RP1
function v = rp1(M,N);
    a = 0.02;
    b = 5;
    Mc = ones(M,1)*b*sin((1:N)*pi/N);
    Ac = a*ones(M,1)*[1:N];
    v = (rand(M,N)-0.5).*Mc + Ac;
end

% RP2
function v = rp2(M,N)
    Ar = rand(M,1)*ones(1,N);
    Mr = rand(M,1)*ones(1,N);
    v = (rand(M,N)-0.5).*Mr + Ar;
end

RP3
function v = rp3(M,N)
    a = 0.5;
    m = 3;
    v = (rand(M,N)-0.5)*m + a;
end


% pdf estimator function
function pdf_est(samples)
    hist(samples,100,'Normalization','pdf');grid on;
    xlabel('Random Variable');
    ylabel('Probability');
    title('PDF Estimate for RV');
end
    