%% Advanced Signal Processing
%% Coursework Assignment 4 - Fixed and adaptive optimal filtering
%% 4.1 Wiener filter

x = randn(1000,1);

b = [1,2,3,2,1];
a = 1;
y = filter(b,a,x);
y = zscore(y);
noise = randn(1,1000);
noise = noise';
z = y + noise;

%% 4.2 - LMS algorithm

order = 5;
mu = 0.5;
[y_est,err,coeffs] = lms(x,z,mu,order);

%% signal adaptation

figure;
plot(y);
hold on;
plot(y_est);
grid on;
title('Actual and LMS-estimated signal');
xlabel('n');
ylabel('y[n]');
legend('Actual system output','Estimate using LMS')

%% filter coefficients

figure('Position',[250,250,800,300]);
for i = 1:5
    hold on;
    plot(coeffs(i,:));
end
grid on;
title("Filter coefficients with adaptation gain \mu = "+mu);
xlabel('n');
ylabel('Coefficient values w[i]');
legend('w[1]','w[2]','w[3]','w[4]','w[5]');

% order determination - use MDL and AIC from section 3

%% 4.4 - Identification of AR processes

% AR(2) coefficient estimation

x2 = 1*randn(10000,1);
b2 = [1];
a2 = [1 0.9 0.2];
y2 = filter(b2,a2,x2);
y2_del = [0; y2(1:10000-1)];
ord2 = 2;
mu2 = 0.005;
[y2_est,err2,coeffs2] = lms(y2_del,y2,mu2,ord2);

% plot coefficients
figure('Position',[250,250,800,300]);
for i = 1:2
    hold on;
    plot(coeffs2(i,:));
end
grid on;
title("Filter coefficients with adaptation gain \mu = "+mu2);
xlabel('n');
ylabel('Coefficient values');
legend('a_{1}','a_{2}');

%% 4.5 - Speech Recognition

% load audio, resample and isolate 1000-sample signal

Fs = 44100; % frequency to be resampled at

% resample audio and isolate appropriate 1000 samples

[E,fs]=audioread('audio/E.m4a');
E = resample(E,44100,fs);
E = E(26000:27000);
E_shifted = [0;E(1:999)];

[A,fs]=audioread('audio/A.m4a');
A = resample(A,44100,fs);
A = A(22000:23000);
A_shifted = [0;A(1:999)];

[S,fs]=audioread('audio/S.m4a');
S = resample(S,44100,fs);
S = S(41000:42000);
S_shifted = [0;S(1:999)];

[T,fs]=audioread('audio/T.m4a');
T = resample(T,44100,fs);
T = T(23000:24000);
T_shifted = [0;T(1:999)];

[X,fs]=audioread('audio/X.m4a');
X = resample(X,44100,fs);
X = X(36000:37000);
X_shifted = [0;X(1:999)];


%% apply LMS routine 

% set order and gain based on MDL, AIC and AICC above
gain = 0.25;
order_E = 3;
order_A = 2;
order_S = 11;
order_T = 6;
order_X = 12;

[E_pred,E_err,E_coeffs] = lms(E_shifted, E, gain, order_E); 
[A_pred,A_err,A_coeffs] = lms(A_shifted, A, gain, order_A); 
[S_pred,S_err,S_coeffs] = lms(S_shifted, S, gain, order_S); 
[T_pred,T_err,T_coeffs] = lms(T_shifted, T, gain, order_T); 
[X_pred,X_err,X_coeffs] = lms(X_shifted, X, gain, order_X); 

figure('Position',[250,250,800,400]);
plot(1:999,E_coeffs);
title("Evolution of Wiener coefficients for letter E with adaptation gain \mu = "+gain);
xlabel('n');
ylabel('w[n]');
legend('w[1]','w[2]','w[3]');
grid on;

figure('Position',[250,250,800,400]);
plot(1:999,A_coeffs);
title("Evolution of Wiener coefficients for letter A with adaptation gain \mu = "+gain);
xlabel('n');
ylabel('w[n]');
legend('w[1]','w[2]');
grid on;

figure('Position',[250,250,800,400]);
plot(1:999,S_coeffs);
title("Evolution of Wiener coefficients for letter S with adaptation gain \mu = "+gain);
xlabel('n');
ylabel('w[n]');
legend('w[1]','w[2]','w[3]','w[4]','w[5]','w[6]','w[7]','w[8]','w[9]','w[10]','w[11]');
grid on;

figure('Position',[250,250,800,400]);
plot(1:999,T_coeffs);
title("Evolution of Wiener coefficients for letter T with adaptation gain \mu = "+gain);xlabel('n');
ylabel('w[n]');
legend('w[1]','w[2]','w[3]','w[4]','w[5]','w[6]');
grid on;

figure('Position',[250,250,800,400]);
plot(1:999,X_coeffs);
title("Evolution of Wiener coefficients for letter X with adaptation gain \mu = "+gain);
xlabel('n');
ylabel('w[n]');
legend('w[1]','w[2]','w[3]','w[4]','w[5]','w[6]','w[7]','w[8]','w[9]','w[10]','w[11]','w[12]');
grid on;

%% determining model order

letter = A_shifted;
p = 50;
letter_norm = zscore(letter);
N = length(letter_norm(:,1));
x = zeros(N,p);

for i = 1:p
    ar = aryule(letter_norm(1:N), i);
    x(:,i) = filter(-1*ar(1:end), 1 ,letter_norm(:,1));
end

D = zeros(N,p);
order = ones(1,p);
D(1,:) = (x(1,:)-order*letter_norm(1)).^2;

for i=2:N 
    D(i,:) = (x(i,:)-order*letter_norm(i)).^2 + D(i-1,:);
end

mdl = log(D(N,(1:p))) + (1:p)*log(N)/N;
aic = log(D(N,(1:p))) + (1:p)*2/N;

figure('Position',[10,10,800,300]);
hold on;
plot(aic,'DisplayName','AIC');
plot(mdl,'DisplayName','MDL');
grid on;
xlabel('Model order p');
ylabel('Magnitude');
title(['MDL and AIC for letter model orders p = 0 to p = ' num2str(p)]);
legend('show');

%% testing predictor

x2 = X;
b2 = [1];
a2 = [1,0.9,0.2];
y2 = filter(b2,a2,x2);
y2_del = [0; y2(1:1000-1)];
ord2 = 2;
mu2 = 0.05;
[y2_est,err2,coeffs2] = lms(y2_del,y2,mu2,ord2);

% plot coefficients
figure('Position',[250,250,800,300]);
for i = 1:2
    hold on;
    plot(coeffs2(i,:));
end
grid on;
title("Filter coefficients with adaptation gain \mu = "+mu2);
xlabel('n');
ylabel('Coefficient values');
legend('a_{1}','a_{2}');



%% function definitions

% LMS routine

function [y_estimate,error,coeffs] = lms(x,z,mu,order)
    N = length(x); % signal length
    coeffs = zeros(order, N-1); % each row represents coefficient for 1 variable (order of filter)
    y_estimate = zeros(N, 1); % must be column vector
    error = zeros(N, 1); % must be column vector
    for i = order+1:N % starts from order+1 because estimate requires number of inputs equal to model order
        a = coeffs(:,i-order); % AR coefficient
        b = x(i:-1:i-order+1); % MA coefficient (remember x is a column vector)
        y_estimate(i) = a'*b;
        error(i) = z(i) - y_estimate(i); % find error between real output z and model-generated estimate
        coeffs(:,i-order+1) = coeffs(:,i-order)+mu*error(i)*b; % fills coefficients matrix with correct coefficients, advancing column each iteration of the for loop
    end
end

% LMS routine with gear shifting

function [y_estimate,error,coeffs] = lms_gear_shifting(x,z,mu,order)
    N = length(x); % signal length
    coeffs = zeros(order, N-1); % each row represents coefficient for 1 variable (order of filter)
    y_estimate = zeros(N, 1); % must be column vector
    error = zeros(N, 1); % must be column vector
    for i = order+1:N % starts from order+1 because estimate requires number of inputs equal to model order
        a = coeffs(:,i-order); % AR coefficient
        b = x(i:-1:i-order+1); % MA coefficient (remember x is a column vector)
        y_estimate(i) = a'*b;
        error(i) = z(i) - y_estimate(i); % find error between real output z and model-generated estimate
        if error(i)>error(i-1)
            mu = 1.1*mu;
        else if error(i)<error(i-1)
                mu = 0.9*mu;
            end
        end
        coeffs(:,i-order+1) = coeffs(:,i-order)+mu*error(i)*b; % fills coefficients matrix with correct coefficients, advancing column each iteration of the for loop
    end
end


