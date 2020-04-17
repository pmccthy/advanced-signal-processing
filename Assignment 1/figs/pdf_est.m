%% 1.3 - Estimation of probability distributions
%% Ex 1
N = 1000;
v = randn(1,N);
pdf_est(v);

% pdf estimator function
function [f,x] = pdf_est(samples)
    [f,x] = hist(samples,100);
    plot(x,f);
end
    