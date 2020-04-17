%% 1.3 - Estimation of probability distributions
%% Ex 1
N = 1000;
v = randn(1,N);
pdf_est(v);

% pdf estimator function
function pdf_est(samples)
    histogram(samples,'BinWidth',0.05,'Normalization','pdf');grid on;
    xlabel('Random Variable');
    ylabel('Probability');
    title('PDF Estimate for RV');
end
    