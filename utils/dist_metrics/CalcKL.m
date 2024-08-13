function [KLdiv, extra] = CalcKL(pd1, pd2)
% calculate Kullback-Leibler divergence for multivariate gaussian distributions,
% namely pd1, pd2, estimated ProbDistributionEst(). 
% Each containing the estimated mean, covariance, and log determinant of 
% the covariance matrix.

% Reference:
% Kullback, S.; Leibler, R.A. (1951). "On information and sufficiency". 
% Annals of Mathematical Statistics. 22 (1): 79–86. 
% doi:10.1214/aoms/1177729694. JSTOR 2236703. MR 0039968.

% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

if rank(pd2.sigma) ~= size(pd2.sigma,1)
    KLdiv = NaN;
else
    sigma1inv = inv(pd2.sigma);
    term1 = trace(sigma1inv * pd1.sigma);
    term2 = (pd2.mu - pd1.mu)' * sigma1inv * (pd2.mu - pd1.mu);
    term3 = pd2.logdetsigma - pd1.logdetsigma;
    KLdiv = 1/2*(term1 + term2 - pd1.p + term3);
    extra.trace = term1;
    extra.meandiffsig = term2;
    extra.logdetsig = term3;
    extra.dim = pd1.p;
end

end