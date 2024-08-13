function Bhattdist = CalcBhattacharyya(pd1, pd2)
% calculate Bhattacharyya distance for multivariate gaussian distributions,
% namely pd1, pd2, estimated ProbDistributionEst(). 
% Each containing the estimated mean, covariance, and log determinant of 
% the covariance matrix.

% Reference:
% T. Kailath, "The Divergence and Bhattacharyya Distance Measures in 
% Signal Selection," in IEEE Transactions on Communication Technology, 
% vol. 15, no. 1, pp. 52-60, February 196

% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

sigma = (pd1.sigma + pd2.sigma)/2;
    [~,S,~] = svd(sigma);
    logdetsigma = sum(log(diag(S)));
    term1 = (pd2.mu - pd1.mu)'* inv(sigma) * (pd2.mu - pd1.mu)/8;
    term2 = (logdetsigma - pd1.logdetsigma/2 - pd2.logdetsigma/2)/2;
    Bhattdist = term1 + term2;
end