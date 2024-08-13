function Jdiv = CalcJKL(pd1, pd2)
% calculate symmetric KL for multivariate gaussian distributions,
% namely pd1, pd2, estimated ProbDistributionEst(). 
% Each containing the estimated mean, covariance, and log determinant of 
% the covariance matrix.

% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

KLdiv1 = CalcKL(pd1, pd2);
KLdiv2 = CalcKL(pd2, pd1);
% Jeffrey's (symmetric)
Jdiv = (KLdiv1 + KLdiv2)/2;
end