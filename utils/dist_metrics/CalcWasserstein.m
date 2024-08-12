function [WDist, meanM, WcovM] = CalcWasserstein(pd1, pd2)
% calculate Wasserstein distance for multivariate gaussian
% distributions

% What is an estimate of the Wasserstein distance between the distributions
% of X and Y using the code linked in the article Random  Matrix-Improved  
% Estimation  of  the  Wasserstein  Distance between  two  Centered  Gaussian  
% Distributions (Malik TIOMOKO & Romain Couillet)
% https://github.com/maliktiomoko/RMTWasserstein
%
% Copyright Matt Harrison, Tsam Kiu Pun, 2024. Brown University
% -------------------------------------------------------------------------

lambda = sort(eig(pd1.sigma * pd2.sigma));
zeta = sort(real(eig(diag(lambda)-(1/pd2.n)*sqrt(lambda)*sqrt(lambda)')));
est = (1/pd1.p)*trace(pd1.sigma+pd2.sigma) ...
     - 2*(sum(sqrt(lambda))-sum(sqrt(zeta)))*(2*pd1.n/pd1.p);
if est < 0
 est = (1/pd1.p)*(trace(pd1.sigma)+trace(pd2.sigma)...
 -2*trace((pd1.sigma^(1/2)*pd2.sigma*pd1.sigma^(1/2))^(1/2)));
end
% include the means
% remove the 1/p scaling
% take the square root
% any imaginary parts are the result of numerical error
m2 = sum((pd1.mu - pd2.mu).^2);
% mean norm diffences
meanM = real(sqrt(m2));
WcovM = real(sqrt(pd1.p*est));
WDist = real(sqrt(pd1.p*est + m2));

end