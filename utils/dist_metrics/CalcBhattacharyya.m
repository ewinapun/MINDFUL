function Bhattdist = CalcBhattacharyya(pd1, pd2)
% calculate Bhattacharyya distance for multivariate gaussian distributions
    sigma = (pd1.sigma + pd2.sigma)/2;
    [~,S,~] = svd(sigma);
    logdetsigma = sum(log(diag(S)));
    term1 = (pd2.mu - pd1.mu)'* inv(sigma) * (pd2.mu - pd1.mu)/8;
    term2 = (logdetsigma - pd1.logdetsigma/2 - pd2.logdetsigma/2)/2;
    Bhattdist = term1 + term2;
end