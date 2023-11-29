function [KLdiv, extra] = CalcKL(pd1, pd2)
% calculate KL divergence for multivariate gaussian
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