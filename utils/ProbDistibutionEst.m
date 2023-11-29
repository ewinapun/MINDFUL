function pd = ProbDistibutionEst(X)
% X is a n x d matrix of n iid d-dim multivariate normal random vectors
% estimate mean and covariance of data X
%             if size(X, 2) > size(X, 1)
%                 X = X';
%             end
    % check data type of X:
    if ~isa(X,'double')
       if isa(X,'single')
            X = double(X);
       else
           error('X should be a floating-point array')
       end
    end
    pd.n = size(X,1);
    pd.p = size(X,2);
    pd.mu = mean(X,1)';
    pd.sigma = X'*X / pd.n;
%             use gmdist to dit a Gaussian mixture distribution to data
%             pgm = fitgmdist(X,1,'Regularize',0.0001);
%             pd.mu = pgm.mu';
%             pd.sigma = double(pgm.Sigma);
    [~,S,~] = svd(pd.sigma);
    pd.logdetsigma = sum(log(diag(S)));
    if rank(pd.sigma) ~= pd.p
        disp('problem with rank')
    end
end