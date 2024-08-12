function [r, ci_diff, xi, bw] = KernelRegression(x, y, pnts, bw, alpha, varargin)
% function [r,xi,bw,CIl,CIr] = KernelRegression(X,Y,pnts,h,...)
%
% 1d Nadaraya-Watson kernel regression for r(x) = E(Y|X=x)
%
% X is a vector of 1d covariates
% Y is the corresponding vector of scalar responses
% (1-alpha) is the confidence level we need
% h is the bandwidth for kernel regression
% pnts is the values x where the density estimate r(x) is to be evaluated.

% r(i) is the estimate of the regression curve at xi(i)
%[CIl(i),CIr(i)] is the confidence interval for r(i) at point xi(i) 
% xi as the output is the same as the input pnts, and output bw is the same as
% input h
%
% The optional arguments are the same as kspdfXhatsity. The main ones of
% interest are 'bandwidth' and 'kernel'. Do not use the arguments
% 'function' or 'plotfcn'.

%   2023 Copyright Mona Khoshnevis, Brown University. All Rights Reserved
%--------------------------------------------------------------------------

% check inputs dimension
if isvector(x)
    x = x(:).';
else
    error('size is wrong should a vector.')
end

if isvector(y)
    y = y(:).';
else
    error('size is wrong should a vector.')
end

if isvector(pnts)
    pnts = pnts(:).';
else
    error('size is wrong. pnts should a vector.')
end

% Compute the pdfXhatominator as a kernel pdfXhatsity estimator of the pdf of x.
n = numel(x);
%bw = n^(-1/5);
[pdfXhat,xi,bw] = ksdensity(x,pnts,varargin{:},'bandwidth',bw,'function','pdf');
% Compute the numerator as a weighted kde of the pdf of y.
% For this to work, y must be positive, so subtract a constant and then add
% it back.
ymin = min(y); ymin = ymin-eps(ymin); ypositive = y-ymin;
r = ksdensity(x,pnts,varargin{:},'weights',ypositive,'bandwidth',bw,'function','pdf')./pdfXhat*mean(ypositive)+ymin;

% Compute the estimator of variance at each point xi(i), as a weighted kde
% with weights (y-r(i)).^2
sigmahat2 = zeros(1,numel(xi));

for i = 1:numel(xi)
    w = (y-r(i)).^2;
    pts = xi(i);
    f = ksdensity(x,pts,varargin{:},'weights',w,'bandwidth',bw,'function','pdf');
    sigmahat2(i) = f*sum(w)*(1/n)/pdfXhat(i);
end

z = norminv(1-alpha/2);

l2normKsq = 1/(2*sqrt(pi));

ci_diff = z .* sqrt( (sigmahat2*l2normKsq)./(n*bw*pdfXhat));

% CIl = r - z .* sqrt( (sigmahat2*l2normKsq)./(n*bw*pdfXhat) );
% CIr = r + z .* sqrt( (sigmahat2*l2normKsq)./(n*bw*pdfXhat) );
end