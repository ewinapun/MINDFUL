function tunings = RegressPerChBootstrap(ND, labels, n, alpha, relativeMD)
% this function estimate the cosine tunings of features
% inputs
%   ND
%       neural data ( timesteps x nfeats )
%   labels
%       intended directions ( timesteps x dim )
%   alpha
%       for calculating confidence level (scalar. Default to 0.05)
% 
% outputs
%   tunings
%       .pd
%           preferred directions
%       .md
%           modulation depth
%       .b
%           tuning parameters
%       .rsq
%           r-square statistic
%       .pval
%           F-statistic p-value
%       .nTimesteps
%           number of time steps in N
%       .iszscore
%           logical if ND seem zscored (inferred)

%
% History:
%   2021.02.19   created by Ewina Pun
%   2023.06.22   last edit by Ewina Pun
%
% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------


if nargin < 5
    % rough estimate if ND is zscored. 
    % If true, MD should be adjusted by baseline firing rate.
    relativeMD = abs(mean(ND,'all')) > 1; 
end

if nargin < 4
    alpha = 0.05;
end

if nargin < 3
    n = 1000; % number of bootstrap
end
   
[nTimesteps, nFeat] = size(ND);

if size(labels,2)==2
    theta = atan2(labels(:,2),labels(:,1));
else
    theta = labels;
end

% cosine tuning model
tuneModel   = [ones(length(theta),1) cos(theta) sin(theta)];

pd   = zeros(nFeat,1);
md   = zeros(nFeat,1);
rsq  = zeros(nFeat,1);
pval = zeros(nFeat,1);
b    = zeros(nFeat,3);
%     b_CI = zeros(nFeat,3);

% preferred direction
getpd = @(b) atan2d(b(:,3),b(:,2));

% modulation depth
getmd = @(b) sqrt(b(:,3).^2+b(:,2).^2);

pds = zeros(nFeat, n);
mds = zeros(nFeat, n);
tic
% Regress through all features
for f = 1:nFeat
    [b(f,:),~,~,~,stats] = regress(ND(:,f), tuneModel, alpha);
    pd(f) = getpd(b(f,:));
    md(f) = getmd(b(f,:)); 
    rsq(f)  = stats(1);
    pval(f) = stats(3);

    % get distributions with bootstrap
    opt = statset('UseParallel',true);
    bs = bootstrp(n, @regress, ND(:,f), tuneModel, alpha,'Options',opt);
    
    % estimating pd distributions
    pds(f,:) = getpd(bs)';

    % estimating md distributions
    if relativeMD
        mds(f,:) = getmd(bs)./bs(:,1); 
    else
        mds(f,:) = getmd(bs);
    end        

%     if md(f) > 0.1
%     pd_ci = prctile(pds(f,:),[100*alpha/2,100*(1-alpha/2)]);
%     md_ci = prctile(mds(f,:),[100*alpha/2,100*(1-alpha/2)]);
%     figure(1001);
%     subplot(1,2,1)
%     histogram(pds(f,:),100);
%     xline(pd(f),'r'); xline(circ_mean(pds(f,:)'*pi/180)*180/pi,'r--')
%     xline(pd_ci(1),'k--');xline(pd_ci(2),'k--');
%     
%     subplot(1,2,2)
%     histogram(mds(f,:));
%     xline(md(f),'r'); xline(mean(mds(f,:)),'r--')
%     xline(md_ci(1),'k--');xline(md_ci(2),'k--');
%     end
end
toc

%b save outputs
tunings.bootstrap.pd      = pds;
tunings.bootstrap.md      = mds;
tunings.pd      = pd;
tunings.md      = md;
tunings.b       = b;
tunings.rsq     = rsq;
tunings.pval    = pval;
tunings.nPoints = nTimesteps;
tunings.iszscore = ~relativeMD;
end