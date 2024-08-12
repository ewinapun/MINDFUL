%% calculate tunings and p-values

% compute the cosine tunings and binned feature responses per session

% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

isBootstrap = 1;
tstart = -5; % 100 ms before go cue
tend = 50; % 1 second from go cue
alltrialss = event.trialStartStop(~event.excludeTrials,:);
inds = alltrialss(event.sessionNumberPerTrial(~event.excludeTrials)==1,1);
inds = RowColon([inds + tstart,inds + tend]);
X1 = NDzc(inds,:);
Y1 = labels(inds,:);
tunings = [];
for day = 1:nday
    disp(day)
    inds = alltrialss(event.sessionNumberPerTrial(~event.excludeTrials)==day,1);
    inds = RowColon([inds + tstart,inds + tend]);
    X2 = NDzc(inds,:);
    Y2 = labels(inds,:);
    if isBootstrap
        tune = RegressPerChBootstrap(X2,Y2);
    else 
        tune = RegressPerCh(X2,Y2);
    end
    tunings = ConcatStruct(tunings,tune,2);
end
tunings.b = reshape(tunings.b, [nfeats,3,nday]);

if isBootstrap
    tunings.bootstrap.pd = reshape(tunings.bootstrap.pd, [nfeats,1000,nday]);
    tunings.bootstrap.md = reshape(tunings.bootstrap.md, [nfeats,1000,nday]);
end

%% bin average feature values with or without Nadaraya-Watson kernel regression
clear binkr binavg
pvalues_sigma = zeros(nfeats,nday);
pvalues_beta = zeros(nfeats,nday);

for day = 1:nday
    disp(day)
    inds = alltrialss(event.sessionNumberPerTrial(~event.excludeTrials)==day,1);
    inds = RowColon([inds + tstart,inds + tend]);
    X2 = NDzc(inds,:);
    Y2 = labels(inds,:);
% uncomment for getting p-values
%     [final_pvalue, p_matrix, F0_vec, F1_vec, b] = Bonferroni_pvalue(X1,Y1,X2,Y2);
%     pvalues_sigma(:,day) = p_matrix(:,1);
%     pvalues_beta(:,day) = p_matrix(:,2);
    [xi, means, ci] = getAvgBinFeatures(X2, Y2,'kernelSmooth',true);
    binkr.means(:,:,day) = means;
    binkr.CIs(:,:,day) = ci;
    binkr.xi = xi;
    [xi, means, ci] = getAvgBinFeatures(X2, Y2,'kernelSmooth',false,'binwidth',45);
    binavg.means(:,:,day) = means;
    binavg.CIs(:,:,day) = ci;
    binavg.xi = xi;
end
tunings.binkr = binkr;
tunings.binavg = binavg;
% save('tunings_NDzc.mat','tunings','xticksday','tstart','tend')

%%
function [xi, means, ci] = getAvgBinFeatures(ND, labels, options)
    arguments
        ND double
        labels double
        options.binwidth = 10 % bin width in degree
        options.bw = 22.5 % bandwidth for regression in degree
        options.kernelSmooth = true
        options.alpha = 0.05; % confidence interval
    end
    nfeats = size(ND,2);
    angles = atan2d(labels(:,2), labels(:,1));
    xi = (-180:options.binwidth:180)'; % range of the bins
    nbins = length(xi);
    means = zeros(nbins, nfeats); % mean
    ci = zeros(nbins, nfeats); % confidence intervals
        
    if options.kernelSmooth
        % repeat in range [x - 2pi, x, x + 2pi] for accounting edge effect
        disp('running Nadaraya-Watson kernel regression')
        x = vertcat(angles-360, angles, angles+360);

        for f = 1:nfeats
            disp(['feature #', num2str(f)])
            y = repmat(ND(:, f),[3,1]);
            [r, ci_diff, xi, bw] = KernelRegression(x, y, xi, options.bw, options.alpha);
            means(:,f) = r;
            ci(:,f) = ci_diff;
% bootstrap to get confidence intervals (not used)
%             r_all = zeros(100, nbins);
%             for b = 1:100
%                 perm = randi(length(x),1,length(x));
%                 r_all(b,:) = KernelRegression(x(perm),y(perm),xi,'bandwidth',options.bw);
%             end
%             means(:,f) = mean(r_all,1);
%             ci(:,f) = 1.96 * std(r_all);
        end
    else
        % Binned neural data
        [cnt, binInds, ~] = PolarHist(angles, nbins, 1);
        err_t = zeros(nbins, 2);
        % Estimate mean
        for f = 1:nfeats
            for t = 1:nbins
                [means(t,f), ~, err_t(t,:)] = normfit(ND(binInds==t, f));
            end
            ci(:,f) = means(:,f) - err_t(:,1); % symmetric error bar
        end
        % shift display order from (0, 360) to (-180, +180)
        means = means([5:8,1:5],:);
        ci = ci([5:8,1:5],:);
    end
end