%% calculate tunings and p-values
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
    tune = RegressPerChBootstrap(X2,Y2);
    tunings = ConcatStruct(tunings,tune,2);
end
tunings.b = reshape(tunings.b, [nfeats,3,nday]);
tunings.bootstrap.pd = reshape(tunings.bootstrap.pd, [nfeats,1000,nday]);
tunings.bootstrap.md = reshape(tunings.bootstrap.md, [nfeats,1000,nday]);
% bin average feature values
% pvalues_sigma = zeros(nfeats,nday);
% pvalues_beta = zeros(nfeats,nday);
clear binkr binavg
for day = 1:nday
    disp(day)
    inds = alltrialss(event.sessionNumberPerTrial(~event.excludeTrials)==day,1);
    inds = RowColon([inds + tstart,inds + tend]);
    X2 = NDzc(inds,:);
    Y2 = labels(inds,:);
%     [final_pvalue, p_matrix, F0_vec, F1_vec, b] = Bonferroni_pvalue(X1,Y1,X2,Y2);
%     pvalues_sigma(:,day) = p_matrix(:,1);pvalues_beta(:,day) = p_matrix(:,2);
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
save('tunings_NDzc.mat','tunings','xticksday','tstart','tend')
