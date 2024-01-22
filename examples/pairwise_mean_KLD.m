%% calculate pairwise KL divergence
params.xlabelFromDay0 = 1;
params.excludeNonTrials = 1;
params.updateHz = 500;
k = MINDFUL(NDzc, event, info, params, extra);

% get a lagged version of Xhat
Xhat = extra.cursorVel;
lag = 1; 
lagXhat = getlagX(Xhat, k.event.blockStartStop, lag);  
Xhatnlag = [Xhat lagXhat];

% setting PCA subspace
subsampleAE_max = 4;
[~, ind] = k.GetSessionData(init_day);
ind = ind(k.extra.angleError(ind) < subsampleAE_max);
refData = k.data(ind,:);
header = ['session ',num2str(init_day),' AE < ', num2str(subsampleAE_max)];

k.CalculatePCA(refData, ind, header);

% uses paarallel pool (parfor) to speed up computation
k.CalcPairwiseDistanceParallel('ND',Xhatnlag(k.params.indSelected,:))

% Supplemental plot pairwise KLD
k.PlotDistance()

% calculate average pairwise KL divergence between sessions
pwmeanKL = NaN(nday,nday);
for day_i = 1:nday
    segInds_i = k.GetSessionSegmentInds(day_i);
    for day_j = 1:nday
        segInds_j = k.GetSessionSegmentInds(day_j);
        pwmeanKL(day_i,day_j) = ...
        mean(k.dist.ND.KLdiv(segInds_i,segInds_j),'all','omitnan');
    end
end
disp('pairwise KLD: ')
disp(pwmeanKL)