%% calculate pairwise KL divergence
params.xlabelFromDay0 = 1;
params.excludeNonTrials = 1;
params.updateHz = 500;
k = MINDFUL(NDzc, event, info, params, extra);
[~, ind] = k.GetSessionData(init_day);
header = ['session ',num2str(init_day)];

% comment out if not selecting time steps
subsampleAE_max = 4;
ind = ind(k.extra.angleError(ind) < subsampleAE_max);
refData = k.data(ind,:);
header = [header,' AE < ', num2str(subsampleAE_max)];

k.CalculatePCA(refData, ind, header);
k.CalcPairwiseDistanceParallel('ND',[Xhat(k.params.indSelected,:) lagXhat(k.params.indSelected,:)])

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