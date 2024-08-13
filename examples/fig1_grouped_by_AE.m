%% Fig 1b. KL Grouped ND by angle error
% Code released with manuscript: Pun et al., "Measuring instability in 
% multi-day human intracortical neural recordings towards stable, 
% long-term brain-computer interfaces".
%
% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

%
clear k params
params.excludeNonTrials = 1;
k = MINDFUL(NDzc, event, info, params, extra);
ae = k.extra.angleError;
ref_ind = find(ae < 4); 
refData = k.data(ref_ind,:); 
header = 'AE < 4';
k.CalculatePCA(refData, ref_ind, header);

sTo = [];
AErange = (0:4:180-4);
AErange = [AErange;AErange+4]';
data_after_pca = k.PCAtransform(k.data);
ref_pd = ProbDistributionEst(data_after_pca(ref_ind,:));
statusprint = fprintf('Running 0 of %d', length(AErange));
for j = 1:length(AErange)-1
% group by even AE interval
    fprintf(repmat('\b',1,statusprint))
    statusprint = fprintf('Running %d of %d\n', j, length(AErange));
    ind = find(ae>=AErange(j,1) & ae<AErange(j,2));
    pd = ProbDistributionEst(data_after_pca(ind,:));
    st = k.CalcDistance(ref_pd, pd);
    st.mAE = mean(ae(ind));
    sTo = ConcatStruct(sTo, st);
end
%% plot KL vs AE
if strcmp(info.participant,'T11')
    color = [0.3725    0.3137    0.6275];
elseif strcmp(info.participant,'T5')
    color = [0.3373    0.7098    0.6627];
end
set(0,'defaultAxesFontSize',20)
figure('Units','inches','Position',[5 5 3 3]);
scatter(sTo.mAE, sTo.KLdiv, 20,'k','filled',...
        "MarkerEdgeColor",color, ...
        "MarkerFaceColor",color)
xlabel('Angle Error'); ylabel('KL divergence')
axis tight; axis square
ylim([0 .333])
yticks([0 .1 .2 .3])
xlim([0 180])
xticks([0 90 180])
xticklabels({['0',char(176)], ['90',char(176)], ['180',char(176)]})
hold on
dlm = fitlm(sTo.mAE,sTo.KLdiv,'Intercept',false);

[r,r_pval] = corr(sTo.mAE,sTo.KLdiv,'type','Pearson','rows','complete');
fprintf('Pearson correlation of mean KL to mAE = %.3f, p = %.1e.\n', r, r_pval)

limAxis = double(axis);
x = linspace(limAxis(1),limAxis(2),20)';
y = double(x*dlm.Coefficients.Estimate);
plot(x, y,'k--','linewidth',2)
text(20, 0.27, ...
    ['r=',num2str(r,2)], 'FontSize',20,'Color',color)
if saveGenFigure;savepdf(gcf,'KLsortedAE');end