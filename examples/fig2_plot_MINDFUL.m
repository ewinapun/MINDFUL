% fig 1c, 2a-c. MINDFUL with KLD
clear k params
close all

params.excludeNonTrials = 1; % change to true to exclude outlier trials
params.winlen = 3000;
params.updateHz = 500;
params.xlabelFromDay0 = 1;
params.pcaDim = 5;
% params.dmetric = {'KLdiv'};

% instantiate the class
k = MINDFUL(NDzc, event, info, params, extra);

% get a lagged version of Xhat
Xhat = extra.cursorVel;
lag = 1; 
lagXhat = k.getlagX(Xhat, k.event.blockStartStop, lag);  
Xhatnlag = [Xhat lagXhat];

[refData, ind] = k.GetSessionData(init_day);
header = ['session ',num2str(init_day)];

% comment out if not selecting time steps
subsampleAE_max = 4;
ind = ind(k.extra.angleError(ind) < subsampleAE_max);
refData = k.data(ind,:);
header = [header,' AE < ', num2str(subsampleAE_max)];

% calculate PCA given refrence data
k.CalculatePCA(refData, ind, header);

% set reference as refrence data
k.SetReference(refData, ind, header);

%%

% calculate MINDFUL given NF only as inputs
k.CalcDistanceFromData('NF');

% calculate MINDFUL with NF + Xhat
k.CalcDistanceFromData('NFXhat', Xhat(k.params.indSelected,:));

% calculate MINDFUL with NF + Xhat + Xhat_lag
k.CalcDistanceFromData('NFXhatlag', Xhatnlag(k.params.indSelected,:));

% no ND just kinematics data Xhat
k.pcaInfo = [];
k.data = Xhat(k.params.indSelected,:);
k.SetReference(k.data(ind,:), ind, header);

% calculate MINDFUL with Xhat only
k.CalcDistanceFromData('Xhat')

% calculate MINDFUL with Xhat + Xhat_lag
k.CalcDistanceFromData('Xhatlag', lagXhat(k.params.indSelected,:))

% specify plot only KLD
k.PlotCompFnvsAE(1,'KLdiv', {'NF','Xhatlag','Xhat','NFXhatlag'}, ...
          {'NF','$$\hat{X} + \hat{X}_{lag}$$','$$\hat{X}$$', ...
          ['NF', '$$ +\hat{X} + \hat{X}_{lag}$$']});
if strcmp(info.participant,'T11')
    subplot(2,2,2);ylim([0,2.25]);subplot(2,2,3);ylim([0,2.25])
end
if saveGenFigure;savepdf(gcf,'KL_vs_AE');end
fprintf("\n")

% save for later
NFXhatlag = k.dist.NFXhatlag;

%% additional info for the table in supplemental
k.PlotCompFnvsAE(1,'KLdiv', ...
    {'NFXhatlag', 'NF', 'NFXhat', 'Xhat', 'Xhatlag'     }, ...
    {['low-D NF', '$$ +\hat{X} + \hat{X}_{lag}$$'], ...
    'low-D NF',['low-D NF', '$$ +\hat{X}$$'], ...
    '$$\hat{X}$$','$$\hat{X} + \hat{X}_{lag}$$'});
% ND_pca = PCAtransform(k, NDzc);
% NDzc_lag = getlagX(ND_pca, k.event.blockStartStop, lag);
% k.CalcDistanceFromData('NFnlag', NDzc_lag(k.params.indSelected,:));

%% plot KLD against outlier for T11 (supplemental)
if all(ismember(info.participant,'T11')) && ~params.excludeNonTrials
    k.PlotCompFnvsOutlier(1,'KLdiv', ...
        {'NF','Xhat','Xhatlag','NFXhatlag'}, ...
        {'NF','$$\hat{X}$$','$$\hat{X} + \hat{X}_{lag}$$', ...
          ['NF', '$$ +\hat{X} + \hat{X}_{lag}$$']});
    if saveGenFigure;savepdf(gcf,'KL_vs_outlier');end
end

%% supplemental hypothesis testing

k.CalcPairwisePvalue()
k.PlotPairwisePval()
if saveGenFigure;savepdf(gcf,'HT');end