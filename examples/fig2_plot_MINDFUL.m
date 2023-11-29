% fig 1c, 2a-c. MINDFUL with KLD
clear k params

params.excludeNonTrials = 0; % change to true to exclude outlier trials
params.winlen = 3000;
params.updateHz = 50;
params.xlabelFromDay0 = 1;
params.pcaDim = 5;

% instantiate the class
k = MINDFUL(NDzc, event, info, params, extra);

% get a lagged version of Xhat
Xhat = extra.cursorVel;
lag = 1; 
lagXhat = getlagX(Xhat, k.event.blockStartStop, lag);  
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

disp(k)

% calculate MINDFUL given ND only as inputs
k.CalcDistanceFromData('ND');

% calculate MINDFUL with ND + Xhat
k.CalcDistanceFromData('NDXhat', Xhat(k.params.indSelected,:));

% calculate MINDFUL with ND + Xhat + Xhat_lag
k.CalcDistanceFromData('NDXhatlag', Xhatnlag(k.params.indSelected,:));

% no ND just kinematics data Xhat
k.pcaInfo = [];
k.data = Xhat(k.params.indSelected,:);
k.SetReference(k.data(ind,:), ind, header);

% calculate MINDFUL with Xhat only
k.CalcDistanceFromData('Xhat')

% calculate MINDFUL with Xhat + Xhat_lag
k.CalcDistanceFromData('Xhatlag', lagXhat(k.params.indSelected,:))

% specify plot only KLD
k.PlotCompFnvsAE(1,'KLdiv', {'ND','Xhat','Xhatlag','NDXhatlag'}, ...
          {'low-D ND','$$\hat{X}$$','$$\hat{X} + \hat{X}_{lag}$$', ...
          ['low-D ND', '$$ +\hat{X} + \hat{X}_{lag}$$']});

% save for later
NDXhatlag = k.dist.NDXhatlag;

% additional info for the table in supplemental
k.PlotCompFnvsAE(1,'KLdiv', ...
    {'NDXhatlag', 'ND', 'NDXhat', 'Xhat', 'Xhatlag'     }, ...
    {['low-D ND', '$$ +\hat{X} + \hat{X}_{lag}$$'], ...
    'low-D ND',['low-D ND', '$$ +\hat{X}$$'], ...
    '$$\hat{X}$$','$$\hat{X} + \hat{X}_{lag}$$'});
% ND_pca = PCAtransform(k, NDzc);
% NDzc_lag = getlagX(ND_pca, k.event.blockStartStop, lag);
% k.CalcDistanceFromData('NDnlag', NDzc_lag(k.params.indSelected,:));

% plot KLD against outlier for T11 (supplemental)
if all(ismember(info.participant,'T11')) && ~params.excludeNonTrials
k.PlotCompFnvsOutlier(1,'KLdiv', {'ND','Xhat','Xhatlag','NDXhatlag'}, ...
          {'low-D ND','$$\hat{X}$$','$$\hat{X} + \hat{X}_{lag}$$', ...
          ['low-D ND', '$$ +\hat{X} + \hat{X}_{lag}$$']});
end
