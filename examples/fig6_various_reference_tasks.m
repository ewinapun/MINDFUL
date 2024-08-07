%% Fig 6 Different tasks for reference

clear k params

% define data paths
rand_target_datapath = [dirpath, filesep, 'T11(additional)', ...
                                 filesep, 'random_targets'];
personal_use_datapath = [dirpath, filesep, 'T11(additional)', ...
                                 filesep, 'personal_use'];

params.excludeNonTrials = 0; % change to true to exclude outlier trials
params.winlen = 3000;
params.updateHz = 500;
params.xlabelFromDay0 = 1;
params.pcaDim = 5;

% instantiate the class
k = MINDFUL(NDzc, event, info, params, extra);

% get a lagged version of Xhat
Xhat = extra.cursorVel;
lag = 1;
lagXhat = k.getlagX(Xhat, k.event.blockStartStop, lag);
Xhatnlag = [Xhat lagXhat];


% Loop between different cases
ref_block_name_str = {'center out', 'fitts', 'personal use', 'mixed tasks'};

p.zscoreFeatures = true; % set false for raw extract features
p.trailingMoving = true; % set true for causal rolling z-score using movmean; set false for block batch z-score

for ii = 1:4
% change which reference block to use
ref_block = ref_block_name_str{ii};
fprintf('loading "%s" as the reference ...\n', ref_block)

switch ref_block
    case 'center out'
        [refData, ind] = k.GetSessionData(1);
%         refData = k.data;
        refXhat = Xhat;
        reflagXhat = lagXhat;
        refXhatnlag = Xhatnlag;
        co_ref = Xhat(ind,:);
        header = 'center out day 0';

    case 'fitts'
        % random target, also named the 'fitts' task
        path = rand_target_datapath;
        [refData, ~, refevent, ~, refextra] = ConcatSavedSessionsData(path, p);
        refXhat = refextra.cursorVel;
        reflagXhat = k.getlagX(refXhat, refevent.blockStartStop, lag);
        refXhatnlag = [refXhat reflagXhat];
        ind = (1:length(refXhat));
        fitt_ref = refXhat(ind,:);
        header = 'Random target on day 7';

    case 'personal use'
        % personal use
        path = personal_use_datapath;
        [refData, ~, refevent, ~, refextra] = ConcatSavedSessionsData(path, p);
        % here trialStartStop refers to start stops of active cursor
        % control, no explicit trials were cued in this dataset
        ind = RowColon(refevent.trialStartStop(1:end,:));
        refXhat = refextra.cursorVel;
        reflagXhat = k.getlagX(refXhat, refevent.blockStartStop, lag);
        refXhatnlag = [refXhat reflagXhat];
        personal_use_ref = refXhat(ind,:);
        header = 'personal use';

    case 'mixed tasks'
        % combine tasks
        % Center out
        [refData, ind] = k.GetSessionData(1);
        refXhat = Xhat(ind,:);
        reflagXhat = lagXhat(ind,:);
        % fitts
        path = rand_target_datapath;
        [refData1, ~, refevent, ~, refextra] = ConcatSavedSessionsData(path, p);
        ind = [ind' length(ind) + (1:length(refData1))];
        refXhat = [refXhat; refextra.cursorVel];
        reflagXhat = [reflagXhat; k.getlagX(refextra.cursorVel, refevent.blockStartStop, lag)];
        % personal use
        path = personal_use_datapath;
        [refData2, ~, refevent, ~, refextra] = ConcatSavedSessionsData(path, p);
        ind = [ind length(ind) + RowColon(refevent.trialStartStop(1:end,:))];
        refXhat = [refXhat ; refextra.cursorVel];
        reflagXhat = [reflagXhat; k.getlagX(refextra.cursorVel, refevent.blockStartStop, lag)];
        % combine all
        refXhatnlag = [refXhat reflagXhat];
        refData = [refData;refData1;refData2];
        combo_ref = refXhat(ind,:);
        header = 'mixed tasks';
end
disp('reference has been set.')
%
% calculate PCA given refrence data
k.CalculatePCA(refData(ind,:), ind, header);

% set reference as refrence data
k.SetReference(refData(ind,:), ind, header);

% calculate MINDFUL with NF + Xhat + Xhat_lag
k.CalcDistanceFromData('NFXhatlag', Xhatnlag(k.params.indSelected,:), refXhatnlag(ind,:) );

% % calculate MINDFUL given NF only as inputs
% k.CalcDistanceFromData('NF');
% 
% % calculate MINDFUL with NF + Xhat
% k.CalcDistanceFromData('NFXhat', Xhat(k.params.indSelected,:), refXhat(ind,:) );

% KLD_temp.pcaInfo = k.pcaInfo;
% 
% % no ND just kinematics data Xhat
% k.pcaInfo = [];
% k.data = Xhat(k.params.indSelected,:);
% k.SetReference(refXhat(ind,:), ind, header);
% 
% % calculate MINDFUL with Xhat only
% k.CalcDistanceFromData('Xhat')
% 
% % calculate MINDFUL with Xhat + Xhat_lag
% k.CalcDistanceFromData('Xhatlag', lagXhat(k.params.indSelected,:), reflagXhat(ind,:) )

% specify plot only KLD
% k.PlotCompFnvsAE(1,'KLdiv', {'NF','Xhat','Xhatlag','NFXhatlag'}, ...
%           {'NF','$$\hat{X}$$','$$\hat{X} + \hat{X}_{lag}$$', ...
%           ['NF', '$$ +\hat{X} + \hat{X}_{lag}$$']});
% text(-1.77,2.65, 'a', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
% text(-0.22,2.65, 'b', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
% text(-1.77,1.21, 'c', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
% text(-0.22,1.21, 'd', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')

KLD_temp.dist = k.dist;

switch ref_block
    case 'fitts'
        KLD_fitts = KLD_temp;
    case 'personal use'
        KLD_PU = KLD_temp;
    case 'center out'
        KLD_co = KLD_temp;
    case 'mixed tasks'
        KLD_combo = KLD_temp;
end

end
%% plot to compare the NFXhatlag in all cases of various reference
k.dist = [];
k.dist.fitts = KLD_fitts.dist.NFXhatlag;
k.dist.pu = KLD_PU.dist.NFXhatlag;
k.dist.co = KLD_co.dist.NFXhatlag;
k.dist.combo = KLD_combo.dist.NFXhatlag;
k.PlotCompFnvsAE(1,'KLdiv',{'co', 'fitts', 'pu', 'combo'}, ...
    {'Center-Out-And-Back (day 0)', 'Random Target (day 7)', ...
    'Personal Use (day 0)', 'Mixed Tasks'})
set(gcf,'Position',[0 0 .65 .65],'Units','normalized')
text(-1.77,2.65, 'a', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-0.22,2.65, 'b', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-1.77,1.21, 'c', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-0.22,1.21, 'd', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
%% plot Xhat polar histogram per segment

p.edges = linspace(-pi,pi,18*4+1);
% p.startStop = event.trialStartStop;
p.rlim = [0 .03]; % probability range

figure('Units','normalized','Position',[0 0 .65 .4])
axesHandle = subplot(1,4,1);
polaraxes('Units',axesHandle.Units,'Position',axesHandle.Position);delete(axesHandle);
discretizePolar(co_ref, 1, p);
title({'Center out task (day 0)', sprintf('(n = %i)',length(co_ref))},'FontSize',16)

axesHandle = subplot(1,4,3);
polaraxes('Units',axesHandle.Units,'Position',axesHandle.Position);delete(axesHandle);
discretizePolar(personal_use_ref, 1, p);
title({'Personal use (day 0)', sprintf('(n = %i)',length(personal_use_ref))},'FontSize',16)

axesHandle = subplot(1,4,2);
polaraxes('Units',axesHandle.Units,'Position',axesHandle.Position);delete(axesHandle);
discretizePolar(fitt_ref, 1, p);
title({'Fitts (day 7)', sprintf('(n = %i)',length(fitt_ref))},'FontSize',16)

axesHandle = subplot(1,4,4);
polaraxes('Units',axesHandle.Units,'Position',axesHandle.Position);delete(axesHandle);
discretizePolar(combo_ref, 1, p);
title({'Mixed tasks (day 0,7)', sprintf('(n = %i)',length(combo_ref))},'FontSize',16)

sgtitle('Decoded kinematics distributions')

