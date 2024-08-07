%% code for figures in the MINDFUL paper
% Copyright Ewina Pun, All Rights Reserved

% add MINDFUL repo to path
filePath  = mfilename('fullpath');
addpath(genpath(fileparts(fileparts(filePath))))

% set up data path were the data is
dirpath = 'Y:\TransferSpace\MINDFUL(2024)';

% specify which participant's data to load
participant = 'T11'; % T11 or T5

% concat path
path = [dirpath, filesep, participant];

% set to true if want to save output figures
saveGenFigure = 1;

% figure settings
set(0,'defaultAxesFontSize',20)
c = get(groot,'DefaultAxesColorOrder');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

% concatenate parsed data
clear p
p.zscoreFeatures = true; % set false for raw extract features
p.trailingMoving = true; % set true for causal rolling z-score using movmean; set false for block batch z-score
[NDzc, labels, event, info, extra] = ConcatSavedSessionsData(path, p);
info.participant = participant;

% set reference day
if strcmp(info.participant,'T5')
    init_day = [1, 2];
elseif strcmp(info.participant,'T11')
    init_day = 1;
end

nfeats = size(NDzc,2);
nday = size(event.sessionStartStop,1);
xticksday = info.trialDay(1:end) - info.trialDay(1);

%% fig 1b. Group ND by angle error
run fig1_grouped_by_AE

%% fig 1c, 2a-c. MINDFUL against AE
run fig2_plot_MINDFUL

%% figure 3: plot tunings
run fig3_plot_tunings.m

%% figure 4: dpca plot
run fig4_dpca_projection.m

%% figure 5: parameter sweep
run fig5_param_sweep.m

%% supplemental: performance plot
run fig_supp_performance_comparison.m

