function plotPerformanceSummaryPerTrial(perf, info, events, params)
% plot performance metrics trial by trial
% 
% Input
%     perf
%         performance metrics
%         .angleErrorPerTrial
%         .orthChanges
%         .pathEfficiency
%         .timeToTarget
%         .trialSuccess
%     info & events
%         can be obtained by
%         [cat_ND, cat_labels, events, info, cat_extra] = ConcatSavedSessionsData(path, p)
%     params 
%         additional parameters
%     
% Output
%   none
%
% Copyright Tsam Kiu Pun, 2023. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

defaultParams.singleSession = 1;
defaultParams.trialDay = [];
defaultParams.xlabelFromDayZero = 1;
defaultParams.skipExcludeTrials = 0;
defaultParams.avgBlocksSuccessRate = 1;
% Optional inputs
if nargin < 4
    params = [];
end
params = MergeBintoA(defaultParams, params);

if isfield(info,'trialDay')
    params.trialDay = info.trialDay;
end
nday = length(info.trialDay);

if nday > 1
    params.singleSession = 0;
end
c = [0    0.4470    0.7410]; % default blue color
if params.avgBlocksSuccessRate
    avgCorrectPerSession = zeros(nday,1);
    csTrial = 1+[0;cumsum(events.trialsPerSession)];
    for day = 1:nday
        tInd = csTrial(day):csTrial(day+1)-1;
        tInd(events.excludeTrials(tInd)) = [];
        avgCorrectPerSession(day) = mean(perf.trialSuccess(tInd))*100;
    end
    perf.percentCorrect = avgCorrectPerSession;
else
    csTrial = 1+[0;cumsum(events.trialsPerBlock)];
    for b = 1:length(events.trialsPerBlock)
        tInd = csTrial(b):csTrial(b+1)-1;
        tInd(events.excludeTrials(tInd)) = [];
        perf.percentCorrect(b) = mean(perf.trialSuccess(tInd))*100;
    end
end

if params.skipExcludeTrials && any(events.excludeTrials)
    perf.angleErrorPerTrial(events.excludeTrials) = [];
    perf.orthChanges(events.excludeTrials) = [];
    perf.pathEfficiency(events.excludeTrials) = [];
    perf.timeToTarget(events.excludeTrials) = [];
    perf.trialSuccess(events.excludeTrials) = [];
end

[xlabelstr, ticklabels, tickInds] = GetPlotLabels(info, events, params);

%%
figure
set(gcf, 'Position',  [100, 100, 1000, 800])
set(0,'defaultAxesFontSize',20)

% trial Success
subplot(5,1,1)
bar(perf.percentCorrect)
if length(perf.percentCorrect) < 20
for i = 1:length(perf.percentCorrect)
    text(i-0.4,perf.percentCorrect(i) - 10,...
        [num2str(perf.percentCorrect(i),3)],...
        'Color','w','FontSize',15)
end
end
if ~params.avgBlocksSuccessRate
    xticks(find(diff([0;events.sessionNumberPerBlock])~=0))
else
    xticks(1:length(perf.percentCorrect))
end
axis tight
box off
title('Success Rate')

xticklabels(ticklabels)
failTrialNum = find(perf.trialSuccess==0);
ylabel('%')

if params.singleSession
    xticks(1:length(ticklabels))
    xticklabels(ticklabels)
end

x = [0;cumsum(perf.timeToTarget(1:end-1))];

% AE
subplot(5,1,2)
if ~params.skipExcludeTrials && any(events.excludeTrials)
    plot(x(events.excludeTrials), perf.angleErrorPerTrial(events.excludeTrials), ...
        'o','Color',[.5 .5 .5],'MarkerFaceColor',[.7 .7 .7])
end
plot(x, perf.angleErrorPerTrial,'.','Color',c)
hold on
plot(x(failTrialNum), perf.angleErrorPerTrial(failTrialNum),'.r')
title('Angle Error Per Trial')
ylabel('Deg')
axis tight

% time to target
subplot(5,1,3)
if ~params.skipExcludeTrials && any(events.excludeTrials)
    plot(x(events.excludeTrials), perf.timeToTarget(events.excludeTrials), ...
        'o','Color',[.5 .5 .5],'MarkerFaceColor',[.7 .7 .7])
end
plot(x, perf.timeToTarget,'.','Color',c)
hold on
plot(x(failTrialNum), perf.timeToTarget(failTrialNum),'.r')
title('Time to target Per Trial')
ylabel('Seconds')

% path efficiency
subplot(5,1,4)
% fill nan values, path efficiency seems wrong
perf.pathEfficiency(isnan(perf.pathEfficiency)) = 0;
if ~params.skipExcludeTrials && any(events.excludeTrials)
    plot(x(events.excludeTrials), perf.pathEfficiency(events.excludeTrials), ...
        'o','Color',[.5 .5 .5],'MarkerFaceColor',[.7 .7 .7])
end
plot(x, perf.pathEfficiency,'.','Color',c)
hold on
plot(x(failTrialNum), perf.pathEfficiency(failTrialNum),'.r')
title('Path Efficiency Per Trial')
ylabel('Best = 1')

% orthogonral directional changes
subplot(5,1,5)
h(1) = plot(x, perf.orthChanges,'.','Color',c, ...
    'DisplayName','successful trials');
hold on
h(2) = plot(x(failTrialNum), perf.orthChanges(failTrialNum),'.r', ...
        'DisplayName','unsuccessful trials');
if ~params.skipExcludeTrials && any(events.excludeTrials)
    h(3) = plot(x(events.excludeTrials), perf.orthChanges(events.excludeTrials), ...
        'o','Color',[.5 .5 .5],'MarkerFaceColor',[.7 .7 .7], ...
        'DisplayName','outlier trials');
end
plot(x, perf.orthChanges,'.','Color',c);
plot(x(failTrialNum), perf.orthChanges(failTrialNum),'.r');
title('Orth. Direction Change Per Trial')
ylabel('Counts')

for splot = 2:5
    subplot(5,1,splot)
    xticks(x(tickInds))
    xticklabels(ticklabels)
    xlim([0, x(end)])
    AddDayTransitionLines(x(tickInds))
end
xlabel(xlabelstr)
leg = legend(h);
leg.Box = 'off';
leg.Position(1) = 0.75;
leg.Position(2) = -0.0056;
leg.FontSize = 14;
end

function [xlabelstr,ticklabels,tickInds] = GetPlotLabels(info, events, p)
    if p.singleSession
    % label multiple blocks instead of days
        xlabelstr = 'Block';
        xtickstr = info.usedBlocks{:};
        ticklabels = xtickstr;
        tickInds = [0;cumsum(events.trialsPerBlock(1:end-1))]+1;
    else
    % label multiple sessions days
        if p.xlabelFromDayZero
            xlabelstr = 'Days since first session';
            xtickstr =  info.trialDay(events.sessionNumberPerBlock) - info.trialDay(1);
            ticklabels = xtickstr(diff([-1 xtickstr])~=0);
        else
            xlabelstr = 'Trial day';
            xtickstr = info.trialDay(events.sessionNumberPerBlock);
            ticklabels = info.trialDay;
        end
        tickInds = [0;cumsum(events.trialsPerSession(1:end-1))]+1;
    end
    if p.skipExcludeTrials && any(events.excludeTrials)
        for i = 2:length(tickInds)
            tickInds(i) = tickInds(i) - sum(events.excludeTrials(1:tickInds(i)-1));
        end
    end
end

function AddDayTransitionLines(loc)
    for ii = 1:length(loc)
        xline(loc(ii),'Color',[.5 .5 .5]);
    end
end