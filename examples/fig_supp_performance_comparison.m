%% supplemental: performance plot
% Code released with manuscript: Pun et al., "Measuring instability in 
% multi-day human intracortical neural recordings towards stable, 
% long-term brain-computer interfaces".
%
% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

%% plot performance summary
plotPerformanceSummaryPerTrial(extra, info, event, params)

%% script to compare performance between early and later sessions
if ismember(info.participant,'T11');n = 11;else;n = 3;end
fieldname = 'trialSuccess';
pe = extra.(fieldname);
tcutoff = sum(event.trialsPerSession(1:n));
peearly = pe(1:tcutoff);peearly = peearly(~isnan(peearly));
pelater = pe(tcutoff+1:end);pelater = pelater(~isnan(pelater));

fprintf('%s: first %i sessions %.2f %c %.2f\n',fieldname,n, mean(peearly)*100,char(177),std(peearly)*100)
fprintf('%s: later sessions %.2f %c %.2f\n',fieldname,mean(pelater)*100,char(177),std(pelater)*100)
fprintf('significant difference: p = %.5f\n',ranksum(peearly,pelater))
sus = extra.angleErrorPerTrial(extra.trialSuccess);
unsus = extra.angleErrorPerTrial(~extra.trialSuccess);
fprintf('\n%s - early: %.2f%c%.2f; later: %.2f%c%.2f; p<0.01 \n',...
        fieldname, mean(sus),char(177),std(sus),...
        mean(unsus),char(177),std(unsus))
fieldnames = {'angleErrorPerTrial','timeToTarget','pathEfficiency','orthChanges'};
for fn = 1:4
    pe = extra.(fieldnames{fn});
    % pe = extra.(fieldname)(~event.excludeTrials);
    tcutoff = sum(event.trialsPerSession(1:n));
    peearly = pe(1:tcutoff);peearly = peearly(~isnan(peearly));
    pelater = pe(tcutoff+1:end);pelater = pelater(~isnan(pelater));
    fprintf('\n%s - early: %.2f%c%.2f; later: %.2f%c%.2f; p<0.01 \n',...
        fieldnames{fn}, mean(peearly),char(177),std(peearly),...
        mean(pelater),char(177),std(pelater))
    fprintf('significant difference: p = %.5f\n',ranksum(peearly,pelater))
end