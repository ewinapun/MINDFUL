%% figure 4: dpca plot

tstart=8;tend=57;
day = 1;
feats = NDzc(RowColon(event.sessionStartStop(day,:)),:);
feats = GaussSmooth2(feats, 50);
clrs = [    
    0.6118    0.1059    0.2706
    0.8196    0.2392    0.3137
    0.9294    0.3725    0.2667
    0.9686    0.6000    0.3412
    0.9922    0.8118    0.4745
    0.3373    0.7098    0.6627
    0.1725    0.4745    0.6667
    0.3725    0.3137    0.6275
];
tnum = event.trialsPerSession(day);
csTPS = [0; cumsum(event.trialsPerSession)];
trialInds = csTPS(day)+1:csTPS(day+1);
trialInds = trialInds(~event.excludeTrials(trialInds));
[catMoveDir, uCatDirs] = ChangeDirsToCategorical(event.movedirVect, 0);
goalLabels = double(catMoveDir);
% gcolors = clrs(goalLabels+1,:);

trialStart = event.trialStartStop(trialInds,1) - event.sessionStartStop(day,1) + 1;
window = (tstart:tend);
margNames = 'Directions';
dpcaParams = [];
dpcaOut = GenericDPCA(feats, goalLabels(trialInds), trialStart, window, margNames, dpcaParams);
dirDepPC = find(dpcaOut.whichMarg==1);
dirInDepPC = find(dpcaOut.whichMarg==2);

%%
pltwin = (tstart:tend);
subplot(2,3,5)
components = [dirDepPC(1),dirDepPC(2)];
for dir = 1:8
    trials = dpcaOut.trialMap(dir,:);
    trials = trials(~isnan(trials));
    trials = trials + csTPS(day);
    inds = event.trialStartStop(trials,1) - event.sessionStartStop(day,1) + 1 + pltwin;
    pltData = feats(inds,:)*dpcaOut.W(:,components);
    pltData = reshape(pltData,[size(inds),size(pltData,2)]);
    meandir = squeeze(mean(pltData,1,'omitnan'));
    plot(meandir(:,1),meandir(:,2),'Color',(clrs(dir,:)),'linewidth',3)
    for t = 1:length(trials)
        plot(pltData(t,:,1),pltData(t,:,2),'Color',[(clrs(dir,:)) 0.15])
    end
end
xlabel(['dPC',num2str(components(1))])
ylabel(['dPC',num2str(components(2))])

subplot(2,3,6)
components = [dirDepPC(1),dirInDepPC(1)];
for dir = 1:8
    trials = dpcaOut.trialMap(dir,:);
    trials = trials(~isnan(trials));
    trials = trials + csTPS(day);
    inds = event.trialStartStop(trials,1) - event.sessionStartStop(day,1) + 1 + pltwin;
    pltData = feats(inds,:)*dpcaOut.W(:,components);
    pltData = reshape(pltData,[size(inds),size(pltData,2)]);
    meandir = squeeze(mean(pltData,1,'omitnan'));
    plot(meandir(:,1),meandir(:,2),'Color',(clrs(dir,:)),'linewidth',3)
    for t = 1:length(trials)
        plot(pltData(t,:,1),pltData(t,:,2),'Color',[(clrs(dir,:)) 0.15])
    end
end
xlabel(['dPC',num2str(components(1))])
ylabel(['dPC',num2str(components(2))])

%% project assessment blocks data onto first session dpca space
varExp = zeros(length(event.pointsPerSession),20); % top 20 PCs
varExpProj = zeros(length(event.pointsPerSession),20);
p.plotFun = [];
%variance explain
for day = 1:length(event.pointsPerSession)
    feats = NDzc(RowColon(event.sessionStartStop(day,:)),:);
    feats = GaussSmooth2(feats, 50);
    trialInds = csTPS(day)+1:csTPS(day+1);
    trialInds = trialInds(~event.excludeTrials(trialInds));
    inds = event.trialStartStop(trialInds,1)+pltwin;
    inds = inds - event.sessionStartStop(day,1) + 1;
    events = event.trialStartStop(trialInds,1)- event.sessionStartStop(day,1) + 1;
    [dpcaFeatures, dpcaLabels] = GetEpochOfData(feats, goalLabels(trialInds), events, window, 'dpca');
    [dpcaOut2, firingRates, firingRatesAverage] = SetupAndComputedPCA(dpcaFeatures, dpcaLabels, p);
    varExp(day,:) = dpcaOut2.explVar.componentVar;
    explVar = dpca_explainedVariance(firingRatesAverage, dpcaOut.W, dpcaOut.V, ...
                        'combinedParams', dpcaOut.combinedParams, ...
                        'X_trial', firingRates, ...
                        'numOfTrials', dpcaOut.numTrials, ...
                        'Cnoise', dpcaOut.Cnoise, ...
                        'centerData', dpcaOut.params.centerData);
    varExpProj(day,:) = explVar.componentVar;
    
%     dpcaOutall{day}.W = dpcaOut2.W;
%     dpcaOutall{day}.V = dpcaOut2.V;
%     [COEFF, ~, ~, ~, EXPLAINED] = pca(feats(inds,:));
end
%% projection plot with VAF 
components = [dirDepPC(1),dirDepPC(2)];
figure('Units','inches','Position',[0 0 14 3])

for day = 1:nday
    subplot(ceil(nday/5),5,day)
    meandir = zeros(length(pltwin),2,8);
    disp(['Projecting day ', num2str(day)])
    trialInds = csTPS(day)+1:csTPS(day+1);
    trialInds = trialInds(~event.excludeTrials(trialInds));
    dirLabelday = goalLabels(trialInds);
    feats = NDzc(RowColon(event.sessionStartStop(day,:)),:);
    feats = GaussSmooth2(feats, 50);
for dir = 1:8
    trials = trialInds(dirLabelday==dir);
    trials = trials(~isnan(trials));
    inds = event.trialStartStop(trials,1)+pltwin;
    inds = inds - event.sessionStartStop(day,1) + 1;
    pltData = feats(inds,:)*dpcaOut.W(:,components);
    pltData = reshape(pltData,[size(inds),size(pltData,2)]);
    meandir(:,:,dir) = squeeze(mean(pltData,1,'omitnan'));
    for t = 1:length(trials)
        plot(pltData(t,:,1),pltData(t,:,2),'Color',[(clrs(dir,:)) 0.1])
    end
for dir = 1:8
    plot(meandir(:,1,dir),meandir(:,2,dir),'Color',(clrs(dir,:)),'linewidth',3)
end
end
axis([-5 5 -2.5 2.5])
axis off
title([{['decoder day ',num2str(info.trialDay(day)-info.trialDay(1))]}...
     ,{['(',num2str(sum(abs(varExpProj(day,1:2))),2),'% VAF)']}],'FontSize',14)
end
% savepdf(gcf,'dpca_VAF')

%% Pearson correlation of session mean KL to VAF (compare to first session)
% pwmeanKL obtained from pairwise_mean_KLD.m
if ~exist('pwmeanKL','var')
    run pairwise_mean_KLD.m
end
[corrout_dPCA,PVAL_dPCA] = corr(pwmeanKL(1,:)',VAF);
fprintf('Pearson correlation of mean KL to VAF = %.3f, p = %.1e.\n', corrout_dPCA, PVAL_dPCA)
figure(11)
for day = 1:nday
    scatter(pwmeanKL(1,day),VAF(day),30,k_cmp(day,:),'filled')
    text(pwmeanKL(1,day),VAF(day)-2, ...
        num2str(info.trialDay(day)-info.trialDay(1)),'FontSize',16)
end
title(['Pearson correlation = ',num2str(corrout_dPCA,2)])
ylabel('VAF (%)');xlabel('mean KLD')
axis([0 3 0 70])

%% Pearson correlation of Angle Error to VAF
mAE = zeros(nday,1);
nPC = 2;
VAF = sum(abs(varExpProj(:,1:nPC)),2);
for day = 1:nday
    % trial indices for each day
    trialInds = csTPS(day)+1:csTPS(day+1);
    trialInds = trialInds(~event.excludeTrials(trialInds));
    mAE(day) = mean(extra.angleErrorPerTrial(trialInds));
end
[corr_AE, PVAL_AE] = corr(mAE,sum(abs(varExpProj(:,1:2)),2));
fprintf('Pearson correlation of mean angle error to VAF = %.3f, p = %.1e.\n', corr_AE, PVAL_AE)
figure(10)

if strcmp(info.participant,'T11');rmv = 20;else;rmv = 5;end
k_cmp = brewermap(xticksday(end)+1+2*rmv,'RdYlBu');
mid = ceil(length(k_cmp)/2);
k_cmp(mid-rmv:mid+rmv,:)=[];
k_cmp = vertcat(k_cmp(1,:)*0.5,k_cmp);

for day = 1:nday
    scatter(mAE(day),VAF(day),30,k_cmp(day,:),'filled')
    text(mAE(day)-2,VAF(day)-2, ...
        num2str(info.trialDay(day)-info.trialDay(1)),'FontSize',16)
end
title(['Pearson correlation = ',num2str(corr_AE,2)])
ylabel('VAF (%)');xlabel('mean angle error(\circ)')
axis([0 180 0 70])
