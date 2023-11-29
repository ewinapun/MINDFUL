%% figure 3: plot tunings

set(0,'defaultAxesfontsize',20)

% get tunings per session
run fit_tunings.m

% get pairwise mean KLD between sessions
run pairwise_mean_KLD.m

% initialize object using a class for PlotTunings
pt = PlotTunings(tunings);

% get MD and PD difference from the reference session day
% mask is 0 if tuning p-value is < 0.05, NaN otherwise
[deltaMD, deltaPD, refDay, tune_feats, mask] = pt.getDeltaTuning('refFirstTune', true);

nanmask = @(x) x+mask; 
nanzscore = @(x) (x-nanmean(x,'all')) ./ nanstd(x,1,'all');
%% plot tunings for all features
sort_data = [nanzscore(nanmask(pt.tunings.md)) nanzscore(nanmask(pt.tunings.pd))];
[~, order] = pt.GetSortedOrder(sort_data, 'distMetric', 'euclidean','linkage', 'ward');
pt.selected_feats = order;

pt.PlotAllTuning('xtickstr', xticksday)
savepdf(gcf,'pdmd_all_feats')
pt.PlotAllTuning('xtickstr', xticksday, 'gray_out_non_sig', true) % gray out the non-significant features
savepdf(gcf,'pdmd_all_feats_sig')

%% now only plot days with more than half days are significantly tuned
sign_tune_feats = find(sum(~isnan(mask),2) > floor(nday/2) & tune_feats);

% hierarchical clustering to order among selected features
sort_data = [nanzscore(deltaMD) nanzscore(deltaPD)];
[~, sig_order] = pt.GetSortedOrder(sort_data(sign_tune_feats,:) , ...
                'distMetric', 'euclidean', 'linkage', 'ward');
sig_order = sign_tune_feats(sig_order);

% possibly additional ordering
% [~,lowerest_md_ind] = min(nansum(pt.tunings.md(sig_order,:).^2,2));
% sig_order = [sig_order(lowerest_md_ind+1:end);sig_order(1:lowerest_md_ind)];

pt.selected_feats = sig_order;

%% gray out the non-significant features
pt.PlotAllTuning('xtickstr', xticksday, 'gray_out_non_sig', true)
savepdf(gcf,'pdmd_per_feats_sig')

hFig = figure(100);
if strcmp(info.participant,'T11')
set(hFig,'Units','inches','Position',[1 2 11.5 6]);
pt.PlotDeltaTuning(deltaMD, deltaPD, 'fig', hFig, 'xticks', (1:2:nday),...
    'xtickstr', xticksday(1:2:end), skip_diff_col=false, gray_out_non_sig=true)
elseif strcmp(info.participant, 'T5')
set(hFig,'Units','inches','Position',[1 2 11 6]);
pt.PlotDeltaTuning(deltaMD, deltaPD, 'fig', hFig, ...
    'xtickstr', xticksday, skip_diff_col=false, gray_out_non_sig=true)
end

savepdf(gcf,'change_pdmd_per_feats_sig')
%% plot tuning curves across sessions
if strcmp(info.participant,'T11')
    channel = [191, 112];
elseif strcmp(info.participant, 'T5')
    channel = [128, 20];
end

feat = find(sig_order==channel(1)); 
[~, cmp] = pt.PlotExampleTuning(sig_order(feat),'fig',figure(101), ...
    'style','cosine', ...
    'ylim', 1.2743, ...
    'title',['Feat ',num2str(feat)]);

feat = find(sig_order==channel(2)); 
pt.PlotExampleTuning(sig_order(feat),'fig',figure(102), ...
    'style','cosine', ...
    'ylim', 1.2743, ...
    'title',['Feat ',num2str(feat)], ...
    'colorbar',xticksday);

[corrR, pVal] = pt.GetPairwiseCorrelation('fig',figure(103), ...
    'xtickstr',xticksday, ...
    'interp',xticksday(end)+1, ...
    'cbXTick',[0,1]);

%% plot between sessions with color days apart
if strcmp(info.participant,'T11');rmv = 20;else;rmv = 5;end
k_cmp = brewermap(xticksday(end)+1+2*rmv,'RdYlBu');
mid = ceil(length(k_cmp)/2);
k_cmp(mid-rmv:mid+rmv,:)=[];
k_cmp = vertcat(k_cmp(1,:)*0.5,k_cmp);

days_apart = abs(repmat(xticksday,nday,1)-xticksday');
sess_apart = abs(repmat((1:nday),nday,1)-(1:nday)');
ind = find(triu(repmat((1:nday),nday,1))~=0); % only upper triangle
% color by days apart
corrFig = figure(104);
% set(corrFig,'Units','inches','Position',[1 1 2.8 2.3])
scatter(pwmeanKL(ind),corrR(ind),[],k_cmp(days_apart(ind)+1,:),'filled')
[corrout,corr_pval] = corr(pwmeanKL(ind),corrR(ind),'type','Pearson','rows','complete');
xmax = max(pwmeanKL(ind));
axis([0 xmax 0 1.02])
axis square
xticks([0 1 2 3])
axh = gca();
text(xmax*0.6, .98,['${\it} r$ = ',num2str(corrout,3)],'Interpreter','Latex')
text(xmax*0.6,.88,['${\it} p$ = ',num2str(corr_pval,1)],'Interpreter','Latex')
ylabel('R'); xlabel('KL divergence')
colormap(axh,k_cmp); 
cbh = colorbar();
ylabel(cbh,'Days apart')
clim = cbh.Limits;
y = linspace(clim(1),clim(2),xticksday(end)+1);
cbh.Ticks = y(1:rmv:end)+ (y(2)-y(1))/2;
cbh.TickLabels = (0:rmv:xticksday(end));
% savepdf(gcf,'KL_vs_corrR_days_apart_triu')

newfig=figure(105);
set(newfig,'Units','inches','Position',[1 1 7.5 6])
tcl = tiledlayout(newfig,'flow','Padding','tight','TileSpacing','tight');
for i = 1:4
    figure(100+i);f = gcf();figure(f)
    ax=gca;ax.Parent=tcl;ax.Layout.Tile=i;
    close(f)
end
figure(newfig)
savepdf(gcf,'tuning_plotOut')

%% tuning examples with kernel regression
figure(106)
set(gcf,'Units','inches','Position',[1 1 3 6])

feat = find(sig_order==channel(1)); 
pt.PlotExampleTuning(sig_order(feat),'fig',subplot(2,1,1), ...
    'style','kernreg', ...
    'ylim', 1.2743, ...
    'title',['Feat ',num2str(feat)], ...
    'colorbar',xticksday);
feat = find(sig_order==channel(2)); 
pt.PlotExampleTuning(sig_order(feat),'fig',subplot(2,1,2), ...
    'style','kernreg', ...
    'ylim', 1.2743, ...
    'title',['Feat ',num2str(feat)], ...
    'colorbar',xticksday);
savepdf(gcf,'tuning_examples_kernel_reg')

%% calculate tuning significance statistics
[sigdeltaPD, sigdeltaMD] = bootstrapDifference(tunings, refDay);

sigdeltaMD = sigdeltaMD(sig_order,:);
sigdeltaPD = sigdeltaPD(sig_order,:);
tuneddMD = deltaMD(sig_order,:);
tuneddPD = deltaPD(sig_order,:);
% pt.PlotSignificantHeatmap(sigdeltaMD, sigdeltaPD)

%  tuned features exhibited significant change in both MD and PD in at least one session
sigdeltaPDMD = floor((sigdeltaPD+sigdeltaMD)/2);
sig_change_feats = find(sum(sigdeltaPDMD,2, 'omitnan') >= 1)';
fprintf(['%i out of %i tuned features that exhibited significant changes' ...
    ' in both MD and PD in at least one session\n'], ...
    length(sig_change_feats), length(sig_order));

f = (1:length(sig_order));
sigTuneddMD = tuneddMD(f,2:end).*sigdeltaMD(f,2:end);
sigTuneddPD = tuneddPD(f,2:end).*sigdeltaPD(f,2:end);

rmv = 1; k_cmp = brewermap(nday + rmv * 2,'RdYlBu');
mid = ceil(length(k_cmp)/2); k_cmp(mid-rmv:mid+rmv,:)=[]; k_cmp = flipud(k_cmp);
k_cmp = vertcat(k_cmp(1,:)*0.5,k_cmp);
figure;
for i = f;scatter(sigTuneddMD(i,:),sigTuneddPD(i,:),20,k_cmp(1:end-1,:),'filled');end
sigTuneddMD((sigTuneddMD)==0) = nan;
sigTuneddPD((sigTuneddPD)==0) = nan;

%% get statistics
avgdMD = nanmean(sigTuneddMD,1)';
avgdPD = nanmean(sigTuneddPD,1)';
[r,r_pval] = corr(avgdMD,avgdPD,'type','Pearson','rows','complete');
% figure;scatter(avgdMD,avgdPD,100,'.');xlabel('average \DeltaMD');ylabel('average \DeltaPD')

% declaring sessions with poor performance
if strcmp(info.participant,'T11');sbad = (12:15)-1;else; sbad = (4:5)-1;end
sgood = setxor((1:size(sigTuneddPD,2)),sbad); % good performance

% PD
v1 = sigTuneddPD(:,sbad); v1 = v1(~isnan(v1));
v2 = sigTuneddPD(:,sgood);v2 = v2(~isnan(v2));
fprintf('sessions with poor performance are %s\n',num2str(sbad+1))
fprintf('avg dPD in poor sessions = %.3g %c %.3g\n', mean(v1),char(177),std(v1))
fprintf('avg dPD in good sessions = %.3g %c %.3g\n:', mean(v2),char(177),std(v2))
fprintf('significant difference: p = %.3g\n\n',ranksum(v1, v2))

% MD
v1 = abs(sigTuneddMD(:,sbad)); v1 = v1(~isnan(v1));
v2 = abs(sigTuneddMD(:,sgood));v2 = v2(~isnan(v2));
fprintf('avg dMD in poor sessions = %.3f %c %.3f\n', mean(v1),char(177),std(v1))
fprintf('avg dMD in good sessions = %.3f %c %.3f\n:', mean(v2),char(177),std(v2))
fprintf('significant difference: p = %.3g\n',ranksum(v1, v2))

%% extra
% plot p values changes of distribution from day 0
% pt.PlotBasicImagesc(pvalues_sigma(order,:), title='p-value (\Sigma)', clim=[0 0.05], xtickstr=xticksday)
% savepdf(gcf,'Bonferroni_pvalue(sigma)')
% pt.PlotBasicImagesc((pvalues_beta(sig_order,1:end)), title='p-value (\beta)', clim=[0 0.05], xtickstr=xticksday)
% savepdf(gcf,'Bonferroni_pvalue(beta)')
% save('tunings_val','offset','t','tunings','xticksday', ...
%     'deltaMD','deltaPD','sig_order','order', ...
%     'pvalues_sigma','pvalues_beta');

% correlations
[corrR, pVal] = pt.GetPairwiseCorrelation('xtickstr',xticksday,'interp',xticksday(end)+1);
% savepdf(gcf,'corr_tuneMap_pw_interp_fine')
% [corrR, pVal] = pt.GetPairwiseCorrelation('xtickstr',xticksday,'interp',round((xticksday(end)+1)/2));
% savepdf(gcf,'corr_tuneMap_pw_interp_coarse')
[corrR, pVal] = pt.GetPairwiseCorrelation('xtickstr',xticksday);
% savepdf(gcf,'corr_tuneMap_pw')

% plot mean KL from day 0 against tuning correlation
figure('Units','inches','Position',[1 1 3 3])
xKL = pwmeanKL(1,:);
scatter(xKL,corrR(1,:),[],k_cmp,'filled')
[corrout,corr_pval] = corr(pwmeanKL(1,:)',corrR(1,:)','type','Pearson','rows','complete');
sprintf('Pearson r = %.3f, p = %.3e', corrout, corr_pval)
title(['${\it} r$ = ',num2str(corrout,3), 'p =',num2str(corr_pval)],'Interpreter','Latex')
ylabel('tuning similarity')
xlabel('mean KL-div')
axis([ min(xKL)*0.8 max(xKL) 0.4 1.1])
% savepdf(gcf,'KL_vs_corrR_day0')
