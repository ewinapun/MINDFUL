%% params sweep
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
lagXhat = k.getlagX(Xhat, k.event.blockStartStop, lag);  

[refData, ind] = k.GetSessionData(init_day);
header = ['session ',num2str(init_day)];

% comment out if not selecting time steps
subsampleAE_max = 4;
ind = ind(k.extra.angleError(ind) < subsampleAE_max);
refData = k.data(ind,:);
header = [header,' AE < ', num2str(subsampleAE_max)];
 
k.CalculatePCA(refData, ind, header);
k.SetReference(refData, ind, header);
disp(k)
%% sweep window length
displabels=[];
k.dist = [];
wlist = [1:5 10 15 20 30 45 60 90 120 180 240]; % window length in second  
startInd = k.params.startStop(:,1);
for wl = wlist
    disp(wl)
    k.params.startStop(:,2) = startInd + wl/0.02 - 1;
    k.params.startStop(k.params.startStop > length(k.data)) = length(k.data);
    fieldname = ['winlen_', num2str(wl),'sec'];
    k.CalcDistanceFromData(fieldname,[Xhat(k.params.indSelected,:) lagXhat(k.params.indSelected,:)])
end

k.GetMedianAEperSeg()
fn = fieldnames(k.dist);
sweep_winlen = k.PlotCompDistvsAE(1,'KLdiv',fn);
save(['sweep_winlen_',info.participant],'-struct','sweep_winlen')
%% plot correlations for both participants
try
    t11 = load('sweep_winlen_T11');
t5 = load('sweep_winlen_T5');
figure('units','inches','Position',[10 1 3.5 3],'defaultAxesfontsize',10.5, ...
    'DefaultLineLineWidth',2)

% h(1) = plot(wlist,t11.Spearman,'-k','DisplayName','Spearman');
% h(2) = plot(wlist,t11.Pearson,'--k','DisplayName','Pearson');
% legend(h);legend('boxoff')

plot(wlist,t11.Spearman,'Marker','x','Color',c(1,:),'MarkerSize',10)
plot(wlist, t5.Spearman,'Marker','.','Color',c(2,:),'MarkerSize',20)
% plot(wlist,t11.Pearson,'--','Marker','x','Color',c(1,:),'MarkerSize',10);
% plot(wlist, t5.Pearson,'--','Marker','.','Color',c(2,:),'MarkerSize',20);

text(200,0.9,'T11','FontSize',16,'Color',c(1,:))
text(210,0.7,'T5','FontSize',16,'Color',c(2,:))

xticks([0 30 60 120 180 240])
grid on
xlabel('window length (seconds)')
ylabel('Spearman correlation of KLD to AE')
axis([-5 240 -0.3 1])
if saveGenFigure;savepdf(gcf,'sweep_WL');end
end
%% AE quantiles
k = MINDFUL(NDzc, event, info, params, extra);

% get a lagged version of Xhat
lag = 1;
Xhat = extra.cursorVel;
lagXhat = k.getlagX(Xhat, k.event.blockStartStop, lag);  

% set reference segment
init_day = 1;
if ismember(info.participant,'T5')
    init_day = [1,2];
end
csSegPerBlk = [0;cumsum(k.params.SegPerBlk)];
loc = find(ismember(k.event.sessionNumberPerBlock,init_day));
segInds = csSegPerBlk(loc(1))+1:csSegPerBlk(loc(end)+1) ...
            - round(k.params.winlen/k.params.updateHz) + 1;
[refData, ind] = k.GetSegmentData(segInds);
header = ['session ',num2str(init_day)];
k.CalculatePCA(refData, ind);
k.SetReference(refData, ind);   
k.CalcDistanceFromData('AE0180',[Xhat(k.params.indSelected,:) lagXhat(k.params.indSelected,:)])
displabels{1} = ['0-180', char(176)];
%
lenInd = length(ind);
ind2 = [];
AEqtile = [0, 4, 8, 13, 20, 30, 50, 100, 180];
for i = 1:length(AEqtile)-1
    ind2{i} = find(k.extra.angleError(ind) >= AEqtile(i) & k.extra.angleError(ind) <= AEqtile(i+1));
    lenInd(i+1) = length(ind2{i});
    displabels{i+1} = [num2str(AEqtile(i)),'-',num2str(AEqtile(i+1)), char(176)];
    k.CalculatePCA(refData(ind2{i},:), ind2{i})
    k.SetReference(refData(ind2{i},:), ind2{i});
    k.CalcDistanceFromData(['AE',num2str(AEqtile(i)),num2str(AEqtile(i+1))], ...
    [Xhat(k.params.indSelected,:) lagXhat(k.params.indSelected,:)]);
end
k.GetMedianAEperSeg()
disp(lenInd)
%% plot correlation for a participant
fn = fieldnames(k.dist);
sweepAE = k.PlotCompDistvsAE(1,'KLdiv',fn,displabels);
figure('Position',[100 100 1200 650],'defaultAxesfontsize',20)
x = categorical(displabels);
x = reordercats(x,displabels);
h = bar(x,[sweepAE.Pearson sweepAE.Spearman],1);
xtickangle(0)
ylim([0 1])
xline(1.5,'--','Color',[.5 .5 .5],'LineWidth',3)
legend(h,'Pearson','Spearman')
legend('boxoff')
title(info.participant)
for i = 1:length(x)
    text(i-0.3, -0.08, ['n = ',num2str(lenInd(i))])
end
xlabel('AE quantiles')
ylabel('Correlation of KLD to AE')
grid on
save(['sweepAE_',info.participant],'-struct','sweepAE')
%% plot Spearman correlation for both participants
t11 = load('sweepAE_T11');
t5 = load('sweepAE_T5');
figure('units','inches','Position',[1 1 10 7],'defaultAxesfontsize',20)

x = categorical(displabels);
x = reordercats(x,displabels);
h = bar(x,[t11.Spearman t5.Spearman],1);
xtickangle(45)
ylim([0.3 1])
xline(1.5,'--','Color',[.5 .5 .5],'LineWidth',3)
legend(h,'T11','T5')
legend('boxoff')
xlabel('AE quantiles')
ylabel('Spearman correlation of KLD to AE')
grid on
savepdf(gcf,'sweepAE_T11_T5(Spearman)')

% plot Pearson correlation for both participants
figure('units','inches','Position',[1 1 10 7],'defaultAxesfontsize',20)
x = categorical(displabels);
x = reordercats(x,displabels);
h = bar(x,[t11.Pearson t5.Pearson],1);
xtickangle(45)
ylim([0.4 1])
xline(1.5,'--','Color',[.5 .5 .5],'LineWidth',3)
legend(h,'T11','T5')
legend('boxoff')
xlabel('AE quantiles')
ylabel('Pearson correlation of KLD to AE')
grid on
savepdf(gcf,'sweepAE_T11_T5(Pearson)')