%% KL Grouped ND by angle error conditioned to different direction as reference
% Code released with manuscript: Pun et al., "Measuring instability in 
% multi-day human intracortical neural recordings towards stable, 
% long-term brain-computer interfaces".
%
% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

clear k params ind ind_dir
params.excludeNonTrials = 1;
k = MINDFUL(NDzc, event, info, params, extra);
ae = k.extra.angleError;
X = labels(k.params.indSelected,:);
theta = atan2(X(:,2),X(:,1));

ind_up = find(theta>6*pi/16 & theta<10*pi/16);
ind_down = find(theta<-6*pi/16 & theta>-10*pi/16);
ind_left = find(abs(theta)>7*pi/8);
ind_right = find(abs(theta)<1*pi/8);

% plot KL vs AE
direction{1} = 'left';  ind_dir{1} = ind_left;
direction{2} = 'up';    ind_dir{2} = ind_up;
direction{3} = 'right'; ind_dir{3} = ind_right;
direction{4} = 'down';  ind_dir{4} = ind_down;
p.edges = linspace(-pi,pi,18*4+1);
if strcmp(info.participant,'T11')
p.rlim = [0 .35]; % probability range
elseif strcmp(info.participant,'T5')
p.rlim = [0 .201]; % probability range
end
set(0,'defaultAxesFontSize',20)
figure('Units','inches','Position',[5 0 12 12]);

for row_dir = 1:4
    opposite_dir = mod(row_dir+2,4); 
    if opposite_dir==0;opposite_dir = 4;end
    ind_row_dir = ind_dir{row_dir};
    ref_ind = ind_row_dir(ae(ind_row_dir) < 4); 

    axesHandle = subplot(4,4,1+4*(row_dir-1));
    polaraxes('Units',axesHandle.Units,'Position',axesHandle.Position);
    delete(axesHandle);
    discretizePolar(X(ref_ind,:), 1, p);
    thetaticks(0:22.5:360);
    set(gca,'fontsize',14)
    Thetalabel = get(gca,'ThetaTickLabel');Thetalabel(2:2:end) = {''};
    set(gca,'ThetaTickLabel',Thetalabel)
    h=gca;
    h.Position(1)=0.10;
    h.Position(2)=h.Position(2)-0.04;
    h.Position(3)=0.13;
    title(['Reference: ' direction{row_dir}],'FontSize',20)
    
    subplot(4,4,2+4*(row_dir-1));
    plot1b(k,(1:length(ae)),ref_ind,'on all',info)
    if row_dir ~= 4; xlabel('');end
    
    subplot(4,4,3+4*(row_dir-1));
    plot1b(k,ind_dir{row_dir},ref_ind,['on ',direction{row_dir}],info)
    if row_dir ~= 4; xlabel('');end
    
    subplot(4,4,4+4*(row_dir-1));
    plot1b(k,ind_dir{opposite_dir},ref_ind,['on ',direction{opposite_dir}],info)
    if row_dir ~= 4; xlabel('');end
end
sgtitle(info.participant, 'Fontsize', 32, 'FontWeight', 'bold')
if strcmp(info.participant,'T11')
text(-6.35,6.70, 'a', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-6.35,4.85, 'b', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-6.35,3.00, 'c', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-6.35,1.15, 'd', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
elseif strcmp(info.participant,'T5')
text(-6.35,6.70, 'e', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-6.35,4.85, 'f', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-6.35,3.00, 'g', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
text(-6.35,1.15, 'h', 'Fontsize', 32, 'FontWeight', 'bold','Units','normalized')
end
savepdf(gcf,sprintf('fig 1b diff ref dir(%s)',info.participant))
%%
params.updateHz = 500;
params.excludeNonTrials = 0;
params.xlabelFromDay0 = 1;

Xhat = extra.cursorVel;
lag = 1;
lagXhat = k.getlagX(Xhat, k.event.blockStartStop, lag);
Xhatnlag = [Xhat lagXhat];

theta = atan2(labels(:,2),labels(:,1));
ae = k.extra.angleError;
dirstr = {'All  ','Left ', 'Down ','Right', 'Up   '};

k = MINDFUL(Xhatnlag, event, info, params, extra);
[~, ind] = k.GetSessionData(init_day);
ind = ind(ae(ind) < 4);
refData = k.data(ind,:);
k.SetReference(refData, ind, header);
k.CalcDistanceFromData('all');

[~, ind] = k.GetSessionData(init_day);
dVel = discretize(theta(ind), linspace(-pi,pi,17)); dVel(dVel==16)=0;
% dirstr = {'All Dirs','Left','Down-Left', 'Down','Down-Right','Right','Up-Right', 'Up','Up-Left'};
header = 'center out day 0; difference directions';
for i=1:4
    ind = find(dVel==4*(i-1) | dVel==4*i-3);
    ind_dir = ind(ae(ind) < 4);zz ind 
    figure;discretizePolar(labels(ind,:), 1, p);
    refData_dir = k.data(ind_dir,:);
    % set reference as refrence data
    k.SetReference(refData_dir, ind, header);
    k.CalcDistanceFromData(['d',num2str(i)]);
end
%
corrs = k.PlotCompDistvsAE(3,'KLdiv',[],dirstr);
set(gcf,'Position',[0 0 .3 .35],'Units','normalized')
title('$$ \hat{X} + \hat{X}_{lag}$$','fontsize',gca().FontSize*1.5,'interpreter','latex')
yyaxis right;ylim([0 1.1])
ylabel('Normalized KLD')
savepdf(gcf,sprintf('fig 2c diff ref dir Xhat(%s)',info.participant))
%%
params.excludeNonTrials = 0;
ae = k.extra.angleError;
k = MINDFUL(NDzc, event, info, params, extra);
[~, ind] = k.GetSessionData(init_day);
ind = ind(ae(ind) < 4);
refData = k.data(ind,:);
k.CalculatePCA(refData, ind, header);
k.SetReference(refData, ind, header);
k.CalcDistanceFromData('all',Xhatnlag(k.params.indSelected,:));
%
[~, ind] = k.GetSessionData(init_day);
dVel = discretize(theta(ind), linspace(-pi,pi,17)); dVel(dVel==16)=0;
% dirstr = {'All Dirs','Left','Down-Left', 'Down','Down-Right','Right','Up-Right', 'Up','Up-Left'};
header = 'center out day 0; difference directions';
for i=1:4
%     ind_dir = find(dVel==2*(i-1) | dVel==2*i-1);
    ind_dir = find(dVel==4*(i-1) | dVel==4*i-3);
    ind_dir = ind_dir(ae(ind_dir) < 4);zz ind 
%     figure;discretizePolar(Xhat(ind_dir,:), 1, p);
    refData_dir = k.data(ind_dir,:);
    k.CalculatePCA(refData_dir, ind, header);
    k.SetReference(refData_dir, ind_dir, header);
    k.CalcDistanceFromData(['d',num2str(i)],Xhatnlag(k.params.indSelected,:));
end
%
corrs = k.PlotCompDistvsAE(3,'KLdiv',[],dirstr);
set(gcf,'Position',[0 0 .3 .35],'Units','normalized')
title('NF + $$ \hat{X} + \hat{X}_{lag}$$','fontsize',gca().FontSize*1.5,'interpreter','latex')
yyaxis right;ylim([0 1.1])
ylabel('Normalized KLD')
savepdf(gcf,sprintf('fig 2c diff ref dir(%s)',info.participant))
