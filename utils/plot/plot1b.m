function plot1b(k, sub_ind,ref_ind,titlestr,info)

if strcmp(info.participant,'T11')
    color = [0.3725    0.3137    0.6275];
elseif strcmp(info.participant,'T5')
    color = [0.3373    0.7098    0.6627];
end
refData = k.data(ref_ind,:); 
header = '';
k.CalculatePCA(refData, ref_ind, header);

sTo = [];
AErange = (0:4:180-4);AErange = [AErange;AErange+4]';
data_after_pca = k.PCAtransform(k.data);
ref_pd = ProbDistibutionEst(data_after_pca(ref_ind,:));
statusprint = fprintf('Running 0 of %d', length(AErange));
data_after_pca = data_after_pca(sub_ind,:);
ae = k.extra.angleError(sub_ind);
for j = 1:length(AErange)-1
% group by even AE interval
    fprintf(repmat('\b',1,statusprint))
    statusprint = fprintf('Running %d of %d\n', j, length(AErange)-1);
    ind = find(ae>=AErange(j,1) & ae<AErange(j,2));
    pd = ProbDistibutionEst(data_after_pca(ind,:));
    st = k.CalcDistance(ref_pd, pd);
    st.mAE = nanmean(ae(ind));
    sTo = ConcatStruct(sTo, st);
end

scatter(sTo.mAE, sTo.KLdiv, 20,'k','filled',...
        "MarkerEdgeColor",color, ...
        "MarkerFaceColor",color)
xlabel('Angle Error'); ylabel('KLD')
axis tight; axis square
ylim([0 max(round(max(sTo.KLdiv)),max(sTo.KLdiv))*1.2])
xlim([0 180])
xticks([0 90 180])
xticklabels({['0',char(176)], ['90',char(176)], ['180',char(176)]})
title(titlestr)

[r,r_pval] = corr(sTo.mAE,sTo.KLdiv,'type','Pearson','rows','complete');
fprintf('Pearson correlation of mean KL to mAE = %.3f, p = %.1e.\n', r, r_pval)

limAxis = double(axis);
text(10, limAxis(4)*0.9, ['r=',num2str(r,2)], 'FontSize',16,'Color',color)

% dlm = fitlm(sTo.mAE,sTo.KLdiv,'Intercept',false);
% x = linspace(limAxis(1),limAxis(2),20)';
% y = double(x*dlm.Coefficients.Estimate);
% plot(x, y,'k--','linewidth',2)
end