function [xi, means, ci] = getAvgBinFeatures(ND, labels, options)
    arguments
        ND double
        labels double
        options.binwidth = 10 % bin width in degree
        options.bw = 22.5 % bandwidth for regression in degree
        options.kernelSmooth = true
        options.alpha = 0.05; % confidence interval
    end
    nfeats = size(ND,2);
    angles = atan2d(labels(:,2), labels(:,1));
    xi = (-180:options.binwidth:180)'; % range of the bins
    nbins = length(xi);
    means = zeros(nbins, nfeats); % mean
    ci = zeros(nbins, nfeats); % confidence intervals
        
    if options.kernelSmooth
    % using Nadaraya-Watson kernel regression
        % repeat in range [x - 2pi, x, x + 2pi for accounting edge effect
        disp('running Nadaraya-Watson kernel regression')
        x = vertcat(angles-360, angles, angles+360);

        for f = 1:nfeats
            disp(['feature #', num2str(f)])
            y = repmat(ND(:, f),[3,1]);
            [r, ci_diff, xi, bw] = KernelRegression(x, y, xi, options.bw, options.alpha);
            means(:,f) = r;
            ci(:,f) = ci_diff;
            % bootstrap to get confidence intervals
%             r_all = zeros(100, nbins);
%             for b = 1:100
%                 perm = randi(length(x),1,length(x));
%                 r_all(b,:) = KernelRegression(x(perm),y(perm),xi,'bandwidth',options.bw);
%             end
%             means(:,f) = mean(r_all,1);
%             ci(:,f) = 1.96 * std(r_all);
        end
    else
        % Binned neural data
        [cnt, binInds, ~] = PolarHist(angles, nbins, 1);
        err_t = zeros(nbins, 2);
        % Estimate mean
        for f = 1:nfeats
            for t = 1:nbins
                [means(t,f), ~, err_t(t,:)] = normfit(ND(binInds==t, f));
            end
            ci(:,f) = means(:,f) - err_t(:,1); % symmetric error bar
        end
        % shift display order from (0, 360) to (-180, +180)
        means = means([5:8,1:5],:);
        ci = ci([5:8,1:5],:);
    end
end