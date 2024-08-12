function [data, extra] = RemoveOutliers( data, opt )
% [data, commonFilled, featFilled] = RemoveOutliers( data, opt )
%
% Removes outliers in a data [nTime x nFeat] matrix using filloutliers.
%
%
% We remove two types of outliers:
%   1) (common) Outliers that are common across all features
%   2) (feat) Outliers per feature
%
%
%
% Inputs
%   data [nTime x nFeat]
%       matrix to be scanned for outliers
%
%   opt (see dfOpt for default options)
%       Any optional parameters for how to compute and threshold outliers
%       for 'common' and 'feat' outlier removal.
%
% Outputs
%   data
%       input data matrix but with outliers replaced as indicated in the
%       opt structure.
%
%   featFilled [nTime x nFeat]
%       A logical matrix of what points in data were considered outliers
%
%   commonFilled [nTime x 1]
%       A logical vector of what time points were considered outliers
%
%
%
% If using matlab < 2017b, can use custom RemoveNoiseOutliers.m
%
%
% History:
%   2019.09.22   Created
%
% Copyright Tommy Hosman, Brown University All Rights Reserved.
%--------------------------------------------------------------------------


% common outliers config
dfOpt.common.enable       = true;
dfOpt.common.interpMethod = 'linear'; % 'center' | 'clip' | 'previous' | 'next' | 'nearest' | 'linear' | 'spline' | 'pchip' | 'makima'
dfOpt.common.threshMethod = 'quartiles'; % Median has slightly tighter bounds
dfOpt.common.threshFactor = 4;


% per feature outlier config
dfOpt.feat.enable         = true;
dfOpt.feat.interpMethod   = 'linear';
dfOpt.feat.threshMethod   = 'customPercentile'; % Use custom percentile function below
% For customPercentile
dfOpt.feat.percentileStride = 2; % How many steps to skip to help speed up percentile calc
dfOpt.feat.percentiles      = [99.5]; % Used if threshMethod is customPercentile. If scalar, then lower threshold is 100-percentiles.
% dfOpt.feat.wirelessOutliers_prct      = []; % If not empty, estimates .percentiles bounds by excluding noise steps > excludeWirelessOutlierThreshold percent
% dfOpt.feat.excludeWirelessOutlierThreshold = 3; % Exclude steps with more than excludeWirelessOutlierThreshold percent noise when estimating percentiles.
% dfOpt.feat.ignoreStepsAroundOutliers = 25; % If threshMethod is customPercentile and this is > 0, we exclude ignoreStepsAroundOutliers steps around exclusion epochs when estimating percentiles. 
dfOpt.feat.excludeCommonOutlierSteps = true; % If threshMethod is customPercentile, then we exclude steps that were common outliers
dfOpt.feat.hardOutlierBounds = []; % 1x2 vector [lowerBound upperBound] to specify min and max feature boundaries. If not empty, will check against.

% threshFactor
% How much to multiply the threshold by. If len 2, first value is for lower
% bound, second is upper. If length 1, lower and upper bound are the same.
dfOpt.feat.threshFactor = 2;

dfOpt.ignoreSteps = logical([]); % Size of data. Specify time steps to ignore.

dfOpt.plot.enable = true;
dfOpt.plot.ch      = [];
dfOpt.plot.saveLoc = [];
dfOpt.plot.plotStr = ''; % Can pass in for plot label, e.g. the feature name
dfOpt.plot.figNum  = 1001;

if nargin < 2, opt = []; end
opt = MergeBintoA(dfOpt, opt);

%--------------------------------------------------------------------------
% First remove common (across array) outliers
%--------------------------------------------------------------------------
if opt.common.enable
    %%
    commonData = ComputeCommonData(data,opt);
    
    commonDataRaw = nanmedian(data,2);
    [~, commonLb, commonUb] = isoutlier(commonData, opt.common.threshMethod, 'thresholdfactor', opt.common.threshFactor);
    commonFilled = commonDataRaw < commonLb | commonDataRaw > commonUb; 
    %%

%%
    if any(commonFilled)
        % Found some outliers
        commonFilled = repmat(commonFilled, 1, size(data,2));
        data = filloutliers(data,opt.common.interpMethod, 'outlierLocation', commonFilled);
    end
else
    commonFilled = false;
    commonLb = nan;
    commonUb = nan;
end


%--------------------------------------------------------------------------
% Now remove outliers per channel
%--------------------------------------------------------------------------
orgData = data;

if opt.feat.enable
    %% Calc per-feature outliers
    opt.feat.commonFilled = commonFilled(:,1);
    
    [data,featFilled, lthresh, uthresh] = FillFeatureOutliers(orgData, opt);
        
    %%
%      opt.plot.ch = 47;
    if ~isempty(opt.plot.ch)
        PlotPlayPlace(orgData, opt, opt.plot.ch)
    end
else
    featFilled = false;
    lthresh = nan;
    uthresh = nan;
end


% Plot
if opt.plot.enable
    PlotOutlierDetection( orgData, data, opt, commonFilled, featFilled, lthresh, uthresh );
    drawnow
end

% Can plot individual features with
% PlotPlayPlace(orgData, opt, 105)

%% Extra output
extra.common.filled     = commonFilled;
extra.common.lowerBound = commonLb;
extra.common.upperBound = commonUb;

extra.feature.filled     = featFilled;
extra.feature.lowerBound = lthresh;
extra.feature.upperBound = uthresh;

end

function [data,outliers, lthresh, uthresh] = FillFeatureOutliers(data, opt)

    if strcmp(opt.feat.threshMethod, 'customPercentile')
        [lthresh, uthresh, featAvg, featPrctl] = RemovePercentileFactor( data, opt );
    else
        [~, lthresh, uthresh, featAvg] = isoutlier(data, opt.feat.threshMethod, 'thresholdfactor', opt.feat.threshFactor);
    end
    
    
    % Handle hard boundaries
    [lthresh,uthresh] = ThresholdTheThresholds(lthresh,uthresh,opt.feat.hardOutlierBounds);
    % Detect outliers
    outliers = data < lthresh | data > uthresh; 
    % Fill outliers
    data = FillOutliersHandleEdgeCase(data, outliers, featAvg, opt.feat.interpMethod);
    
end

function [lthresh,uthresh] = ThresholdTheThresholds(lthresh,uthresh,hardOutlierBounds)
if ~isempty(hardOutlierBounds)
    lthresh(lthresh<hardOutlierBounds(1)) = hardOutlierBounds(1);
    uthresh(uthresh>hardOutlierBounds(2)) = hardOutlierBounds(2);
end
end

function [lthresh, uthresh, featAvg, featPrctl] = RemovePercentileFactor( data, opt )
% Thresholds by a percentile * a factor
%
% todo: allow for a lower percentile to be applied
%       this would be lower percentile - a factor




featPrctl   = EstFeatThreshold( data, opt )'; % Note the transpose
featAvg     = featPrctl(2,:); % Median
featBounds  = diff(featPrctl);

nBounds = size(featBounds,1);
if isscalar(opt.feat.threshFactor)
    threshFactor = repmat(opt.feat.threshFactor, nBounds,1);
end
if any(threshFactor < 0)
    disp(threshFactor)
    error('We assume that the upper threshFactor is greater than or equal to 0, but we found %.f', upperThreshFactor)
end


%% Set thresholds

lthresh  = featAvg - featBounds(1,:).*threshFactor(1);
uthresh  = featAvg + featBounds(2,:).*threshFactor(2);


end

function featPrctl = EstFeatThreshold(data, opt)
    [nPoints, nFeat] = size(data);
    percentilePerFeat = GetPercentilePerFeature(opt.feat.percentiles,nFeat);
    nPercentiles = size(percentilePerFeat,2);
    featPrctl = zeros(nFeat,nPercentiles);
    
    % Select steps for estimating percentile
    if ~isempty(opt.ignoreSteps)
        selectSteps = ~opt.ignoreSteps;
    else
        selectSteps = true(size(data));
    end
    
    % Ignore common outliers
    if opt.feat.excludeCommonOutlierSteps && ~isscalar(opt.feat.commonFilled) 
        selectSteps(opt.feat.commonFilled,:) = false;
    end
    
    % Stride percentile estimate?
    if opt.feat.percentileStride > 1
        strideInds = 1:opt.feat.percentileStride:nPoints;
        data = data(strideInds,:);
        selectSteps = selectSteps(strideInds,:);
    end
    
    
    % Calc percentile per feature
    for featI = 1:nFeat
        selI = selectSteps(:,featI);
        % Estimate percentile
        featPrctl(featI,:) = prctile( data(selI,featI), percentilePerFeat(featI,:) );
    end
end

function data = FillOutliersHandleEdgeCase(data, outliers, centerVal, interpMethod, fillArgs)
% data = FillOutliersHandleEdgeCase(data, outliers, centerVal, interpMethod, fillArgs)
% Handles edge case when an outlier is the first or last index in data
% 
% data: [nStep x nDim]
% outliers: logical [nStep x nDim]
% centerVal: [1 x nDim]
% interpMethod: see variable fill in filloutliers.m
% fillArgs: optional additional arguments to go into filloutliers.m
if nargin < 5
    fillArgs = {};
end

% Edge case!
% Given interpMethods like linear, pchip, spline and outliers on the first
% or last point of data, append the center to allow interpolation to be
% bounded by the median.
% For example:
%     data = [ones(1,5)*1000 400:-200:0];
%     outliers = false(size(data));
%     outliers(1:5) = true;
%     filledData = filloutliers(data,'linear', 'outlierLocation', outliers);

isOutlierOnEdge = any(outliers(1,:)) || any(outliers(end,:));
isSuseptableInterpMethod = ismember(interpMethod, {'linear', 'pchip','spline', 'previous', 'next'});
appendEdges = isSuseptableInterpMethod && isOutlierOnEdge;

if appendEdges
    data = cat(1, centerVal, data, centerVal);
    notOutlier = false(size(centerVal));
    outliers = cat(1, notOutlier, outliers, notOutlier);
end


%% Fill outliers
data = filloutliers(data,interpMethod,'outlierLocation', outliers,fillArgs{:});


%% Edge case clean up
if appendEdges
    data([1 end],:) = [];
end

end

function percentilePerFeat = GetPercentilePerFeature(percentilePerFeat, nFeat)

    if isscalar(percentilePerFeat)
        percentilePerFeat = repmat(percentilePerFeat,nFeat,1);
    end

    % Append lower bound and median
    medianPercentile = repmat(50,nFeat,1);
    if isvector(percentilePerFeat)
        lowerBound = 100-percentilePerFeat;
        percentilePerFeat = cat(2, lowerBound, medianPercentile, percentilePerFeat);
    end

end

function commonData = ComputeCommonData(data,opt)
    nFeat = size(data,2);
    
    if ~isempty(opt.ignoreSteps) && isequal(size(data), size(opt.ignoreSteps))
        data(opt.ignoreSteps) = nan;
    end
       
    commonData = nanmedian(data,2);
end


function PlotOutlierDetection( orgData, data, opt, commonFilled,  featFilled, lthresh, uthresh )
nPlt = 1;


    [nTime, nFeat] = size(data);
        nFilledPerFeat = sum(featFilled);
        [sv, si] = sort( nFilledPerFeat, 'descend' );

        uCh = si(sv(1:nPlt)>0);
%         % Append an average feature to the end
%         uF = mean(data);
%         [~, suf] = sort(uF);
%         avgInd = suf(round(length(suf)/2));
%         uCh(end+1) = avgInd;
%         uCh = unique(uCh, 'stable');



        h = [];
        % Get the figure, but do not necessarily bring into focus
        fh = findobj('Number', opt.plot.figNum);
        if length(fh) == 1
            set(0, 'CurrentFigure', fh)
        else
            fh = figure(opt.plot.figNum);
        end
        set(gcf, 'Position',  [100, 0, 1000, 1000])
        fh.Name = 'Outlier summary';
        if ~isempty(opt.plot.plotStr)
            fh.Name = sprintf('Outlier summary for %s',opt.plot.plotStr);
        end
        clf
        
        
        nStableSp = 3;
        nSpRow = nStableSp + nPlt;
        
        
        % Plot number of outliers per feature
        subplot(nSpRow,1,1)
        barData = cat(1, repmat(sum(commonFilled(:,1)),1,nFeat), sum(featFilled))';
        bh = bar(barData, 'stacked');
        legend(bh, 'Common outliers', 'Per-feature outliers')
        
        xlabel('Feature number')
        ylabel('counts')
        titleStr = 'Number of outlier steps removed';
        if ~isempty(opt.plot.plotStr)
            loc = [0.005 0.7 0.15 0.3];
            ha = annotation('textbox', loc, 'string', opt.plot.plotStr, ...
                'edgecolor', 'none', 'Interpreter', 'none', 'FontSize', 16 );
        end
        title(titleStr)
        
        subplot(nSpRow,1,2);
        plot(uthresh)
        plot(lthresh)
        legend('Upper threshold', 'Lower threshold')
        
        xlabel('Feature number')
        ylabel('Threshold value')
        title('Outlier Thresholds')
        
        
        % Plot percent of features with outliers
        ax_nt(1) = subplot(nSpRow,1,3);
        plotT = (1:nTime).*0.02;
        percentFeatFilled = mean(featFilled,2).*100;
        percentCommonFilled = mean(commonFilled,2).*100;
        
        
        plot(plotT,percentFeatFilled, 'DisplayName', 'per-feature detected outliers')
        plot(plotT,percentCommonFilled, 'DisplayName', 'common detected outliers')
%         if ~isempty(opt.feat.wirelessOutliers_prct)
%             percentNs5Outliers = mean( opt.feat.wirelessOutliers_prct, 2 );
%             plot(plotT, percentNs5Outliers, 'DisplayName', 'NS5 outliers')
%         end
        legend()
        xlabel('Time (sec)')
        ylabel('%')
        title('Percent features with outliers')

        
        
        % Plot Features and their detected outliers
        if ~isempty(uCh)
            for uc = 1:length( uCh )

                ax_nt(uc+1) = subplot(nSpRow,1,uc+nStableSp);
                f = uCh(uc);

                [h, legendStr] = PlotOutlierFeature(f, orgData, data, featFilled, lthresh, uthresh );
                % If we decide we want to plot multiple features again...
    %             titleStr = sprintf('Outliers for feature %d (%d worst feat, %d outliers)', f, uc, nFilledPerFeat(f));

                titleStr = sprintf('Feature with most outliers\nFeature %d (%d outliers)', f, nFilledPerFeat(f));
                title(titleStr);
                if uc ~= length( uCh )
                    xlabel('')
                    ylabel('')
                end
            end
        
            warning off
            legend(h, legendStr{:} )
            warning on
        end
        % Last plotted feature is the average
%         titleStr = sprintf('Feature with average number of outliers (%d outliers)\nOutliers for feature %d', nFilledPerFeat(f), f);
%         title(titleStr);


        linkaxes(ax_nt, 'x')
        axis(ax_nt, 'tight')
        
        if ~isempty(opt.plot.saveLoc)
            drawnow
            print(opt.plot.saveLoc, '-dpng', '-r200')
        end
end

function [h, legendStr] = PlotOutlierFeature(f, orgData, data, featFilled, lthresh, uthresh )

[nTime, nFeat] = size(data);

    samplePeriod = 0.02;
    t = (1:nTime).*samplePeriod;

    if isvector(lthresh)
        lt = [lthresh(f) lthresh(f)];
        ut = [uthresh(f) uthresh(f)];
        threshT = [1 nTime];
    else
        lt = lthresh(:,f);
        ut = uthresh(:,f);
        threshT = 1:length(uthresh);
    end

    hold on;
    h(1) = plot(t,  data(:, f) );
    h(2) = plot( t(threshT), lt, 'k--');
    plot( t(threshT), ut, 'k--')

    legendStr{1} = 'Updated data';
    legendStr{2} = 'Outlier threshold';
    outlierLoc = find(featFilled(:,f));
    if ~isempty(outlierLoc)
        h(3) = plot(t(outlierLoc), orgData(outlierLoc,f), 'rx');
        legendStr(3:4) = {'Detected outliers'};
    end
%     xlim([min(t),max(t)]);

    xlabel('Time (sec)')
    ylabel('Feature value')
    title(sprintf('Outliers for feature %d', f))


end


function PlotPlayPlace(orgData, opt, f)
if ~exist('f','var') || isempty(f)
    f = 112;
end
%%
%     f = mod(324-1,192)+1;
%     ch = f;
%     f = 112;
    ch = mod(f-1,192)+1;
%     opt.feat.threshFactor = 2; %3;

    [data,featFilled, lthresh, uthresh] = FillFeatureOutliers(orgData, opt);
%     [data,featFilled, lthresh, uthresh]  = filloutliers(data,opt.feat.interpMethod, 'movmean', [9e3 0], 'thresholdfactor', opt.feat.threshFactor);
    
%    SetFigureFontSize(20,555)
%     figure(555)
    clf
    sp(1) = subplot(211);
    t = (1:length(orgData))./50;
    hh = plot(t, orgData(:,f));


    [h, legendStr] = PlotOutlierFeature(f, orgData, data, featFilled, lthresh, uthresh );

    legend(h, legendStr{:});
    
    sp(2) = subplot(212);
    plot(t, opt.feat.wirelessOutliers_prct(:,ch), 'DisplayName', 'Percent wireless outliers')
    overThresh = opt.feat.wirelessOutliers_prct(:,ch);
    overThresh(opt.feat.wirelessOutliers_prct(:,ch)<opt.feat.excludeWirelessOutlierThreshold) = nan;
    plot(t, overThresh, '-x', 'DisplayName', 'Over wireless outlier threshold')
    xlabel('Time (sec)')
    ylabel('Percent wirelses outliers')
    legend();
    linkaxes(sp,'x')

end

