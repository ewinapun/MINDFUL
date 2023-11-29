function dpcaOut = GenericDPCA(features, labels, events, window, margNames, dpcaParams, warpAlign, maxNs5Outliers, maxNs5OutlierThreshold)
% Wrapper to call dPCA given events, winInds and other common parameters
% Note, for features and labels, the size does matter
% 
% Required
%   features [nSteps x nFeatDim]
%   labels   [nSteps x nLabelDim] Or [nEvents x nLabelDim]
%   events   [nEvents x 1]
%   window   [nWin x 1]
if nargin < 5
    margNames = 'Stimulus';
end
if nargin < 6
    dpcaParams = [];
end
if nargin < 7
    warpAlign = false; % Requires additional library
end
if nargin < 8
    maxNs5Outliers = [];
end
if nargin < 9
    maxNs5OutlierThreshold = 5;
end



%% Handle trial exclusion and features, labels extraction

% Labels should be [nStep x nDim], or [nTrials x nDim]
if isrow(labels)
    labels = labels(:);
end


% Exclude events
excludeEvents = false(size(events));
if ndims(features) == 2
    excludeOutOfBounds = (events+max(window)) > size(features,1) |  (events+min(window)) < 1;
end
excludeEvents(excludeOutOfBounds) = true;
if ~isempty(maxNs5Outliers)
    excludeOutlierEvents = ExcludeOutlierEvents(maxNs5Outliers, maxNs5OutlierThreshold, events, window);
    excludeEvents(excludeOutlierEvents) = true;
end





if ndims(features) == 3
    features = features(~excludeEvents,:,:);
end
if size(labels,1) == length(events)
    labels = labels(~excludeEvents,:);
end
excludeEventVals = events(excludeEvents);
events(excludeEvents) = [];

if ndims(features) == 3
    % nTrials x nTime x nFeat
    dpcaFeatures = features;
    dpcaLabels = labels;
else
    [dpcaFeatures, dpcaLabels] = GetEpochOfData(features, labels, events, window, 'dpca');
end




%% Remove trials based on average feature response

[~, rmvExtra] = RemoveOutliers(squeeze(mean(dpcaFeatures,2)));
isAvgSigOutlier = any( rmvExtra.feature.filled, 2 );


[~, Xpc] = pca(squeeze(mean(dpcaFeatures,2)));
[tf, lthresh, uthresh, center] = isoutlier(Xpc(:,1:5),'median', 'ThresholdFactor', 4.);
remvOutlierTrials = any(tf,2);


% 
% anyOut = find(remvOutlierTrials | isAvgSigOutlier);
% bp = Xpc(anyOut,1:3);
% figure(234234)
% clf
% bar(bp)
% px = find(ismember(anyOut,find(isAvgSigOutlier)));
% py = max(bp(px,:),[],2);
% plot(px,py, 'rx')
% 
% 
% px = find(ismember(anyOut,find(remvOutlierTrials)));
% py = max(bp(px,:),[],2);
% plot(px,py, 'ks')

isAvgSigOutlier = isAvgSigOutlier | remvOutlierTrials;
%%
% Do not exclude if more than 1% was excluded
nRmvRsp = sum(isAvgSigOutlier);
if nRmvRsp> 0 && nRmvRsp/length(isAvgSigOutlier) < 0.05
    fprintf('\nRemoving avg sig trials! (%d)\n', sum(isAvgSigOutlier))
    excludeEventVals = sort(cat(1, excludeEventVals, events(isAvgSigOutlier)));
    events(isAvgSigOutlier) = [];
    dpcaLabels(isAvgSigOutlier,:) = [];
    dpcaFeatures(isAvgSigOutlier,:,:) = [];
end

if ~isempty(excludeEventVals)
    nExcludeEvents = length(excludeEventVals);
    prctExcludeEvents = length(excludeEventVals) / length(excludeEvents)*100;
    fprintf('Excluding %d (%.1f%%) events\n', nExcludeEvents, prctExcludeEvents);
end

%% Setup and run dPCA
if ~iscategorical(dpcaLabels)
    dpcaLabels = categorical(dpcaLabels);
end
dpcaLabels = removecats(dpcaLabels);

dpcaP = [];
dpcaP.plotPCA = false; %true;
dpcaP.pltTime = window./50;
dpcaP.marginalizationNames = AddMarginalization(dpcaLabels, margNames, 'Condition-independent');
if length(dpcaP.marginalizationNames) == 3
    addInteraction = 1;
    for ii = 1:size(dpcaLabels,2)
        if length(unique(dpcaLabels(:,ii))) == 1
            addInteraction = 0;
        end
    end
    if addInteraction
        dpcaP.marginalizationNames = AddMarginalization(dpcaLabels, dpcaP.marginalizationNames, 'Interaction');
    end
end
% dpcaP.lambda = 0.0025;
dpcaP = MergeBintoA(dpcaP,dpcaParams);


%% run dPCA
dpcaOut = SetupAndComputedPCA(dpcaFeatures, dpcaLabels, dpcaP);



%% Run Warp alignment if requested
if ~isempty(warpAlign) && (isstruct( warpAlign ) || warpAlign)
    dpcaOut = WarpAndRunDPCA(warpAlign, dpcaOut, dpcaFeatures, dpcaLabels, dpcaP);
end



%% Save info about what epoch and label we used
taskInfo.labels = dpcaLabels;
taskInfo.events = events;
taskInfo.window = window;
taskInfo.excludeEvents = excludeEventVals;

dpcaOut.taskInfo = taskInfo;

end

function excludeEvents = ExcludeOutlierEvents(maxNs5Outliers, maxNs5OutlierThreshold, eventInds, winInds)
    outlierEpoch = GetEpochOfData(maxNs5Outliers(:), [], eventInds, winInds, 'dpca');
    prctOutliersPerEvent = mean(outlierEpoch,2);
    excludeEvents = prctOutliersPerEvent>maxNs5OutlierThreshold;
end

function marginalizationNames = AddMarginalization(dpcaLabels, marginalizationNames, addMarg)

if ~iscell(marginalizationNames)
    marginalizationNames = {marginalizationNames};
end
if ~iscell(addMarg)
    addMarg = {addMarg};
end

if length(marginalizationNames) < size(dpcaLabels,2)
    error('Found %d label columns but only %d marginalization names. Please add a marginalization name per label column.', ...
        size(dpcaLabels,2),length(marginalizationNames));
end

if ~any(ismember(marginalizationNames,addMarg))
    marginalizationNames = cat(2, marginalizationNames, addMarg);
end

end




function dpcaOut = WarpAndRunDPCA(warpAlign, dpcaOut, dpcaFeatures, dpcaLabels, dpcaP)

if isequal(warpAlign,2)
    sWarped = CallAffineWarp(); % Reload previous affine warp results
elseif isstruct(warpAlign)
    sWarped = CallAffineWarp(dpcaFeatures,warpAlign);
else
    sWarped = CallAffineWarp(dpcaFeatures);
end

if ~isempty(sWarped)
    fn = fieldnames(sWarped.warped_data);
    for ii = 1:length(fn)
        dpcaP.titleStr = sprintf('%s warped', fn{ii});
        dpcaOut.(['warped_' fn{ii}]) = SetupAndComputedPCA(sWarped.warped_data.(fn{ii}), dpcaLabels, dpcaP);
        drawnow
    end
end

end