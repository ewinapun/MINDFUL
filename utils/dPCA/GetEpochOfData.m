function [epochFeatures, epochLabels, epoch] = GetEpochOfData(features, labels, eventInds, window, outputType, epoch)
% [epochFeatures, epochLabels, epoch] = GetEpochOfData(features, labels, eventInds, window, outputType)
% 
% Given event inds and a window, indexes into features and labels and outputs
% in standard shapes.
% 
% NOTE! This function will exclude trials whose event + window goes outside
%       the features/labels (e.g. negative inds, or > size(features,1).
% 
% Inputs:
% 
%  features [nPoints x nDim_feat]
%    Features to be indexed into.
% 
%  labels [nPoints x nDim_lbl]
%    Label values for the respective features.
% 
%  eventInds
%    Vector of events to align the window vector, window. E.g. trial start
%    inds.
% 
%  window
%    Vector of offset inds from eventInds used to index into features and
%    labels. E.g. 0:49
% 
%  outputType (default false)
%    How the output is reshaped (F: features, L: labels):
%     false (0) F: [nWindowTrial x nFeat],      L: [nWindowTrial x nDim_lbl]
%     true  (1) F: [nWindow x nEvents x nFeat], L: [nWindow x nEvents x nDim_lbl]
%     'dPCA'    F: [nEvents x nWindow x nFeat], L: [nEvents x nDim_lbl]
% 
% 
% Output:
%   epochFeatures
%       Indexed and reshaped features. Shape dependent on outputType.
%   epochLabels
%       Indexed and reshaped labels. Shape dependent on outputType.
%   epoch
%       Inds used to index into features, labels. Shape [nWindow x nEvents]
% 
% History:
% Copyright Tommy Hosman, All Rights Reserved
%   2019   

if nargin < 4
    window = 0:49;
end
if nargin < 5
    outputType = 0;
end
if nargin < 6
    epoch = [];
end
nPoints = size(features,1);
lblSize = size(labels); % Is updated below after possible label shape adjustments
nEvents = length(eventInds);

if isvector(window)
    nWindow = length(window);
elseif ismatrix(window)
    nWindow = size(window,1);
else
    error('Unexpected window shape: %s', num2str(size(window)));
end
    

% Transpose the labels (if labels are nDim x nPoints)
if length(lblSize) == 2 && lblSize(2) == nPoints
    labels = labels';
end
isemptyLabel = all(lblSize==0);


% Get the requested epoch from eventInds and window
if isempty(epoch)
    epoch = GetWindowAroundEvents(eventInds, window);
end
excludeTrls = (epoch(end,:) > nPoints);
epoch(:,excludeTrls) = [];


nWindowTrial = numel(epoch);
lblSize      = size(labels);

% Is the label events x nDim (as opposed to standard nPoints x nDim)
if lblSize(1) == nEvents
    labelIsPerEvent = true;
    labels(excludeTrls,:) = [];
else
    labelIsPerEvent = false;
end


% Reshape features to be nWindow x nEvents x nFeat
epochFeatures = reshape(features(epoch,:), [size(epoch) size(features,2)]);

%% Output shape
% dPCA  F: [nEvents x nWindow x nFeat], L: [nEvents x nDim_lbl]
% true  F: [nWindowTrial x nFeat],      L: [nWindowTrial x nDim_lbl]
% false F: [nWindow x nEvents x nFeat], L: [nWindow x nEvents x nDim_lbl]
if ischar(outputType) && strcmpi(outputType, 'dpca')
    % outputType dPCA  F: [nEvents x nWindow x nFeat], L: [nEvents x nDim_lbl]
    
    % One label per trial. 
    %
    % Try to find the label at the event index (e.g. when windInds==0).
    % Otherwise take the first label value.
    labelInd = find(window==0);
    if isempty(labelInd)
        labelInd = 1;
    end
    
    % nEvents x nWindow x nFeat
    epochFeatures = permute(epochFeatures,[2,1,3]);
    % nEvents x nDim_lbl
    if ~isemptyLabel
        if labelIsPerEvent
            epochLabels = labels;
        else
            epochLabels = reshape(labels(epoch(labelInd,:),:), [nEvents, size(labels,2)]);
        end
    end
    
    
elseif ~outputType %  && ndims(epochFeatures) >= 3
    % outputType false  F: [nWindowTrial x nFeat], L: [nWindowTrial x nDim_lbl]    
    %
    %
    % From: nWindow x nEvents x nFeat -> (nWindow,nEvents) x nFeat
    epochFeatures = reshape(epochFeatures, nWindowTrial, size(epochFeatures,3));    
    
    if ~isemptyLabel
        if labelIsPerEvent
            % Repeat the event label nWindow (nWindow x nEvents x nDim_lbl)
            epochLabels = repmat( reshape(labels, [1 size(labels)]), nWindow, 1, 1 );
            % Reshape to be nWindowTrial x nDim_lbl
            epochLabels = reshape(epochLabels, nWindowTrial,lblSize(2:end));
            
        else
            epochLabels = labels(epoch(:),:);
        end
    end


else
    % outputType true F: [nWindow x nEvents x nFeat], L: [nWindow x nEvents x nDim_lbl]
    % nWindow x nEvents x nFeat
    if ~isemptyLabel
        if labelIsPerEvent
            % Repeat the event label nWindow (nWindow x nEvents x nDim_lbl)
            epochLabels = repmat( reshape(labels, [1 size(labels)]), nWindow, 1, 1 );
        else
            epochLabels = reshape(labels(epoch,:), [size(epoch), size(labels,2)]);
        end
    end
end

if isemptyLabel
    epochLabels = [];
end


end