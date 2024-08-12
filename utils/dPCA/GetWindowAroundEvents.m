function outputInds = GetWindowAroundEvents(events, window, removeLessThanEqZero)
% This function takes in a vector of events (e.g. [10 5030]), a window
% (e.g. -2:2), and returns a matrix (nWin x nEvent) where each event has
% the window applied.
%
% This can be useful for looking at epochs aligned at events in time series
% data.
%
%
% Inputs:
%   events
%       Vector array where each element is to be surrounded by 'window'.
%
%   window
%       A vector to be applied to each element in 'events'
%
%   removeLessThanEqZero (default true)
%       Flag if true, will remove columns that have any points less than
%       equal to 0. This would be used if you want valid matlab indices.
%
% Output
%   A matrix of size nWindow x nEvents where window is applied to each
%   value in events.
%
%
% For example: GetWindowAroundEvents([10 16 200 5030], -2:2) will
% return the following matrix:
% outputInds =
%
%            8          14         198        5028
%            9          15         199        5029
%           10          16         200        5030
%           11          17         201        5031
%           12          18         202        5032
%
% Example 2: Boolean events (do not remove columns with <= 0 points)
% removeLessThanEqZero = 0;
% outputInds = GetWindowAroundEvents([ 1 0 1 0 1 0 0 0 0 0 0 0 1 ], [-2:2 5], removeLessThanEqZero)
% Function Output:
%
% outputInds =
%
%     -1     1     3    11
%      0     2     4    12
%      1     3     5    13
%      2     4     6    14
%      3     5     7    15
%      6     8    10    18
%
% History:
%   2016   Copyright Tommy Hosman, Brown University. All Rights Reserved
%--------------------------------------------------------------------------



if ~exist('removeLessThanEqZero','var') || isempty(removeLessThanEqZero)
    removeLessThanEqZero = 1; % Removes columns with values <= zero
end



%% Handle logical input

% If we sent in logicals or only 1's and 0's, then perform a find on these events.
if islogical( events ) || ( numel(events) == sum( events == 0 | events== 1 ) )
    events = find(events);
end

%% Handle an empty window
if isempty(window)
    outputInds = events;
    return
end


%% Handle events shape
% event shape is 1 x nEvents
if ~isrow(events)
    events = events(:)'; % Want events to be a row
end


%% Handle window shape
% repWindow shape is [nSteps x nEvents]
if isvector(window)
    % Force window to be a column and repeat for each event
    repWindow = repmat(window(:),1, length(events));

elseif ismatrix(window) && size(window,2) == length(events)
    % If window is a matrix it must already be the right shape 
    % [nSteps x nEvents]
    repWindow = window;
else
    error('Unexpected window input. Expecting a vector or matrix')
end


%% Apply the repeated window to each event
repEvents  = repmat(events,size(repWindow,1),1);
outputInds = repEvents + repWindow;




%% Remove indices that are less than or equal to zero if requested
if removeLessThanEqZero
    removeZeroInds = any(outputInds <= 0);
    if any(removeZeroInds(:))
        fprintf('Removing %d columns that are <=0\n', sum(removeZeroInds(:)));
        outputInds( :, removeZeroInds) = [];
    end
end

end