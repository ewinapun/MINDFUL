function inds = RowColon( startStops, vectorize )
% inds = RowColon( startStops, vectorize )
% 
% Creates a [cell] array of indices given a matrix of [n x 2] start stop inds.
% 
% 
% If the second column (the stop index) is lower than the first column,
% inds will decrease, i.e. [startIndex:-1:stopIndex].
% 
% 
% Input:
%   startStops 
%       [n x 2] maxtrix.
%       A matrix whose columns represent the respective start and stops inds.
%       For example, [2 5; 3 6; 10 100];
%   
%   vectorize (default true)
%       If true, inds is a [1 x N] vector of all concatenated inds.
% 
% Output:
%   inds
%       [1 x N] vector of the inds if vectorize is true.
%       Or
%       [n x 1] cell array (when vectorize is false)
% 
% 
% Example: Output a cell array of inds from a 3 x 2 start stop matrix, 'a'
% 
% a = [2 5; 3 6; 10 100];
% b = RowColon( a, 0 )
%
% b =
% 
%   3×1 cell array
% 
%     {1×4  double}
%     {1×4  double}
%     {1×91 double}
% 
% 
% Example: Output vector with a decreasing start-stop set.
% 
% a = [5 10; 2 -2];
% b = RowColon( a )
% 
% b =
%      5     6     7     8     9    10     2     1     0    -1    -2
%
%--------------------------------------------------------------------------
% History:
%   2019   Copyright Tommy Hosman, All Rights Reserved
%--------------------------------------------------------------------------

%% Handle inputs
    if nargin == 1
        vectorize = true;
    end

    [nRows, nCols] = size(startStops);
    
    % Error check and handle possible transposed input
    if nCols ~= 2
        if nRows == 2
            % Suspected transpose of input
            warning('Suspected mis-formatted matrix of size [%d x %d]. %s expects start stop index matrix of size [N x 2]. Transposing.', ...
                nRows, nCols, mfilename);
            startStops = startStops';
            nRows = nCols;
        else
            error('Suspected mis-formatted matrix of size [%d x %d]. %s expects start stop index matrix of size [N x 2].', ...
                nRows, nCols, mfilename);
        end
    end
    
    
%% Map start stops

    % Create a cell array of expanded inds
    inds = cell(nRows,1);
    for ii = 1:nRows
        if startStops(ii,1) < startStops(ii,2)
            % Increment indexing, e.g. [0, 1, 2]
            inds{ii} = startStops(ii,1):startStops(ii,2);
        else
            % Decrement indexing, e.g. [0, -1, -2]
            inds{ii} = startStops(ii,1):-1:startStops(ii,2);
        end
    end

    
    if vectorize
        inds = cat(2,inds{:});
    end

end