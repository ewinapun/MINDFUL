function lagX = getlagX(X, startStop, lag)
% Define the function getlagX which needs three inputs: matrix X, matrix startStop and lag value
% Note that it will pad the initial part of each specified segment with zeros due to the lag.

% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

%     check input shape
    if size(X,1)==1 || size(X,2)==1
        X = X(:);
    end
    if size(X,1) < size(X,2)
        X = X';
    end
    % run loop over the rows of startStop (block based)
    if isempty(startStop)
        startStop = [1 size(X, 1)];
    else
        if size(startStop,2) ~= 2
            error('Input startStop should contain start and stop for each row [n x 2].')
        end
    end

    if ~isscalar(lag)
        error('Input lag should be a scalar.')
    end
    % Create a zero matrix, lagX, with the same size as X
    lagX = zeros(size(X));

    for t = 1:size(startStop,1)
        % assign the first element of row as start
        start = startStop(t,1);

        % assign the second element of row as stop
        stop = startStop(t,2);
        
        % set the relevant entries in lagX using the start and stop indices from startStop and the specified lag
        % This assigns, in range starting from start index + lag upto stop index, the values of X starting from start upto stop - lag. 
        lagX(start+lag:stop,:) = X(start:stop-lag,:);
    end
end