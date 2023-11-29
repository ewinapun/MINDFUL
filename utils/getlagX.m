%% return a lagged X with multiple startStop
function lagX = getlagX(X, startStop, lag)
    % pad with zeros for lag
    lagX = zeros(size(X));
    for t = 1:length(startStop)
        start = startStop(t,1);
        stop = startStop(t,2);
        lagX(start+lag:stop,:)= X(start:stop-lag,:);
    end
end