function [data, bgzscoreExtra] = BGzscoreNew( data, opt )
% [data, bgzscoreExtra] = BGzscoreNew( data, opt )
% 
% Function that zscores the data
%
%
% Inputs
%   data [nTime x nFeat]
%       A data matrix of time x features. z-scoring is performed on the
%       first dimension (time).
%
%   opt (see dfOpt for default options)
%       Structure to specify options which override the default options.
%
% Outputs
%   data
%       z-scored data
%
%   bgzscoreExtra
%       structure that specifies additional information about the zscore.
%       .opt - what options were used
%       .avgData - The average for each feature
%       .normData - The norm for each feature (e.g. std)
%
%
% Copyright Tommy Hosman, Brown University. All Rights Reserved
%--------------------------------------------------------------------------


%% Default options

% Z score
dfOpt.avgSubEnable      = true; % Set true to subtract the avg from the data
dfOpt.normDivideEnable  = true; % Set true to divide by the calculated norm
dfOpt.avgFun            = 'movmean'; % Function name or handle, e.g. 'mean' or @mean
dfOpt.normFun           = 'movstd'; % Function name or handle, e.g. @std. More robust methods: @mad, @iqr
dfOpt.avgFunArgs        = {'omitnan'};
dfOpt.normFunArgs       = {'omitnan'};


% Moving avg and norm params are used when avgFun, normFun are mov<fun>, e.g. movmean.
dfOpt.isMoving  = 'auto'; % If true, we use the avgWin_s, normWin_s params. If auto, then look if avgFun,normFun have 'mov' in their name.
dfOpt.avgWin_s  = 180;    % seconds to average over if using a rolling average function
dfOpt.normWin_s = 180;    % second to compute std over if using a rolling norm function
dfOpt.winFs     = 50;     % Hz. Converts avgWin_s and normWin_s to steps
dfOpt.trailingMoving = false; % If true, use trailing moving average [K 0]. If false, use centered moving average [K/2 K/2].


% Ignore steps when estimating mean/norm
dfOpt.ignoreSteps = false(0); % Logical matrix same size as data to ignore mean, norm estimates
dfOpt.fillIgnored = 'linear'; % How to set ignored means and variances. See fillmissing


if nargin < 2, opt = []; end
opt = MergeBintoA(dfOpt, opt);

if any(isfield(opt, {'norm','interpEstsOnNoise'}))
    error('Updated opt struct. Please see this functions dfOpt for the updated opt struct.')
end



%% Z-Score the data

[data, avgData, normData] = ZScoreData( data, opt );


%% Return additional information
if nargout > 1
    bgzscoreExtra.opt = opt;
    bgzscoreExtra.avgData = avgData;
    bgzscoreExtra.normData = normData;
end

end



function [data, avgData, normData] = ZScoreData( data, opt )
% Zscores the data with some norm checks
%
% If interested in fancier methods, look at normalize
% e.g. normalize(data, 'zscore', 'robust')
%
%
% Another method used for zscoring (not implemented here)
%
% min-max normalizing, while discounting low (firing rate) features with the offset
% See Kaufman, M. T. et al. The largest response component in the motor cortex reflects movement timing but not movement type. eNeuro 3, (2016)


[avgData, normData] = deal([]); % Init outputs


modData = SetIgnoreStepsToNan(data,opt);

opt = AddMoveWindowArgs(opt);


%--------------------------------------------------------------------------
% Compute average and norm based on the respective functions



% Subtract the estimated average
if opt.avgSubEnable
    if ischar(opt.avgFun)
        opt.avgFun = str2func(opt.avgFun);
    end
    
    % Estimate average
    avgData = opt.avgFun( modData, opt.avgFunArgs{:} );
    avgData = FillNansWithEdgeCase(avgData, opt.fillIgnored);
    
    data = bsxfun(@minus, data, avgData);
end


% Divide by the estimated norm
if opt.normDivideEnable
    if ischar(opt.normFun)
        opt.normFun = str2func(opt.normFun);
    end
    
    % Estimate norm
    normData = opt.normFun( modData, opt.normFunArgs{:} );
    normData = FillNansWithEdgeCase(normData, opt.fillIgnored);
    normData(normData<1e-5) = 1; % Don't want to divide by (near) 0
    
    data = bsxfun(@rdivide, data, normData);
end



end

function modData = SetIgnoreStepsToNan(modData,opt)
    % Set ignore steps to nan in modData to ignore outliers

    isIgnoreSameSize  = isequal(size(modData), size(opt.ignoreSteps));
    isIgnoreSameSteps = isequal(size(modData,1), size(opt.ignoreSteps,1));
    isIgnore = isIgnoreSameSize || isIgnoreSameSteps;
    if isIgnore
        if isIgnoreSameSize
            modData(opt.ignoreSteps) = nan;
        elseif isIgnoreSameSteps
            modData(opt.ignoreSteps,:) = nan;
        end
    end

end

function data = FillNansWithEdgeCase(data, interpMethod)
    
    % Check for outliers/nans edge case
    outliers = isnan(data);
    if ~any(outliers(:))
        return
    end
    
    % If outliers on first or last time step, some interpolation methods
    % will misbehave. Append the median (per row, or feature) to the
    % front/back of data.
    % Possibly could append the first/last non-nan instance?
    
    isOutlierOnEdge = any(outliers(1,:)) || any(outliers(end,:));
    isSuseptableInterpMethod = ismember(interpMethod, {'linear', 'pchip','spline', 'previous', 'next'});
    appendMedian = isOutlierOnEdge && isSuseptableInterpMethod;
    
    if appendMedian
        medianPerFeat = nanmedian(data);
        data = cat(1, medianPerFeat, data, medianPerFeat);    
    end
    
    % Fill nans
    data = fillmissing(data, interpMethod);
    
    
    if appendMedian
        data([1 end], :) = [];
    end
    
end

function opt = AddMoveWindowArgs(opt)
    if opt.isMoving
        % Check if true, or if auto look for mov keywords
        if islogical(opt.isMoving) || isnumeric(opt.isMoving)
            isMoveAvg  = true;
            isMoveNorm = true;
        elseif ischar(opt.isMoving) && strcmp(opt.isMoving,'auto')
            % Auto, look at avg and norm function names
            if ischar( opt.avgFun )
                avgFunStr = opt.avgFun;
            else
                % Assume a function
                avgFunStr = func2str( opt.avgFun );
            end
            
            if ischar( opt.normFun )
                normFunStr = opt.normFun;
            else
                % Assume a function
                normFunStr = func2str( opt.normFun );
            end
            
            moveKeywords = {'mov'};
            isMoveAvg  = contains(avgFunStr, moveKeywords);
            isMoveNorm = contains(normFunStr, moveKeywords);
        else
            opt.isMoving
            error('Unexpected opt.isMoving.')
        end

        % If moving functions, add the window to the start of the function
        % arguments.
        if isMoveAvg
            avgWin = opt.avgWin_s*opt.winFs;
            if opt.trailingMoving
                % trailing moving average
                opt.avgFunArgs  = [{[avgWin 0]}, opt.avgFunArgs];
            else
                % centered moving average
                opt.avgFunArgs  = [{avgWin}, opt.avgFunArgs];
            end
        end

        if isMoveNorm
            normWin = opt.normWin_s*opt.winFs;
            if opt.trailingMoving
                % trailing moving average
                opt.normFunArgs = [{[normWin 0]}, opt.normFunArgs];
            else
                % centered moving average
                opt.normFunArgs  = [{normWin}, opt.normFunArgs];
            end
        end
    end

end