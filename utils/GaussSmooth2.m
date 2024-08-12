function data = GaussSmooth2( data, standardDev_ms, dim, samplePeriod_ms, alpha )
% Smooths data with a gaussian kernel with a standard deviation of
% standardDev_ms along dimension dim.
% 
% standardDev_ms == 0 will just return the original data.
% 
% alpha is a measure of the density of the kernel. Larger alpha has better
% resolution but takes longer for filtering.
% 
% Note, it is faster to smooth the first dim (dim==1), but faster to
% specify the dimension than to transpose. 
% 
% History:
%   2018   Update to use filter instead of conv which is faster
%   2016   Copyright Tommy Hosman, Brown University. All Rights Reserved
%--------------------------------------------------------------------------


if ~exist('standardDev_ms','var') || isempty(standardDev_ms)
    standardDev_ms  = 40; % ms
end

% Special case checking
if standardDev_ms <= 0
    % Assume that no smoothing was actually requested
    data = data;
    return;
end


defaultDim = 1;
if ~exist('dim','var')
    % What dimension do we apply the kernel to?
    % Set to defaultDim later in the code.
    dim  = []; 
end

if ~exist('alpha','var') || isempty(alpha)
    % Higher alpha has longer gaussian tails, but takes longer. 3 is a
    % reasonable tradeoff
    alpha = 3;
end

if ~exist('samplePeriod_ms','var') || isempty(samplePeriod_ms)
    samplePeriod_ms = 20; % ms
end

if standardDev_ms < samplePeriod_ms
    error('Your standardDev_ms (%.3f) must be >= samplePeriod_ms (%.3f)', standardDev_ms, samplePeriod_ms)
end

% Make sure samplePeriod_ms is >= 1
if samplePeriod_ms < 1
    mult = 1/samplePeriod_ms;
    samplePeriod_ms = samplePeriod_ms*mult;
    standardDev_ms = standardDev_ms*mult;
end


dataDims = size(data);


% If a row is passed and we are running with the default dim (dim isempty).
% Set defaultDim to be 2. To smooth along the row.
useDefaultDim = isempty(dim);
if isrow(data) && useDefaultDim
    defaultDim = 2;
end

if useDefaultDim
    dim = defaultDim;
end
    



%% Build gaussian filter

% See https://www.mathworks.com/help/signal/ref/gausswin.html for details.
% alpha       = (win-1)/(2*standardDev_ms);
winSize     = round( (alpha * standardDev_ms * 2 ) + 1 );
gaussFilter = gausswin(winSize, alpha);


% Down sample gaussian filter to our sample rate
inds = 1:samplePeriod_ms:winSize; % sampled points of our gaussian curve.
convGaussFilter = gaussFilter(inds) / sum(gaussFilter(inds)); % Normalize.
convWin = length(convGaussFilter);
halfWin = floor(convWin/2);

% If odd, make even
if bitand(convWin,1)
    convWin = convWin-1;
end



%% Create index variables and filter

% Variables for indexing
% [padIndsFront, padIndsBack, padIndsRmv] = deal(arrayfun(@(d) 1:d, dataDims, 'UniformOutput', false));
[padIndsFront, padIndsBack, padIndsRmv] = deal(repmat({':'}, 1, length(dataDims)));
padIndsFront{dim}   = 1:halfWin;
padIndsBack{dim}    = dataDims(dim)-halfWin+1:dataDims(dim);
padIndsRmv{dim}     = convWin+1:dataDims(dim)+convWin;


% Padded variable for filtering
tempVect = cat(dim,flip(data(padIndsFront{:}),dim),  data, flip(data(padIndsBack{:}),dim)); % Pad
data   = filter(convGaussFilter, 1, tempVect, [], dim); % Filter
data   = data(padIndsRmv{:}); % Remove the padded data






end