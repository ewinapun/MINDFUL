function [cnt, inds, outEdges] = PolarHist(ang, N, isDeg, logicalInds)
% [cnt, inds, edges] = PolarHist(ang, N, isDeg)
% 
% This function takes a vector of angles (assumed to be in radians) and
% a number of bins and computes the polar histogram for equally separated
% bins around a circle. 
% 
% Additionally, it provides a logical index into ang for each bin in the
% form of an N x length(ang) logical matrix.
% 
% Inputs:
%   ang         (required) A vector of angles (assumed in radians, unless isDeg is set).
%   N           (8)     Number of bins in the histogram.
%   isDeg       (false) If the ang vector is in degrees.
%   logicalInds (false) If true, output inds will be a one hot vector for each angle passed.
%
% Outputs:
%   cnt         The number of angles in each bin.
%   inds        The index into which bin each ang sample falls into
%               Or
%               An N x length(ang) logical matrix set where an angle in ang
%               is inside the bin.
%               
%   edges       The edges in radians for each bin.
% 
% 
% Written by: Tommy Hosman 12/21/2017

%% Handle Inputs
if ~exist('isDeg','var') || isempty(isDeg)
    isDeg = 0;
end
if ~exist('logicalInds','var') || isempty(logicalInds)
    logicalInds = 0;
end
if ~exist('N','var') || isempty(N)
    N = 8;
end

isCart = 0;
if all(size(ang) >1)
    isCart = 1;
    % Assume cartesian coords
    if size(ang,1) == 2
        ang = ang';
    end
    ang = atan2(ang(:,2), ang(:,1));    
end

% Convert to radians. 
% But don't do it if we converted to from catesian to rads already (isCart).
if isDeg && ~isCart
    ang = deg2rad(ang);
end


%% Initialize outputs, edges   

inds     = false(N,length(ang));
outEdges = zeros(N,2);
edges    = linspace(0,2*pi,N*2+1);


%% Make sure ang is between 0 and 2pi and is a row vector

% Force ang to be a row vector
ang = ang(:)';

% Force all angles to be between 0 and 2pi
calcAng = mod(ang,2*pi);


%% Calc Polar Histogram
% edges = linspace(-pi/N, 2*pi-pi/N, 2*N+1);%(0:N)*2*pi/N - pi/N;
% [cnt, outEdges, angInds] = histcounts(rem(calcAng,edges(end)),edges);

% Set up the bins
outEdges(1,:)     = [edges(end-1) edges(2)];               % Set up the border between 0 and 2pi
outEdges(2:end,:) = [edges(2:2:end-2)' edges(4:2:end)'];   % handle all other angles

% Calculate the logical indexing
inds(1,:)     = calcAng > outEdges(1,1) | calcAng <= outEdges(1,2);         % Handle the border case between 0 and 2pi
inds(2:end,:) = calcAng > outEdges(2:end,1) & calcAng <= outEdges(2:end,2); % Handle all other angles

% Calculate the hist count
cnt = sum(inds,2);


if isDeg
    outEdges = rad2deg(outEdges);
end
    
if ~logicalInds
    [inds,~] = find(inds); % row index is the bin location
end

end