function colors = cmap(nColors, colorStyle, lightColorThresh)
% Returns colors from a particular colorStyle that are not too close to
% white (see lightColorThresh).
% 
% Uses brewermap for creating and referencing different colormaps, but will
% discard lighter colors (where all rgb > lightColorThresh).
% 
% See also brewermap, brewermap_view
%
%--------------------------------------------------------------------------
% History:
%   2021.06   Copyright Tommy Hosman, Brown University. All Rights Reserved
%--------------------------------------------------------------------------


persistent isBrewerOnPath

% Set threshold so our colors are not too close to white.
if nargin < 1
    nColors = [];
end
if nargin < 2
    colorStyle = 'spectral';
end
if nargin < 3
    lightColorThresh = 0.55; %0.525;
end


if isempty(isBrewerOnPath)
    isBrewerOnPath = ~isempty(which('brewermap'));
end
if ~isBrewerOnPath
    error('This function requires brewermap.m which can be found online. Alternatively, ''lines'' should work in place of cmap.');
end
    

% Handle special cases
isSpecial = ismember( lower(colorStyle), {'paired', 'set1', 'set2', 'set3'});
if isSpecial
    colors = brewermap(nColors,colorStyle);
    return;
end
    
%% Get 4x of the requested number of colors
nOver  = 8*nColors;
colors = brewermap(nOver,colorStyle);


%% Exclude the colors that are too light
badColors = all(colors>lightColorThresh,2);

if all(badColors)
    warning('All colors are over threshold. Not excluding any.')
    colors = brewermap(nColors,colorStyle);
    return
end


colors(badColors,:) = [];
newN = size(colors,1);

% If nColors is empty, use all good colors
if isempty(nColors) % empty nColors is used if setting colormap
    nColors = newN;
end

%% Compute indexing into new non-light colors
if strcmpi(colorStyle, 'spectral')
    clrInds = round( linspace(1,newN,nColors) );
else
    
    % Split in half (in cases of rdbu)
    hlf = floor(nColors/2);
    hlfNew = round(newN/2);
    if rem(nColors,2) % If odd
        firstHalf = round( linspace(1,hlfNew,hlf) );
        secondHalf = round( linspace(hlfNew+1,newN,hlf+1) );

    else
        if hlf == 1
            firstHalf = 1;
        else
            firstHalf = round( linspace(1,hlfNew,hlf) );
        end
        secondHalf = round( linspace(hlfNew+1,newN,hlf) );

    end
    clrInds = [firstHalf secondHalf];
end


%% ReIndex into colors
colors = colors(clrInds,:);

if size(colors,1) ~= nColors
    error('%d colors were requested but we found %d colors... We messed up.', nColors, size(colors,1));
end

end
