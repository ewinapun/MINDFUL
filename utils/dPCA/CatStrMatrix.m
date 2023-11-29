function catStr = CatStrMatrix( m, betweenVal, padStr,varargin)
% catStr = CatStrMatrix( m, betweenVal, padStr, padStrArgs)
% 
% Concatenates cell str by row adding betweenVal between each element.
% 
% betweenVal: str that separates each concatenated element
% padStr: pads elements in m. Passed as 2nd input to pad.m, can be a
% positive integer, true, or char side, e.g. 'left'
% 
% Note: If passed betweenVal that is empty, betweenVal uses the default.
% If you want to concatenate strs with no space, use strcat.

if nargin < 2 || isempty(betweenVal)
    betweenVal = {' | '};
end
if nargin < 3
    padStr = false;
end

[r,c] = size(m);

if iscategorical(m)
    m = cellstr(m);
end

isEmpty = cellfun('isempty', m);

if padStr
    if ischar(padStr) || (isnumeric(padStr) && padStr>0)
        m = pad(m, padStr, varargin{:});
    elseif islogical( padStr ) % True
        m = pad(m, 'both');
    end
end

%%
if ~iscell(betweenVal)
    betweenVal = {betweenVal};
end
if size(betweenVal) ~= size(m);
    betweenVal = repmat(betweenVal,r,1);
end

betweenVal(isEmpty) = {''};
%%
mm = cell(c+1,1); % Might need to handle odd/even c
for ii = 1:c
    mi = (ii*2)-1;
    % Entry
    mm{mi} = m(:,ii);
    
    % Set empty
    emptyEntries = cellfun('isempty', mm{mi});
    
    % Seperator
    if ii < c
        mm{mi+1} = cell( r, 1 );
        mm{mi+1}(~emptyEntries) = betweenVal(~emptyEntries);
    end
end
%%
catStr = strcat(mm{:});
end