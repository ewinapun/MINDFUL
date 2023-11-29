function [outMoveDir, uMoveDirs] = ChangeDirsToCategorical(moveDir, exactMatch)
% Map Directional vectors to categorical values.
% Input:
%     moveDir: [N x 2] direction vector
%     exactMatch: (default to True) angles must be within tolerance
%
% Output:
%     outMoveDir
%     uMoveDirs
%--------------------------------------------------------------------------
% History:
%   2019   Copyright Tommy Hosman, All Rights Reserved
%--------------------------------------------------------------------------

    if nargin < 2
        exactMatch = true;
    end
    
    tol = 1e-5;
    dirStrMap = {
            [0,             0], 'No Direction'  % No Direction is used below. If wording changes. Change it below as well!
            [0.5000         0], 'Right'
            [0.3500    0.3500], 'Up-Right'
            [0         0.5000], 'Up'
            [-0.35000    0.35], 'Up-Left'
            [-0.5000        0], 'Left'
            [-0.35      -0.35], 'Down-Left'         
            [0.00        -0.5], 'Down'
            [0.3500   -0.3500], 'Down-Right'
            };

    moveMapStr  = dirStrMap(:,2);
    moveMapVect = cat(1,dirStrMap{:,1});
    
%     if isa(moveDir, 'single')
%         moveMapVect = single(moveMapVect);
%     end
    
    
    %% Convert to angle
    moveMapRad = atan2(moveMapVect(:,2),moveMapVect(:,1));
    moveDirRad = atan2(moveDir(:,2),moveDir(:,1));
    
    
    %% Map to direction

    moveDirInd = SelectByClosest(moveDirRad,moveMapRad, exactMatch, tol);
    
    
    % Look for no direction (literal [0,0] vector, not by angle)
    % Set no-dir afterwards because 0,0 and 0.5,0 can be confused as 0 rad
    noDirInd = (ismember(moveMapStr,'No Direction'));
    isNoDir  = ismembertol(moveDir, moveMapVect(noDirInd,:), tol, 'ByRows',1);
    if any(isNoDir)
        moveDirInd(isNoDir) = 1;
    end
    
    
    % Check for not set inds (nans)
    nanDirInds = find(isnan(moveDirInd));
    if ~SanityCheck(moveDir, nanDirInds)
        keyboard; % NaNs found (see print in SanityCheck)
    end
    
    %% Map to direction for output
    uMoveDirInds = unique(moveDirInd(~isnan(moveDirInd)));
    outMoveDir   = categorical(moveDirInd(:), uMoveDirInds, moveMapStr(uMoveDirInds));
    
    uMoveDirs = moveMapVect(uMoveDirInds,:);
end

function isValid = SanityCheck(moveDir, nanDirInds)
    
    isValid = isempty(nanDirInds);
    if ~isValid
        invalidInfoStr = sprintf('\n\nDid not catch all move directions (%d out of %d)\nPausing (keyboard)...\n\n', length(nanDirInds), length(moveDir));
        fprintf(2,invalidInfoStr);
        for ii = 1:length(nanDirInds)
            fprintf(2, 'Trial index: %03d   Direction vector:  %20s\n', nanDirInds(ii), num2str(moveDir(nanDirInds(ii),:)))
        end
        fprintf(2,invalidInfoStr);
    end
end

function moveDirInd = SelectByClosest(moveDirRad,moveMapRad, exactMatch, tol)
% Use tolerance formula from ismembertol
% TOL*max(abs([A(:);B(:)]))
scaledTol = tol*max(abs([moveDirRad(:); moveMapRad(:)]));

moveMapRad(1) = []; % Remove the no direction (it gets confused with right)
[minV, minI] = min(abs(loc_circ_dist2(moveDirRad,moveMapRad)),[],2);
moveDirInd = minI + 1; % Shift back to no direction aware index

%% If require exact match, then scale based on 
if exactMatch
    moveDirInd(minV > scaledTol) = nan;
end

end



function r =  loc_circ_dist2(x,y)
% Local function copied to not require dependency...
% 
% r = circ_dist(alpha, beta)
%   All pairwise difference x_i-y_j around the circle computed efficiently.
%
% Input:
%     alpha       sample of linear random variable
%     beta       sample of linear random variable
%
% Output:
%     r       matrix with pairwise differences
%
% References:
%     Biostatistical Analysis, J. H. Zar, p. 651
%
% PHB 3/19/2009
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 2
  y = x;
end

if size(x,2)>size(x,1)
  x = x';
end

if size(y,2)>size(y,1)
  y = y';
end

r = angle(repmat(exp(1i*x),1,length(y)) ...
       ./ repmat(exp(1i*y'),length(x),1));
end
