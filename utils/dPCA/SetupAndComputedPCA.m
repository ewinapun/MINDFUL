function [dpcaOut, firingRates, firingRatesAverage] = SetupAndComputedPCA(feat, labelGrp, p)

% Given features [nTrial x nTime x nFeatures] and labelGrp [nTrial x nGroups] 
% (sets of unique label index per group), this function arranges the data
% in preparation for dPCA calls.

% History:
%   2016   Copyright Tommy Hosman, Brown University. All Rights Reserved
%--------------------------------------------------------------------------

HandleDependencies();
if nargin < 1 && nargout == 0
    return
end
    
if nargin < 3
    p = [];
end

% Default params
dfParams.plotFun = @dpca_plot_standard; % @dpca_plot_default;
dfParams.plotPCA = false;
dfParams.pltTime = []; % time value to be plotted for each time step
dfParams.marginalizationNames = {'Stimulus', 'Condition-independent', 'Interaction'};
dfParams.marginalizationColours = [];
dfParams.conditionNames = {}; % Stim conditions
dfParams.lambda = 4.9879e-04; %'optimal'; % 'optimal'
dfParams.Cnoise = []; % 'optimal'
dfParams.centerData = 1;
dfParams.titleStr = '';
dfParams.runDPCA = true;
dfParams.returnMargs = true;
dfParams.printTrialConditionCounts = true;
dfParams.nDPCs = 20;
dfParams.nanOutsideOfTrial = [];
dfParams.breakTimeOut = []; % Marg number here
dfParams.dpcFigNum = [];
dfParams.pcaFigNum = [];
dfParams.sessionIndex = [];
p = MergeBintoA(dfParams, p);

p.nDPCs = min(p.nDPCs, size(feat,3));

%% Organize data by conditions

[firingRates, firingRatesAverage, dpcaOut, p] = OrganizeFeaturesByCondition(feat, labelGrp, p);




%% Compute dPCA
pca_h = gobjects(0);
dpca_h = gobjects(0);
if p.runDPCA
    
    % If the number of name entries does not match the number of
    % marginalizations, and they are the default marginalization entries,
    % clear p.marginalizationNames.
    if length(p.marginalizationNames) < ndims(firingRatesAverage)-1 && ...
            isequal( p.marginalizationNames, dfParams.marginalizationNames )
        p.marginalizationNames = {};
    end
    
%     [accuracy, brier] = dpca_classificationAccuracy(firingRatesAverage, firingRates, numTrials, 'simultaneous', 1, ...
%         'combinedParams', combinedParams, 'lambda', 0.001);
    if strcmp(p.lambda, 'optimal')
        [optimalLambda, optimalLambdas] = dpca_optimizeLambda(firingRatesAverage, ...
                    firingRates, dpcaOut.numTrials, ...
                    'simultaneous', 1, ...
                    'combinedParams', dpcaOut.combinedParams,  ...
                    'display', 'no', ...
                    'lambdas',        1e-07 * 1.5.^[0:3:35]);
        p.lambda = optimalLambda;
    end
%%
    if isempty(p.Cnoise)
        Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
            firingRates, dpcaOut.numTrials, 'simultaneous', 1);
    else
        Cnoise = p.Cnoise;
    end

    [W,V,whichMarg] = dpca(firingRatesAverage, p.nDPCs, ...
        'combinedParams', dpcaOut.combinedParams, ...
        'Cnoise', Cnoise, ...
        'lambda', p.lambda, ...
        'centerData', p.centerData, ...
        'marginalizationNames', p.marginalizationNames);


    explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
                    'combinedParams', dpcaOut.combinedParams, ...
                    'X_trial', firingRates, ...
                    'numOfTrials', dpcaOut.numTrials, ...
                    'Cnoise', Cnoise, ...
                    'centerData', p.centerData);


    %% Plot dPCA
    if ~isempty(p.plotFun)
        if isempty(p.marginalizationColours)
            dims = size(firingRatesAverage);
            nMargColors = length(p.marginalizationNames); % max(dims(2:end-1));
            p.marginalizationColours = lines(nMargColors);
%             nCondColors = length(combinedParams);
%             p.conditionColours = cmap(nCondColors, 'spectral');
        end
        

        [~, dpca_h] = dpca_plot(firingRatesAverage, W, V, p.plotFun, ...
            'explainedVar', explVar, ...
            'marginalizationNames', p.marginalizationNames, ...
            'marginalizationColours', p.marginalizationColours, ...
            'whichMarg', whichMarg,                 ...
            'timeMarginalization', dpcaOut.timeMarginalization,           ...
            'legendSubplot', 16, ...
            'time', p.pltTime, ...
            'titleStr', p.titleStr, ...
            'conditionNames', p.conditionNames, ...
            'fig',          p.dpcFigNum);
        
        if p.plotPCA
            %%
            pca_h = dpca_perMarginalization(firingRatesAverage, p.plotFun, ...
                        'combinedParams', dpcaOut.combinedParams, ...
                        'marginalizationNames',p.marginalizationNames, ...
                        'timeMarginalization', dpcaOut.timeMarginalization, ...
                        'time', p.pltTime, ...
                        'titleStr', p.titleStr, ...
                        'conditionNames', p.conditionNames, ...
                        'centerData', p.centerData, ...
                        'fig',          p.pcaFigNum);
        end
        

    end
    


    dpcaOut.W = W;
    dpcaOut.V = V;

    dpcaOut.Cnoise    = Cnoise;
    dpcaOut.whichMarg = whichMarg;
    dpcaOut.explVar   = explVar;
end

if p.returnMargs
    [Xmargs, margNums] = dpca_marginalize(firingRatesAverage, 'combinedParams', dpcaOut.combinedParams, ...
                        'ifFull', 'yes', ...
                        'ifFlat', 'no');

    dpcaOut.margs.Xmargs = Xmargs;
    dpcaOut.margs.margNums = margNums;
end

dpcaOut.params = p;
dpcaOut.handles.dpca = dpca_h;
dpcaOut.handles.pca = pca_h;
end

function HandleDependencies()
% See if dpca repo is on path. 
% 
% If not, try to find them in the analysis repo and add them to the matlab
% path:  Analysis\Shareware\dPCA-master\matlab\
    if isempty(which('dpca_getNoiseCovariance'))
        fileLoc = mfilename('fullpath');
        userInd = regexp(fileLoc, 'UserAnalysis', 'once');
        if isempty(userInd)
            error('Cannot find the dpca repo on path. Expected this file to be in a subfolder of the UserAnalysis folder. Please add dPCA-master in Shareware to path.')
        end
        rootAnalysisRepo = fileLoc(1:userInd-1);
        dpcaSharewareDir = fullfile(rootAnalysisRepo, 'Shareware', 'dPCA-master', 'matlab');

        if ~isfolder(dpcaSharewareDir)
            error('Cannot find the dpca repo on path. Expected dPCA repo to be in Analysis repo''s Shareware folder. Please add dPCA-master in Shareware to path.')
        end

%         if contains(userpath, 'Tommy')
%             ag 'C:\Users\Tommy\Documents\MATLAB\libraries\dPCA-master\matlab';
%         else
        fprintf('Adding %s and subfolders to path\n', dpcaSharewareDir);
        addpath(genpath(dpcaSharewareDir));
%         end
    end
end

function [firingRates, firingRatesAverage, dpcaOut, p] = OrganizeFeaturesByCondition(feat, labelGrp, p)
    [nTrials, nDims] = size(labelGrp);
cvIn = {};
nGroups = {};
grpClassNames = {};
grps = zeros(size(labelGrp));
% Find the number of groups for each label dimension
for ii = 1:nDims
    if iscategorical(labelGrp(:,ii))
        tmpLbl = removecats(labelGrp(:,ii));
        
    else
        tmpLbl = labelGrp(:,ii);
    end
    [grps(:,ii), uGrp]  = grp2idx(tmpLbl);
    nGroups{ii} = length( uGrp );
    cvIn{ii}    = 1:length( uGrp );
    grpClassNames{ii} = uGrp;
end

if isempty(p.conditionNames) && iscategorical(labelGrp)
    p.conditionNames = cat(1,grpClassNames{:});
end
%%
% Create a matrix of all combinations of each group, dim
cv = combvec(cvIn{:})';
tmpNumTrials = zeros(nGroups{:});

% Simple expand into cv group name sets
grpSet = cell(size(cv));
for ii = 1:size(cv,2)
    grpSet(:,ii) = cat(1,grpClassNames{ii}(cv(:,ii)));
end
grpSet = CatStrMatrix(grpSet,[],true);
% Count the number of trials in each group, dim combination
for ii = 1:size(cv,1)
    % Find the trials whose label inds match this combination of label inds
    selTrials = ismember(grps, cv(ii,:), 'rows');
    dimInd2Cell = num2cell(cv(ii,:));
    tmpNumTrials( dimInd2Cell{:} ) = sum(selTrials);
end

% This is used so we can init our firingRates matrix
maxTrials = max(tmpNumTrials(:));


%% Organize aligned features groups tasks

ndPerm = [3, 2, 1];
[nTr, nTime, nFeat] = size(feat);
firingRates = nan(nFeat, nGroups{:}, nTime, maxTrials);
numTrials   = nan(nFeat, nGroups{:} );
trialMap    = nan(nGroups{:}, maxTrials );
toExclude   = cell(1,size(cv,2));
for ii = 1:size(cv,1)
    
    
    % Find the trials whose label inds match this combination of label inds
    selTrials  = find(ismember(grps, cv(ii,:), 'rows'));
    nSelTrials = length(selTrials);
    frTrInds   = 1:nSelTrials;
    
    if p.printTrialConditionCounts
        fprintf('%s %03d trials\n', grpSet{ii}, nSelTrials)
    end

    
    if ~isempty( p.nanOutsideOfTrial )
        for jj = 1:nSelTrials
            nanInds = p.nanOutsideOfTrial(selTrials(jj))+1:nTime;
            feat(selTrials(jj),nanInds,:) = nan;
        end
    end
    
    dimInd2Cell = num2cell(cv(ii,:));
    numTrials( :, dimInd2Cell{:} ) = nSelTrials;    
    firingRates(:,dimInd2Cell{:},:,frTrInds) = permute(feat(selTrials,:,:), ndPerm);
    trialMap(dimInd2Cell{:}, frTrInds) = selTrials;
    
    if ~nSelTrials
        for jj = 1:size(cv,2)
            toExclude{jj}(end+1) = cv(ii,jj);
        end
    end
end


%% Sanity check nans

nFRDim = ndims(firingRates);
if nFRDim == 4
    % -- For only stim (target direction) --
    % firingRates array has [N S T E] size; here we ignore the 1st dimension 
    % (neurons), i.e. we have the following parameters:
    %    1 - stimulus 
    %    2 - time
    % There is one pairwise interaction:
    %    [1 2] - stimulus/time interaction
    
    combinedParams = {{1, [1 2]}, {2}};  
    timeMarginalization = 2;
    [combinedParams, p] = BreakTimeOut(combinedParams,p, timeMarginalization);
    
%     firingRates = BootstrapData(firingRates);
    firingRatesAverage = nanmean(firingRates,nFRDim);
    %%
    

    %%
    [N,S,T,E] = size(firingRates);

    for n = 1:size(firingRates,1)
        for stmI = 1:size(firingRates,2)
            assert(isempty(find(isnan(firingRatesAverage(n,stmI,:)), 1)), 'Something is wrong!')
        end
    end

elseif nFRDim == 5
 
    % -- For stim/decision --
    % 
    % firingRates array has [N S D T E] size; here we ignore the 1st dimension 
    % (neurons), i.e. we have the following parameters:
    %    1 - stimulus 
    %    2 - decision
    %    3 - time
    % There are three pairwise interactions:
    %    [1 3] - stimulus/time interaction
    %    [2 3] - decision/time interaction
    %    [1 2] - stimulus/decision interaction
    % And one three-way interaction:
    %    [1 2 3] - rest
    if length(p.marginalizationNames) == 3
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}};
    else
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    end
    % 
    timeMarginalization = 3;
    [combinedParams, p] = BreakTimeOut(combinedParams,p,timeMarginalization);

%%

    firingRatesAverage = nanmean(firingRates,nFRDim);
    [N,S,D,T,E] = size(firingRates);
    for n = 1:size(firingRates,1)
        for stmI = 1:size(firingRates,2)
            for d = 1:size(firingRates,3)
                assert(isempty(find(isnan(firingRatesAverage(n,stmI,d,:)), 1)), 'Something is wrong!')
            end
        end
    end
elseif nFRDim == 6
    % -- For stim/decision --
    % 
    % firingRates array has [N S D1 D2 T E] size; here we ignore the 1st dimension 
    % (neurons), i.e. we have the following parameters:
    %    1 - stimulus 
    %    2 - decision1 (laterality)
    %    3 - decision2 (proximality)
    %    4 - time
    % There are three pairwise interactions:
    %    [1 4] - stimulus/time interaction
    %    [2 4] - decision1/time interaction
    %    [3 4] - decision2/time interaction
    %    [1 2 3] - stimulus/decision interaction
    % And one three-way interaction:
    %    [1 2 3 4] - rest
    if length(p.marginalizationNames) == 5
        combinedParams = {{1, [1 4]}, {2, [2 4]}, {3, [3 4]}, {4},...
                          {[1 2], [1 3], [2 3], [1 2 3], [1 2 3 4]}};
    else
        combinedParams = {{1, [1 4]}, {2, [2 4]}, {3, [3 4]}, {4}};
    end
    timeMarginalization = 4;




    firingRatesAverage = nanmean(firingRates,nFRDim);
    [N,S,D1, D2,T,E] = size(firingRates);
    for n = 1:size(firingRates,1)
        for stmI = 1:size(firingRates,2)
            for d1 = 1:size(firingRates,3)
                for d2 = 1:size(firingRates,4)
                    assert(isempty(find(isnan(firingRatesAverage(n,stmI,d1,d2,:)), 1)), 'Something is wrong!')
                end
            end
        end
    end
elseif nFRDim == 7
    % -- For stim/decision --
    % 
    % firingRates array has [N S D1 D2 D3 T E] size; here we ignore the 1st dimension 
    % (neurons), i.e. we have the following parameters:
    %    1 - stimulus 
    %    2 - decision1 (laterality)
    %    3 - decision2 (proximality)
    %    4 - decision2 (proximality)
    %    5 - time
    % There are three pairwise interactions:
    %    [1 4] - stimulus/time interaction
    %    [2 4] - decision1/time interaction
    %    [3 4] - decision2/time interaction
    %    [1 2 3] - stimulus/decision interaction
    % And one three-way interaction:
    %    [1 2 3 4 5] - rest
    nSubsets = nFRDim-2;
    
    %%
    
    S = GetSubsets(1:nSubsets);
    usedInds = [];
    combinedParams = {};
    for ii = 1:nSubsets
        if ii < nSubsets+1
            addInds = [ ii, 2*nSubsets + ii + 1];            
        elseif ii == nSubsets
            addInds = ii;            
        end
        combinedParams{ii} = S(addInds);        
        usedInds = [usedInds addInds];
    end
    combinedParams{end+1} = S(setdiff(1:length(S)-1, usedInds));
    combinedParams{end+1} = S(end);
    
    timeMarginalization = nSubsets;




    firingRatesAverage = nanmean(firingRates,nFRDim);
    [N,S,D1, D2,D3,T,E] = size(firingRates);
    for n = 1:size(firingRates,1)
        for stmI = 1:size(firingRates,2)
            for d1 = 1:size(firingRates,3)
                for d2 = 1:size(firingRates,4)
                    for d3 = 1:size(firingRates,5)
                    assert(isempty(find(isnan(firingRatesAverage(n,stmI,d1,d2,d3,:)), 1)), 'Something is wrong!')
                    end
                end
            end
        end
    end
        
else
    error('unsupported Firing rate dim size %d', nFRDim)
end



dpcaOut.numTrials = numTrials;
dpcaOut.combinedParams = combinedParams;
dpcaOut.timeMarginalization = timeMarginalization;
dpcaOut.trialMap = trialMap;


end



function f = BootstrapData(f)
    permInds = [4 3 2 1];
    f = permute(f,permInds);
    dims = size(f);
    f = reshape(f, [size(f,1:2), prod(size(f,3:length(dims)))]);
    nTrials = size(f,1);
    for ii = 1:size(f,3)
        f(:,:,ii) = bootstrp(nTrials,@nanmean, f(:,:,ii));
    end
    f = reshape(f, dims);
    f = ipermute(f,permInds);
end
function [combinedParams, p] = BreakTimeOut( combinedParams,p, timeMarginalization)

for ii = 1:length(p.breakTimeOut)
    bo = p.breakTimeOut(ii);
    combMarg = combinedParams{bo};
    timeInds = find(cellfun(@(x) any(ismember(x, timeMarginalization)), combMarg));
    if ~isempty(timeInds)
        timeMargs = combMarg(timeInds);
        combMarg(timeInds) = [];
        combinedParams{bo} = combMarg;
        combinedParams{end+1} = timeMargs;
        
        p.marginalizationNames{end+1} = [p.marginalizationNames{bo} ' + Time'];
    end
    
%     switch nDims
%         case 4
%             
%             combinedParams = {{1}, {2}, {[1 2]}};
%             
%             dpcaOut.params.marginalizationNames
%             p.marginalizationNames{end+1} = 'Gesture + Time';
%         case 5
%             if length(p.marginalizationNames) == 3
%                 combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}};
%             else
%                 combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
%             end
%     end
end



end

