function [cat_ND, cat_labels, events, info, cat_extra] = ConcatSavedSessionsData(path, p)

% inputs
%   path 
%       load extracted sessions from this path
%       expected to contain [session name]/[block number]
%       with /data, /info and /task
%
%   p
%    - structure to specify options which override the default options. 
%       - oneCalBlock
%           only load first block per session
%       - skipSessionName 
%           cell struct containing the session names to be skipped
%       - requestSessionNames 
%           cell struct containing the requested session names
%
% History:
%   2021.02.19   created by Ewina Pun
%   2022.08.25   last edit by Ewina Pun
%
%   Copyright Ewina Pun, All Rights Reserved

if nargin < 2
    p = []; 
end
opt.oneCalBlock = 0;
opt.skipSessionNames = [];
opt.requestSessionNames = [];
p = MergeBintoA(opt, p);

addpath(genpath(fullfile(path)))
UDIR = dir(path);
% Extract only those that are directories.
subFolders = UDIR([UDIR.isdir]);
subFolders = subFolders(~ismember({subFolders(:).name},{'.','..'}));
% exclude skip sessions
if ~isempty(p.skipSessionNames)
    subFolders = subFolders(~ismember({subFolders(:).name}, p.skipSessionNames));
end
% only load requested sessions
if ~isempty(p.requestSessionNames)
    subFolders = subFolders(ismember({subFolders(:).name}, p.requestSessionNames));
end
cat_ND      = [];
cat_labels  = [];
events      = {};
info        = {};
cat_extra   = {};
usedSessionList = {};
sessionStartInd         = 0;
events.trialsPerBlock = [];
events.trialsPerSession = [];
events.pointsPerSession = [];
events.excludeTrials = true(0);
events.sessionStartStop = [];
events.trialStartStop   = [];
events.blockStartStop   = [];
events.sessionNumberPerTrial = [];
events.sessionNumberPerBlock = [];
events.uniqueGoalStateLabels = [];
events.goalState = [];
events.goalPosition = [];
events.movedirVect = [];
events.gestGoalState = [];
usedBlocks = {};

for k = 1:length(subFolders)
    try
    fprintf('\n[%d/%d] Concatenating %s: ', k, length(subFolders), subFolders(k).name);
    usedSessionList = [usedSessionList; subFolders(k).name];
    SDIR = dir(fullfile(subFolders(k).folder,subFolders(k).name));
    % sort by ascending time/block
    SDIR = SDIR(~ismember({SDIR(:).name},{'.','..'}));
    [~,idx] = sort([SDIR.datenum]);
    SDIR = SDIR(idx);
    [neuralData, labels, tmpEvents, extra, blocknum] = concatBlocksData(SDIR, p);
    
    cat_ND = vertcat(cat_ND, neuralData);
    cat_labels = vertcat(cat_labels, labels);
    cat_extra = ConcatStruct(cat_extra, extra);
    usedBlocks = vertcat(usedBlocks, {blocknum});
    
    %Concatenate events data
    
    events.excludeTrials = vertcat(events.excludeTrials, tmpEvents.excludeTrials);
    events.movedirVect = vertcat(events.movedirVect, tmpEvents.movedirVect);
    events.gestGoalState = vertcat(events.gestGoalState, tmpEvents.gestGoalState);

%     events.goalPosition = vertcat(events.goalPosition, tmpEvents.goalPosition);
    
    events.trialStartStop = vertcat(events.trialStartStop, ...
        double(tmpEvents.trialStartStop) + sum(events.pointsPerSession));

    events.blockStartStop = vertcat(events.blockStartStop, ...
        double(tmpEvents.blockStartStop) + sum(events.pointsPerSession));

    events.sessionStartStop(k, 1) = sessionStartInd + 1;
    sessionStartInd = sessionStartInd + tmpEvents.pointsPerSession;
    events.sessionStartStop(k, 2) = sessionStartInd;
    
    events.sessionNumberPerTrial = vertcat(events.sessionNumberPerTrial, ...
        repmat(length(usedSessionList), size(tmpEvents.trialStartStop, 1), 1));

    events.sessionNumberPerBlock = vertcat(events.sessionNumberPerBlock, ...
        repmat(length(usedSessionList), size(tmpEvents.blockStartStop, 1), 1));
    
    events.trialsPerSession = [events.trialsPerSession; tmpEvents.trialsPerSession];
    events.trialsPerBlock = [events.trialsPerBlock; tmpEvents.trialsPerBlock];
    events.pointsPerSession = [events.pointsPerSession; tmpEvents.pointsPerSession];
    catch err2
        UnrollError(err2);
        return;
    end
end
info.usedSessionList = usedSessionList;
info.usedBlocks = usedBlocks;
info.sessionDates = [];
info.trialDay = [];
for i = 1:length(info.usedSessionList)
    [info.trialDay(i), sdate] = GetTrialDay(info.usedSessionList{i});
    info.sessionDates = [info.sessionDates;sdate];
end
fprintf('\nDone\n')
end


function [cat_ND, cat_labels, events, cat_extra, blocknum] = concatBlocksData(SDIR, p)
% concat all wanted blocks from the same session first
% dont care about extra for now
if ~exist('p','var')
    p =[];
end
opt.zscoreFeatures = 0;
p = MergeBintoA( opt, p );

cat_ND      = [];
cat_labels  = [];
cat_cursorVel  = [];
cat_extra   = {};
events      = {};
events.trialStartStop   = [];
events.blockStartStop   = [];
events.trialsPerSession = 0;
events.pointsPerSession = 0;
events.trialsPerBlock = [];
events.excludeTrials = true(0);
events.goalPosition = [];
events.movedirVect = [];
events.gestGoalState = [];
blocknum    = [];

clear goalState uniqueGoalStateLabels

if isfield(p,'selectBlock')
    blocks = p.selectBlock;
else
    blocks = 1:length(SDIR);
end

for b = blocks
    if ~SDIR(b).isdir
        continue
    end
    block_folder = SDIR(b).name;
    load(fullfile(SDIR(b).folder, block_folder, 'task.mat'))
    
    % extract only the wanted gametype, extract all available blocks if no
    % gametype specified
    if isfield(p, 'useGameTypes') && ~strcmp(p.useGameTypes, name)
        disp('Game type of requested blocks do not match found blocks')
        continue
    end
    
    extra = load(fullfile(SDIR(b).folder, block_folder, 'info.mat'));
    if isfield(extra,'angleErrorPerTrial')
        extra.meanAEPerBlock = nanmean(extra.angleErrorPerTrial);
    end
    
    load(fullfile(SDIR(b).folder, block_folder, 'data.mat'));
    
    % account for different names for neural data
    if ~exist('data','var')
        if exist('feats','var')
            data = feats;
        end
        if exist('nctx','var')
            data = nctx;
        end
        if exist('nctx','var') && exist('spikePower','var')
            data = [nctx spikePower];
        end
    end
    data = double(data);
    if p.zscoreFeatures
        data = BGzscoreNew(data, p);
%         data = zscore(data);
    end
    cat_ND = vertcat(cat_ND, data);
    cat_labels = vertcat(cat_labels, labels);
    cat_cursorVel = vertcat(cat_cursorVel, cursorVel);
    
    extra.normVal = {};
    if exist('normVal','var')
        extra.normVal = ConcatStruct(extra.normVal, normVal);
    end

    cat_extra = ConcatStruct(cat_extra, extra);

    clear data   
    
    %Concatenate events data
    if exist('excludeTrials','var')
        events.excludeTrials = vertcat(events.excludeTrials, excludeTrials);
    end
%     events.goalPosition = vertcat(events.goalPosition, goalPosition);    
    if exist('moveDirVect','var')
        events.movedirVect = vertcat(events.movedirVect, moveDirVect);    
    end
    events.trialStartStop = vertcat(events.trialStartStop, ...
        startStops + events.pointsPerSession);
    trialLen = length(startStops);
    blockStartStopInds = [1 size(labels,1)];
    events.blockStartStop = vertcat(events.blockStartStop, ...
        blockStartStopInds + events.pointsPerSession);
    events.trialsPerSession = events.trialsPerSession + trialLen;
    events.trialsPerBlock = vertcat(events.trialsPerBlock, trialLen);
    nPointsPerBlock = size(labels, 1);
    events.pointsPerSession = events.pointsPerSession + nPointsPerBlock;
    if exist('gestGoalState','var')
        events.gestGoalState = vertcat(events.gestGoalState , gestGoalState);
    end
    blocknum = [blocknum; str2num(block_folder(7:end))]; 
    if isfield(p,'oneCalBlock') && p.oneCalBlock
        % just load one block from the session
        return
    end
    fprintf('\n[%s]', block_folder)
end

cat_extra.cursorVel = cat_cursorVel;
fprintf('\nPoints Per Session: %d', length(cat_ND))

end