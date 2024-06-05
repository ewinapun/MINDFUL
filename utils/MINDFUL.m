classdef MINDFUL < handle
    
% Monitoring Instability in Neural Data For Useful Long-term iBCI
% a class to calculate and plot distribution distances with various metrics
% addpath(genpath(fullfile('D:\Code\analysis\UserAnalysis\Ewina\DistMetric')))
% example:
%     event.blockStartStop = slc.blockStartStops;
%     event.blkNum = blk;
%     params.winlen = 1000;
%     params.updateHz = 50;
%     params.xlabelFromDay0 = 1;
%     params.pcaDim = 5;
%     info.participant = 'T11';
%     k = MINDFUL(data, event, info, params);
%     % optional to select features
%     k.SelectFeatures(feats);
%     % optional user-defined reference data to compare to other data
%     [refData, refInd] = k.GetBlockData(1);
%     header = sprintf('block %i',refblocknum);
%     k.CalculatePCA(refData, refInd, header)
%     k.SetReference(refData, refInd, header)
%     k.CalcDistanceFromData('filename')
%     k.PlotDistance()

% % % alternate outside for loop to calc distance with user input data
% % sTo = [];
% % tic
% % for j = 1:k.n
% %     data2 = GetSegmentData(k,j);
% %     st = k.CalcDistance(data2);
% %     sTo = ConcatStruct(sTo, st);
% % end
% % k.dist = sTo;
% % toc
% 
% History:
%   2021.02.19   created by Ewina Pun
%   2024.01.21   last edit by Ewina Pun
%
%   Copyright Ewina Pun, All Rights Reserved

    properties
        data            % neural data, [nSteps x nFeat]
        labels          % directional velocities,  [nSteps x nDim]
        event           % struct containing important events
        info            % struct contatining sessions info
        params          % Defines how to calc distances
        extra           % struct containing extra info about performance
        mAE = [];       % struct containing median angle error
        n = [];         % number of segments for comparisons
        reference = []; % reference data to calc distance
        dist            % struct contatining distance metric
        pd = [];        % probability distribution info of data
        pcaInfo = [];   % PCA transform info
        pval = [];      % Bonferroni pvalue
    end
    
    methods
        function obj = MINDFUL(data, events, info, params, extra)
            
            defaultParams.winlen = 3000; % window length (in time steps)
            defaultParams.updateHz = 1000; % frequency of updating a distance (in time steps) 
            defaultParams.excludeNonTrials = 0; % true if only include trial data
            defaultParams.gsmooth = 0; % window length for gaussian smooth data (in time steps)
            defaultParams.selected_sessions = []; % select some sessions
            defaultParams.dmetric = {'KLdiv', 'Bhatt', 'Wass'};
            defaultParams.pcaDim = 5; % number of dimension for PCA dim reduction
            defaultParams.xlabelFromDay0 = 0; % default to display trial day, or 1 = day from first session
            defaultParams.setBtwDays2NaN = 1; % set between days indices to NaN
            defaultParams.NanInds = [];
            defaultInfo.usedSessionList = 1;
            
            if nargin < 1; data = []; end
            if nargin < 2; events = []; end
            if nargin < 3; info = []; end
            if nargin < 4; params = []; end
            if nargin < 5; extra = []; end
            
            obj.data   = double(data);
            obj.event  = events;
            obj.info   = MergeBintoA(defaultInfo, info);
            obj.params = MergeBintoA(defaultParams, params);
            obj.extra  = extra; % angular error
            obj.n      = floor(length(data)/obj.params.updateHz);
            obj.pcaInfo.nDim = obj.params.pcaDim; % number of dimension using PCA
            
            if ~iscell(obj.params.dmetric)
                obj.params.dmetric = cellstr(obj.params.dmetric);
            end
            if ~isempty(obj.event)
                if isfield(obj.event,'blockStartStops') || ...
                   ~isfield(obj.event,'blockStartStop')
                    obj.event.blockStartStop = obj.event.blockStartStops;
                    obj.event = rmfield(obj.event,'blockStartStops');
                end
                % check if only some sessions are selected
                if isempty(obj.params.selected_sessions)
                   obj.params.selected_sessions = (1:length(obj.info.usedSessionList));
                end
                % default blk indices
                if ~isfield(obj.info,'blkInds')
                    obj.info.blkInds = (1:size(obj.event.blockStartStop,1));
                end
                % check if only one session
                if ~isfield(obj.event,'sessionNumberPerBlock')
                    obj.event.sessionNumberPerBlock = ones(size(obj.event.blockStartStop,1),1);
                end
                if ~isempty(obj.event.sessionNumberPerBlock)
                    blkInds = ismember(obj.event.sessionNumberPerBlock, ...
                                       obj.params.selected_sessions);
                    obj.info.blkInds = find(blkInds');
                else
                    obj.info.blkInds = (1:size(obj.event.blockStartStop,1));
                end
                obj = SelectEpoch(obj);
            end
            obj = UpdatePlotLabels(obj);
        end
        
        function obj = SelectEpoch(obj)
            SegPerBlk = zeros(length(obj.info.blkInds), 1);
            ind = []; 
            dataSelected = []; 
            indSelected = [];
            selectedStart = [];
            if obj.params.excludeNonTrials
                obj.event.trialStartStop(obj.event.excludeTrials,:) = [];
                for i = 1: length(obj.info.blkInds)
                    % exclude non-trials time steps
                    bt = find(obj.event.trialStartStop(:,1)>=...
                        obj.event.blockStartStop(obj.info.blkInds(i),1),1,'first');
                    et = find(obj.event.trialStartStop(:,2)<=...
                        obj.event.blockStartStop(obj.info.blkInds(i),2),1,'last');
                    trialInds = RowColon(obj.event.trialStartStop(bt:et,:));
                    selectedStart = [selectedStart; ...
                        trialInds(1:obj.params.updateHz:end-obj.params.updateHz)'];
                    SegPerBlk(i) = floor(length(trialInds)/obj.params.updateHz);
                    ind = [ind length(dataSelected) + ...
                        (0:SegPerBlk(i)-1)*obj.params.updateHz + 1];
                    NDTemp = obj.data(trialInds,:);
                    if obj.params.gsmooth > 0
                        NDTemp = GaussSmooth2(NDTemp, obj.params.gsmooth);
                    end
                    dataSelected = [dataSelected; NDTemp];
                    indSelected = [indSelected trialInds];
                end
                obj.data = dataSelected;
            else
                for i = 1: length(obj.info.blkInds)
                    startBlkInd = obj.event.blockStartStop(obj.info.blkInds(i), 1);
                    endBlkInd   = obj.event.blockStartStop(obj.info.blkInds(i), 2);
                    SegPerBlk(i) = floor((endBlkInd - startBlkInd)/obj.params.updateHz);
                    ind = [ind startBlkInd + (0:SegPerBlk(i)-1)*obj.params.updateHz];
                end
                indSelected = 1:size(obj.data,1);
            end
            obj.n = length(ind);
            obj.params.indSelected = indSelected;
            obj.params.SegPerBlk = SegPerBlk;
            obj.params.startStop = [ind; ind + obj.params.winlen - 1]';
            obj.params.startStop(obj.params.startStop > ...
                            length(obj.data)) = length(obj.data);
%             obj.params.startStop = [ind - obj.params.winlen;ind - 1]';
%             obj.params.startStop(obj.params.startStop < 1) = 1;
            if isfield(obj.extra,'angleError')
                obj.extra.angleError = obj.extra.angleError(obj.params.indSelected);
            end
            if isfield(obj.extra,'cursorVel')
                obj.labels = obj.extra.cursorVel(obj.params.indSelected,:);
                obj.extra = rmfield(obj.extra, 'cursorVel');
            end
        end
        
%         function obj = EvenKinSampling(obj, Xhat)
%         % even kinematic sampling
%             Xhat = Xhat(obj.params.indSelected,:);
%             [indSelected, winlen] = discretizePolar(Xhat, obj.n, obj.params);
%             obj.data = obj.data(indSelected,:);
%             obj.params.winlen = winlen;
%             obj.params.indSelected = indSelected;
%             ind = (0:obj.n-1)*obj.params.updateHz + 1;
%             obj.params.startStop = [ind; ind + obj.params.winlen - 1]';
%         end

        function obj = CalculatePCA(obj, refData, refInd, header)
        % transform data by top n PCs from refData
        % refData: [nSteps x nFeat], nFeat has to equal to that of obj.data
        % refInd : optional indices if refData is from obj.data
        % header : information about data used in pca
            refData = double(refData);
            [coeff,~,~,~,explained,mu] = pca(refData);
            obj.pcaInfo.coeff = coeff; % [nFeat x nFeat]
            obj.pcaInfo.explained = explained; % variance explained
            obj.pcaInfo.mu = mu'; % mean [nFeat x 1]
            if nargin < 4; header = []; end
            if nargin < 3; refInd = []; end
            obj.pcaInfo.ind = refInd;
            obj.pcaInfo.header = header;
            fprintf('\nPCA Info:\n')
            disp(obj.pcaInfo)
%             figure;plot(cumsum(explained),'k');ylim([0 100])
%             xlabel('dimensions');ylabel('variance explained (%)');
        end

        function obj = SetReference(obj, refData, refInd, header)
        % set reference for comparison in obj
        % refData: [nSteps x nFeat] nFeat has to equal to that of obj.data
        % refInd : [1 x nSteps] optional indices if refData is from obj.data
        % header : information about data used in reference
            if nargin < 4; header = []; end
            if nargin < 3; refInd = []; end
            obj.reference.data = refData;
            obj.reference.ind = refInd;
            obj.reference.header = header;
            fprintf('\nReference Info:\n')
            disp(obj.reference)
        end
        
        function newdata = PCAtransform(obj, varargin)
        % transform obj.data with PCA coefficients 
            if length(varargin) == 1
                if isa(varargin{1},'char')
                    indata = obj.data;
                elseif isa(varargin{1},'float')
                    indata = varargin{1};
                end
            else
                indata = obj.data;
            end
            if isa(indata,'single')
                indata = double(indata);
            end
            if isfield(obj.pcaInfo,'coeff')
                newdata = (indata - obj.pcaInfo.mu' ) * ...
                        obj.pcaInfo.coeff(:,1:obj.pcaInfo.nDim);
                obj.pcaInfo.done = 1;
            else
                newdata = indata;
                obj.pcaInfo.done = 0;
            end
        end
        
        function obj = SelectFeatures(obj, feats)
            obj.data = obj.data(:,feats);
        end
        
        function obj = CalcDistanceFromData(obj, filename, indata, refindata)        
            % perform PCA transform
            newdata = PCAtransform(obj);
            
            % check if previously reference data has been specified
            if ~isempty(obj.reference.data)
                if ~isfield(obj.pcaInfo,'nDim')
                    obj.pcaInfo.nDim = obj.params.pcaDim;
                end
                if size(obj.reference.data,2) > obj.pcaInfo.nDim
                    data1 = PCAtransform(obj, obj.reference.data);
                else
                    data1 = obj.reference.data;
                end
            else
                [~, ind] = GetSegmentData(obj,1); 
                data1 = newdata(ind, :);
            end
            
            % optional input data
            if nargin > 2
            % check if a input data was specified (not in pca)
                newdata =  [newdata indata];
                % concat in reference
                if nargin > 3
                    % concat input data in reference 
                    data1 = [data1 refindata];
                else
                    data1 = [data1 indata(obj.reference.ind,:)];
                end
            end
            
            % check if a field name for this comparison was specified
            if nargin < 2; filename = 'dist'; end
            
            % estimate distribution for data1
            pd1 = ProbDistibutionEst(data1);
            % dcat = [];
            tic
            statusprint = fprintf('Running 0 of %d', obj.n);
            for i = 1:obj.n
                fprintf(repmat('\b',1,statusprint))
                statusprint = fprintf('Running %d of %d\n', i, obj.n);
                [~, ind] = GetSegmentData(obj,i); % updateHz x #feats  
                data2 = newdata(ind, :);
                pd2 = ProbDistibutionEst(data2);
                d = CalcDistance(obj, pd1, pd2);
                fn = fieldnames(d);
                for f = 1:length(fn)
                    dcat(i).(fn{f}) = d.(fn{f});
                end
            end
            fn = fieldnames(dcat);
            for i = 1:obj.n
                for f = 1:length(fn)
                    dcat2.(fn{f})(i) = dcat(i).(fn{f});
                end
            end
                toc
            obj.dist.(filename) = dcat2;
        end
        
        function obj = CalcPairwiseDistanceParallel(obj, filename, indata)
            
            % pairwise comparison n x n for each distance metrics
            newdata = PCAtransform(obj); % perform PCA transform
            
            % check if a input data was specified and concat
            if nargin == 3
                newdata = [newdata indata];
            end
            
            % check if a field name for this comparison was specified
            if nargin < 2; filename = 'dist'; end

            % estimate prob. distribution
            for i = 1:obj.n
                [~, ind] = GetSegmentData(obj,i);
                dataSeg = newdata(ind, :);
                pd{i} = ProbDistibutionEst(dataSeg);
            end
            
            tic
            obj.pd = pd;
            % pairwise comparison
            N = obj.n;
            parfor i = 1:N
%                 statusprint = fprintf('Running %d of %d\n', i, N);
%                 fprintf(repmat('\b',1,statusprint))
                for j = 1:N
                    d = CalcDistance(obj, pd{i}, pd{j});
                    fn = fieldnames(d);
                    for f = 1:length(fn)
                        dcatdim(i,j).(fn{f}) = d.(fn{f});
                    end
                    
                end
            end
            toc
            
            fn = fieldnames(dcatdim);
            for i = 1:N
                for j = 1:N
                    for f = 1:length(fn)
                        dcatdim2.(fn{f})(i,j) = dcatdim(i,j).(fn{f});
                    end
                end
            end
            obj.dist.(filename) = dcatdim2;
        end
        
        function obj = CalcPairwiseDistance(obj, filename, indata)
            % pairwise comparison n x n for each distance metrics
            
            newdata = PCAtransform(obj); % perform PCA transform
            if nargin > 2
            % check if a input data was specified and concat
                newdata = [newdata indata];
            end
            
            % check if a field name for this comparison was specified
            if nargin < 2; filename = 'dist'; end

            % estimate prob. distribution
            for i = 1:obj.n
                [~, ind] = GetSegmentData(obj,i);
                dataSeg = newdata(ind, :);
                obj.pd{i} = ProbDistibutionEst(dataSeg);
            end
            tic
            dcatdim2 = [];
            statusprint = fprintf('Running 0 of %d', obj.n);
            % pairwise comparison
            for i = 1:obj.n
                fprintf(repmat('\b',1,statusprint))
                statusprint = fprintf('Running %d of %d\n', i, obj.n);
                dcatdim1 = [];
                for j = 1:obj.n
                    d = CalcDistance(obj, obj.pd{i}, obj.pd{j});
                    dcatdim1 = ConcatStruct(dcatdim1, d);
                end
                dcatdim2 = ConcatStruct(dcatdim2, dcatdim1, 2);
            end
            toc
            obj.dist.(filename) = dcatdim2;
        end
        
        function obj = CalcPairwisePvalue(obj)
            % pairwise comparison n x n for each distance metrics
            newdata = PCAtransform(obj);
           
            tic
            pvalue = zeros(obj.n);
            statusprint = fprintf('Running 0 of %d', obj.n);
            % pairwise comparison
            for i = 1:obj.n
                fprintf(repmat('\b',1,statusprint))
                statusprint = fprintf('Running %d of %d\n', i, obj.n);
                [~, ind1] = GetSegmentData(obj,i);
                for j = 1:obj.n
                    [~, ind2] = GetSegmentData(obj,j);
                    pvalue(i,j) = Bonferroni_pvalue(newdata(ind1,:),...
                        obj.labels(ind1,:),newdata(ind2,:),obj.labels(ind2,:));
                end
            end
            toc
            obj.pval = pvalue;
        end
        
        function d = CalcDistance(obj, data1, data2)
            % data1 is a n x d matrix of n iid d-dim multivariate normal random vectors
            % data2 is a m x d matrix of m iid d-dim multivariate normal random vectors
            
            % can also take in pd struct as data1 and data2
            % pd = ProbDistibutionEst(data)
            % d = CalcDistance(obj, pd1, pd2)
            if nargin < 3
                data2 = data1;
                data1 = obj.reference.data;
            end 
            if isa(data1, 'float')
                % flip dimension assuming n , m > d
                if size(data1, 1) < size(data1, 2)
                    data1 = data1';
                end
                pd1 = ProbDistibutionEst(data1);
            elseif isstruct(data1)
                pd1 = data1;
            else
                error('data1 should be either a double or a pd struct')
            end
            if isa(data2,'float')
                % flip dimension assuming n , m > d
                if size(data2, 1) < size(data2, 2)
                    data2 = data2';
                end
                pd2 = ProbDistibutionEst(data2);
            elseif isstruct(data2)
                pd2 = data2;
            else
                error('data2 should be either a double or a pd struct')
            end
            
            % calculate distances requested
            d = [];
            try
                if ismember('Wass',obj.params.dmetric)
                    [d.Wass, d.meanM, d.WcovM] = ...
                        CalcWasserstein(pd1, pd2);
                end
                if ismember('KLdiv',obj.params.dmetric)
                    d.KLdiv = CalcKL(pd1, pd2);
                    d.Jdiv = CalcJKL(pd1, pd2);
                end
                if ismember('Bhatt',obj.params.dmetric)
                    d.Bhatt = CalcBhattacharyya(pd1, pd2);
                end
            catch err
                disp(err.message)
            end
        end
        
        function mAE = GetMedianAEperSeg(obj, AE)
        % calculate median angle error for each segments
            mAE = [];
            if nargin == 2
                obj.extra.angleError = AE(obj.params.indSelected);
            end
            if isempty(obj.extra.angleError)
                error('median cannot be calculated. AE is not provided.')
            end
            statusprint = fprintf('Running 0 of %d', obj.n);
            for i = 1:obj.n
                fprintf(repmat('\b',1,statusprint))
                statusprint = fprintf('Running %d of %d\n', i, obj.n);
                ind = RowColon(obj.params.startStop(i,:));
                mAE = [mAE; median(obj.extra.angleError(ind),'omitnan')];
            end
            obj.mAE = mAE;
        end
        
        function [r, rho, r_pval, rho_pval] = CalcDistCorr2mAE(obj, dist, nsmooth, AE)
        % calculate Distances' Correlation to AE
        % rp: Pearson's correlation
        % rhos: Spearman's correlation
            if nargin < 3
                nsmooth = 1; % no smoothing
            end
            if nargin < 4 || isempty(AE)
                x = movmean(obj.mAE, nsmooth,'omitnan');
            else
                x = movmean(AE, nsmooth,'omitnan');
            end
            y = movmean(dist, nsmooth,'omitnan');
            if size(x, 2) > size(x, 1)
                x = x';
            end
            if size(y, 2) > size(y, 1)
                y = y';
            end
            [r,r_pval] = corr(x,y,'type','Pearson','rows','complete');
            [rho, rho_pval] = corr(x,y,'type','Spearman','rows','complete');
        end
        
        function [data, ind] = GetSegmentData(obj, segInd)
        % get data and the indices given segment indices
            ind = RowColon(obj.params.startStop(segInd,:));
            ind = unique(ind)'; % get rid of repeated indices
            data = obj.data(ind,:);
        end

        function segInds = GetSessionSegmentInds(obj, day)
        % get segment indices on a given session day
            obj.params.csSegPerBlk = [0;cumsum(obj.params.SegPerBlk)];
            loc = find(ismember(obj.event.sessionNumberPerBlock,day));
            segInds = obj.params.csSegPerBlk(loc(1))+1: ...
                      obj.params.csSegPerBlk(loc(end)+1) ...
                      - round(obj.params.winlen/obj.params.updateHz) + 1;
        end 

        function [data, ind] = GetSessionData(obj, day)
        % get data and the indices given session day
            segInds = GetSessionSegmentInds(obj, day);
            [data, ind] = GetSegmentData(obj, segInds);
        end 

        function [data, ind] = GetBlockData(obj, blocks)
        % get data and the indices given a block
            obj.params.csSegPerBlk = [0;cumsum(obj.params.SegPerBlk)];
            segInds = [];
            for loc = blocks
                inds = obj.params.csSegPerBlk(loc)+1: ...
                          obj.params.csSegPerBlk(loc+1) ...
                          - round(obj.params.winlen/obj.params.updateHz) + 1;
                segInds = [segInds inds];
            end
            [data, ind] = GetSegmentData(obj, segInds);
        end 

%         function [trialsPerSeg, trialNumStartStop] = GetTrialsPerSeg(obj)
%         % get the number of trials per segments 
%         % to do: when updateHz is small
%             trialsPerSeg = []; trialNumStartStop = [];
%             for i = 1:length(obj.info.blkInds)
%                 btnow = find(obj.event.trialStartStop(:,1)>= ...
%                             obj.event.blockStartStop(obj.info.blkInds(i),1),...
%                             1, 'first');
%                 for j = 1:obj.params.SegPerBlk(i)
%                     etnow = find(obj.event.trialStartStop(:,2)>= ...
%                                 obj.event.trialStartStop(btnow,2)...
%                                 + obj.params.updateHz, ...
%                                 1, 'first');
%                     if j == obj.params.SegPerBlk(i)
%                         trialsPerSeg = [trialsPerSeg;etnow - btnow + 1];
%                     else
%                         trialsPerSeg = [trialsPerSeg;etnow - btnow];
%                     end
%                     trialNumStartStop = [trialNumStartStop;[btnow etnow-1]];
%                     btnow = etnow;
%                 end
%             end
%         end
        
        function obj = SetIndstoNaN(obj, lastn)
        % towards the end of day not enough datapoints to estimate distribution
            if isempty(obj.params.NanInds)
                NanInds = [];
                if  nargin < 2 || isempty(lastn) 
                    lastn = round(obj.params.winlen/obj.params.updateHz);
                end
                csSegPerBlk = cumsum(obj.params.SegPerBlk);
                for day = obj.params.selected_sessions
                    loc = find(obj.event.sessionNumberPerBlock(obj.info.blkInds) == day);
                    segs = csSegPerBlk(loc(end)) - lastn + 1:csSegPerBlk(loc(end));
                    NanInds = [NanInds segs];
                end
                obj.params.NanInds = NanInds;
            end
            obj.mAE(obj.params.NanInds) = NaN;
            fn = fieldnames(obj.dist);
            for ii = 1:length(fn)
                dfn = fieldnames(obj.dist.(fn{ii}));
                for jj = 1:length(dfn)
                    % set indices to NaN if provided (in between sessions)
                    if size(obj.dist.(fn{ii}).(dfn{jj})) == [obj.n, obj.n]
                        obj.dist.(fn{ii}).(dfn{jj})(:,obj.params.NanInds) = NaN;
                        obj.dist.(fn{ii}).(dfn{jj})(obj.params.NanInds,:) = NaN;
                    else
                        obj.dist.(fn{ii}).(dfn{jj})(obj.params.NanInds) = NaN;
                    end
                end
            end
        end
        
        function obj = RemoveLastn(obj, lastn)
        %% towards the end of day not enough datapoints to estimate distribution
            if isempty(obj.params.NanInds)
                NanInds = [];
                if  nargin < 2 || isempty(lastn) 
                    lastn = round(obj.params.winlen/obj.params.updateHz);
                end
                csSegPerBlk = cumsum(obj.params.SegPerBlk);
                for day = obj.params.selected_sessions
                    loc = find(obj.event.sessionNumberPerBlock(obj.info.blkInds) == day);
                    segs = csSegPerBlk(loc(end)) - lastn + 1:csSegPerBlk(loc(end));
                    NanInds = [NanInds segs];
                end
                obj.params.NanInds = NanInds;
            end
            obj.mAE(obj.params.NanInds) = NaN;
            fn = fieldnames(obj.dist);
            for ii = 1:length(fn)
                dfn = fieldnames(obj.dist.(fn{ii}));
                for jj = 1:length(dfn)
                    % set indices to NaN if provided (in between sessions)
                    obj.dist.(fn{ii}).(dfn{jj})(obj.params.NanInds) = NaN;
                end
            end
        end
    
        function lagX = getlagX(obj, X, startStop, lag)
            % pad with zeros for lag
            lagX = zeros(size(X));
            for t = 1:size(startStop,1)
                start = startStop(t,1);
                stop = startStop(t,2);
                lagX(start+lag:stop,:)= X(start:stop-lag,:);
            end
        end
    end
    
    methods
        function dispcorr = PlotCompDistvsAE(obj, nsmooth, dfn, fn, fnDispName)
        % PlotCompDistvsAE(obj, nsmooth, dfn, fn, fnDispName)
        % generate subplots that compare different distance measures against
        % angular error, display their pearson correlation after smoothing with 
        % moving averaging
        % add AE by GetMedianAEperSeg(obj, AE) first
        % inputs:
        %     nsmooth     - moving average with a sliding window of length k
        %     dfn         - distance metric function names
        %     fn          - fieldname in m to be plot; if not provided, plot all fields in m
        %     fnDispName  - fieldname display name
        
        % set which fieldnames to be plot

            if nargin < 4 || isempty(fn) 
                fn = fieldnames(obj.dist);
            end
            
            if nargin < 5; fnDispName = fn; end
            if ~iscell(fn); fn = cellstr(fn); end                        
            if ~iscell(fnDispName); fnDispName = cellstr(fnDispName); end
            
            if any(~isfield(obj.dist, fn))
                error('One or more fieldname does not exist in struct') 
            end
            
            if nargin < 3; dfn = []; end
            
            % check dfn: distance metric function names
            [dfn, dfnStr] = CheckDistanceName(dfn);
            
            if nargin < 2; nsmooth = 1; end
            
            figure('Units','normalized','Position',[0 0 1 1],'DefaultLineLineWidth',1)
            
            % ha = ttight_subplot(length(dfn)/2,2,0.1,[],0.1);
            % set line colors
            % bcolor = brewermap(length(fn)+5,'OrRd')*0.9;
%             ha = AnnotateDetails(obj, nsmooth);
            c = get(groot,'DefaultAxesColorOrder'); 
            if length(fn) < 7
                c = [c;c*0.6];
                c(3,:) = [];
                yaxiscolor = [c(1,:);c(2,:)];
            else
                c = [c(1,:); c(2,:); cmap(length(fn)+3)*0.8];
                set(gcf, 'DefaultLineLineWidth',1.3);
                yaxiscolor = [c(1,:);c(4,:)];
            end
            
            if isempty(obj.mAE) && ~isempty(obj.extra.angleError)
                obj.GetMedianAEperSeg()
            end
            
            % set indices to NaN if provided (in between sessions)
            if obj.params.setBtwDays2NaN
                obj.SetIndstoNaN()
            end
            
            % looping over different measures
            for jj = 1:length(dfn)
                if length(dfn) > 1
                    subplot(2, 3, jj)
                end
%                 axes(ha(jj)); 
                colororder(yaxiscolor)
                yyaxis left
                % plot median AE
                y1 = movmean(obj.mAE, nsmooth,'omitnan');
                y1(obj.params.NanInds) = NaN;
                h(1) = plot(y1,'Color',c(1,:),'LineWidth',2);
                ylim([15, 150]) % to match between participants
                axis tight
                ylabel('Median AE')

                % looping over different fieldname to compare
                dispcorr = [];
                for ii = 1:length(fn)
                    yyaxis left
                    ylim([10 150])
                    y2 = obj.dist.(fn{ii}).(dfn{jj})(1,1:obj.n);
                    [rp, rhos] = CalcDistCorr2mAE(obj, y2, 1);
                    dispcorr.Pearson(ii,1) = rp;
                    dispcorr.Spearman(ii,1) = rhos;
                    if length(fn) < 10
                        limAxis = double(axis);
                        str = [(fnDispName{ii}),' - {\it r} = ', num2str(rp,3),...
                                        '; \rho = ', num2str(rhos,3)];
                        text(limAxis(2)*0.01,limAxis(3)+(limAxis(4)-limAxis(3))*(1-0.06*ii),...
                            str,'fontsize',gca().FontSize*0.8,'color',c(1+ii,:))
                    else
                        legend(h(2:end),'FontSize',15,'Location','northwest')
                        legend('boxoff')
                    end
                    hold on
                    yyaxis right
                    y2 = movmean(y2, nsmooth,'omitnan');
                    fprintf('%.2f\n',prctile(y2,99.5))
                    y2 = y2/prctile(y2,99.5);
                    y2(obj.params.NanInds) = NaN;
                    h(ii+1)= plot(y2,'Color',c(ii+1,:),'linestyle','-', ...
                        'Marker','none','DisplayName',fnDispName{ii});
                    axis tight
                    if obj.params.setBtwDays2NaN
                        obj.AddTransitionLines
                    end
                end
                disp(dfnStr{jj})
                disp(struct2table(dispcorr))
                % plotting details
                xlabel(obj.info.plot.xlabelstr)
                xticks(obj.info.plot.tickInds)
                xticklabels(obj.info.plot.tickNames)
                xlim([1 obj.n])
                title(dfnStr{jj})
                ylabel(dfnStr{jj})
            end
        end
        
        function PlotDistCorr2mAE(obj, nsmooth, dfn, fn, fnDispName)
        % generate subplots that compare different distance measures against
        % angular error, display their pearson correlation after smoothing with 
        % moving averaging
        % add AE by GetMedianAEperSeg(obj, AE) first

        % inputs:
        %     nsmooth     - moving average with a sliding window of length k
        %     dfn         - distance metric function names
        %     fn          - fieldname in m to be plot; if not provided, plot all fields in m
        %     fnDispName  - fieldname display name

            % set which fieldnames to be plot
            if nargin < 4 || isempty(fn) 
                fn = fieldnames(obj.dist);
            end            
            
            if ~iscell(fn);fn = cellstr(fn);end
            if ~iscell(fnDispName); fnDispName = cellstr(fnDispName); end

            if any(~isfield(obj.dist, fn))
                error('One or more fieldname does not exist') 
            end
            
            if nargin < 5; fnDispName = fn; end
            if nargin < 3; dfn = []; end
            % check dfn: distance metric function names
            [dfn, dfnStr] = CheckDistanceName(dfn);
            if nargin < 2; nsmooth = 1; end
            
            figure
            set(gcf, 'Position',  [0, 100, 666*length(dfn), 800],'DefaultLineLineWidth',1.5)
            % set line colors
            %bcolor = brewermap(length(fn)+5,'OrRd')*0.9;
            c = get(groot,'DefaultAxesColorOrder'); 
            if length(fn) < 8
                c = [c;c*0.6];
                c(3,:) = [];
                yaxiscolor = [c(1,:);c(2,:)];
            else
                c = [c(1,:); cmap(length(fn)+3)*0.85];
                set(gcf, 'DefaultLineLineWidth',1.3);
                yaxiscolor = [c(1,:);c(4,:)];
            end
            
            % set indices to NaN if provided (in between sessions)
            if obj.params.setBtwDays2NaN
                obj.SetIndstoNaN()
            end
            
            % looping over different measures
            for jj = 1:length(dfn)
                disp(dfn{jj})
                subplot(2, 3, jj)
                % looping over different matrices to compare
                dispcorr = [];
                for ii = 1:length(fn)
                    hold on
                    y = obj.dist.(fn{ii}).(dfn{jj});                    
                    [rp, rhos] = CalcDistCorr2mAE(obj, y, nsmooth);
                    x = movmean(obj.mAE, nsmooth,'omitnan');
                    y = movmean(y, nsmooth,'omitnan');
                    x(obj.params.NanInds) = NaN;
                    y(obj.params.NanInds) = NaN;
                    h(ii)= plot(x,y, '.','Color',c(ii+1,:),'DisplayName',fnDispName{ii});
                    [B,BINT,R,RINT,STATS] = regress(y,[ones(length(x),1) x ]);
                    axis tight
                    limAxis = double(axis);
                    shadedErrorBar([limAxis(1) limAxis(2)],...
                        [1 limAxis(1); 1 limAxis(2)]*B, ...
                        [1 limAxis(1); 1 limAxis(2)]*(B-BINT(:,1)), ...
                        'Color',c(ii+1,:))
                    str = [(fnDispName{ii}),' - {\it r} = ', num2str(rp,2),...
                                    '; \rho = ', num2str(rhos,2)];
                    text(limAxis(1),limAxis(3)+(limAxis(4)-limAxis(3))*(1-0.08*ii),...
                        str,'fontsize',12,'color',c(1+ii,:))
%                     ylim([min(obj.mAE) max(obj.mAE)])
                    xlabel('median AE')
                    dispcorr = [dispcorr;rp rhos];
                end
                title(dfnStr{jj})
                ylabel(dfnStr{jj})
            end
        end
        
        function PlotCompFnvsAE(obj, nsmooth, dfn, fn, fnDispName, axesFontsize, same_ax_color)
        % generate subplots that compare different distance measures against
        % angular error, display their pearson correlation after smoothing with 
        % moving averaging
        % add AE by GetMedianAEperSeg(obj, AE) first
        % inputs:
        %     nsmooth     - moving average with a sliding window of length k
        %     dfn         - distance metric function names
        %     fn          - fieldname in m to be plot; if not provided, plot all fields in m
        %     fnDispName  - fieldname display name
  
            % set which fieldnames to be plot
            if nargin < 4 || isempty(fn) 
                fn = fieldnames(obj.dist);
            end
            if ~iscell(fn); fn = cellstr(fn); end
            if ~iscell(fnDispName); fnDispName = cellstr(fnDispName); end
            
            if any(~isfield(obj.dist, fn))
                error('One or more fieldname does not exist in struct') 
            end
            if nargin < 7; same_ax_color = true; end
            if nargin < 6; axesFontsize = 20; end
            if nargin < 5; fnDispName = fn; end
            if nargin < 3; dfn = []; end
            [dfn, dfnStr] = CheckDistanceName(dfn);
            if nargin < 2; nsmooth = 1; end
            
            figure('Units','normalized',...
                    'Position',[0 0 .65 .65], ...
                    'DefaultLineLineWidth',1, ...
                    'DefaultAxesFontsize',axesFontsize)
            
            AnnotateDetails(obj, nsmooth);
            c = get(groot,'DefaultAxesColorOrder');
            if same_ax_color
                c = [c(1,:);repmat(c(2,:),length(fn),1)];
            else
                c(6,:) = [];c(3,:) = [];
                c = [c;c*0.6]; 
            end

            if isempty(obj.mAE) && ~isempty(obj.extra.angleError)
                obj.GetMedianAEperSeg();
            end
            
            % set indices to NaN if provided (in between sessions)
            if obj.params.setBtwDays2NaN
                obj.SetIndstoNaN();
            end
            
            % looping over different measures
            ax = zeros(length(fn),1);
            for ii = 1:length(fn)
                if length(fn) == 4
                    ax(ii) = subplot(2, 2, ii);
                else
                    ax(ii) = subplot(length(fn), 1, ii);
                end
                yaxiscolor = [c(1,:);c(1+ii,:)];
                colororder(yaxiscolor)
                yyaxis left
                % plot median AE
                y1 = movmean(obj.mAE, nsmooth,'omitnan');
                y1(obj.params.NanInds) = NaN;
                plot(y1,'Color',c(1,:),'LineWidth',1);
                axis tight
                ylabel('Median AE')

                % looping over different fieldname to compare
                y2 = obj.dist.(fn{ii}).(dfn{1})(1,1:obj.n);
                [rp, rhos, r_pval, rho_pval] = CalcDistCorr2mAE(obj, y2, nsmooth);
                
                yyaxis left
                ylim([15, 150]) % to match between participants
                if length(fn) < 8
                    limAxis = double(axis);
                    str = ['$r$ = ', num2str(rp,2),';$\rho$ = ', num2str(rhos,2)];
                    text(limAxis(2)*0.01,limAxis(4)*0.95,...
                        str,'fontsize',gca().FontSize,'color',c(1+ii,:), ...
                        'Fontweight','Bold','interpreter','latex')
                end
                hold on                
                yyaxis right
                y2 = movmean(y2, nsmooth,'omitnan');
                y2(obj.params.NanInds) = NaN;
                plot(y2,'Color',c(ii+1,:),'DisplayName',fnDispName{ii});
                axis tight
                ylim([0 prctile(y2, 99.8)])
                if obj.params.setBtwDays2NaN
                    obj.AddTransitionLines
                end
                fprintf('%.3f %.3f\n',rp, rhos)
%                 fprintf('Pearson %.3f, p = %.3e\n',rp, r_pval)
%                 fprintf('Spearson %.3f, p = %.3e\n',rhos, rho_pval)
%                 disp(dfn{1})
%                 disp(dispcorr)
                % plotting details
                xlabel(obj.info.plot.xlabelstr)
                xticks(obj.info.plot.tickInds)
                xticklabels(obj.info.plot.tickNames)
                xtickangle(0)
                xlim([1 obj.n])
                title(fnDispName{ii},'fontsize',gca().FontSize*1.5,'interpreter','latex')
                ylabel(dfnStr{1})
            end
            linkaxes(ax,'x')
        end

        function PlotCompFnvsOutlier(obj, nsmooth, dfn, fn, fnDispName, axesFontsize)
        % generate subplots that compare different distance measures against
        % angular error, display their pearson correlation after smoothing with 
        % moving averaging
        % add AE by GetMedianAEperSeg(obj, AE) first
        % inputs:
        %     nsmooth     - moving average with a sliding window of length k
        %     dfn         - distance metric function names
        %     fn          - fieldname in m to be plot; if not provided, plot all fields in m
        %     fnDispName  - fieldname display name
  
            % set which fieldnames to be plot
            if nargin < 4 || isempty(fn) 
                fn = fieldnames(obj.dist);
            end
            if ~iscell(fn); fn = cellstr(fn); end
            if ~iscell(fnDispName); fnDispName = cellstr(fnDispName); end
            if any(~isfield(obj.dist, fn))
                error('One or more fieldname does not exist in struct') 
            end
            if nargin < 6; axesFontsize = 20; end
            if nargin < 5; fnDispName = fn; end
            if nargin < 3; dfn = []; end
            [dfn, dfnStr] = CheckDistanceName(dfn);
            if nargin < 2; nsmooth = 1; end
            
            figure('Units','normalized','Position',[0 0 1 1], ...
                'DefaultLineLineWidth',1, ...
                'DefaultAxesFontsize',axesFontsize)
            
            AnnotateDetails(obj, nsmooth);
            c = get(groot,'DefaultAxesColorOrder'); 
            c = [c;c*0.6];
            c(6,:) = [];
            c(3,:) = [];
            
            % mean outliers 
            mOutlier = zeros(obj.n,1);
            for i = 1:obj.n
                ind = RowColon(obj.params.startStop(i,:));
                mOutlier(i) = mean(obj.extra.avgOutliers(ind),'omitnan');
            end        

            % set indices to NaN if provided (in between sessions)
            if obj.params.setBtwDays2NaN
                obj.SetIndstoNaN()
            end
            
            % looping over different measures
            ax = zeros(length(fn),1);
            for ii = 1:length(fn)
                if length(fn) == 4
                    ax(ii) = subplot(2, 2, ii);
                else
                    ax(ii) = subplot(length(fn), 1, ii);
                end
                yaxiscolor = [c(1,:);c(1+ii,:)];
                colororder(yaxiscolor)
                yyaxis left
                % plot mean outlier
                y1 = movmean(mOutlier, nsmooth,'omitnan');
                y1(obj.params.NanInds) = NaN;
                plot(y1,'Color',c(1,:),'LineWidth',1);
                axis tight
                ylabel('Mean Outlier (%)')
                % looping over different fieldname to compare
                y2 = obj.dist.(fn{ii}).(dfn{1})(1,1:obj.n);
                ylim([-.3, 100])
                hold on
                
                yyaxis right
                y2 = movmean(y2, nsmooth,'omitnan');
                y2(obj.params.NanInds) = NaN;
                plot(y2,'Color',c(ii+1,:),...
                            'DisplayName',fnDispName{ii});
                ylim([0 prctile(y2, 100)])
                if obj.params.setBtwDays2NaN
                    obj.AddTransitionLines
                end

                % plotting details
                xlabel(obj.info.plot.xlabelstr)
                xticks(obj.info.plot.tickInds)
                xticklabels(obj.info.plot.tickNames)
                xlim([1 obj.n])
                title(fnDispName{ii},'fontsize',gca().FontSize*1.5,'interpreter','latex')
                ylabel(dfnStr{1})
            end
            linkaxes(ax,'x')
        end        
        
        function PlotDistance(obj, dfn, fn, fnDispName, nsmooth)
        % generate subplots that compare different distance measures
        % inputs:
        %     dfn         - distance metric function names
        %     fn          - fieldname in m to be plot; if not provided, plot all fields in m
        %     fnDispName  - fieldname display name
        %     nsmooth     - moving average with a sliding window of length k
            
            % set which fieldnames to be plot
            if nargin < 5; nsmooth = 1; end
            if nargin < 3 || isempty(fn) 
                fn = fieldnames(obj.dist);
            end
            if ~iscell(fn);fn = cellstr(fn);end
            if any(~isfield(obj.dist, fn))
                error('One or more fieldname does not exist in struct m') 
            end
            if nargin < 4 || isempty(fnDispName); fnDispName = fn; end
            if ~iscell(fnDispName); fnDispName = cellstr(fnDispName); end
            if nargin < 2; dfn = []; end
            [dfn, dfnStr] = CheckDistanceName(dfn);
            
            % if obj.dist does not have a substruct with a filename
            % add a dummy substruct fn
            if any(ismember(fn,dfn)) || isempty(fn) 
                fn = {'dist'};
                obj.dist.(fn{1}) = obj.dist;
            end
            
            % set indices to NaN if provided (in between sessions)
            if obj.params.setBtwDays2NaN
                obj.SetIndstoNaN()
            end
            
            figure('Units','normalized','Position',[0 0 1 1], ...
                'DefaultLineLineWidth',1.5, ...
                'defaultAxesFontSize',20)
            AnnotateDetails(obj, nsmooth);
            c = get(groot,'DefaultAxesColorOrder');
            c = [c;c*0.6];
            
            % looping over different measures
            for jj = 1:length(dfn)
                if jj >= 1
                    ax(jj) = subplot(2, 3, jj);
                end
                colororder([c(1,:);c(2,:)])
                % looping over different matrices to compare
                y = obj.dist.(fn{1}).(dfn{jj});
                % y = obj.dist.(fn{1}).(dfn{jj})(1:obj.params.NanInds(1)-1,1:obj.params.NanInds(1)-1);
                if length(fn) == 1 && size(y,1) == size(y,2)
                % plot pairwise comparison with imagesc
                    plotPairwiseFlag = 1;
                    % set NaN for values not used for calculated meanKLD
                    mask = tril(NaN(obj.n), -1);
                    y = y+mask;
%                     y(isnan(y)) = -1;
%                     imagesc(y)
                    imagesc(y,'AlphaData',~isnan(y))
                else
                    plotPairwiseFlag = 0;
                    for ii = 1:length(fn)
                        y = obj.dist.(fn{ii}).(dfn{jj});
                        if size(y,1) == 1 || size(y,2) == 1
                            y = movmean(y, nsmooth,'omitnan');
                            y(obj.params.NanInds) = NaN;
                            plot(y,'Color',c(ii+1,:),...
                                'DisplayName',fnDispName{ii})
                        end
                        if obj.params.setBtwDays2NaN
                            obj.AddTransitionLines
                        end
                        hold on
                    end
                end
              % plotting details
                title(dfnStr{jj})
                ylabel(dfnStr{jj})
                xlabel(obj.info.plot.xlabelstr)
                xticks(obj.info.plot.tickInds)
                xtickangle(45)
                xticklabels(obj.info.plot.tickNames)
                if plotPairwiseFlag
                    hold on
                    yticks(obj.info.plot.tickInds)
                    yticklabels(obj.info.plot.tickNames)
                    ylabel(obj.info.plot.ylabelstr)
                    % lineInds = obj.info.plot.lineInds;
                    % range = (lineInds(1)+1:lineInds(end));
                    clim([0 round(prctile(y, 99 ,'all')*1.1,1,'significant')])
                    cbar = colorbar;
%                     colormap([.25 .25 .25;colormap]);
                    set(cbar,'YLim',[0 4])%cbar.YLim(2)]);
                    set(gca,'YDir','reverse', 'Layer','top','Color',[1 1 1]*.25)
                end
                axis tight
            end
            if length(fn) > 1
                legend()
            end
            linkaxes(ax,'xy')
        end
        
        function PlotPairwisePval(obj)
            % plot pairwise comparison with imagesc
            figure('Position',[0 0 570 500],'DefaultLineLineWidth',1.5)
            imagesc(obj.pval)
            cb = colorbar();
            set(gca,'YDir','reverse', 'Layer','top')
            % plotting details
            xlabel(obj.info.plot.xlabelstr)
            xticks(obj.info.plot.tickInds)
            xticklabels(obj.info.plot.tickNames)
            ylabel(obj.info.plot.ylabelstr)
            yticks(obj.info.plot.tickInds)
            yticklabels(obj.info.plot.tickNames)
            clim([0 0.05])
            ylabel(cb,'P-value')
            obj.AddTransitionLines('xy')
            axis tight
        end
        
        function obj = UpdatePlotLabels(obj, skip)
            if length(obj.info.usedSessionList) == 1
                % label multiple blocks instead of days
                obj.info.plot.xlabelstr = 'block';
                if isfield(obj.event, 'blkNum')
                    xtickstr = obj.event.blkNum;
                else
                    xtickstr = obj.info.blkInds;
                end
                ticklabels = xtickstr;
                xlabelstr = 'block';
            else
                % label multiple sessions days
                if isfield(obj.info, 'blkInds')
                xtickstr = obj.info.trialDay(...
                    obj.event.sessionNumberPerBlock(obj.info.blkInds));
                xlabelstr = 'Trial day';
                if obj.params.xlabelFromDay0
                    xtickstr = xtickstr - xtickstr(1);
                    xlabelstr = 'Days since first session';
                end
                ticklabels = xtickstr(diff([-1 xtickstr])~=0);
                else
                    xtickstr = [];
                    xlabelstr = [];
                    ticklabels = [];
                end
            end
            obj.info.plot.xtickstr = xtickstr;
            obj.info.plot.xlabelstr = xlabelstr;
            obj.info.plot.ylabelstr = [];
            obj.info.plot.lineInds = [0;cumsum(obj.params.SegPerBlk)];
            tickInds = obj.info.plot.lineInds(diff([-1 xtickstr])~=0)+1;
            
            % decide if display every block/day or skip some
            totalChar = numel(num2str(xtickstr'));
            if nargin < 2 || isempty(skip)
                skip = ceil(totalChar/45);
            end
            obj.info.plot.tickNames = ticklabels(1:skip:end);
            obj.info.plot.tickInds = tickInds(1:skip:end);
            
        end
        
        function ha = AnnotateDetails(obj, nsmooth)
        % add annotation text at the corner of the plot
            if nargin < 2 || isempty(nsmooth) || nsmooth == 1
                nsmooth = 0; 
            end
            % info string to display
            disstr = sprintf(['%i blocks over %i sessions\n\n',...
                           'sampling Hz: %gs\n',...
                           'window len: %gs\n',...
                           'smooth steps: %g\n\n'],...
                            length(obj.info.blkInds),...
                            length(obj.info.usedSessionList),...
                            obj.params.updateHz*0.02, ...
                            obj.params.winlen*0.02, ...
                            nsmooth);
            % print header of reference data
            if isfield(obj.reference,'header') && ~isempty(obj.reference.header)
            disstr = [disstr, sprintf('reference: \n%s\n',obj.reference.header)];
            end
            % print size of reference data
            if isfield(obj.reference,'data')
                if  ~isempty(obj.pcaInfo) && obj.pcaInfo.done
                    dataDimension = obj.pcaInfo.nDim;
                else
                    dataDimension = size(obj.reference.data,2);
                end
                disstr = [disstr, sprintf('size~(%d x %d)\n\n',...
                                size(obj.reference.data,1), dataDimension)];
            end
            % print header of pca data
            if isfield(obj.pcaInfo,'header') && ~isempty(obj.pcaInfo.header)
                disstr = [disstr, sprintf('PCA: \n%s',obj.pcaInfo.header)];
            end
            % add annotation
            ha = annotation('textbox',[0.005 0.7 0.075 0.3], 'string', disstr, ...
            'edgecolor', 'none', 'Interpreter', 'none', 'FontSize', 12);
        end

        function AddTransitionLines(obj, plotparam)
            if nargin < 2
                plotparam = 'x';
            end
            csSegPerBlk = cumsum(obj.params.SegPerBlk)+1;
            for day = obj.params.selected_sessions
                loc = find(obj.event.sessionNumberPerBlock(obj.info.blkInds) == day);
                % if ismember(csSegPerBlk(loc(end)),obj.info.plot.tickInds)
                    switch plotparam
                        case 'x'
                            xline(csSegPerBlk(loc(end)),'Color',[.5 .5 .5]);
                        case 'xy'
                            xline(csSegPerBlk(loc(end))-0.5,'Color',[.5 .5 .5]);
                            yline(csSegPerBlk(loc(end))-0.5,'Color',[.5 .5 .5]);
                    end
                % end
            end
        end


    end

end