function fh = dpca_perMarginalization(Xfull, plotFunction, varargin)

% dpca_perMarginalization(X, plotFunction, ...) performs PCA in each
% marginalization of X and plots the components using plotFunction, a
% pointer to the function that plots one component (see dpca_plot_default() for
% the template).

% dpca_perMarginalization(..., 'PARAM1',val1, 'PARAM2',val2, ...) 
% specifies optional parameter name/value pairs:
%
%  'combinedParams' - cell array of cell arrays specifying 
%                     which marginalizations should be added up together,
%                     e.g. for the three-parameter case with parameters
%                           1: stimulus
%                           2: decision
%                           3: time
%                     one could use the following value:
%                     {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}.
%
% 'timeEvents'      - time-points that should be marked on each subplot
%  'marginalizationNames'   - names of each marginalization
%  'time'                   - time axis
%
% 'timeSplits'      - an array of K integer numbers specifying time splits
%                     for time period splitting. All marginalizations will
%                     be additionally split into K+1 marginalizations,
%                     apart from the one corresponding to the last
%                     parameter (which is assumed to be time).
%
% 'timeParameter'   - is only used together with 'timeSplits', and must be
%                     provided. Specifies the time parameter. In the
%                     example above it is equal to 3.
%
% 'notToSplit'      - is only used together with 'timeSplits'. A cell array
%                     of cell arrays specifying which marginalizations
%                     should NOT be split. If not provided, all
%                     marginalizations will be split.


% default input parameters
options = struct('combinedParams', [],       ...   
                 'timeEvents',     [],       ...
                 'time',           [], ...   
                 'marginalizationNames', [], ...
                 'timeSplits',     [],       ...
                 'timeParameter',  [],       ...
                 'notToSplit',     [], ...
                 'marginalizationColours', [], ...
                 'timeMarginalization', [], ...
                 'conditionNames', [], ...
                 'titleStr', '', ...
                 'centerData', true, ...
                 'fig',             []);

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end
options.legendSubplot = 16;
% options.conditionNames = [];
 
% centering
X = Xfull(:,:);
if options.centerData
    % centering
    X = bsxfun(@minus, X, nanmean(X,2));
end
XfullCen = reshape(X, size(Xfull));
xNanInds = any(isnan(X));
nonNanX = X(:,~xNanInds);
% total variance
totalVar = nansum(X(:).^2);

%%
% marginalize
[Xmargs, margNums] = dpca_marginalize(XfullCen, 'combinedParams', options.combinedParams, ...
                    'timeSplits', options.timeSplits, ...
                    'timeParameter', options.timeParameter, ...
                    'notToSplit', options.notToSplit, ...
                    'ifFlat', 'yes');

if ~isempty(options.timeMarginalization)
    % Shuffle time marg to be first
    margOrder = [options.timeMarginalization setdiff(1:length(Xmargs),options.timeMarginalization)];
%     Xmargs = Xmargs(margOrder);
    margNums = margNums(margOrder);
end

PCs = [];
vars = [];
margs = [];

ncompsPerMarg = 3;

for ii=1:length(Xmargs)
    %[~,S,V] = svd(Xmargs{m});      % this is very slow!
    m = margNums(ii);
    xm = Xmargs{m}(:,~xNanInds);
    margVar(m) = sum(xm(:).^2)/totalVar*100;
    
    if size(xm,1)<size(xm,2)
        XX = xm*xm';
        [U,S] = eig(XX);
        [~,ind] = sort(abs(diag(S)), 'descend');
        S = sqrt(S(ind,ind));
        U = U(:,ind);
        SV = U'*Xmargs{m};
    else
        XX = xm'*xm;
        [U,S] = eig(XX);
        [~,ind] = sort(abs(diag(S)), 'descend');
        S = sqrt(S(ind,ind));
        U = U(:,ind);
        SV = (U*S)';
    end        
    
    %PCs = [PCs; S(1:10,1:10)*V(:,1:10)'];
    PCs = [PCs; SV(1:ncompsPerMarg,:)];
    vars = [vars; diag(S(1:ncompsPerMarg,1:ncompsPerMarg)).^2];
    margs = [margs repmat(m, [1 ncompsPerMarg])];
end
[~,componentRank] = sort(vars,'descend');
% PCs = PCs(ind,:);
% margs = margs(ind);
%PCs = PCs(1:15,:);
vars = vars / totalVar * 100;

dims = size(Xfull);
Z = reshape(PCs, [length(vars) dims(2:end)]);

yspan = max(abs(Z(:)));

if ~isempty(plotFunction)
    plotFontSize = 16;
    
    nMargs = length(Xmargs);
    if isempty(options.marginalizationNames)
        options.marginalizationNames = sprintfc('Marg %d', margNums);
    end
    if isempty(options.fig)
        fh = figure;
    else
        fh = figure(options.fig);
        clf
    end
    figName = 'pca on margs';
    if ~contains(fh.Name,figName)
        fh.Name = strtrim( [fh.Name ' ' figName] );
    end
    set(fh,'units','normalized','outerposition',[0 0 1 1])

    for i=1:prod(nMargs*ncompsPerMarg)

        subplot(nMargs,ncompsPerMarg,i)

        cln = {i};
        for j=2:ndims(Z)
            cln{j} = ':';
        end
        
        plotFunction(Z(cln{:}), options.time, [-yspan yspan]*1.1, vars(i), find(componentRank==i), options.timeEvents, [], 1, options)
        

        showXLabel = ~mod(i-1,ncompsPerMarg);
        if showXLabel
            xlbl_h = ylabel( options.marginalizationNames{margs(i)} );
            xlbl_h.FontSize = plotFontSize;
            xlbl_h.FontWeight = 'bold';
        end
    
    end


% legend
if ~isempty(options.legendSubplot)   
    s = subplot(4,4,options.legendSubplot);
    delete(s)
    subplot(4,4,options.legendSubplot)
    plotFunction('legend', size(Xfull), options)
end


end
loc = [0.005 0.7 0.5 0.3];
if ~isempty(options.titleStr)
    additionalStr = sprintf('\n%s', options.titleStr);
else
    additionalStr = '';
end
str = sprintf('PCA on marginalizations%s',additionalStr);
h = annotation('textbox', loc, 'string', str, 'edgecolor', 'none', 'fontSize', plotFontSize, 'FontWeight', 'bold' );
%%
if ~nargout
    clear fh
end
end