function varargout=shadedErrorBar(x,y,errBar,varargin)
% generate continuous error bar area around a line plot
%
% function H=shadedErrorBar(x,y,errBar, ...
%
% H=shadedErrorBar(x,y,errBar, ..., 'lineProps', {...}, ...
% 
% 
% Purpose 
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
%
% Inputs (required)
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
%          (with respect to y -- see demo). 
%          
%          ** alternatively ** 
%          errBar can be a char to estimate a confidence interval (95%)
%            'boot' - use bootci to estimate the confidence interval
%            'norm' - use normfit to estimate the confidence interval
% 
%          errBar can be a cellArray of two function handles. The 
%          first defines statistic the line should be and the second 
%          defines the error bar.
%
% Inputs (optional, param/value pairs)
% 'lineProps' - ['-k' by default] defines the properties of
%             the data line. e.g.:    
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% 'transparent' - [true  by default] if true, the shaded error
%               bar is made transparent. However, for a transparent
%               vector image you will need to save as PDF, not EPS,
%               and set the figure renderer to "painters". An EPS 
%               will only be transparent if you set the renderer 
%               to OpenGL, however this makes a raster image.
% 'patchSaturation'- [0.2 by default] The saturation of the patch color.
%
% 'axes'        - specify the axes handle to plot to
%
%
% Outputs
% H - a structure of handles to the generated plot objects.
%
%
% Examples:
% y=randn(30,80); 
% x=1:size(y,2);
%
% 1)
% shadedErrorBar(x,mean(y,1),std(y),'lineprops','g');
%
% 2)
% shadedErrorBar(x,y,{@median,@std},'lineprops',{'r-o','markerfacecolor','r'});
%
% 3)
% shadedErrorBar([],y,{@median,@(x) std(x)*1.96},'lineprops',{'r-o','markerfacecolor','k'});
%
% 4)
% Overlay two transparent lines:
% clf
% y=randn(30,80)*10; 
% x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'lineprops','-r','transparent',1);
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'lineprops','-b','transparent',1);
% hold off
%
%
% Rob Campbell - November 2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input arguments
if nargin < 3, errBar='norm'; end

params = inputParser;
params.CaseSensitive = false;
params.KeepUnmatched = true;
params.addParameter('lineProps', {}, @(x) ischar(x) | iscell(x));
params.addParameter('transparent', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('edgeLineWidthFactor', 0.1, @isnumeric); % if 0 then do not plot edge lines
params.addParameter('patchSaturation', 0.2, @(x) isnumeric(x) && x>=0 && x<=1);
params.addParameter('axes', [] );

lineProps = {};
try
    params.parse(varargin{:});
catch
    % On error, use defaults and assume any inputs are lineProps
    if ~isempty(varargin) && isempty(fieldnames(params.Unmatched))
        lineProps = varargin;
    end
end

if isempty(lineProps)
    lineProps =  params.Results.lineProps;
end


%Extract values from the inputParser

transparent =  params.Results.transparent;
patchSaturation = params.Results.patchSaturation;
edgeLineWidthFactor = params.Results.edgeLineWidthFactor;
ax = params.Results.axes;


% Assume unmatched fields are str-value line properties
if isempty(lineProps) && ~isempty(params.Unmatched)
    unmatchedFields = fieldnames(params.Unmatched);
    for ii = 1:length(unmatchedFields)
        jj = (ii-1)*2 + 1;
        lineProps{jj} = unmatchedFields{ii};
        lineProps{jj+1} = params.Unmatched.(unmatchedFields{ii});
    end    
end

if isempty(ax); ax = gca; end
if ~iscell(lineProps), lineProps={lineProps}; end


%Process y using function handles if needed to make the error bar dynamically
if ischar(errBar)
    
    confidenceIntervalAlpha = 0.05; % Compute the 95% confidence interval
    
    % Confidence intervals are averaged on the first dim. Transpose if we
    % are oriented incorrectly
    if ~isempty(x)
        if ~iscell(y) && length(x) == size(y,1) && length(x) ~= size(y,2)
            y = y';
        end
    end
    
    if isvector(y)
        % If just a vector, no way to compute CI
        % Should we assume the user meant the shape?
        % Currently, forcing each point to be a time point.
        y=y(:).';
        errBar=zeros(size(y,2),1);
    else
    [y, errBar] = ComputeCI(errBar, y, confidenceIntervalAlpha );
    end
             
elseif iscell(errBar) 
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
elseif isvector(y)
    y=y(:).';
end

if isempty(x)
    x=1:length(y);
elseif isvector(x)
    x=x(:).';
end


%Make upper and lower error bars if only one was specified
if length(errBar)==numel(errBar)
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar);
    f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('length(x) must equal length(errBar)')
end


%Log the hold status so we don't change
initialHoldStatus=ishold(ax);
if ~initialHoldStatus, hold(ax,'on'),  end

H = makePlot(ax,x,y,errBar,lineProps,transparent,patchSaturation, edgeLineWidthFactor);

if ~initialHoldStatus, hold(ax,'off'), end

if nargout
    varargout{1}=H;
    varargout{2} = H.mainLine;
end

end

function H = makePlot(ax,x,y,errBar,lineProps,transparent,patchSaturation, edgeLineWidthFactor)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot to get the parameters of the line
    
    H.mainLine=plot(ax,x,y,lineProps{:});


    % Work out the color of the shaded region and associated lines.
    % Here we have the option of choosing alpha or a de-saturated
    % solid colour for the patch surface.
    mainLineColor=get(H.mainLine,'MarkerFaceColor');
    if any(strcmp( mainLineColor, {'none', 'flat'}))
        mainLineColor=get(H.mainLine,'color');
    end
    edgeColor=mainLineColor+(1-mainLineColor)*0.55;

    if transparent
        faceAlpha=patchSaturation;
        patchColor=mainLineColor;
    else
        faceAlpha=1;
        patchColor=mainLineColor+(1-mainLineColor)*(1-patchSaturation);
    end


    %Calculate the error bars
    uE=y+errBar(1,:);
    lE=y-errBar(2,:);


    %Add the patch error bar



    %Make the patch
    yP=[lE,fliplr(uE)];
    xP=[x,fliplr(x)];

    %remove nans otherwise patch won't work
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];


    patchProps = {  'facecolor',patchColor, ...
                    'edgecolor','none', ...
                    'facealpha',faceAlpha};
    if(isdatetime(x))
        H.patch=fill(ax,(xP),yP,1, patchProps{:});
    else
        H.patch=patch(ax,xP,yP,1, patchProps{:});
    end


    %Make pretty edges around the patch. 
    if edgeLineWidthFactor > 0
        edgeLineWidth = H.mainLine.LineWidth * edgeLineWidthFactor;
        H.edge(1)=plot(ax,x,lE,'-','color',edgeColor, 'LineWidth', edgeLineWidth);
        H.edge(2)=plot(ax,x,uE,'-','color',edgeColor, 'LineWidth', edgeLineWidth);
    end



    uistack(H.mainLine,'top') % Bring the main line to the top


end

function [y, errBar] = ComputeCI( computeType, y, alpha )
if ~iscell(y)
    y = {y};
end
%%
switch computeType
    case 'boot'
        %%
        nBoot = 1000;
        ii = 1;
        [u,ci] = BootstrapCi(y{ii}, nBoot, alpha);
        for ii = 2:length( y )
            [u(ii),ci(:,ii)] = BootstrapCi(y{ii}, nBoot, alpha);
        end

        
    case 'norm'
        [u, ~, ci] = normfit(y{1}, alpha);
        for ii = 2:length( y )
            [u(ii), ~, ci(:,ii)] = normfit(y{ii}, alpha);
        end

        
    otherwise
        error('Unrecognized error bar command: %s', computeType)
end

errBar = [u-ci(1,:); ...
          ci(2,:)-u];
y      = u;
end

function [u,ci] = BootstrapCi(y, nBoot, alpha)
    try
        ci = bootci(nBoot, {@nanmean, y}, 'alpha', alpha);
    catch
        ci = bootci(nBoot, {@nanmean, y}, 'alpha', alpha,'type', 'cper');
    end
    u    = nanmean(y);
end