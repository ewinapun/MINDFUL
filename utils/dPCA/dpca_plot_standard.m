function dpca_plot_standard(data, time, yspan, explVar, compNum, events, signif, marg, options)

% Modify this function to adjust how components are plotted.
%
% Parameters are as follows:
%   data      - data matrix, size(data,1)=1 because it's only one component
%   time      - time axis
%   yspan     - y-axis spab
%   explVar   - variance of this component
%   compNum   - component number
%   events    - time events to be marked on the time axis
%   signif    - marks time-point where component is significant
%   marg      - marginalization number


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displaying legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if nargin < 9
    options = [];
end
allMarkers = {'-', '--', '.-', ':', '-o', '-s', '-^', '-p', '-h', '-*'};

if strcmp(data, 'legend')
    
    % Assuming legend is the last step, let's combine all axes that plotted
    % dPC projections
    LinkPlotProjectinAxes()
    
    
    % if there is only time and no other parameter - do nothing
    if length(time) == 2
        return
    end
         
    isSpecial5 = length(time) == 5 && time(3) == 2 && time(4) == 2;
    
    % if there is one parameter
    if length(time) == 3
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        
        options = yspan;
        colors = GetColors(numOfStimuli, options);
        if isempty(options.conditionNames)
            options.conditionNames = sprintfc('Stimulus %d', 1:numOfStimuli);
        end
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, options.conditionNames{f})
        end
        axis([0 3 -1 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return

    % two parameters: stimulus and decision
    elseif length(time) == 4
        options = yspan;
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        nCond2 = time(3);
        if isempty(options.conditionNames)
            options.conditionNames = cat(2, sprintfc(' %d', 1:numOfStimuli), sprintfc(' %d', 1:nCond2) );
        end
        if isempty(options.marginalizationNames)
            options.marginalizationNames = sprintfc('Condition %d', 1:2);
        end
        colors = GetColors(numOfStimuli, options);
        hold on
        %%
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, sprintf('%s: %s',options.marginalizationNames{1}, options.conditionNames{f}));
        end
        
        if numOfStimuli == nCond2 && nCond2 == length(options.conditionNames)
            cond2Offset = 0;
        else
            cond2Offset = numOfStimuli;
        end
        for ii = 1:nCond2
            yLoc = -(ii+1);
            markerI = mod(ii,length(allMarkers)) + 1;
            plot([0.5 1], [yLoc yLoc], ['k' allMarkers{markerI}], 'LineWidth', 2)
            text(1.2, yLoc, sprintf('%s: %s', options.marginalizationNames{2}, options.conditionNames{cond2Offset+ii}) )
        end
                
        
        axis([0 3 yLoc-0.75 0.75+numOfStimuli])
        
        
        nTotalConditions = numOfStimuli + nCond2;
        UpdateLegendAxesSize(nTotalConditions)
        
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return
    elseif isSpecial5
        %%
        % different stimuli in different colours and binary condition as
        % solid/dashed
        options = yspan;
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        if isempty(options.conditionNames)
            options.conditionNames = sprintfc('Stimulus %d', 1:numOfStimuli);
        end
        colors = GetColors(numOfStimuli, options);
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, [options.marginalizationNames{1} '  ' num2str(f)])
        end
        plot([0.5 1], [-2 -2], 'k-', 'LineWidth', 2)
        plot([0.5 1], [-3 -3], 'k-+', 'LineWidth', 2)
        
        plot([0.5 1], [-4 -4], 'k--', 'LineWidth', 2)
        plot([0.5 1], [-5 -5], 'k:', 'LineWidth', 2)

        if isfield(options, 'conditionNames') && length(options.conditionNames) == 4
            text(1.2, -2, options.conditionNames{1})
            text(1.2, -3, options.conditionNames{2})

            text(1.2, -4, options.conditionNames{3})
            text(1.2, -5, options.conditionNames{4})
            
        else
            text(1.2, -2, sprintf('%s 1 %s 1', options.marginalizationNames{2:3}))
            text(1.2, -3, sprintf('%s 2 %s 1', options.marginalizationNames{2:3}))

            text(1.2, -4, sprintf('%s 1 %s 2', options.marginalizationNames{2:3}))
            text(1.2, -5, sprintf('%s 2 %s 2', options.marginalizationNames{2:3}))
        end
        
        axis([0 3 -5.5 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        ax = gca;
        ax.Position(4)  = ax.Position(4) + ax.Position(4)*.25;
        return;
    % other cases - do nothing
    else
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up the subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(time)
    time = 1:size(data, ndims(data));
end
axis([time(1) time(end) yspan])
hold on

if ~isempty(explVar)
    title(['Component #' num2str(compNum) ' [' num2str(explVar,'%.1f') '%]'])
else
    title(['Component #' num2str(compNum)])
end

if ~isempty(events)
    plot([events; events], yspan, 'Color', [0.6 0.6 0.6])
end

if ~isempty(signif)
    signif(signif==0) = nan;
    plot(time, signif + yspan(1) + (yspan(2)-yspan(1))*0.05, 'k', 'LineWidth', 3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if ndims(data) == 2
    % only time - plot it
    plot(time, squeeze(data(1, :)), 'k', 'LineWidth', 2)

elseif ndims(data) == 3
    % different stimuli in different colours
    numOfStimuli = size(data, 2);
    colors = GetColors(numOfStimuli, options);
    for f=1:numOfStimuli
        plot(time, squeeze(data(1,f,:)),'color', colors(f,:), 'LineWidth', 2)    
    end

elseif ndims(data) == 4
    % different stimuli in different colours and different markers for each condition
    numOfStimuli = size(data, 2);
    numCond2 = size(data, 3);
    colors = GetColors(numOfStimuli, options);

    for f=1:numOfStimuli
        for ci = 1:numCond2
            markerI = mod(ci,length(allMarkers)) + 1;
            m = allMarkers{markerI};
            plot(time, squeeze(data(1, f, ci, :)), m, 'color', colors(f,:), 'LineWidth', 2, 'MarkerSize', 3)        
        end
    end

elseif ndims(data) == 5 && size(data,3)==2 && size(data,4)==2
    % different stimuli in different colours and binary condition as
    % solid/dashed
    numOfStimuli = size(data, 2);
    colors = GetColors(numOfStimuli, options);

    for f=1:numOfStimuli 
        plot(time, squeeze(data(1, f, 1, 1, :)), 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 2, 1, :)), '-.', 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 1, 2, :)), '--', 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 2, 2, :)), ':', 'color', colors(f,:), 'LineWidth', 2)
    end

else
    % in all other cases pool all conditions and plot them in different
    % colours
    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    data = reshape(data, size(data,1), []);
    data = data';
    
    plot(time, data, 'LineWidth', 2)    
end

ax = gca;
ax.Tag = 'Projection plot';


function colors = GetColors(nColors, options)
    try
        colors = cmap(nColors);
    catch
        colors = lines(nColors);
    end


function UpdateLegendAxesSize(nTotalConditions)
    axYPerCond = 0.02;
    
    ax = gca;
    if ax.Position(4)/nTotalConditions < axYPerCond
        newAxY = axYPerCond*nTotalConditions;
        shiftAxY = newAxY - ax.Position(4);
%         ax.Position(2) = max(ax.Position(2) - shiftAxY, 0); % Do not shift below 0
%         ax.Position(4) = newAxY;
    end

    
function LinkPlotProjectinAxes()
    % We set the axes to have a Tag at bottom of plotting projections
    % Line ~232
    f = gcf;
    axs = findobj( f.Children, 'Tag', 'Projection plot', '-depth', 1);
    linkaxes(axs,'xy');