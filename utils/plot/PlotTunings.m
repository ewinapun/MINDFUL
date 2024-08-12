classdef PlotTunings
% This class is a collection of plot functions for visualizing cosine tuning curves
% it is optimized to compare across days tuning changes
%
% tunings
%       pd: [nfeat × nday double] preferred directions
%       md: [nfeat × nday double] modulation depth
%        b: [nfeat × 3 x nday double] estimated parameters to 
%                                     [1 cos(theta) sin(theta)]
%      rsq: [nfeat × nday double] R-square to the regression fit 
%     pval: [nfeat × nday double] p-value of the regression fit
% Example:
% pt = PlotTunings(tunings)
% pt.GetPairwiseCorrelation
% GetSortedOrder          
% PlotAllTuning           
% PlotBasicHeatmap        
% PlotDeltaTuning         
% PlotExampleTuning       
% PlotTunings             
% PlotTuningsCompare      

% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

    properties
        tunings (1,1) struct  % tunings struct
        nfeat (1,1) double    % number of features in tunings
        nday (1,1) double     % number of of days
        selected_feats double % which features to be displayed; can be ordered
    end

    methods
        function obj = PlotTunings(tunings)
            obj.tunings = tunings;
            obj.nfeat = size(tunings.pd, 1);
            obj.nday = size(tunings.pd, 2);
            obj.selected_feats = (1:obj.nfeat);
        end

        function [hFig, clrs] = PlotExampleTuning(obj, ch, options)
        % plot tuning curves across sessions
            arguments
                obj PlotTunings
                ch double % channel or feature number
                options.tunings = [] % specify tuning
                options.style = 'cosine' % plot the cosine tuning (default)
                % if other options, cosine tunings become dotted line
                % other style options include
                % 'binavg': bin average firing rate 
                % 'kernreg' kernreg
                options.colorbar = [] % TickLabels for colorbar
                options.fig = []
                options.day = 1:obj.nday % specify which day to plot
                options.ylim = []
                options.title = [];
            end     
            assert(ch > 0, 'Channel number is invalid. Has to be an positive integer')
            
            if ~isempty(options.tunings)
                tuning = options.tunings;
            elseif ~isempty(obj.tunings)
                tuning = obj.tunings;
            end
    
            if isempty(options.fig)
                hFig = figure('Units','inches','Position',[1 1 4 3]);
            else
                hFig = options.fig;
            end
            
            rmv = 1;
            c = brewermap(obj.nday + rmv * 2,'RdYlBu');
            mid = ceil(length(c)/2);clrs = c;
            clrs(mid-rmv:mid+rmv,:)=[]; clrs = flipud(clrs);
            clrs = vertcat(c(end,:)*0.5,clrs);
            
            x = (-180:180); % discretize angles
            y = squeeze(tuning.b(ch,1,options.day)) + ...
                squeeze(tuning.b(ch,2,options.day)) * cosd(x) + ...
                squeeze(tuning.b(ch,3,options.day)) * sind(x);
            hLine = plot(x,y,'Tag','Signal');
            set(hLine,{'Color'},num2cell(clrs(options.day,:),2),'LineWidth',1.5);
            
            bin_ymax = 0;
            if ~strcmp(options.style,'cosine') % add binned feature average
                set(hLine,'LineStyle','--');
                if strcmp(options.style,'binavg')
                    binstyle = tuning.binavg;
                elseif strcmp(options.style,'kernreg')
                    binstyle = tuning.binkr;
                end
                biny = binstyle.means;
                for day = options.day
                    shadedErrorBar(binstyle.xi, ...
                    biny(:,ch,day), binstyle.CIs(:,ch,day) ...
                    ,'lineprops',{ ...
                    'Color', clrs(day,:), ...
                    'markerfacecolor',clrs(day,:), ...
                    'LineWidth',1.5, ...
                    },'patchSaturation',0.1);
                end   
                bin_ymax = max(max(max(biny(:),[],2)), -min(min(biny(:),[],2))) * 1.1;
            end
            
            % set ylim maximum
            if ~isempty(options.ylim)
                ymax = options.ylim;
            else
                est_ymax = max(max(y,[],'all'), - min(y,[],'all')) * 1.1;
                ymax = max([est_ymax; bin_ymax]);
            end
            ylabel('firing rate (Hz)')
            axis([-200 200 -ymax ymax])
            axis square
            xlabel('Direction')
            xticks([-180 0 180])
            xtickangle(0)
            xticklabels({'-180\circ','0\circ', '180\circ'})
            
            if isempty(options.title)
                title(['Feat ', num2str(ch)])
            else
                title(options.title)
            end

            % set ylabel
            if obj.tunings.iszscore
                ylabel('z-scored FR (Hz)')
            else
                ylabel('FR (Hz)')
            end
    
            % set colorbar
            if ~isempty(options.colorbar)
                axh = gca();
                colormap(axh,clrs); 
                cbh = colorbar();
                ylabel(cbh,'Days')
                clim = cbh.Limits;
                ty = linspace(clim(1),clim(2),obj.nday+1);
                if length(options.day)==obj.nday
                    if obj.nday > 6
                        ind = (1:2:obj.nday);
                    else
                        ind = (1:obj.nday);
                    end
                else
                    ind = options.day;
                end
                cbh.Ticks = ty(ind)+ (ty(2)-ty(1))/2;
                cbh.TickLabels = options.colorbar(ind);
            end
            
            % add marker to significantly tuned units when plotting cosine
            % tunings
            if strcmp(options.style,'cosine')
                plotm = (tuning.pval(ch,options.day) < 0.05);
                for d = 1:length(options.day)
                    if plotm(d)
                        if length(options.day)==1
                            hl = hLine;
                        else
                            hl = hLine(d);
                        end
                        [~,iPk] = max(y(d,:));
                        hAxes = ancestor(hl, 'axes');
                        % use the color of the line
                        % color = get(hl,'Color');
                        hl = line(hl.XData(iPk),y(d,iPk), ...
                            'Parent',hAxes, ...
                            'Marker','o', ...
                            'MarkerEdgeColor','k', ...
                            'LineStyle','none', ...
                            'Color',get(hl,'Color'), ...
                            'tag','Peak');
                        % using MATLAB use offset inverted triangular marker
                        signal.internal.findpeaks.plotpkmarkers(hl,y(d,iPk));
                    else
                        % mark curve to dotted if not significant
                        hl.LineStyle = '--';
                    end
                end
            end
        end
        
        function PlotTuningsCompare(obj, ch, tunings2)
            % used for comparing two sets of tunings 
            % (e.g. z-score and no-zscore)
            newfig = figure('Units','inches','Position',[0 0 10 5]);
            tcl = tiledlayout(newfig,1,2);

            PlotExampleTuning(obj, ch);
            curfig = get(groot,'CurrentFigure');figure(curfig);
            ax=gca; ax.Parent=tcl; ax.Layout.Tile=1;close(curfig);
            
            PlotExampleTuning(obj, ch, tunings2);
            curfig = get(groot,'CurrentFigure');figure(curfig);
            ax=gca; ax.Parent=tcl; ax.Layout.Tile=2;close(curfig);
            
%             figlist = get(groot,'Children');
%             newfig = figure('Units','inches','Position',[0 0 7 3]);
%             tcl = tiledlayout(newfig,1,2);
%             for i = 1:numel(figlist)
%                 figure(figlist(i));
%                 ax=gca; ax.Parent=tcl;ax.Layout.Tile=i;
%             end
%             close(figlist);
        end

        function hFig = PlotExampleTuningInteractive(obj, ch, varargin)
        % plot tuning curves across sessions
        % can use left/right arrow key to change between ch
            hFig = PlotExampleTuning(obj, ch, varargin{:});
%             set(hFig,'Units','inches','Position',[0 0 4.5 4]);
            set(hFig,'KeyPressFcn', @OnKeyPress);
            set(hFig,'UserData',0);
            drawnow();
            while isvalid(hFig)
                pause(0.005);
                if ~isvalid(hFig)
                    break
                end
                next = get(hFig,'UserData');
                if next~=0
                    clf(hFig)
                    ch = ch + next;
                    PlotExampleTuning(obj, ch, 'fig', hFig, varargin{:});
                    set(hFig,'UserData',0);    
                end
            end
            
            function OnKeyPress(obj, event)
                switch event.Key
                    case 'leftarrow'
                        next = - 1;
                    case 'rightarrow'
                        next = 1;
                end
                set(obj,'UserData', next)
            end
        end

        function PlotBasicHeatmap(obj, data, varargin)
            figure('Units','inches','Position',[0 0 6.5 6])
            himg = imagesc(data);
            ConfigPlotDetails(obj, himg, varargin{:})
        end

        function PlotSignificantHeatmap(obj, sigdMD, sigdPD, varargin)
            % only works for matrix that only contains (NaN, 0 and 1)
            figure('Units','inches','Position',[0 0 12.5 6])
            % gray out NaNs
            ax(1) = subplot(1,2,1);
            sigdMD(isnan(sigdMD)) = -1;
            himg = imagesc(sigdMD);
            ConfigPlotDetails(obj, himg, varargin{:})
            title('Significant MD change')
            colormap(ax(1), [.75 .75 .75; parula(2)])
            ax(2) = subplot(1,2,2);
            sigdPD(isnan(sigdPD)) = -1;
            himg = imagesc(sigdPD);
            ConfigPlotDetails(obj, himg, varargin{:})
            title('Significant PD change')
            colormap(ax(2), [.75 .75 .75; parula(2)])
        end

        function PlotAllTuning(obj, options)
            % plot prefered direction and modulation depth

            arguments
                obj PlotTunings
                options.fig = [];
                options.xtickstr = (1:obj.nday)
                options.gray_out_non_sig = false % based on p-value
                options.title = ''
            end
            
            pd = obj.tunings.pd;
            md = obj.tunings.md;
            
            if options.gray_out_non_sig
                pd(obj.tunings.pval > 0.05) = 1000;
                md(obj.tunings.pval > 0.05) = 1000;
            end
            
            if isempty(options.fig)
                hFig = figure('Units','inches','Position',[0 0 12 6]);
            else
                hFig = options.fig;
            end

            ax(1) = subplot(1,2,1);
            himg = imagesc(md(obj.selected_feats,:));
            obj.ConfigPlotDetails(himg, 'cbtitle','MD', ...
                'xtickstr',options.xtickstr, ...
                'clim',[0 1.05], 'cbXTick',[0, 0.5, 1]);
            k = brewermap(200,'OrRd');
            colormap(ax(1),[k;.75 .75 .75]);
            
            ax(2) = subplot(1,2,2);
            himg = imagesc(pd(obj.selected_feats,:));
            kb = brewermap(180,'Blues');
            kg = brewermap(180,'Greens');
            colormap(ax(2),[kb; flipud(kg);.75 .75 .75]);
            obj.ConfigPlotDetails(himg, 'cbtitle','PD', ...
                'xtickstr',options.xtickstr, ...
                'clim', [-180 180], ...
                'cbXTickLabel', {'-180\circ','0\circ','180\circ'}, ...
                'cbXTick', [-180, 0, 180]);
            linkaxes([ax(1),ax(2)])
            sgtitle(options.title,'FontWeight','Bold')
        end

        function hFig = PlotDeltaTuning(obj, deltaMD, deltaPD, options)
            % plot prefered direction and modulation depth changes
            arguments
                obj PlotTunings
                deltaMD double
                deltaPD double
                options.fig = []
                options.xtickstr = (1:obj.nday) % xtick label strings
                options.xticks = (1:obj.nday) % x tickx 
                options.gray_out_non_sig = false % based on p-value
                options.skip_diff_col = true % skip plotting the diff column
                options.diff_col = 1 % default to first column to calculate difference
                options.title = '' % plot title
            end
            
            clim_l = -0.5;%round(prctile(deltaMD(:),1),1,'significant');
            clim_u = -clim_l;
            crange = [clim_l clim_u];
            deltaMD(deltaMD <= clim_l) = clim_l+1;
            deltaMD(deltaMD >= clim_u) = clim_u-1;

            if options.gray_out_non_sig
                deltaMD(obj.tunings.pval > 0.05) = 1000;
                deltaPD(obj.tunings.pval > 0.05) = 1000;
            end

            if options.skip_diff_col
                options.xticks(end) = [];
                options.xtickstr(:,options.diff_col) = [];
                deltaPD(:,options.diff_col) = [];
                deltaMD(:,options.diff_col) = [];
            end
            if isempty(options.fig)
                hFig = figure('Units','inches','Position',[1 2 12 6]);
            else
                hFig = options.fig;
            end
            ax(1) = subplot(1,2,1);
            himg = imagesc(deltaMD(obj.selected_feats,:));
            k = brewermap(200,'RdBu');
            colormap(ax(1), [k;.75 .75 .75])
            dispclim = linspace(clim_l,clim_u,3); 
            displ = dispclim;
%             displ = cell(3,1);
%             for i=1:length(dispclim);displ{i} = [num2str(dispclim(i)),'%'];end
            obj.ConfigPlotDetails(himg, 'xticks', options.xticks, ...
                                    'xtickstr', options.xtickstr, ...
                                    'clim', crange, ...
                                    'cbtitle', '\DeltaMD from day 0', ...
                                    'cbXTick', dispclim, cbXTickLabel=displ);
            
            ax(2) = subplot(1,2,2);
            himg = imagesc(deltaPD(obj.selected_feats,:));
            k = brewermap(179,'Blues');
            colormap(ax(2),[k;.75 .75 .75]);
            obj.ConfigPlotDetails(himg, 'xticks', options.xticks, ...
                                    'xtickstr', options.xtickstr, ...
                                    'clim', [0 181], ...
                                    'cbtitle', '|\DeltaPD| from day 0', ...
                                    'cbXTickLabel', {'0\circ','180\circ'}, ...
                                    'cbXTick', [0,180]);
            linkaxes([ax(1),ax(2)])
        end
        
        function [corrR, pVal, hFig] = GetPairwiseCorrelation(obj, options)
            % pairwise correlation plot
            arguments
                obj PlotTunings
                options.fig = []
                options.xtickstr double = (1:obj.nday)
                options.interp double = [] % to interpolate on both axis
                options.clim (1,2) double = [0 1]
                options.cbXTick = []
            end
            assert(length(options.xtickstr)==obj.nday, ...
                'length of xtickstr does not match with nday')
            % compute pairwise correlation and p-values using tuning weights 
            obj.nfeat = size(obj.tunings.pval,1);
            obj.nday = size(obj.tunings.pval,2);
            tuneMap = obj.tunings.b; % obj.nfeat x 3 x obj.nday
            corrR = zeros(obj.nday);
            pVal = zeros(obj.nday);
            disp(['selected_feats:', num2str(length(obj.selected_feats))])
            for i = 1:obj.nday
                tuneMap(obj.tunings.pval(:,i)>0.05,:,i) = NaN;
                for j = 1:obj.nday
                    tuneMap(obj.tunings.pval(:,j)>0.05,:,j) = NaN;
                    % 'Rows','complete' ignore rows with NaN
                    [R,P] = corrcoef(tuneMap(obj.selected_feats,:,i), ...
                                     tuneMap(obj.selected_feats,:,j), ...
                                     'Rows','complete');
                    corrR(i,j) = R(1,2); % correlation
                    pVal(i,j) = P(1,2); % p-values
                end
            end
            
            if isempty(options.fig)
                hFig = figure('Units','inches','Position',[0 0 6.5 6]);
            else
                hFig = options.fig;
            end
            % interpolation of axis to account for days apart between
            % session
            if ~isempty(options.interp)
                x = repmat(options.xtickstr, length(options.xtickstr), 1);
                xi = linspace(0, max(options.xtickstr), options.interp);
                x2 = repmat(xi, length(xi), 1);
                V = interp2(x, x', corrR, x2, x2', 'linear');
                imagesc(V)
                loc = find(diff([0 discretize(round(xi),options.xtickstr)])~=0);
                xtickloc = [1 loc(2:end)-1 options.interp];
            else
                imagesc(corrR)
                xtickloc = (1:obj.nday);
            end
            xtickstr = options.xtickstr;
            if length(xtickloc) > 7
                xtickloc = xtickloc(1:2:end);
                xtickstr = options.xtickstr(1:2:end);
            end
            xlabel('Days');xticks(xtickloc);xticklabels(xtickstr)
            ylabel('Days');yticks(xtickloc);yticklabels(xtickstr)
            xtickangle(45)
            hcb = colorbar;
            clim(options.clim)
            if ~isempty(options.cbXTick)
                hcb.XTick = options.cbXTick;
            end
            Rpos = hcb.Position;
            hcb.Title.String = 'R';
            hcb.Title.Units = 'normalized';
            hcb.Title.Position = [.5,0.42,0];
            axis tight
            axis square
            set(gca,'YDir','reverse', 'TickDir', 'out')
            rectangle('position', [.5, .5, xtickloc(end),xtickloc(end)], 'edgecolor', [0 0 0])
        end

        function hcb = ConfigPlotDetails(obj, himg, options)
            % add figure configs
            arguments
                obj (1,1) PlotTunings
                himg 
                options.xticks double = (1:obj.nday)
                options.xtickstr double = (1:obj.nday)
                options.title (1,1) string = ''
                options.cbtitle (1,1) string = ''
                options.clim  double {mustBeNumeric} = []
                options.cbXTick double = []
                options.cbXTickLabel string = []
            end
            himg_nfeat = size(himg.CData,1);
            rectangle('position', [.5, .5, size(himg.CData,2) size(himg.CData,1)], 'edgecolor', [0 0 0])
            axis tight
            xticks(options.xticks)
            xticklabels(options.xtickstr)
            xtickangle(45)
            xlabel('Days since first session')

            if mod(himg_nfeat, 96) == 0
                yticks(linspace(0, obj.nfeat, ceil(obj.nfeat/96)+1))
                if himg_nfeat ==384
                    ylabel('Features (ncTX, SpkW)')
                    text(-.25, 336,'M')
                    text(-.25, 240,'L')
                else
                    ylabel('Features (ncTX)')
                end
                text(-.25, 144,'M')
                text(-.25, 48,'L')
            else
                ylabel('Features')
            end

            % indicate two feature types
            if himg_nfeat == 384;line([.5 obj.nday+0.5], [192,192], 'Color', 'k','lineWidth',2);end

            title(options.title)
            
            % colorbar
            hcb = colorbar;
            hcb.Title.String = options.cbtitle;
            hcb.Title.Rotation = 90;
            hcb.Title.Units = 'normalized';
            hcb.Title.Position = [4.25,0.5,0];
            if ~isempty(options.clim)
                clim(options.clim);
            end
            if ~isempty(options.cbXTickLabel)
                hcb.XTickLabel = options.cbXTickLabel;
            end
            if ~isempty(options.cbXTick)
                hcb.XTick = options.cbXTick;
            end
        end
        
        function [deltaMD, deltaPD, refDay, tuned_feats, mask] = getDeltaTuning(obj, options)
            
            % return angle difference between two PDs (from 0 - 180)
            % modulation depth (MD) difference is the absolute
            % difference from the reference day, not relative.
            
            arguments
                obj PlotTunings
                options.refDay (1,1) double = 1 % set ref day to a specific day
                options.refFirstTune logical = true % set ref day to subtract from
                options.mask_non_sig logical = true % set value to NaN if p-value is < 0.05
            end

            mask = zeros(size(obj.tunings.pval)); 
            mask(obj.tunings.pval>=0.05) = NaN;
            deltaMD = NaN(size(mask));
            deltaPD = NaN(size(mask));

            % compare delta from first session with significant tuning
            % instead of specify one specific session
            if options.refFirstTune

                % find feats that has at least one session with tunings
                tuned_feats = any(~isnan(mask),2);
                    
                % first first session with tuning per feature
                refDay = ones(obj.nfeat,1); % default reference day to day 0
                for f = 1:obj.nfeat
                    if tuned_feats(f)
                        refDay(f) = find(~isnan(mask(f,:)),1,'first');
                    end
                end

                fprintf('\n%i features have significant tuning on at least one session.\n\n', ...
                    sum(tuned_feats))
            else
                % default reference day to refDay
                refDay = repmat(options.refDay, obj.nfeat, 1);

                % find feats where refDay has tunings
                tuned_feats = any(~isnan(mask(:,options.refDay)),2);
                fprintf('\n%i features have significant tuning on session number %i.\n\n', ...
                    sum(tuned_feats), options.refDay)

            end

            % get delta
            for f = 1:obj.nfeat
                deltaMD(f,:) = (obj.tunings.md(f,:) - obj.tunings.md(f,refDay(f)));
                deltaPD(f,:)  = angdiffdeg(obj.tunings.pd(f,:), obj.tunings.pd(f,refDay(f)));
            end

            % apply mask to set value to NaN if p-value is < 0.05
            if options.mask_non_sig
                fprintf('Value is set to NaN if p-value is < 0.05.\n\n')
                deltaMD = deltaMD + mask;
                deltaPD = deltaPD + mask;
            end

        end

        function [data, order] = GetSortedOrder(obj, data, varargin)
            % sort channels either by columns or hierarchical clustering
            
            order = zeros(size(data,1),1);
            % sort per array per feature type
            if mod(size(data,1), 96) == 0
                for i = 1:size(data,1)/96
                    chs = (1:96) + 96*(i-1);
                    tmp = getOrderWithinRange(obj, data(chs,:), i, varargin{:});
                    order(chs) = tmp + 96*(i-1);
                end
            else
            % sort all channels
                order = getOrderWithinRange(obj, data, 1, varargin{:});
            end
            data = data(order,:);
        end

        function order = getOrderWithinRange(obj, data, i, options)
            arguments
                obj (1,1) PlotTunings
                data double
                i 
                options.distMetric = 'squaredeuclidean';
                options.linkage = 'ward';
                options.sortByFirstCol = false;
            end

            if ~options.sortByFirstCol
                % hierarchical clustering
%               Pairwise distance between observations
                data(isnan(data)) = 0;
                try
                    Y = pdist(data,options.distMetric);
                    disp([options.distMetric,' distance was used for pdist.'])
                catch
                    Y = pdist(data,'squaredeuclidean');
                    disp('Squared Euclidean distance was used for pdist.')
                end
                Z = linkage(Y,options.linkage);
                disp([options.linkage,' method was used for linkage.'])
                order = optimalleaforder(Z,Y);%,'CRITERIA','group');
            else 
                % first column is either 0 or NaN (non significant)
                [~, tmp] = sort(data(:, 1));
                ind = tmp(data(tmp,1)==0); [~, tmp1] = sort(data(ind, end));
                ind = tmp(data(tmp,1)~=0); [~, tmp2] = sort(data(ind, end));
                order = [tmp(tmp1); tmp(tmp2 + length(tmp1))];
                if mod(i,2); order = flip(order); end
            end
        end     
    end
end

function delta = angdiffdeg(x, y)
    % absolute difference in angles in degree
    normDeg  = mod(x - y,360);
    delta = min(360-normDeg , normDeg);
end