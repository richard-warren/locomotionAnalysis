function mouseAvgs = plotDvPsth(data, dv, plotVar, opts)


% used to plot PSTHs of continuous data signals stored in data struct //
% data has one row per trial // dv is the variable to be plotted as a
% function of x (eg plotting velocity vs time) // plotVar is the field in data containing conditions that
% should be plotted as different colored lines // rowVar is plotted as
% different subplot, stacked one on top of another // this function
% averages all trials for given mouse, then averages across mice // assumes
% that dv is matrix with 2 rows, first row containing y values and second
% row containing x values // opts is cell array containing name value pairs
% of settings to adjust // note that all settings are in structure 's'


% settings
s.errorFcn = @(x) nanstd(x)/sqrt(size(x,1)); % function for error bars
% s.errorFcn = @(x) nanstd(x);
s.plotMouseAvgs = false;
s.showLegend = true;
s.mouseColors = 'lines';
s.conditionColors = 'hsv'; % this can also be set directly to a matrix of colors by the user if the number of conditions are known in advance
s.errorAlpha = .2;
s.xlim = [];
s.rowVar = []; % plot levels of this var as different subplot rows
s.plotConditions = {}; % levels of plotVar to show (otherwise determined as all unique elements in plotVar)
s.lineWidth = 3;
s.mouseAlpha = .4;
s.errorEdges = .4;
s.showErrorBars = true;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% create dummy plotvar and rowVar if not provided by user (all ones)
temp = num2cell(ones(size(data)));
if isempty(s.rowVar); [data.rowVar] = temp{:}; rowVar = 'rowVar'; end % creates dummy row containing all ones
if ~exist('plotVar', 'var') || isempty(plotVar); [data.plotVar] = temp{:}; plotVar = 'plotVar'; end

% initializations
if isempty(s.plotConditions)
    if isa(data(1).(plotVar), 'char')
        s.plotConditions = unique({data.(plotVar)});
    else
        s.plotConditions = num2cell(unique([data.(plotVar)]));
    end
end

if ~iscell(s.plotConditions)
    s.plotConditions = num2cell(unique([data.(plotVar)]));
end

if isa(data(1).(rowVar), 'char'); rowConditions = unique({data.(rowVar)}); else; rowConditions = num2cell(unique([data.(rowVar)])); end
mice = unique({data.mouse});
xGrid = data(1).([dv 'X']);
s.mouseColors = eval([s.mouseColors '(length(mice))']);
if ischar(s.conditionColors); s.conditionColors = eval([s.conditionColors '(length(s.plotConditions))']); end % set colorspace if color is specified as a string

% collect data for each mouse in each condition
mouseAvgs = nan(length(s.plotConditions), length(rowConditions), length(mice), length(xGrid)); % plt condition X row condition X mouse X position
dvData = reshape([data.(dv)],[],length(data))'; % turn dv into matrix
for i = 1:length(s.plotConditions)
    for j = 1:length(rowConditions)
        for k = 1:length(mice)
            bins = cellfun(@(x) isequal(x, s.plotConditions{i}), {data.(plotVar)}) & ...
                   cellfun(@(x) isequal(x, rowConditions{j}), {data.(rowVar)}) & ...
                   strcmp({data.mouse}, mice{k});
            mouseAvgs(i,j,k,:) = nanmean(dvData(bins,:),1);
        end
    end
end


% plot 
for i = 1:length(rowConditions)
    if length(rowConditions)>1; subplot(length(rowConditions),1,i); end
    
    for j = 1:length(s.plotConditions)
        
        if s.plotMouseAvgs
            for k = 1:length(mice)
                plot(xGrid, squeeze(mouseAvgs(j,i,k,:)), ...
                    'LineWidth', 1, 'Color', [s.mouseColors(k,:) s.mouseAlpha]); hold on
            end
        end
        
        if s.showErrorBars
            shadedErrorBar(xGrid, squeeze(mouseAvgs(j,i,:,:)), {@nanmean, s.errorFcn}, ...
                'lineprops', {'linewidth', s.lineWidth, 'color', s.conditionColors(j,:)}, 'patchSaturation', s.errorAlpha); hold on;
        else
%             keyboard
            conditionAvgs = squeeze(mouseAvgs(j,i,:,:));
            plot(xGrid, nanmean(conditionAvgs,1), 'linewidth', s.lineWidth, 'color', s.conditionColors(j,:)); hold on;
        end
    end
    
    % pimp fig
    set(gca, 'Box', 'off', 'TickDir', 'out')
    if ~isempty(s.xlim); set(gca, 'XLim', s.xlim); end
end

% add legend
if s.showLegend
    
    for i = 1:length(s.plotConditions); lines(i) = plot([nan nan], 'color', s.conditionColors(i,:), 'LineWidth', s.lineWidth); end % create dummy lines
    if islogical(s.plotConditions{1})
        conditionNames = {[plotVar ': false'], [plotVar ': true']};
        if s.plotConditions{1}; conditionNames = fliplr(conditionNames); end % if [1 0] rather than [0 1], flip condition names
    elseif ischar(s.plotConditions{1})
        conditionNames = s.plotConditions;
    elseif isnumeric(s.plotConditions{1})
        conditionNames = cellfun(@(x) {num2str(x)}, plotConditions);
    end
    
    legend(lines, conditionNames, 'Location', 'best', 'Box', 'off', 'AutoUpdate', 'off');
end






