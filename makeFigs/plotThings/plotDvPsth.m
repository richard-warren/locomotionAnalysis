function plotDvPsth(data, dv, xLims, plotVar, rowVar, opts)


% used to plot PSTHs of continuous data signals stored on data struct //
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

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% create dummy plotvar and rowVar if not provided by user (all ones)
temp = num2cell(ones(size(data)));
if ~exist('rowVar', 'var') || isempty(rowVar); [data.rowVar] = temp{:}; rowVar = 'rowVar'; end % creates dummy row containing all ones
if ~exist('plotVar', 'var') || isempty(plotVar); [data.plotVar] = temp{:}; plotVar = 'plotVar'; end

% initializations
if isa(data(1).(plotVar), 'char'); plotConditions = unique({data.(plotVar)}); else; plotConditions = num2cell(unique([data.(plotVar)])); end
if isa(data(1).(rowVar), 'char'); rowConditions = unique({data.(rowVar)}); else; rowConditions = num2cell(unique([data.(rowVar)])); end
mice = unique({data.mouse});
xGrid = data(1).([dv 'X']);
s.mouseColors = eval([s.mouseColors '(length(mice))']);
if ischar(s.conditionColors); s.conditionColors = eval([s.conditionColors '(length(plotConditions))']); end % set colorspace if color is specified as a string


% collect data for each mouse in each condition
mouseAvgs = nan(length(plotConditions), length(rowConditions), length(mice), length(xGrid)); % condition X mouse X light off/on X position
dvData = reshape([data.(dv)],[],length(data))'; % turn dv into matrix
for i = 1:length(plotConditions)
    for j = 1:length(rowConditions)
        for k = 1:length(mice)
            bins = cellfun(@(x) isequal(x, plotConditions{i}), {data.(plotVar)}) & ...
                   cellfun(@(x) isequal(x, rowConditions{j}), {data.(rowVar)}) & ...
                   strcmp({data.mouse}, mice{k});
            mouseAvgs(i,j,k,:) = nanmean(dvData(bins,:),1);
        end
    end
end


% plot 
for i = 1:length(rowConditions)
    if length(rowConditions)>1; subplot(length(rowConditions),1,i); end
    
    for j = 1:length(plotConditions)
        if length(mice)>1
            shadedErrorBar(xGrid, squeeze(mouseAvgs(j,i,:,:)), {@nanmean, s.errorFcn}, ...
                'lineprops', {'linewidth', 3, 'color', s.conditionColors(j,:)}, 'patchSaturation', s.errorAlpha); hold on;
        else
            plot(xGrid, squeeze(mouseAvgs(j,i,:,:)), 'linewidth', 3, 'color', s.conditionColors(j,:)); hold on;
        end
        
        if s.plotMouseAvgs
            for k = 1:length(mice)
                plot(xGrid, squeeze(mouseAvgs(j,i,k,:)), ...
                    'LineWidth', 1, 'Color', [s.mouseColors(k,:) .4]); hold on
            end
        end
    end
    
    % pimp fig
    set(gca, 'XLim', xLims, 'Box', 'off', 'TickDir', 'out')
end

% add legend
if s.showLegend
    
    for i = 1:length(plotConditions); lines(i) = plot([nan nan], 'color', colors(i,:), 'LineWidth', 2); end % create dummy lines
    if islogical(plotConditions{1})
        conditionNames = {[plotVar ': false'], [plotVar ': true']};
        if plotConditions{1}; conditionNames = fliplr(conditionNames); end % if [1 0] rather than [0 1], flip condition names
    elseif ischar(plotConditions{1})
        conditionNames = plotConditions;
    elseif isnumeric(plotConditions{1})
        conditionNames = cellfun(@(x) {num2str(x)}, plotConditions);
    end
    
    legend(lines, conditionNames, 'Location', 'best', 'Box', 'off', 'AutoUpdate', 'off');
end






