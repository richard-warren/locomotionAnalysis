function plotDvPsth(data, dv, xLims, plotVar, rowVar)


% used to plot PSTHs of continuous data signals stored on data struct //
% data has one row per trial // dv is the variable to be plotted as a
% function of x (eg plotting velocity vs time) // plotVar is the field in data containing conditions that
% should be plotted as different colored lines // rowVar is plotted as
% different subplot, stacked one on top of another // this function
% averages all trials for given mouse, then averages across mice // assumes
% that dv is matrix with 2 rows, first row containing y values and second
% row containing x values

% temp
% data = getNestedStructFields(dataTemp, ...
%     {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', 'isBigStep', 'obsOnPositions', ...
%     'condition', 'velContinuousAtContact', 'velVsPosition', 'isTrialSuccess', 'isWheelBreak', 'wiskContactPosition'});
% dv = 'velVsPosition';
% plotVar = 'isLightOn';
% rowVar = 'condition';
% xLims = [-.5 .2];

% settings
errorFcn = @(x) nanstd(x)/sqrt(size(x,1)); % function for error bars
plotMouseAvgs = false;
% errorFcn = @(x) nanstd(x);


% initializations
if isa(data(1).(plotVar), 'char'); plotConditions = unique({data.(plotVar)}); else; plotConditions = num2cell(unique([data.(plotVar)])); end
if isa(data(1).(rowVar), 'char'); rowConditions = unique({data.(rowVar)}); else; rowConditions = num2cell(unique([data.(rowVar)])); end
colors = hsv(length(plotConditions));
mice = unique({data.mouse});
xGrid = data(1).(dv)(2,:); % this will fail if first trial has not been


% collect data for each mouse in each condition for both light on and light off trials
mouseAvgs = nan(length(plotConditions), length(rowConditions), length(mice), length(xGrid)); % condition X mouse X light off/on X position
for i = 1:length(plotConditions)
    for j = 1:length(rowConditions)
        for k = 1:length(mice)
            inds = find(cellfun(@(x) isequal(x, plotConditions{i}), {data.(plotVar)}) & ...
                        cellfun(@(x) isequal(x, rowConditions{j}), {data.(rowVar)}) & ...
                        strcmp({data.mouse}, mice{k}));
            
            mouseData = nan(length(inds), length(xGrid));
            for m = 1:length(inds)
                if ~isempty(data(inds(m)).(dv)); mouseData(m,:) = data(inds(m)).(dv)(1,:); end
            end
            mouseAvgs(i,j,k,:) = nanmean(mouseData,1);
        end
    end
end


% plot for light on/off
for i = 1:length(rowConditions)
    subplot(length(rowConditions),1,i)
    
    for j = 1:length(plotConditions)
        shadedErrorBar(xGrid, squeeze(mouseAvgs(j,i,:,:)), {@nanmean, errorFcn}, ...
            'lineprops', {'linewidth', 3, 'color', colors(j,:)}, 'patchSaturation', .1); hold on;
        
        if plotMouseAvgs
            for k = 1:length(mice)
                plot(xGrid, squeeze(mouseAvgs(j,i,k,:)), ...
                    'LineWidth', 1, 'Color', [colors(j,:) .4]); hold on
            end
        end
    end
    yLims = get(gca, 'YLim');
    line([0 0], yLims, 'color', [.5 .5 .5])
    set(gca, 'YLim', yLims)
    
    % pimp fig
    set(gca, 'XLim', xLims, 'Box', 'off')
end

for i = 1:length(plotConditions); lines(i) = plot([nan nan], 'color', colors(i,:), 'LineWidth', 2); end % create dummy lines
try; legend(lines, plotConditions, 'Location', 'northeast', 'Box', 'off', 'AutoUpdate', 'off'); catch; end






