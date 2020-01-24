function plotDecisionHeatmaps(data, varargin)


% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot

s.colors = [];  % color(s0 of line showing probability of big step
s.xLims = [-.03 .015]*1000;
s.yLims = [-.03 .03]*1000;
s.binWidth = 5;  % (mm) width for sliding average of big ste probability
s.binNum = 100;  % number of bins for sliding average of big step probability
s.saveLocation = '';  % if provided, save figure automatically to this location

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.colors); s.colors = jet(length(s.levels)); end

flat = flattenData(data, {'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'isBigStep', ...
    'isLightOn', 'isTrialSuccess', 'modPawOnlySwing', s.condition});

% restrict to desired trials
if s.successOnly; flat = flat([flat.isTrialSuccess]); end
if s.modPawOnlySwing; flat = flat([flat.modPawOnlySwing]==1); end
if s.lightOffOnly; flat = flat(~[flat.isLightOn]); end

[~, condition] = ismember({flat.(s.condition)}, s.levels);  % turn the 'condition' into numbers
figure('Color', 'white', 'Position', [2006 540 300*(length(s.levels)+1) 250], 'MenuBar', 'none');




% make one plot per condition
binCenters = linspace(s.xLims(1), s.xLims(2), s.binNum);
bigStepProbs = nan(length(s.levels), length(binCenters));

for i = 1:length(s.levels)
    
    subplot(1, length(s.levels)+1, i)
    bins = condition==i;
    
    if i==1; ylab = 'landing distance (mm)'; else; ylab = ''; end
    heatmapRick([flat(bins).modPawPredictedDistanceToObs]*1000, [flat(bins).modPawDistanceToObs]*1000, ...
        'xLims', s.xLims, 'yLims', s.yLims, 'colormap', 'hot', ...
        'xlabel', 'predicted landing distance (mm)', 'ylabel', ylab)
    set(gca, 'DataAspectRatio', [1 1 1])
    line(s.xLims+[-1 1]*range(s.xLims), s.xLims+[-1 1]*range(s.xLims), 'color', [0 0 0 .5], 'lineWidth', 3)

    % add moving average of big step probability
    for j = 1:length(binCenters)
        binLims = binCenters(j) + [-s.binWidth/2 s.binWidth/2];
        binsSub = bins & ...
               [flat.modPawPredictedDistanceToObs]*1000 > binLims(1) & ...
               [flat.modPawPredictedDistanceToObs]*1000 <= binLims(2);
        bigStepProbs(i,j) = nanmean([flat(binsSub).isBigStep]);
    end

    yyaxis right
    plot(binCenters, bigStepProbs(i,:), 'LineWidth', 3, 'Color', s.colors(i,:))
    set(gca, 'YColor', s.colors(i,:), 'YTick', 0:.5:1, 'box', 'off', 'ylim', [0 1])
    if i==length(s.levels); ylabel('big step probability'); end
    
    title(s.levels{i});
end

% overlay log plots
subplot(1, length(s.levels)+1, length(s.levels)+1); hold on
for i = 1:length(s.levels)
    plot(bigStepProbs(i,:), 'Color', s.colors(i,:), 'LineWidth', 2)
end


% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end





