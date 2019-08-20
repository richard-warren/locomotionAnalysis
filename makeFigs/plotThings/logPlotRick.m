function logPlotRick(x, y, opts)

% plots probability of logical variable being true (y axis) as function of
% scalar (x axis) // x and y are the scalar and logical variable,
% respectively // varNames is a cell array containing names of x, then y
% value // conditions is of length x, indicating the condition identity of
% values is x and y // each condition is plotted as a line of a different
% color // the names of the conditions are stored in conditionNames, a cell
% array of length max(conditions)

% settings
s.xPercentileLims = [5 95]; % exclude x values outside these percentile limits
s.binWidth = []; % expressed as fraction of x limits (determined by xPercentileLims)
s.binNum = 100;
s.lineWidth = 3;
s.conditionNames = {}; % names of conditions for legend
s.xlabel = []; % two element cell array containing names of x and y axes
s.ylabel = [];
s.colors = 'hsv'; % color scheme // can be specified either as a string, or as an nX3 color map, where n is number of conditions
s.conditions = ones(1,length(x));
s.xlim = []; % if specified, only includes bins in this range
s.ylim = []; % if specified, set y range to this
s.smoothing = [];
s.varyTransparency = false; % if true, then transparency reflects density of sampling at a given x position
s.mouseNames = {}; % if this exists (is set in opts), kinematics are first averaged within, then across mice // cell array of mouse name corresponding to each trial
s.plotMice = false; % if true, plot thin lines showing individual mice, in addition to mean across mice (only can do this if mouseNames are provided)
s.mouseTransparency = .4; % transparency of lines for individual mice
s.plotMouseErrors = true; % if mouseNames provided, plot std across mice
s.errorFcn = @(x) nanstd(x); % error function used the error across mice
s.errorAlpha = .2;
s.computeVariance = false; % if true, computes variance of y instead of y


% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end
% keyboard

% initializations
if ischar(s.colors); s.colors = eval([s.colors '(max(s.conditions))']); end % set colorspace if color is specified as a string
if isempty(s.xlim); s.xlim = prctile(x, s.xPercentileLims); end
if isempty(s.binWidth); s.binWidth = .05*range(s.xlim); end % if bin width not specified, just use 5% of x range
binCenters = linspace(s.xlim(1), s.xlim(2), s.binNum);
if ~isempty(s.smoothing); smoothingSmps = s.smoothing/diff(binCenters(1:2)); end
if ~isempty(s.mouseNames)
    [~, rowConditions] = ismember(s.mouseNames, unique(s.mouseNames));
else
    rowConditions = ones(length(x),1);
end
x = x(:); y = y(:); s.conditions = s.conditions(:); rowConditions = rowConditions(:); % ensure orientations are matched
hold on



for i = 1:max(s.conditions)
    
    % collect data for condition, looping across bins (rows correspond to
    % the different mice, if s.mouseNames is provided)
    conditionData = nan(max(rowConditions),s.binNum);
    for r = 1:max(rowConditions)
        for j = 1:s.binNum
            bins = x > binCenters(j)-s.binWidth*.5 & ...
                   x <= binCenters(j)+s.binWidth*.5 & ...
                   s.conditions==i & ...
                   rowConditions==r;
           if s.computeVariance 
                conditionData(r,j) = nanstd(y(bins))^2;
           else
               conditionData(r,j) = nanmean(y(bins));
           end
        end
        if ~isempty(s.smoothing); conditionData(r,:) = smooth(conditionData(r,:), smoothingSmps); end
    end
    
    % plot lines for individual mice
    if ~isempty(s.mouseNames) && s.plotMice
        for r = 1:max(rowConditions)
            plot(binCenters, conditionData(r,:), 'LineWidth', 1, 'Color', [s.colors(i,:) s.mouseTransparency])
        end
    end
    
    % plot condition mean
    if ~isempty(s.mouseNames) && s.plotMouseErrors
        shadedErrorBar(binCenters, conditionData, {@nanmean, s.errorFcn}, ...
            'lineprops', {'linewidth', s.lineWidth, 'color', s.colors(i,:)}, 'patchSaturation', s.errorAlpha); hold on;
    else
        plot(binCenters, nanmean(conditionData,1), 'LineWidth', s.lineWidth, 'Color', s.colors(i,:));
    end
end

set(gca, 'XLim', s.xlim, 'Box', 'off')
if ~isempty(s.xlabel); xlabel(s.xlabel); end
if ~isempty(s.ylabel); ylabel(s.ylabel); end
if ~isempty(s.ylim); set(gca, 'YLim', s.ylim); end
if ~isempty(s.conditionNames)
    for i = 1:length(s.conditionNames); lines(i) = plot([nan nan], 'color', s.colors(i,:), 'LineWidth', 2); end % create dummy lines
    legend(lines, s.conditionNames, 'Box', 'off', 'Location', 'best');
end

