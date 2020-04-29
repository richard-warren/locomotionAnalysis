function heatmaps = plotDecisionHeatmaps(flat, varargin)


% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.plotMice = true;  % whether to show plots for individual mice
s.avgMice = true;  % whether to compute DVs for each mouse, then average // otherwise, pool across mice

s.colormap = 'hot';  % heat map colormap
s.colors = [];  % color(s0 of line showing probability of big step
s.xLims = [-.03 .015]*1000;
s.yLims = [-.03 .03]*1000;
s.binWidth = 10;  % (mm) width for sliding average of big ste probability
s.binNum = 100;  % number of bins for sliding average of big step probability
s.saveLocation = '';  % if provided, save figure automatically to this location
s.outcome = 'isBigStep';  % plot the probability of 'outcome' as a function of 'modPawPredictedDistanceToObs'
s.normalize = '';  % 'row' or 'col' normalize the sum of each row or column to 1
s.lineThickness = .05;  % expressed as fraction of y axis range
s.lineAlhpa = .8;
s.plotProbs = true;  % whether to plot probability of outcome being true
s.subplotDims = [];  % user can overwrite the dimensions of the subplots, although this should only not wreak havock under certain conditions (e.g. when there is only one row or column)

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials
s.deltaMin = 0;  % exclude trials where little step length is adjusted less than deltaMin


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.colors); s.colors = jet(length(s.levels)); end
if isstruct(flat); flat = struct2table(flat); end

x = linspace(s.xLims(1), s.xLims(2), s.binNum);
y = linspace(s.yLims(1), s.yLims(2), s.binNum);


% restrict to desired trials
if s.successOnly; flat = flat(flat.isTrialSuccess,:); end
if s.modPawOnlySwing; flat = flat(flat.modPawOnlySwing==1,:); end
if s.lightOffOnly; flat = flat(~flat.isLightOn,:); end
if s.modSwingContactsMax; flat = flat(flat.modSwingContacts<=s.modSwingContactsMax, :); end
% if s.deltaMin
%     minDif = std(flat.preModPawDeltaLength) * s.deltaMin;
%     flat = flat(abs(flat.modPawDeltaLength)>minDif,:);
% end
if s.deltaMin; flat = flat( ~(abs(zscore(flat.modPawDeltaLength))<s.deltaMin & [flat.isBigStep]==0), :); end



if ~isempty(s.condition)  % if no condition provided, create dummy variable
    [~, condition] = ismember(flat.(s.condition), s.levels);  % turn the 'condition' into numbers
else
    condition = ones(height(flat),1);
end



if ~s.avgMice; flat.mouse = repmat({'temp'}, height(flat), 1); end  % if pooling across mice, rename all mice with dummy string
mice = unique(flat.mouse);


% compute heatmaps and big step probabilities
binCenters = linspace(s.xLims(1), s.xLims(2), s.binNum);
[heatmaps, probs, thicknesses] = deal(cell(length(s.levels), length(mice)));

for i = 1:length(s.levels)
    conditionBins = condition==i;
    fprintf('%s: ', s.levels{i})
    
    for j = 1:length(mice)
        fprintf('%i/%i ', j, length(mice))
        bins = conditionBins(:) & strcmp(flat.mouse, mice{j});
    
        if any(bins)
            heatmaps{i,j} = heatmapRick(flat.modPawPredictedDistanceToObs(bins)*1000, flat.modPawDistanceToObs(bins)*1000, ...
                'xLims', s.xLims, 'yLims', s.yLims, 'showPlot', false, 'binNum', s.binNum, 'normalize', s.normalize);

            % compute moving average for big step probability
            probs{i,j} = zeros(1, length(binCenters));
            for k = 1:length(binCenters)
                binLims = binCenters(k) + [-s.binWidth/2 s.binWidth/2];
                binsSub = bins & ...
                       flat.modPawPredictedDistanceToObs*1000 > binLims(1) & ...
                       flat.modPawPredictedDistanceToObs*1000 <= binLims(2);
                if any(binsSub); probs{i,j}(k) = nanmean(flat.(s.outcome)(binsSub)); end
            end
            thicknesses{i,j} = ksdensity(flat.modPawPredictedDistanceToObs*1000, x);  % compute marginals for x
        else
            heatmaps{i,j} = zeros(s.binNum, s.binNum);
            probs{i,j} = zeros(1, length(binCenters));
        end
    end
    fprintf('\n')
end

    
% MAKE PLOTS
if (s.plotMice && s.avgMice); rows = 1+length(mice); else; rows = 1; end
if isempty(s.subplotDims); s.subplotDims = [rows, length(s.levels)+s.plotProbs]; end
figure('Color', 'white', 'Position', [200 10 200*s.subplotDims(2) 200*s.subplotDims(1)], 'MenuBar', 'none');
colormap(s.colormap)


% average heatmaps
for i = 1:length(s.levels)
    
    % heatmap
    subplot(s.subplotDims(1), s.subplotDims(2), i); hold on
    imagesc(x, y, nanmean(cat(3, heatmaps{i,:}), 3))
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'off')
    line(s.xLims+[-1 1]*range(s.xLims), s.xLims+[-1 1]*range(s.xLims), 'color', [0 0 0 .5], 'lineWidth', 3)
    if i==1; ylabel('average'); end
    set(gca, 'XLim', s.xLims, 'YLim', s.yLims)
    
    % probabilities
    if s.plotProbs
        yyaxis right; hold on
        t = mean(cat(1,thicknesses{i,:}),1);  % line thickness
        t = t*s.lineThickness / max(t);
        p = nanmean(cat(1, probs{i,:}), 1);
        patch([x fliplr(x)], [p+t, fliplr(p-t)], s.colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', s.lineAlhpa)
        set(gca, 'YColor', s.colors(i,:), 'ytick', [0 1], 'ylim', [0 1])
    end
    
    title(s.levels{i})
end

% average log plots
if s.plotProbs
    subplot(s.subplotDims(1), s.subplotDims(2), length(s.levels)+s.plotProbs); hold on
    for i = 1:length(s.levels)
        t = mean(cat(1,thicknesses{i,:}),1);
        t = t*s.lineThickness / max(t);
        p = nanmean(cat(1, probs{i,:}), 1);

        patch([x fliplr(x)], [p+t, fliplr(p-t)], s.colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', s.lineAlhpa)
        set(gca, 'YColor', s.colors(i,:), 'ytick', [0 1], 'ylim', [0 1])
        title([s.outcome ' probability'])
    end
end


% mouse plots
if (s.plotMice && s.avgMice)
    for m = 1:length(mice)

        % heatmaps
        for i = 1:length(s.levels)
            subplot(s.subplotDims(1), s.subplotDims(2), i + m*(length(s.levels)+s.plotProbs))
            imagesc(x, y, heatmaps{i,m})
            set(gca, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'off')
            line(s.xLims+[-1 1]*range(s.xLims), s.xLims+[-1 1]*range(s.xLims), 'color', [0 0 0 .5], 'lineWidth', 3)
            if i==1; ylabel(mice{m}); end

            if s.plotProbs
                yyaxis right
                t = thicknesses{i,m};
                t = t*s.lineThickness / max(t);
                p = probs{i,m};

                patch([x fliplr(x)], [p+t, fliplr(p-t)], s.colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', s.lineAlhpa)
                set(gca, 'YColor', s.colors(i,:), 'ytick', [0 1], 'ylim', [0 1])
            end
        end

        % log plot
        if s.plotProbs
            subplot(s.subplotDims(1), s.subplotDims(2), (m+1)*(length(s.levels)+s.plotProbs)); hold on
            for i = 1:length(s.levels)
                t = thicknesses{i,m};
                t = t*s.lineThickness / max(t);
                p = probs{i,m};

                patch([x fliplr(x)], [p+t, fliplr(p-t)], s.colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', s.lineAlhpa)
                set(gca, 'YColor', s.colors(i,:), 'ytick', [0 1], 'ylim', [0 1])
            end
        end
    end
end

% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end





