function plotDecisionTrials(flat, varargin)



% plot trial kinematics with landing position distribution below

% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.outcome = 'isBigStep';  % variable used to color trials // isBigStep or isModPawLengthened

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials
s.deltaMin = 0;  % exclude little step trials where modPawDeltaLength is less than deltaLim standard deviations

s.colors = flipud(colorme(2, 'offset', .2, 'showSamples', false));  % colors for little and big steps
s.ctlStepColor = [.5 .5 .5];
s.obsColor = [188 125 181] / 255;
s.rowColors = [];  % if provided, overwrites s.colors and plots all traces in the same color per row

s.trialsToShow = 50;
s.histoFillAlpha = .2;
s.xLims = [-.11 .08];

s.histoOffset = .005;  % (m)
s.histoHgt = .015;  % (m)

s.saveLocation = '';  % if provided, save figure automatically to this location




% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isstruct(flat); flat = struct2table(flat); end

% restrict to desired trials
if s.successOnly; flat = flat(flat.isTrialSuccess,:); end
if s.modPawOnlySwing; flat = flat(flat.modPawOnlySwing==1,:); end
if s.lightOffOnly; flat = flat(~flat.isLightOn,:); end
if s.deltaMin
    minDif = std(flat.preModPawDeltaLength) * s.deltaMin;
    flat = flat(abs(flat.modPawDeltaLength)>minDif,:);
end


xGrid = linspace(s.xLims(1), s.xLims(2), 500);  % grid for histograms

kinData = permute(cat(3, flat.modPawKinInterp{:}), [3,1,2]);
kinDataCtl = permute(cat(3, flat.preModPawKinInterp{:}), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1);
kinData(flat.firstModPaw==3,2,:) = -kinData(flat.firstModPaw==3,2,:); % flip st mod paw is always paw 2
kinDataCtl(flat.firstModPaw==3,2,:) = -kinDataCtl(flat.firstModPaw==3,2,:); % flip st mod paw is always paw 2

figure('Color', 'white', 'Position', [200 100 900.00 250*length(s.levels)], 'MenuBar', 'none');
if ~isempty(s.condition)  % if no condition provided, create dummy variable
    [~, condition] = ismember(flat.(s.condition), s.levels);  % turn the 'condition' into numbers
else
    condition = ones(height(flat), 1);
end


for i = 1:length(s.levels)
    subplot(length(s.levels), 1, i)
    bins = condition == i;
    
    if ~isempty(s.rowColors)
        s.colors = repmat(s.rowColors(i,:),2,1);
    end
    
    % plot kinematics
    plotKinematics(kinData(bins,[1 3],:), flat.obsHgt(bins), flat.(s.outcome)(bins) + 1, ...
        'colors', s.colors, 'trialsToOverlay', s.trialsToShow, 'trialAlpha', .4, 'lineAlpha', 0, 'yLimZero', false, 'plotObs', false)
    plotKinematics(kinDataCtl(bins,[1 3],:), flat.obsHgt(bins), ones(1,sum(bins)), ...
        'colors', s.ctlStepColor, 'lineWidth', 5, 'yLimZero', false, 'obsColors', s.obsColor)
    set(gca, 'XLim', s.xLims)

    % plot pdfs
    longShortRatio = nanmean(flat.(s.outcome)(bins));
    kdCtl = ksdensity(kinDataCtl(bins,1,end), xGrid);
    kdLong = ksdensity(kinData(bins & flat.(s.outcome)==1,1,end), xGrid) * longShortRatio;
    kdShort = ksdensity(kinData(bins & flat.(s.outcome)~=1,1,end), xGrid) * (1-longShortRatio);

    % scale y axis to fit in same subplot as kinematics
    pdfMax = max([kdLong kdCtl kdShort]);
    kdCtl = -kdCtl * (s.histoHgt/pdfMax) - s.histoOffset;
    kdLong = -kdLong * (s.histoHgt/pdfMax) - s.histoOffset;
    kdShort = -kdShort * (s.histoHgt/pdfMax) - s.histoOffset;


    % plot that shit
    fill([xGrid xGrid(1)], [kdCtl kdCtl(1)], s.ctlStepColor, 'FaceAlpha', s.histoFillAlpha, 'EdgeColor', 'none')
    plot(xGrid, kdCtl, 'Color', s.ctlStepColor, 'LineWidth', 2)

    fill([xGrid xGrid(1)], [kdLong kdLong(1)], s.colors(2,:), 'FaceAlpha', s.histoFillAlpha, 'EdgeColor', 'none')
    plot(xGrid, kdLong, 'Color', s.colors(2,:), 'LineWidth', 2)

    fill([xGrid xGrid(1)], [kdShort kdShort(1)], s.colors(1,:), 'FaceAlpha', s.histoFillAlpha, 'EdgeColor', 'none')
    plot(xGrid, kdShort, 'Color', s.colors(1,:), 'LineWidth', 2)
    
    title(s.levels{i})
end


% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end
