function plotDecisionTrials(flat, varargin)



% plot trial kinematics with landing position distribution below

% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.outcome = 'isBigStep';  % variable used to color trials // isBigStep or isModPawLengthened
s.view = 'top';  % 'top' or 'bot'

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
s.showTitles = true;

s.showHistos = true;
s.poolHistos = false;  % whether to pool the big step and little step histos
s.histoOffset = .005;  % (m)
s.histoHgt = .015;  % (m)

s.saveLocation = '';  % if provided, save figure automatically to this location




% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isstruct(flat); flat = struct2table(flat); end

% set view specific parameters
switch s.view
    case 'top'
        isBotView = false;
        dims = [1 3];
        yLims = [-s.histoOffset .02];
    case 'bot'
        isBotView = true;
        dims = [1 2];
        yLims = [-.02 .02];
end

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


xGrid = linspace(s.xLims(1), s.xLims(2), 500);  % grid for histograms

kinData = permute(cat(3, flat.modPawKinInterp{:}), [3,1,2]);
kinDataCtl = permute(cat(3, flat.preModPawKinInterp{:}), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1);
kinData(flat.firstModPaw==3,2,:) = -kinData(flat.firstModPaw==3,2,:); % flip st mod paw is always paw 2
kinDataCtl(flat.firstModPaw==3,2,:) = -kinDataCtl(flat.firstModPaw==3,2,:); % flip st mod paw is always paw 2

figure('Color', 'white', 'Position', [200 0 900.00 250*length(s.levels)], 'MenuBar', 'none');
if ~isempty(s.condition)  % if no condition provided, create dummy variable
    [~, condition] = ismember(flat.(s.condition), s.levels);  % turn the 'condition' into numbers
else
    condition = ones(height(flat), 1);
end


for i = 1:length(s.levels)
    subplot(length(s.levels), 1, i)
    bins = condition == i;
    
%     if ~isempty(s.rowColors); s.ctlStepColor = s.rowColors(i,:); end
    if ~isempty(s.rowColors); s.colors = repmat(s.rowColors(i,:),2,1); end
    
    % plot kinematics
    plotKinematics(kinData(bins,dims,:), flat.obsHgt(bins), flat.(s.outcome)(bins) + 1, ...
        'colors', s.colors, 'trialsToOverlay', s.trialsToShow, 'trialAlpha', .4, 'lineAlpha', 0, ...
        'yLimZero', false, 'plotObs', false, 'isBotView', isBotView)
    plotKinematics(kinDataCtl(bins,dims,:), flat.obsHgt(bins), ones(1,sum(bins)), ...
        'colors', s.ctlStepColor, 'lineWidth', 5, 'yLimZero', false, 'obsColors', s.obsColor, ...
        'isBotView', isBotView)
    set(gca, 'XLim', s.xLims, 'YLim', yLims)

    % plot pdfs
    if s.showHistos
        longShortRatio = nanmean(flat.(s.outcome)(bins));
        kdCtl = ksdensity(kinDataCtl(bins,1,end), xGrid);
        if ~s.poolHistos
            kdLong = ksdensity(kinData(bins & flat.(s.outcome)==1,1,end), xGrid) * longShortRatio;
            kdShort = ksdensity(kinData(bins & flat.(s.outcome)~=1,1,end), xGrid) * (1-longShortRatio);
        else
            kd = ksdensity(kinData(bins,1,end), xGrid);
        end

        % scale y axis to fit in same subplot as kinematics
        if ~s.poolHistos
            pdfMax = max([kdLong kdCtl kdShort]);
        else
            pdfMax = max([kd kdCtl]);
        end

        kdCtl = yLims(1) - kdCtl * (s.histoHgt/pdfMax);
        if ~s.poolHistos
            kdLong = yLims(1) - kdLong * (s.histoHgt/pdfMax);
            kdShort = yLims(1) - kdShort * (s.histoHgt/pdfMax);
        else
            kd = yLims(1) - kd * (s.histoHgt/pdfMax);
        end


        % plot that shit
        fill([xGrid xGrid(1)], [kdCtl kdCtl(1)], s.ctlStepColor, 'FaceAlpha', s.histoFillAlpha, 'EdgeColor', 'none')
        plot(xGrid, kdCtl, 'Color', s.ctlStepColor, 'LineWidth', 2)

        if ~s.poolHistos
            fill([xGrid xGrid(1)], [kdLong kdLong(1)], s.colors(2,:), 'FaceAlpha', s.histoFillAlpha, 'EdgeColor', 'none')
            plot(xGrid, kdLong, 'Color', s.colors(2,:), 'LineWidth', 2)

            fill([xGrid xGrid(1)], [kdShort kdShort(1)], s.colors(1,:), 'FaceAlpha', s.histoFillAlpha, 'EdgeColor', 'none')
            plot(xGrid, kdShort, 'Color', s.colors(1,:), 'LineWidth', 2)
        else
            fill([xGrid xGrid(1)], [kd kd(1)], s.colors(1,:), 'FaceAlpha', s.histoFillAlpha, 'EdgeColor', 'none')
            plot(xGrid, kd, 'Color', s.colors(1,:), 'LineWidth', 2)
        end

        if s.showTitles; title(s.levels{i}); end

        set(gca, 'ylim', [yLims(1)-s.histoHgt yLims(2)])
    end
end


% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end
