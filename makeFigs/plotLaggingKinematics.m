function plotLaggingKinematics(flat, varargin)



% plot trial kinematics of final control step for lagging forepaw // idea
% is to see how distance of planting paw to obs is distributed across
% trials and how/whether it varies from what you'd expect from chance...

% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials

s.colors = flipud(colorme(2, 'offset', .2, 'showSamples', false));  % colors for little and big steps
s.ctlStepColor = [.5 .5 .5];
s.obsColor = [188 125 181] / 255;

s.trialsToShow = 50;
s.histoFillAlpha = .2;
s.xLims = [-.15 0];

s.histoOffset = .005;  % (m)
s.histoHgt = .015;  % (m)

s.randSeed = [];  % for selecting the same random trials to show

s.saveLocation = '';  % if provided, save figure automatically to this location




% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isstruct(flat); flat = struct2table(flat); end
if ~isempty(s.randSeed); rng(s.randSeed); end  % initialize random seed for reproduceibility

% restrict to desired trials
if s.successOnly; flat = flat(flat.isTrialSuccess,:); end
if s.modPawOnlySwing; flat = flat(flat.modPawOnlySwing==1,:); end
if s.lightOffOnly; flat = flat(~flat.isLightOn,:); end


xGrid = linspace(s.xLims(1), s.xLims(2), 500);  % grid for histograms
kinData = permute(cat(3, flat.laggingPenultKin{:}), [3,1,2]);

figure('Color', 'white', 'Position', [2001.00 10 900.00 180*(length(s.levels)+1)], 'MenuBar', 'none');
if ~isempty(s.condition)  % if no condition provided, create dummy variable
    [~, condition] = ismember(flat.(s.condition), s.levels);  % turn the 'condition' into numbers
else
    condition = ones(height(flat), 1);
end


pdfs = nan(length(s.levels), length(xGrid));
landingPositions = cellfun(@(x) x(1,end), flat.laggingPenultKin);

for i = 1:length(s.levels)
    subplot(length(s.levels)+1, 1, i)
    bins = condition == i;
    
    % plot kinematics
    plotKinematics(kinData(bins,[1 3],:), flat.obsHgt(bins), ones(1,sum(bins)), ...
        'colors', s.colors(i,:), 'trialsToOverlay', s.trialsToShow, 'trialAlpha', .4, 'lineAlpha', 0, 'yLimZero', false, 'obsColors', s.obsColor)
    set(gca, 'XLim', s.xLims)
    
    % compute and plot pdf
    subplot(length(s.levels)+1, 1, length(s.levels)+1); hold on
    pdfs(i,:) = ksdensity(landingPositions(bins), xGrid);
    fill([xGrid(1) xGrid xGrid(end)], [0 pdfs(i,:) 0], s.colors(i,:), 'FaceAlpha', s.histoFillAlpha, 'EdgeColor', 'none')
    plot(xGrid, pdfs(i,:), 'Color', s.colors(i,:), 'LineWidth', 2)
end
set(gca, 'XLim', s.xLims, 'YColor', 'none', 'XTick', [])



% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end
