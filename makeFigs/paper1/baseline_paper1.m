%% load experiment data
fprintf('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')


%% ----------
% PLOT THINGS
%  ----------

%% speed vs. position

% settings
yLims = [0 .8];
obsOnColor = [0 0 0];
obsOnAlpha = .05;
meanColor = [0 0 0];

% initializations
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition'});
flat = flat([flat.isLightOn]);
flat = flat(~[flat.isWheelBreak]);
close all; figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off')

% add obstacle rectangle and lines
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
    'FaceColor', [obsOnColor obsOnAlpha], 'EdgeColor', 'none');
line([0 0], yLims, 'linewidth', 2, 'color', obsOnColor)
    
% plot
plotDvPsth(flat, 'velVsPosition', [-.5 .2], [], [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [0 0 0], 'errorFcn', @(x) nanstd(x)})

% pimp fig
set(gca, 'YLim', yLims);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineVel');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% kinematics

% settings
colNames = {'hind limb', 'fore limb'};
rowNames = {'lagging', 'leading'};
isLeading = [false false true true];
isFore = [false true false true];

obsHgtBins = 4; % discretize obstacle heights into this many bins
xLims = [-.05 .05];
yLims = [0 .016];
colors = hsv(4);
fading = .25; % within a plot, each paw's color fades from fading*color -> color

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'isFore', 'isLeading', 'stepOverKinInterp'});

obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLightOn]);


% initializations
conditions = [flat.obsHgtDiscretized];
kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);
kinData = kinData(:,[1,3],:); % keep only x and z dimensions        

%
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [[2017 708 1608 242]])

for i = 1:4
    subplot(2,2,i)
    plotColor = repmat(colors(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    bins = [flat.isLeading]==isLeading(i) & ...
           [flat.isFore]==isFore(i);
    plotKinematics(kinData(bins,:,:), [flat(bins).obsHgt], conditions(bins), ...
        {'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8})
    set(gca, 'XLim', xLims, 'YLim', yLims)
end

% add row, col labels
subplot(2,2,1);
text(0, yLims(2), colNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(xLims(1), mean(yLims), rowNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90);
subplot(2,2,2);
text(0, yLims(2), colNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subplot(2,2,3);
text(xLims(1), mean(yLims), rowNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90);

% !!! reposition figures
% subplot(2,2,1); pos=get(gca,'OuterPosition'); 

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineKinematics');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');















