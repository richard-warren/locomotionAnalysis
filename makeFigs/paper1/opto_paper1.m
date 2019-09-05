%% OPTO ANALYSIS

% first select dataset, then set global settings, then load data, then add powerCondition, then plots things!


%% select dataset

clear all;
exp = 'sen2mm';

switch exp
    
    % mtc with light 2mm above skull
    case 'mtc2mm'
        sessions = {'190822_000', '190822_001', '190822_002', ...
            '190827_000', '190827_001', '190827_002'};
        
    % sen with light 2mm above skull
    case 'sen2mm'
%         sessions = {'190823_000', '190823_001', '190823_002', ...
%             '190828_000', '190828_001', '190828_002'};
        sessions = {'190828_000', '190828_001', '190828_002'};
        
    % olf with light 2mm above skull
    case 'olf2mm'
        sessions = {'190826_000', '190826_001', '190826_002', ...
            '190829_002', '190829_003', '190829_004'};
end

%% global settings
groupPowers = false;
colors = [0 0 0; winter(3)]; % assumes 3 light on conditions
colorsGrouped = [.2 .2 .2; .2 .6 .9];  % light off, light on
varsToAvg = {'mouse'};


% initializations
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'optoNotes');
sessionInfo = sessionInfo(ismember(sessionInfo.session, sessions),:);
powers = cellfun(@str2num, strsplit(sessionInfo.power___{1}, ', '), 'UniformOutput', false);  % NOTE: assumes all sessions share the same power!
powers = [0, cat(2, powers{:})];

if groupPowers
    vars.condition = struct('name', 'isOptoOn', 'levels', [0 1], 'levelNames', {{'no opto', 'opto'}});
    dvName='isOptoOn';
else
    vars.condition = struct('name', 'powerCondition', 'levels', 1:length(powers), 'levelNames', {sprintfc('%.2f', powers)});
    dvName='powerCondition';
end
vars.conditionGrouped = struct('name', 'isOptoOn', 'levels', [0 1], 'levelNames', {{'no opto', 'opto'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
conditionNames = vars.condition.levelNames;

conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isHind = struct('name', 'isFore', 'condition', @(x) x==0);
conditionals.none = struct('name', '', 'condition', @(x) x);  % what does this do???
conditionals.isFast = struct('name', 'velAtWiskContact', 'condition', @(x) x>.3);
     
%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', ['opto_' exp '.mat']), 'data'); disp([exp ' opto data loaded!'])

%% compute experiment data
addToLoadedData = true;
if exist('data', 'var') && addToLoadedData; data = getExperimentData(sessionInfo, 'all', data); else; data = getExperimentData(sessionInfo, 'all'); end
save(fullfile(getenv('OBSDATADIR'), 'matlabData', ['opto_' exp '.mat']), 'data'); disp('data saved')

%% add powerCondition (binned light power)

for i = 1:length(data.data)
    for j = 1:length(data.data(i).sessions)
        sesPowers = [data.data(i).sessions(j).trials.optoPower];
        powerCondition = num2cell(knnsearch(powers', sesPowers'));
        [data.data(i).sessions(j).trials.powerCondition] = deal(powerCondition{:});
    end
end


%% ----------
% PLOT THINGS
%  ----------


%% BARS

% initializations
scatAlpha = .8;
addStats = true;
rows = 4;
cols = 3;
figure('name', exp, 'color', 'white', 'menubar', 'none', 'position', [1984 113 916 871])
figConditionals = [conditionals.none];

% success
subplot(rows, cols, 1);
conditions = [vars.condition];
dv = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'success rate', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'ylim', [0 1], 'ytick', 0:.5:1, ...
    'compareToFirstOnly', true, 'smpNames', {data.data.mouse}, 'numVariables', length(conditions)})

% vel
subplot(rows, cols, 2);
conditions = [vars.condition];
dv = getDvMatrix(data, 'velAtWiskContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'velocity (m/s)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), 'ylim', [0 .8]...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats})

% tail height
subplot(rows, cols, 3);
conditions = [vars.condition];
dv = getDvMatrix(data, 'tailHgtAtWiskContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'tail height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), 'ylim', [], ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'numVariables', length(conditions)})

% pre obs height
subplot(rows, cols, 4);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'pre paw height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'numVariables', length(conditions)})

% max obs height
subplot(rows, cols, 5);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'stepOverMaxHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'max paw height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'numVariables', length(conditions)})

% control step height
subplot(rows, cols, 6);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'controlStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'control step height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'numVariables', length(conditions)})

% height shaping (correlations)
subplot(rows, cols, 7);
conditions = [vars.conditionGrouped];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, ...
    [conditionals.isFore; conditionals.isLeading; figConditionals], 'corr');
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'paw shaping (correlation)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', colorsGrouped, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'numVariables', length(conditions)})


% compute decision determinism and threshold
flat = flattenData(data, {'mouse', 'modPawPredictedDistanceToObs', 'conditionNew', 'isBigStep', 'isOptoOn', 'powerCondition'});
mice = unique({flat.mouse});
lightConditions = [false, true];
[accuracies, thresholds, f1s] = deal(nan(length(conditionNames), length(mice)));
conditions = [vars.condition];

for m = 1:length(conditions.levels)
    for j = 1:length(mice)
        bins = [flat.(dvName)]==conditions.levels(m) & ...
                strcmp({flat.mouse}, mice{j});
        glm = fitglm([flat(bins).modPawPredictedDistanceToObs]', [flat(bins).isBigStep]', 'Distribution', 'binomial');
        coeffs = glm.Coefficients.Estimate;
        predictions = round(predict(glm, [flat(bins).modPawPredictedDistanceToObs]'));
        confusion = confusionmat([flat(bins).isBigStep], logical(predictions));
        precision = confusion(2,2)/sum(confusion(:,2));
        recall = confusion(2,2)/sum(confusion(2,:));

        f1s(m,j) = harmmean([precision, recall]);
        accuracies(m,j) = mean(predictions==[flat(bins).isBigStep]');
        thresholds(m,j) = (.5-coeffs(1)) / coeffs(2); % solve for prediction = .5
    end
end

% decision determinism (glm accuracy)
% !!! how to apply fig conditionals to flattened data?
subplot(rows, cols, 8);
barPlotRick(f1s, {'conditionNames', {conditionNames}, 'ylabel', 'f1 score', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'ylim', [.5 1], 'numVariables', length(conditions)})

% decision threshold
% !!! how to apply fig conditionals to flattened data?
subplot(rows, cols, 9);
barPlotRick(thresholds*1000, {'conditionNames', {conditionNames}, 'ylabel', 'big step threshold (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'numVariables', length(conditions)})

% ventral touches
subplot(rows, cols, 10);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'isVentralContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'ventral touch probability', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'numVariables', length(conditions)})

% dorsal touches
subplot(rows, cols, 11);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'isDorsalContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'dorsal touch probability', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', scatAlpha, 'showStats', addStats, 'numVariables', length(conditions)})


% save
file = fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_bars']);
savefig(gcf, file);


%% SPEED VS. POSITION

% settings
yLims = [0 1];
obsOnColor = [0 0 0];
obsOnAlpha = .05;

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition', 'isOptoOn', 'optoOnPositions', 'optoPower', 'powerCondition'});
figure('name', exp, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300]); hold on

% add light on, obs on, and obs at nose markers
optoOnPos = nanmedian([flat.optoOnPositions]);
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
patch([x(1) x(1) x(2) x(2)], [yLims(1) yLims(2) yLims(2) yLims(1)], ...
    [obsOnColor], 'FaceAlpha', obsOnAlpha, 'EdgeColor', 'none')
line([0 0], yLims, 'linewidth', 2, 'color', obsOnColor)
line([optoOnPos optoOnPos], yLims, 'linewidth', 2, 'color', colors(2,:))
    
% plot
plotDvPsth(flat, 'velVsPosition', dvName, ...
    {'showLegend', false, 'errorFcn', @(x) nanstd(x), 'conditionColors', colors, ...
    'plotMouseAvgs', false})

% pimp fig
set(gca, 'YLim', yLims);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

% save
file = fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_speed']);
savefig(gcf, file);

% INDIVIDUAL MICE PLOTS
mice = unique({flat.mouse});
figure('name', exp, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 350 200*length(mice) 150], 'inverthardcopy', 'off'); hold on

for i = 1:length(mice)
    subplot(1,length(mice),i)
    flatSub = flat(strcmp({flat.mouse}, mice{i}));
    
    plotDvPsth(flatSub, 'velVsPosition', dvName, ...
        {'showLegend', false, 'errorFcn', @(x) nanstd(x), 'conditionColors', colors, ...
        'plotMouseAvgs', false})
    set(gca, 'ylim', yLims)
    
    optoOnPos = nanmedian([flatSub.optoOnPositions]);
    x = [nanmean([flatSub.obsOnPositions]) nanmean([flatSub.obsOffPositions])];
    rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
        'FaceColor', [obsOnColor obsOnAlpha], 'EdgeColor', 'none');
    line([0 0], yLims, 'linewidth', 2, 'color', obsOnColor)
    line([optoOnPos optoOnPos], yLims, 'linewidth', 2, 'color', colors(2,:))
    title(mice{i})
end


%% HEIGHT SHAPING


% SCATTER ACROSS ALL ANIMALS

% settings
xLims = [3 10];
isHgtPreObs = true; % do you measure the height at peak (false) or before the paw reaches the obstacle (true)

isLeading = [true false true false]; % sequence of conditions for plots
isFore = [true true false false];
conditionNames = {'leading fore', 'trailing fore', 'leading hind', 'trailing hind'};

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'powerCondition', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isLeading', 'stepOverMaxHgt', 'isPawSuccess', 'isOptoOn'});
flat = flat(~isnan([flat.obsHgt]) & ...
            ~isnan([flat.stepOverMaxHgt]) & ...
            ~isnan([flat.preObsHgt])); % add conditionals here

% get condition numbers, where each condition is unique combination of isLeading and isFore
a = [isLeading; isFore]'; % each row is a condition (00,01,10,11)
b = [[flat.isLeading]; [flat.isFore]]';
[~, stepTypeConditions] = ismember(b, a, 'rows');

optoCondition = [flat.isOptoOn]+1;
% optoCondition = [flat.powerCondition];
obsHgts = [flat.obsHgt]*1000;
pawHgts = [flat.preObsHgt]*1000;
% pawHgts = [flat.stepOverMaxHgt]*1000;


% OBS HGT PAW HGT MOVING AVG
figure('name', exp, 'Color', 'white', 'Position', [2000 94 786 706], 'MenuBar', 'none');

for i = 1:4
    
    subplot(2,2,i);
    bins = i==stepTypeConditions;
    
    plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
    logPlotRick(obsHgts(bins), pawHgts(bins), ...
        {'colors', colorsGrouped, 'conditions', optoCondition(bins), 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)', ...
        'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat(bins).mouse}, 'plotMice', false, ...
        'errorFcn', @(x) std(x)/sqrt(size(x,1))})
    
    set(gca, 'xlim', [4 10])    
    title(conditionNames{i})
end

% save
file = fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_shaping']);
savefig(gcf, file);

%% KINEMATICS

% settings
obsHgtBins = 3; % discretize obstacle heights into this many bins
fading = .25; % within a plot, each paw's color fades from fading*color -> color
xLims = [-.06 0];
yLims = [0 .016];

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'powerCondition', 'isWheelBreak', ...
    'obsHgt', 'isFore', 'isLeading', 'preObsKin', 'controlStepKinInterp', 'powerCondition', 'isOptoOn'});
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLeading] & ...
            [flat.isFore]);

kinData = permute(cat(3, flat.preObsKin), [3,1,2]);
kinDataCtl = permute(cat(3, flat.controlStepKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1); % change the x starting x position of ctl steps to match steps over


figure('name', exp, 'position', [2000 200 400 800], 'color', 'white', 'menubar', 'none'); hold on;
for i = 1:4
    subplot(5,1,i)
    bins = [flat.powerCondition]==i;
    plotColor = repmat(colors(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], [flat(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims)
end

subplot(5,1,5)
plotKinematics(kinData(:,[1,3],:), [flat.obsHgt], [flat.powerCondition], ...
    {'colors', colors, 'obsAlpha', .25, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims, 'YLim', yLims)

% save
file = fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_kinematics']);
savefig(gcf, file);







