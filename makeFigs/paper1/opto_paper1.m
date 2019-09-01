%% INITIALIZATIONS


% global settings
% note: assumes all sesssions have same brain region and light powers

% sessions = {'190822_000', '190822_001', '190822_002', ...
%             '190827_000', '190827_001', '190827_002', 190827_003};  % mtc, 2mm dorsal
% sessions = {'190823_000', '190823_001', '190823_002', ...
%             '190828_000', '190828_001', '190828_002', '190828_003'};  % sen, 2mm dorsal
sessions = {'190826_000', '190826_001', '190826_002', ...
            '190829_001', '190829_002', '190829_003', '190829_004'};  % olf, 2mm dorsal
sessions = '190829_001';  % single session test

groupPowers = false;

colors = [.2 .2 .2; .2 .6 .9];  % light off, light on
varsToAvg = {'mouse'};
matchPropensities = false;
varsToMatch = {'velAtWiskContact', 'tailHgtAtWiskContact'};
manipPercent = 25; % take manipPercent percent of best matched manip trials


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
colorsPowers = [0 0 0; winter(length(vars.condition.levels)-1)]; % scale linearly from black to color
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
if matchPropensities; fileSuffix = '_matched'; else; fileSuffix = ''; end
conditionNames = vars.condition.levelNames;

conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isHind = struct('name', 'isFore', 'condition', @(x) x==0);
conditionals.none = struct('name', '', 'condition', @(x) x);  % what does this do???
conditionals.isFast = struct('name', 'velAtWiskContact', 'condition', @(x) x>.3);

exp = unique(sessionInfo.brainRegion);
exp = exp{:};
     
%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', [exp '_opto_data.mat']), 'data'); disp([exp ' opto data loaded!'])
% data = data(strcmp({data.mouse}, 'vgt6')); varsToAvg = {'session'}; % run this line to show bars for all sessions of a given mouse!

%% compute experiment data
loadOldData = false;
if exist('data', 'var') && loadOldData; data = getExperimentData(sessionInfo, 'all', data); else; data = getExperimentData(sessionInfo, 'all'); end
save(fullfile(getenv('OBSDATADIR'), 'matlabData', [exp '_opto_data.mat']), 'data'); disp('data saved')

%% compute experiment from scratch, in parallel
data = cell(1,length(mice));    
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data = cat(2,data{:});
save(fullfile(getenv('OBSDATADIR'), 'matlabData', [exp '_opto_data.mat']), 'data'); disp('all done!')

%% add powerCondition (binned light power)

for i = 1:length(data.data)
    for j = 1:length(data.data(i).sessions)
        sesPowers = [data.data(i).sessions(j).trials.optoPower];
        powerCondition = num2cell(knnsearch(powers', sesPowers'));
        [data.data(i).sessions(j).trials.powerCondition] = deal(powerCondition{:});
    end
end


%% propensity score matching
if matchPropensities

    mice = {data.data.mouse};
    flat = struct2table(flattenData(data, [{'mouse', 'session', 'trial', 'isOptoOn'} varsToMatch]));
    varBins = ismember(flat.Properties.VariableNames, varsToMatch);
    metaBins = ismember(flat.Properties.VariableNames, {'mouse', 'session', 'trial'});

    % find matched trials
    matchedTrials = cell2table(cell(0,3), 'VariableNames', {'mouse', 'session', 'trial'});
    for mouse = mice
        flatSub = flat(strcmp(flat.mouse, mouse{1}),:);
        X = table2array(flatSub(:, varBins));
        y = flatSub.isOptoOn; % is trial in the manip condition
        matchedPairs = propensityMatching(X, y, ...
            {'percentileThresh', manipPercent, 'predictorNames', varsToMatch, 'verbose', false});
        matchedTrials = [matchedTrials; flatSub(matchedPairs(:), metaBins)];
    end

    % get rid of non-matched trials!
    for i = 1:length(data.data)
        for j = 1:length(data.data(i).sessions)
            bins = strcmp(matchedTrials.mouse, data.data(i).mouse) & ...
                   strcmp(matchedTrials.session, data.data(i).sessions(j).session);
            sesTrials = matchedTrials.trial(bins);
            data.data(i).sessions(j).trials = data.data(i).sessions(j).trials(sesTrials);
        end

        % remove unused sessions
        isSessionUsed = ismember({data.data(i).sessions.session}, unique(matchedTrials.session));
        data.data(i).sessions = data.data(i).sessions(isSessionUsed);
    end
end



%% ----------
% PLOT THINGS
%  ----------


%% BARS

% initializations
close all
addStats = true;
rows = 4;
cols = 3;
figure('name', exp, 'color', 'white', 'menubar', 'none', 'position', [1984 393 1019 591])
figConditionals = [conditionals.none];

% success
subplot(rows, cols, 1);
conditions = [vars.condition];
dv = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'success rate', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'ylim', [0 1], 'ytick', 0:.5:1, ...
    'compareToFirstOnly', true, 'smpNames', {data.data.mouse}, 'numVariables', length(conditions)})

% vel
subplot(rows, cols, 2);
conditions = [vars.condition];
dv = getDvMatrix(data, 'velAtWiskContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'velocity (m/s)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,2,1), 'ylim', [0 .6]...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats})

% angle
% subplot(rows, cols, 3);
% conditions = [vars.condition];
% dv = getDvMatrix(data, 'angleAtWiskContactContra', conditions, varsToAvg, figConditionals);
% barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', ['body angle contra (' char(176) ')'], ...
%     'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,2,1), 'ylim', [-15 15], ...
%     'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% tail height
subplot(rows, cols, 3);
conditions = [vars.condition];
dv = getDvMatrix(data, 'tailHgtAtWiskContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'tail height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,2,1), 'ylim', [], ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'numVariables', length(conditions)})

% pre obs height
subplot(rows, cols, 4);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'pre paw height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'numVariables', length(conditions)})

% max obs height
subplot(rows, cols, 5);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'stepOverMaxHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'max paw height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'numVariables', length(conditions)})

% control step height
subplot(rows, cols, 6);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'controlStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'control step height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'numVariables', length(conditions)})

% height shaping (correlations)
subplot(rows, cols, 7);
conditions = [vars.condition];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, ...
    [conditionals.isFore; conditionals.isLeading; figConditionals], 'corr');
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'paw shaping (correlation)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,8,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'numVariables', length(conditions)})


% compute decision determinism and threshold
flat = flattenData(data, {'mouse', 'modPawPredictedDistanceToObs', 'conditionNew', 'isBigStep', 'isOptoOn', 'powerCondition'});
mice = unique({flat.mouse});
lightConditions = [false, true];
[accuracies, thresholds, f1s] = deal(nan(length(conditionNames), length(mice)));

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
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'ylim', [.5 1], 'numVariables', length(conditions)})

% decision threshold
% !!! how to apply fig conditionals to flattened data?
subplot(rows, cols, 9);
barPlotRick(thresholds*1000, {'conditionNames', {conditionNames}, 'ylabel', 'big step threshold (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'numVariables', length(conditions)})

% ventral touches
subplot(rows, cols, 10);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'isVentralContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'ventral touch probability', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'numVariables', length(conditions)})

% dorsal touches
subplot(rows, cols, 11);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'isDorsalContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'dorsal touch probability', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colorsPowers,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', addStats, 'numVariables', length(conditions)})


% save
file = fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_bars']);
savefig(gcf, [file fileSuffix]);


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
rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
    'FaceColor', [obsOnColor obsOnAlpha], 'EdgeColor', 'none');
line([0 0], yLims, 'linewidth', 2, 'color', obsOnColor)
line([optoOnPos optoOnPos], yLims, 'linewidth', 2, 'color', colors(2,:))
    
% plot
plotDvPsth(flat, 'velVsPosition', dvName, ...
    {'showLegend', false, 'errorFcn', @(x) nanstd(x), 'conditionColors', colorsPowers, ...
    'plotMouseAvgs', false})

% pimp fig
set(gca, 'YLim', yLims);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

% save
file = fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_speed']);
saveas(gcf, [file fileSuffix '.fig']);

% INDIVIDUAL MICE PLOTS
mice = unique({flat.mouse});
figure('name', exp, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 350 200*length(mice) 150], 'inverthardcopy', 'off'); hold on

for i = 1:length(mice)
    subplot(1,length(mice),i)
    flatSub = flat(strcmp({flat.mouse}, mice{i}));
    
    plotDvPsth(flatSub, 'velVsPosition', dvName, ...
        {'showLegend', false, 'errorFcn', @(x) nanstd(x), 'conditionColors', colorsPowers, ...
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


%% KINEMATICS

% settings
xLims = [-.05 .05];
yLims = [0 .016];


% initializations
colNames = {'contra', 'ipsi'};
rowNames = {'muscimol', 'lesion'};
conditions = {[1,2], [1,2], [3,4,5], [3,4,5]};
if matchPropensities; conditions = {[1,2], [1,2], [3,4], [3,4]}; end
isContraFirst = [true false true false];

flat = flattenData(data, {'mouse', 'session', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'isContraFirst', 'isFore', 'isLeading', 'stepOverKinInterp', 'paw', 'isOptoOn', 'preObsKin'});
flat = flat(~[flat.isWheelBreak] & [flat.isLeading] & [flat.isFore]); % add conditionals here
kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);



figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 700 1600 250])

for i = 1:4
    subplot(2,2,i)
    bins = ismember([flat.conditionNew], conditions{i}) & ...
           [flat.isContraFirst]==isContraFirst(i);
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], [flat(bins).conditionNew]-min(conditions{i})+1 , ...
        {'colors', colors(conditions{i},:), 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
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

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'motorCortex', 'kinematics');
saveas(gcf, [file fileSuffix], 'svg');









