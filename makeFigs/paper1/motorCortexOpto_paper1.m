%% INITIALIZATIONS


% global settings
sessions = {'190814_000', '190814_001', '190814_002'}; % ASSUMES ALL SESSIONS PROVIDED ONLY HAVE DATA FROM ONE BRAIN REGION

varsToAvg = {'mouse'};
colors = hsv(2);
matchPropensities = false;
varsToMatch = {'velAtWiskContact', 'tailHgtAtWiskContact'};
manipPercent = 25; % take manipPercent percent of best matched manip trials


% initializations
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'optoNotes');
sessionInfo = sessionInfo(ismember(sessionInfo.session, sessions),:);

vars.condition = struct('name', 'isOptoOn', 'levels', [0 1], 'levelNames', {{'no opto', 'opto'}});
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
rows = 4;
cols = 3;
figure('name', 'motorCortex', 'color', 'white', 'menubar', 'none', 'position', [1984 393 1019 591], 'Renderer', 'painters')
figConditionals = [conditionals.none];

% success
subplot(rows, cols, 1);
conditions = [vars.condition];
dv = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'success rate', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [0 1], 'ytick', 0:.5:1, ...
    'compareToFirstOnly', true})

% vel
subplot(rows, cols, 2);
conditions = [vars.condition];
dv = getDvMatrix(data, 'velAtWiskContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'velocity (m/s)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), 'ylim', [0 .6]...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% angle
% subplot(rows, cols, 3);
% conditions = [vars.condition];
% dv = getDvMatrix(data, 'angleAtWiskContactContra', conditions, varsToAvg, figConditionals);
% barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', ['body angle contra (' char(176) ')'], ...
%     'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), 'ylim', [-15 15], ...
%     'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% tail height
subplot(rows, cols, 4);
conditions = [vars.condition];
dv = getDvMatrix(data, 'tailHgtAtWiskContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'tail height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), 'ylim', [], ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% pre obs height
subplot(rows, cols, 5);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'pre paw height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% max obs height
subplot(rows, cols, 6);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'stepOverMaxHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'max paw height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% control step height
subplot(rows, cols, 7);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'controlStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'control step height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% height shaping (correlations)
subplot(rows, cols, 8);
conditions = [vars.condition];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, ...
    [conditionals.isFore; conditionals.isLeading; figConditionals], 'corr');
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'paw shaping (correlation)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,8,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})


% compute decision determinism and threshold
flat = flattenData(data, {'mouse', 'modPawPredictedDistanceToObs', 'conditionNew', 'isBigStep', 'isOptoOn'});
mice = unique({flat.mouse});
lightConditions = [false, true];
[accuracies, thresholds, f1s] = deal(nan(length(conditionNames), length(mice)));

for m = 1:length(conditions.levels)
    for j = 1:length(mice)
        bins = [flat.isOptoOn]==conditions.levels(m) & ...
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
subplot(rows, cols, 9);
barPlotRick(f1s, {'conditionNames', {conditionNames}, 'ylabel', 'f1 score', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [.5 1]})

% decision threshold
% !!! how to apply fig conditionals to flattened data?
subplot(rows, cols, 10);
barPlotRick(thresholds*1000, {'conditionNames', {conditionNames}, 'ylabel', 'big step threshold (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% ventral touches
subplot(rows, cols, 11);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'isVentralContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'ventral touch probability', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% dorsal touches
subplot(rows, cols, 3);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'isDorsalContact', conditions, varsToAvg, figConditionals);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'dorsal touch probability', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})


% save
file = fullfile(getenv('OBSDATADIR'), 'figures', 'opto', 'bars');
saveas(gcf, [file fileSuffix], 'svg');


%% KINEMATICS

% settings
xLims = [-.05 .05];
yLims = [0 .016];


% initializations
colNames = {'contra', 'ipsi'};
rowNames = {'muscimol', 'lesion'};
conditions = {[1,2], [1,2], [3,4,5], [3,4,5]};
if matchPropensities; conditions = {[1,2], [1,2], [3,4], [3,4]};; end
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









