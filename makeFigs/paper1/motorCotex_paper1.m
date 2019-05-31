%% INITIALIZATIONS


% global settings
varsToAvg = {'mouse'};
earlySessions = 1:3;
lateSessions = 5:7;
colors = [217, 65, 244; 244, 149, 66] / 255; % mus, lesion
darkening = .25; % how much to darken control condition relative to manipulated condition
matchPropensities = false;
manipPercent = 25; % take manipPercent percent of best matched manip trials



% initializations
vars.conditionMus = struct('name', 'condition', 'levels', {{'saline', 'muscimol'}}, 'levelNames', {{'saline', 'muscimol'}});
vars.conditionLes = struct('name', 'condition', 'levels', {{'pre', 'post'}}, 'levelNames', {{'pre', 'post'}});
vars.condition = struct('name', 'conditionNew', 'levels', [1 2 3 4 5], 'levelNames', {{'sal', 'mus', 'pre', 'postE', 'postL'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.conditionNum = struct('name', 'conditionNum', 'levels', 1:8);
fileSuffix = '';
if matchPropensities
    fileSuffix = '_matched'; 
    vars.condition = struct('name', 'conditionNew', 'levels', [1 2 3 4], 'levelNames', {{'sal', 'mus', 'pre', 'postE'}});
end
conditionNames = vars.condition.levelNames;

conditionals.isEarly = struct('name', 'conditionNum', 'condition', @(x) ismember(x,earlySessions));
conditionals.isLate = struct('name', 'conditionNum', 'condition', @(x) ismember(x,lateSessions));
conditionals.isPre = struct('name', 'condition', 'condition', @(x) strcmp(x, 'pre'));
conditionals.isPost = struct('name', 'condition', 'condition', @(x) strcmp(x, 'post'));
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isHind = struct('name', 'isFore', 'condition', @(x) x==0);
conditionals.none = struct('name', '', 'condition', @(x) x);

colors = [repmat(colors(1,:),2,1); repmat(colors(2,:),3,1)] .* ...
         [linspace(darkening,1,2), linspace(darkening,1,3)]';
     
% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'mtc_muscimol_data.mat'), 'data');
dataMus = data;
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'mtc_lesion_data.mat'), 'data');
dataLes = data;
clear data
disp('data loaded!')

% put together in one giant data structure
for i = 1:length(dataMus.data)
    
    data(i,1).mouse = dataMus.data(i).mouse;
    data(i,1).sessions = [dataMus.data(i).sessions; dataLes.data(i).sessions];
    
    bins = strcmp({data(i).sessions.condition}, 'saline');
    for j = 1:length(data(i).sessions)
        
        % get condition
        if strcmp({data(i).sessions(j).condition}, 'saline'); c = 1;
        elseif strcmp({data(i).sessions(j).condition}, 'muscimol'); c = 2;
        elseif strcmp({data(i).sessions(j).condition}, 'pre'); c = 3;
        elseif strcmp({data(i).sessions(j).condition}, 'post') && ismember(data(i).sessions(j).conditionNum, earlySessions); c = 4;
        elseif strcmp({data(i).sessions(j).condition}, 'post') && ismember(data(i).sessions(j).conditionNum, lateSessions); c = 5;
        else; c = nan; end
        data(i).sessions(j).conditionNew = c;
    end
end




% propensity score matching
if matchPropensities
    varsToMatch = {'velAtWiskContact', 'angleAtWiskContactContra'};
    experiments = {[1,2], [3,4]};

    mice = {data.mouse};
    flat = struct2table(flattenData(data, [{'mouse', 'session', 'trial', 'conditionNew', 'isLightOn'} varsToMatch]));
    varBins = ismember(flat.Properties.VariableNames, varsToMatch);
    metaBins = ismember(flat.Properties.VariableNames, {'mouse', 'session', 'trial'});

    % find matched trials
    matchedTrials = cell2table(cell(0,3), 'VariableNames', {'mouse', 'session', 'trial'});
    for mouse = mice
        for exp = experiments
            for light = [false true]
                bins = strcmp(flat.mouse, mouse{1}) & ...
                       ismember(flat.conditionNew, exp{1}) & ...
                       flat.isLightOn==light;
                flatSub = flat(bins,:);
                X = table2array(flatSub(:, varBins));
                y = flatSub.conditionNew==exp{1}(2); % is trial in the manip condition
                matchedPairs = propensityMatching(X, y, ...
                    {'percentileThresh', manipPercent, 'predictorNames', varsToMatch, 'verbose', false});
                matchedTrials = [matchedTrials; flatSub(matchedPairs(:), metaBins)];
            end
        end
    end

    % get rid ofnon-matched trials!
    for i = 1:length(data)
        for j = 1:length(data(i).sessions)
            bins = strcmp(matchedTrials.mouse, data(i).mouse) & ...
                   strcmp(matchedTrials.session, data(i).sessions(j).session);
            sesTrials = matchedTrials.trial(bins);
            data(i).sessions(j).trials = data(i).sessions(j).trials(sesTrials);
        end

        % remove unused sessions
        isSessionUsed = ismember({data(i).sessions.session}, unique(matchedTrials.session));
        data(i).sessions = data(i).sessions(isSessionUsed);
    end
end



%% ----------
% PLOT THINGS
%  ----------


%% BARS


% initializations
rows = 4;
cols = 3;
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1600 225*rows])

% success
subplot(rows, cols, 1);
conditions = [vars.isLightOn; vars.condition];
dv = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'success rate', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [0 1], 'ytick', 0:.5:1, ...
    'compareToFirstOnly', true})

% vel
subplot(rows, cols, 2);
conditions = [vars.isLightOn; vars.condition];
dv = getDvMatrix(data, 'trialVel', conditions, varsToAvg);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'velocity (m/s)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), 'ylim', [0 .6]...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% angle
subplot(rows, cols, 3);
conditions = [vars.isLightOn; vars.condition];
dv = getDvMatrix(data, 'trialAngleContra', conditions, varsToAvg);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', ['body angle contra (' char(176) ')'], ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), 'ylim', [-15 15], ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% tail height
subplot(rows, cols, 4);
conditions = [vars.isLightOn; vars.condition];
dv = getDvMatrix(data, 'tailHgt', conditions, varsToAvg);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'tail height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), 'ylim', [], ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% pre obs height
subplot(rows, cols, 5);
conditions = [vars.isFore; vars.isContra; vars.condition];
dv = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'pre paw height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% max obs height
subplot(rows, cols, 6);
conditions = [vars.isFore; vars.isContra; vars.condition];
dv = getDvMatrix(data, 'stepOverMaxHgt', conditions, varsToAvg) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'max paw height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% control step height
subplot(rows, cols, 7);
conditions = [vars.isFore; vars.isContra; vars.condition];
dv = getDvMatrix(data, 'controlStepHgt', conditions, varsToAvg) * 1000;
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'control step height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% height shaping (correlations)
subplot(rows, cols, 8);
conditions = [vars.isLightOn; vars.condition];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, [conditionals.isFore; conditionals.isLeading], 'corr');
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'paw shaping (correlation)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,8,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})


% compute decision determinism and threshold
flat = flattenData(data, {'mouse', 'modPawPredictedDistanceToObs', 'conditionNew', 'isBigStep', 'isLightOn'});
mice = unique({flat.mouse});
lightConditions = [false, true];
[accuracies, thresholds] = deal(nan(2, length(conditionNames), length(mice)));

for i = 1:2
    for m = 1:length(conditionNames)
        for j = 1:length(mice)
            bins = [flat.conditionNew]==m & ...
                    [flat.isLightOn]==lightConditions(i) & ...
                    strcmp({flat.mouse}, mice{j});
            glm = fitglm([flat(bins).modPawPredictedDistanceToObs]', [flat(bins).isBigStep]', 'Distribution', 'binomial');
            accuracies(i,m,j) = mean(round(predict(glm, [flat(bins).modPawPredictedDistanceToObs]'))==[flat(bins).isBigStep]');
            coeffs = glm.Coefficients.Estimate;
            thresholds(i,m,j) = (.5-coeffs(1)) / coeffs(2); % solve for prediction = .5
        end
    end
end

% decision determinism (glm accuracy)
subplot(rows, cols, 9);
barPlotRick(accuracies, {'conditionNames', {{'light off', 'light on'}, conditionNames}, 'ylabel', 'glm accuracy', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [.5 1]})

% decision threshold
subplot(rows, cols, 10);
barPlotRick(thresholds*1000, {'conditionNames', {{'light off', 'light on'}, conditionNames}, 'ylabel', 'big step threshold (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,2,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% ventral touches
subplot(rows, cols, 11);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'isVentralContact', conditions, varsToAvg);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'ventral touch probability', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})

% dorsal touches
subplot(rows, cols, 12);
conditions = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'isDorsalContact', conditions, varsToAvg);
barPlotRick(dv, {'conditionNames', {conditions.levelNames}, 'ylabel', 'dorsal touch probability', ...
    'showViolins', false, 'lineThickness', 2, 'addBars', true, 'conditionColors', repmat(colors,4,1), ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false})


% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'motorCortex', 'bars');
saveas(gcf, [file fileSuffix], 'svg');

%% VARS OVER TIME

% settings
rows = 3;
cols = 3;
sessionsPostLesion = -2:8;


% initializations
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1200 200*rows])
vars.sessionsPostLesion = struct('name', 'sessionsPostLesion', 'levels', sessionsPostLesion);

for i = 1:length(dataMus) % compute session number relative to first lesion session
    firstLesSession = find(strcmp({dataLes(i).data.sessions.condition}, 'post'), 1, 'first');
    sessionsPostLesion = num2cell([dataLes(i).data.sessions.sessionNum] - firstLesSession);
    [dataLes(i).data.sessions.sessionsPostLesion] = sessionsPostLesion{:};
end


% success

% success
subplot(rows, cols, 1);
dv = getDvMatrix(dataLes, 'isTrialSuccess', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'success rate'});

% velocity
subplot(rows, cols, 2);
dv = getDvMatrix(dataLes, 'trialVel', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'velocity (m/s)'});

% body angle
subplot(rows, cols, 3);
dv = getDvMatrix(dataLes, 'trialAngleContra', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', ['body angle contra (' char(176) ')']});

% tail height
subplot(rows, cols, 4);
dv = getDvMatrix(dataLes, 'tailHgt', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'tail height'});

% pre obs height
subplot(rows, cols, 5);
dv = getDvMatrix(dataLes, 'preObsHgt', vars.sessionsPostLesion, {'mouse'}, [conditionals.isFore; conditionals.isLeading]);
sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'pre obs height'});

% pre obs height
subplot(rows, cols, 6);
dv = getDvMatrix(dataLes, 'stepOverMaxHgt', vars.sessionsPostLesion, {'mouse'}, [conditionals.isFore; conditionals.isLeading]);
sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'max paw height'});

% control step height
subplot(rows, cols, 7);
dv = getDvMatrix(dataLes, 'controlStepHgt', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'control step height'});

% height shaping (correlations)
subplot(rows, cols, 8);
dv = getSlopeMatrix(dataLes, {'obsHgt', 'preObsHgt'}, vars.sessionsPostLesion, ...
    {'mouse'}, {'session'}, [conditionals.isFore; conditionals.isLeading], 'corr');
sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'height shaping (correlation)'});

% add lines at time of lesion
for i = 1:8; subplot(rows,cols,i); line([-.5 -.5], get(gca,'YLim'), 'color', 'red'); end

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'motorCortex', 'sessionsOverTime');
saveas(gcf, [file fileSuffix], 'svg');

%% KINEMATICS

% settings
xLims = [-.05 .05];
yLims = [0 .016];


% initializations
colNames = {'contra', 'ipsi'};
rowNames = {'muscimol', 'lesion'};
conditions = {[1,2], [1,2], [3,4,5], [3,4,5]};
isContraFirst = [true false true false];

flat = flattenData(data, {'mouse', 'session', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'isContraFirst', 'isFore', 'isLeading', 'stepOverKinInterp', 'paw', 'conditionNew', 'preObsKin'});
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


%%











