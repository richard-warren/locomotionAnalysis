%% compute experiment data from scratch (only need to do once)

% settings
dataset = 'mtc_muscimol';  % senLesion, mtc_lesion, mtc_muscimol

if strcmp(dataset,'senLesion'); sheet='senLesionNotes'; elseif strcmp(dataset,'mtc_lesion'); sheet='mtcLesionNotes'; elseif strcmp(dataset,'mtc_muscimol'); sheet='muscimolNotes'; end
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', sheet);
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', [dataset '_data.mat']), 'data', '-v7.3'); disp('data saved!')
clear all


%% initializations
clear all;  % best to clear workspace before loading these super large datasets

% settings
dataset = 'mtc_lesion';
poolSenLesionConditions = true;  % whether to use all conditions or pool postBi and postContra
splitEarlyLate = true;  % whether to split early and late post-lesion sessions
earlySessions = [1 3];  % min and max sessions to include in 'early' lesion sessions
lateSessions = [5 7];  % min and max sessions to include in 'late' lesion sessions

matchTrials = false;  % whether to use propensity score matching to control for baseline characteristics of locomotion (varsToMatch)
varsToMatch = {'velAtWiskContact', 'angleAtWiskContactContra', 'tailHgtAtWiskContact'};
manipPercent = 25;  % take manipPercent percent of best matched manip trials



% initializations
global_config;
if matchTrials; suffix1='_matched'; else; suffix1=''; end
if ~poolSenLesionConditions && strcmp(dataset,'senLesion'); suffix2='_allConditions'; else; suffix2=''; end
if ~poolSenLesionConditions && splitEarlyLate; disp('WARNING! Splitting early and late sessions only makes sense when pooling senLesionConditions'); keyboard; end


% motor cortex muscimol
if strcmp(dataset, 'mtc_muscimol')
    colors = [ctlStepColor; muscimolColor];
    vars.condition = struct('name', 'condition', 'levels', {{'saline', 'muscimol'}}, 'levelNames', {{'saline', 'muscimol'}});
    figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
    matchConditions = {'saline' 'muscimol'};  % for propensity score matching, math these two conditions
    
% motor cortex lesion
elseif strcmp(dataset, 'mtc_lesion')
    if splitEarlyLate
        colors = [ctlStepColor; lesionColor; mean([ctlStepColor; lesionColor],1)];
        vars.condition = struct('name', 'condition', 'levels', {{'pre', 'postEarly', 'postLate'}}, 'levelNames', {{'pre', 'postEarly', 'postLate'}});
        figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
        matchConditions = {'pre' 'postEarly'};  % for propensity score matching, math these two conditions
    else
        colors = [ctlStepColor; lesionColor];
        vars.condition = struct('name', 'condition', 'levels', {{'pre', 'post'}}, 'levelNames', {{'pre', 'post'}});
        figConditionals = struct('name', 'conditionNum', 'condition', @(x) x>=earlySessions(1) & x<=earlySessions(2));
        matchConditions = {'pre' 'post'};  % for propensity score matching, math these two conditions
    end

% barrel cortex lesion
elseif strcmp(dataset, 'senLesion')
    
    if poolSenLesionConditions
        if splitEarlyLate
            colors = [ctlStepColor; lesionColor; mean([ctlStepColor; lesionColor],1); lesionColor*.25];
            vars.condition = struct('name', 'condition', 'levels', {{'pre', 'postEarly', 'postLate', 'noWisk'}}, 'levelNames', {{'pre', 'postEarly', 'postLate', 'no whiskers'}});
            matchConditions = {'pre' 'postEarly'};  % for propensity score matching, math these two conditions
            figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
        else 
            colors = [ctlStepColor; lesionColor; lesionColor*.5];
            vars.condition = struct('name', 'condition', 'levels', {{'pre', 'post', 'noWisk'}}, 'levelNames', {{'pre', 'post', 'no whiskers'}});
            matchConditions = {'pre' 'post'};  % for propensity score matching, math these two conditions
            figConditionals = struct('name', 'conditionNum', 'condition', @(x) x>=earlySessions(1) & x<=earlySessions(2));
        end
    else
        colors = winter(6);
        vars.condition = struct('name', 'condition', ...
            'levels', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}}, ...
            'levelNames', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}});
        matchConditions = {'pre' 'postContra'};  % for propensity score matching, math these two conditions
    end
    
    extraVars = {'condition'};  % vars to extract using flattenData that are specific to this experiment
    
end


% load data
fprintf(['loading ' dataset ' data... ']);
load(fullfile(getenv('OBSDATADIR'), 'matlabData', [dataset '_data.mat']), 'data');
mice = {data.data.mouse};


% pool postContra and postBi for senLesion experiments
if strcmp(dataset, 'senLesion') && poolSenLesionConditions
    fprintf('pooling postContra and postBi conditions... ');
    for i = 1:length(data.data)  % loop across mice
        for j = 1:length(data.data(i).sessions)
            condition = data.data(i).sessions(j).condition;
            if strcmp(condition, 'postContra') || strcmp(condition, 'postBi')
                data.data(i).sessions(j).condition = 'post';
            end
        end
    end
end


% split early and last postLesion lessions
if splitEarlyLate && ~strcmp(dataset, 'mtc_muscimol')  % only do this for lesion experiments
    fprintf('splitting early and late post lesion epochs... ');
    
    for i = 1:length(data.data)
        postSessions = find(contains({data.data(i).sessions.condition}, 'post'));
        
        % early sessions
        inds = earlySessions(1):earlySessions(2);
        inds = inds(ismember(inds, 1:length(postSessions)));  % keep only those inds that can be found in postSessions
        if ~isempty(inds)
            for j = postSessions(inds); data.data(i).sessions(j).condition = 'postEarly'; end
        end
        
        % late sessions
        inds = lateSessions(1):lateSessions(2);
        inds = inds(ismember(inds, 1:length(postSessions)));  % keep only those inds that can be found in postSessions
        if ~isempty(inds)
            for j = postSessions(inds); data.data(i).sessions(j).condition = 'postLate'; end
        end
    end    
end


% propensity score matching
if matchTrials
    fprintf('matching propensities... ');
    
    flat = struct2table(flattenData(data, [{'mouse', 'session', 'trial', 'condition', 'isLightOn', 'conditionNum'} varsToMatch]));
    varBins = ismember(flat.Properties.VariableNames, varsToMatch);
    metaBins = ismember(flat.Properties.VariableNames, {'mouse', 'session', 'trial'});  % metadata associated with each trial

    % find matched trials
    matchedTrials = cell2table(cell(0,3), 'VariableNames', {'mouse', 'session', 'trial'});  % struct containing the mouse, session, and trial number for propensity score matched trials
    for mouse = mice

        bins = strcmp(flat.mouse, mouse{1}) & ismember(flat.condition, matchConditions);
        if contains(dataset, {'lesion', 'Lesion'}); bins = bins & [flat.conditionNum]>=earlySessions(1) & [flat.conditionNum]<=earlySessions(2); end
        flat_sub = flat(bins,:);

        X = table2array(flat_sub(:, varBins));
        y = ismember(flat_sub.condition, matchConditions{2}); % is trial in the manip condition
        matchedPairs = propensityMatching(X, y, ...
            {'percentileThresh', manipPercent, 'predictorNames', varsToMatch, 'verbose', true});
        matchedTrials = [matchedTrials; flat_sub(matchedPairs(:), metaBins)];
    end

    % get rid of non-matched trials
    data_matched = data;
    for i = 1:length(data_matched.data)
        for j = 1:length(data_matched.data(i).sessions)
            bins = strcmp(matchedTrials.mouse, data_matched.data(i).mouse) & ...
                   strcmp(matchedTrials.session, data_matched.data(i).sessions(j).session);
            sesTrials = sort(matchedTrials.trial(bins));
            data_matched.data(i).sessions(j).trials = data_matched.data(i).sessions(j).trials(sesTrials);  % !!! this assumes all trials are present in sessions(j).trials
            if ~isequal(sesTrials, [data_matched.data(i).sessions(j).trials.trial]'); disp('WARNING! Trials should be the same in matchedTrials and data_matched'); keyboard; end
        end

        % remove unused sessions
        isSessionUsed = ismember({data_matched.data(i).sessions.session}, unique(matchedTrials.session));
        data_matched.data(i).sessions = data_matched.data(i).sessions(isSessionUsed);
    end
    data = data_matched; clear data_matched;
end


% define variables and conditionals to be used in bar plots
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore', 'hind'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isContra = struct('name', 'isContra', 'condition', @(x) x==1);

fprintf('data loaded!\n');


%% bars

% success rate
figure('position', [2000.00 472.00 320 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'isTrialSuccess', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
set(gca, 'YTick', 0:.5:1)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_success' suffix1 suffix2]), 'svg');

% velocity
figure('position', [2300.00 472.00 320 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'trialVel', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'velocity (m/s)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
set(gca, 'YTick', 0:.2:.6)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_velocity' suffix1 suffix2]), 'svg');

% tail height
figure('position', [2600.00 472.00 320 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'tailHgt', vars.condition, {'mouse'}, [figConditionals])*1000;
barFancy(dv, 'ylabel', 'tail height (mm)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
set(gca, 'YTick', 0:10:20)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_tailHeight' suffix1 suffix2]), 'svg');

% body angle
figure('position', [2900.00 472.00 320 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'trialAngleContra', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'body angle ({\circ})', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
set(gca, 'YTick', -15:15:15)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_bodyAngle' suffix1 suffix2]), 'svg');

% paw height
figure('position', [2000.00 100 320 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isContra; vars.isFore; vars.condition];
dv = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, [figConditionals]) * 1000;
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'paw height (mm)', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:})
set(gca, 'YTick', 0:7:14)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_pawHeight' suffix1 suffix2]), 'svg');

% paw correlations
figure('position', [2300.00 100 320 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.condition];
tempConditionals = [figConditionals; conditionals.isLeading; conditionals.isFore];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, figVars, {'mouse'}, {'session'}, tempConditionals, 'corr'); % perform regression for each session, then average slopes across sessions for each mouse
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'paw-hurdle correlation', 'levelNames', {figVars.levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'textRotation', 0)
% set(gca, 'YTick', 0:7:14)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_correlations' suffix1 suffix2]), 'svg');

% paw contacts
figure('position', [2600.00 100 700 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isContra; vars.isFore; vars.isLeading; vars.condition];
dv = 1-getDvMatrix(data, 'anyTouchFrames', figVars, {'mouse'}, [figConditionals]);
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:})
% set(gca, 'YTick', 0:7:14)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_pawContacts' suffix1 suffix2]), 'svg');

% baseline step heights
figure('position', [2800.00 100 700 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'controlStepHgt', figVars, {'mouse'}, [figConditionals])*1000;
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'control step height (mm)', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:})
% set(gca, 'YTick', 0:7:14)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_baselineHeights' suffix1 suffix2]), 'svg');

if matchTrials
    figure('name', 'matched variables', 'position', [1997.00 838.00 659.00 195.00], 'color', 'white', 'menubar', 'none');
    for i = 1:length(varsToMatch)
        subplot(1,length(varsToMatch),i)
        dv = getDvMatrix(data, varsToMatch{i}, vars.condition, {'mouse'}, [figConditionals]);
        barFancy(dv, 'ylabel', varsToMatch{i}, 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
    end
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_controlledVars' suffix1 suffix2]), 'svg');
end


%% sessions over time


if strcmp(dataset, 'mtc_lesion')
    sessionsToShow = -1:8;
elseif strcmp(dataset, 'senLesion')
    sessionsToShow = -2:4;
end
manipInd = find(sessionsToShow==0);

    
% add sessionsPostLesion to data structure
vars.sessionsPostLesion = struct('name', 'sessionsPostLesion', 'levels', sessionsToShow);
for i = 1:length(data.data) % compute session number relative to first lesion session
    sessionNums = [data.data(i).sessions.sessionNum];
    firstLesSession = sessionNums(find(contains({data.data(i).sessions.condition}, 'post'), 1, 'first'));
    sessionsPostLesion = num2cell([data.data(i).sessions.sessionNum] - firstLesSession+1);
    [data.data(i).sessions.sessionsPostLesion] = sessionsPostLesion{:};

    if strcmp(dataset, 'senLesion')  % make sure preTrim not included in pre, and postTrim not included in post
        for j = 1:length(data.data(i).sessions)
            if ismember(data.data(i).sessions(j).condition, {'preTrim', 'noWisk'})
                data.data(i).sessions(j).sessionsPostLesion = nan;
            end
        end
    end
end


% success
figure('position', [2018.00 100 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'isTrialSuccess', vars.sessionsPostLesion, {'mouse'});
if strcmp(dataset, 'senLesion'); dv = dv(:,[1,3:end]); end  % !!! this is a hack to remove the mouse who only has one post lesion sessions
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'success rate', 'xlabel', 'sessions post lesion', ...
    'meanColor', axisColor, 'colors', 'lines', 'alpha', .5);
set(gca, 'YLim', [0 1], 'YTick', 0:.5:1, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', lesionColor, 'linewidth', 3); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_succesOverSessions' suffix1 suffix2]), 'svg');

% success stats
fprintf('\nSUCCESS VS. BASELINE\n')
bl = mean(dv(1:manipInd-1,:),1);
for i = 1:sessionsToShow(end)+1
    [sig, p] = ttest(bl, dv(manipInd+i-1,:));
    fprintf('session %i p: %.5f\n', sessionsToShow(manipInd+i-1), p);
end


% velocity
figure('position', [2018.00 400 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'velAtWiskContact', vars.sessionsPostLesion, {'mouse'});
if strcmp(dataset, 'senLesion'); dv = dv(:,[1,3:end]); end  % !!! this is a hack to remove the mouse who only has one post lesion sessions
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'velocity (m/s)', 'xlabel', 'sessions post lesion', ...
    'meanColor', axisColor, 'colors', 'lines', 'alpha', .5);
set(gca, 'YTick', [.2 .35 .5], 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', lesionColor, 'linewidth', 3); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_velOverSessions' suffix1 suffix2]), 'svg');

% stats
fprintf('\nVELOCITY VS. BASELINE\n')
bl = mean(dv(1:manipInd-1,:),1);
for i = 1:sessionsToShow(end)+1
    [sig, p] = ttest(bl, dv(manipInd+i-1,:));
    fprintf('session %i p: %.5f\n', sessionsToShow(manipInd+i-1), p);
end

% correlations
figure('position', [2018.00 700 521.00 239.00], 'color', 'white', 'menubar', 'none')
tempConditionals = [figConditionals; conditionals.isLeading; conditionals.isFore];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.sessionsPostLesion, {'mouse'}, {'session'}, tempConditionals, 'corr'); % perform regression for each session, then average slopes across sessions for each mouse
if strcmp(dataset, 'senLesion'); dv = dv(:,[1,3:end]); end  % !!! this is a hack to remove the mouse who only has one post lesion sessions
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'correlation', 'xlabel', 'sessions post lesion', ...
    'meanColor', axisColor, 'colors', 'lines', 'alpha', .5);
set(gca, 'YTick', [.1 .35 .6], 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', lesionColor, 'linewidth', 3); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_corrOverSessions' suffix1 suffix2]), 'svg');

% stats
fprintf('\nCORRS VS. BASELINE\n')
bl = mean(dv(1:manipInd-1,:),1);
for i = 1:sessionsToShow(end)+1
    [sig, p] = ttest(bl, dv(manipInd+i-1,:));
    fprintf('session %i p: %.5f\n', sessionsToShow(manipInd+i-1), p);
end


%% kinematics of leading forepaw

% settings
xLims = [-.05 0];

flat = flattenData(data, {'mouse', 'session', 'trial', ...
    'preObsKin', 'conditionNum', 'condition', 'obsHgt', 'isLeading', 'isFore', 'isContra'});
if contains(dataset, {'lesion', 'Lesion'})
    flat = flat([flat.conditionNum]>=earlySessions(1) & [flat.conditionNum]<=earlySessions(2));
end
flat = flat([flat.isLeading] & [flat.isFore]);  % add conditionals as required
[~,condition] = ismember({flat.condition}, matchConditions);
kinData = permute(cat(3, flat.preObsKin), [3,1,2]);

figure('position', [2018.00 700 521.00 239.00], 'color', 'white', 'menubar', 'none')
plotKinematics(kinData(:,[1,3],:), [flat.obsHgt], condition, ...
        'colors', colors(1:2,:), 'obsAlpha', 1, 'lineAlpha', .8, ...
        'obsColors', repelem(obsColor,2,1), 'mouseNames', {flat.mouse}, ...
        'lineWidth', 3, 'errorFcn', @std) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims)

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_kinematics' suffix1 suffix2]), 'svg');


%% decision making

% settings
kFolds = 10;
% predictors = {'modPawX', 'obsHgt', 'velAtWiskContact', 'wiskContactPosition'};
modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
lightOffOnly = false;  % whether to restrict to light on trials


flat = flattenData(data, {'mouse', 'session', 'trial', 'conditionNum', 'modPawOnlySwing', 'condition', ...
    'isBigStep', 'obsHgt', 'velAtWiskContact', 'wiskContactPosition', 'modPawKin', 'contactInd', 'isLightOn'});
if contains(dataset, {'lesion', 'Lesion'})
    flat = flat([flat.conditionNum]>=earlySessions(1) & [flat.conditionNum]<=earlySessions(2));
end
if lightOffOnly; flat = flat([flat.isLightOn]~=1); end
if modPawOnlySwing; flat = flat([flat.modPawOnlySwing]==1); end
flat = flat(~isnan([flat.isBigStep]));

% compute x position of first mod paw at moment of whisker contact
for i = 1:length(flat)
    ind = max(flat(i).contactInd,1);
    ind = min(ind, size(flat(i).modPawKin,2));
    flat(i).modPawX = flat(i).modPawKin(1,ind);
end


%% compute models
accuracies = nan(3, length(mice));  % rows are: pre accuracy, post accuracy, pre evaluated on post accuracy
X = [[flat.modPawX]; [flat.obsHgt]; [flat.velAtWiskContact]; [flat.wiskContactPosition]]';
y = [flat.isBigStep];

for i = 1:length(mice)
    
    mouseBins = strcmp({flat.mouse}, mice{i});
    preBins = strcmp({flat.condition}, matchConditions(1));
    postBins = strcmp({flat.condition}, matchConditions(2));
    
    % pre model
    X_sub = X(mouseBins & preBins,:);
    y_sub = y(mouseBins & preBins);
    crossVals = cvpartition(length(y_sub), 'kfold', kFolds);  % corss validation splits
    mouseAccuracies = nan(1, kFolds);
    
    for j = 1:kFolds
        glm = fitglm(X_sub(crossVals.training(j),:), y_sub(crossVals.training(j)), 'Distribution', 'binomial');
        predictions = round(predict(glm, X_sub(crossVals.test(j),:)));
        mouseAccuracies(j) = mean(y_sub(crossVals.test(j))==predictions');
    end
    glmPre = fitglm(X_sub, y_sub, 'Distribution', 'binomial');
    accuracies(1,i) = mean(mouseAccuracies);
    
    
    % post model
    X_sub = X(mouseBins & postBins,:);
    y_sub = y(mouseBins & postBins);
    crossVals = cvpartition(length(y_sub), 'kfold', kFolds);  % corss validation splits
    mouseAccuracies = nan(1, kFolds);
    
    for j = 1:kFolds
        glm = fitglm(X_sub(crossVals.training(j),:), y_sub(crossVals.training(j)), 'Distribution', 'binomial');
        predictions = round(predict(glm, X_sub(crossVals.test(j),:)));
        mouseAccuracies(j) = mean(y_sub(crossVals.test(j))==predictions');
    end
    accuracies(2,i) = mean(mouseAccuracies);
    
    % pre model evaluated on post data
    accuracies(3,i) = mean(round(predict(glmPre, X_sub)) == y_sub'); 
end






