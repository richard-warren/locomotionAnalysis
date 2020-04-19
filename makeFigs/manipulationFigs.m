%% compute experiment data from scratch (only need to do once)

% settings
dataset = 'senLesion';  % senLesion, mtc_lesion, mtc_muscimol

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
clear all; close all  % best to clear workspace before loading these super large datasets

% settings
dataset = 'senLesion';
poolSenLesionConditions = true;  % whether to use all conditions or pool postBi and postContra
splitEarlyLate = true;  % whether to split early and late post-lesion sessions
earlySessions = [1 1];  % min and max sessions to include in 'early' lesion sessions
lateSessions = [4 4];  % min and max sessions to include in 'late' lesion sessions
preSessions = 2;  % only include the most recent 'preSessions' in the 'pre' condition

matchTrials = false;  % whether to use propensity score matching to control for baseline characteristics of locomotion (varsToMatch)
varsToMatch = {'velAtWiskContact', 'angleAtWiskContactContra', 'tailHgtAtWiskContact'};
manipPercent = 20;  % take manipPercent percent of best matched manip trials
miceToExclude = {'sen11'};




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
%         colors = winter(6);
        colors = [ctlStepColor; ctlStepColor; repmat(lesionColor,4,1).*linspace(1,.25,4)'];
        vars.condition = struct('name', 'condition', ...
            'levels', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}}, ...
            'levelNames', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}});
        matchConditions = {'pre' 'postContra'};  % for propensity score matching, math these two conditions
        figConditionals = struct('name', 'conditionNum', 'condition', @(x) x>=earlySessions(1) & x<=earlySessions(2));
    end
    
    extraVars = {'condition'};  % vars to extract using flattenData that are specific to this experiment
    
end


% load data
fprintf(['loading ' dataset ' data... ']);
% tic; load(fullfile(getenv('OBSDATADIR'), 'matlabData', [dataset '_data.mat']), 'data'); toc
tic; load(fullfile('C:\Users\richa\Desktop\', 'matlabData', [dataset '_data.mat']), 'data'); toc  % when working on home computer
data.data = data.data(~ismember({data.data.mouse}, miceToExclude));
mice = {data.data.mouse};


%% pool postContra and postBi for senLesion experiments
if strcmp(dataset, 'senLesion') && poolSenLesionConditions
    fprintf('pooling postContra and postBi conditions... ');
    for i = 1:length(data.data)  % loop across mice    
        for j = 1:length(data.data(i).sessions)
            condition = data.data(i).sessions(j).condition;
            if strcmp(condition, 'postContra') || strcmp(condition, 'postBi')
                data.data(i).sessions(j).condition = 'post';
            elseif strcmp(condition, 'postIpsi') % !!! should i be putting together postIpsi data and pre data?
                data.data(i).sessions(j).condition = 'pre';
            end
        end
    end
end

% remove excess pre sessions
if contains(dataset, 'esion')  % only do this for the lesion data
    fprintf('restricting to %i pre sessions... ', preSessions);
    for i = 1:length(data.data)  % loop across mice    
        firstLesSession = find(contains({data.data(i).sessions.condition}, 'post'), 1, 'first');
        firstPreSession = find(strcmp({data.data(i).sessions.condition}, 'pre'), 1, 'first');
        [data.data(i).sessions(firstPreSession:firstLesSession-preSessions-1).condition] = deal('');
    end
end


% split early and last postLesion lessions
if splitEarlyLate && ~strcmp(dataset, 'mtc_muscimol')  % only do this for lesion experiments
    fprintf('splitting early and late post lesion epochs... ');
    
    for i = 1:length(data.data)
        postSessions = find(strcmp({data.data(i).sessions.condition}, 'post'));
        
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

        bins = strcmp(flat.mouse, mouse{1}) & ...
            ismember(flat.condition, matchConditions);
        if contains(dataset, {'lesion', 'Lesion'})
            bins = bins & [flat.conditionNum]>=earlySessions(1) & [flat.conditionNum]<=earlySessions(2);
            fprintf('restricting matched sessions to ''early'' data...');
        end
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
            if ~isempty(sesTrials)
                data_matched.data(i).sessions(j).trials = data_matched.data(i).sessions(j).trials(sesTrials);  % !!! this assumes all trials are present in sessions(j).trials
                if ~isequal(sesTrials, [data_matched.data(i).sessions(j).trials.trial]'); disp('WARNING! Trials should be the same in matchedTrials and data_matched'); keyboard; end
            end
        end

        % remove unused sessions
        isSessionUsed = ismember({data_matched.data(i).sessions.session}, unique(matchedTrials.session));
        data_matched.data(i).sessions = data_matched.data(i).sessions(isSessionUsed);
    end
    data_unmatched = data;
    data = data_matched; clear data_matched;
end


% define variables and conditionals to be used in bar plots
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore', 'hind'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isHind = struct('name', 'isFore', 'condition', @(x) x==0);
conditionals.isContra = struct('name', 'isContra', 'condition', @(x) x==1);

colorsRaw = colors;  % keep original colors before applying following transformation
if strcmp(dataset, 'mtc_muscimol') && matchTrials; colors(2,:) = mean(colors,1); end  % if matching trials, split difference btwn control and manipulated colors

fprintf('data loaded!\n');


%% plot propensity score matching

% settings
binNum = 100;  % resolution of x axis
percentileLims = [1 99];  % x axis limits
mouseAlpha = .25;

flat = flattenData(data_unmatched, [{'mouse', 'condition'} varsToMatch]);
flat_matched = flattenData(data, [{'mouse', 'condition'} varsToMatch]);

colors_unmatched = colorsRaw;
colors_matched = [colorsRaw(1,:); mean(colorsRaw,1)];

% pdfs
close all
figure('color', 'white', 'menubar', 'none', 'position', [1997.00 595.00 600.00 345.00]);

for i = 1:length(varsToMatch)
    xLims = prctile([flat.(varsToMatch{i})], percentileLims);
    xLims = xLims + [-1 1]*diff(xLims)*.15;
    x = linspace(xLims(1), xLims(2), binNum);
    conditionMeans = nan(length(mice), 2);  % mouse means for two conditions
    
    [pdfs, pdfs_matched] = deal(nan(length(mice), binNum));
    for k = 1:length(matchConditions)
        for j = 1:length(mice)
            
            % unmatched
            subplot(length(varsToMatch),2,(i-1)*2+1); hold on
            bins = strcmp({flat.mouse}, mice{j}) & strcmp({flat.condition}, matchConditions{k});
            pdfs(j,:) = ksdensity([flat(bins).(varsToMatch{i})], x);
            plot(x, pdfs(j,:), 'Color', [colors_unmatched(k,:) mouseAlpha])
            
            % matched
            subplot(length(varsToMatch),2,(i-1)*2+2); hold on
            bins = strcmp({flat_matched.mouse}, mice{j}) & strcmp({flat_matched.condition}, matchConditions{k});
            pdfs_matched(j,:) = ksdensity([flat_matched(bins).(varsToMatch{i})], x);
            plot(x, pdfs_matched(j,:), 'Color', [colors_matched(k,:) mouseAlpha])
        end
        
        % unmatched means
        subplot(length(varsToMatch),2,(i-1)*2+1); hold on
        if i==1; title('all trials'); end
        plot(x, mean(pdfs,1), 'Color', colors_unmatched(k,:), 'LineWidth', 3)
        fill([x(1) x x(end)], [0 mean(pdfs,1) 0], colors_unmatched(k,:), 'FaceAlpha', .2, 'EdgeColor', 'none')
        set(gca, 'XLim', xLims, 'YTick', [], 'YColor', 'none')
        xlabel(varsToMatch{i})
        
        % matched means
        subplot(length(varsToMatch),2,(i-1)*2+2); hold on
        if i==1; title('matched trials'); end
        plot(x, mean(pdfs_matched,1), 'Color', colors_matched(k,:), 'LineWidth', 3)
        fill([x(1) x x(end)], [0 mean(pdfs_matched,1) 0], colors_unmatched(k,:), 'FaceAlpha', .2, 'EdgeColor', 'none')
        set(gca, 'XLim', xLims, 'YTick', [], 'YColor', 'none')
    end
end

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_matchedHistos' suffix1 suffix2]), 'svg');
fprintf('%.3f of trials used in matched sub-population\n', size(flat_matched) / size(flat))


%% bars

% success rate
if strcmp(dataset, 'senLesion'); props = {'YLim', [.2 1], 'YTick', [.2 .6 1]}; else; props = {}; end
figure('position', [200 472.00 300 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'isTrialSuccess', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4], 'test', 'ttest', props{:})
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_success' suffix1 suffix2]), 'svg');

%% velocity
figure('position', [200 472.00 300 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'trialVel', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'velocity (m/s)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_velocity' suffix1 suffix2]), 'svg');

% tail height
figure('position', [400 472.00 300 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'tailHgt', vars.condition, {'mouse'}, [figConditionals])*1000;
barFancy(dv, 'ylabel', 'tail height (mm)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_tailHeight' suffix1 suffix2]), 'svg');

% body angle
figure('position', [600 472.00 300 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'trialAngleContra', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'body angle ({\circ})', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_bodyAngle' suffix1 suffix2]), 'svg');

%% paw height
figure('position', [200 100 300 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isContra; vars.isFore; vars.condition];
if strcmp(dataset, 'senLesion'); figVars = [vars.isFore; vars.condition]; end
dv = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, [figConditionals]) * 1000;
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'paw height (mm)', 'levelNames', {figVars(1:end).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'constantEdgeColor', [.15 .15 .15], ...
    'comparisons', [1 2; 1 3; 1 4; 5 6; 5 7; 5 8], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_pawHeight' suffix1 suffix2]), 'svg');

%% paw correlations
figure('position', [200 100 300 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.condition];
tempConditionals = [figConditionals; conditionals.isLeading; conditionals.isFore];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, figVars, {'mouse'}, {'session'}, tempConditionals, 'corr'); % perform regression for each session, then average slopes across sessions for each mouse
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'paw:obstacle correlation', 'levelNames', {figVars.levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'textRotation', 0, ...
    'comparisons', [1 2; 1 3; 1 4], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_correlations' suffix1 suffix2]), 'svg');

%% paw contacts
figure('position', [2600.00 100 700 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isContra; vars.isFore; vars.isLeading; vars.condition];
if strcmp(dataset, 'senLesion'); figVars = [vars.isFore; vars.condition]; end
dv = 1-getDvMatrix(data, 'anyTouchFrames', figVars, {'mouse'}, [figConditionals]);
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:})
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_pawContacts' suffix1 suffix2]), 'svg');

% baseline step heights
figure('position', [2800.00 100 700 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isFore; vars.condition];
dv = getDvMatrix(data, 'controlStepHgt', figVars, {'mouse'}, [figConditionals])*1000;
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'control step height (mm)', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:})
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_baselineHeights' suffix1 suffix2]), 'svg');

%% forepaw height
figure('position', [200 100 300 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isContra; vars.condition];
if strcmp(dataset, 'senLesion'); figVars = [vars.condition]; end
dv = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, [figConditionals; conditionals.isFore]) * 1000;
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'paw height (mm)', 'levelNames', {figVars(1:end).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_forepawHeight' suffix1 suffix2]), 'svg');

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

sessionsToShow = -(preSessions-1):lateSessions(end);
manipInd = find(sessionsToShow==0);

    
% add sessionsPostLesion to data structure
vars.sessionsPostLesion = struct('name', 'sessionsPostLesion', 'levels', sessionsToShow);
for i = 1:length(data.data) % compute session number relative to first lesion session
    sessionNums = [data.data(i).sessions.sessionNum];
    firstLesSession = sessionNums(find(ismember({data.data(i).sessions.condition}, {'post', 'postEarly', 'postLate'}), 1, 'first'));  % we don't want 'postIpsi' getting in here!!!
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


%% success
figure('position', [200 400 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'isTrialSuccess', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'success rate', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
set(gca, 'YLim', [0 1], 'YTick', 0:.5:1, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_succesOverSessions' suffix1 suffix2]), 'svg');

%% body angle
figure('position', [2018.00 200 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'trialAngleContra', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'body angle', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
set(gca, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_angleOverSessions' suffix1 suffix2]), 'svg');

%% tail height
figure('position', [2018.00 300 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'tailHgt', vars.sessionsPostLesion, {'mouse'});
% if strcmp(dataset, 'senLesion'); dv = dv(:,[1,3:end]); end  % !!! this is a hack to remove the mouse who only has one post lesion sessions
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'tail height', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
set(gca, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_tailHgtOverSessions' suffix1 suffix2]), 'svg');

%% velocity
figure('position', [200 400 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'velAtWiskContact', vars.sessionsPostLesion, {'mouse'});
% if strcmp(dataset, 'senLesion'); dv = dv(:,[1,3:end]); end  % !!! this is a hack to remove the mouse who only has one post lesion sessions
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'velocity (m/s)', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
if strcmp(dataset, 'senLesion'); set(gca, 'YLim', [0 .8]); end
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_velOverSessions' suffix1 suffix2]), 'svg');

%% forepaw height
figure('position', [200 500 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'preObsHgt', vars.sessionsPostLesion, {'mouse'}, conditionals.isFore) * 1000;
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'forepaw height (mm)', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
if strcmp(dataset, 'senLesion'); set(gca, 'YLim', [6 10.5]); end
set(gca, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_forepawHgtOverSessions' suffix1 suffix2]), 'svg');

%% hindpaw height
figure('position', [2018.00 500 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'preObsHgt', vars.sessionsPostLesion, {'mouse'}, conditionals.isHind) * 1000;
% if strcmp(dataset, 'senLesion'); dv = dv(:,[1,3:end]); end  % !!! this is a hack to remove the mouse who only has one post lesion sessions
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'hindpaw height (mm)', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
set(gca, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_hindpawHgtOverSessions' suffix1 suffix2]), 'svg');

%% correlations
figure('position', [200 600 521.00 239.00], 'color', 'white', 'menubar', 'none')
tempConditionals = [figConditionals; conditionals.isLeading; conditionals.isFore];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.sessionsPostLesion, {'mouse'}, {'session'}, tempConditionals, 'corr'); % perform regression for each session, then average slopes across sessions for each mouse
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'correlation', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
yLims = get(gca, 'ylim'); ln = line([.5 .5],yLims, 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom'); set(gca, 'ylim', yLims);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_corrOverSessions' suffix1 suffix2]), 'svg');

%% success over days post wisk trim

sessionsToShow = -(preSessions+1):3;

% add sessionsPostTim to data structure
vars.sessionsPostTrim = struct('name', 'sessionsPostTrim', 'levels', sessionsToShow);
for i = 1:length(data.data) % compute session number relative to first whisker trim session
    sessionNums = [data.data(i).sessions.sessionNum];
    firstTrimSession = sessionNums(find(strcmp({data.data(i).sessions.condition}, 'noWisk'), 1, 'first'));  % we don't want 'postIpsi' getting in here!!!
    sessionsPostTrim = num2cell([data.data(i).sessions.sessionNum] - firstTrimSession+1);
    [data.data(i).sessions.sessionsPostTrim] = sessionsPostTrim{:};

    % make sure 'pre' not included
    for j = 1:length(data.data(i).sessions)
        if ismember(data.data(i).sessions(j).condition, {'preTrim', 'pre'})
            data.data(i).sessions(j).sessionsPostTrim = nan;
        end
    end
end

% success post trim
figure('position', [2018.00 100 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'isTrialSuccess', vars.sessionsPostTrim, {'mouse'});
if strcmp(dataset, 'senLesion'); dv = dv(:,[1,3:end]); end  % !!! this is a hack to remove the mouse who only has one post lesion sessions
sesPlotRick(dv', 'xvals', vars.sessionsPostTrim.levels, 'ylabel', 'success rate', 'xlabel', 'sessions post whisker trim', ...
    'meanColor', axisColor, 'colors', 'lines', 'alpha', .5);
set(gca, 'YLim', [0 1], 'YTick', 0:.5:1, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', lesionColor, 'linewidth', 3); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_succesOverSessionsTrim' suffix1 suffix2]), 'svg');


%% kinematics of leading forepaw

% settings
xLims = [-.05 0];
earlyLate = 'early';  % 'early' or 'late'

flat = flattenData(data, {'mouse', 'session', 'trial', ...
    'preObsKin', 'conditionNum', 'condition', 'obsHgt', 'isLeading', 'isFore', 'isContra'});
if contains(dataset, {'lesion', 'Lesion'})
    if strcmp(earlyLate, 'early')
        flat = flat([flat.conditionNum]>=earlySessions(1) & [flat.conditionNum]<=earlySessions(2));  % for early sessions
        colors_temp = colors;
    elseif strcmp(earlyLate, 'late')
        postBins = (strcmp({flat.condition}, 'post') & [flat.conditionNum]>=lateSessions(1) & [flat.conditionNum]<=lateSessions(2));
        flat = flat(postBins | strcmp({flat.condition}, 'pre'));  % for late sessions
        colors_temp = [ctlStepColor; mean([ctlStepColor; lesionColor],1)];
    end
end
flat = flat([flat.isLeading] & [flat.isFore]);  % add conditionals as required
[~,condition] = ismember({flat.condition}, matchConditions);
kinData = permute(cat(3, flat.preObsKin), [3,1,2]);

figure('position', [2018.00 700 521.00 239.00], 'color', 'white', 'menubar', 'none')
plotKinematics(kinData(:,[1,3],:), [flat.obsHgt], condition, ...
        'colors', colors_temp, 'obsAlpha', 1, 'lineAlpha', .8, ...
        'obsColors', repelem(obsColor,2,1), 'mouseNames', {flat.mouse}, ...
        'lineWidth', 3, 'errorFcn', @std) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims)

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_kinematics' suffix1 suffix2 '_' earlyLate]), 'svg');


%% decision making
flat = flattenData(data, ...
    [m.predictorsAll, {'mouse', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', 'modPawOnlySwing', 'isTrialSuccess', 'condition', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'modPawKinInterp', 'preModPawKinInterp', 'firstModPaw', 'preModPawDeltaLength'}]);

%% heatmaps
plotDecisionHeatmaps(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'normalize', 'col', 'xLims', [-20 15], ...
    'deltaMin', 0, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', false, 'plotMice', false, 'colors', colors, 'outcome', 'isModPawLengthened', ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_heatMaps' suffix1 suffix2]));

%% entropy
heatmaps = plotDecisionHeatmaps(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'normalize', 'col', 'xLims', [-20 15], ...
    'deltaMin', 0, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', true, 'plotMice', false, 'colors', colors, 'outcome', 'isModPawLengthened');
ent = cellfun(@entropy, heatmaps);

figure('position', [2000 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
barFancy(ent, 'levelNames', {vars.condition.levelNames}, 'ylabel', 'landing position entropy', ...
    'colors', sensColors, barProperties{:}, 'showBars', false, 'lineThickness', 4)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceEntropy'), 'svg');

%% trials scatters
plotDecisionTrials(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'outcome', 'isBigStep', ...
    'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'rowColors', colors, 'xLims', [-.11 .06], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_decisionKin' suffix1 suffix2]));

%% model accuracies
accuracies = plotModelAccuracies(flat, m.predictors, 'isModPawLengthened', 'modelTransfers', [1 2], ...
    'weightClasses', true, 'condition', 'condition', 'levels', vars.condition.levels, 'kFolds', 5, ...
    'deltaMin', .5, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', colors, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_models' suffix1 suffix2]));

%% decision threshold
plotDecisionThresholds(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'outcome', 'isModPawLengthened', ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', colors, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_thresholds' suffix1 suffix2]));

%% model predictors
plotPredictors(flat, m.predictors, 'isModPawLengthenedshj', 'colors', colors, 'avgMice', true, ...
    'condition', 'condition', 'levels', vars.condition.levels, ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', true, 'lightOffOnly', false);

%% BASELINE LOCOMOTION SCATTERS

close all
vars.conditionSub = struct('name', 'condition', 'levels', {matchConditions}, 'levelNames', {matchConditions});  % only first two conditions
figure('position', [200 525 757 275], 'color', 'white', 'menubar', 'none');

% velocity
subplot(1,3,1)
dv = getDvMatrix(data, 'trialVel', vars.conditionSub, {'mouse'}, [figConditionals]);
pairedScatter(dv, 'xlabel', matchConditions{1}, 'ylabel', matchConditions{2}, 'title', 'velocity (m/s)', ...
    'colors', colors(1:2,:))

% tail height
subplot(1,3,2)
dv = getDvMatrix(data, 'tailHgt', vars.conditionSub, {'mouse'}, [figConditionals])*1000;
pairedScatter(dv, 'xlabel', matchConditions{1}, 'ylabel', matchConditions{2}, 'title', 'tail height (mm)', ...
    'colors', colors(1:2,:))

% body angle
subplot(1,3,3)
dv = getDvMatrix(data, 'trialAngleContra', vars.conditionSub, {'mouse'}, [figConditionals]);
pairedScatter(dv, 'xlabel', matchConditions{1}, 'ylabel', matchConditions{2}, 'title', 'body angle (degrees)', ...
    'colors', colors(1:2,:))
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_locoScatters' suffix1 suffix2]), 'svg');

% %% body angle (temp, should give same answer as body angle above!)
% figure('position', [2018.00 200 521.00 239.00], 'color', 'white', 'menubar', 'none')
% dv = getDvMatrix(data, 'trialAngleContra', vars.sessionsPostLesion, {'mouse'});
% sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'body angle', 'xlabel', 'sessions post lesion', ...
%     'compareTo', 1:manipInd, sessionPlotProperties{:});
% set(gca, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
% ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
% saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_angleOverSessions' suffix1 suffix2]), 'svg');
