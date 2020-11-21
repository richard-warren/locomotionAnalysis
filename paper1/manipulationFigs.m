%% compute experiment data from scratch (only need to do once)

% settings
dataset = 'mtc_lesion';  % senLesion, mtc_lesion, mtc_muscimol

if strcmp(dataset,'senLesion'); sheet='senLesionNotes'; elseif strcmp(dataset,'mtc_lesion'); sheet='mtcLesionNotes'; elseif strcmp(dataset,'mtc_muscimol'); sheet='muscimolNotes'; end
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', sheet);
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('SSD'), 'paper1', [dataset '_data.mat']), 'data', '-v7.3'); disp('data saved!')
clear all


%% initializations
% clear all; close all  % best to clear workspace before loading these super large datasets

% settings
dataset = 'senLesion';  % senLesion, mtc_muscimol, or mtc_lesion
poolSenLesionConditions = false;  % whether to use all conditions or pool postBi and postContra
splitEarlyLate = false;  % whether to split early and late post-lesion sessions
if strcmp(dataset, 'senLesion')
    earlySessions = [1 1];  % min and max sessions to include in 'early' lesion sessions ([1 1] for sen lesions, [1 3] for mtc lesions)
    lateSessions = [4 4];  % min and max sessions to include in 'late' lesion sessions
else
    earlySessions = [1 3];  % min and max sessions to include in 'early' lesion sessions ([1 1] for sen lesions, [1 3] for mtc lesions)
    lateSessions = [5 7];  % min and max sessions to include in 'late' lesion sessions
end
preSessions = 2;  % only include the most recent 'preSessions' in the 'pre' condition

matchTrials = false;  % whether to use propensity score matching to control for baseline characteristics of locomotion (varsToMatch)
varsToMatch = {'velAtWiskContact', 'angleAtWiskContactContra', 'tailHgtAtWiskContact'};
manipPercent = 20;  % take manipPercent percent of best matched manip trials
miceToExclude = {'sen11'};
% miceToExclude = {'sen11', 'sen13', 'sen15', 'sen16'};




% initializations
paper1_config;
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
        figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
        matchConditions = {'pre' 'post'};  % for propensity score matching, math these two conditions
    end

% barrel cortex lesion
elseif strcmp(dataset, 'senLesion')
    
    if poolSenLesionConditions
        if splitEarlyLate
            colors = [ctlStepColor; lesionColor; mean([ctlStepColor; lesionColor],1); lesionColor*.25];
            vars.condition = struct('name', 'condition', 'levels', {{'pre', 'postEarly', 'postLate', 'noWisk'}}, 'levelNames', {{'pre', 'postEarly', 'postLate', 'no whiskers'}});
%             vars.condition = struct('name', 'condition', 'levels', {{'preTrim', 'pre', 'postEarly', 'postLate', 'noWisk'}}, 'levelNames', {{'preTrim', 'pre', 'postEarly', 'postLate', 'no whiskers'}});
            matchConditions = {'pre' 'postEarly'};  % for propensity score matching, math these two conditions
            figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
        else 
            colors = [ctlStepColor; lesionColor; lesionColor*.5];
            vars.condition = struct('name', 'condition', 'levels', {{'pre', 'post', 'noWisk'}}, 'levelNames', {{'pre', 'post', 'no whiskers'}});
%             vars.condition = struct('name', 'condition', 'levels', {{'preTrim', 'pre', 'post', 'noWisk'}}, 'levelNames', {{'preTrim', 'pre', 'post', 'no whiskers'}});
            matchConditions = {'pre' 'post'};  % for propensity score matching, math these two conditions
            figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
        end
    else
        colors = [ctlStepColor; ctlStepColor; repmat(lesionColor,4,1).*linspace(1,.25,4)'];
%         vars.condition = struct('name', 'condition', ...
%             'levels', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}}, ...
%             'levelNames', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}});
        vars.condition = struct('name', 'condition', ...  % this version has no preTrim
            'levels', {{'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}}, ...
            'levelNames', {{'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}});
        matchConditions = {'pre', 'postContra'};  % for propensity score matching, math these two conditions
    end
    
    extraVars = {'condition'};  % vars to extract using flattenData that are specific to this experiment
    
end


% load data
fprintf(['loading ' dataset ' data... ']);
tic; load(fullfile(getenv('SSD'), 'paper1', [dataset '_data.mat']), 'data'); toc
data_backup = data;  % store unaltered data so it doesn't need to be loaded again


% exclude mice
data.data = data.data(~ismember({data.data.mouse}, miceToExclude));
mice = {data.data.mouse};

% pool postContra and postBi for senLesion experiments
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
        [data.data(i).sessions(firstPreSession:firstLesSession-preSessions-1).condition] = deal(' ');
    end
end

% restrict to early lesion sessions
if contains(dataset, 'esion') && ~splitEarlyLate  % only do this for the lesion data
    fprintf('restricting post lesion sessions to %i->%i... ', earlySessions(1), earlySessions(2));
    for i = 1:length(data.data)  % loop across mice
        bins = contains({data.data(i).sessions.condition}, 'post') & ...  % if a lesion session
            ([data.data(i).sessions.conditionNum]<earlySessions(1) | ...  % and if not an early session
            [data.data(i).sessions.conditionNum]>earlySessions(2));
        [data.data(i).sessions(bins).condition] = deal(' ');              % get rid of them!      
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
vars.isBigStep = struct('name', 'isBigStep', 'levels', [0 1], 'levelNames', {{'little step', 'big step'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore', 'hind'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isModPawContra = struct('name', 'isModPawContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isHind = struct('name', 'isFore', 'condition', @(x) x==0);
conditionals.isContra = struct('name', 'isContra', 'condition', @(x) x==1);
conditionals.modPawOnlySwing = struct('name', 'modPawOnlySwing', 'condition', @(x) x==1);
conditionals.isPre = struct('name', 'condition', 'condition', @(x) strcmp(x, 'pre'));  % for restricting analyses to 'pre' condition

colorsRaw = colors;  % keep original colors before applying following transformation
if strcmp(dataset, 'mtc_muscimol') && matchTrials; colors(2,:) = mean(colors,1); end  % if matching trials, split difference btwn control and manipulated colors

fprintf('data loaded!\n');


%% plot propensity score matching

% settings
binNum = 100;  % resolution of x axis
percentileLims = [1 99];  % x axis limits
mouseAlpha = .3;
pThresh = [.05 .01 .001];       % !!! CURRENTLY MUST BE ORDERED FROM LARGEST TO SMALLEST!
symbols = {'*', '**', '***'};   % symbols associated with the pThresh values above (they will appear above the brackets connecting the conditions to be compared)

flat = flattenData(data_unmatched, [{'mouse', 'condition'} varsToMatch]);
flat_matched = flattenData(data, [{'mouse', 'condition'} varsToMatch]);

colors_unmatched = colorsRaw;
colors_matched = [colorsRaw(1,:); mean(colorsRaw,1)];

% pdfs
close all
figure('color', 'white', 'menubar', 'none', 'position', [200 595.00 600.00 345.00]);
fprintf('\nttests----------------\n')

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
        xlabel(varsToMatch{i})
        set(gca, 'XLim', xLims, 'YTick', [], 'YColor', 'none', 'tickdir', 'out')
        xTicks = get(gca, 'XTick');
        set(gca, 'XTick', xTicks([1,end]))
        
        % matched means
        subplot(length(varsToMatch),2,(i-1)*2+2); hold on
        if i==1; title('matched trials'); end
        plot(x, mean(pdfs_matched,1), 'Color', colors_matched(k,:), 'LineWidth', 3)
        fill([x(1) x x(end)], [0 mean(pdfs_matched,1) 0], colors_unmatched(k,:), 'FaceAlpha', .2, 'EdgeColor', 'none')
        set(gca, 'XLim', xLims, 'YTick', [], 'YColor', 'none', 'tickdir', 'out')
        xTicks = get(gca, 'XTick');
        set(gca, 'XTick', xTicks([1,end]))
    end
    
    % stats
    subplot(length(varsToMatch),2,(i-1)*2+1); hold on
    mat = getDvMatrix(data_unmatched, varsToMatch{i}, vars.condition, {'mouse'}, figConditionals);
    [~, p] = ttest(mat(1,:), mat(2,:));
    pInd = find(p<pThresh, 1, 'last');
    if isempty(pInd); t = 'n.s.'; else; t = symbols{pInd}; end
    text(xLims(1)+range(xLims)*.8, max(get(gca,'YLim'))*.9, t, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    fprintf('%s (unmatched) -> %.4f %s\n', varsToMatch{i}, p, t)
    
    % stats, matched
    subplot(length(varsToMatch),2,(i-1)*2+2); hold on
    mat = getDvMatrix(data, varsToMatch{i}, vars.condition, {'mouse'}, figConditionals);
    [~, p] = ttest(mat(1,:), mat(2,:));
    pInd = find(p<pThresh, 1, 'last');
    if isempty(pInd); t = 'n.s.'; else; t = symbols{pInd}; end
    text(xLims(1)+range(xLims)*.8, max(get(gca,'YLim'))*.9, t, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    fprintf('%s (matched) -> %.4f %s\n', varsToMatch{i}, p, t)
end

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_matchedHistos' suffix1 suffix2]), 'svg');
fprintf('%.3f of trials used in matched sub-population\n', size(flat_matched) / size(flat))


%% success rate
if strcmp(dataset, 'senLesion'); props = {'YLim', [.2 1], 'YTick', [.2 .6 1]}; else; props = {}; end
figure('position', [200 472.00 300 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'isTrialSuccess', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
    'comparisons', [ones(length(vars.condition.levels)-1,1), [2:length(vars.condition.levels)]'], 'test', 'ttest', props{:})
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_success' suffix1 suffix2]), 'svg');

%% velocity
figure('position', [200 472.00 300 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'trialVel', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'velocity (m/s)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0, ...
    'comparisons', [ones(length(vars.condition.levels)-1,1), [2:length(vars.condition.levels)]'], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_velocity' suffix1 suffix2]), 'svg');

%% tail height
figure('position', [400 472.00 300 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'tailHgt', vars.condition, {'mouse'}, [figConditionals])*1000;
barFancy(dv, 'ylabel', 'tail height (mm)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_tailHeight' suffix1 suffix2]), 'svg');

%% body angle
figure('position', [600 472.00 300 328.00], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'trialAngleContra', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'body angle ({\circ})', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, 'textRotation', 0)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_bodyAngle' suffix1 suffix2]), 'svg');

%% paw height (figures s5d)
figure('position', [200 100 300 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isFore; vars.isContra; vars.condition];
if strcmp(dataset, 'senLesion'); figVars = [vars.isFore; vars.condition]; end
dv = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, [figConditionals]) * 1000;
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'paw height (mm)', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'constantEdgeColor', [.15 .15 .15], ...
    'comparisons', [1 2; 3 4; 5 6; 7 8], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_pawHeight' suffix1 suffix2]), 'svg');

%% paw correlations
figure('position', [200 100 300 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.condition];
tempConditionals = [figConditionals; conditionals.isLeading; conditionals.isFore];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, figVars, {'mouse'}, {'session'}, tempConditionals, 'corr'); % perform regression for each session, then average slopes across sessions for each mouse
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'paw-obstacle correlation', 'levelNames', {figVars.levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'textRotation', 0, ...
    'comparisons', [1 2; 1 3; 1 4], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_correlations' suffix1 suffix2]), 'svg');

%% paw contacts (figures s5b)
figure('position', [200 100 300 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isFore; vars.isContra; vars.condition];
if strcmp(dataset, 'senLesion'); figVars = [vars.condition; vars.isFore]; end
dv = 1-getDvMatrix(data, 'anyTouchFrames', figVars, {'mouse'}, [figConditionals]);
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'constantEdgeColor', [.15 .15 .15], ...
    'comparisons', [1 2; 3 4; 5 6; 7 8], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_pawContacts' suffix1 suffix2]), 'svg');

%% baseline step heights (figures s5a)
figure('position', [200 100 300 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isFore; vars.isContra; vars.condition];
dv = getDvMatrix(data, 'controlStepHgt', figVars, {'mouse'}, [figConditionals])*1000;
barFancy(dv, 'ylabel', 'control step height (mm)', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,4,1), barProperties{:}, ...
    'comparisons', [1 2; 3 4; 5 6; 7 8], 'test', 'ttest', 'addLabelLines', false)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_baselineHeights' suffix1 suffix2]), 'svg');

%% ventral contacts,  'grabs' (figures s5f)
figure('position', [200 100 300 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isFore; vars.isContra; vars.condition];
dv = getDvMatrix(data, 'isVentralContact', figVars, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'ventral contact rate', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,4,1), barProperties{:}, ...
    'comparisons', [1 2; 3 4; 5 6; 7 8], 'test', 'ttest', 'addLabelLines', false)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_ventralContacts' suffix1 suffix2]), 'svg');

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

sessionsToShow = -(preSessions-1):7;
manipInd = find(sessionsToShow==0);
plotIpsi = true;  % whether to align everything to the ipsilateral lesion

    
% add sessionsPostLesion to data structure
mouseBins = true(1,length(mice));
vars.sessionsPostLesion = struct('name', 'sessionsPostLesion', 'levels', sessionsToShow);
for i = 1:length(data.data) % compute session number relative to first lesion session
    sessionNums = [data.data(i).sessions.sessionNum];
    
    if plotIpsi
        firstLesSession = sessionNums(find(strcmp({data.data(i).sessions.condition}, 'postIpsi'), 1, 'first'));  % we don't want 'postIpsi' getting in here!!!
    else
        firstLesSession = sessionNums(find(ismember({data.data(i).sessions.condition}, {'post', 'postEarly', 'postLate', 'postContra', 'postBi'}), 1, 'first'));  % we don't want 'postIpsi' getting in here!!!
    end
    
    if ~isempty(firstLesSession)
        sessionsPostLesion = num2cell([data.data(i).sessions.sessionNum] - firstLesSession+1);
        [data.data(i).sessions.sessionsPostLesion] = sessionsPostLesion{:};

        if strcmp(dataset, 'senLesion')  % make sure preTrim not included in pre, and postTrim not included in post
            for j = 1:length(data.data(i).sessions)
                if ismember(data.data(i).sessions(j).condition, {'preTrim', 'noWisk'})
                    data.data(i).sessions(j).sessionsPostLesion = nan;
                end
            end
        end
    else
        [data.data(i).sessions.sessionsPostLesion] = deal(nan);
        mouseBins(i) = false;
    end
end


%% success
figure('position', [200 400 521.00 239.00], 'color', 'white', 'menubar', 'none')
if plotIpsi; xLims = [sessionsToShow(1) 5]; else; xLims = [sessionsToShow(1) sessionsToShow(end)]; end
dv = getDvMatrix(data.data(mouseBins), 'isTrialSuccess', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'success rate', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:}, 'xlim', xLims);
set(gca, 'YLim', [0 1], 'YTick', 0:.5:1)
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
if plotIpsi
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_succesOverSessionsIpsi' suffix1 suffix2]), 'svg');
else
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_succesOverSessions' suffix1 suffix2]), 'svg');
end

%% body angle
figure('position', [200 200 521.00 239.00], 'color', 'white', 'menubar', 'none')
dv = getDvMatrix(data, 'trialAngleContra', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'body angle', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
set(gca, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_angleOverSessions' suffix1 suffix2]), 'svg');

%% tail height
figure('position', [200 300 521.00 239.00], 'color', 'white', 'menubar', 'none')
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
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'leading forepaw height (mm)', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:});
if strcmp(dataset, 'senLesion'); set(gca, 'YLim', [6 10.5]); end
set(gca, 'XLim', [sessionsToShow(1) sessionsToShow(end)])
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_forepawHgtOverSessions' suffix1 suffix2]), 'svg');

%% hindpaw height
figure('position', [200 500 521.00 239.00], 'color', 'white', 'menubar', 'none')
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

%% success over days post wisk trim (senLesion only)

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
else
    colors_temp = colors;
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


%% decision making initialization
flat = flattenData(data, ...
    [m.predictorsAll, {'mouse', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', 'isModPawContra', ...
    'modPawOnlySwing', 'isTrialSuccess', 'condition', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'isWheelBreak', ...
    'modPawKinInterp', 'preModPawKinInterp', 'firstModPaw', 'preModPawDeltaLength', 'modSwingContacts', 'session'}]);

%% heatmaps (figures f6b f6e s6f)
if strcmp(dataset, 'senLesion'); dims = [2 2]; else; dims = []; end
plotDecisionHeatmaps(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'normalize', 'col', 'xLims', [-20 15], ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', false, 'plotMice', false, 'colors', colors, 'outcome', 'isModPawLengthened', 'plotProbs', false, 'subplotDims', dims, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_heatMaps' suffix1 suffix2]));

%% (temp; to see individual mice) heatmaps
plotDecisionHeatmaps(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'normalize', 'col', 'xLims', [-20 15], ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', true, 'plotMice', true, 'colors', colors, 'outcome', 'isModPawLengthened', 'plotProbs', false, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_heatMaps' suffix1 suffix2]));

%% (temp; to compare ispi and contra paw first) heatmaps
plotDecisionHeatmaps(flat(strcmp({flat.condition}, 'pre')), 'condition', 'isModPawContra', 'levels', {0,1}, 'normalize', 'col', 'xLims', [-20 15], ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', true, 'plotMice', true, 'colors', colors, 'outcome', 'isModPawLengthened', 'plotProbs', false, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_heatMaps' suffix1 suffix2]));

%% landing position entropy
plotEntropies(flat, 'condition', 'condition', 'levels', vars.condition.levels, ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', colors, 'barProps', [barProperties, {'comparisons', [1 2]}]);

%% trial scatters (figures f6a f6d s6g)
plotDecisionTrials(flat([flat.isWheelBreak]), 'condition', 'condition', 'levels', vars.condition.levels, 'outcome', 'isBigStep', ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'rowColors', colors, 'xLims', [-.11 .06], 'obsColor', obsColor, 'poolHistos', true, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_decisionKin' suffix1 suffix2]));

%% model accuracies (figures s6a s6c s6h)
rng(1)
if strcmp(dataset, 'senLesion'); comparisons = [2 4; 2 6; 2 8]; xfers = []; else; comparisons = [2 4]; xfers = []; end
[accuracies, ~, temp] = plotModelAccuracies(flat, m.predictors, 'isModPawLengthened', 'modelTransfers', xfers, ...
    'weightClasses', true, 'condition', 'condition', 'levels', vars.condition.levels, 'kFolds', 15, ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', colors, 'barProps', [barProperties {'comparisons', comparisons, 'constantEdgeColor', [.15 .15 .15]}], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_models' suffix1 suffix2]));

%% (temp, for including preTrim condition) model accuracies
plotModelAccuracies(flat, m.predictorsAll, 'isModPawLengthened', 'modelTransfers', [], ...
    'weightClasses', true, 'condition', 'condition', 'levels', vars.condition.levels, 'kFolds', 10, ...
    'modSwingContactsMax', false, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'barProps', [barProperties {'YLim', [.2 1], 'comparisons', [2 4; 2 6; 2 8], 'constantEdgeColor', [.15 .15 .15]}]);

%% landing position entropy
plotEntropies(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'outcome', 'isModPawLengthened', ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'modSwingContactsMax', m.modSwingContactsMax, 'colors', colors, 'barProps', barProperties);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_landingHistos' suffix1 suffix2]), 'svg');

%% model predictors (figures f6c s6e s6j)
% predictorsToShow = m.predictors;
predictorsToShow = m.predictors(1:4);

[~, inds] = ismember(predictorsToShow, m.predictors);
plotPredictors(flat, m.predictors(inds), 'isModPawLengthened', 'colors', colors, 'avgMice', true, ...
    'condition', 'condition', 'levels', vars.condition.levels(1:2), 'overlayConditions', true, 'mouseAlpha', .15, 'subplotDims', [2 2], 'names', m.predictorsNamed(inds), ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', true, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_predictors' suffix1 suffix2]));

%% distribution bimodality
close all
plotBimodalities(flat, 'bic', true, 'condition', 'condition', 'levels', vars.condition.levels, 'outcome', 'isModPawLengthened', ...
    'modSwingContactsMax', 0, 'deltaMin', m.deltaMin, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', colors, 'barProps', [barProperties, 'comparisons', [2 4], 'test', 'ttest'], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceBimodality'));

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


%% -------------------------------
%  FIGURES WITH SUBPLOTS FOR PAPER
%  -------------------------------

%% success rate and correlation scatters (stacked on top of one another for final figure)

% settings
limsSuccess = [0 1];
limsCorr = [-.1 .8];

close all
figure('position', [173 115 241 489], 'color', 'white', 'menubar', 'none');

% success
subplot(2,1,1)
dv = getDvMatrix(data, 'isTrialSuccess', vars.condition, {'mouse'}, [figConditionals]);
pairedScatter(dv, 'ylabel', matchConditions{2}, 'title', 'success rate', ...
    'colors', colors(1:2,:), 'lims', limsSuccess)

% correlation
subplot(2,1,2)
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.condition, {'mouse'}, {'session'}, [conditionals.isFore; conditionals.isLeading], 'corr'); % perform regression for each session, then average slopes across sessions for each mouse
pairedScatter(dv, 'xlabel', matchConditions{1}, 'ylabel', matchConditions{2}, 'title', 'paw-obstacle correlation', ...
    'colors', colors(1:2,:), 'lims', limsCorr)

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_dvScatters' suffix1 suffix2]), 'svg');


%% paw height + paw contacts

figure('position', [200 100 700 328.00], 'color', 'white', 'menubar', 'none');

% paw height
subplot(1,2,1)
figVars = [vars.isContra; vars.isFore; vars.condition];
if strcmp(dataset, 'senLesion'); figVars = [vars.isFore; vars.condition]; end
dv = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, [figConditionals]) * 1000;
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'paw height (mm)', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'constantEdgeColor', [.15 .15 .15], ...
    'YTick', [6 10 14], 'YLim', [6 14], ...
    'comparisons', [1 2; 3 4; 5 6; 7 8], 'test', 'ttest')

% paw contacts
subplot(1,2,2)
figVars = [vars.isContra; vars.isFore; vars.condition];
if strcmp(dataset, 'senLesion'); figVars = [vars.isFore; vars.condition]; end
dv = 1-getDvMatrix(data, 'anyTouchFrames', figVars, {'mouse'}, [figConditionals]);
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {figVars(1:end-1).levelNames}, ...
    'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'constantEdgeColor', [.15 .15 .15], ...
    'comparisons', [1 2; 3 4; 5 6; 7 8], 'test', 'ttest', 'YTick', [0 .5 1], 'YLim', [0 1])
subplot(1,2,2)
% % the following code also breaks things down by leading / lagging
% figVars = [vars.isContra; vars.isFore; vars.isLeading; vars.condition];
% if strcmp(dataset, 'senLesion'); figVars = [vars.isFore; vars.condition]; end
% dv = 1-getDvMatrix(data, 'anyTouchFrames', figVars, {'mouse'}, [figConditionals]);
% colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
% barFancy(dv, 'ylabel', 'success rate', 'levelNames', {figVars(1:end-1).levelNames}, ...
%     'colors', repmat(colors,colorRepeats,1), barProperties{:}, 'constantEdgeColor', [.15 .15 .15], ...
%     'comparisons', reshape(1:16, 2, 8)', 'test', 'ttest', 'YTick', [0 .5 1], 'YLim', [0 1])

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_heightAndSuccess' suffix1 suffix2]), 'svg');

%% success rate + velocity over time (need to run 'sessions over time' with plotIpsi=false code above to initialize)


close all
figure('position', [1444 224 538 457], 'color', 'white', 'menubar', 'none')

% success
subplot(2,1,1)
dv = getDvMatrix(data, 'isTrialSuccess', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'success rate', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:}, 'xlim', vars.sessionsPostLesion.levels([1 end]));
set(gca, 'YLim', [0 1], 'YTick', 0:.5:1)
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')

% velocity
subplot(2,1,2)
dv = getDvMatrix(data, 'velAtWiskContact', vars.sessionsPostLesion, {'mouse'});
sesPlotRick(dv', 'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'velocity (m/s)', 'xlabel', 'sessions post lesion', ...
    'compareTo', 1:manipInd, sessionPlotProperties{:}, 'xlim', vars.sessionsPostLesion.levels([1 end]));
ln = line([.5 .5], get(gca, 'ylim'), 'color', [lesionColor .9], 'linewidth', 2); uistack(ln, 'bottom')

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_successAndVelOverSessions' suffix1 suffix2]), 'svg');

%% landing distance bar plots (figures s6b s6d s6i)

% landing distance
if strcmp(dataset, 'senLesion'); cmp = [1 2; 1 3; 1 4; 5 6; 5 7; 5 8]; edgeColor=[.15 .15 .15]; else; cmp = [1 2; 3 4]; edgeColor = []; end
figure('position', [200 472.00 300 255], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'modPawDistanceToObsAbs', [vars.isBigStep; vars.condition], {'mouse'}, [figConditionals; conditionals.modPawOnlySwing])*1000;
barFancy(dv, 'ylabel', 'landing distance (mm)', 'levelNames', {vars.isBigStep.levelNames}, 'colors', repmat(colors,2,1), barProperties{:}, ...
    'comparisons', cmp, 'test', 'ttest', 'constantEdgeColor', edgeColor)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_landingDistance' suffix1 suffix2]), 'svg');

%% (temp) landing distance bar plots broken down by isModPawContra

% landing distance
figure('position', [200 472.00 300 255], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'modPawDistanceToObsAbs', [vars.isModPawContra], {'mouse'}, [figConditionals; conditionals.isPre])*1000;
barFancy(dv, 'ylabel', 'landing distance (mm)', 'levelNames', {vars.isModPawContra.levelNames}, barProperties{:}, ...
    'comparisons', [1 2], 'test', 'ttest')

%% success and paw height as fcn of obs height
close all
flat = struct2table(flattenData(data, {'mouse', 'session', 'trial', 'condition', 'conditionNum', 'isTrialSuccess', 'isLeading', 'isFore', 'obsHgt', 'preObsHgt'}));
[~, conditions] = ismember(flat.condition, vars.condition.levels);

% obs height vs paw height
figure('Color', 'white', 'Position', [200 400 400 300], 'MenuBar', 'none');
plot([0 10], [0 10], 'Color', [obsColor .4], 'LineWidth', 3) % add unity line
lfBins = [flat.isLeading] & [flat.isFore];  % leading forepaw bins
logPlotRick(flat.obsHgt(lfBins)*1000, flat.preObsHgt(lfBins)*1000, ...
    'colors', colors, 'conditions', conditions(lfBins), 'xlabel', 'obstcle height (mm)', 'ylabel', 'paw height (mm)', 'plotMice', false, ...
    'xlim', [4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', flat.mouse(lfBins), ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)))
set(gca, 'xlim', [4 10], 'YTick', 4:2:12, 'YLim', [4 12])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_heightShapingMovingAverages' suffix1 suffix2]), 'svg');

% obs height vs success rate
figure('Color', 'white', 'Position', [800 400 400 300], 'MenuBar', 'none');
lfBins = [flat.isLeading] & [flat.isFore];  % leading forepaw bins
logPlotRick(flat.obsHgt(lfBins)*1000, flat.isTrialSuccess(lfBins), ...
    'colors', colors, 'conditions', conditions(lfBins), 'xlabel', 'obstcle height (mm)', 'ylabel', 'success rate', 'plotMice', false, ...
    'xlim', [4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', flat.mouse(lfBins), ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)))
set(gca, 'xlim', [4 10], 'YLim', [0 1])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_successVsHeight' suffix1 suffix2]), 'svg');


%% compare ipsi and contra lesion effect sizes
% (need to set poolSenLesionConditions=false in inits)
poolSenLesionConditions










