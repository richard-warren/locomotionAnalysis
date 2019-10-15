%% INITIALIZATIONS

% load global settings
clear all; tcm190910_config;


% settings
dataset = 'senLesion';
poolSenLesionConditions = true;  % whether to use all conditions or pool postBi and postContra
maxEarlySessions = 3;  % only include this many days post lesion in lesions figs

% motor cortex muscimol
if strcmp(dataset, 'mtc_muscimol')
    colors = manipColors;
    vars.condition = struct('name', 'condition', 'levels', {{'saline', 'muscimol'}}, 'levelNames', {{'saline', 'muscimol'}});
    figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
    extraVars = {'condition'};  % vars to extract using flattenData that are specific to this experiment
    
% motor cortex lesion
elseif strcmp(dataset, 'mtc_lesion')
    colors = manipColors;
    vars.condition = struct('name', 'condition', 'levels', {{'pre', 'post'}}, 'levelNames', {{'pre', 'post'}});
	figConditionals = struct('name', 'conditionNum', 'condition', @(x) x<=maxEarlySessions);
    extraVars = {'condition'};  % vars to extract using flattenData that are specific to this experiment

% barrel cortex lesion
elseif strcmp(dataset, 'senLesion')
    
    if poolSenLesionConditions
        colors = [manipColors; manipColors(2,:)*.5];
        vars.condition = struct('name', 'condition', ...
            'levels', {{'pre', 'postBi', 'noWisk'}}, ...
            'levelNames', {{'pre', 'post', 'no whiskers'}});
    else
        colors = [axisColor; axisColor; winter(4)];
        vars.condition = struct('name', 'condition', ...
            'levels', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}}, ...
            'levelNames', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}});
    end
    figConditionals = struct('name', 'conditionNum', 'condition', @(x) x<=maxEarlySessions);
    extraVars = {'condition'};  % vars to extract using flattenData that are specific to this experiment

% opto (settings the same for all brain regions)
elseif contains(dataset, 'opto')
    colors = optoColors;
    vars.condition = struct('name', 'powerCondition', 'levels', [1 2 3 4], 'levelNames', {{'0', '1.4', '4.5', '15mW'}});
    figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
    extraVars = {'powerCondition', 'optoOnPositions', 'isOptoOn'};  % vars to extract using flattenData that are specific to this experiment
end

% load data
fprintf(['loading ' dataset ' data... ']);
load(fullfile(getenv('OBSDATADIR'), 'matlabData', [dataset '_data.mat']), 'data');

% pool postContra and postBi for senLesion experiments
if strcmp(dataset, 'senLesion') && poolSenLesionConditions
    fprintf('pooling postContra and postBi conditions... ');
    for i = 1:length(data.data)  % loop across mice
        for j = 1:length(data.data(i).sessions)
            if strcmp(data.data(i).sessions(j).condition, 'postContra')
                data.data(i).sessions(j).condition = 'postBi';
            end
        end
    end
end

flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'obsHgt', 'isPawSuccess', ...
    'isLeading', 'isFore', 'stepOverKinInterp', 'paw', 'controlStepKinInterp', 'preObsHgt', 'stepType', ...
    'isBigStep', 'preObsKin', 'conditionNum', extraVars{:}});

if contains(dataset, {'lesion', 'Lesion'})
    fprintf('restricting flat to early sessions... ');
    flat = flat([flat.conditionNum]<=maxEarlySessions);
end
disp([dataset ' data loaded!'])

% define variables and conditionals to be used in bar plots
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
figPadding = (length(vars.condition.levels)-2) * .25;
conditionals.isObsHigh = struct('name', 'obsHgt', 'condition', @(x) x>.008);
conditionals.isObsLow = struct('name', 'obsHgt', 'condition', @(x) x<.005);


%% MAKE FIGURE WITH SUCCESS, PAW HEIGHT, SHAPING, AND VEL

% SUCCESS
if contains(dataset, 'opto')
    cols = 3;
    figWidth = 9;
else
    cols = 2;
    figWidth = 6.28;
end

figure('name', dataset, 'units', 'inches', 'position', [21.18 3.21 figWidth 7.5], ...
    'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on
subplot(2,cols,1)
% figure('units', 'inches', 'position', [6.98 5.07 2.44+figPadding figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
dv = getDvMatrix(data, 'isTrialSuccess', vars.condition, {'mouse'}, [figConditionals]);
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
    'YLim', [0 1], 'YTick', [0 .5 1], 'barWidth', .8)
% set(gca, 'position', [.2 .11 .78 .82])
% print -clipboard -dmeta

% PAW HEIGHT
subplot(2,cols,2)
% figure('units', 'inches', 'position', [6.98 5.07 2.44+figPadding figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
dv = getDvMatrix(data, 'preObsHgt', vars.condition, {'mouse'}, [figConditionals; conditionals.isFore])*1000;
barFancy(dv, 'ylabel', 'paw height (mm)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
    'YTick', 0:2:12, 'barWidth', .8, 'YLim', [0 12])
% set(gca, 'position', [.2 .11 .78 .82])
% print -clipboard -dmeta

% PAW SHAPING
if contains(dataset, 'opto')
    subplot(2,cols,3)
    dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.condition, {'mouse'}, {'session'}, ...
        [figConditionals; conditionals.isFore; conditionals.isLeading], 'corr');
    barFancy(dv, 'levelNames', {vars.condition.levelNames}, 'ylabel', 'paw-obstacle correlation', 'colors', colors, barProperties{:}, ...
        'YLim', [-.4 .8])
    print -clipboard -dmeta
end


% SPEED VS POSITION
subplot(2,cols,cols+1:cols*2)
yLims = [.2 .7];
% figure('units', 'inches', 'position', [21.18 3.21 6.28 figHgt*.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on

% add obstacle rectangle and lines
hold on
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
line([0 0], yLims, 'linewidth', 2, 'color', axisColor)
line([x(1) x(1)], yLims, 'linewidth', 2, 'color', axisColor)
text(x(1), yLims(2), '\itobstacle on', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font)
text(0, yLims(2), '\itobstacle reached', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font)

if ~contains(dataset, 'opto')  % plots for non-opto experiments
    plotDvPsth(flat, 'velVsPosition', vars.condition.name, ...
        {'showLegend', false, 'conditionColors', colors, 'xlim', [-.5 .2], 'showErrorBars', true, ... 
         'plotConditions', vars.condition.levels, 'errorAlpha', .3, 'lineWidth', 3})

else  % plots for opto experiments
    
    % add light on, obs on, and obs at nose markers
    optoOnPos = nanmedian([flat.optoOnPositions]);
    x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
    line([optoOnPos optoOnPos], yLims, 'linewidth', 4, 'color', mean(colors(2:end,:),1))
    text(optoOnPos, yLims(2), '\itlaser on', 'color', axisColor, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font, 'color', mean(colors(2:end,:),1))

    plotDvPsth(flat, 'velVsPosition', vars.condition.name, ...
        {'showLegend', false, 'conditionColors', colors, 'xlim', [-.9 .2], 'showErrorBars', false, ... 
         'plotConditions', vars.condition.levels, 'errorAlpha', .3, 'lineWidth', 3})
end

set(gca, 'YLim', yLims, 'YTick', linspace(yLims(1),yLims(2),3));
set(gca, 'YLim', yLims);
xlabel('distance (m)')
ylabel('velocity (m/s)')
set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black', 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta


%% PAW SHAPING

figure('units', 'inches', 'position', [6.98 5.07 2.44+figPadding figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.condition, {'mouse'}, {'session'}, ...
    [figConditionals; conditionals.isFore; conditionals.isLeading], 'corr');
barFancy(dv, 'levelNames', {vars.condition.levelNames}, 'ylabel', 'paw-obstacle correlation', 'colors', colors, barProperties{:})
set(gca, 'position', [.2 .11 .78 .82])
print -clipboard -dmeta


% %% SHAPING MOVING AVGS
% 
% flat_sub = flat(~isnan([flat.obsHgt]) & ...
%                 [flat.isLeading] & ...
%                 [flat.isFore]);
% if ~contains(dataset, 'opto')
%     [~, conditions] = ismember({flat_sub.(vars.condition.name)}, vars.condition.levels);  % convert condition name into condition number
%     colorsTemp = colors;
% else
%     conditions = [flat_sub.isOptoOn]+1;
%     colorsTemp = manipColors;
% end
% 
% figure('units', 'inches', 'position', [20.95 1.88 3.80 3.7], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
% plot([0 100], [0 100], 'Color', axisColor*.5, 'LineWidth', 3) % add unity line
% logPlotRick([flat_sub.obsHgt]*1000, [flat_sub.preObsHgt]*1000, ...
%     {'colors', colorsTemp, 'conditions', conditions, 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)',  ...
%     'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat_sub.mouse}, ...
%     'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'lineWidth', 2})
% set(gca, 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor, 'FontSize', fontSize, 'FontName', font)
% print -clipboard -dmeta

%% KINEMATICS

% settings
obsHgtBins = 4; % discretize obstacle heights into this many bins
fading = .5; % within a plot, each paw's color fades from fading*color -> color
xLims = [-.058 0];
yLims = [0 .016];

% initializations
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat_sub = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLeading] & ...
            [flat.isFore]);
if ~contains(dataset, 'opto')
    [~, conditions] = ismember({flat_sub.(vars.condition.name)}, vars.condition.levels);  % convert condition name into condition number
    colorsTemp = colors;
else
    conditions = [flat_sub.isOptoOn]+1;
    colorsTemp = manipColors;
end

figure('units', 'inches', 'position', [10.67 3.16 4.51 1.44*max(conditions)], ...
    'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on
kinData = permute(cat(3, flat_sub.preObsKin), [3,1,2]);
kinDataCtl = permute(cat(3, flat_sub.controlStepKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1); % change the x starting x position of ctl steps to match steps over

rows = max(conditions);
set(gca, 'visible', 'off')
for i = 1:rows
    ax = axes('position', [0 1-i*(1/rows) 1 1/rows ], 'color', 'black');
    bins = conditions==i;
    plotColor = repmat(colorsTemp(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    plotKinematics(kinData(bins,[1,3],:), [flat_sub(bins).obsHgt], [flat_sub(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'mouseNames', {flat_sub(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims, 'color', 'black')
end
print -clipboard -dmeta


%% SESSIONS OVER TIME

close all

if strcmp(dataset, 'mtc_lesion')
    sessionsToShow = -2:8;
    postCondition = {'post'};
elseif strcmp(dataset, 'senLesion')
    sessionsToShow = -3:6;
    postCondition = {'postBi', 'postContra'};
end

if strcmp(dataset, 'mtc_lesion') || (strcmp(dataset, 'senLesion') && poolSenLesionConditions)
    
    axProps = {'color', 'black', 'xcolor', axisColor', 'ycolor', axisColor, 'FontName', font, 'FontSize', fontSize};
    
    % initializations
    figure('units', 'inches', 'position', [2 2 6 6], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
    vars.sessionsPostLesion = struct('name', 'sessionsPostLesion', 'levels', sessionsToShow);

    for i = 1:length(data.data) % compute session number relative to first lesion session
        sessionNums = [data.data(i).sessions.sessionNum];
        firstLesSession = sessionNums(find(ismember({data.data(i).sessions.condition}, postCondition), 1, 'first'));
        sessionsPostLesion = num2cell([data.data(i).sessions.sessionNum] - firstLesSession);
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
    subplot(2,1,1);
    dv = getDvMatrix(data, 'isTrialSuccess', vars.sessionsPostLesion, {'mouse'});
    sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'success rate', ...
        'meanColor', axisColor, 'colors', mouseColors, 'alpha', .6});
    set(gca, axProps{:}, 'ylim', [0 1])
    ln = line(-[.5 .5], get(gca, 'ylim'), 'color', manipColors(2,:), 'linewidth', 2); uistack(ln, 'bottom')

    % velocity
    subplot(2,1,2);
    dv = getDvMatrix(data, 'velAtWiskContact', vars.sessionsPostLesion, {'mouse'});
    sesPlotRick(dv', {'xvals', vars.sessionsPostLesion.levels, 'ylabel', 'velocity (m/s)', 'xlabel', 'days post lesion', ...
        'meanColor', axisColor, 'colors', mouseColors, 'alpha', .6});
    ln = line(-[.5 .5], get(gca, 'ylim'), 'color', manipColors(2,:), 'linewidth', 2); uistack(ln, 'bottom')
    set(gca, axProps{:})
    
    print -clipboard -dmeta
end



%% PROPENSITY SCORE MATCHING


% settings
varsToMatch = {'velAtWiskContact', 'angleAtWiskContactContra', 'tailHgtAtWiskContact'};
manipPercent = 25; % take manipPercent percent of best matched manip trials
manipConditions = {'post', 'muscimol'};  % the names of the conditions corresponding to the manipulation


% propensity score matching
mice = {data.data.mouse};
flat_sub = struct2table(flattenData(data, [{'mouse', 'session', 'trial', 'condition', 'isLightOn', 'conditionNum'} varsToMatch]));
varBins = ismember(flat_sub.Properties.VariableNames, varsToMatch);
metaBins = ismember(flat_sub.Properties.VariableNames, {'mouse', 'session', 'trial'});  % metadata associated with each trial

% find matched trials
matchedTrials = cell2table(cell(0,3), 'VariableNames', {'mouse', 'session', 'trial'});
for mouse = mice
    
    bins = strcmp(flat_sub.mouse, mouse{1}) & ...
           ismember(flat_sub.condition, vars.condition.levels) & ...
           [flat_sub.conditionNum]<=maxEarlySessions;
    flat_sub_sub = flat_sub(bins,:);
    
    
	X = table2array(flat_sub_sub(:, varBins));
    y = ismember(flat_sub_sub.condition, manipConditions); % is trial in the manip condition
    matchedPairs = propensityMatching(X, y, ...
        {'percentileThresh', manipPercent, 'predictorNames', varsToMatch, 'verbose', false});
    matchedTrials = [matchedTrials; flat_sub_sub(matchedPairs(:), metaBins)];
    
end

% get rid of non-matched trials!
data_matched = data;

for i = 1:length(data_matched.data)
    for j = 1:length(data_matched.data(i).sessions)
        bins = strcmp(matchedTrials.mouse, data_matched.data(i).mouse) & ...
               strcmp(matchedTrials.session, data_matched.data(i).sessions(j).session);
        sesTrials = matchedTrials.trial(bins);
        data_matched.data(i).sessions(j).trials = data_matched.data(i).sessions(j).trials(sesTrials);
    end

    % remove unused sessions
    isSessionUsed = ismember({data_matched.data(i).sessions.session}, unique(matchedTrials.session));
    data_matched.data(i).sessions = data_matched.data(i).sessions(isSessionUsed);
end

% plot!
figure('name', dataset, 'units', 'inches', 'position', [2 2 12 6], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)

datasets = {data; data_matched};
plotInd = 1;
for i = 1:length(datasets)
    subplot(2,5,1+(i-1)*5)
    dv = getDvMatrix(datasets{i}, 'isTrialSuccess', vars.condition, {'mouse'}, figConditionals);
    barFancy(dv, 'ylabel', 'success rate', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
        'barWidth', .75, 'YLim', [0 1], 'YTick', [0 .5 1], 'ylabelPosX', -.7);

    subplot(2,5,3+(i-1)*5)
    dv = getDvMatrix(datasets{i}, 'velAtWiskContact', vars.condition, {'mouse'}, figConditionals);
    barFancy(dv, 'ylabel', 'velocity', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
        'barWidth', .75, 'YLim', [0 .6], 'YTick', [0 .3 .6], 'ylabelPosX', -.7);

    subplot(2,5,4+(i-1)*5)
    dv = getDvMatrix(datasets{i}, 'angleAtWiskContactContra', vars.condition, {'mouse'}, figConditionals);
    barFancy(dv, 'ylabel', 'body angle', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
        'barWidth', .75, 'YLim', [-12 12], 'YTick', [-12 0 12], 'ylabelPosX', -.7);

    subplot(2,5,5+(i-1)*5)
    dv = getDvMatrix(datasets{i}, 'tailHgtAtWiskContact', vars.condition, {'mouse'}, figConditionals) * 1000;
    barFancy(dv, 'ylabel', 'tail height (mm)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
        'barWidth', .75, 'YLim', [0 20], 'YTick', 0:10:20, 'ylabelPosX', -.7); 
end

print -clipboard -dmeta





