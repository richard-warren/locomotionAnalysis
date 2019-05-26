%% INITIALIZATIONS


% global settings
varsToAvg = {'mouse'};
earlySessions = 1:3;
lateSessions = 5:7;
colors = [217, 65, 244; 244, 149, 66] / 255; % mus, lesion
darkening = .25; % how much to darken control condition relative to manipulated condition



% initializations
vars.conditionMus = struct('name', 'condition', 'levels', {{'saline', 'muscimol'}}, 'levelNames', {{'saline', 'muscimol'}});
vars.conditionLes = struct('name', 'condition', 'levels', {{'pre', 'post'}}, 'levelNames', {{'pre', 'post'}});
vars.condition = struct('name', 'conditionNew', 'levels', [1 2 3 4 5], 'levelNames', {{'saline', 'muscimol', 'pre lesion', 'post early', 'post late'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.conditionnum = struct('name', 'conditionNum', 'levels', 1:8);

conditionNames = cat(2, vars.conditionMus.levelNames, {'pre lesion', 'post early', 'post late'});

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
for i = 1:length(dataMus)
    
    data(i,1).mouse = dataMus(i).data.mouse;
    data(i,1).sessions = [dataMus(i).data.sessions; dataLes(i).data.sessions];
    
    bins = strcmp({data(i).sessions.condition}, 'saline');
    for j = 1:length(data(i).sessions)
        
        % get condition
        if strcmp({data(i).sessions(j).condition}, 'saline'); c = 1;
        elseif strcmp({data(i).sessions(j).condition}, 'muscimol'); c = 2;
        elseif strcmp({data(i).sessions(j).condition}, 'pre'); c = 3;
        elseif strcmp({data(i).sessions(j).condition}, 'post') & ismember(data(i).sessions(j).conditionNum, earlySessions); c = 4;
        elseif strcmp({data(i).sessions(j).condition}, 'post') & ismember(data(i).sessions(j).conditionNum, lateSessions); c = 5;
        else; c = nan; end
        data(i).sessions(j).conditionNew = c;
    end
end




%% ----------
% PLOT THINGS
%  ----------


%% BARS


% initializations
rows = 5;
cols = 2;
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1200 200*rows])

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

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'motorCortex', 'bars');
saveas(gcf, file, 'svg');

%% VARS OF TIME

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
saveas(gcf, file, 'svg');



