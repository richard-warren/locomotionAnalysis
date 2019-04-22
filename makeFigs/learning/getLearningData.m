% manually creating metadata spreadsheet for learning is too difficult due
% to number of sessions, so this automatically creates table with learning
% metadata that can then be passed to getExperimentData to analyze the
% sessions // user enters names of mice to include in the analysis // user
% needs to manually check that only mice are included that have continuous
% days of records during the learning phase

% TO DO: could get quite a few more mice if i allow single skipped
% sessions, but would need to determine whether these are sessions where
% there was practice, but no recording, or sessions where there was water
% in cage, or no obstacles, for example... // fix 180702 sesions with no
% foil on obstacle... fuck!


% settings
noBrSessions = 2;
brSessions = 5;
wildTypeOnly = true;
vars = {'condition', 'conditionNum', 'sessionNum', 'isTrialSuccess', 'trialVel', ...
    'velVsPosition', 'isWheelBreak'}; % only compute these variables for each session




% compute learning metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'sessions');
mouseInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'mice');
mice = unique(sessionInfo.mouse);
conditions = [repmat({'noWheelBreaks'},1,noBrSessions), repmat({'wheelBreaks'},1,brSessions)]';
learningInfo = table('Size', [0,4], ...
    'VariableNames', {'session', 'mouse', 'condition', 'include'}, ...
    'VariableTypes', {'cell', 'cell', 'cell', 'logical'});


for i = 1:length(mice)
    
    % only perform analysis on WT
    
    genotype = mouseInfo.line_Genotype(strcmp(mouseInfo.mouse, mice{i}));
    if isequal(genotype, {'C57/B6'}) && wildTypeOnly
        
        % find last no break sessions and first break sessions
        bins = strcmp(sessionInfo.mouse, mice{i}) & logical([sessionInfo.include]);
        breakInds = find(strcmp(sessionInfo.experiment, 'obsBrBar') & bins, brSessions, 'first'); % inds in sessionInfo of rows containing first break sessions
        if any(breakInds)
            noBreakInds = find(strcmp(sessionInfo.experiment, 'obsNoBrBar') & bins & ...
                               [1:height(sessionInfo)]'<breakInds(1), noBrSessions, 'last'); % inds in sessionInfo of rows containing last no break sessions (last conditional ensures we don't take no break sessions occuring after initial break sessions)
        else
            noBreakInds = [];
        end
        
        % make sure all sessions occur on consecutive days
        sessions = sessionInfo.session([noBreakInds; breakInds]);
        if length(sessions)==(noBrSessions+brSessions)
            dates = cellfun(@(x) datetime(str2num(['20' x(1:2)]), str2num(x(3:4)), str2num(x(5:6))), sessions); % this parses the names of the sessions into matlab datetime objects
            diffs = diff(dates);
            diffs.Format = 'd';
            areDaysConsecutive = all(diffs==days(1));

            if areDaysConsecutive
                inds = height(learningInfo)+1:height(learningInfo)+noBrSessions+brSessions;
                learningInfo.session(inds) = sessions;
                learningInfo.mouse(inds) = {mice{i}};
                learningInfo.condition(inds) = conditions;
                learningInfo.include(inds) = 1;
            end
        end
    end
end

sessionInfo = learningInfo;
mice = unique(sessionInfo.mouse);

%% compute kinData for all sessions (only need to do once)
overwriteKindata = false;
sessions = unique(sessionInfo(logical([sessionInfo.include]),:).session);
for i = 1:length(sessions)
    if ~exist(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'kinData.mat'), 'file') || overwriteKindata
        getKinematicData5(sessions{i});
    end
end

%% overwrite kinData for a given mice
mice = {'sen2', 'sen3', 'sen6'};
sessions = unique(sessionInfo(logical([sessionInfo.include]) & ismember(sessionInfo.mouse, mice),:).session);
for i = 1:length(sessions)
    try; getKinematicData5(sessions{i}); catch; fprintf('%s: error!', sessions{i}); end
end


%% load experiment data
disp('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'learning_data.mat'), 'data'); disp('learning data loaded!')

%% compute new data and append to loaded data
loadOldData = true;
if exist('data', 'var') && loadOldData; data = getExperimentData(sessionInfo, vars, data); else; data = getExperimentData(sessionInfo, 'all'); end
disp('saving data...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'learning_data.mat'), 'data', '-v7.3'); disp('data saved');

%% compute experiment from scratch, in parallel
data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), vars); end
data = cat(2,data{:});
save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'learning_data.mat'), 'data');















