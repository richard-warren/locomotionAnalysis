%% general 'scratch pad' for modelling paper

sessions = getEphysSessions();

%% perform various analyses on all sessions

sessions = getEphysSessions;
sessions = sessions(1:33);  % temp

overwrite = false;
tic

% parpool('local', 4);  % set number of workers
parfor i = 1:length(sessions)
    folderSes = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
    folder = fullfile(getenv('SSD'), 'modelling');
    
    try
        % format ephys data
        if overwrite || ~exist(fullfile(folderSes, 'neuralData.mat'), 'file')
            formatEphysData(sessions{i})
        end
        
        % predictors
        if overwrite || ~exist(fullfile(folder, 'predictors', [sessions{i} '_predictors.mat']), 'file')
            getPredictors(sessions{i}, 'plotPredictors', true, 'visible', 'off')
        end
        
        % neural responses
        if overwrite || ~exist(fullfile(folder, 'responses', [sessions{i} '_responses.mat']), 'file')
            getNeuralResponses(sessions{i})
        end
        
        % feature importance
        if overwrite || ~exist(fullfile(folder, 'importance', [sessions{i} '_importance.mat']), 'file')
            getFeatureImportance(sessions{i})
        end
        
        % plot neural responses
%         plotNeuralResponses(sessions{i}, 'visible', false, 'showImportance', true)
    
    catch exception
        fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
    end
end
toc


%% look into MI

session = '180922_001';  % reward_all
unitLow = 84;

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'responses.mat'), 'responses');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), 'predictors');
response = responses{'reward_all', 'response'}{1};












