%% general 'scratch pad' for modelling paper

sessions = getEphysSessions();

%% perform various analyses on all sessions

sessions = getEphysSessions;
sessions = sessions(1:33);  % temp

overwrite = true;
tic

% parpool('local', 2);  % set number of workers
parfor i = 1:length(sessions)
    folder = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
    try
        % format ephys data
%         if overwrite || ~exist(fullfile(folder, 'neuralData.mat'), 'file')
%             formatEphysData(sessions{i})
%         end
        
        % predictors
%         if overwrite || ~exist(fullfile(folder, 'modelling', 'predictors.mat'), 'file')
%             getPredictors(sessions{i})
%         end
        
        % neural responses
        if overwrite || ~exist(fullfile(folder, 'modelling', 'responses.mat'), 'file')
            getNeuralResponses(sessions{i})
        end
        
        % feature importance
        if overwrite || ~exist(fullfile(folder, 'modelling', 'importance.mat'), 'file')
            getFeatureImportance(sessions{i})
        end
        
        % plot neural responses
        plotNeuralResponses(sessions{i}, 'visible', false, 'showImportance', true)
        
        % close figures
%         close all
        
        fprintf('----------- %s complete! -----------\n', sessions{i})
    
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












