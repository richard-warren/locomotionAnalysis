%% general 'scratch pad' for modelling paper

sessions = getEphysSessions();

%% perform various analyses on all sessions

sessions = getEphysSessions;
sessions = sessions(1:33);  % temp

overwrite = true;
tic

% parpool('local', 4);  % set number of workers
for i = 1:length(sessions)
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
        
%         % neural responses
%         if overwrite || ~exist(fullfile(folder, 'modelling', 'responses.mat'), 'file')
%             getNeuralResponses(sessions{i})
%             plotNeuralResponses(sessions{i}, 'visible', false)
%         end
%         
%         % feature importance
        if overwrite || ~exist(fullfile(folder, 'modelling', 'importance.mat'), 'file')
            getFeatureImportance(sessions{i})
            plotNeuralResponses(sessions{i}, 'visible', false)
        end
%         
%         % plot neural responses
%         plotNeuralResponses(sessions{i}, 'visible', false, 'showImportance', false)
%         
%         % close figures
%         close all
    catch
        fprintf('%s: problem with analysis\n', sessions{i})
    end
end
toc


%% look into MI

session = '200703_000';  % reward_all
unitLow = 170;
unitHigh = 275;

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'responses.mat'), 'responses');
responseLow = responses([responses.cell]==unitLow).responses{'reward_all', 'response'}{1};
responseHigh = responses([responses.cell]==unitHigh).responses{'reward_all', 'response'}{1};
%%

contact = responses([responses.cell]==261).responses{'paw1LH_contact_ventral', 'response'}{1};













