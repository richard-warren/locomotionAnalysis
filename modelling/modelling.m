%% general 'scratch pad' for modelling paper

sessions = getEphysSessions();

%% perform various analyses on all sessions

sessions = getEphysSessions;
sessions = sessions(1:33);  % temp

overwrite = true;
tic

parpool('local', 4);  % set number of workers
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


%% test PSTH function

session = '181030_000';
unit = 132;
predictor = 'paw1LH_stride';

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), 'predictors');
events = table2array(predictors(predictor, 'data')); events = events{1};
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), 'spkTimes', 'unit_ids');
spkTimes = spkTimes{unit_ids==unit};

%%
plotPSTH(spkTimes, events, 'removeNoSpikeTrials', true, ...
    'eventLims', [-1 1], 'epochLims', [0 1], 'xlabel', predictor, 'epochDurationLims', [20 80]);





