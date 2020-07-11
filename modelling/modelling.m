%% general 'scratch pad' for modelling paper

sessions = getEphysSessions();

%% perform various analyses on all sessions

overwrite = true;

for i = 1:length(sessions)
    folder = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
    try
        % format ephys data
%         if overwrite || ~exist(fullfile(folder, 'neuralData.mat'), 'file')
%             formatEphysData(sessions{i})
%         end
        
        % predictors
        if overwrite || ~exist(fullfile(folder, 'modelling', 'predictors.mat'), 'file')
            getPredictors(sessions{i})
        end
        
        % neural responses
%         if overwrite || ~exist(fullfile(folder, 'modelling', 'responses.mat'), 'file')
%             getNeuralResponses(sessions{i})
%         end
        
        % feature importance
%         if overwrite || ~exist(fullfile(folder, 'modelling', 'importance.mat'), 'file')
%             getFeatureImportance(sessions{i})
%         end
        
        % plot neural responses
%         plotNeuralResponses(sessions{i}, 'visible', false)
        
        % close figures
%         close all
    catch
        fprintf('%s: problem with analysis\n', sessions{i})
    end
end


%% create missing .dat files

for i = 1:length(ephysSessions)
    try
        ephysFolder = dir(fullfile(getenv('OBSDATADIR'), 'sessions', ephysSessions{i}, 'ephys_*'));
        datFile = dir(fullfile(ephysFolder.folder, ephysFolder.name, '*.dat'));
        if isempty(datFile)
            packContFiles(ephysSessions{i});
        end
    catch
        fprintf('%s: problem with analysis\n', ephysSessions{i})
    end
end


