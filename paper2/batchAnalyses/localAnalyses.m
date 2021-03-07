% perform analyses that don't require access to engram

sessions = getEphysSessions();
% sessions = sessions(1:33);  % temp
% sessions = {'181020_001'};

overwrite = true;


tic
parfor i = 1:length(sessions)
    folder = fullfile(getenv('SSD'), 'paper2', 'modelling');
    try
        % neural responses
        filename = fullfile(folder, 'responses', [sessions{i} '_responses.mat']);
        if overwrite || ~exist(filename, 'file')
            getNeuralResponses(sessions{i})
        end
        
        % plot neural responses
        existingPlots = dir(['E:\lab_files\paper2\plots\responses\' sessions{i} '*']);
        if overwrite || isempty(existingPlots)
            plotNeuralResponses(sessions{i}, 'visible', false, 'showImportance', false)
        end

    catch exception
        fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
    end
end
fprintf('finished in %.1f minutes\n', toc/60)


%% histo

% settings
overwrite = false;


load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), 'ephysHistoTable')
mice = unique(ephysHistoTable.mouseID);

for i = 1:length(mice)
    filename = ['E:\lab_files\paper2\histo\registration\' mice{i} '_registration.mat'];
    if ~exist(filename, 'file') || overwrite
        registerBrain(mice{i});
        close all
    end
end

%% copy runAnalyzed.mat from engram to local SSD

% copy runAnalyzed.mat from engram to local folder for faster reading

% settings
overwrite = false;


data = getUnitInfo();
sessions = unique(data.session);

for i = 1:length(sessions)
    engramfile = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat');
    localfile = fullfile(getenv('SSD'), 'paper2', 'modelling', 'runAnalyzed', ...
        [sessions{i} '_runAnalyzed.mat']);
    
    if ~exist(localfile, 'file')
        fprintf('(%3i/%i) copying %s from engram to local SSD\n', ...
            i, length(sessions), sessions{i})
        copyfile(engramfile, localfile);
    end
end
disp('all done!')




%% todo (automatically copy files from engram to local, skipping files that are already there...)


