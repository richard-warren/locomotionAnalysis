% perform analyses that don't require access to engram

sessions = getEphysSessions();
% sessions = sessions(1:33);  % temp
% sessions = {'181020_001'};

overwrite = false;


tic
for i = 1:length(sessions)
    folder = fullfile(getenv('SSD'), 'paper2', 'modelling');
%     try
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

        % design matrices (default)
        filename = fullfile(folder, 'designMatrices', [sessions{i} '_designMatrix.mat']);
        if overwrite || ~exist(filename, 'file')
            makeDesignMatrix(sessions{i}, 'timeDegrees', 3, 'outputFileName', filename);
        end
        
        % design matrices (epochs)
        filename = fullfile(folder, 'designMatrices', 'epochs', [sessions{i} '_designMatrix.mat']);
        if overwrite || ~exist(filename, 'file')
            makeDesignMatrix(sessions{i}, 'timeDegrees', 3, 'outputFileName', filename, ...
                'predictorSpreadsheet', 'C:\Users\richa\Desktop\github\locomotionAnalysis\paper2\glm\epoch_predictorSettings.xlsx');
        end
        
%     catch exception
%         fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
%     end
end
fprintf('finished in %.1f minutes\n', toc/60)


%% histo

load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), 'ephysHistoTable')
mice = unique(ephysHistoTable.mouseID);

overwrite = false;

for i = 1:length(mice)
    filename = ['E:\lab_files\paper2\histo\registration\' mice{i} '_registration.mat'];
    if ~exist(filename, 'file') || overwrite
        registerBrain(mice{i});
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





%% todo (automatically copy files from engram to local, skipping files that are already there...)


