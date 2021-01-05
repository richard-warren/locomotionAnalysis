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

        % design matrices
        filename = fullfile(folder, 'designMatrices', [sessions{i} '_designMatrix.mat']);
        if overwrite || ~exist(filename, 'file')
            makeDesignMatrix(sessions{i}, 'timeDegrees', 3, 'outputFileName', filename);
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

%% todo (automatically copy files from engram to local, skipping files that are already there...)


