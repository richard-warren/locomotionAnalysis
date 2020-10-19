% perform analyses that don't require access to engram

sessions = getEphysSessions();
% sessions = sessions(1:33);  % temp
% sessions = {'181004_003'};

overwrite = true;


tic
parfor i = 1:length(sessions)
%     try
        % neural responses
%         filename = fullfile(folder, 'responses', [sessions{i} '_responses.mat']);
%         if overwrite || ~exist(filename, 'file')
%             getNeuralResponses(sessions{i})
%         end
        
        % plot neural responses (local)
%         plotNeuralResponses(sessions{i}, 'visible', false, 'showImportance', false)

        % design matrices (requires predictors only...)
        filename = fullfile(folder, 'designMatrices', [sessions{i} '_designMatrix.mat']);
        if overwrite || ~exist(filename, 'file')
            makeDesignMatrix(sessions{i}, 'timeDegrees', 3, 'outputFileName', filename);
        end
        
%     catch exception
%         fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
%     end
end
toc


%% todo (automatically copy files from engram to local, skipping files that are already there...)


