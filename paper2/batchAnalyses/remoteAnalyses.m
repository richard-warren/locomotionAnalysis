%% perform analyses that require access to engram, which is slow from home ethernet

sessions = getEphysSessions();
% sessions = sessions(1:33);  % temp
overwrite = true;


tic
parfor i = 1:length(sessions)
    try
        % format ephys data
        filename = fullfile(getenv('OBSDATADIR'), 'data_transfer', 'neuralData', [sessions{i} '_neuralData.mat']);
        if overwrite || ~exist(filename, 'file')
            formatEphysData(sessions{i}, 'outputFileName', filename, 'kernel', 'gauss', 'kernelSig', .02, 'plot', false)
        end
        
        % predictors
        filename = fullfile(getenv('OBSDATADIR'), 'data_transfer', 'predictors', [sessions{i} '_predictors.mat']);
%         if overwrite || ~exist(filename, 'file')
%             getPredictors(sessions{i}, 'outputFileName', filename, 'plot', true, 'visible', 'off')
%         end
            
    catch exception
        fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
    end
end
toc