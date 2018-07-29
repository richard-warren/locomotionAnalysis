% to do: make sure it doesnt get caught trying to reanalyze problematic
% sessions // have separate memory for whether dlc and spike analyses are
% done


% settings
checkFrequency = 60; % seconds
sessionsDir = [getenv('OBSDATADIR') 'sessions\'];
schedulerDir = [getenv('GITDIR') 'locomotionAnalysis\tracking\deepLabCut\dlcScheduler\'];
dlcPath = 'C:\Users\rick\Desktop\github\DeepLabCut';

% initializations
load([schedulerDir 'analyzedSessions.mat'], 'analyzedSessions');


% start scheduling
fprintf('waiting for new files to analyze...')
while true
    
    % find new sessions
    sessions = dir(sessionsDir);
    sessions = {sessions([false false sessions(3:end).isdir]).name};
    newSessions = sessions(~ismember(sessions, analyzedSessions));
    validSessionNameBins = cellfun(@(x) length(x)==10 && any(strfind(x,'_')), newSessions); % this is hacky, but checks that the length of the file name is 10 and thhat it also contains an underscore, which conforms to naming concentions for sessions
    newSessions = newSessions(validSessionNameBins);
    
    if ~isempty(newSessions)
        
        % wait for files to be transferred
        while ~(exist([sessionsDir newSessions{1} '\runTop.mp4'], 'file') && ...
                exist([sessionsDir newSessions{1} '\runTop.mp4'], 'file'));
            pause(5);
        end
        fprintf('\n%s: deteceted new session\n', newSessions{1});
        
        
        % DeepLabCut analysis
        try
            dlcAnalysisSuccessful = false;
            currentTime = clock;
            fprintf('%s: starting DeepLabCut analysis at %i:%i...\n', newSessions{1}, currentTime(4), currentTime(5))
            tic; [~,~] = system(['cd ' dlcPath ' && batchDLC.bat ' newSessions{1}]);
            fprintf('%s: DeepLabCut analysis finished in %.1f hours\n', newSessions{1}, toc/60/60)
            dlcAnalysisSuccessful = true;
        catch
            fprintf('%s: problem with DeepLabCut analysis!\n', newSessions{1})
        end
        
        
        % spike analysis
        disp('starting to analyze sessions...')
        try
            spikeAnalysis2(sessions{i});
        catch
            fprintf('%s: problem with spike analysis!\n', newSessions{1})
        end
        
        
        % save that session has been analyzed
        if dlcAnalysisSuccessful
            analyzedSessions{end+1} = newSessions{1};
            save([schedulerDir 'analyzedSessions.mat'], 'analyzedSessions');
            fprintf('%s: dlcScheduler analysis completed\n', newSessions{1});
            newSessions = {};
        end
        fprintf('waiting for new files to analyze...')
    end
    
    pause(checkFrequency)
end