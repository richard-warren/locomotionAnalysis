function autoAnalyze(varargin)

% this script keeps track of which sessions have been analyzed using
% analyzeSession() // it automatically checks the sessions folder and
% analyzes sessions that have not yet been analyzed


% settings
s.plotDiagnostics = true;
s.checkFrequency = 60; % (seconds)

% initializations
analyzedFile = fullfile(getenv('GITDIR'), 'locomotionAnalysis', 'tracking', 'autoAnalyze', 'analyzedSessions.mat');
load(analyzedFile, 'analyzedSessions');


fprintf('\nwaiting for new sessions to analyze...')
failedSessions = {};
while true
    
    % find new sessions
    sessions = dir(fullfile(getenv('OBSDATADIR'), 'sessions'));
    sessions = {sessions([sessions.isdir]).name};
    sessions = sessions(cellfun(@(x) ~isempty(regexp(x, '\d\d\d\d\d\d_\d\d\d', 'once')), sessions));  % only keep valid session names
    newSessions = sessions(~ismember(sessions, [analyzedSessions, failedSessions]));
    
    if ~isempty(newSessions)
        for i = 1:length(newSessions)
            fprintf('\n\n---------%s: deteceted new session---------\n', newSessions{i});
            pause(120);  % wait to make sure data are fully transfered to engram
            
            try
                analyzeSession(newSessions{i}, 'plotDiagnostics', s.plotDiagnostics)
                analyzedSessions{end+1} = newSessions{i};
                save(analyzedFile, 'analyzedSessions');
            catch
                fprintf('%s: WARNING! autoAnalyze failed!\n', newSessions{i})
                failedSessions{end+1} = newSessions{i};
            end
        end
        fprintf('waiting for new files to analyze...')
    end
    
    pause(s.checkFrequency)
end

