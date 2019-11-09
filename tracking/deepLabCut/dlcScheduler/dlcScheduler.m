% todo: make sure it doesnt get caught trying to reanalyze problematic
% sessions // have separate memory for whether dlc and spike analyses are
% done


% settings
checkFrequency = 60; % seconds
sessionsDir = [getenv('OBSDATADIR') 'sessions\'];
schedulerDir = [getenv('GITDIR') 'locomotionAnalysis\tracking\deepLabCut\dlcScheduler\'];
dlcPath = [getenv('GITDIR') 'DeepLabCutBatch'];

% initializations
load([schedulerDir 'analyzedSessions.mat'], 'analyzedSessions');


% start scheduling
fprintf('waiting for new files to analyze...')
while true
    
    % find new sessions
    sessions = dir(sessionsDir);
    sessions = {sessions([false false sessions(3:end).isdir]).name};
    newSessions = sessions(~ismember(sessions, analyzedSessions));
    validSessionNameBins = cellfun(@(x) length(x)==10 && any(strfind(x,'_')), newSessions); % this is hacky, but checks that the length of the file name is 10 and that it also contains an underscore, which conforms to naming concentions for sessions
    newSessions = newSessions(validSessionNameBins);
    
    if ~isempty(newSessions)
        
        % wait for files to be transferred
        fprintf('\n\n---------%s: deteceted new session---------\n', newSessions{1});
        pause(120); % a hack to make sure files transfer completely before starting analysis... sad hack face:   :(
        
        
        % DeepLabCut analysis
        try
            rootDir = fullfile(getenv('OBSDATADIR'), 'sessions', newSessions{1});
            
            if (exist(fullfile(rootDir, 'run.mp4'), 'file') || exist(fullfile(rootDir, 'runTop.mp4'), 'file'))  % check if video was recorded for this session
                
                dlcAnalysisSuccessful = false;
                currentTime = clock;
                if ~exist(fullfile(rootDir, 'run.mp4')); concatTopBotVids(newSessions{1}); end  % old sessions were recorded with separate top and bot views, which need to be concatenated
                
                % if using new dimensions crop to old dimensions
                vid = VideoReader(fullfile(rootDir, 'run.mp4'));
                dims = [vid.Height, vid.Width];
                clear vid
                
                if ~isequal(dims, [406 396])    
                    fprintf('%s: cropping videos to match old video dimensions...\n', newSessions{1});
                    copyfile(fullfile(rootDir, 'run.mp4'), fullfile(rootDir, 'run_originalDimensions.mp4'))  % copy and rename original dimension files
                    system(['ffmpeg -y -loglevel panic -r 250 -i ' fullfile(rootDir, 'run_originalDimensions.mp4') ...
                        ' -filter:v "crop=396:406:44:52" -vb 10M -vcodec mpeg4 ' fullfile(rootDir, 'run.mp4')]);
                end
                
                fprintf('%s: starting DeepLabCut analysis at %i:%i...\n', newSessions{1}, currentTime(4), currentTime(5))
                    tic; [~,~] = system([dlcPath(1:2) ' && cd ' dlcPath ' && batchDLC.bat ' newSessions{1}]); % first move to correct drive, then execute script
                fprintf('%s: DeepLabCut analysis finished in %.1f hours\n', newSessions{1}, toc/60/60)
                dlcAnalysisSuccessful = true;
            else
                dlcAnalysisSuccessful = true;
            end
        catch
            fprintf('%s: problem with DeepLabCut analysis!\n', newSessions{1})
        end
        
        
        % spike analysis
        try 
            spikeAnalysis2(newSessions{1});
            showWiskContactFrames(newSessions{1});
%             checkObsLight(sessions{1});
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