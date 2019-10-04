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
            if (exist([sessionsDir newSessions{1} '\runTop.mp4'], 'file') && ...
                    exist([sessionsDir newSessions{1} '\runBot.mp4'], 'file'))
                dlcAnalysisSuccessful = false;
                currentTime = clock;
                
                % check if video has new or old cropping
                rootDir = fullfile(getenv('OBSDATADIR'), 'sessions', newSessions{1});
                vidTop = VideoReader(rootDir, 'runTop.mp4');
                vidBot = VideoReader(rootDir, 'runBot.mp4'));
                dims = [vidTop.Height, vidTop.Width, vidBot.Height, vidBot.Width];
                clear vidTop vidBot
                
                % if using original dimensions
                if isequal(dims, [168 396 238 396])
                    fprintf('%s: starting DeepLabCut analysis at %i:%i...\n', newSessions{1}, currentTime(4), currentTime(5))
                    tic; [~,~] = system([dlcPath(1:2) ' && cd ' dlcPath ' && batchDLC.bat ' newSessions{1}]); % first move to correct drive, then execute script
                    
                % if using new dimensions
                elseif isequal(dims, [214 448 238 448])
                    
                    fprintf('%s: cropping videos to match old video dimensions...', newSessions{1});
                    
                    % copy and rename original dimension files
                    copyfile(fullfile(rootDir, 'runBot.mp4'), fullfile(rootDir, 'runBot_originalDimensions.mp4'))
                    copyfile(fullfile(rootDir, 'runTop.mp4'), fullfile(rootDir, 'runTop_originalDimensions.mp4'))
                    
                    % crop so dimensions match old dimensions
                    system(['ffmpeg -y -loglevel panic -r 250 -i ' fullfile(rootDir, 'runTop.mp4') ...
                        ' -filter:v "crop=396:168:68:52" -vb 10M -vcodec mpeg4 ' fullfile(rootDir, 'runTop.mp4')]);
                    system(['ffmpeg -y -loglevel panic -r 250 -i ' fullfile(rootDir, 'runBot.mp4') ...
                        ' -filter:v "crop=396:238:68:0" -vb 10M -vcodec mpeg4 ' fullfile(rootDir, 'runBot.mp4')]);
                    
                    fprintf('%s: starting DeepLabCut analysis at %i:%i...\n', newSessions{1}, currentTime(4), currentTime(5))
                    tic; [~,~] = system([dlcPath(1:2) ' && cd ' dlcPath ' && batchDLC.bat ' newSessions{1}]); % first move to correct drive, then execute script
                    
                % if dimensions not recognized
                else
                    fprintf('%s: WARNING! Unexpected frame dimensions!\n', newSessions{1})
                end
                
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
            checkObsLight(sessions{1});
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