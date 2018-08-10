

% settings
dir = 'Z:\RAW\sharing\experimentVids\whiskerTrimVids\';
trialPortion = .1;
overWriteVids = false;

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);


mice = unique(sessionInfo.mouse);
    
for j = 1:length(mice)

    mouseBins = strcmp(sessionInfo.mouse, mice{j});
    sessions = sessionInfo.session(mouseBins);
    mouseDir = [dir mice{j} '\'];
    mkdir(mouseDir);

    for k = 1:length(sessions)
        sessionBin = strcmp(sessionInfo.session, sessions{k});
        condition = sessionInfo.preOrPost{sessionBin};

        fileName = sprintf('%sday%i-%s-%s.avi', mouseDir, k, condition, sessions{k});
        if overWriteVids || ~exist(fileName, 'file') % only overwrite existing vid if overWriteVids is true
            load([getenv('OBSDATADIR') 'sessions\' sessions{k} '\runAnalyzed.mat'], 'isLightOn');
            makeVidWisk(fileName, sessions{k}, [-.05 .1], .15, trialPortion, {'OFF', 'ON'}, isLightOn+1);
        end
    end     
end

disp('all done!')


