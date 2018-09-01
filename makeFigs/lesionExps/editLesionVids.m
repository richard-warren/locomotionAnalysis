

% settings
dir = 'Z:\RAW\sharing\experimentVids\lesionVids\';
trialPortion = .1;
overWriteVids = false;

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'lesionNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);

brainRegions = unique(sessionInfo.brainRegion);

for i = 1:length(brainRegions)
    
    % make directory
    mkdir([dir brainRegions{i}]);
    
    brainRegionBins = strcmp(sessionInfo.brainRegion, brainRegions{i});
    mice = unique(sessionInfo.mouse(brainRegionBins));
    
    for j = 1:length(mice)
        
        mouseBins = strcmp(sessionInfo.mouse, mice{j});
        sessions = sessionInfo.session(mouseBins);
        mouseDir = [dir brainRegions{i} '\' mice{j} '\'];
        mkdir(mouseDir);
        
        for k = 1:length(sessions)
            sessionBin = strcmp(sessionInfo.session, sessions{k});
            condition = sessionInfo.preOrPost{sessionBin};
            side = sessionInfo.side{sessionBin};
            
            fileName = sprintf('%sday%i-%s-%s-%s.avi', mouseDir, k, condition, side, sessions{k});
            if overWriteVids || ~exist(fileName, 'file') % only overwrite existing vid if overWriteVids is true
                try
                    load([getenv('OBSDATADIR') 'sessions\' sessions{k} '\runAnalyzed.mat'], 'isLightOn');
                    makeVidWisk(fileName, sessions{k}, [-.05 .1], .15, trialPortion, {'OFF', 'ON'}, isLightOn+1);
                catch
                    fprintf('%s: failed to edit video!\n', sessions{k})
                end
            end
        end     
    end
end
disp('all done!')


