%% find out which sessions have 'originalDimensions' files

files = dir(fullfile(getenv('OBSDATADIR'), 'sessions'));
sessions = {files([files.isdir]).name};
origSessions = {};

fileExists = false(1,length(sessions));
fprintf('\n\n--------------------looking for originalDimensions--------------------\n')
for i = 1:length(sessions)
    dirSub = dir(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, '*.mp4'));
    bins = contains({dirSub.name}, 'originalDimensions');
    if any(bins)
        fprintf('%s: ', sessions{i})
        fprintf('%s ', dirSub(bins).name)
        fprintf('\n')
        origSessions{end+1} = sessions{i};
    end
end
disp('all done!')

%% find sessions with ephys folders

files = dir(fullfile(getenv('OBSDATADIR'), 'sessions'));
sessions = {files([files.isdir]).name};
ephysSessions = {};

hasEphysFolder = false(1,length(sessions));
fprintf('\n\n--------------------looking for ephys folders--------------------\n')
for i = 1:length(sessions)
    dirSub = dir(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}));
    dirSub = dirSub([dirSub.isdir]);
    bins = contains({dirSub.name}, 'ephys_');
    if any(bins)
        fprintf('%s: ', sessions{i})
        fprintf('%s ', dirSub(bins).name)
        fprintf('\n')
        ephysSessions{end+1} = sessions{i};
    end
end
disp('all done!')


% %% get rid of cropped views and concat top and bot FOR EPHYS SESSIONS ONLY
% % note: this does not delete the un-concatenated versions, which could optionally be done later to save disk space
% 
% for i = 1:length(ephysSessions)
%     
%     folder = fullfile(getenv('OBSDATADIR'), 'sessions', ephysSessions{i});
%     dirSub = dir(fullfile(folder, '*.mp4'));
%     
%     % rename originalDimensions
%     origInds = find(contains({dirSub.name}, '_originalDimensions'));
%     if ~isempty(origInds)
%         fprintf('%s: renaming files ', ephysSessions{i})
%         for j = 1:length(origInds)
%             fprintf('%s ', dirSub(origInds(j)).name)
%             movefile(fullfile(folder, dirSub(origInds(j)).name), ...
%                      fullfile(folder, erase(dirSub(origInds(j)).name, '_originalDimensions')));
%         end
%         fprintf('\n')
%         
%         % concatenate views if originally recorded un-concatenated
%         if exist(fullfile(folder, 'runTop.mp4')); concatTopBotVids(ephysSessions{i}); end
%         fprintf('\n')
%     end
% end

%% reanalyze ephys sessions

for i = 1%:length(ephysSessions)
    fprintf('\n___________%i/%i___________\n', i, length(ephysSessions))
    try
        analyzeSession(ephysSessions{i}, ...
            'overwriteVars', 'all', ...
            'verbose', true, ...
            'superVerbose', false, ...
            'rerunRunNetwork', true, ...
            'rerunWiskNetwork', true, ...
            'rerunPawContactNetwork', true, ...
            'rerunWiskContactNetwork', true);
    catch
        fprintf('%s: problem with analysis!\n', sessions{i})
    end
end
disp('all done!')


%% to clear up disk space we could:


% get rid of runTop, runBot when run exists

% get rid of originalDimensions for old sessions OR ...

% get rid of cropped for old sessions and reanalyze: could potentially
% analyze on cropped vids with old network, then shift the coordinates to
% accomodate the cropping, and throw away cropped vid, performing the rest
% of the analysis on the uncropped video // alternatively, see if old DLC
% can handle uncropped vids, and reanalyze like that...







