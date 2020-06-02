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


%% get rid of cropped versions for ephys experiments

% rename _orig files to overwrite cropped // if runTop and runBot, concat to run // get rid of runTop and runBot?


%% to clear up disk space we could:


% get rid of runTop, runBot when run exists

% get rid of originalDimensions for old sessions OR ...

% get rid of cropped for old sessions and reanalyze: could potentially
% analyze on cropped vids with old network, then shift the coordinates to
% accomodate the cropping, and throw away cropped vid, performing the rest
% of the analysis on the uncropped video // alternatively, see if old DLC
% can handle uncropped vids, and reanalyze like that...







