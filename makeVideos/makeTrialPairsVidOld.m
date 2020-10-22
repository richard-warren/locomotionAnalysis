function makeTrialPairsVidOld(filename, sessions, trials, trialText)

% makes vid with different trials show at same time as different rows in
% same video, with all pausing at moment of wisk contact // number of rows
% in video = number of rows in sessions and trials // both session and
% trials are n rows by n trials matrices with name of session and trial
% number in each entry

% settings
prePostTime = [-.4 .4]; % (s) time before and after wisk contact to show
contrastLims = [.2 1]; % pixels at these proportional values are mapped to 0 and 255
playBackSpeed = .1;
wiskPause = 1; % pause video for this many seconds at wisk contact
showTrialsTwice = true;

% initializations
numTrials = size(trials,2);
uniqueSessions = unique(sessions(:));
rows = size(sessions,1);

[sessionsData, vidTops, vidBots] = deal(cell(1, length(uniqueSessions)));
for i = 1:length(uniqueSessions)
    sessionsData{i} = load([getenv('OBSDATADIR') 'sessions\' uniqueSessions{i} '\runAnalyzed.mat'], ... 
        'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'wiskContactTimes');
    vidTops{i} = VideoReader([getenv('OBSDATADIR') 'sessions\' uniqueSessions{i} '\runTop.mp4']);
    vidBots{i} = VideoReader([getenv('OBSDATADIR') 'sessions\' uniqueSessions{i} '\runBot.mp4']);
end

fps = median(cellfun(@(x) x.FrameRate, cat(2,vidTops,vidBots)));
times = prePostTime(1) : (1/fps) : prePostTime(2);
contactInd = find(times>=0, 1, 'first');
targetFps = round(fps * playBackSpeed);
vidWriter = VideoWriter(filename, 'MPEG-4');
set(vidWriter, 'FrameRate', targetFps);
open(vidWriter);


% iterate through trials
for i = repelem(1:numTrials,showTrialsTwice+1)
    
%     disp(i/numTrials)
    
    % get frame inds for each row of video
    trialInds = nan(rows, length(times));
    for row = 1:rows
        sesBin = strcmp(uniqueSessions, sessions(row,i));
        contactTime = sessionsData{sesBin}.wiskContactTimes(trials(row,i));
        trialTimes = times + contactTime; % desired times
        trialInds(row,:) = knnsearch(sessionsData{sesBin}.frameTimeStamps, trialTimes')'; % inds of frames with times closest to desired times
    end
    
    for j = 1:size(trialInds, 2)
        
        % get frame
        frames = cell(1,rows);
        for row = 1:rows
            sesBin = strcmp(uniqueSessions, sessions(row,i));
            frame = cat(1, rgb2gray(read(vidTops{sesBin}, trialInds(row,j))), rgb2gray(read(vidBots{sesBin},trialInds(row,j))));
            frame = insertText(frame, [size(frame,2) size(frame,1)], trialText{row,i},...
                                   'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
            frames{row} = frame;
        end
        frame = cat(1, frames{:});
        frame = imadjust(frame, contrastLims, [0 1]);
        
        % write to video
        writeVideo(vidWriter, frame);
        
        % pause video at momemt of contact
        if j==contactInd
            for k = 1:(wiskPause*targetFps); writeVideo(vidWriter, frame); end
        end
    end
end

close(vidWriter)
disp('all done!')