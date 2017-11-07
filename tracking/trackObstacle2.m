% get col sums for many obstacle frames

% load data
session = '171023_002';
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
vid = VideoReader([dataDir session '\runBot.mp4']);
load([dataDir session '\runAnalyzed.mat'], 'obsOnTimes', 'obsOffTimes',...
                                           'obsPositions', 'obsTimes',...
                                           'frameTimeStamps');
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

% settings
frameEdges = [.336 .444];
statFrames = 100;
minObsPos = 80; % pixels... tracking is unreliable to the left because of wheel corner... this is temporary hack bro!


%% !!! statistically determine pixel threshold
% for each obs trial, select n frames, compute col sums
% aplly k means to all sums and use highest mean as threshold for detecting object in column


colSumsStats = nan(statFrames, vid.Width);

for i = 1:statFrames
    
    % select random trial and obs position
    trial = round(rand(1) * (length(obsOnTimes)-1))+1;
    framePosit = rand(1) * diff(frameEdges) + min(frameEdges);
    
    % find frame index
    isInTrial = obsTimes>obsOnTimes(trial) & obsTimes<obsOffTimes(trial);              % binary vector of obsTimes that are in trial
    frameTime = obsTimes( find(isInTrial & obsPositions>=framePosit, 1, 'first') );    % time at which obstacle reaches position
    frameInd = find(frameTimeStamps>=frameTime, 1, 'first');                           % frame index at which obstacles reaches position
    
    % get frame
    frame = rgb2gray(read(vid, frameInd));
    colSumsStats(i,:) = sum(frame,1) / size(frame,1);
    disp(i)
    
end


colSumsStats = colSumsStats(:);
guesses = kmeans(colSumsStats, 2);
pixThresh = mean(colSumsStats(guesses==2));

%% !!! find obstacles
% go through all frames in obsOn epochs, getting column sums and thresholding
% eliminate non-adjacent colums
% get mean of remaining columns

close all; figure;
im = imshow(rgb2gray(read(vid,1)));
pimpFig
obsLine = line([0 0], [1 size(frame,1)], 'Visible', 'off', 'lineWidth', 15, 'color', 'red');

% colSums = nan(vid.NumberOfFrames, vid.Width);
obsPosits = nan(1,length(frameTimeStamps));
delay = .001;

for i = 1:length(obsOnTimes)

    frameInds = find(frameTimeStamps>obsOnTimes(i) & frameTimeStamps<obsOffTimes(i));

    for j = 1:length(frameInds)
        
        frame = rgb2gray(read(vid, frameInds(j)));
        
        colSum = sum(frame,1) / size(frame,1);
        colSum = colSum > pixThresh;
        inds = find(colSum);
        
        pos = median(inds);
        if pos < minObsPos; pos = 0; end
        obsPosits(frameInds(j)) = pos;
                     
        
        % update image preview
        set(im, 'CData', frame);
        
        % update lines
        if ~isnan(pos)
            set(obsLine, 'XData', [pos pos], 'Visible', 'on');
            delay = .03;
        else
            set(obsLine, 'Visible', 'off');
            delay = .001;
        end
        
        pause(delay);
        
    end
end




%% verification
% ???


%% plot
% prepare figure
close all; figure;

hist(colSumsStats(:), size(frame,2)/2)

pimpFig;



