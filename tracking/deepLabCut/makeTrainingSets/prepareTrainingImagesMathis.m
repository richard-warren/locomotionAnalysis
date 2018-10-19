function prepareTrainingImagesMathis(writeDir, trainingData, view, features, targetSize)

% same as prepareTrainingImages but accepts a target image size rather than
% a scaling factor


% settings
fileType = '.png';

% get original image size
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' trainingData(1).session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' trainingData(1).session '\runBot.mp4']);
originalSize = [vidTop.Height+vidBot.Height, vidTop.Width];

% get xy values for all labelled features in all frames
isLabeled = false(length(trainingData), length(features)); % keeps track of whether features were labelled in given frames
positions = nan(2, length(features), length(trainingData)); % stores all features in all frames

for i = 1:length(features)
    
    % extract x and y values
    x = cellfun(@(j) j(1), {trainingData.(features{i})}) * (targetSize(2)/originalSize(2));
    y = cellfun(@(j) j(2), {trainingData.(features{i})}) * (targetSize(1)/originalSize(1));
    positions(1,i,:) = x;
    positions(2,i,:) = y;
    
    realBins = ~isnan(x);
    isLabeled(realBins, i) = 1;

end

% structInds = find(all(isLabeled, 2)' & [trainingData.includeFrame]); % only use frames where all frames are labelled
structInds = find([trainingData.includeFrame]); % only use frames where all features are labelled (some can bbe labelled as occluded)


% create images (only for frames in which everything is labelled
currentSession = '';

for i = 1:length(structInds)
    
    % load new video if reached a new session
    if ~strcmp(currentSession, trainingData(structInds(i)).session)
        currentSession = trainingData(structInds(i)).session;
        
        if strcmp(view, 'top')
            vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runTop.mp4']);
        elseif strcmp(view, 'bot')
            vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
        elseif strcmp(view, 'both')
            vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runTop.mp4']);
            vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
        end
        
        fprintf('processing session: %s\n', currentSession)
    end
    
    % get frame
    currentFrame = trainingData(structInds(i)).frameNum;
    if strcmp(view, 'both')
        frame = cat(1, rgb2gray(read(vid, currentFrame)), rgb2gray(read(vidBot, currentFrame)));
    else
        frame = rgb2gray(read(vid, currentFrame));
    end
    
    % rescale frame
    if ~all(targetSize==originalSize); frame = imresize(frame, targetSize); end
    
    imwrite(frame, [writeDir 'img' num2str(i) fileType])
    
end


% save spreadsheet for each feature, which are used by deepLabCut
for i = 1:length(features)
    
    X = squeeze(positions(1,i,structInds));
    Y = squeeze(positions(2,i,structInds));
    
    % set to 0 all NaN values (this is how deepLabCut expects to see occluded features)
    X(isnan(X)) = 0;
    Y(isnan(Y)) = 0;
    
    featureTable = table(X, Y);
    writetable(featureTable, [writeDir features{i} '.csv'], 'delimiter', ' ')
end

disp('all done!')






