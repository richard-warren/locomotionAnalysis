function prepareTrainingImages(writeDir, trainingData, features)


% settings
fileType = '.png';


% get xy values for all labelled features in all frames
isLabeled = false(length(trainingData), length(features)); % keeps track of whether features were labelled in given frames
positions = nan(2, length(features), length(trainingData)); % stores all features in all frames

for i = 1:length(features)
    
    % extract x and y values
    x = cellfun(@(j) j(1), {trainingData.(features{i})});
    realBins = ~isnan(x);
    isLabeled(realBins, i) = 1;
    y = nan(size(x));
    y(realBins) = cellfun(@(j) j(2), {trainingData(realBins).(features{i})});
    
    positions(1,i,:) = x;
    positions(2,i,:) = y;
    
end
structInds = find(all(isLabeled, 2)); % only use frames where all frames are labelled



% create images (only for frames in which everything is labelled
currentSession = trainingData(structInds(1)).session;
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);

for i = 1:length(structInds)
    
    if ~strcmp(currentSession, trainingData(structInds(i)).session)
        currentSession = trainingData(structInds(i)).session;
        vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
    end
    frame = rgb2gray(read(vid, trainingData(structInds(i)).frameNum));
    imwrite(frame, [writeDir 'img' num2str(i) fileType])
    
end


% save spreadsheet for each feature, which are used by deepLabCut
for i = 1:length(features)
    
    X = squeeze(positions(1,i,structInds));
    Y = squeeze(positions(2,i,structInds));
    featureTable = table(X, Y);
    
    writetable(featureTable, [writeDir features{i} '.csv'], 'delimiter', ' ')
end

disp('all done!')






