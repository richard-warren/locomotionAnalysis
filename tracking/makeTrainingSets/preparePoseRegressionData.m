function preparePoseRegressionData(sessions, totalEgs)

% % temp
% sessions = {'180122_000', '180122_001', '180122_002'};
% totalEgs = 500;

% settings
writeDir = 'C:\Users\rick\Desktop\trainingExamples\poseRegression\';
targetSize = [227 227];

% initializations
if ~exist(writeDir, 'dir'); mkdir(writeDir); end
if ~exist([writeDir '\imgs'], 'dir'); mkdir([writeDir '\imgs']); end
lastSessionInd = 0;
imgInd = 1;
% features = nan(targetSize(1), targetSize(2), 3, totalEgs);
features = table({'img', 'x1','y1','x2','y2','x3','y3','x4','y4'});

% concatinate all labeled data sets
sessionInds = []; % stores the session identity for each saved location
sessionFrameInds = []; % stores the frame inds
locationsAll = [];

for i = 1:length(sessions)
    
    % get labeled locations for single session
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\locationsBotCorrected.mat'], 'locations');
    
    % remove nan entries
    frameInds = find(~isnan(sum(squeeze(locations.locationsCorrected(:,1,:)),2)))'; % only keep locations where all paws have a non-nan entry
    locations = locations.locationsCorrected(frameInds,:,:);
    
    % store
    locationsAll = cat(1, locationsAll, locations);
    sessionInds = [sessionInds i*ones(1,size(locationsAll,1))];
    sessionFrameInds = [sessionFrameInds frameInds];
    
end


% select random frames
locationInds = randperm(size(locationsAll,1), totalEgs);
locationInds = sort(locationInds);

% restructure locations into feature matrix
locations = locationsAll(locationInds,:,:);
locations = reshape(locations, totalEgs, 8);

for i = locationInds
    
    % load new video if you have reached the next session
    if sessionInds(i) ~= lastSessionInd
        fprintf('loading session %s', sessions{sessionInds(i)})
        vid = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{sessionInds(i)} '\runBot.mp4']);
        bg = getBgImage(vid, 1000, 120, 2*10e-4, false);
        load([getenv('OBSDATADIR') 'sessions\' sessions{sessionInds(i)} '\runAnalyzed.mat'], 'obsPixPositions')
        lastSessionInd = sessionInds(i);
    end
    
    % get frame
    frame = rgb2gray(read(vid, sessionFrameInds(i)));
    frame = frame - bg;
    
    % mask obstacle
    if ~isnan(obsPixPositions(sessionFrameInds(i)))
        frame = maskObs(frame, obsPixPositions(sessionFrameInds(i)));
    end
    
    % save image
    img = uint8(imresize(frame, 'outputsize', targetSize));
    img = repmat(img, 1, 1, 3);
    imwrite(img, [writeDir 'imgs\img' num2str(imgInd) '.tif'])
    
%     features(:,:,:,imgInd) = img;
    
    % report progress
    disp(imgInd/totalEgs)
    imgInd = imgInd +  1;
end

save([writeDir 'pawLocations.mat'], 'features', 'locations')


% % sanity check pltting to see ifthings worked correctly
% ind = 100;
% ind = locationInds(ind);
% 
% vid = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{sessionInds(ind)} '\runBot.mp4']);
% frame = rgb2gray(read(vid, sessionFrameInds(ind)));
% close all; figure;
% imshow(frame);
% hold on; scatter(locationsAll(ind,[1 3 5 7]), locationsAll(ind,[2 4 6 8]))








