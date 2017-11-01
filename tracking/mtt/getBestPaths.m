
% load tracking data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\tracked.mat', 'locations')

% user settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\botTest.mp4';

unaryWeight = 1;
pairwiseWeight = 0;
occludedWeight = 0;
occlusionGridSpacing = 15;

maxVelocity = 20;
pawsToShow = 1;
vidDelay = .05;

% initializations
frameHeight = 242; % !!! hacky temp
frameWidth = 398;
gridX = 1 : occlusionGridSpacing : frameWidth;
gridY = 1 : occlusionGridSpacing : frameHeight;
numOccluded = length(gridX) * length(gridY);


% compute unary potentials
unary = cell(1,length(locations));

for i = 1:length(locations)
    
    frameUnaries = nan(4, length(locations(i).x) + numOccluded);
    
    frameUnaries(1,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 0, numOccluded, unaryWeight); % LH
    frameUnaries(2,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 1, numOccluded, unaryWeight); % LF
    frameUnaries(3,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 0, numOccluded, unaryWeight); % RH
    frameUnaries(4,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 1, numOccluded, unaryWeight); % RF
    
    
    unary{i} = frameUnaries;
    
end



%% compute pairwise potentials
pairwise = cell(1,length(locations)-1);

for i = 1:length(locations)-1
    
    pairwise{i} = getPairwisePotentials([locations(i+1).x, locations(i+1).y], [locations(i).x, locations(i).y], maxVelocity, pairwiseWeight, occludedWeight);
    
end




%% find most probable paths!

objectNum = size(unary{1}, 1);
labels = nan(length(unary), objectNum);
nodeScores = cell(1, length(unary));
backPointers = cell(1, length(unary));

% initializations

nodeScores{1} = unary{1};


for i = 2:length(unary)
    
    nodeScores{i}   = nan(objectNum, size(unary{i},2));
    backPointers{i} = nan(objectNum, size(unary{i},2));
    
    for j = 1:objectNum
        
        previousScores = nodeScores{i-1}(j,:);
        currentUnary   = unary{i}(j,:)';
        
        allTransitionScores = repmat(previousScores, length(currentUnary), 1) .* ...
                              (pairwise{i-1} + repmat(currentUnary, 1, length(previousScores)));
        
        
        [nodeScores{i}(j,:), backPointers{i}(j,:)] = max(allTransitionScores, [], 2);
%         [~, labels(i,j)] = max(nodeScores{i}(j,:));
%         disp(backPointers{i}(j,:)); keyboard
%         keyboard
        
    end
end


% back trace

[pathScores, labels(end,:)] = max(nodeScores{end}, [], 2);
% 
for i = fliplr(1:length(unary)-1)
    for j = 1:objectNum
        
        labels(i,j) = backPointers{i+1}(labels(i+1,j));

    end
end




% VISUALIZE TRACKING


% initializations
startFrame = 1;
vid = VideoReader(vidFile);
sampleFrame = rgb2gray(read(vid,startFrame));
totalFrames = vid.NumberOfFrames;
cmap = winter(length(pawsToShow));

% prepare figure
close all; figure('position', [567 383 698 400], 'color', 'black'); colormap gray


rawIm = image(sampleFrame, 'CDataMapping', 'scaled');
rawAxis = gca;
set(rawAxis, 'visible', 'off')
hold on;
scatterPts =    scatter(rawAxis, zeros(1,length(pawsToShow)), zeros(1,length(pawsToShow)), 200, cmap, 'filled'); hold on
scatterPtsAll = scatter(rawAxis, 0, 0, 200, [1 1 1]);


for i = startFrame:totalFrames
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % update figure
    set(rawIm, 'CData', frame);
    ind = backPointers{i+1}(pawsToShow,1);
    
    xs = [locations(i).x; 0]; % add zero for occluded state
    ys = [locations(i).y; 0]; % add zero for occluded state
    
    set(scatterPts, 'XData', xs(backPointers{i+1}(pawsToShow,1)), 'YData', ys(backPointers{i+1}(pawsToShow,1)));
    set(scatterPtsAll, 'XData', locations(i).x, 'YData', locations(i).y);
    
    % pause to reflcet on the little things...
    pause(vidDelay);
end










