% to try: ***way of regularizing total values... // smarter reappearance function?


% load tracking data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\tracked.mat', 'locations')

% user settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\botTest.mp4';

unaryWeight = 1;
pairwiseWeight = .1;
occludedWeight = .1;
occlusionGridSpacing = 20;

maxVelocity = 30;
paws = 1:4;
vidDelay = .05;

% initializations
frameHeight = 242; % !!! hacky temp
frameWidth = 398;
[gridX, gridY] = meshgrid(1:occlusionGridSpacing:frameWidth,...
                          1:occlusionGridSpacing:frameHeight);
gridPts = [gridX(:), gridY(:)];
numOccluded = size(gridPts,1);


% compute unary potentials
unary = cell(1,length(locations));

for i = 1:length(locations)

    % ensure first unaries are not zero (this will cause all paths to be zero for all time!)
    if (i==1 && unaryWeight==0)
        unaryWeightTemp = 1;
    else
        unaryWeightTemp = unaryWeight;
    end

    frameUnaries = nan(4, length(locations(i).x) + numOccluded);
    
    frameUnaries(1,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 0, numOccluded, unaryWeightTemp); % RH
    frameUnaries(2,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 1, numOccluded, unaryWeightTemp); % RF
    frameUnaries(3,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 0, numOccluded, unaryWeightTemp); % LH
    frameUnaries(4,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 1, numOccluded, unaryWeightTemp); % LF
    
    
    unary{i} = frameUnaries;
    
end



% compute pairwise potentials
pairwise = cell(1,length(locations)-1);

for i = 1:length(locations)-1
    
    pairwise{i} = getPairwisePotentials([locations(i+1).x, locations(i+1).y], [locations(i).x, locations(i).y],...
                                        maxVelocity, pairwiseWeight, occludedWeight, gridPts);
    
end




% find most probable paths!

objectNum = size(unary{1}, 1);
labels = nan(length(unary), objectNum);
nodeScores = cell(1, length(unary));
backPointers = cell(1, length(unary)-1);

% initializations

nodeScores{1} = unary{1};


for i = 2:length(unary)
    
    nodeScores{i}   = nan(objectNum, size(pairwise{i-1},1));
    backPointers{i-1} = nan(objectNum, size(pairwise{i-1},1));
    
    for j = 1:objectNum
        
        previousScores = nodeScores{i-1}(j,:);
        currentUnary   = unary{i}(j,:)';
        
        allTransitionScores = repmat(previousScores, length(currentUnary), 1) .* ...
                              (pairwise{i-1} + repmat(currentUnary, 1, length(previousScores)));
        allTransitionScores(pairwise{i-1}==0) = 0; % make invalid transitions impossible
        
        [nodeScores{i}(j,:), backPointers{i-1}(j,:)] = max(allTransitionScores, [], 2);
        backPointers{i-1}(j,:) = backPointers{i-1}(j,:);% .* double(nodeScores{i}(j,:)>0);
    end
end


% back trace

[pathScores, labels(end,:)] = max(nodeScores{end}, [], 2);
%
for i = fliplr(1:length(unary)-1)
    for j = 1:objectNum
        
        labels(i,j) = backPointers{i}(j, labels(i+1,j));

    end
end




%% VISUALIZE TRACKING


% initializations
startFrame = 1;
vid = VideoReader(vidFile);
sampleFrame = rgb2gray(read(vid,startFrame));
totalFrames = vid.NumberOfFrames;
cmap = winter(length(paws));

% prepare figure
close all; figure('position', [567 383 698 400], 'color', 'black'); colormap gray


rawIm = image(sampleFrame, 'CDataMapping', 'scaled');
rawAxis = gca;
set(rawAxis, 'visible', 'off')
hold on;
scatterPts =    scatter(rawAxis, zeros(1,length(paws)), zeros(1,length(paws)), 200, cmap, 'filled'); hold on
scatterPtsAll = scatter(rawAxis, 0, 0, 200, 'green', 'linewidth', 2);


for i = startFrame:totalFrames
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % update figure
    set(rawIm, 'CData', frame);
    
    xs = [locations(i).x; gridPts(:,1)];
    ys = [locations(i).y; gridPts(:,2)];
    
%     [~, maxInds] = max(nodeScores{i}(paws,:), [], 2);
%     set(scatterPts, 'XData', xs(maxInds), 'YData', ys(maxInds));

    set(scatterPts, 'XData', xs(labels(i,paws)), 'YData', ys(labels(i,paws)));

    set(scatterPtsAll, 'XData', locations(i).x, 'YData', locations(i).y);
    
    % pause to reflcet on the little things...
    pause(vidDelay);
end










