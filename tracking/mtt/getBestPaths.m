
% load tracking data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\tracked.mat', 'locations')

% user settings
occludedCost = .1;
occludedStates = 1;
maxVelocity = 15;
velocityWeight = 0;
anchorWeight = 1;
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\botTest.mp4';
paws = 1:4;
delay = .05;

% initializations
frameHeight = 242; % !!! hacky temp
frameWidth = 398;


% compute unary potentials
unary = cell(1,length(locations));

for i = 1:length(locations)
    
    frameUnaries = nan(4, length(locations(i).x) + occludedStates);
    
    frameUnaries(1,1:end-occludedStates) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 0); % LH
    frameUnaries(2,1:end-occludedStates) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 1); % LF
    frameUnaries(3,1:end-occludedStates) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 0); % RH
    frameUnaries(4,1:end-occludedStates) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 1); % RF
    
    frameUnaries(:,end-occludedStates+1:end) = 0;%occludedCost;
    
    unary{i} = frameUnaries * anchorWeight; % documentation on match2nd appears to be incorrect... this matrix should be flipped i think...
    
end



% compute pairwise potentials
pairwise = cell(1,length(locations)-1);

for i = 1:length(locations)-1
    
    pairwise{i} = getPairwisePotentials([locations(i+1).x, locations(i+1).y], [locations(i).x, locations(i).y], maxVelocity, velocityWeight, occludedCost);
    
end



% function getBestPaths(unary, pairwise)

% CURRENTLY HAS ONLY ONE OCCLUDED STATE, AND NO EXCLUSION CONSTRAINTS!!!


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
cmap = winter(length(paws));

% prepare figure
close all; figure('position', [680 144 698 400]); colormap gray


rawIm = image(sampleFrame, 'CDataMapping', 'scaled');
rawAxis = gca;
hold on;
scatterPts =    scatter(rawAxis, zeros(1,length(paws)), zeros(1,length(paws)), 200, cmap, 'filled'); hold on
scatterPtsAll = scatter(rawAxis, 0, 0, 200, [1 1 1]);


for i = startFrame:totalFrames
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % update figure
    set(rawIm, 'CData', frame);
    ind = backPointers{i+1}(paws,1);
    
    xs = [locations(i).x; 0]; % add zero for occluded state
    ys = [locations(i).y; 0]; % add zero for occluded state
    
%     set(scatterPts, 'XData', xs(backPointers{i+1}(paws,1)), 'YData', ys(backPointers{i+1}(paws,1)));
    set(scatterPtsAll, 'XData', locations(i).x, 'YData', locations(i).y);
    
    % pause to reflcet on the little things...
    pause(delay);
end










