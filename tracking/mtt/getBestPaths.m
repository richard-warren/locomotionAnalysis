% function getBestPaths(unary, pairwise)

% CURRENTLY HAS ONLY ONE OCCLUDED STATE, AND NO EXCLUSION CONSTRAINTS!!!


objectNum = size(unary{1}, 1);
labels = nan(length(unary), objectNum);
nodeScores = cell(1, length(unary));
backPointers = cell(1, length(unary)-1);

% initializations
% (pick most probable label based solely on unary potentials)
nodeScores{1} = unary{1};


for i = 2:length(unary)
    
    nodeScores{i}   = nan(objectNum, size(unary{i},2));
    backPointers{i-1} = nan(objectNum, size(unary{i},2));
    
    for j = 1:objectNum
        
        previousScores = nodeScores{i-1}(j,:);
        currentUnary   = unary{i}(j,:)';
        
        allTransitionScores = repmat(previousScores, length(currentUnary), 1) + ...
                              pairwise{i-1}*0 + ...
                              repmat(currentUnary, 1, length(previousScores));
        
        [nodeScores{i}(j,:), backPointers{i-1}(j,:)] = max(allTransitionScores, [], 2);
        
    end
end


% back trace

[pathScores, labels(end,:)] = max(nodeScores{end}, [], 2);

for i = fliplr(1:length(unary)-1)
    
    labels(i,:) = backPointers{i}(labels(i+1,:));
    
    
end

%% visualize tracking



% USER SETTINGS

% settings
vidFile = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\testVideo\botTest.mp4';
load('C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\trackedData\tracked.mat', 'locations')
paw = 1;

% initializations
startFrame = 1;
vid = VideoReader(vidFile);
sampleFrame = rgb2gray(read(vid,startFrame));
totalFrames = vid.NumberOfFrames;
cmap = winter(4);

% prepare figure
close all; figure('position', [680 144 698 400]); colormap gray


rawIm = image(sampleFrame, 'CDataMapping', 'scaled');
rawAxis = gca;
hold on; scatterPts = scatter(rawAxis, 0, 0, 200, 'filled');


for i = startFrame:totalFrames
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % update figure
    set(rawIm, 'CData', frame);
    set(scatterPts, 'XData', locations(i).x(labels(i,paw)), 'YData', locations(i).y(labels(i,paw)));
    
    % pause to reflcet on the little things...
    pause(.001);
end










