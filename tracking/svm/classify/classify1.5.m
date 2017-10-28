
% USER SETTINGS

% settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\botTest.mp4';
classifier = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawBot.mat';
startFrame = 1;
boxSize = [37 37];
overlapThresh = .2;
scoreThresh = .75;

% load sample frame and sample sub-frame
vid = VideoReader(vidFile);
sampleFrame = rgb2gray(read(vid,startFrame));
frmHgt = size(sampleFrame,1);
frmWid = size(sampleFrame,2);
xMin = 20; % x and yMin are a temporary hack until i crop the videso properly
yMin = 15;

% load classifier
load(classifier)

% prepare figure
close all; figure('position', [680 144 698 834]); colormap gray

rawAxis = subaxis(2,1,1, 'spacing', 0.01, 'margin', .01);
rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
hold on; scatterPts = scatter(rawAxis, 0, 0, 200, 'filled');

predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');

%
for i = startFrame:vid.numberofframes
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % filter with svm
    frameFiltered = - (conv2(frame, reshape(model.w, boxSize(1), boxSize(2))', 'same') - model.rho);
    frameFiltered(frameFiltered < scoreThresh) = 0;
    frameFiltered(1:yMin,:) = 0;
    frameFiltered(:,1:xMin) = 0;
    [trackedPts, ~] = nonMaximumSupress(frameFiltered, boxSize, overlapThresh);
    
    
    
    % update figure
    set(rawIm, 'CData', frame);
    set(predictIm, 'CData', frameFiltered)
    set(scatterPts, 'XData', trackedPts(2,:), 'YData', trackedPts(1,:));

    % pause to reflcet on the little things...
    pause(.001);
end





