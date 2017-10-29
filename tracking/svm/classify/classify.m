
% USER SETTINGS

% settings
vidFile = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\testVideo\botTest.mp4';
classifier = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\classifiers\pawBot.mat';
dataDir = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\trackedData\';

startFrame = 1;
overlapThresh = .5;
scoreThresh = .5;
simpleThresh = 150;

% initializations
vid = VideoReader(vidFile);
sampleFrame = rgb2gray(read(vid,startFrame));
xMin = 20; % x and yMin are a temporary hack until i crop the videso properly
yMin = 15;
totalFrames = vid.NumberOfFrames;
cmap = winter(4);

% load classifier
load(classifier, 'model', 'subHgt', 'subWid')

% prepare figure
close all; figure('position', [680 144 698 834]); colormap gray

rawAxis = subaxis(2,1,1, 'spacing', 0.01, 'margin', .01);
rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
hold on; scatterPts = scatter(rawAxis, 0, 0, 200, 'filled');

predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');


locations = struct();

for i = startFrame:totalFrames
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % filter with svm
    frameFiltered = - (conv2(double(frame), reshape(model.w, subHgt, subWid), 'same') - model.rho);
    
    frameFiltered(frameFiltered < scoreThresh) = 0;
    frameFiltered(1:yMin,:) = 0;
    frameFiltered(:,1:xMin) = 0;
    [x, y] = nonMaximumSupress(frameFiltered, [subHgt subWid], overlapThresh);  
    
    % update figure
    set(rawIm, 'CData', frame);
    set(predictIm, 'CData', frameFiltered)
    
    
    % compute proximity to corner
    
    set(scatterPts, 'XData', x(labels(i,:)), 'YData', y(labels(i,:)));

    % store data
    locations(i).x = x;
    locations(i).y = y;
    
    % pause to reflcet on the little things...
    pause(.001);
end

% save([dataDir 'tracked.mat'], 'locations');





