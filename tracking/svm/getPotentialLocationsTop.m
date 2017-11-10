function getPotentialLocationsTop(showTracking)

% settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runTop.mp4';
classifier = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawTop.mat';
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\';
objectNum = 4;
circRoiPts = [55 146; 209 108; 364 135];


startFrame = 1;
overlapThresh = .7; % .5 for bottom // 
scoreThresh = 0;


% load classifier
load(classifier, 'model', 'subHgt', 'subWid')


% initializations
vid = VideoReader(vidFile);
sampleFrame = rgb2gray(read(vid,startFrame));
totalFrames = vid.NumberOfFrames;
cmap = winter(4);
kernel = reshape(model.w, subHgt, subWid);
wheelMask = getWheelMask(circRoiPts, [vid.Height vid.Width]);
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\locationsBot.mat', 'locationsBot');

% prepare figure
if showTracking

    close all;
    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(2,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
    hold on; scatterPtsAll = scatter(rawAxis, 0, 0, 100, 'filled', 'red');
    hold on; scatterPts = scatter(rawAxis, 0, 0, 200, 'filled', 'blue');

    predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
end


potentialLocationsTop = struct();

for i = startFrame:totalFrames
    
    disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % mask with circle
    frame = frame .* wheelMask;
    
    % filter with svm and apply non-maxima suppression
    frameFiltered = - (conv2(double(frame), kernel, 'same') - model.rho);
    frameFiltered(frameFiltered < scoreThresh) = 0;
    [x, y, scores] = nonMaximumSupress(frameFiltered, [subHgt subWid], overlapThresh);
    
    
    
    % store data
    potentialLocationsTop(i).x = x;
    potentialLocationsTop(i).y = y;
    potentialLocationsTop(i).scores = scores;
    
    
    if showTracking
        
        % update figure
        set(rawIm, 'CData', frame);
        set(predictIm, 'CData', frameFiltered)
        set(scatterPtsAll, 'XData', x, 'YData', y);
        
        % pause to reflcet on the little things...
        pause(.05);
    end
end

save([dataDir 'potentialLocationsTop.mat'], 'potentialLocationsTop');
close all




