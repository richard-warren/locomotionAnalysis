
view = 'Bot';

% settings
vidFile = ['C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\run' view '.mp4'];
classifier = ['C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\paw' view '.mat'];
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\';
showTracking = true;

startFrame = 1;
overlapThresh = .5;
scoreThresh = 0;


% load classifier
load(classifier, 'model', 'subHgt', 'subWid')


% initializations
vid = VideoReader(vidFile);
sampleFrame = rgb2gray(read(vid,startFrame));
xMin = 20; % x and yMin are a temporary hack until i crop the videso properly
yMin = 15;
totalFrames = vid.NumberOfFrames;
cmap = winter(4);
kernel = reshape(model.w, subHgt, subWid);


% prepare figure
if showTracking

    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(2,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
    hold on; scatterPts = scatter(rawAxis, 0, 0, 200, 'filled');

    predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
end


locations = struct();

for i = startFrame:totalFrames
    
%     disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % filter with svm
    frameFiltered = - (conv2(double(frame), kernel, 'same') - model.rho);
    
    frameFiltered(frameFiltered < scoreThresh) = 0;
    frameFiltered(1:yMin,:) = 0;
    frameFiltered(:,1:xMin) = 0;
    [x, y, scores] = nonMaximumSupress(frameFiltered, [subHgt subWid], overlapThresh);  
    
    
    % store data
    locations(i).x = x;
    locations(i).y = y;
    locations(i).scores = scores;
    
    if showTracking
        
        % update figure
        set(rawIm, 'CData', frame);
        set(predictIm, 'CData', frameFiltered)
        set(scatterPts, 'XData', x, 'YData', y);
        
        % pause to reflcet on the little things...
        pause(.001);
    end
end

save([dataDir 'tracked' view '.mat'], 'locations');
close all




