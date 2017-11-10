

% settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runTop.mp4';
classifier = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawTop.mat';
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\';
showTracking = true;
objectNum = 4;
circRoiPts = [55 146; 209 108; 364 135];
xBuffer = 20;


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
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\trackedBot', 'locationsBot');

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


locationsTop = struct();

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
    
    for j = 1:objectNum
        
        % check if paw is occluded
        occludedByBins = abs(locationsBot(i).x(j) - locationsBot(i).x) < xBuffer & ...
                      (locationsBot(i).y > locationsBot(i).y(j));
        
        if any(occludedByBins)
            locationsTop(i).x(j) = nan;
            locationsTop(i).y(j) = nan;
        else
            validXBins = abs(locationsBot(i).x(j) - x) < xBuffer;
            [~, ind] = max(validXBins .* scores);
            locationsTop(i).x(j) = x(ind);
            locationsTop(i).y(j) = y(ind);
        end
        
    end
    
    % store data
%     locationsTop(i).x = x;
%     locationsTop(i).y = y;
%     
    if showTracking
        
        % update figure
        set(rawIm, 'CData', frame);
        set(predictIm, 'CData', frameFiltered)
%         set(scatterPtsAll, 'XData', x, 'YData', y);
        set(scatterPts, 'XData', locationsTop(i).x, 'YData', locationsTop(i).y);
        
        % pause to reflcet on the little things...
        pause(.05);
    end
end

save([dataDir 'trackedTop.mat'], 'locationsTop');
close all




