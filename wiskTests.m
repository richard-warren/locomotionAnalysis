%% use crude thresholding algorithm to git whisk polygon

session = 'wiskTest2';

vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runWisk.mp4']);
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions');
bgRaw = getBgImage(vid, 500, false); % !!! need to ensure obstacle is not in these frames
%%

% settings
frameInds = find(obsPixPositions>300 & obsPixPositions<500);
thresh = 10;
xMinMax = [80 305];
yMinMax = [70 250];
xMin = 80;
edgeMin = 80;
bgThresh = 40;
dilation = 5;
obsThresh = 200;

% initializations
dilateKernel = strel('diamond', dilation);
bg = imdilate(bgRaw, dilateKernel);
bgThreshed = bg>bgThresh;
bgThreshed = imdilate(bgThreshed, dilateKernel);
bgThreshed = ~bgThreshed;
close all; figure('position', [2000 2 561 994]);
subplot(2,1,1)
imRaw = imagesc(rgb2gray(read(vid,1))); hold on;

contourPlot = plot(0, 0, 'color', 'red', 'linewidth', 3);
subplot(2,1,2)
imThresh = imagesc(rgb2gray(read(vid,1))>thresh); hold on;

for i = frameInds(5000:end)
    
    % get frame and thresholded frame
    frameRaw = rgb2gray(read(vid,i));
%     frame = frameRaw - bg;
    obsMask = uint8(frameRaw < obsThresh);
    obsMask = imdilate(obsMask, dilateKernel);
    frameThreshed = (frameRaw >= thresh);
    frameThreshed = frameThreshed .* bgThreshed;
%     frameThreshed([1:yMinMax(1), yMinMax(2):end],:) = 0;
%     frameThreshed(:, [1:xMinMax(1), xMinMax(2):end]) = 0;
    
%     
%     frameThreshed(:, 1:xMin) = 0;
    
    % get boundary points
    [yTip, xTip] = ind2sub(size(frameThreshed), find(frameThreshed));
    
    try
        pts = convhull(xTip, yTip);
        xTip = xTip(pts);
        yTip = yTip(pts);

        validInds = xTip>xMinMax(1) & yTip>yMinMax(1);
        xTip = xTip(validInds);
        yTip = yTip(validInds);
    catch
        xTip = [];
    end

    
    
    set(imRaw, 'CData', frameRaw);
%     set(scat, 'XData', x, 'YData', y)
    set(contourPlot, 'XData', xTip, 'YData', yTip)
    set(imThresh, 'CData', frameThreshed);
    pause(.05);
    
    
end

%% prep file for whisker tracking

dir = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest2\';

vid = VideoReader([dir 'runWisk.mp4']);
vidWrite = VideoWriter([dir 'runWiskEdited.mp4']);

for i = 1:vid.NumberOfFrames
    
    
    
end



%% show whiski whisker tracking

% settings
vidFile = [getenv('OBSDATADIR') 'sessions\' session '\runWisk.mp4'];
wiskData = [getenv('OBSDATADIR') 'sessions\' session '\runWisk.measures'];
minLength = 0;

% initializations
vid = VideoReader(vidFile);
measurements = LoadMeasurements(wiskData);


close all; figure('position', [2000 600 550 400]);
im = imshow(rgb2gray(read(vid,1))); hold on;
scatBase = scatter(10, 10, 50, 'blue', 'filled');
scatTip = scatter(10, 10, 50, 'red', 'filled');
polyPlot = plot(0, 0, 'color', 'red', 'linewidth', 3);


for i = 10000:vid.NumberOfFrames
    
    % get frame
    frame = read(vid,i);
    
    % get wisk base and tips
    inds = ([measurements.fid]==i) & [measurements.length]>minLength;
    xBase = [measurements(inds).follicle_x];
    yBase = [measurements(inds).follicle_y];
    xTip = [measurements(inds).tip_x];
    yTip = [measurements(inds).tip_y];
    
    % get convex hull around tips
    try
        pts = convhull(xTip, yTip);
    catch
        pts = [];
    end
    
    % update preview
    set(im, 'CData', frame);
    set(scatTip, 'XData', xTip, 'YData', yTip);
    set(scatBase, 'XData', xBase, 'YData', yBase);
    set(polyPlot, 'XData', xTip(pts), 'YData', yTip(pts));
    pause(.1);
    
end




















