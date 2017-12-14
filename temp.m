

vid = VideoReader('C:\Users\rick\Desktop\wiskTest\wisk.mp4');
bgRaw = getBgImage(vid, 1000, false); % !!! need to ensure obstacle is not in these frames
%%

% settings
thresh = 10;
xMinMax = [80 305];
yMinMax = [70 250];
xMin = 80;
edgeMin = 80;
bgThresh = 40;
dilation = 10;

% initializations
bgThreshed = bg>bgThresh;
dilateKernel = strel('diamond', dilation);
bg = imdilate(bgRaw, dilateKernel);
bgThreshed = imdilate(bgThreshed, dilateKernel);
bgThreshed = ~bgThreshed;
close all; figure('position', [2000 2 561 994]);
subplot(2,1,1)
imRaw = imagesc(rgb2gray(read(vid,1))); hold on;

contourPlot = plot(0, 0, 'color', 'red', 'linewidth', 3);
subplot(2,1,2)
imThresh = imagesc(rgb2gray(read(vid,1))>thresh); hold on;

for i = 1:vid.NumberofFrames
    
    % get frame and thresholded frame
    frameRaw = rgb2gray(read(vid,i));
    frame = frameRaw - bg;
    frameThreshed = (frame >= thresh);
    frameThreshed = frameThreshed .* bgThreshed;
    frameThreshed([1:yMinMax(1), yMinMax(2):end],:) = 0;
    frameThreshed(:, [1:xMinMax(1), xMinMax(2):end]) = 0;
    
%     
%     frameThreshed(:, 1:xMin) = 0;
    
    % get boundary points
    [y, x] = ind2sub(size(frameThreshed), find(frameThreshed));
    
    try
        pts = convhull(x, y);
        x = x(pts);
        y = y(pts);

        validInds = x>xMinMax(1) & y>yMinMax(1);
        x = x(validInds);
        y = y(validInds);
    catch
        x = [];
    end

    
    
    set(imRaw, 'CData', frameRaw);
%     set(scat, 'XData', x, 'YData', y)
    set(contourPlot, 'XData', x, 'YData', y)
    set(imThresh, 'CData', frameThreshed);
    pause(.05);
    
    
end

%%



close all; imshow(frame); hold on; scatter(x, y)