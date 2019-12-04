function showSingleFrameTracking(session, trial, varargin)


% settings
s.addWiskCam = false;  % whether to overlay view from whisker camera
s.trailingSizes = [10 50];  % range of sizes for trailing scatter points
s.mainSize = 100;  % size of main scatter points
s.obsFramePixels = 300;  % where is the obstacle within the frame
s.contrastLims = [0 .8];
s.obsWidth = 5;  % (pixels)
s.wiskBorder = 4;  % white border to draw around whisker camera
s.numPoints = 10;  % number of scatter points to show
s.deltaFrames = 2;  % how many frames apart should each scatter point be
s.pawColors = jet(4);
s.featuresToShow = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', ...
                    'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot', ...
                    'nose_top', 'nose_bot', 'tailMid_top', 'tailMid_bot', 'tailBase_top', 'tailBase_bot'};


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions', 'pixelsPerM', 'frameTimeStampsWisk');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);


% show image
ind = find(frameTimeStamps > obsOnTimes(trial) & ...
           frameTimeStamps < obsOffTimes(trial) & ...
           obsPixPositions' < s.obsFramePixels, 1, 'first');  % ind of frame where mouse is getting over obstacle
img = double(cat(1, rgb2gray(read(vidTop, ind)), rgb2gray(read(vidBot, ind))));
img = uint8(img * (255/max(img(:))));
img = imadjust(img, s.contrastLims, [0 1]);

figure('name', sprintf('%s, trial %i', session, trial), ...
    'Color', 'white', 'Position', [100 100 size(img,2) size(img,1)], 'MenuBar', 'none')
colormap gray
image(img, 'cdatamapping', 'scaled'); hold on;
set(gca, 'XLim', [1, size(img,2)], 'YLim', [1, size(img,2)], ...
    'visible', 'off', 'Units', 'pixels', 'Position', [1 1 size(img,2) size(img,1)]);


% add whisker frame
if s.addWiskCam
    
    % get whisker frame and frame alignment
    vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
    indWisk = knnsearch(frameTimeStampsWisk, frameTimeStamps(ind));
    frameWisk = rgb2gray(read(vidWisk, indWisk));
    frameTop = rgb2gray(read(vidTop, ind));
    [yWiskPos, xWiskPos, wiskScaling] = getSubFramePosition(frameTop(:,:), frameWisk(:,:), .35:.005:.45);
    
    % adjust whisker frame
    frameWisk = imresize(frameWisk, wiskScaling);
    frameWisk = imadjust(frameWisk, [.5 1], [0 1]);
    frameWisk = 255 - frameWisk;
    
    % add border to frame
    frameWisk([1:s.wiskBorder, end-s.wiskBorder:end], :) = 255;
    frameWisk(:, [1:s.wiskBorder, end-s.wiskBorder:end]) = 255;
    frameWisk = repmat(frameWisk,1,1,3);
    
    xInds = xWiskPos:xWiskPos+size(frameWisk,2)-1;
    yInds = yWiskPos:yWiskPos+size(frameWisk,1)-1;
    image(xInds, yInds, frameWisk);
    
    xLims = [1, max(size(img,2), xInds(end))];
    yLims = [min(1,yWiskPos), size(img,1)];
    set(gcf, 'position', [100 100 range(xLims) range(yLims)])
    set(gca, 'xlim', xLims, 'ylim', yLims, 'units', 'normalized', 'position', [0 0 1 1])
end


% add kinematic traces or scatters
sizes = [linspace(s.trailingSizes(1), s.trailingSizes(2), (s.numPoints-1)) s.mainSize];
alphas = linspace(0, 1, s.numPoints);
for i = 1:length(features)
    if ismember(features{i}, s.featuresToShow)
        if contains(features{i}, 'paw')
            pawNum = str2double(features{i}(4));
            color = s.pawColors(pawNum,:);
        else
            color = [.8 .8 .8];
        end

        % traces
        inds = ind-s.deltaFrames*(s.numPoints-1) : s.deltaFrames : ind;
        x = locations(inds, 1, i);
        y = locations(inds, 2, i);
        for j = 1:s.numPoints
            scatter(x(j), y(j), sizes(j), color, 'filled', 'MarkerFaceAlpha', alphas(j)); hold on
        end
    end
end


% add obstacle
obsPos = obsPixPositions(ind);
rectangle('Position', [obsPos-s.obsWidth, vidTop.Height, s.obsWidth*2, vidBot.Height], ...
    'EdgeColor', 'none', 'FaceColor', [1 1 1 .6]);