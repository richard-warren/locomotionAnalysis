function makeTrackingVidDLC(session, frameInds)

% makes video of top and bot views with tracked positions overlayed
% if duration==0 the video is made for the whole file

% user settings
showInterpedValues = true;
showText = false;
show3D = false;

% 3d plot settings
pawTraceLength = 6;
azimuth = -10;
elevation = 10;
xLims = [0 .1];
yLims = [-0.0381 0.0381]; % 3 inches, which is width of wheel
zLims = [-.01 .02];

fps = 25; % playback fps
circSize = 75;
botPawInds = 1:4;
topPawInds = 8:11;
colorMap = 'hsv';
ommitFeatures = {'gen'};

% load tracking and vid data
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'frameTimeStamps', 'nosePos', 'mToPixMapping');
mToPixFactor = abs(median(mToPixMapping(:,1)));
locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw2.csv']); % get raw tracking data
locationsTable = locationsTable(:,2:end); % remove index column
[locations, features, featurePairInds, isInterped] = fixTrackingDLC(locationsTable, frameTimeStamps);
vidTop = VideoReader([getenv('OBSDATADIR') '\sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') '\sessions\' session '\runBot.mp4']);


% initializations
featureBins = ~ismember(features, ommitFeatures);
% set colors s.t. matching features in top and bot view have same color
cmap = eval(sprintf('%s(%i);', colorMap, length(features)));
for i = 1:size(featurePairInds,1)
    cmap(featurePairInds(i,2),:) = cmap(featurePairInds(i,1),:);
end
cmap = cmap(featureBins,:);
vidWrite = VideoWriter([getenv('OBSDATADIR') '\editedVid\' session 'trackingSample.mp4'], 'MPEG-4');
set(vidWrite, 'FrameRate', fps, 'Quality', 50);
open(vidWrite);
height = vidBot.Height+vidTop.Height;
width = vidBot.Width;

% put together xyz for paws only
if show3D
    locations3D = nan(size(locations,1), 3, 4);
    locations3D(:,1:2,:) = locations(:,:,botPawInds);
    locations3D(:,3,:) = locations(:,2,topPawInds);
    locations3D(:,2,:) = locations3D(:,2,:) - (nosePos(2)+vidTop.Height); % subtract midline from all y values
    locations3D(:,3,:) = -locations3D(:,3,:) + nanmean(locations3D(:,3,:)); % flip z and subtract mean
    locations3D = locations3D / abs(mToPixFactor); % convert to meters
    
    width = width*2;
end



% set up figure
mainFig = figure('name', session, 'color', [0 0 0], 'position', [1925, 50, width, height], 'menubar', 'none');
colormap gray
mainAx = subplot(1,1+show3D,1);
imshow = image(zeros(height, vidTop.Width), 'CDataMapping', 'scaled'); hold on;
set(mainAx, 'visible', 'off', 'CLim', [0 255], 'units', 'pixels', ...
    'position', [0 0 vidBot.Width height])
circs = scatter(mainAx, zeros(1,sum(featureBins)), zeros(1,sum(featureBins)), circSize, cmap, 'filled'); hold on
if showText; txt = text(mainAx, 5 , 5, '', 'color', [1 1 1]); end

if show3D
    plotAx = subplot(1,2,2);
    
    for i = 1:size(locations3D,3)
        plotLines{i} = plot3(0,0,0, 'color', cmap(botPawInds(i),:), 'linewidth', 2); hold on;
    end
    set(gca, 'units', 'pixels', 'position', [vidBot.Width+1 0 vidBot.Width height], ...
        'color', 'black', 'view', [azimuth elevation], 'xlim', xLims, 'ylim', yLims, 'zlim', zLims)
    daspect([1 1 1])
end



% write video
for i = frameInds'
    
    % plot frame
    frameTop = read(vidTop,i);
    frameBot = read(vidBot,i);
    frame = cat(1, frameTop, frameBot);
    frame = imadjust(squeeze(frame(:,:,1)), [.1 1], [0 1]);
    set(imshow, 'CData', frame);
    if showText; set(txt, 'String', ['frame ' num2str(i)]); end

    % update circles
    x = locations(i,1,featureBins);
    y = locations(i,2,featureBins);
    
    if ~showInterpedValues
        x(isInterped(i,featureBins)) = nan;
        y(isInterped(i,featureBins)) = nan;
    end
    
    % update 3D plot
    if show3D
        for j = 1:size(locations3D,3)
            set(plotLines{j}, ...
                'XData', locations3D(i-pawTraceLength:i,1,j), ...
                'YData', locations3D(i-pawTraceLength:i,2,j), ...
                'ZData', locations3D(i-pawTraceLength:i,3,j));
        end
    end
    
    set(circs, 'XData', x, 'YData', y);

    writeVideo(vidWrite, getframe(gcf));
end

close(vidWrite);
close(mainFig);


