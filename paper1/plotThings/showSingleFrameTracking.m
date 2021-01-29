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
s.otherColors = [.8 .8 .8];
s.ind = nan;
s.removeRows = [];  % start and end rows to remove from the image // use to cut out the bright line separating top and bottom views
s.featuresToShow = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', ...
                    'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot', ...
                    'nose_top', 'nose_bot', 'tailMid_top', 'tailMid_bot', 'tailBase_top', 'tailBase_bot'};


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions', 'pixelsPerM', 'frameTimeStampsWisk');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM, 'scoreThresh', scoreThresh);


% find frame number
if isnan(s.ind)
    s.ind = find(frameTimeStamps > obsOnTimes(trial) & ...
                 frameTimeStamps < obsOffTimes(trial) & ...
                 obsPixPositions' < s.obsFramePixels, 1, 'first');  % ind of frame where mouse is getting over obstacle
end


% load frame
if s.addWiskCam
    vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
    img = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, s.ind, 'runContrast', s.contrastLims);
else
    img = imadjust(read(vid, s.ind), s.contrastLims);
end

if ~isempty(s.removeRows)
    img(s.removeRows(1):s.removeRows(2),:) = [];
    locations(:,2,:) = locations(:,2,:) - ...
        (locations(:,2,:)>s.removeRows(2))*diff(s.removeRows);  % shift tracking to compenstate for taking out rows in the middle
end



figure('name', sprintf('%s, trial %i', session, trial), ...
    'Color', 'white', 'Position', [100 100 size(img,2) size(img,1)], 'MenuBar', 'none')
colormap gray
image(img, 'cdatamapping', 'scaled'); hold on;
set(gca, 'XLim', [1, size(img,2)], 'YLim', [1, size(img,1)], ...
    'visible', 'off', 'Units', 'normalized', 'Position', [0 0 1 1]);


% add kinematic traces or scatters
sizes = [linspace(s.trailingSizes(1), s.trailingSizes(2), (s.numPoints-1)) s.mainSize];
alphas = linspace(0, 1, s.numPoints);
for i = 1:length(features)
    if ismember(features{i}, s.featuresToShow)
        if contains(features{i}, 'paw')
            pawNum = str2double(features{i}(4));
            color = s.pawColors(pawNum,:);
        else
            color = s.otherColors;
        end

        % traces
        inds = s.ind-s.deltaFrames*(s.numPoints-1) : s.deltaFrames : s.ind;
        x = locations(inds, 1, i);
        y = locations(inds, 2, i);
        for j = 1:s.numPoints
            scatter(x(j), y(j), sizes(j), color, 'filled', 'MarkerFaceAlpha', alphas(j)); hold on
        end
    end
end


% % add obstacle
% obsPos = obsPixPositions(s.ind);
% rectangle('Position', [obsPos-s.obsWidth, vidTop.Height, s.obsWidth*2, vidBot.Height], ...
%     'EdgeColor', 'none', 'FaceColor', [1 1 1 .6]);

