function imgs = showTrackingOverFrames(session, trials, imgNum, varargin)

% generates figures with raw tracking data for a single trial show, along
% with images overlaid on top of the traces that show snapshots of
% the behavior over time! so cool // imgs is number of imgs to show on the
% left and right of center img, which contains obstacle


% settings
s.showFig = false;  % whether to show figure (otherwise is only written to disk)
s.writeToDisk = false;
s.topOnly = false;  % whether to show only top view
s.kinematicsOnly = false;  % if true, don't show frames or overwrite obstalce
s.scatLines = false;  % if true, use scatter instead of lines to show kinematics
s.imgSpacing = 0; % spacing between image tiles
s.obsFramePixels = 300; % in the middle frame, at what pixel is the obstacle
s.featuresToShow = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', ...
                    'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot'};

s.edgeFading = 60; % fading at the edges of images
s.contrastLims = [.2 .8]; % pixels at these proportional values are mapped to 0 and 255
s.pawColors = jet(4);
s.scatSize = 80;
s.obsWidth = 5;
s.alpha = .8;


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions', 'wheelPositions', 'pixelsPerM', ...
    'wheelTimes', 'wheelToObsPixPosMappings', 'obsPixPositionsUninterped');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);
obsPixPositionsContinuous = getObsPixPositionsContinuous(...
    wheelToObsPixPosMappings, wheelTimes, wheelPositions, frameTimeStamps, ...
    obsPixPositions, obsPixPositionsUninterped, obsOnTimes, obsOffTimes);
fade = repmat(linspace(0,1,s.edgeFading), vidTop.Height+vidBot.Height, 1);

imgs = cell(1, length(trials));
for t = 1:length(trials)
    
    trial = trials(t);
    trialLocations = locations;
    trialLocations(:,1,:) = trialLocations(:,1,:) - obsPixPositionsContinuous(trial,:)'; % unravel coordinates to account for movement of wheel
    
    
    middleFrameInd = find(obsPixPositionsContinuous(trial,:)<s.obsFramePixels, 1, 'first'); % ind of frame where mouse is getting over obstacle
    middleFrameObsPos = round(obsPixPositionsContinuous(trial,middleFrameInd));
    frameObsPosits = middleFrameObsPos + fliplr(-imgNum:imgNum)*(vidTop.Width+s.imgSpacing); % obsPixPositions for frames, from left to right
    imgInds = knnsearch(obsPixPositionsContinuous(trial,:)', frameObsPosits'); % inds of imgs corresponding to the desired frame positions (framePosits)
    imgOffsets = -frameObsPosits + max(frameObsPosits) + 1; % positions in the image of all of the frames
    obsPos = middleFrameObsPos + imgOffsets(imgNum+1); % position of obstacle in final frame
    imgDims = [vidBot.Height+vidTop.Height imgOffsets(end)+vidTop.Width];
    allImgs = double(zeros(imgNum, imgDims(1), imgDims(2)));
    if s.showFig; vis = 'on'; else; vis = 'off'; end
    figure('name', sprintf('%s, trial %i', session, trial), ...
        'Color', 'white', 'Position', [100 442 imgDims(2) imgDims(1)], 'MenuBar', 'none', 'visible', vis)

    if ~s.kinematicsOnly
        
        % create image montage
        for i = 1:length(imgInds)
            img = double(cat(1, rgb2gray(read(vidTop, imgInds(i))), rgb2gray(read(vidBot, imgInds(i)))));
            img(:,1:s.edgeFading) = img(:,1:s.edgeFading) .* fade;  % left fade
            img(:,end+1-s.edgeFading:end) = img(:,end+1-s.edgeFading:end) .* fliplr(fade);  % right fade
            allImgs(i,:,imgOffsets(i)+1:imgOffsets(i)+vidTop.Width) = double(img)*imgNum;
        end


        % pimp and show image
        colormap gray
        montage = squeeze(mean(allImgs,1));
        montage = uint8(montage * (255/max(montage(:))));
        montage = imadjust(montage, s.contrastLims, [0 1]);
        image(montage, 'cdatamapping', 'scaled'); hold on;
    end
    
    
    % and kinematic traces or scatters
    for i = 1:length(features)
        if ismember(features{i}, s.featuresToShow)
            if contains(features{i}, 'paw')
                pawNum = str2double(features{i}(4));
                color = s.pawColors(pawNum,:);
            else
                color = [.8 .8 .8];
            end
            
            % traces
            x = trialLocations(imgInds(1):imgInds(end),1,i) + obsPos;
            y = trialLocations(imgInds(1):imgInds(end),2,i);
            if s.scatLines
                scatter(x, y, 20, color, 'filled', 'MarkerFaceAlpha', s.alpha); hold on
            else
                plot(x, y, 'LineWidth', 2, 'Color', [color s.alpha]); hold on
            end
        end
    end
    
    % add scatters on paws at times of frames
    for i = 1:length(features)
        if ismember(features{i}, s.featuresToShow)
            if contains(features{i}, 'paw')
                pawNum = str2double(features{i}(4));
                color = s.pawColors(pawNum,:);
            else
                color = [.8 .8 .8];
            end
            
            % scatters
            x = trialLocations(imgInds,1,i) + obsPos;
            y = trialLocations(imgInds,2,i);
            scatter(x, y, s.scatSize, color, 'filled')
        end
    end
    
    % add obstacle
    if ~s.kinematicsOnly
        rectangle('Position', [obsPos-s.obsWidth, vidTop.Height, s.obsWidth*2, vidBot.Height], ...
            'EdgeColor', 'none', 'FaceColor', [1 1 1 .8]);
    end


    % pimp fig
    set(gca, 'XLim', [1,imgDims(2)], 'YLim', [1,imgDims(1)], ...
        'visible', 'off', 'Units', 'pixels', 'Position', [1 1 imgDims(2) imgDims(1)]);
    if s.topOnly
        pos = get(gcf, 'position');
        set(gcf, 'position', [pos(1) pos(2) pos(3) vidTop.Height])
        set(gca, 'YLim', [1 vidTop.Height], 'Units', 'normalized', 'Position', [0 0 1 1])
    end

    % save
    if s.writeToDisk
        file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'trackingOverFrames', ...
            [session '_imgs' num2str(imgNum*2+1) '_trial' num2str(trial) '.png']);
        fprintf('writing %s to disk...\n', file)
        saveas(gcf, file)
    end
    
    imgs{t} = frame2im(getframe(gcf));
end
