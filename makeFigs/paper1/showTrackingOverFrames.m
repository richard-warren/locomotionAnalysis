function showTrackingOverFrames(session, trials, imgs)

% generates figures with raw tracking data for a single trial show, along
% with images overlaid on top of the traces that show snapshots of
% the behavior over time! so cool // imgs is number of imgs to show on the
% left and right of center img, which contains obstacle


% settings
% session = '180703_000';
% trials = 10;
% imgs = 1; % how many images to add to the left and right
imgSpacing = 0; % spacing between image tiles
obsFramePixels = 300; % in the middle frame, at what pixel is the obstacle
featuresToShow = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', ...
                  'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot'};

edgeFading = 60; % fading at the edges of images
contrastLims = [.2 .8]; % pixels at these proportional values are mapped to 0 and 255
pawColors = jet(4);
scatSize = 80;
obsWidth = 5;
lineAlpha = .8;


% initializations
vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions', 'wheelPositions', ...
    'wheelTimes', 'wheelCenter', 'wheelRadius', 'mToPixMapping', 'obsPosToWheelPosMappings', 'obsPixPositionsUninterped');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
% topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
% stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, ...
%     wheelTimes, wheelCenter, wheelRadius, 250, mToPixMapping(1));
obsPixPositionsContinuous = getObsPixPositionsContinuous(...
    obsPosToWheelPosMappings, wheelTimes, wheelPositions, frameTimeStamps, ...
    obsPixPositions, obsPixPositionsUninterped, obsOnTimes, obsOffTimes);
fade = repmat(linspace(0,1,edgeFading), vidTop.Height+vidBot.Height, 1);






for trial = trials
    
    trialLocations = locations;
    trialLocations(:,1,:) = trialLocations(:,1,:) - obsPixPositionsContinuous(trial,:)'; % unravel coordinates to account for movement of wheel
    
    middleFrameInd = find(obsPixPositionsContinuous(trial,:)<obsFramePixels, 1, 'first'); % ind of frame where mouse is getting over obstacle
    middleFrameObsPos = round(obsPixPositionsContinuous(trial,middleFrameInd));
    frameObsPosits = middleFrameObsPos + fliplr(-imgs:imgs)*(vidTop.Width+imgSpacing); % obsPixPositions for frames, from left to right
    imgInds = knnsearch(obsPixPositionsContinuous(trial,:)', frameObsPosits'); % inds of imgs corresponding to the desired frame positions (framePosits)
    imgOffsets = -frameObsPosits + max(frameObsPosits) + 1; % positions in the image of all of the frames
    obsPos = middleFrameObsPos + imgOffsets(imgs+1); % position of obstacle in final frame
    
    imgDims = [vidBot.Height+vidTop.Height imgOffsets(end)+vidTop.Width];
    allImgs = double(zeros(imgs, imgDims(1), imgDims(2)));
    
    figure('Color', 'white', 'Position', [100 442 imgDims(2) imgDims(1)], 'MenuBar', 'none', 'visible', 'off')


    % create image montage
    for i = 1:length(imgInds)
        img = double(cat(1, rgb2gray(read(vidTop, imgInds(i))), rgb2gray(read(vidBot, imgInds(i)))));
        img(:,1:edgeFading) = img(:,1:edgeFading) .* fade; % left fade
        img(:,end+1-edgeFading:end) = img(:,end+1-edgeFading:end) .* fliplr(fade); % right fade
        allImgs(i,:,imgOffsets(i)+1:imgOffsets(i)+vidTop.Width) = double(img)*imgs;
    end
    
    
    % pimp and show image
    colormap gray
    montage = squeeze(mean(allImgs,1));
    montage = uint8(montage * (255/max(montage(:))));
    montage = imadjust(montage, contrastLims, [0 1]);
    image(montage, 'cdatamapping', 'scaled'); hold on;
    
    
    % and kinematic traces and scatters
    for i = 1:length(features)
        if ismember(features{i}, featuresToShow)
            if contains(features{i}, 'paw')
                pawNum = str2double(features{i}(4));
                color = pawColors(pawNum,:);
            else
                color = [.8 .8 .8];
            end
            
            % traces
            x = trialLocations(imgInds(1):imgInds(end),1,i) + obsPos;
            y = trialLocations(imgInds(1):imgInds(end),2,i);
            plot(x, y, 'LineWidth', 2, 'Color', [color lineAlpha]); hold on
            
            % scatters
            x = trialLocations(imgInds,1,i) + obsPos;
            y = trialLocations(imgInds,2,i);
            scatter(x, y, scatSize, color, 'filled')
        end
    end
    
    
    % add obstacle
    rectangle('Position', [obsPos-obsWidth, vidTop.Height, obsWidth*2, vidBot.Height], ...
        'EdgeColor', 'none', 'FaceColor', [1 1 1 .8]);


    % pimp fig
    set(gca, 'XLim', [1,imgDims(2)], 'YLim', [1,imgDims(1)], ...
        'visible', 'off', 'Units', 'pixels', 'Position', [1 1 imgDims(2) imgDims(1)]);

    % save
    file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'trackingOverFrames', ...
        [session '_imgs' num2str(imgs*2+1) '_trial' num2str(trial) '.png']);
    fprintf('writing %s to disk...\n', file)
    saveas(gcf, file)
end
