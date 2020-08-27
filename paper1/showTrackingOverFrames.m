function imgs = showTrackingOverFrames(session, trials, imgNum, varargin)

% generates figures with raw tracking data for a single trial show, along
% with images overlaid on top of the traces that show snapshots of
% the behavior over time! so cool // imgs is number of imgs to show on the
% left and right of center img, which contains obstacle


% settings
s.showFig = false;  % whether to show figure (otherwise is only written to disk)
s.topOnly = false;  % whether to show only top view
s.kinematicsOnly = false;  % if true, don't show frames or overwrite obstalce
s.scatLines = false;  % if true, use scatter instead of lines to show kinematics
s.imgSpacing = 0; % spacing between image tiles
s.overObsPixels = 10; % in the middle frame, how many pixels beyond the obstacle should the leading paw be

s.edgeFading = 60; % fading at the edges of images
s.contrastLims = [.2 .8]; % pixels at these proportional values are mapped to 0 and 255
s.pawColors = jet(4);
s.scatSize = 80;
s.obsWidth = 5;
s.alpha = .8;
s.highlightStepOver = true;  % whether to trace over step over obstacles with thicker line


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), 'kinData');
fade = repmat(linspace(0,1,s.edgeFading), vidTop.Height+vidBot.Height, 1);

imgs = cell(1, length(trials));
for t = 1:length(trials)
    
    trial = trials(t);
    trialLocations = kinData(trial).locationsPix;  % convert back to pixels
    trialLocations(:,1,:) = trialLocations(:,1,:) - kinData(trial).obsPixPos';
    
    middleFrameInd = find(any(squeeze(trialLocations(:,1,:)) > s.overObsPixels, 2), 1, 'first');  % ind of frame where mouse is getting over obstacle
    middleFrameObsPos = round(kinData(trial).obsPixPos(middleFrameInd));
    frameObsPosits = middleFrameObsPos + fliplr(-imgNum:imgNum)*(vidTop.Width+s.imgSpacing);  % obsPixPositions for frames, from left to right
    inds = knnsearch(kinData(trial).obsPixPos', frameObsPosits');  % inds of imgs corresponding to the desired frame positions (framePosits)
    imgInds = kinData(trial).trialInds(inds);
    imgOffsets = -frameObsPosits + max(frameObsPosits) + 1; % positions in the image of all of the frames
    obsPos = middleFrameObsPos + imgOffsets(imgNum+1); % position of obstacle in final frame
    imgDims = [vidBot.Height+vidTop.Height imgOffsets(end)+vidTop.Width];
    allImgs = double(zeros(imgNum, imgDims(1), imgDims(2)));
    if s.showFig; vis = 'on'; else; vis = 'off'; end
    figure('name', sprintf('%s, trial %i', session, trial), ...
        'Color', 'white', 'Position', [2000 442 imgDims(2) imgDims(1)], 'MenuBar', 'none', 'visible', vis)

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
    scats = cell(1,4);
    for i = 1:4
        color = s.pawColors(kinData(trial).pawOverSequence==i, :);  % arrange colors to correspond to sequence of paws over obstalce

        % traces
        x = repmat(trialLocations(inds(1):inds(end),1,i) + obsPos, 2, 1);
        y_bot = trialLocations(inds(1):inds(end),2,i);
        y_top = trialLocations(inds(1):inds(end),3,i);
        y = [y_bot; y_top];
        if s.scatLines
            scatter(x, y, 20, color, 'filled', 'MarkerFaceAlpha', s.alpha); hold on
        else
            plot(x, y, 'LineWidth', 2, 'Color', [color s.alpha]); hold on
        end
        
        
        % highlight step over obstacle
        if s.highlightStepOver
            
            numModSteps = size(kinData(trial).modifiedLocations{i}, 1);
            bins =kinData(trial).modifiedStepIdentities(:,i)==numModSteps;
            
            % plot bot
            x = trialLocations(bins,1,i) + obsPos;
            y_bot = trialLocations(bins,2,i);
            y_top = trialLocations(bins,3,i);
            plot(x, y_bot, 'LineWidth', 4, 'Color', [color .8])
            plot(x, y_top, 'LineWidth', 4, 'Color', [color .8])
        end
        
        
        % scatter over paws
        x = repmat(trialLocations(inds,1,i) + obsPos, 2, 1);
        y_bot = trialLocations(inds,2,i);
        y_top = trialLocations(inds,3,i);
        y = [y_bot; y_top];
        scats{i} = scatter(x, y, s.scatSize, color, 'filled', 'MarkerEdgeColor', 'white', 'LineWidth', 2);
    end
    uistack(cat(2, scats{:}), 'top')  % move scatter to top
    
    
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

    imgs{t} = frame2im(getframe(gcf));
end
