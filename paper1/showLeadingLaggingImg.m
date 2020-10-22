function showLeadingLaggingImg(session, trial, varargin)

% make img showing four pics of mouse getting over obstacle, show leading
% fore, lagging fore, leading hind, then lagging hind // overlay kinematic
% traces for one paw per img, as well as randomly selected trials
% transparently overlaid in background // randomly selects valid trial
% for the main image, and randomly selects among valid trials the
% background traces // the file name reflects the randomly selected trial
% for the main image // if imgTrial is supplied by user, it forces this trial to be the
% background image

% settings
s.pawPos = .005; % (m) each frame is selected when the paw of interest is pawPos in front of obs (can be negative too)
s.imgSpacing = 0; % pixels for separating images
s.edgeFading = 60; % fading at the edges of images
s.contrastLims = [.05 .8]; % pixels at these proportional values are mapped to 0 and 255
s.colors = hsv(4);
s.vertical = false;  % whether to lay images on top of one another rather than side by side

s.overlays = 8; % how many trial kinematics to overlay
s.overlayAlpha = .5;
s.overlayWidth = 2;
s.mainWidth = 2;
s.mainAlpha = .5;
s.randSeed = [];  % for selecting the same random trials to show
s.scatter = false;  % whether to scatter rather than plot kinematic overlays



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
pawSequence = [3 2 4 1]; % only include trials with this sequence of paws going over
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), 'kinData')
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'frameTimeStamps', 'pixelsPerM')
fade = repmat(linspace(0,1,s.edgeFading), vid.Height, 1);

% get kinematics in original pixel coordinates (kinData kinematics have
% been transformed st they cannot be directly overlaid on frames)
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM, 'scoreThresh', scoreThresh);
clear locationsTable

% choose trial(s)
if s.vertical
    dims = [vid.Height*4, vid.Width];
else    
    dims = [vid.Height, vid.Width*4 + s.imgSpacing*3];
end
pawOverSequence = nan(length(kinData),4);
pawOverSequence([kinData.isTrialAnalyzed],:) = [kinData.pawOverSequence]';
trials = find(ismember(pawOverSequence, pawSequence, 'rows') & ...
              [kinData.isTrialAnalyzed]'); % all trials with correct sequence of paws over obstacle

if ~ismember(trial, trials)
    disp('WARNING! Aborting because requested trial does not have the correct paw over sequence...');
else
    
    % select overlay trials
    trials = trials(trials~=trial); % ensure that imgTrial is not represented in the population of background trials
    if ~isempty(s.randSeed); rng(s.randSeed); end  % initialize random seed for reproduceibility
    trials = trials(randperm(length(trials), s.overlays));


    figure('name', sprintf('leadLagImg_%s_trial%i.png', session, trial), ...
        'menubar', 'none', 'Position', [2000 200 dims(2) dims(1)], 'InvertHardcopy', 'off');
    colormap gray;
    montage = uint8(zeros(dims));
    img = image(montage, 'cdatamapping', 'scaled'); hold on


    % create image montage
    scats = cell(1,4);
    for i = 1:4

        % add frame to image montage
        if s.vertical
            offset_x = 0;
            offset_y = (i-1)*vid.Height+1;
        else
            offset_x = (i-1)*(vid.Width+s.imgSpacing)+1;
            offset_y = 0;
        end
        featureBin = contains(features, ['paw' num2str(pawSequence(i))]) & contains(features, 'top');
        frameInd = find(kinData(trial).locations(:,1,pawSequence(i))<s.pawPos, 1, 'last'); % ind of frame within trial
        frameInd = kinData(trial).trialInds(frameInd); % ind of frame within entire video
        scats{i} = scatter(locations(frameInd,1,featureBin)+offset_x, locations(frameInd,2,featureBin)+offset_y, ...
            100, s.colors(i,:), 'filled', 'MarkerEdgeColor', 'white', 'LineWidth', 2); hold on

        frame = double(rgb2gray(read(vid, frameInd)));
        frame(:,1:s.edgeFading) = frame(:,1:s.edgeFading) .* fade; % left fade
        frame(:,end+1-s.edgeFading:end) = frame(:,end+1-s.edgeFading:end) .* fliplr(fade); % right fade
        
        montage(offset_y+1:offset_y+vid.Height, offset_x+1:offset_x+vid.Width) = uint8(frame);

        for j = [trial trials']
            % get inds of mod step start, mid mod step where the frame is taken, and mod step end
            stepBins = kinData(j).modifiedStepIdentities(:,pawSequence(i)); % 
            stepStart = find(stepBins==max(stepBins),1,'first');
            stepEnd = find(stepBins==max(stepBins),1,'last');

            % convert inds st they are relative to all frames, not just trial
            stepStart = kinData(j).trialInds(stepStart);
            stepEnd = kinData(j).trialInds(stepEnd);

            % plot kinematics
            x = locations(stepStart:stepEnd,1,featureBin)+offset_x;
            y = locations(stepStart:stepEnd,2,featureBin)+offset_y;
            if j==trial
                plot(x, y, 'LineWidth', s.mainWidth, 'Color', [s.colors(i,:) s.mainAlpha])
            else
                if s.scatter
                    scatter(x, y, 25, s.colors(i,:), 'filled', 'MarkerFaceAlpha', s.overlayAlpha)
                else
                    plot(x, y, 'LineWidth', s.overlayWidth, 'Color', [s.colors(i,:) s.overlayAlpha])
                end
            end
        end
    end
    uistack(cat(2, scats{:}), 'top')  % move scatter to top
    
    % adjust contast and update montage
    montage = imadjust(montage, s.contrastLims, [0 1]);
    set(img, 'CData', montage)

    % pimp fig
    set(gca, 'Visible', 'off', 'Position', [0 0 1 1])

    % save
    if s.vertical
        file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'leadLagImgs_vertical.png');
    else
        file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'leadLagImgs.png');
    end
    fprintf('writing %s to disk...\n', file)
    saveas(gcf, file)
end