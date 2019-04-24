function showLeadingLaggingImg(session, imgTrial)

% make img showing four pics of mouse getting over obstacle, show leading
% fore, lagging fore, leading hind, then lagging hind // overlay kinematic
% traces for one paw per img, as well as randomly selected trials
% transparently overlaid in background // randomly selects valid trial
% for the main image, and randomly selects among valid trials the
% background traces // the file name reflects the randomly selected trial
% for the main image // if imgTrial is supplied by user, it forces this trial to be the
% background image

% settings
% session = '190318_000';
% trial = 74;
pawPos = .005; % (m) each frame is selected when the paw of interest is pawPos in front of obs (can be negative too)
imgSpacing = 0; % pixels for separating images
edgeFading = 60; % fading at the edges of images
contrastLims = [.05 .8]; % pixels at these proportional values are mapped to 0 and 255
colors = hsv(4);

overlays = 8; % how many trial kinematics to overlay
overlayAlpha = .2;
overlayWidth = 2;
mainWidth = 4;
mainAlpha = 1;


% initializations
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
pawSequence = [3 2 4 1]; % only include trials with this sequence of paws going over
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), 'kinData')
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'frameTimeStamps')
fade = repmat(linspace(0,1,edgeFading), vid.Height, 1);

% get kinematics in original pixel coordinates (kinData kinematics have
% been transformed st they cannot be directly overlaid on frames)
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
clear locationsTable

% choose trial(s)
dims = [vid.Height, vid.Width*4 + imgSpacing*3];
pawOverSequence = nan(length(kinData),4);
pawOverSequence([kinData.isTrialAnalyzed],:) = [kinData.pawOverSequence]';
trials = find(ismember(pawOverSequence, pawSequence, 'rows') & ...
              [kinData.isTrialAnalyzed]'); % all trials with correct sequence of paws over obstacle
if exist('imgTrial', 'var')
    trials = trials(trials~=imgTrial); % ensure that imgTrial is not represented in the population of background trials
else
    imgTrial = trials(randperm(length(trials),1));
end
trials = trials(trials~=imgTrial); % ensure that imgTrial is not represented in the population of background trials
trials = trials(randperm(length(trials), overlays));





figure('menubar', 'none', 'Position', [2000 600 dims(2) dims(1)]);
colormap gray;
montage = uint8(zeros(dims));
img = image(montage, 'cdatamapping', 'scaled'); hold on


% create image montage
for i = 1:4

    % add frame to image montage
    offset = (i-1)*(vid.Width+imgSpacing)+1;
    featureBin = contains(features, ['paw' num2str(pawSequence(i))]) & contains(features, 'top');
    frameInd = find(kinData(imgTrial).locations(:,1,pawSequence(i))<pawPos, 1, 'last'); % ind of frame within trial
    frameInd = kinData(imgTrial).trialInds(frameInd); % ind of frame within entire video
    scatter(locations(frameInd,1,featureBin)+offset, locations(frameInd,2,featureBin), ...
        150, colors(i,:), 'filled')

    frame = double(rgb2gray(read(vid, frameInd)));
    frame(:,1:edgeFading) = frame(:,1:edgeFading) .* fade; % left fade
    frame(:,end+1-edgeFading:end) = frame(:,end+1-edgeFading:end) .* fliplr(fade); % right fade
    montage(:,offset:offset+vid.Width-1) = uint8(frame);

%         keyboard
    for trial = [imgTrial trials']
        % get inds of mod step start, mid mod step where the frame is taken, and mod step end
        stepBins = kinData(trial).modifiedStepIdentities(:,pawSequence(i)); % 
        stepStart = find(stepBins==max(stepBins),1,'first');
        stepEnd = find(stepBins==max(stepBins),1,'last');

        % convert inds st they are relative to all frames, not just trial
        stepStart = kinData(trial).trialInds(stepStart);
        stepEnd = kinData(trial).trialInds(stepEnd);

        % plot kinematics
        if trial==imgTrial    
            plot(locations(stepStart:stepEnd,1,featureBin)+offset, locations(stepStart:stepEnd,2,featureBin), ...
                'LineWidth', mainWidth, 'Color', [colors(i,:) mainAlpha])
        else
            plot(locations(stepStart:stepEnd,1,featureBin)+offset, locations(stepStart:stepEnd,2,featureBin), ...
                'LineWidth', overlayWidth, 'Color', [colors(i,:) overlayAlpha])
        end
    end
end

% adjust contast and update montage
montage = imadjust(montage, contrastLims, [0 1]);
set(img, 'CData', montage)

% pimp fig
set(gca, 'Visible', 'off', 'Position', [0 0 1 1])

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', ...
    sprintf('leadLagImg_%s_trial%i.png', session, imgTrial));
fprintf('writing %s to disk...\n', file)
saveas(gcf, file)