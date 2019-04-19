% generates figures with raw tracking data for a single trial show, along
% with images overlaid on top of the traces that show snapshots of
% the behavior over time! so cool

% settings
session = '180703_000';
trialNum = 10;
pixStart = 700; % what pixel to start at
edgeFading = 60;
contrastLims = [.2 .8]; % pixels at these proportional values are mapped to 0 and 255
featuresToShow = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', ...
                  'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot'};
pawColors = jet(4);
imgs = 3; % how many images to overlay on figures
scatSize = 80;


%% initializations
vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions', 'wheelPositions', ...
    'wheelTimes', 'wheelCenter', 'wheelRadius', 'mToPixMapping');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data

%%
[locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, ...
    wheelTimes, wheelCenter, wheelRadius, 250, mToPixMapping(1));
fade = repmat(linspace(0,1,edgeFading), vidTop.Height+vidBot.Height, 1);
frameOffsets = round(-obsPixPositions' + max(obsPixPositions));
locations(:,1,:) = locations(:,1,:) + frameOffsets; % unravel coordinates to account for movement of wheel

% trials = sort(randperm(length(obsOnTimes), trialNum));
trials = 99;

for trial = trials
    startInd = find(frameTimeStamps>obsOnTimes(trial) & obsPixPositions'<pixStart, 1, 'first');
    endInd = find(frameTimeStamps<obsOffTimes(trial), 1, 'last');
    trialInds = startInd:endInd;
    imgOffsets = linspace(frameOffsets(trialInds(1)), frameOffsets(trialInds(end)), imgs); % evenly spaced image offsets
    imgInds = trialInds(knnsearch(frameOffsets(trialInds), imgOffsets'));
    dims = [vidBot.Height+vidTop.Height range(frameOffsets(trialInds))+vidTop.Width];
    allImgs = double(zeros(imgs, dims(1), dims(2)));



    % plot image
    for i = 1:length(imgInds)
        img = double(cat(1, rgb2gray(read(vidTop, imgInds(i))), rgb2gray(read(vidBot, imgInds(i)))));
        img(:,1:edgeFading) = img(:,1:edgeFading) .* fade; % left fade
        img(:,end+1-edgeFading:end) = img(:,end+1-edgeFading:end) .* fliplr(fade); % right fade
        allImgs(i,:,frameOffsets(imgInds(i))+1:frameOffsets(imgInds(i))+vidTop.Width) = double(img)*imgs;
    end
    
    close all;
    figure('Color', 'white', 'Position', [2000 442 dims(2) dims(1)], 'MenuBar', 'none')
    colormap gray
    montage = squeeze(mean(allImgs,1));
    montage = uint8(montage * (255/max(montage(:))));
    montage = imadjust(montage, contrastLims, [0 1]);
    
    % add obstacle
    obsInd = round(frameOffsets(imgInds(2)) + obsPixPositions(imgInds(2)));
    obsWidth = 5;
    montage(vidTop.Height+1:end, obsInd-obsWidth:obsInd+obsWidth) = 255;
    
    image(montage, 'cdatamapping', 'scaled'); hold on;
    
    
    
    

    % plot traces
    for i = 1:length(features)
        if ismember(features{i}, featuresToShow)
            if contains(features{i}, 'paw')
                pawNum = str2double(features{i}(4));
    %             inds = trialInds(~stanceBins(trialInds,pawNum)); % use this instead of subsequent line to hide kinematics during stance
                inds = trialInds;
                color = pawColors(pawNum,:);
            else
                inds = trialInds;
                color = [.8 .8 .8];
            end
            plot(locations(inds,1,i), locations(inds,2,i), ...
                'LineWidth', 2, 'Color', color); hold on
        end
    end

    % add scatters
    for i = 1:length(features)
        if ismember(features{i}, featuresToShow)
            if contains(features{i}, 'paw')
                color = pawColors(str2double(features{i}(4)),:);
            else
                color = [.8 .8 .8];
            end
            x = squeeze(locations(imgInds,1,i));
            y = squeeze(locations(imgInds,2,i));
            scatter(x,y, scatSize, color, 'filled')
        end
    end

    % pimp fig
    set(gca, 'DataAspectRatio', [1 1 1], ...
        'XLim', [frameOffsets(trialInds(1)) frameOffsets(trialInds(end))+vidTop.Width], ...
        'YLim', [0 (vidBot.Height+vidTop.Height)], ...
        'visible', 'off', 'color', 'black', ...
        'Position', [0 0 1 1]);

    % save
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'trackingOverFrames', [session '_trial_' num2str(trial) '.png']))
end
