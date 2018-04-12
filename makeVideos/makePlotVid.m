% function makePlotVid

% get trial
load([getenv('OBSDATADIR') 'kinematicData.mat'], 'data');
data = data([data.oneSwingOneStance] & ~[data.isFlipped]);
trialInd = randperm(length(data), 1);
session = data(trialInd).session;
trials = data(trialInd).trial;


%%

% settings 
contrastLims = [.1 .9];
paws = [2 3];

% initializations
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\stepSegmentation.mat'], 'modifiedStepIdentities')
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBotCorrected.mat'], 'locations')
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions')
locations = locations.locationsCorrected;
posRange = round(range(obsPixPositions));
dims = [vid.Height posRange+vid.Width];
obsPos = posRange;


%% prepare figure and objects
close all;
fig = figure('menubar', 'none', 'position', [1600 0 dims(2) dims(1)], 'color', 'black');

% frame
colormap gray
frame = zeros(vid.Height, vid.Width);
frameShow = image(1:vid.Width, 1:vid.Height, frame, 'cdatamapping', 'scaled'); hold on;

% obstacle
line([obsPos obsPos], [1 vid.Height], 'color', 'white', 'linewidth', 8);

% kinematic plots
rightPlot = plot(0,0, 'linewidth', 3);
leftPlot = plot(0,0, 'linewidth', 3);

% scatter points for end of kinematic trajectories
rightScatter = scatter(0,0);
leftScatter = scatter(0,0);



ax = gca;
set(ax, 'color', 'black', 'position', [0 0 1 1], 'xlim', [1 dims(2)], 'ylim', [1 vid.Height], 'visible', 'off', 'clim', [0 1]);



for i = 1:length(trials)
    
    frameBins = frameTimeStamps>=obsOnTimes(trials(i)) & frameTimeStamps<=obsOffTimes(trials(i));
    frameInds = find(frameBins);
    
    for j = 1:length(frameInds)
        
        % update frame
        frame = rgb2gray(read(vid, frameInds(j)));
        frame = double(frame) / 255;
        frame = imadjust(frame, contrastLims, [0 1]);
        leftInd = round(obsPos - obsPixPositions(frameInds(j)));
        set(frameShow, 'XData', (1:vid.Width)+leftInd, 'CData', frame);
        
        % update trajectories
        validBins = frameBins & (1:length(frameTimeStamps))'<=frameInds(j); % bins in trial up to and including current frame
        rightBins = modifiedStepIdentities(:,3)==1 & validBins;
        leftBins = modifiedStepIdentities(:,2)==1 & validBins;
        
        
        
        if any(rightBins)
            rightLocations = locations(rightBins,:,3);
            rightLocations(:,1) = rightLocations(:,1) - obsPixPositions(rightBins)' + obsPos;
            set(rightPlot, 'XData', rightLocations(:,1), 'YData', rightLocations(:,2))
        end
        
        if any(leftBins)
            leftLocations = locations(leftBins,:,2);
            leftLocations(:,1) = leftLocations(:,1) - obsPixPositions(leftBins)' + obsPos;
            set(leftPlot, 'XData', leftLocations(:,1), 'YData', leftLocations(:,2))
        end
        
        pause(.05)
        
        
    end
    
end


% close(fig)


















