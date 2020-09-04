function makePlotVid


% settings
trialNum = 1;
pixelsLeftOfPos = 1000;
contrastLims = [.1 .9];
paws = [2 3];
circSize = 100;
obsPixPosBuffer = 40; % replace the first obsPixPosBuffer frames with positions inferred from wheel positions - this is to avoid the duration of veloivty ramp
fps = 250;
playbackSpeed = .1;
postWiskSlowDown = 2; % slow down by this factor after whisk touches obs
wiskPause = 3;



% find trials
load([getenv('OBSDATADIR') 'kinematicData.mat'], 'data');

numModSteps = cellfun(@(x) x(1,3), {data.modStepNum});
predictedDistances = [data.swingStartDistance] + [data.predictedLengths]; % predicted distance to obs
deltaLengths = cellfun(@(x) x(1,3), {data.modifiedSwingLengths}) - [data.predictedLengths];

inds = find([data.oneSwingOneStance] & ~[data.isFlipped] & numModSteps~=1 & ...
    predictedDistances>-.005 & predictedDistances<-.0025 & deltaLengths<-.01);
% inds = find([data.oneSwingOneStance] & ~[data.isFlipped] & numModSteps~=1 & deltaLengths>-.002 & deltaLengths<.002);
% inds = find([data.oneSwingOneStance] & ~[data.isFlipped] & numModSteps==1 & deltaLengths>.02 & predictedDistances>.002);

trialInds = inds(randperm(length(inds), trialNum));
sessions = {data(trialInds).session};
[sessions, sortInds] = sort(sessions);
trialInds = trialInds(sortInds);





% initializations
vidWriter = VideoWriter([getenv('OBSDATADIR') 'editedVid\plotVid.mp4'], 'MPEG-4');
frameRate = round(playbackSpeed*fps);
set(vidWriter, 'FrameRate', frameRate);
open(vidWriter);
colors = [.25 1 1; .25 1 .25]; % use these if you want different colors per step type (lengthened or shortened)
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{1} '\runBot.mp4']);
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{1} '\runTop.mp4']);
dims = [vid.Height+vidTop.Height pixelsLeftOfPos+vid.Width];


% prepare figure and objects
fig = figure('menubar', 'none', 'position', [2000 100 dims(2) dims(1)], 'color', 'black');

% frame
colormap gray
frame = zeros(vid.Height, vid.Width);
frameTop = zeros(vidTop.Height, vidTop.Width);
frameShow = image(1:vid.Width, 1:vid.Height+vidTop.Height, cat(1, frame, frameTop), 'cdatamapping', 'scaled'); hold on;


% kinematic plots
plots = cell(1,length(paws));
for i = 1:length(plots)
    plots{i} = plot(0,0, 'linewidth', 3);
end

% scatter points for end of kinematic trajectories
scatters = cell(1,length(paws));
for i = 1:length(scatters)
    scatters{i} = scatter(0,0, circSize, colors(i,:), 'filled');
end

% obstacle
line([pixelsLeftOfPos pixelsLeftOfPos], [vidTop.Height vidTop.Height+vid.Height], 'color', 'white', 'linewidth', 8);

% predicted distance line
predLine = line([0 0], [vidTop.Height vidTop.Height+vid.Height], 'color', 'white', 'linewidth', 3, 'linestyle', ':', 'visible', 'off');


ax = gca;
set(ax, 'color', 'black', 'position', [0 0 1 1], 'xlim', [1 dims(2)], 'ylim', [1 vid.Height+vidTop.Height], 'visible', 'off', 'clim', [0 1]);



for i = 1:length(trialInds)
    
    session = data(trialInds(i)).session;
    trial = data(trialInds(i)).trial;
    
    % load session data if it is not already loaded
    if i==1 || ~strcmp(sessions{i-1}, sessions{i})
        vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
        vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
        load([getenv('OBSDATADIR') 'sessions\' session '\tracking\stepSegmentation.mat'], 'modifiedStepIdentities')
        load([getenv('OBSDATADIR') 'sessions\' session '\wiskContactData.mat'], 'contactTimes')
        load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBotCorrected.mat'], 'locations')
        load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
            'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions', 'wheelPositions', 'wheelTimes', 'mToPixMapping')
        mToPixMapping = median(mToPixMapping,1);
        wheelPosInterp = interp1(wheelTimes, wheelPositions, frameTimeStamps);
        locations = locations.locationsCorrected;
    end
    
    % reset graphics objects
    set(predLine, 'visible', 'off')
    for j=1:length(paws)
        set(plots{j}, 'XData', 0, 'YData', 0)
        set(scatters{j}, 'XData', 0, 'YData', 0)
    end

    
    % get frame info
    frameBins = frameTimeStamps>=obsOnTimes(trial) & frameTimeStamps<=obsOffTimes(trial);
    frameInds = find(frameBins);
    numModSteps = max(modifiedStepIdentities(frameInds,:), [], 1); % number of mod steps for each of 4 paws
    contactInd = find(frameTimeStamps>=contactTimes(trial), 1, 'first');
    
    % get final swing inds
    lastSwingInds = nan(1,length(paws));
    for k = 1:length(paws)
        lastSwingInds(k) = find(modifiedStepIdentities(:,paws(k))==1 & frameBins, 1, 'last');
    end
    
    % replace nans at edges of obsPixPositions by figuring out linear mapping for wheelPositions to obsPixPos
    startInd = find(isnan(obsPixPositions) & frameTimeStamps'>obsOffTimes(trial-1), 1, 'first');
    endInd = find(isnan(obsPixPositions) & frameTimeStamps'<obsOnTimes(trial+1), 1, 'last');
    pixPosStart = find(frameBins & ~isnan(obsPixPositions)', 1, 'first') + obsPixPosBuffer;
    pixPosEnd = find(frameBins & ~isnan(obsPixPositions)', 1, 'last');
    fit = polyfit(wheelPosInterp(pixPosStart:pixPosEnd), obsPixPositions(pixPosStart:pixPosEnd)', 1);
    predictedObsPixPositions = wheelPosInterp*fit(1) + fit(2);
    
    
    % recompute frame inds to restrict from pixelsLeftOfPos to moment mouse exits screen
    startInd = find((pixelsLeftOfPos - predictedObsPixPositions+vid.Width)>0 & frameTimeStamps>obsOffTimes(trial-1), 1, 'first') - 1;
    endInd = find(pixelsLeftOfPos - predictedObsPixPositions>dims(2),1,'first'); % extend frameInds so they go until mouse is off the screen
    frameInds = startInd:endInd;
    
    % iterate through all frames
    for j = 1:length(frameInds)
        
        % update frame
        frame = rgb2gray(read(vid, frameInds(j)));
        frame = double(frame) / 255;
        frame = imadjust(frame, contrastLims, [0 1]);
        
        frameTop = rgb2gray(read(vidTop, frameInds(j)));
        frameTop = double(frameTop) / 255;
        frameTop = imadjust(frameTop, contrastLims, [0 1]);
        
        
        frameLeftInd = round(pixelsLeftOfPos - predictedObsPixPositions(frameInds(j)));
        set(frameShow, 'XData', (1:vid.Width)+frameLeftInd, 'CData', cat(1, frameTop, frame));
        
        
        for k = 1:length(paws)
            
            locationBins = modifiedStepIdentities(:,paws(k))==1 & frameBins & (1:length(frameTimeStamps))'<=frameInds(j);

            if any(locationBins)
                
                % update plot
                pawLocations = locations(locationBins,:,paws(k));
                pawLocations(:,2) = pawLocations(:,2) + vidTop.Height; % add y offset so traces don't end up on top view
                pawLocations(:,1) = pawLocations(:,1) - obsPixPositions(locationBins)' + pixelsLeftOfPos;
                if numModSteps(paws(k))==1; trialColor = colors(2,:); else; trialColor = colors(1,:); end
                set(plots{k}, 'XData', pawLocations(:,1), 'YData', pawLocations(:,2), 'color', trialColor)
                
                % update scatter
                if frameInds(j)==lastSwingInds(k)
                    set(scatters{k}, 'XData', pawLocations(end,1), ...
                        'YData', pawLocations(end,2), 'CData', trialColor)
                end
            end
        end
        
        if frameInds(j)==contactInd
            predictedPos = pixelsLeftOfPos + round(predictedDistances(trialInds(i)) * abs(mToPixMapping(1)));
            set(predLine, 'XData', [predictedPos predictedPos], 'visible', 'on')
            for k = 1:wiskPause*frameRate-1; writeVideo(vidWriter, getframe(gcf)); end
        end
        
        % write frames
        if frameInds(j)>=contactInd
            for k = 1:postWiskSlowDown
                writeVideo(vidWriter, getframe(gcf));
            end
        else
            writeVideo(vidWriter, getframe(gcf));
        end
    end
end


close(fig)
close(vidWriter)






