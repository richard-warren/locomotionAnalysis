function getLocations(session, steps, minVel, showTracking)

% performs paw tracking for a single session!
% if makeVid is true, it also automatically makes a video of the tracking



% settings
% obsPrePost = [.2 .2]; % include this many meters before and after obs turns on
obsPrePostBot = [.55 .25];
obsPrePostTop = [.55 .25];
velPrePost = [.1 .1]; % compute trials velocity between these obstacle positions (relative to tip of mouse's nose)
anchorPtsBot = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner
colors = hsv(4); % red green blue purple


% initializations
trackingDir = [getenv('OBSDATADIR') 'sessions\' session '\tracking'];
if ~exist(trackingDir, 'dir'); mkdir(trackingDir); end % make tracking directory if it doesn't already exist

xMapping = [getenv('GITDIR') 'locomotionAnalysis\xAlignment\xLinearMapping.mat'];
load(xMapping, 'xLinearMapping');

load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps', 'nosePos', ...
    'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'obsOnTimes', 'obsOffTimes', 'mToPixMapping', 'targetFs')
obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));

[frameIndsBot, trialIdentities] = getTrialFrameInds(minVel, obsPrePostBot, velPrePost, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);
frameIndsTop = getTrialFrameInds(minVel, obsPrePostTop, velPrePost, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);
% fprintf('%s: analyzing %i sessions', session, length(unique(trialIdentities(~isnan(trialIdentities)))));

mToPixFactor = median(mToPixMapping(:,1)); % get mapping from meters to pixels

vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
wheelPoints = getWheelPoints(vidTop);


% potential locations bot
if any(strcmp('potBot', steps))
    
    % settings
    scoreThresh = 0;
    
    % svm1
    load([getenv('OBSDATADIR') 'tracking\classifiers\pawBot1'], 'model', 'subFrameSize');
    model1 = model; clear model; subFrameSize1 = subFrameSize; clear subFrameSize;
    load([getenv('OBSDATADIR') 'tracking\classifiers\pawBot2AlexNet (fore vs hind).mat'], 'convNetwork', 'subFrameSize');
    model2 = convNetwork; clear convNetwork; subFrameSize2 = subFrameSize; clear subFrameSize;
    classNum = model2.Layers(end).OutputSize - 1; % not including NOT PAW class
    
    % find paws
    tic; potentialLocationsBot = getPotentialLocationsBot(vidBot, model1, model2, classNum, ...
        subFrameSize1, subFrameSize2, scoreThresh, obsPixPositions, frameIndsBot, trialIdentities, showTracking);

    % save results
    save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');
    save([getenv('OBSDATADIR') 'sessions\' session '\tracking\velocityInfo.mat'], 'trialVels', 'minVel', 'velPositions');
    fprintf('%s: potentialLocationsBot analyzed in %.1f minutes\n', session, (toc/60))
end


% locations bot
if any(strcmp('bot', steps))
    
    % get bot locations
    showPotentialLocations = 1;
    
    % initializations
    if ~exist('potentialLocationsBot', 'var') % load potentialLocations if it was not already computed above
        load([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsBot.mat'])
    end
    

    locations = getLocationsBot(potentialLocationsBot, anchorPtsBot, frameTimeStamps, vidBot.Width, vidBot.Height);
    save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'], 'locations');
    fprintf('%s: locationsBot analyzed\n', session)
    
    if showTracking
        showLocations(vidBot, find(locations.isAnalyzed), potentialLocationsBot, ...
            locations.locationsRaw, locations.trialIdentities, 0, .02, anchorPtsBot, hsv(4));
    end
end


% stance analysis
if any(strcmp('stance', steps))
    
    % initializations
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBotCorrected.mat'], 'locations')
    fs = 250;
    
    % fix x alignment for bottom view
    locationsBot = locations.locationsCorrected;
    locationsBot(:,1,:) = locationsBot(:,1,:)*xLinearMapping(1) + xLinearMapping(2);
    
    % get stance bins
    stanceBins = getStanceBins(vidTop, wheelPoints, squeeze(locationsBot(:,1,:)), locations.trialIdentities', ...
        fs, mToPixFactor, wheelPositions, wheelTimes, targetFs, frameTimeStamps);
    
    save([getenv('OBSDATADIR') 'sessions\' session '\tracking\stanceBins.mat'], 'stanceBins')

end


% potential locations top
if any(strcmp('potTop', steps))
    
    % get potential locations for top (svm)
    scoreThresh = 0;

    % initializations
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBotCorrected.mat'], 'locations')
    load([getenv('OBSDATADIR') 'tracking\classifiers\pawTop1'], 'model', 'subFrameSize');
    model1 = model; clear model; subFrameSize1 = subFrameSize; clear subFrameSize;
    load([getenv('OBSDATADIR') 'tracking\classifiers\pawTop2AlexNet'], 'convNetwork', 'subFrameSize')
    model2 = convNetwork; clear convNetwork; subFrameSize2 = subFrameSize;
    classNum = model2.Layers(end).OutputSize - 1; % not including NOT PAW class


    tic;
    potentialLocationsTop = getPotentialLocationsTop(vidTop, locations, model1, model2, ...
        classNum, subFrameSize1, subFrameSize2, scoreThresh, frameIndsTop, showTracking);
    save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');
    fprintf('%s: potentialLocationsTop analyzed in %.1f minutes\n', session, toc/60)
end


% locations top
if any(strcmp('top', steps))
    
    % settings
    showPotentialLocations = 1;
    fs = 250;
    
    % initializations
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsTop.mat'], 'potentialLocationsTop')
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBotCorrected.mat'], 'locations')
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\stanceBins.mat'], 'stanceBins')

    % fix x alignment for bottom view
    locationsBot = locations.locationsCorrected;
    locationsBot(:,1,:) = locationsBot(:,1,:)*xLinearMapping(1) + xLinearMapping(2);

    locations = getLocationsTop(potentialLocationsTop, locationsBot, stanceBins, ...
        wheelPositions, wheelTimes, targetFs, mToPixFactor, obsPixPositions, frameTimeStamps, fs, wheelPoints);
    save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsTop.mat'], 'locations');
    fprintf('%s: locationsTop analyzed\n', session)
    
    if showTracking
        showLocations(vidTop, find(locations.isAnalyzed), potentialLocationsTop, ...
            locations.locationsRaw, locations.trialIdentities, showPotentialLocations, ...
            .04, anchorPtsBot, hsv(4), locationsBot);
    end
end



