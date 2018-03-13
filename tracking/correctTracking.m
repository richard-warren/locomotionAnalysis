function correctTracking(outputFile, vid, locations, frameInds, vidDelay, anchorPts, paws)


% settings
circSize = 200;
vidSizeScaling = 1.25;
colors = hsv(4);
maxCorrectionInterpolation = 20; % max number of frames that can be interpolated between locationCorrections


% initializations
pawPair = [4 3];
playing = true;
paused = false;
frameUpdating = true;
currentFrameInd = 1;
lastCorrectedPaw = [];
sampleFrame = rgb2gray(read(vid,currentFrameInd));
w = waitbar(0, 'correction progress...', 'position', [1500 50 270 56.2500]);


% initialize fields for manual corrections (if they do not already exist)
if ~isfield(locations, 'locationCorrections'); locations.locationCorrections = nan(size(locations.locationsRaw)); end
if ~isfield(locations, 'locationsCorrected'); locations.locationsCorrected = nan(size(locations.locationsRaw)); end
if ~isfield(locations, 'isCorrected'); locations.isCorrected = false(length(locations.isAnalyzed),1); end

% initialize corrected data (for every trial, replace x and y entries with non-nan entries in xCorrections and yCorrections // then interpolate on a trial by trial basis)
for i = unique(locations.trialIdentities(~isnan(locations.trialIdentities)))'
    mergeCorrectionsAndInterpolate(i, 1:4);
end



% prepare figure
fig = figure('name', outputFile, 'units', 'pixels', 'position', [600 400 vid.Width*vidSizeScaling vid.Height*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);

colormap gray
preview = image(sampleFrame, 'CDataMapping', 'scaled'); hold on;
rawAxis = gca;
set(rawAxis, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vid.Width*vidSizeScaling vid.Height*vidSizeScaling]);
circSizes = circSize * ones(1,length(anchorPts));


% % prepare lines showing x locations of bottom tracked paws
% lines = cell(4,1);
% if exist('lineLocations', 'var')
%     for i = 1:4
%         lines{i} = line([0 0], [vid.Height vid.Height-50], 'color', colors(i,:));
%     end
% end


% prepare circles showing locationsRaw
scatterRaw = scatter(rawAxis, zeros(1,length(anchorPts)), zeros(1,length(anchorPts)),...
    circSizes, colors, 'linewidth', 2.5); hold on


% prepare colored circles in the corner to tell you which paw goes where
hold on; scatter(rawAxis, [anchorPts{1}(1) anchorPts{2}(1) anchorPts{3}(1) anchorPts{4}(1)] .* (vid.Width-1) + 1,...
                 [anchorPts{1}(2) anchorPts{2}(2) anchorPts{3}(2) anchorPts{4}(2)] .* (vid.Height-1) + 1,...
                 circSizes, colors, 'filled', 'linewidth', 3);     % show anchor points

% prepare impoints, draggable markers used to show / adjust tracking
for i = 1:4
    impoints{i} = impoint(gca, [10 10]*i);
    setColor(impoints{i}, colors(i,:));
    addNewPositionCallback(impoints{i}, @(x) addCorrection(i));
end


% main loop
while playing
    while paused; pause(.001); end
    updateFrame(1);
end
close(fig)
close(w)

% ---------
% FUNCTIONS
% ---------

% keypress controls
function changeFrames(~,~)
    
    key = double(get(fig, 'currentcharacter'));
    
    switch key
        
        % LEFT: move frame backward
        case 28
            pause(.001);
            paused = true;
            updateFrame(-1);
        
        % RIGHT: move frame forward
        case 29
            pause(.001);
            paused = true;
            updateFrame(1);
        
        % 'f': select frame
        case 102
            pause(.001);
            paused = true;
            input = inputdlg('enter frame number');
            currentFrameInd = find(frameInds == str2double(input{1}));
            updateFrame(1);
        
        % ESCAPE: close window
        case 27
            playing = false;
            paused = false;
            
        % 's': save current progress
        case 115
            save(outputFile, 'locations')
            m = msgbox('saving!!!'); pause(.5); close(m)
            
        % 'c': go to latest corrected frame
        case 99
            pause(.001);
            paused = true;
            currentFrameInd = find(frameInds == find(locations.isCorrected,1,'last'));
            msgbox('moved to latest corrected frame');
            updateFrame(0);
            
        % 'i': for last corrected paw, interpolate between last two corrected values
        case 105
            
            % find last two corrected inds for lastCorrectedPaw
            lastTwoCorrectedInds = find(~isnan(locations.locationCorrections(1:frameInds(currentFrameInd),1,lastCorrectedPaw)), 2, 'last');
            
            interpolatedTrials = locations.trialIdentities(lastTwoCorrectedInds);
            
            if diff(interpolatedTrials)==0 && diff(lastTwoCorrectedInds)<maxCorrectionInterpolation % if inds belong to the same trial and they arenot too far apart
            
                inds = lastTwoCorrectedInds(1):lastTwoCorrectedInds(2);
                
                % interpolate x and y valuesinds = lastTwoCorrectedInds(1):lastTwoCorrectedInds(2); between last two corrected points and store in locationsCorrected
                for j = 1:2
                    xInterp = fillmissing((locations.locationCorrections(inds,j,lastCorrectedPaw)), 'linear');
                    locations.locationCorrections(inds,j,lastCorrectedPaw) = xInterp;
                end

                % recorrect and interpolate the trial for lastCorrectedPaw
                mergeCorrectionsAndInterpolate(interpolatedTrials(1), lastCorrectedPaw)
                
                % rewind to first interpolated frame
                for j = 1 : ((frameInds(currentFrameInd) - lastTwoCorrectedInds(1)) + 5)
                    updateFrame(-1)
                    pause(vidDelay)
                end
            end
            
        % 'p': user enters the paw pair to be swapped
        case 112
            pause(.001);
            paused = true;
            input = inputdlg('enter paw pair number');
            pawPair(1) = str2num(input{1}(1));
            pawPair(2) = str2num(input{1}(2));
            
        % open bracket: flip pawPair
        case 91
            pause(.001);
            paused = true;
            locations.locationCorrections(frameInds(currentFrameInd),:,pawPair) = ...
                locations.locationsCorrected(frameInds(currentFrameInd), :, fliplr(pawPair));
            mergeCorrectionsAndInterpolate(locations.trialIdentities(frameInds(currentFrameInd)), pawPair);
            updateFrame(0);
            
        
        % SPACEBAR: pause playback
        case 32
            paused = ~paused;
    end
end



% update frame preview
function updateFrame(frameStep)
    
    currentFrameInd = currentFrameInd + frameStep;
    if currentFrameInd > length(frameInds); currentFrameInd = length(frameInds);
    elseif currentFrameInd < 1; currentFrameInd = 1; end
    
    % record that this frame has been corrected (somebody looked at it and verified it was correct or corrected it manually)
    locations.isCorrected(frameInds(currentFrameInd)) = true;
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid, frameInds(currentFrameInd)));
    
    
    % add vertical lines
%     if exist('lineLocations', 'var')
%         inds = lineLocations.x(frameInds(currentFrameInd),:);
%         for j = 1:4
%             set(lines{j}, 'XData', [inds(j) inds(j)])
%         end
%     end
    
    
    % add frame number
    frame = insertText(frame, [size(frame,2) size(frame,1)], ...
        sprintf('trial %i, frame %i', locations.trialIdentities(frameInds(currentFrameInd)), frameInds(currentFrameInd)),...
	    'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
    
    % update figure
    set(preview, 'CData', frame);
    
    
    % update point locations to locationsCorrected
    frameUpdating = true;
    for j = 1:4
        setPosition(impoints{j}, locations.locationsCorrected(frameInds(currentFrameInd),1,j), ...
            locations.locationsCorrected(frameInds(currentFrameInd),2,j));
    end
    frameUpdating = false;
    
    % update circles to locationsRaw (change this so it updates based on locationsCorrected but only when there is nan both in locationsCorrections and locationsRaw)
    frameLocations = locations.locationsRaw(frameInds(currentFrameInd),:,:);
    frameLocationCorrections = locations.locationCorrections(frameInds(currentFrameInd),:,:);
    correctedBins = ~isnan(frameLocationCorrections);
    frameLocations(correctedBins) = frameLocationCorrections(correctedBins);
    frameLocations = squeeze(frameLocations);
    
    set(scatterRaw, 'XData', frameLocations(1,:), 'YData', frameLocations(2,:), 'visible', 'on');
    
    % pause to reflcet on the little things...
    pause(vidDelay);
    
    waitbar(currentFrameInd / length(frameInds))
end





% adds corrected values to raw values for a trial, and interpolates missing values
function mergeCorrectionsAndInterpolate(trial, paws)
    
    trialInds = ([locations.trialIdentities]==trial);
    
    for dimension = 1:2
        for paw = paws
            
            % extract locations for a single trial, paw, and dimension (x, y, or z)
            trialLocations = squeeze(locations.locationsRaw(trialInds, dimension, paw));
            trialCorrections = squeeze(locations.locationCorrections(trialInds, dimension, paw));
            
            % if there are any manually corrected values, incorporate them
            correctedInds = ~isnan(trialCorrections);
            if ~isempty(correctedInds); trialLocations(correctedInds) = trialCorrections(correctedInds); end
            
            % fill in missing values and store in locationsCorrected field
            trialLocations = fillmissing(trialLocations, 'spline', 'endvalues', 'none');
%             if exist('smoothing', 'var'), trialLocations = smooth(trialLocations, smoothing); end
            locations.locationsCorrected(trialInds, dimension, paw) = trialLocations;
        end
    end
end


function addCorrection(paw)
    
    if ~frameUpdating
        trial = locations.trialIdentities(frameInds(currentFrameInd));
        locations.locationCorrections(frameInds(currentFrameInd), :, paw) = getPosition(impoints{paw});
        mergeCorrectionsAndInterpolate(trial, paw);
        lastCorrectedPaw = paw;
        updateFrame(0);
    end
end
    
    
end













