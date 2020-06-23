function analyzeSession(session, varargin)

    % performs preliminary analyses on a session, including turning raw
    % data in run.mat into useful data in runAnalyzed.mat // also has
    % wrappers that execute neural network algorithms // variable names are
    % listed as arguments to the analyzeVar() function in the code blocks
    % below // to recompute a variable, pass 'overwriteVars' as follows:
    %
    % analyzeSession('191118_001', 'overwriteVars', 'rewardTimes')
    % analyzeSession('191118_001', 'overwriteVars', {'rewardTimes, isLightOn'}) % for multiple variables
    
    
    % settings
    s.verbose = true;
    s.superVerbose = false;  % whether to also display 
    s.targetFs = 1000; % frequency that positional data will be resampled to
    s.overwriteVars = '';
    s.plotObsTracking = false;  % whether to check obstacle tracking of wheel velocity by plotting them on top of one another
    
    s.rerunRunNetwork = false;
    s.rerunWiskNetwork = false;
    s.rerunWiskContactNetwork = false;
    s.rerunPawContactNetwork = false;
    
    s.showLickFig = false;  % whether to show figure for lickTimes analysis

    % rig characteristics
    s.whEncoderSteps = 2880; % 720cpr * 4t
    s.wheelRad = 95.25; % mm
    s.obEncoderSteps = 1000; % 250cpr * 4
    s.obsPulleyRad = 96 / (2*pi); % radius of timing pulley driving belt of obstacles platform
    s.obsDiameter = 3.175; % (mm)

    % initializations
    if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
    if ischar(s.overwriteVars); s.overwriteVars = {s.overwriteVars}; end  % make sure overwriteVars is in cell format
    anythingAnalyzed = false; % results are only saved if something was found that wasn't already analyzed
    sessionDir = fullfile(getenv('OBSDATADIR'), 'sessions', session);
    isRunVid = exist(fullfile(sessionDir, 'run.mp4'), 'file');
    isWiskVid = exist(fullfile(sessionDir, 'runWisk.mp4'), 'file');

    % load or initialize data structure
    if exist(fullfile(sessionDir, 'runAnalyzed.mat'), 'file') && ~isequal(s.overwriteVars, {'all'})
        data = load(fullfile(sessionDir, 'runAnalyzed.mat'));
    else
        data = struct();
    end
    computedVars = fieldnames(data);


    
    
    % run tracking
    if (~exist(fullfile(sessionDir, 'trackedFeaturesRaw.csv'), 'file') || s.rerunRunNetwork) && isRunVid
        dpkAnalysis(session, 'verbose', s.superVerbose, ...
            'vid', 'run.mp4', ...
            'model', 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_StackedDenseNet.h5', ...
            'skeleton', 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv', ...
            'output', 'trackedFeaturesRaw.csv')
    end
    
    
    
    
    % face tracking
    if (~exist(fullfile(sessionDir, 'trackedFeaturesRaw_wisk.csv'), 'file') || s.rerunWiskNetwork) && isWiskVid
        dpkAnalysis(session, 'verbose', s.superVerbose, ...
            'vid', 'runWisk.mp4', ...
            'model', 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_wisk_StackedDenseNet.h5', ...
            'skeleton', 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv', ...
            'output', 'trackedFeaturesRaw_wisk.csv')
    end
    
    
    
    
    % whisker angle
    if analyzeVar('whiskerAngle') && exist(fullfile(sessionDir, 'trackedFeaturesRaw_wisk.csv'), 'file')
        
        if s.verbose; fprintf('%s: getting whisker angles\n', session); end
        if ~exist('locationsWisk', 'var'); locationsWisk = readtable(fullfile(sessionDir, 'trackedFeaturesRaw_wisk.csv')); end
        
        % settings
        confidenceThresh = .5;
        medianFiltering = 5;
        smoothing = 5;
        angleLims = [-120 -60];

        pad = [median(locationsWisk.wisk_pad) median(locationsWisk.wisk_pad_1)];  % location of whisker pad
        
        valid = locationsWisk.wisk_caudal_2>confidenceThresh | locationsWisk.wisk_rostral_2>confidenceThresh;

        % get [x y] matrices for rostral and caudal whiskers
        c = [medfilt1(locationsWisk.wisk_caudal,medianFiltering), medfilt1(locationsWisk.wisk_caudal_1,medianFiltering)];
        r = [medfilt1(locationsWisk.wisk_rostral,medianFiltering), medfilt1(locationsWisk.wisk_rostral_1,medianFiltering)];

        % average position of rostral and caudal whisker
        avg = mean(cat(3,c,r), 3);
        avg = avg - pad;  % normalize to location of whisker pad
        avg(:,2) = -avg(:,2);  % flip y axis (because y is reversed in tracking)
        
        % compute angles
        whiskerAngle = rad2deg(atan2(avg(:,2), avg(:,1)));
        whiskerAngle = smooth(whiskerAngle, smoothing);
        whiskerAngle(whiskerAngle<angleLims(1) | whiskerAngle>angleLims(2) | ~valid) = nan;
        
        whiskerAngle = fillmissing(whiskerAngle, 'linear');
        saveVars('whiskerAngle', whiskerAngle)
    end
    
    
    
    % decode wheel position
    if analyzeVar('wheelPositions', 'wheelTimes')

        if s.verbose; fprintf('%s: decoding wheel position\n', session); end
        load(fullfile(sessionDir, 'run.mat'), 'whEncodA', 'whEncodB')

        [wheelPositions, wheelTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                     whEncodB.times, whEncodB.level,...
                                                     s.whEncoderSteps, s.wheelRad, s.targetFs, session);
        
         saveVars('wheelPositions', wheelPositions, 'wheelTimes', wheelTimes, 'targetFs', s.targetFs);
    end
    
    
    
    % analyze reward times
    if analyzeVar('rewardTimes', 'isRewardSurprise', 'omissionTimes')
        
        if s.verbose; fprintf('%s: getting reward times\n', session); end
        load(fullfile(sessionDir, 'run.mat'), 'reward', 'Cue')
        
        % settings
        minRewardInteveral = 1;  % remove rewards detected within minRewardInterval of eachother
        
        % reward times
        if isfield(reward, 'values')  % if recorded as analog input (sessions prior to 191523)
            rewardInds = find(diff(reward.values>2)==1) + 1;
            rewardTimes = reward.times(rewardInds);
        else
            rewardTimes = reward.times(logical(reward.level));  % keep only transitions from low to high
        end
        rewardTimes = rewardTimes(logical([1; diff(rewardTimes)>minRewardInteveral])); % remove reward times occuring within minRewardInteveral seconds of eachother
        
        % determine whether reward was surprise (delivered earlier than expected)
        rewardDistances = diff(interp1(data.wheelTimes, data.wheelPositions, rewardTimes));  % distance (m) between rewards
        isRewardSurprise = [false; rewardDistances < nanmedian(rewardDistances) * .75];
        
        % find omissions (cue without reward)
        if exist('Cue', 'var')
            cueInds = find(diff(Cue.values>2)==1)+1;
            cueTimes  = Cue.times(cueInds);
            cueTimes = cueTimes(logical([1; diff(cueTimes)>minRewardInteveral])); % remove reward times occuring within minRewardInteveral seconds of eachother
            [~, closestRewardDistance] = knnsearch(rewardTimes, cueTimes);
            omissionTimes = cueTimes(closestRewardDistance>1);
        else
            omissionTimes = [];
        end

        saveVars('rewardTimes', rewardTimes, 'isRewardSurprise', isRewardSurprise, 'omissionTimes', omissionTimes);
    end
    
    
    
    
    % led indices (run camera)
    if analyzeVar('ledInds') && isRunVid && isWiskVid
        
        if s.verbose; fprintf('%s: getting LED indices (for run camera)... ', session); end
        if ~exist('locations', 'var'); locations = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        
        % settings
        confidenceThresh = .8;
        minDuration = 4;  % (frames)
        minDistanceTime = 250;  % (frames)
        maxDistanceSpace = 4;  % (pixels)
        
        if any(contains(locations.Properties.VariableNames, 'LED'))
            ledInds = getLedInds(locations.LED_2, locations.LED, locations.LED_1, ...
                confidenceThresh, minDuration, minDistanceTime, maxDistanceSpace);
        else
            fprintf(' LED not tracked in trackedfeaturesRaw.csv\n')
            ledInds = [];
        end
        
        saveVars('ledInds', ledInds)
    end
    
    
    
    
    % led indices (wisk camera)
    if analyzeVar('ledIndsWisk') && exist(fullfile(sessionDir, 'trackedFeaturesRaw_wisk.csv'), 'file')
        
        if s.verbose; fprintf('%s: getting LED indices (for whisker camera)... ', session); end
        if ~exist('locationsWisk', 'var'); locationsWisk = readtable(fullfile(sessionDir, 'trackedFeaturesRaw_wisk.csv')); end
        
        if ~isempty(data.ledInds)
        
            % settings
            confidenceThresh = .7;
            minDuration = 4;  % (frames)
            minDistanceTime = 250;  % (frames)
            maxDistanceSpace = 2;  % (pixels)

            ledIndsWisk = getLedInds(locationsWisk.LED_2, locationsWisk.LED, locationsWisk.LED_1, ...
                confidenceThresh, minDuration, minDistanceTime, maxDistanceSpace);
        else
            ledIndsWisk = [];
            fprintf('no LED detected in run camera, so assuming no LED in whisker camera.\n');
        end
        
        saveVars('ledIndsWisk', ledIndsWisk)
    end
    
    
    
    
    % decode obstacle position (based on obstacle track rotary encoder)
    if analyzeVar('obsPositions', 'obsTimes')
        
        if s.verbose; fprintf('%s: decoding obstacle position\n', session); end
        load(fullfile(sessionDir, 'run.mat'), 'obEncodA', 'obEncodB')
        
        if ~isempty(obEncodA.times) && ~isempty(obEncodB.times)
            [obsPositions, obsTimes] = rotaryDecoder(obEncodA.times, obEncodA.level,...
                                                     obEncodB.times, obEncodB.level,...
                                                     s.obEncoderSteps, s.obsPulleyRad, s.targetFs, session);
        else
            obsPositions = [];
            obsTimes = [];
        end
        
        saveVars('obsPositions', obsPositions, 'obsTimes', obsTimes, 'targetFs', s.targetFs);
    end




    % get obstacle on and off times
    % (ensuring that first event is obs turning ON and last is obs turning OFF)
    if analyzeVar('obsOnTimes', 'obsOffTimes')

        if s.verbose; fprintf('%s: getting obstacle on and off times\n', session); end
        load(fullfile(sessionDir, 'run.mat'), 'obsOn')

        if exist('obsOn', 'var')
            firstOnInd = find(obsOn.level, 1, 'first');
            lastOffInd = find(~obsOn.level, 1, 'last');

            obsOn.level = obsOn.level(firstOnInd:lastOffInd);
            obsOn.times = obsOn.times(firstOnInd:lastOffInd);

            obsOnTimes  =  obsOn.times(logical(obsOn.level));
            obsOffTimes = obsOn.times(logical(~obsOn.level));
        else
            obsOnTimes = [];
            obsOffTimes = [];
        end
        
        % get rid of really close obsOnTimes (debounce)
        validBins = (obsOffTimes-obsOnTimes) > .1; % not too close together
        obsOnTimes = obsOnTimes(validBins);
        obsOffTimes = obsOffTimes(validBins);
        
        saveVars('obsOnTimes', obsOnTimes, 'obsOffTimes', obsOffTimes)
    end
    
    
    
    
    % compute whether each trial had a wheel break or not
    if analyzeVar('isWheelBreak')
        
        if s.verbose; fprintf('%s: finding wheel break trials\n', session); end
        load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mat'), 'breaks')

        isWheelBreak = false(size(data.obsOnTimes));
        for i = 1:length(data.obsOnTimes)
            isWheelBreak(i) = any(breaks.times>data.obsOnTimes(i) & breaks.times<data.obsOffTimes(i));
        end
        
        saveVars('isWheelBreak', isWheelBreak);
    end
    
    
    
    
    % check that obstacle tracked wheel position well
    if analyzeVar('obsTracking') && ~isempty(data.obsOnTimes)
        
        if s.verbose; fprintf('%s: checking obstacle tracking of wheel velocity\n', session); end
        
        % settings
        velTime = .01;  % (s) compute velocity over this time window
        obsOnBuffer = .25;  % (s) don't consider obstacle tracking within obsOnBuffer of obstacle first turning on, because this is the time period where the speed is still ramping up
        velTolerance = .02;  % (m/s) tracking is considered good when obstacle velocity is within velTolerance of wheel velocity
        warningThresh = .1;  % if warningThresh of trial has poor obstacle tracking, a warning message is thrown

        % initializations
        wheelVel = getVelocity(data.wheelPositions, velTime, s.targetFs);
        obsVel = getVelocity(data.obsPositions, velTime, s.targetFs);

        % get trial data
        obsTracking = struct();
        rowInd = 1;
        anyPoorTrackingTrials = false;  % keeps tracking of whether any poorly tracked trials have been detected
        for i = 1:length(data.obsOnTimes)

            wheelBins = data.wheelTimes > (data.obsOnTimes(i)+obsOnBuffer) & ...
                        data.wheelTimes < data.obsOffTimes(i);
            obsBins = data.obsTimes > (data.obsOnTimes(i)+obsOnBuffer) & ...
                      data.obsTimes < data.obsOffTimes(i);
            
            wheelVelTrial = wheelVel(wheelBins);
            obsVelTrial = obsVel(obsBins);
            times = data.wheelTimes(wheelBins);
            if ~isequal(times, data.obsTimes(obsBins))
                obsVelTrial = interp1(data.obsTimes(obsBins), obsVelTrial, times);
            end

            obsTracking(rowInd).wheelVel = wheelVelTrial;
            obsTracking(rowInd).obsVel = obsVelTrial;
            obsTracking(rowInd).times = times;
            obsTracking(rowInd).percentBadTracking = nanmean(abs(wheelVelTrial-obsVelTrial) > velTolerance);
            if obsTracking(rowInd).percentBadTracking > warningThresh
                if ~anyPoorTrackingTrials
                    fprintf('%s: WARNING! poor obstacle tracking in trial(s): %i', session, i);
                    anyPoorTrackingTrials = true;
                else
                    fprintf(' %i', i);
                end
            end
            rowInd = rowInd + 1; 
        end
        if anyPoorTrackingTrials; fprintf('\n'); end
        if s.plotObsTracking; plotObsTracking(session, obsTracking); end
        
        saveVars('obsTracking', obsTracking)
    end
    


    
    % get obstacle light on and off times
    % (ensuring that first event is obs turning ON and last is obs turning OFF)
    if analyzeVar('obsLightOnTimes', 'obsLightOffTimes')

        if s.verbose; fprintf('%s: getting obstacle light on and off times\n', session); end
        load(fullfile(sessionDir, 'run.mat'), 'obsLight')

        % settings
        minObsLightInterval = .1;
        obsLightOnTimes = [];
        obsLightOffTimes = [];
        
        if exist('obsLight', 'var')
            
            if isfield(obsLight, 'values') % if recorded as analog input (sessions prior to 191523)
                if any(obsLight.values>2)
                    obsLightOnInds  = find(diff(obsLight.values>2)==1) + 1;
                    obsLightOffInds = find(diff(obsLight.values>2)==-1) + 1;
                    obsLightOnTimes = obsLight.times(obsLightOnInds);
                    obsLightOffTimes = obsLight.times(obsLightOffInds);
                end
            else
                if ~isempty(obsLight.times)
                    obsLightOnTimes = obsLight.times(logical(obsLight.level));
                    obsLightOffTimes = obsLight.times(~logical(obsLight.level));
                end
            end
            
            if ~isempty(obsLightOnTimes)
                
                % remove times occuring within minObsLightInterval seconds of each other
                obsLightOnTimes  = obsLightOnTimes(logical([1; diff(obsLightOnTimes)>minObsLightInterval]));
                obsLightOffTimes = obsLightOffTimes(logical([1; diff(obsLightOffTimes)>minObsLightInterval]));

                % ensure first time is on and last time is off
                obsLightOnTimes =  obsLightOnTimes(obsLightOnTimes < obsLightOffTimes(end));
                obsLightOffTimes = obsLightOffTimes(obsLightOffTimes > obsLightOnTimes(1));
                
            end
        end

        saveVars('obsLightOnTimes', obsLightOnTimes, 'obsLightOffTimes', obsLightOffTimes)
    end
    
    
    
    
    % determine whether light was on or off for every trial
    if analyzeVar('isLightOn')
        
        if isfield(data, 'obsLightOnTimes')
            
            if s.verbose; fprintf('%s: determing whether each trial is light on or light off\n', session); end
            isLightOn = false(size(data.obsOnTimes));

            if ~isempty(data.obsLightOnTimes) 
                for i = 1:length(data.obsOnTimes)
                    isLightOn(i) = min(abs(data.obsOnTimes(i) - data.obsLightOnTimes)) < 1; % did the light turn on near when the obstacle turned on
                end
            end
        else
            isLightOn = [];
        end

        saveVars('isLightOn', isLightOn)
    end




    % get frame timeStamps
    if analyzeVar('frameTimeStamps') && exist(fullfile(sessionDir, 'run.csv'), 'file')

        if s.verbose; fprintf('%s: getting frame time stamps\n', session); end
        load(fullfile(sessionDir, 'run.mat'), 'exposure')

        camMetadata = dlmread(fullfile(sessionDir, 'run.csv')); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
        frameCounts = camMetadata(:,2);
        timeStampsFlir = timeStampDecoderFLIR(camMetadata(:,3));
        
        if ~exist(fullfile(sessionDir, 'alignmentFrames.csv'), 'file')
            frameTimeStamps = getFrameTimes(exposure.times, timeStampsFlir, frameCounts, session);
        else
            alignmentFrames = readtable(fullfile(sessionDir, 'alignmentFrames.csv'));
            frameTimeStamps = getFrameTimes(exposure.times, timeStampsFlir, frameCounts, session, ...
                [alignmentFrames.ttlNumberRun, alignmentFrames.frameNumberRun]);
        end

        saveVars('frameTimeStamps', frameTimeStamps)
    end
    
    
    
    
    % get wisk frame timeStamps
    if analyzeVar('frameTimeStampsWisk') && exist(fullfile(sessionDir, 'wisk.csv'), 'file')
            
        if s.verbose; fprintf('%s: getting wisk frame time stamps\n', session); end
        load(fullfile(sessionDir, 'run.mat'), 'exposure')

        camMetadataWisk = dlmread(fullfile(sessionDir, 'wisk.csv')); % columns: point grey counter, point grey timestamps (uninterpretted)
        frameCountsWisk = camMetadataWisk(:,1);
        timeStampsFlirWisk = timeStampDecoderFLIR(camMetadataWisk(:,2));
        
        if ~exist(fullfile(sessionDir, 'alignmentFrames.csv'), 'file')
            frameTimeStampsWisk = getFrameTimes(exposure.times, timeStampsFlirWisk, frameCountsWisk, session);
        else
            alignmentFrames = readtable(fullfile(sessionDir, 'alignmentFrames.csv'));
            frameTimeStampsWisk = getFrameTimes(exposure.times, timeStampsFlirWisk, frameCountsWisk, session, ...
                [alignmentFrames.ttlNumberWisk, alignmentFrames.frameNumberWisk]);
        end

        saveVars('frameTimeStampsWisk', frameTimeStampsWisk)
    end
    
    
    
    
    % licking
    if analyzeVar('lickTimes') && exist(fullfile(sessionDir, 'trackedFeaturesRaw_wisk.csv'), 'file')
        
        if s.verbose; fprintf('%s: getting lick times\n', session); end
        if ~exist('locationsWisk', 'var'); locationsWisk = readtable(fullfile(sessionDir, 'trackedFeaturesRaw_wisk.csv')); end
        
        % setttings
        confidenceThresh = .5;
        smoothing = 5;  % (frames) mean smoothing for x, y, and confidence
        minTimeDiff = 20;  % (frames) licks must be at least this many frames apart in time
        maxDistance = 40;  % (pixels) maximum vertical (y) distance of tongue to median tongue position
        minDistance = 2;  % (pixels) tongue must be minDistance pixels away (vertically) from median tongue position for lick to count
        
        % get x, y, and confidence
        y = smooth(locationsWisk.tongue_1, smoothing);
        conf = smooth(locationsWisk.tongue_2, smoothing);
        
        % find distance to median tongue position
        tonguePos = median(y(conf>confidenceThresh));
        dist = y-tonguePos;  % vertical distance to median tongue position
        
        % remove low confidence and tracking outlier frames
        valid = conf>confidenceThresh & dist<maxDistance;
        sig = y;  % y position of tongue
        sig = medfilt1(sig, 3);
        sig(~valid) = nan;
        
        sig = sig-tonguePos;  % express relative to median tongue position
        
        % get lick times
        [~, lickInds] = findpeaks(sig, 'MinPeakDistance', minTimeDiff, 'MinPeakHeight', minDistance);
        lickTimes = data.frameTimeStampsWisk(lickInds);
        
        if s.showLickFig
            figure('name', sprintf('%s: lickTimes', session), 'color', 'white', 'position', [108.00 53.00 1664.00 921.00]);
            rows = 6;
            cols = 3;
            xLims = [-1 4];  % (seconds) time pre and post reward to plot
            yLims = [-20 maxDistance]; % (pixels) relative to median tongue position
            
            for i = 1:(rows*cols)
                subplot(rows, cols, i); hold on
                plot(data.frameTimeStampsWisk - data.rewardTimes(i), sig);
                scatter(lickTimes - data.rewardTimes(i), sig(lickInds), 20, 'filled');
                plot([0 0], get(gca, 'ylim'), 'color', [0 0 .6])
                plot(xLims, [minDistance minDistance], 'color', [.6 .6 .6])
                set(gca, 'xlim', xLims, 'ylim', yLims);
            end
            pause(.1)
            xlabel('time from water (s)')
            ylabel('tongue pos - median tongue pos (pixels)')
        end
        
        saveVars('lickTimes', lickTimes);
    end




    % get webCam timeStamps if webCam data exist
    if analyzeVar('webCamTimeStamps') && ...
            exist(fullfile(sessionDir, 'webCam.csv'), 'file') &&...
            exist(fullfile(sessionDir, 'run.csv'), 'file') &&...
            any(strcmp(fieldnames(data), 'frameTimeStamps'))

        if s.verbose; fprintf('%s: getting webcam time stamps\n', session); end

        % get main camera frame times wrt computer clock
        camMetadataRun = dlmread(fullfile(sessionDir, 'run.csv'));
        camSysClock = camMetadataRun(:,1) / 1000;
        camTimeSteps = cumsum([0; diff(camSysClock)<0]);
        camSysClock = camSysClock + camTimeSteps;  % remove discontinuities
        camSysClock = camSysClock - camSysClock(1); % set first time to zero
        
        % get web camera frame times wrt computer clock
        webCamSysClock = dlmread(fullfile(sessionDir, 'webCam.csv')) / 1000; % convert from ms to s
        webCamTimeSteps = cumsum([0; diff(webCamSysClock)<0]);
        webCamSysClock = webCamSysClock + webCamTimeSteps;  % remove discontinuities
        webCamSysClock = webCamSysClock - webCamSysClock(1);  % set first time to zero
        
        % get main camera frame times wrt to spike clock
        camSpikeClock = data.frameTimeStamps;
        
        % determine web camera frame times with respect to spike clock
        bins = ~isnan(camSpikeClock);
        mapping = polyfit(camSysClock(bins), camSpikeClock(bins), 1);
        webCamTimeStamps = polyval(mapping, webCamSysClock);

        saveVars('webCamTimeStamps', webCamTimeStamps)
    end
    
    
    
    
    % get nose position
    if analyzeVar('nosePos') && exist(fullfile(sessionDir, 'trackedFeaturesRaw.csv'), 'file')

        if s.verbose; fprintf('%s: getting nose position\n', session); end

        if ~exist('locations', 'var'); locations = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        noseBotX = median(locations.nose_bot);
        noseBotY = median(locations.nose_bot_1);

        saveVars('nosePos', [noseBotX noseBotY])
    end
    
    
    
    
    % get mToPixMapping and obstacle pixel positions in bottom view
    if analyzeVar('obsPixPositions', 'obsPixPositionsUninterped', 'obsToObsPixPosMappings', ...
            'wheelToObsPixPosMappings', 'pixelsPerM', 'obsPositionsFixed') && ...
            exist(fullfile(sessionDir, 'trackedFeaturesRaw.csv'), 'file') && ...
            ~isempty(data.obsOnTimes)
   
        if s.verbose; fprintf('%s: tracking obstacles in bottom view\n', session); end
        if ~exist('locations', 'var'); locations = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        
        % settings
        confidenceThresh = .8;
        
        % get obs pix positions
        obsHighX = locations.obsHigh_bot;
        obsLowX = locations.obsLow_bot;
        obsHighScores = locations.obsHigh_bot_2;
        obsLowScores = locations.obsLow_bot_2;
        
        % ensure that confidences are high for top and bottom points of obstacle, and ensure that both have x values that are close to one another
        validInds = abs(obsHighX - obsLowX) < 10 & ...
                    obsHighScores>confidenceThresh & ...
                    obsLowScores>confidenceThresh;
        obsPixPositions = mean([obsHighX, obsLowX], 2);
        obsPixPositions(~validInds) = nan;
        
        % store version of obsPixPositions without interpolations=
        obsPixPositionsUninterped = obsPixPositions; % these are the obsPixPositions wihtout any interpolation (only where the obs is in the frame and successfully tracked)
        
        % get positions of obstacle and wheel at the time of each frame
        obsPositionsInterp = interp1(data.obsTimes, data.obsPositions, data.frameTimeStamps); 
        wheelPositionsInterp = interp1(data.wheelTimes, data.wheelPositions, data.frameTimeStamps);
        
        % determine mappings from wheel and obs positions to obsPixPositions
        obsToObsPixPosMappings = nan(length(data.obsOnTimes), 2); % mapping between obs positions (rotary encoder) and obs pix positions (video tracking)
        wheelToObsPixPosMappings = nan(length(data.obsOnTimes), 2); % mapping between wheel positions (rotary encoder) and obs pix positions (video tracking)
        for i = 1:length(data.obsOnTimes)
            bins = data.frameTimeStamps>data.obsOnTimes(i) & ...
                   data.frameTimeStamps<=data.obsOffTimes(i) & ...
                   ~isnan(obsPixPositions);
            
            if any(bins)
                obsToObsPixPosMappings(i,:) = polyfit(obsPositionsInterp(bins), obsPixPositions(bins), 1);
                wheelToObsPixPosMappings(i,:) = polyfit(wheelPositionsInterp(bins), obsPixPositions(bins), 1);
            end          
        end
        
        % use obs position from rotary encoder to infer pix positions when obs is out of frame
        epochTimes = [data.obsOnTimes; data.frameTimeStamps(end)];
        for i = 1:(length(epochTimes)-1)
            epochBins = data.frameTimeStamps>epochTimes(i) & ...
                        data.frameTimeStamps<=epochTimes(i+1);
            interpBins = epochBins & isnan(obsPixPositions); % bins that should be interpolated over
            if any(epochBins)
                obsPixPositions(interpBins) = polyval(obsToObsPixPosMappings(i,:), obsPositionsInterp(interpBins));
            end
        end
        
        % compute obsPositionsFixed, which is realigned st 0 corresponds to the position at which obsetacle is beneath the nose
        obsPositionsFixed = fixObsPositions(data.obsPositions, data.obsTimes, obsPixPositions, ...
            data.frameTimeStamps, data.obsOnTimes, data.obsOffTimes, data.nosePos(1));
        
        % compute pixel per mm ratio
        pixelsPerM = abs(nanmedian(obsToObsPixPosMappings(:,1)));

        saveVars('obsPixPositions', obsPixPositions', ...  % positions of obstacle along track in pixel coordinates // positions when out of frame are inferred based on obstacle rotary encoder and linear mapping between pixel and encoder coordinate systems
                 'obsPixPositionsUninterped', obsPixPositionsUninterped, ...  % positions of obstacle in pixel coordinates only when it is in the camera view and accurately tracked
                 'obsToObsPixPosMappings', obsToObsPixPosMappings, ...  % linear mapping from obstacle position (meters) to obstacle position (pixels)
                 'wheelToObsPixPosMappings', wheelToObsPixPosMappings, ...  % linear mapping from wheel position (meters) to obstacle position (pixels)
                 'pixelsPerM', pixelsPerM, ...
                 'obsPositionsFixed', obsPositionsFixed);
    end
    
    
    
    
    % get wheel points
    if analyzeVar('wheelCenter', 'wheelRadius') && exist(fullfile(sessionDir, 'run.csv'), 'file')
        
        if s.verbose; fprintf('%s: getting wheel center and radius ', session); end
        
        if ~exist('locations', 'var'); locations = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        
        % if wheel was tracked with the neural network
        if any(contains(locations.Properties.VariableNames, 'wheel'))
            fprintf('(using neural network tracking)\n')
            wheelPoints = [median(locations.wheel_left), median(locations.wheel_left_1)
                           median(locations.wheel_mid), median(locations.wheel_mid_1)
                           median(locations.wheel_right), median(locations.wheel_right_1)];
        
        % otherwise use hacky algorithm
        else
            fprintf('(using hacky algorithm)\n')
            wheelPoints = getWheelPoints(session);
        end
        
        [wheelRadius, wheelCenter] = fitCircle(wheelPoints);
        
        saveVars('wheelCenter', wheelCenter, 'wheelRadius', wheelRadius);
    end
    
    
    
    
    % get body angle
    if analyzeVar('bodyAngles') && exist(fullfile(sessionDir, 'trackedFeaturesRaw.csv'), 'file')
        
        if s.verbose; fprintf('%s: getting body angle\n', session); end
        
        if ~exist('locations', 'var'); locations = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        
        % settings
        confidenceThresh = .8;
        percentileLims = [1 99];
        
        tailXY = [locations.tailBase_bot, locations.tailBase_bot_1];
        conf = locations.tailBase_bot_2;
        tailXY = data.nosePos - tailXY;  % set X to number of pixels behind nose
        angles = rad2deg(atan2(tailXY(:,2), tailXY(:,1)));
        
        % get rid of low confidence and outlier time points
        angles(conf<confidenceThresh) = nan;
        lims = prctile(angles, percentileLims);
        angles(angles<lims(1) | angles>lims(2)) = nan;
        
        % interpolate
        angles = fillmissing(angles, 'linear');
        
        saveVars('bodyAngles', angles)
    end
    
    
    
    
    % get height of obs for each trial
    if analyzeVar('obsHeights') && ~isempty(data.obsOnTimes) && isRunVid
        
        if s.verbose; fprintf('%s: getting obstacle heights\n', session); end
        
        % load tracking data if not already open
        if ~exist('locations', 'var'); locations = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        obsTopY = locations.obs_top_1;
        obsTopScores = locations.obs_top_2;
        
        wheelTopPix = data.wheelCenter(2) - data.wheelRadius;
        bins = ~isnan(data.obsPixPositions)' & ...
                      obsTopScores>.99 & ...
                      obsTopY<wheelTopPix;  % obstacle can't be below the top of the wheel
        
        obsHeights = nan(1,length(data.obsOnTimes));
        for i = 1:length(data.obsOnTimes)
            trialBins = data.frameTimeStamps>data.obsOnTimes(i) & ...
                        data.frameTimeStamps<data.obsOffTimes(i) & ...
                        bins;
            medianObsY = median(obsTopY(trialBins));
            obsHeightPix = wheelTopPix - medianObsY;
            obsHeight = obsHeightPix / (data.pixelsPerM/1000) + (s.obsDiameter/2);  % second term accounts for the fact that center of obs is tracked, but height is the topmost part of the obstacle
            obsHeights(i) = obsHeight;
        end
        
        saveVars('obsHeights', obsHeights)
    end
    
    
    
    
    % neural network classifier to determine whether paw is touching obs
    if (analyzeVar('touches', 'touchesPerPaw', 'touchConfidences', 'touchClassNames') || ...
            ~exist(fullfile(sessionDir, 'pawAnalyzed.csv'), 'file') || ...
            isempty(readtable(fullfile(sessionDir, 'pawAnalyzed.csv')))) && ...
            exist(fullfile(sessionDir, 'trackedFeaturesRaw.csv'), 'file') && ...
            ~isempty(data.obsOnTimes)
        
        if s.verbose; fprintf('%s: getting paw contacts\n', session); end
        
        % settings
        pythonPath = 'C:\Users\rick\Anaconda3\envs\fastai\python.exe';
        confidenceThresh = .5;
        confidenceThreshForeDorsal = .6;  % fore dorsal is prone to false positives // emperically .9 results in good sensitivity/specificity tradeoff
        proximityThresh = 20;  % how close does a paw have to be to the obstacle to be assigned to it for a touch
        classesToAssignToPaw = {'fore_dorsal', 'fore_ventral', 'hind_dorsal', 'hind_ventral_low'};  % other touch types will be ignored

        % run neural network classifier
        if ~exist(fullfile(sessionDir, 'pawAnalyzed.csv'), 'file') || ...
                s.rerunPawContactNetwork || ...
                isempty(readtable(fullfile(sessionDir, 'pawAnalyzed.csv')))
            
            fprintf('  %s: running paw contact neural network\n', session);
            save(fullfile(sessionDir, 'runAnalyzed.mat'), '-struct', 'data');  % first save the file so analyzeVideo.py can access it
            [~,~] = system([pythonPath ' ' fullfile('tracking', 'pawContact', 'expandanalyze.py') ...
                            ' ' fullfile(getenv('OBSDATADIR'), 'sessions') ' ' session]);
        end
        
        % get touch class for each frame
        pawAnalyzed = readtable(fullfile(sessionDir, 'pawAnalyzed.csv'));
        touchClassNames = pawAnalyzed.Properties.VariableNames(2:end); % first column is frame number
        allScores = table2array(pawAnalyzed(:,2:end));
        [touchConfidences, classInds] = max(allScores, [], 2);
        noTouchInd = find(strcmp(touchClassNames, 'no_touch'));
        
        % remove low confidence touches
        foreDorsalInd = find(strcmp(touchClassNames, 'fore_dorsal'));
        classInds(touchConfidences<confidenceThresh & classInds~=foreDorsalInd) = noTouchInd;
        classInds(touchConfidences<confidenceThreshForeDorsal & classInds==foreDorsalInd) = noTouchInd;
        
        touches = ones(size(data.frameTimeStamps)) * noTouchInd;
        touches(pawAnalyzed.framenum) = classInds; % not all frames are analyzed // only those whethere paws are close to obstacle
        
        
        % figure out which paws are touching obs in each touch frame
        
        % get xz positions for paws
        if ~exist('locations', 'var'); locations = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
        [locations, features] = fixTracking(locations, data.frameTimeStamps, data.pixelsPerM, 'scoreThresh', scoreThresh);
        pawXZ = nan(size(locations,1), 2, 4);
        for i = 1:4
            pawXBin = contains(features, ['paw' num2str(i)]) & contains(features, '_bot');
            pawZBins = contains(features, ['paw' num2str(i)]) & contains(features, '_top');
            pawXZ(:,1,i) = locations(:,1,pawXBin);
            pawXZ(:,2,i) = locations(:,2,pawZBins);
        end
        
        % get obs height in pixels for each trial
        % (replace missing values with median of tracked values for each trial)
        obsBin = ismember(features, 'obs_top');
        obsHeights = nan(size(locations,1),1);
        for i = 1:length(data.obsOnTimes)
            trialBins = data.frameTimeStamps>data.obsOnTimes(i) & ...
                        data.frameTimeStamps<data.obsOffTimes(i);
            medianHgt = nanmedian(locations(trialBins,2,obsBin));
            obsHeights(trialBins) = locations(trialBins,2,obsBin);
            obsHeights(trialBins & isnan(locations(:,2,obsBin))) = medianHgt;
        end
        obsHeights = medfilt1(obsHeights,5);

        % get xz distance of paws to obs at all times
        dx = squeeze(pawXZ(:,1,:)) - repmat(data.obsPixPositions',1,4);
        dz = squeeze(pawXZ(:,2,:)) - repmat(obsHeights,1,4);
        pawDistances = sqrt(dx.^2 + dz.^2);

        % determine which paws are touching obs
        classesToAssignToPawInds = find(ismember(touchClassNames, classesToAssignToPaw));
        touchesToAssignToPaws = touches .* ismember(touches, classesToAssignToPawInds);
        touchesPerPaw = repmat(touchesToAssignToPaws,1,4) .* double(pawDistances<proximityThresh);
        
        % only fore paw classes for forepaws and hind paw classes for hind paws
        foreClassInds = find(contains(touchClassNames, 'fore'));
        hindClassInds = find(contains(touchClassNames, 'hind'));
        touchesPerPaw(:,[2,3]) = touchesPerPaw(:,[2,3]) .* double(ismember(touchesPerPaw(:,[2,3]), foreClassInds));
        touchesPerPaw(:,[1,4]) = touchesPerPaw(:,[1,4]) .* double(ismember(touchesPerPaw(:,[1,4]), hindClassInds));
        
        % touchConfidences are only for analyzed frames // asign confidence of 1 to unanalyzed frames
        temp = touchConfidences;
        touchConfidences = ones(1,length(data.frameTimeStamps));
        touchConfidences(pawAnalyzed.framenum) = temp;
        
        saveVars('touches', touches, ...  % categorical vector of the types of touches in each frame, with categories listed in touchClassNames (noTouchInd for no touch)
                 'touchesPerPaw', touchesPerPaw, ...  % categorical matrix encoding the type of touch for each paw at every frame (zeros for no touch)
                 'touchClassNames', touchClassNames, ...  % names of touch classes
                 'touchConfidences', touchConfidences)  % condifence in touch classification
    end
    
    
    
    
    % run whisker contact network
    if analyzeVar('wiskContactFrames', 'wiskContactFramesConfidences', 'wiskContactPositions', 'wiskContactTimes') && ...
            ~isempty(data.obsOnTimes) && ...
            exist(fullfile(sessionDir, 'runWisk.mp4'), 'file') && ...
            exist(fullfile(sessionDir, 'trackedFeaturesRaw.csv'), 'file')
        
        if s.verbose; fprintf('%s: getting whisker contacts\n', session); end
        
        % settings
        pythonPath = 'C:\Users\rick\Anaconda3\envs\deepLabCut\python.exe';
        
        % run neural network classifier
        if ~exist(fullfile(sessionDir, 'whiskerAnalyzed.csv'), 'file') ||  s.rerunWiskContactNetwork
            fprintf('  %s: running wisk contact network\n', session)
            if anythingAnalyzed; save(fullfile(sessionDir, 'runAnalyzed.mat'), '-struct', 'data'); end % first save the file so analyzeVideo.py can access it
            [~,~] = system([pythonPath ' tracking\whiskerContact\cropanalyzevideo.py ' getenv('OBSDATADIR') 'sessions ' session]);
        end
        wiskContactData = readtable(fullfile(sessionDir, 'whiskerAnalyzed.csv'));
        
        % extract contact positions and times
        contactTimes = nan(1,length(data.obsOnTimes));
        contactPositions = nan(1,length(data.obsOnTimes));
        
        % set to nan -1 contact frames
        notFoundBins = wiskContactData.framenum==-1;
        wiskContactData(notFoundBins,:) = num2cell(nan(sum(notFoundBins),2));
        
        for i = 1:length(data.obsOnTimes)
            if wiskContactData.framenum(i)>0
                time = data.frameTimeStampsWisk(wiskContactData.framenum(i));
                if time>data.obsOnTimes(i) && time<data.obsOffTimes(i)
                    contactTimes(i) = time;
                    contactPositions(i) = interp1(data.obsTimes, data.obsPositionsFixed, contactTimes(i));
                else
                    fprintf('%s: WARNING - wisk contact detected when obstacle not on!\n', session)
                end
            end
        end
        
        % todo: check for invalid contact positions here?
        saveVars('wiskContactFrames', wiskContactData.framenum, ...
                 'wiskContactFramesConfidences', wiskContactData.confidence, ...
                 'wiskContactTimes', contactTimes, ...
                 'wiskContactPositions', contactPositions);
    end
    
    
    
    
    % save results
    if anythingAnalyzed
        save(fullfile(sessionDir, 'runAnalyzed.mat'), '-struct', 'data')
        fprintf('%s: data analyzed and saved\n', session)
    else
        fprintf('%s: no new variables computed\n', session)
    end
    
    
    
    
    % ---------
    % FUNCTIONS
    % ---------
    
    % determines whether to analyze a variable based on whether it has
    % already been analyzed or if overwrite is requested
    function analyze = analyzeVar(varargin)
        % varargin: vars that will be checked for whether they should be computed
        
        if ~strcmp(s.overwriteVars, 'all')
            analyze = any(~ismember(varargin, computedVars)) || any(ismember(s.overwriteVars, varargin));
        else
            analyze = true;
        end
    end
    
    
    % adds variables to 'data' struct
    function saveVars(varargin)
        % varargin: name-value pairs where the name is the name of the field to be added to 'data', and value is the value to be assigned to that field
        
        for v = 1:2:length(varargin); data.(varargin{v}) = varargin{v+1}; end
        anythingAnalyzed = true;
    end
    
    
    % get led indices for run or whisker cam
    function inds = getLedInds(conf, x, y, confidenceThresh, minDuration, minDistanceTime, maxDistanceSpace)
        
        conf = conf > confidenceThresh;
        inds = find(diff(conf)==1)+1;
        
        if ~isempty(inds)
            offInds = find(diff(conf)==-1)+1;

            offInds = offInds(offInds>inds(1));  % ensure first switch is ON
            inds = inds(inds<offInds(end));  % ensure last switch is OFF

            % remove brief inds
            isBrief = offInds - inds >= minDuration;
            inds = inds(isBrief);
            
            % remove inds too close together in time
            isTooClose = [Inf; diff(inds)] < minDistanceTime;
            inds = inds(~isTooClose);
            
            % remove inds too far from the median LED location (the LED should never move)
            distances = sqrt(sum(([x(inds), y(inds)] - [median(x(inds)), median(y(inds))]).^2,2));  % distances of each detected LED to the median LED position
            inds = inds(distances < maxDistanceSpace);
            
            if length(data.rewardTimes) ~= length(inds)
                fprintf('WARNING! %i rewards detected in Spike2 and %i in camera! Saving empty vector.\n', length(data.rewardTimes), length(inds));
            else
                fprintf('\n')
            end
        else
            inds = [];
            fprintf('no LED detected\n');
        end
    end
end





