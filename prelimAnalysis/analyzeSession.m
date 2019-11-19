function analyzeSession(session, varargin)

    % performs preliminary analyses on a session, including turning raw
    % data in run.mat into useful data in runAnalyzed.mat // also has
    % wrappers that execute neural network algorithms // variable names are
    % listed as arguments to the analyzeVar() function in the code blocks
    % below // to recompute a variable, pass 'overwriteVars' as follows:
    %
    % analyzeSession('191118_001', 'overwriteVars', 'reawrdTimes')
    % analyzeSession('191118_001', 'overwriteVars', {'rewardTimes, isLightOn'}) % for multiple variables
    


    % settings
    s.targetFs = 1000; % frequency that positional data will be resampled to
    s.overwriteVars = '';
    s.plotObsTracking = true;  % whether to check obstacle tracking of wheel velocity by plotting them on top of one another
    
    s.rerunRunNetwork = false;
    s.rerunFaceNetwork = false;
    s.rerunWiskContactNetwork = false;
    s.rerunPawContactNetwork = false;

    % rig characteristics
    s.whEncoderSteps = 2880; % 720cpr * 4
    s.wheelRad = 95.25; % mm
    s.obEncoderSteps = 1000; % 250cpr * 4
    s.obsPulleyRad = 96 / (2*pi); % radius of timing pulley driving belt of obstacles platform
    s.obsDiameter = 3.175; % (mm)

    % initializations
    if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
    if ischar(s.overwriteVars); s.overwriteVars = {s.overwriteVars}; end  % make sure overwriteVars is in cell format
    anythingAnalyzed = false; % results are only saved if something was found that wasn't already analyzed

    % load or initialize data structure
    sessionDir = fullfile(getenv('OBSDATADIR'), 'sessions', session);
    if exist(fullfile(sessionDir, 'runAnalyzed.mat'), 'file') && ~isequal(s.overwriteVars, {'all'})
        data = load(fullfile(sessionDir, 'runAnalyzed.mat'));
    else
        data = struct();
    end
    computedVars = fieldnames(data);


    
    
    % analyze reward times
    if analyzeVar('rewardTimes')
        
        fprintf('%s: getting reward times\n', session)
        load(fullfile(sessionDir, 'run.mat'), 'reward')
        
        % settings
        minRewardInteveral = 1;  % remove rewards detected within minRewardInterval of eachother
        
        % find reward times
        if isfield(reward, 'values')  % if recorded as analog input (sessions prior to 191523)
            rewardInds = find(diff(reward.values>2)==1) + 1;
            rewardTimes = reward.times(rewardInds);
        else
            rewardTimes = reward.times(logical(reward.level));  % keep only transitions from low to high
        end

        % remove reward times occuring within minRewardInteveral seconds of eachother
        rewardTimes = rewardTimes(logical([1; diff(rewardTimes)>minRewardInteveral]));

        saveVars('rewardTimes', rewardTimes);
    end


    

    % decode obstacle position (based on obstacle track rotary encoder)
    if analyzeVar('obsPositions', 'obsTimes')
        load(fullfile(sessionDir, 'run.mat'), 'obEncodA', 'obEncodB')
        fprintf('%s: decoding obstacle position\n', session)
        
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




    % decode wheel position
    if analyzeVar('wheelPositions', 'wheelTimes')

        fprintf('%s: decoding wheel position\n', session)
        load(fullfile(sessionDir, 'run.mat'), 'whEncodA', 'whEncodB')

        [wheelPositions, wheelTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                     whEncodB.times, whEncodB.level,...
                                                     s.whEncoderSteps, s.wheelRad, s.targetFs, session);
        
         saveVars('wheelPositions', wheelPositions, 'wheelTimes', wheelTimes, 'targetFs', s.targetFs);
    end




    % get obstacle on and off times
    % (ensuring that first event is obs turning ON and last is obs turning OFF)
    if analyzeVar('obsOnTimes', 'obsOffTimes')

        fprintf('%s: getting obstacle on and off times\n', session)
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
    
    
    
    
    % check that obstacle tracked wheel position well
    if analyzeVar('obsTracking') && ~isempty(data.obsOnTimes)
        
        fprintf('%s: checking obstacle tracking of wheel velocity\n', session)
        
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
                    fprintf('  %s: WARNING! poor obstacle tracking in trial(s): %i', session, i);
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

        fprintf('%s: getting obstacle light on and off times\n', session)
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
            
            fprintf('%s: determing whether each trial is light on or light off\n', session)
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

        fprintf('%s: getting frame time stamps\n', session)
        load(fullfile(sessionDir, 'run.mat'), 'exposure')

        camMetadata = dlmread(fullfile(sessionDir, 'run.csv')); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
        frameCounts = camMetadata(:,2);
        timeStampsFlir = timeStampDecoderFLIR(camMetadata(:,3));
        frameTimeStamps = getFrameTimes(exposure.times, timeStampsFlir, frameCounts, session);

        saveVars('frameTimeStamps', frameTimeStamps)
    end
    
    
    
    
    % get wisk frame timeStamps
    if analyzeVar('frameTimeStampsWisk') && exist(fullfile(sessionDir, 'wisk.csv'), 'file')
            
        fprintf('%s: getting wisk frame time stamps\n', session)
        load(fullfile(sessionDir, 'run.mat'), 'exposure')

        camMetadataWisk = dlmread(fullfile(sessionDir, 'wisk.csv')); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
        frameCountsWisk = camMetadataWisk(:,1);
        timeStampsFlirWisk = timeStampDecoderFLIR(camMetadataWisk(:,2));
        frameTimeStampsWisk = getFrameTimes(exposure.times, timeStampsFlirWisk, frameCountsWisk, session);

        saveVars('frameTimeStampsWisk', frameTimeStampsWisk)
    end




    % get webCam timeStamps if webCam data exist
    if analyzeVar('webCamTimeStamps') && ...
            exist(fullfile(sessionDir, 'webCam.csv'), 'file') &&...
            exist(fullfile(sessionDir, 'run.csv'), 'file') &&...
            any(strcmp(fieldnames(data), 'frameTimeStamps'))

        fprintf('%s: getting webcam time stamps\n', session)

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

        fprintf('%s: getting nose position\n', session)

        locationsTable = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); % get raw tracking data
        noseBotX = median(locationsTable.nose_bot);
        noseBotY = median(locationsTable.nose_bot_1);

        saveVars('nosePos', [noseBotX noseBotY])
    end
    
    
    
    
    % get mToPixMapping and obstacle pixel positions in bottom view
    if analyzeVar('obsPixPositions', 'obsPixPositionsUninterped', 'obsToObsPixPosMappings', ...
            'wheelToObsPixPosMappings', 'pixelsPerM', 'obsPositionsFixed') && ...
            exist(fullfile(sessionDir, 'trackedFeaturesRaw.csv'), 'file') && ...
            ~isempty(data.obsOnTimes)
   
        fprintf('%s: tracking obstacles in bottom view\n', session)
        
        % load tracking data if not already open
        if ~exist('locationsTable', 'var'); locationsTable = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        
        % get obs pix positions
        obsHighX = locationsTable.obsHigh_bot;
        obsLowX = locationsTable.obsLow_bot;
        obsHighScores = locationsTable.obsHigh_bot_2;
        obsLowScores = locationsTable.obsLow_bot_2;
        
        % ensure that confidences are high for top and bottom points of obstacle, and ensure that both have x values that are close to one another
        validInds = abs(obsHighX - obsLowX) < 10 & ...
                    obsHighScores>0.99 & ...
                    obsLowScores>0.99;
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
        
        fprintf('%s: getting wheel center and radius\n', session)
        
        wheelPoints = getWheelPoints(session);
        [wheelRadius, wheelCenter] = fitCircle(wheelPoints);
        
        saveVars('wheelCenter', wheelCenter, 'wheelRadius', wheelRadius);
    end
    
    
    
    
    % get body angle
    if analyzeVar('bodyAngles') && exist(fullfile(sessionDir, 'trackedFeaturesRaw.csv'), 'file')
        
        fprintf('%s: getting body angle\n', session)
        
        if ~exist('locationsTable', 'var'); locationsTable = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        bodyAngles = getSessionBodyAngles(locationsTable, data.nosePos);
        
        saveVars('bodyAngles', bodyAngles)
    end
    
    
    
    
    % get height of obs for each trial
    if analyzeVar('obsHeights') && ~isempty(data.obsOnTimes)
        
        fprintf('%s: getting obstacle heights\n', session)
        
        % load tracking data if not already open
        if ~exist('locationsTable', 'var'); locationsTable = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        obsTopY = locationsTable.obs_top_1;
        obsTopScores = locationsTable.obs_top_2;
        
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
            obsHeight = obsHeightPix / (data.pixelsPerM/1000) + (s.obsDiameter/2); % second term accounts for the fact that center of obs is tracked, but height is the topmost part of the obstacle
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
        
        fprintf('%s: getting paw contacts\n', session)
        
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
            
            fprintf('%s: running paw contact neural network\n', session);
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
        if ~exist('locationsTable', 'var'); locationsTable = readtable(fullfile(sessionDir, 'trackedFeaturesRaw.csv')); end
        [locations, features] = fixTracking(locationsTable, data.frameTimeStamps);
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
        
        fprintf('%s: getting whisker contacts\n', session)
        
        % settings
        pythonPath = 'C:\Users\rick\Anaconda3\envs\deepLabCut\python.exe';
        
        % run neural network classifier
        if ~exist(fullfile(sessionDir, 'whiskerAnalyzed.csv'), 'file') ||  s.rerunWiskContactNetwork
            fprintf('%s: running wisk contact network\n', session)
            if anythingAnalyzed; save(fullfile(sessionDir, 'runAnalyzed.mat'), '-struct', 'data'); end % first save the file so analyzeVideo.py can access it
            [~,~] = system([pythonPath ' tracking\whiskerContact\cropanalyzevideo.py ' getenv('OBSDATADIR') 'sessions ' session]);
        end
        wiskContactData = readtable(fullfile(sessionDir, 'whiskerAnalyzed.csv'));
        
        % extract contact positions and times
        contactTimes = nan(1,length(data.obsOnTimes));
        contactPositions = nan(1,length(data.obsOnTimes));
        
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
end