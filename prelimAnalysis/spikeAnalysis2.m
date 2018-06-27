function spikeAnalysis2(session, varsToOverWrite)

    % performs preliminary analyses on spike data and saves results in runAnalyzed.mat
    %
    % computes the following:
    %         reward times
    %         debounced touch signal and touch on/off times
    %         decoding stepper motor commands
    %         decoding obstacle position
    %         decoding wheel position
    %         getting obstacle on and off times
    %         getting obstacle light on and off times
    %         whether a trial has light on or not
    %         getting frame time stamps
    %         getting wisk frame time stamps
    %         getting webcam time stamps
    %         get position of tip of nose from bot view
    %         linear mapping between meters and pixels
    %
    % for each session, loads existing runAnalyzed.mat
    % if a computed variable is not already stored in runAnalyzed.mat AND the files necessary to compute it exist, it computes the variable
    % if a variable name is included in cell array varsToOverWrite, it re-computes the variable even if it already exists, so long as
    % the necessary files exist in the directory to compute it
    %
    % input:     sessions          session to analyze
    %            varsToOverwrite   cell array of of variables that should be re-computed



    % settings
    targetFs = 1000; % frequency that positional data will be resampled to
    minRewardInteveral = 1;
    minObsLightInterval = .1;

    % rig characteristics
    whEncoderSteps = 2880; % 720cpr * 4
    wheelRad = 95.25; % mm
    obEncoderSteps = 1000; % 250cpr * 4
    obsRad = 96 / (2*pi); % radius of timing pulley driving belt of obstacles platform

    % if no variables to overwrite are specified, set to default
    if ~exist('varsToOverWrite', 'var')
        varsToOverWrite = {' '};
    end
    
    


    anythingAnalyzed = false; % results are only saved if something was found that wasn't already analyzed

    % load or initialize data structure
    sessionDir = [getenv('OBSDATADIR') 'sessions\' session '\'];

    if exist([sessionDir 'runAnalyzed.mat'], 'file')
        varStruct = load([sessionDir 'runAnalyzed.mat']);
    else
        varStruct = struct();
    end
    varNames = fieldnames(varStruct);


    
    
    % analyze reward times
    if analyzeVar({'rewardTimes'}, varNames, varsToOverWrite)

        fprintf('%s: getting reward times\n', session)
        load([sessionDir 'run.mat'], 'reward')

        % find reward times
        rewardInds = find(diff(reward.values>2)==1) + 1;
        rewardTimes = reward.times(rewardInds);

        % remove reward times occuring within minRewardInteveral seconds of eachother
        rewardTimes = rewardTimes(logical([1; diff(rewardTimes)>minRewardInteveral]));

        % save values
        varStruct.rewardTimes = rewardTimes;
        anythingAnalyzed = true;
    end




    % decode stepper motor commands
    if analyzeVar({'motorPositions', 'motorTimes'}, varNames, varsToOverWrite)

        load([sessionDir 'run.mat'], 'step', 'stepDir')

        % decode stepper motor
        if exist('stepDir', 'var') && ~isempty(stepDir.times)
            fprintf('%s: decoding stepper motor commands\n', session)
            [motorPositions, motorTimes] = motorDecoder(stepDir.level, stepDir.times, step.times, targetFs);
        else
            motorPositions = [];
            motorTimes = [];
        end

        % save values
        varStruct.motorPositions = motorPositions;
        varStruct.motorTimes = motorTimes;
        varStruct.targetFs = targetFs;
        anythingAnalyzed = true;
    end




    % decode obstacle position (based on obstacle track rotary encoder)
    if analyzeVar({'obsPositions', 'obsTimes'}, varNames, varsToOverWrite)

        load([sessionDir 'run.mat'], 'obEncodA', 'obEncodB')

        if exist('obEncodA', 'var') && ~isempty(obEncodA.times)
            fprintf('%s: decoding obstacle position\n', session)
            [obsPositions, obsTimes] = rotaryDecoder(obEncodA.times, obEncodA.level,...
                                                     obEncodB.times, obEncodB.level,...
                                                     obEncoderSteps, obsRad, targetFs, session);
        else
            obsPositions = [];
            obsTimes = [];
        end

        % save values
        varStruct.obsPositions = obsPositions;
        varStruct.obsTimes = obsTimes;
        varStruct.targetFs = targetFs;
        anythingAnalyzed = true;
    end




    % decode wheel position
    if analyzeVar({'wheelPositions', 'wheelTimes'}, varNames, varsToOverWrite)

        fprintf('%s: decoding wheel position\n', session)
        load([sessionDir 'run.mat'], 'whEncodA', 'whEncodB')

        [wheelPositions, wheelTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                     whEncodB.times, whEncodB.level,...
                                                     whEncoderSteps, wheelRad, targetFs, session);
        % save values
        varStruct.wheelPositions = wheelPositions;
        varStruct.wheelTimes = wheelTimes;
        varStruct.targetFs = targetFs;
        anythingAnalyzed = true;
    end




    % get obstacle on and off times
    % (ensuring that first event is obs turning ON and last is obs turning OFF)
    if analyzeVar({'obsOnTimes', 'obsOffTimes'}, varNames, varsToOverWrite)

        fprintf('%s: getting obstacle on and off times\n', session)
        load([sessionDir 'run.mat'], 'obsOn')

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
        

        % save values
        varStruct.obsOnTimes = obsOnTimes;
        varStruct.obsOffTimes = obsOffTimes; 
        anythingAnalyzed = true;
    end
    
    
    
    
    % determine obstacle height for every trial
    if analyzeVar({'obsHeights'}, varNames, varsToOverWrite) && ...
       exist([sessionDir '\trialInfo.csv'], 'file')
       
        fprintf('%s: getting obstacle heights\n', session)
        obsHeightData = dlmread([sessionDir '\trialInfo.csv']); % columns: obsHeight (mm), timestamp (s), timestamps (ms)
        obsTimesSeconds = obsHeightData(:,2);
        validInds = [abs(diff(obsTimesSeconds))>0; true]; %if two times are very close to eachother only believe the second of the two times
        obsHeights = obsHeightData(validInds,1);
        try
        obsHeights = obsHeights(1:length(varStruct.obsOnTimes)); % this gets rid of last entry, which we shouldn't need
        catch; keyboard; end
        
        % there should be one more obsHeight than there all trials because an obs height is randomly chosen
        if sum(validInds)-1 > length(obsHeights)
            fprintf('%s: WARNING: problem assigning obstacle height metadata to correct trial\n', session)
        end

        % save values
        varStruct.obsHeights = obsHeights;
        anythingAnalyzed = true;
    end




    % get obstacle light on and off times
    % (ensuring that first event is obs turning ON and last is obs turning OFF)
    if analyzeVar({'obsLightOnTimes', 'obsLightOffTimes'}, varNames, varsToOverWrite)

        fprintf('%s: getting obstacle light on and off times\n', session)
        load([sessionDir 'run.mat'], 'obsLight')

        if exist('obsLight', 'var') && any(obsLight.values>2)

            % find reward times
            obsLightOnInds  = find(diff(obsLight.values>2)==1) + 1;
            obsLightOffInds = find(diff(obsLight.values>2)==-1) + 1;
            obsLightOnTimes = obsLight.times(obsLightOnInds);
            obsLightOffTimes = obsLight.times(obsLightOffInds);

            % remove reward times occuring within minRewardInteveral seconds of each other
            obsLightOnTimes  = obsLightOnTimes(logical([1; diff(obsLightOnTimes)>minObsLightInterval]));
            obsLightOffTimes = obsLightOffTimes(logical([1; diff(obsLightOffTimes)>minObsLightInterval]));

            % ensure first time is on and last time is off
            obsLightOnTimes =  obsLightOnTimes(obsLightOnTimes < obsLightOffTimes(end));
            obsLightOffTimes = obsLightOffTimes(obsLightOffTimes > obsLightOnTimes(1));

        else
            obsLightOnTimes = [];
            obsLightOffTimes = [];
        end

        % save values
        varStruct.obsLightOnTimes = obsLightOnTimes;
        varStruct.obsLightOffTimes = obsLightOffTimes; 
        anythingAnalyzed = true;
    end
    
    
    
    
    % determine whether light was on or off for every trial
    if analyzeVar({'isLightOn'}, varNames, varsToOverWrite)
        
        if isfield(varStruct, 'obsLightOnTimes') && isfield(varStruct, 'obsOnTimes')
            
            fprintf('%s: determing whether each trial is light on or light off\n', session)
            isLightOn = false(size(varStruct.obsOnTimes));

            if ~isempty(varStruct.obsLightOnTimes) 
                for i = 1:length(varStruct.obsOnTimes)
                    isLightOn(i) = min(abs(varStruct.obsOnTimes(i) - varStruct.obsLightOnTimes)) < 1; % did the light turn on near when the obstacle turned on
                end
            end
        else
            isLightOn = [];
        end

        % save values
        varStruct.isLightOn = isLightOn;
        anythingAnalyzed = true;
    end



    % debounce touch signal and get touch on/off times
    if analyzeVar({'touchSig', 'touchSigTimes', 'touchOnTimes', 'touchOffTimes'}, varNames, varsToOverWrite)

            fprintf('%s: getting touch signal\n', session)
            load([sessionDir 'run.mat'], 'touch', 'breaks')

            % debounce touch signal
%             if exist('breaks', 'var')
%                 [touchSig, touchOnTimes, touchOffTimes] = debounceTouch(touch.values, touch.times, varStruct.obsOffTimes, breaks.times);
%             else
%                 [touchSig, touchOnTimes, touchOffTimes] = debounceTouch(touch.values, touch.times, varStruct.obsOffTimes);
%             end
            touchSig = double(touch.values>2);
            touchOnTimes = touch.times([false; diff(touch.values>2)==1]);
            touchOffTimes = touch.times([diff(touch.values>2)==-1; false]);

            % save values
            varStruct.touchSig = touchSig;
            varStruct.touchSigTimes = touch.times;
            varStruct.touchOnTimes = touchOnTimes;
            varStruct.touchOffTimes = touchOffTimes;
            anythingAnalyzed = true;
    end




    % get frame timeStamps
    if analyzeVar({'frameTimeStamps'}, varNames, varsToOverWrite)

        if exist([sessionDir 'run.csv'], 'file')

            fprintf('%s: getting frame time stamps\n', session)
            load([sessionDir '\run.mat'], 'exposure')

            % get camera metadata and spike timestamps
            camMetadataWisk = dlmread([sessionDir '\run.csv']); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
            frameCountsWisk = camMetadataWisk(:,2);
            timeStampsFlirWisk = timeStampDecoderFLIR(camMetadataWisk(:,3));

            if length(exposure.times) >= length(frameCountsWisk)
                frameTimeStamps = getFrameTimes(exposure.times, timeStampsFlirWisk, frameCountsWisk, session);
            else
                disp('  there are more frames than exposure TTLs... saving frameTimeStamps as empty vector')
                frameTimeStamps = [];
            end

            % save values
            varStruct.frameTimeStamps = frameTimeStamps;
            anythingAnalyzed = true;
        end
    end
    
    
    
    
    % get wisk frame timeStamps
    if analyzeVar({'frameTimeStampsWisk'}, varNames, varsToOverWrite)

        if exist([sessionDir 'wisk.csv'], 'file')
            
            fprintf('%s: getting wisk frame time stamps\n', session)
            load([sessionDir '\run.mat'], 'exposure')

            % get camera metadata and spike timestamps
            camMetadataWisk = dlmread([sessionDir '\wisk.csv']); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
            frameCountsWisk = camMetadataWisk(:,1);
            timeStampsFlirWisk = timeStampDecoderFLIR(camMetadataWisk(:,2));

            if length(exposure.times) >= length(frameCountsWisk)
                frameTimeStampsWisk = getFrameTimes(exposure.times, timeStampsFlirWisk, frameCountsWisk, session);
            else
                disp('  there are more frames than exposure TTLs... saving frameTimeStamps as empty vector')
                frameTimeStampsWisk = [];
            end

            % save values
            varStruct.frameTimeStampsWisk = frameTimeStampsWisk;
            anythingAnalyzed = true;
        end
    end




    % get webCam timeStamps if webCam data exist
    if analyzeVar({'webCamTimeStamps'}, varNames, varsToOverWrite)

        if exist([sessionDir 'webCam.csv'], 'file') &&...
           exist([sessionDir 'run.csv'], 'file') &&...
           any(strcmp(fieldnames(varStruct), 'frameTimeStamps'))

            fprintf('%s: getting webcam time stamps\n', session)

            % load data
            camMetadataWisk = dlmread([sessionDir '\run.csv']);
            camSysClock = camMetadataWisk(:,1) / 1000;
            camSpikeClock = varStruct.frameTimeStamps;
            webCamSysClock = dlmread([sessionDir '\webCam.csv']) / 1000; % convert from ms to s

            % remove discontinuities
            webCamTimeSteps = cumsum([0; diff(webCamSysClock)<0]);
            webCamSysClock = webCamSysClock + webCamTimeSteps;
            webCamSysClock = webCamSysClock - webCamSysClock(1); % set first time to zero

            camTimeSteps = cumsum([0; diff(camSysClock)<0]);
            camSysClock = camSysClock + camTimeSteps;
            camSysClock = camSysClock - camSysClock(1); % set first time to zero

            % determine spike clock times from system clock times
            validInds = ~isnan(camSpikeClock);
            sysToSpike = polyfit(camSysClock(validInds), camSpikeClock(validInds), 1);
            webCamSpikeClock = webCamSysClock * sysToSpike(1) + sysToSpike(2);

            % save
            varStruct.webCamTimeStamps = webCamSpikeClock;
            anythingAnalyzed = true;
        end
    end
    
    
    
    
    % get nose position
    if analyzeVar({'nosePos'}, varNames, varsToOverWrite) && ...
       exist([sessionDir 'trackedFeaturesRaw.csv'], 'file')

        fprintf('%s: getting nose position\n', session)

        % load data
        locationsTable = readtable([sessionDir 'trackedFeaturesRaw.csv']); % get raw tracking data
        noseBotX = median(locationsTable.nose_bot);
        noseBotY = median(locationsTable.nose_bot_1);

        % save
        varStruct.nosePos = [noseBotX noseBotY];
        anythingAnalyzed = true;
    end
    
    
    
    % get mToPixMapping and obstacle pixel positions in bottom view
    if analyzeVar({'obsPixPositions', 'mToPixMapping'}, varNames, varsToOverWrite) && ...
       exist([sessionDir 'trackedFeaturesRaw.csv'], 'file') && ...
       ~isempty(varStruct.obsOnTimes)

        fprintf('%s: tracking obstacles in bottom view\n', session)
        
        % load tracking data if not already open
        if ~exist('locationsTable', 'var'); locationsTable = readtable([sessionDir 'trackedFeaturesRaw.csv']); end
        
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
        obsPositionsInterp = interp1(varStruct.obsTimes, varStruct.obsPositions, varStruct.frameTimeStamps); % get position of obstacle for all frames
        
        % use obs position from rotary encoder to infer pix positions when obs is out of frame
        epochTimes = [0; varStruct.obsOnTimes; varStruct.frameTimeStamps(end)];
        trialMappings = nan(length(varStruct.obsOnTimes), 2);
        for i = 1:(length(epochTimes)-1)
            epochBins = varStruct.frameTimeStamps>epochTimes(i) & ...
                        varStruct.frameTimeStamps<=epochTimes(i+1);
            getMappingBins = epochBins & ~isnan(obsPixPositions);
            interpBins =  epochBins & isnan(obsPixPositions);
            if any(epochBins)
                trialMappings(i,:) = polyfit(obsPositionsInterp(getMappingBins), obsPixPositions(getMappingBins), 1);
                obsPixPositions(interpBins) = obsPositionsInterp(interpBins)*varStruct.mToPixMapping(1) +varStruct.mToPixMapping(2);
            end
        end

        % save
        varStruct.obsPixPositions = obsPixPositions';
        varStruct.mToPixMapping = nanmedian(trialMappings,1);
        anythingAnalyzed = true;
    end
    
    
    
    
    % get wheel points
    if analyzeVar({'wheelCenter', 'wheelRadius'}, varNames, varsToOverWrite)
        
        fprintf('%s: getting wheel center and radius\n', session)
        
        vidTop = VideoReader([sessionDir 'runTop.mp4']);
        wheelPoints = getWheelPoints(vidTop);
        [wheelRadius, wheelCenter] = fitCircle(wheelPoints);
        
        % save
        varStruct.wheelCenter = wheelCenter;
        varStruct.wheelRadius = wheelRadius;
        anythingAnalyzed = true;
    end
    
    
    
    
    % get height of obs for each trial
    if analyzeVar({'obsHeightsVid'}, varNames, varsToOverWrite) && ...
       any(strcmp(fieldnames(varStruct), 'obsPixPositions'))
        
        fprintf('%s: getting obstacle heights\n', session)
        
        % settings
        obsDiameter = 3.175; % (mm)
        
        % load tracking data if not already open
        if ~exist('locationsTable', 'var'); locationsTable = readtable([sessionDir 'trackedFeaturesRaw.csv']); end
        obsTopY = locationsTable.obs_top_1;
        obsTopScores = locationsTable.obs_top_2;
        
        obsHeightsVid = nan(1,length(varStruct.obsOnTimes));
        for i = 1:length(varStruct.obsOnTimes)
            trialBins = varStruct.frameTimeStamps>varStruct.obsOnTimes(i) & ...
                        varStruct.frameTimeStamps<varStruct.obsOffTimes(i) & ...
                        ~isnan(varStruct.obsPixPositions)' & ...
                        obsTopScores>.99;
            medianObsY = median(obsTopY(trialBins));
            obsHeightPix = varStruct.wheelCenter(2) - varStruct.wheelRadius - medianObsY;
            obsHeightsVid(i) = (obsHeightPix / abs(varStruct.mToPixMapping(1)))*1000 + (obsDiameter/2); % second term accounts for the fact that center of obs is tracked, but height is the topmost part of the obstacle
        end
        
        % temp
        figure; scatter(varStruct.obsHeights, obsHeightsVid); hold on; line(get(gca,'xlim'), get(gca,'xlim'))
        
        % save
        varStruct.obsHeightsVid = obsHeightsVid;
        anythingAnalyzed = true;
    end
    
    
    
    
    % save results
    if anythingAnalyzed
        save([sessionDir 'runAnalyzed.mat'], '-struct', 'varStruct')
        fprintf('%s: data analyzed and saved\n', session)
    end
    
    
    
    
    % ---------
    % FUNCTIONS
    % ---------
    
    function analyze = analyzeVar(vars, varNames, varsToOverWrite)
        if ~strcmp(varsToOverWrite,'all')
            analyze = ~any(ismember(varNames, vars)) || any(ismember(varsToOverWrite, vars));
        else
            analyze = true;
        end
    end
end


