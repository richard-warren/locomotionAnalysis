function spikeAnalysis(dataDir, dataFolder, varsToOverWrite)

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
    %         video tracking of obstacle in bottom view
    %         for each track, gets linear mapping of obsPosition (m) to obsPixPosition (pixels), which can be used to subsequently determine mapping between pixels and meters
    %         get times of wisk contact, and pixels at which wisk contact occurs within wisk vid
    %
    % for each session, loads existing runAnalyzed.mat
    % if a computed variable is not already stored in runAnalyzed.mat AND the files necessary to compute it exist, it computes the variable
    % if a variable name is included in cell array varsToOverWrite, it re-computes the variable even if it already exists, so long as
    % the necessary files exist in the directory to compute it
    %
    % input:     dataDir           directory containing session folders
    %            dataFolder        !!!
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
    
    


    anythingAnalyzed = false;

    % load or initialize data structure
    sessionDir = [dataDir dataFolder '\'];

    if exist([sessionDir 'runAnalyzed.mat'], 'file')
        varStruct = load([sessionDir 'runAnalyzed.mat']);
    else
        varStruct = struct();
    end
    varNames = fieldnames(varStruct);



    % analyze reward times
    if analyzeVar('rewardTimes', varNames, varsToOverWrite)

        fprintf('%s: getting reward times\n', dataFolder)
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
    if analyzeVar('motorPositions', varNames, varsToOverWrite) ||...
       analyzeVar('motorTimes', varNames, varsToOverWrite)

        load([sessionDir 'run.mat'], 'step', 'stepDir')

        % decode stepper motor
        if exist('stepDir', 'var') && ~isempty(stepDir.times)
            fprintf('%s: decoding stepper motor commands\n', dataFolder)
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
    if analyzeVar('obsPositions', varNames, varsToOverWrite) ||...
       analyzeVar('obsTimes', varNames, varsToOverWrite)

        load([sessionDir 'run.mat'], 'obEncodA', 'obEncodB')

        if exist('obEncodA', 'var') && ~isempty(obEncodA.times)
            fprintf('%s: decoding obstacle position\n', dataFolder)
            [obsPositions, obsTimes] = rotaryDecoder(obEncodA.times, obEncodA.level,...
                                                         obEncodB.times, obEncodB.level,...
                                                         obEncoderSteps, obsRad, targetFs, dataFolder);
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
    if analyzeVar('wheelPositions', varNames, varsToOverWrite) ||...
       analyzeVar('wheelTimes', varNames, varsToOverWrite)

        fprintf('%s: decoding wheel position\n', dataFolder)
        load([sessionDir 'run.mat'], 'whEncodA', 'whEncodB')

        [wheelPositions, wheelTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                     whEncodB.times, whEncodB.level,...
                                                     whEncoderSteps, wheelRad, targetFs, dataFolder);
        % save values
        varStruct.wheelPositions = wheelPositions;
        varStruct.wheelTimes = wheelTimes;
        varStruct.targetFs = targetFs;
        anythingAnalyzed = true;
    end




    % get obstacle on and off times
    % (ensuring that first event is obs turning ON and last is obs turning OFF)
    if analyzeVar('obsOnTimes', varNames, varsToOverWrite) ||...
       analyzeVar('obsOffTimes', varNames, varsToOverWrite)

        fprintf('%s: getting obstacle on and off times\n', dataFolder)
        load([sessionDir 'run.mat'], 'obsOn')

        if exist('obsOn', 'var')
            firstOnInd  = find(obsOn.level, 1, 'first');
            lastOffInd  = find(~obsOn.level, 1, 'last');

            obsOn.level = obsOn.level(firstOnInd:lastOffInd);
            obsOn.times = obsOn.times(firstOnInd:lastOffInd);

            obsOnTimes  =  obsOn.times(logical(obsOn.level));
            obsOffTimes = obsOn.times(logical(~obsOn.level));
        else
            obsOnTimes = [];
            obsOffTimes = [];
        end

        % save values
        varStruct.obsOnTimes = obsOnTimes;
        varStruct.obsOffTimes = obsOffTimes; 
        anythingAnalyzed = true;
    end




    % get obstacle light on and off times
    % (ensuring that first event is obs turning ON and last is obs turning OFF)
    if analyzeVar('obsLightOnTimes', varNames, varsToOverWrite) ||...
       analyzeVar('obsLightOffTimes', varNames, varsToOverWrite)

        fprintf('%s: getting obstacle light on and off times\n', dataFolder)
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
    if analyzeVar('isLightOn', varNames, varsToOverWrite)
        if isfield(varStruct, 'obsLightOnTimes') && isfield(varStruct, 'obsOnTimes')
            
            fprintf('%s: determing whether each trial is light on or light off\n', dataFolder)
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
    if exist('touch', 'var')
        if analyzeVar('touchSig', varNames, varsToOverWrite) ||...
           analyzeVar('touchOnTimes', varNames, varsToOverWrite) ||...
           analyzeVar('touchSigTimes', varNames, varsToOverWrite) ||...
           analyzeVar('touchOffTimes', varNames, varsToOverWrite)

            fprintf('%s: debouncing touch signal\n', dataFolder)
            load([sessionDir 'run.mat'], 'touch', 'breaks')

            % decode stepper motor
            if exist('breaks', 'var')
                [touchSig, touchOnTimes, touchOffTimes] = debounceTouch(touch.values, touch.times, varStruct.obsOffTimes, breaks.times);
            else
                [touchSig, touchOnTimes, touchOffTimes] = debounceTouch(touch.values, touch.times, varStruct.obsOffTimes);
            end

            % save values
            varStruct.touchSig = touchSig;
            varStruct.touchSigTimes = touch.times;
            varStruct.touchOnTimes = touchOnTimes;
            varStruct.touchOffTimes = touchOffTimes;
            anythingAnalyzed = true;
        end
    end




    % get frame timeStamps
    if analyzeVar('frameTimeStamps', varNames, varsToOverWrite)

        if exist([sessionDir 'run.csv'], 'file')

            fprintf('%s: getting frame time stamps\n', dataFolder)
            load([sessionDir '\run.mat'], 'exposure')

            % get camera metadata and spike timestamps
            camMetadataWisk = dlmread([sessionDir '\run.csv']); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
            frameCountsWisk = camMetadataWisk(:,2);
            timeStampsFlirWisk = timeStampDecoderFLIR(camMetadataWisk(:,3));

            if length(exposure.times) >= length(frameCountsWisk)
                frameTimeStamps = getFrameTimes(exposure.times, timeStampsFlirWisk, frameCountsWisk, dataFolder);
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
    if analyzeVar('frameTimeStampsWisk', varNames, varsToOverWrite)

        if exist([sessionDir 'wisk.csv'], 'file')
            
            fprintf('%s: getting wisk frame time stamps\n', dataFolder)
            load([sessionDir '\run.mat'], 'exposure')

            % get camera metadata and spike timestamps
            camMetadataWisk = dlmread([sessionDir '\wisk.csv']); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
            frameCountsWisk = camMetadataWisk(:,1);
            timeStampsFlirWisk = timeStampDecoderFLIR(camMetadataWisk(:,2));

            if length(exposure.times) >= length(frameCountsWisk)
                frameTimeStampsWisk = getFrameTimes(exposure.times, timeStampsFlirWisk, frameCountsWisk, dataFolder);
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
    if analyzeVar('webCamTimeStamps', varNames, varsToOverWrite)

        if exist([sessionDir 'webCam.csv'], 'file') &&...
           exist([sessionDir 'run.csv'], 'file') &&...
           any(strcmp(fieldnames(varStruct), 'frameTimeStamps'))

            fprintf('%s: getting webcam time stamps\n', dataFolder)

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
    
    
    
    
    % get position of tip of nose
    if analyzeVar('nosePos', varNames, varsToOverWrite)
        
        if exist([dataDir dataFolder '\runBot.mp4'], 'file')

            fprintf('%s: getting nose position in bottom view\n', dataFolder)

            % get nose position
            vidBot = VideoReader([dataDir dataFolder '\runBot.mp4']);
            [noseX, noseY] = getNosePos(vidBot);

            % save
            varStruct.nosePos = [noseX noseY];
            anythingAnalyzed = true;
        end
    end




    % track the pixel positions of the obstacle in bot
    if analyzeVar('obsPixPositions', varNames, varsToOverWrite) || ...
       analyzeVar('mToPixMapping', varNames, varsToOverWrite)

        if ~isempty('obsOntimes') && ...
           any(strcmp(fieldnames(varStruct), 'frameTimeStamps'))

            fprintf('%s: tracking obstacles in bottom view\n', dataFolder)
            
            % settings
            vidBot = VideoReader([dataDir dataFolder '\runBot.mp4']);
            xLims = [60 vidBot.Width-40];
            yLims = [1 vidBot.Height];
            pixThreshFactor = 2.25;
            invertColors = false;
            showTracking = true;
            obsMinThickness = 10;

            % track obstacle in bottom view
            [obsPixPositions, mappings] = trackObstacles(vidBot, varStruct.obsOnTimes, varStruct.obsOffTimes,...
                varStruct.frameTimeStamps, varStruct.obsPositions, varStruct.obsTimes,...
                xLims, yLims, pixThreshFactor, obsMinThickness, invertColors, showTracking);

            % save
            varStruct.obsPixPositions = obsPixPositions;
            varStruct.mToPixMapping = mappings;
            anythingAnalyzed = true;
        end
    end
    


    
    % get wisk contact times and pixels at which contacts occur
    if analyzeVar('wiskTouchPixels', varNames, varsToOverWrite) ||...
       analyzeVar('wiskTouchSignal', varNames, varsToOverWrite)
        
        if ~isempty('obsOntimes') && ...
           exist([sessionDir 'runWisk.mp4'], 'file') && ...
           isfield(varStruct, 'frameTimeStampsWisk')  %!!! probably needs more checks here
            
            fprintf('%s: analyzing wisk contacts\n', dataFolder)

            % settings
            showTracking = false;

            % initializations
            vidWisk = VideoReader([dataDir dataFolder '\runWisk.mp4']);

            % track obstacle in bottom view
            [wiskTouchSignal, wiskTouchPixels] = getWiskContacts(vidWisk, showTracking,...
                varStruct.frameTimeStampsWisk, varStruct.obsPixPositions);

            % save
            varStruct.wiskTouchSignal = wiskTouchSignal;
            varStruct.wiskTouchPixels = wiskTouchPixels;
            anythingAnalyzed = true;
        end
    end




    % save results
    if anythingAnalyzed
        save([sessionDir 'runAnalyzed.mat'], '-struct', 'varStruct')
        fprintf('%s: data analyzed and saved\n', dataFolder)
    end
    
    
    
    
    % ---------
    % FUNCTIONS
    % ---------
    
    function analyze = analyzeVar(var, varNames, varsToOverWrite)
        analyze = ~any(strcmp(varNames, var)) || any(strcmp(varsToOverWrite, var));
    end
end


