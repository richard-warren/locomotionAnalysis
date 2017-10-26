function spikeAnalysis(dataDir, varsToOverWrite)

    % performs preliminary analyses on spike data and saves results in runAnalyzed.mat
    %
    % computes the following:
    %         reward times
    %         debounced touch signal and touch on/off times
    %         decoding stepper motor commands
    %         decoding obstacle position
    %         decoding wheel position
    %         getting obstacle on and off times
    %         getting frame time stamps
    %         getting webcam time stamps
    %
    % for each session, loads existing runAnalyzed.mat
    % if a computed variable is not already stored in runAnalyzed.mat AND the files necessary to compute it exist, it computes the variable
    % if a variable name is included in cell array varsToOverWrite, it re-computes the variable even if it already exists, so long as
    % the necessary files exist in the directory to compute it
    %
    % input:     dataDir           directory containing session folders
    %            varsToOverwrite   cell array of of variables that should be re-computed



    % settings
    targetFs = 1000; % frequency that positional data will be resampled to
    minRewardInteveral = 1;

    % rig characteristics
    whEncoderSteps = 2880; % 720cpr * 4
    wheelRad = 95.25; % mm
    obEncoderSteps = 1000; % 250cpr * 4
    obsRad = 96 / (2*pi); % radius of timing pulley driving belt of obstacles platform

    % if no variables to overwrite are specified, set to default
    if nargin==1
        varsToOverWrite = {' '};
    end
    
    % initialize folder list
    dataFolders = uigetdir2(dataDir, 'select folders to analyze');


    
    
    % iterate over data folders and analyze those that have not been analyzed
    for i = 1:length(dataFolders)
        
        anythingAnalyzed = false;

        % load or initialize data structure
        sessionDir = [dataFolders{i} '\'];
        nameStartInd = find(dataFolders{1}=='\',1,'last') + 1;
        
        if exist([sessionDir 'runAnalyzed.mat'], 'file')
            varStruct = load([sessionDir 'runAnalyzed.mat']);
        else
            varStruct = struct();
        end
        varNames = fieldnames(varStruct);
        


        % analyze reward times
        if analyzeVar('rewardTimes', varNames, varsToOverWrite)
            
            fprintf('%s: getting reward times\n', dataFolders{i}(nameStartInd:end))
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
                fprintf('%s: decoding stepper motor commands\n', dataFolders{i}(nameStartInd:end))
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
                fprintf('%s: decoding obstacle position\n', dataFolders{i}(nameStartInd:end))
                [obsPositions, obsTimes] = rotaryDecoder(obEncodA.times, obEncodA.level,...
                                                             obEncodB.times, obEncodB.level,...
                                                             obEncoderSteps, obsRad, targetFs);
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
            
            fprintf('%s: decoding wheel position\n', dataFolders{i}(nameStartInd:end))
            load([sessionDir 'run.mat'], 'whEncodA', 'whEncodB')
            
            [wheelPositions, wheelTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                         whEncodB.times, whEncodB.level,...
                                                         whEncoderSteps, wheelRad, targetFs);
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
       
            fprintf('%s: getting obstacle on and off times\n', dataFolders{i}(nameStartInd:end))
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
        
        
        
        
        % debounce touch signal and get touch on/off times
        if analyzeVar('touchSig', varNames, varsToOverWrite) ||...
           analyzeVar('touchOnTimes', varNames, varsToOverWrite) ||...
           analyzeVar('touchOffTimes', varNames, varsToOverWrite)
            
            fprintf('%s: debouncing touch signal\n', dataFolders{i}(nameStartInd:end))
            load([sessionDir 'run.mat'], 'touch', 'breaks')
                        
            % decode stepper motor
            if exist('breaks', 'var')
                [touchSig, touchOnTimes, touchOffTimes] = debounceTouch(touch.values, touch.times, varStruct.obsOffTimes, breaks.times);
            else
                [touchSig, touchOnTimes, touchOffTimes] = debounceTouch(touch.values, touch.times, varStruct.obsOffTimes);
            end
            
            % save values
            varStruct.touchSig = touchSig;
            varStruct.touchOnTimes = touchOnTimes;
            varStruct.touchOffTimes = touchOffTimes;
            anythingAnalyzed = true;
        end
        
        
        
        
        % get frame timeStamps
        if analyzeVar('frameTimeStamps', varNames, varsToOverWrite)
            
            if exist([sessionDir 'run.csv'], 'file')
                
                fprintf('%s: getting frame time stamps\n', dataFolders{i}(nameStartInd:end))
                load([sessionDir '\run.mat'], 'exposure')

                % get camera metadata and spike timestamps
                camMetadata = dlmread([sessionDir '\run.csv']); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
                frameCounts = camMetadata(:,2);
                timeStampsFlir = timeStampDecoderFLIR(camMetadata(:,3));

                if length(exposure.times) >= length(frameCounts)
                    frameTimeStamps = getFrameTimes(exposure.times, timeStampsFlir, frameCounts);
                else
                    disp('  there are more frames than exposure TTLs... saving frameTimeStamps as empty vector')
                    frameTimeStamps = [];
                end
                
                % save values
                varStruct.frameTimeStamps = frameTimeStamps;
                anythingAnalyzed = true;
            end
        end
        
        
        
        
        % get webCam timeStamps if webCam data exist
        if analyzeVar('webCamTimeStamps', varNames, varsToOverWrite)
            
            if exist([sessionDir 'webCam.csv'], 'file') &&...
               exist([sessionDir 'run.csv'], 'file') &&...
               any(strcmp(fieldnames(varStruct), 'frameTimeStamps'))
                
                fprintf('%s: getting webcam time stamps\n', dataFolders{i}(nameStartInd:end))
                
                % load data
                camMetadata = dlmread([sessionDir '\run.csv']);
                camSysClock = camMetadata(:,1) / 1000;
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
        
        
        
        
        % save results
        if anythingAnalyzed
            save([sessionDir 'runAnalyzed.mat'], '-struct', 'varStruct')
            fprintf('----------\n')
        end
    end
    
    
    
    
    % ---------
    % FUNCTIONS
    % ---------
    
    function analyze = analyzeVar(var, varNames, varsToOverWrite)
        analyze = ~any(strcmp(varNames, var)) || any(strcmp(varsToOverWrite, var));
    end
end

