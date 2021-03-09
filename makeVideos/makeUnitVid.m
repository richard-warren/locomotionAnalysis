function makeUnitVid(session, unit_id, fileName, varargin)

% creates videos of mouse running with the raw neural trace plotted in real time
% along with the sound of the neuron // tell it the session, the unit_id
% of the neuron to plot, and a matrix with the start and stop times
% (columns) of trials to plot (rows) // automatically selects trialsToShow
% trials from timeEpochs matrix // if timeEpochs not provided, assumes
% timeEpochs correspond to obs on and off times


% todo: gate sound of spikes



% settings
s.trialsToShow = 5; % how many random trials to show if specific trial numbers are not indicated
s.vidType = 'showObsEvents'; % choose from 'showObsEvents' or 'showRewardEvents' or 'showSpecificTimeWindows' // todo: should select automatically based on arguments
s.specificObsTrials = []; % pick specific trials to show. Time is from obsOn to obsOff. Trial number refers to obsOn and obsOff times.
s.specificRewardTrials = []; % pick specific trials to show. Time is reward delivery +- timeBuffer. Trial number refers to reward times.
s.specificTimeWindows = []; % put in specific start times and end times for generating videos.
s.timeBuffer = [2 2]; % how many seconds before and after reward delivery to show. Default is 2s before and 2s after.

s.playbackSpeed = 0.25;
s.voltageWindow = .75;
s.audioGain = 15;
s.yLims = [-800 800];
s.lowPassFreq = 6000; % 6000 // set to false to turn off lowpass
s.maskVoltage = true;  % whether to mask out the voltage when the unit is not firing

s.includeWiskCam = true;
s.compressVideo = true;

s.spkScatterColor = [1 1 0];
s.lineColors = [.8 .4 1];

% reassign settings contained in opts
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end
% fileName = strrep(fileName, ' ', '_');  % a hack to avoid errors with spaces when using ffmpeg


% set up video readers / writer
disp('initializing...')
% vidName = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run_originalDimensions.mp4');
vidName = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4');
if ~exist(vidName, 'file'); concatTopBotVids(session); end
vid = VideoReader(vidName);

if s.includeWiskCam; vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4')); end
initialFs = vid.FrameRate;
vidWriter = vision.VideoFileWriter(fileName, ...
    'AudioInputPort', true, ...
    'FrameRate', round(initialFs*s.playbackSpeed));
% vidWriter.VideoCompressor = 'MJPEG Compressor';

% load spike data
disp('loading spike data...');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'obsOnTimes', 'obsOffTimes', 'rewardTimes', 'wiskContactFrames', 'isLightOn', ...
    'obsPixPositions', 'frameTimeStampsWisk', 'wiskContactTimes')

% get position where wisk frame should overlap with run frame
if s.includeWiskCam    
    disp('figuring out overlapping b/w runWisk and run vid...');
    obsInWiskCamInds = find(obsPixPositions>vid.Width-50 & obsPixPositions<vid.Width);    
    % find first time point at which both wisk and run cams have a frame and obs is in wisk cam
    for i = obsInWiskCamInds
        wiskInd = find(frameTimeStampsWisk==frameTimeStamps(i), 1);
        if ~isempty(wiskInd); topInd = i; break; end
    end
    
    [wiskFrame, yWiskPos, xWiskPos, wiskScaling] = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, topInd);
end


% determine frame dimensions
if s.includeWiskCam
    frameDim = size(wiskFrame);
else
    frameDim = [vid.Height, vid.Width];
end


% get neural data
disp('getting neural data...');
ephysInfo = getSessionEphysInfo(session);
[~, unit_ids] = getGoodSpkInds(session);
bestChannels = getBestChannels(session, 'returnPhysicolLayout', true);  % should this be true or false?
bestChannel = bestChannels(unit_id==unit_ids);
getVoltage = @(data, channel, inds) data.Data.Data(channel,inds);
data = memmapfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysInfo.ephysFolder, [ephysInfo.fileNameBase '_CHs.dat']), ...
    'Format', {'int16', [ephysInfo.channelNum, ephysInfo.smps], 'Data'}, 'Writable', false);

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'openEphysToSpikeMapping', 'spkRates', 'timeStamps', 'unit_ids', 'spkTimes')
timeStampsMapped = polyval(openEphysToSpikeMapping, ephysInfo.timeStamps);
audioSmpsPerFrame = round((1/initialFs) * ephysInfo.fs);


% create low pass filter
disp('creating low pass filter...');
if s.lowPassFreq
    lp = s.lowPassFreq * 2 / ephysInfo.fs;
    ls = s.lowPassFreq * 4 / ephysInfo.fs; % one octave above pass band
    ls = min(ls, .999); % make sure stopband freq doesn't exceed nyquist freqency
    [NN, Wn] = buttord(lp, ls, 3, 60);
    [B1,A1] = butter(NN, Wn, 'low');
end


% set up figure
disp('setting up fig for vid...');
fig = figure('color', [0 0 0], 'position', [1400, 400, frameDim(2), frameDim(1)*(1/.8)], 'menubar', 'none');
traceLength = s.voltageWindow*ephysInfo.fs;

axes('position', [0 .2 1 .8], 'CLim', [0 255]); colormap gray
im = image(uint8(zeros(frameDim)), 'CDataMapping', 'scaled'); hold on;
set(gca, 'Visible', 'off')

plotAxis = axes('position', [0 0 1 .2], 'Color', 'black');
tracePlot = plot(plotAxis, 1:traceLength, nan(1,traceLength), 'color', 'white'); set(gca, 'color', 'black'); hold on
set(gca, 'Visible', 'off', 'YLimMode', 'manual', 'YLim', s.yLims)

unitInd = find(unit_ids==unit_id);
obsOnLine = line(plotAxis, [0 0], s.yLims, 'linewidth', 2, 'color', s.lineColors);
obsOnText = text(plotAxis, 0, s.yLims(2), 'obstacle on', 'Color', 'white');
wiskLine = line(plotAxis, [0 0], s.yLims, 'linewidth', 2, 'color', s.lineColors);
wiskText = text(plotAxis, 0, s.yLims(2), 'whisker contact', 'Color', 'white');
rewardLine = line(plotAxis, [0 0], s.yLims, 'linewidth', 2, 'color', s.lineColors);
rewardText = text(plotAxis, 0, s.yLims(2), 'reward', 'Color', 'white');


% get timeEpochs for trials to display in the vid
% current setup only supports 'showObsEvents' and 'showRewardEvents'
minTime = timeStamps(find(~isnan(spkRates(unitInd,:)),1,'first'));
maxTime = timeStamps(find(~isnan(spkRates(unitInd,:)),1,'last'));

switch s.vidType
    case 'showObsEvents'
        timeEpochs = cat(2, obsOnTimes, obsOffTimes);
        validTrials = find(timeEpochs(:,1)>minTime & timeEpochs(:,2)<maxTime);  % make sure trials aren't too long
        
        if isempty(s.specificObsTrials) && isempty(s.specificRewardTrials)  % specific trials are not indicated 
            trialsToShow = validTrials(round(linspace(2,length(validTrials)-1,s.trialsToShow)));
            timeEpochs = timeEpochs(trialsToShow, :);
        else
            trialNum = s.specificObsTrials;
            trialsToShow = trialNum(ismember(trialNum, validTrials));
            timeEpochs = timeEpochs(trialsToShow, :);
        end
        
    case 'showRewardEvents'
        timeEpochs = cat(2, rewardTimes - s.timeBuffer(1), rewardTimes + s.timeBuffer(2));
        validTrials = find(timeEpochs(:,1)>minTime & timeEpochs(:,2)<maxTime); 

        % make sure trials aren't too long
        if isempty(s.specificObsTrials) && isempty(s.specificRewardTrials)  % specific trials are not indicated 
            trialsToShow = validTrials(round(linspace(1,length(validTrials),s.trialsToShow)));
            timeEpochs = timeEpochs(trialsToShow, :);
        else
            trialNum = s.specificRewardTrials;
            trialsToShow = trialNum(ismember(trialNum, validTrials));
            timeEpochs = timeEpochs(trialsToShow, :);
        end
        
    case 'showSpecificTimeWindows'
        bins = s.specificTimeWindows(:,1)>minTime & s.specificTimeWindows(:,2)<maxTime;  % bins for epochs that fall within the valid times for the unit
        timeEpochs = s.specificTimeWindows(bins,:); 
        if any(~bins)
            fprintf('WARNING! Trials ')
            fprintf('%i ', find(~bins))
            fprintf('ommitted because they are out of the valid time range for the unit!\n')
        end
        trialsToShow = 1:size(timeEpochs, 1);        
end


% scatter only those points appearning within time epochs
spkBins = any(spkTimes{unitInd} > (timeEpochs(:,1)-s.voltageWindow)' & ...
              spkTimes{unitInd} < (timeEpochs(:,2)+s.voltageWindow)', 2);
scatter(spkTimes{unitInd}(spkBins), ...  % todo: should really do this on a trial by trial basis... this is super inefficeint
    repmat(s.yLims(1),1,sum(spkBins)) + range(s.yLims)*.1, ... % y values
    10, s.spkScatterColor, 'filled');


% create video
fprintf('\nwriting video... ')
for i = 1:length(trialsToShow)
    
    % get trial frame indices
    trialInds = find(frameTimeStamps>timeEpochs(i,1) & ...
                     frameTimeStamps<timeEpochs(i,2));
                 
    % get voltage for entire trial
    voltageBins = timeStampsMapped>timeEpochs(i,1)-s.voltageWindow-1 & ...
                  timeStampsMapped<timeEpochs(i,2)+1; % add and subtract 1 as a buffer
    timeStampsSub = timeStampsMapped(voltageBins);
    voltageRaw = getVoltage(data, bestChannel, voltageBins); % maybe replace this by reading data within trial prior to writing each trials data...
    
    if s.lowPassFreq
        voltageRaw = filter(B1, A1, voltageRaw); % filter forwards then backwards to avoid phase shifts
        voltageRaw = int16(fliplr(filter(B1, A1, fliplr(voltageRaw))));
    end
    
    
    if s.maskVoltage
        % create mask to silence recording when unit isn't firing
    
        % first determine height of kernel for single spike (a hack // should be done in getFiringRate())
        temp = getFiringRate(mean(timeStampsSub), 'times', timeStampsSub(:)', ...
            'kernel', 'doubleExp', 'kernelRise', 1e-6, 'kernelFall', 1e-3);
        maskMax = max(temp);

        % now compute and normalize mask
        mask = getFiringRate(spkTimes{unitInd}, 'times', timeStampsSub(:)', ...
            'kernel', 'doubleExp', 'kernelRise', 1e-6, 'kernelFall', 1e-3);
        mask(mask>maskMax) = maskMax;
        mask = mask / maskMax;

        % apply the mask! omg!
        voltageRaw = int16(double(voltageRaw) .* mask);
    end
    
    voltage = double(voltageRaw) * ephysInfo.bitVolts;  % convert to microvolts
    
    % update obstacle and whisker contact lines
    if isLightOn(trialsToShow(i)); obsOnString = 'obstacle (light on)'; else; obsOnString = 'obstacle (light off)'; end
    updateTextAndLine(obsOnText, obsOnLine, obsOnTimes, obsOnString)
    updateTextAndLine(wiskText, wiskLine, wiskContactTimes)
    updateTextAndLine(rewardText, rewardLine, rewardTimes)
    
    % get frames for trials
    for j = trialInds'
        
        % get run frame with wisk frame matched to it
        if s.includeWiskCam
            [frame, ~, ~, ~] = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, j, ...
                'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, 'isPaddingWhite', false);
        else
            frame = rgb2gray(read(vid, j));
        end
        
        % add trial number onto frames
        position = [10, 10];
        if strcmp(s.vidType, 'showObsEvents')
            textString = ['trial ', num2str(trialsToShow(i))];
            RGB = insertText(frame, position, textString, 'TextColor','white');
            frame = rgb2gray(RGB);
        elseif strcmp(s.vidType, 'showRewardEvents')
            textString = ['Reward Trial ', num2str(trialsToShow(i))];
            RGB = insertText(frame, position, textString, 'TextColor','white');
            frame = rgb2gray(RGB);
        end
        
        % update frame
        set(im, 'CData', frame);
        
        % get voltage
        traceStartInd = find(timeStampsSub>(frameTimeStamps(j)-s.voltageWindow), 1, 'first');
        traceInds = traceStartInd:traceStartInd+traceLength-1;
        trace = voltage(traceInds);
        times = timeStampsSub(traceInds)';
        set(tracePlot, 'xdata', times, 'ydata', trace);
        set(gca, 'xlim', [times(1) times(end)])
        
        % get audio
        audioStartInd = find(timeStampsSub>=frameTimeStamps(j), 1, 'first');
        audio = voltageRaw(audioStartInd:audioStartInd+audioSmpsPerFrame-1)';
        audio = audio * s.audioGain;
        
        % add fade in/out for first and last sample
        if i==1 && j==trialInds(1)
            audio = int16(double(audio) .* linspace(0,1,length(audio))'); % fade in on first sample
        elseif i==length(trialsToShow) && j==trialInds(end)
            audio = int16(double(audio) .* linspace(1,0,length(audio))'); % fade out on last sample
        end
        
        % write to file
        frame = getframe(fig);
        vidWriter(frame.cdata, audio);

    end
end
release(vidWriter);
close(fig)


% compress video
if s.compressVideo
    baseDir = fileparts(fileName);
    [~,~] = system([fileName(1) ': & ffmpeg -i "' fileName '" -vcodec mpeg4 -vb 10M -y ' baseDir '\temp.avi']); % run ffmpeg to compress file
    delete(fileName)
    movefile([baseDir '\temp.avi'], fileName)
end

disp(' all done!')





% update position of lines and text marking events of interest
function updateTextAndLine(text, line, eventTimes, textString)
    eventInd = find(eventTimes>=timeEpochs(i,1)-s.voltageWindow & ...
                    eventTimes<=timeEpochs(i,2),1,'first');
    if ~isempty(eventInd)
        eventTime = eventTimes(eventInd);
        set(text, 'Position', [eventTime+s.voltageWindow*.01 s.yLims(2)])
        set(line, 'XData', [eventTime eventTime])
    end
    if exist('textString', 'var'); set(text, 'String', textString); end
end

end