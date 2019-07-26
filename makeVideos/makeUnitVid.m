function makeUnitVid(session, unit_id, fileName, timeEpochs, opts)

% creates videos of mouse running with the raw neural trace plotted in real time
% along with the sound of the neuron // tell it the session, the unit_id
% of the neuron to plot, and a matrix with the start and stop times
% (columns) of trials to plot (rows) // automatically selects trialsToShow
% trials from timeEpochs matrix // if timeEpochs not provided, assumes
% timeEpochs correspond to obs on and off times


% TO DO:
% gate sound of spikes



% settings
s.trialsToShow = 5;
s.contrastLims = [.1 .9]; % pixels at these proportional values are mapped to 0 and 255
s.playbackSpeed = 0.25;
s.voltageWindow = .75;
s.audioGain = 15;
s.yLims = [-800 800];
s.wiskBorder = 2;
s.lowPassFreq = 6000; % 6000 // set to false to turn off lowpass

s.includeWiskCam = true;
s.compressVideo = true;

s.spkScatterColor = [1 1 0];
s.lineColors = [.8 .4 1];

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end



% set up video readers / writer
disp('initializing...')
vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
if s.includeWiskCam; vidWisk = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runWisk.mp4']); end
initialFs = vidTop.FrameRate;
vidWriter = vision.VideoFileWriter(fileName, ...
    'AudioInputPort', true, ...
    'FrameRate', round(initialFs*s.playbackSpeed));
vidWriter.VideoCompressor = 'MJPEG Compressor';

% load spike data
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'obsOnTimes', 'obsOffTimes', 'rewardTimes', 'wiskContactFrames', 'isLightOn', ...
    'obsPixPositions', 'frameTimeStampsWisk', 'wiskContactTimes')

% get position where wisk frame should overlap with runTop frame
if s.includeWiskCam    
    obsInWiskCamInds = find(obsPixPositions>vidTop.Width-50 & obsPixPositions<vidTop.Width);
    
    % find first time point at which both wisk and run cams have a frame and obs is in wisk cam
    for i = obsInWiskCamInds
        wiskInd = find(frameTimeStampsWisk==frameTimeStamps(i));
        if ~isempty(wiskInd); topInd = i; break; end
    end

    frameTop = rgb2gray(read(vidTop, topInd));
    frameWisk = rgb2gray(read(vidWisk, wiskInd));
    [yWiskPos, xWiskPos, wiskScaling] = getSubFramePosition(frameTop, frameWisk, .35:.005 :.45);
    smpWiskFrame = imresize(frameWisk, wiskScaling);
end


% determine frame dimensions
if s.includeWiskCam
    frameDim = round([vidTop.Height + vidBot.Height, xWiskPos+size(smpWiskFrame,2)]);
else
    frameDim = [vidTop.Height + vidBot.Height, vidBot.Width];
end


% get neural data
ephysInfo = getSessionEphysInfo(session);
[bestChannels, unit_ids_all] = getBestChannels(session, ephysInfo);
bestChannel = bestChannels(unit_id==unit_ids_all);
getVoltage = @(data, channel, inds) data.Data.Data(channel,inds);
data = memmapfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysInfo.ephysFolder, [ephysInfo.fileNameBase '_CHs.dat']), ...
    'Format', {'int16', [ephysInfo.channelNum, ephysInfo.smps], 'Data'}, 'Writable', false);


% create low pass filter
if s.lowPassFreq
    lp = s.lowPassFreq * 2 / ephysInfo.fs;
    ls = s.lowPassFreq * 4 / ephysInfo.fs; % one octave above pass band
    ls = min(ls, .999); % make sure stopband freq doesn't exceed nyquist freqency
    [NN, Wn] = buttord(lp, ls, 3, 60);
    [B1,A1] = butter(NN, Wn, 'low');
end


% other initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'openEphysToSpikeMapping', 'spkRates', 'timeStamps', 'unit_ids', 'spkTimes')
timeStampsMapped = polyval(openEphysToSpikeMapping, ephysInfo.timeStamps);
audioSmpsPerFrame = round((1/initialFs) * ephysInfo.fs);


% set up figure
fig = figure('color', [0 0 0], 'position', [1925, 50, frameDim(2), frameDim(1)], 'menubar', 'none');
traceLength = s.voltageWindow*ephysInfo.fs;

axes('position', [0 .2 1 .8], 'CLim', [0 255]); colormap gray
im = image(uint8(zeros(frameDim)), 'CDataMapping', 'scaled'); hold on;
set(gca, 'Visible', 'off')

plotAxis = axes('position', [0 0 1 .2], 'Color', 'black');
tracePlot = plot(plotAxis, 1:traceLength, nan(1,traceLength), 'color', 'white'); set(gca, 'color', 'black'); hold on
set(gca, 'Visible', 'off', 'YLimMode', 'manual', 'YLim', s.yLims)

unitInd = find(unit_ids==unit_id);
scatter(spkTimes{unitInd}, ...
    repmat(s.yLims(1),1,length(spkTimes{unitInd})) + range(s.yLims)*.1, ... % y values
    10, s.spkScatterColor, 'filled');
obsOnLine = line(plotAxis, [0 0], s.yLims, 'linewidth', 2, 'color', s.lineColors);
obsOnText = text(plotAxis, 0, s.yLims(2), 'obstacle on', 'Color', 'white');
wiskLine = line(plotAxis, [0 0], s.yLims, 'linewidth', 2, 'color', s.lineColors);
wiskText = text(plotAxis, 0, s.yLims(2), 'whisker contact', 'Color', 'white');
rewardLine = line(plotAxis, [0 0], s.yLims, 'linewidth', 2, 'color', s.lineColors);
rewardText = text(plotAxis, 0, s.yLims(2), 'reward', 'Color', 'white');

if ~exist('timeEpochs', 'var'); timeEpochs = cat(2, obsOnTimes, obsOffTimes); end



% get min and max time for unit
minTime = timeStamps(find(~isnan(spkRates(unitInd,:)),1,'first'));
maxTime = timeStamps(find(~isnan(spkRates(unitInd,:)),1,'last'));
validTrials = find(timeEpochs(:,1)>minTime & timeEpochs(:,2)<maxTime); % make sure trials aren't too long
trials = validTrials(round(linspace(2,length(validTrials)-1,s.trialsToShow))); % remove first and last inds to avoid catching beginning and end of usable spikes period

% create video
fprintf('\nwriting video... ')
for i = 1:length(trials)
    
    % get trial frame indices
    trialInds = find(frameTimeStamps>timeEpochs(trials(i),1) & ...
                     frameTimeStamps<timeEpochs(trials(i),2));
                 
    % get voltage for entire trial
    voltageBins = timeStampsMapped>timeEpochs(trials(i),1)-s.voltageWindow-1 & ...
                  timeStampsMapped<timeEpochs(trials(i),2)+1; % add and subtract 1 as a buffer
    timeStampsSub = timeStampsMapped(voltageBins);
    voltageRaw = getVoltage(data, bestChannel, voltageBins); % maybe replace this by reading data within trial prior to writing each trials data...
    if s.lowPassFreq
        voltageRaw = filter(B1, A1, voltageRaw); % filter forwards then backwards to avoid phase shifts
        voltageRaw = int16(fliplr(filter(B1, A1, fliplr(voltageRaw))));
    end
    voltage = double(voltageRaw) * ephysInfo.bitVolts; % convert to microvolts
    
    
    % update obstacle and whisker contact lines
    if isLightOn(trials(i)); obsOnString = 'obstacle (light on)'; else; obsOnString = 'obstacle (light off)'; end
	updateTextAndLine(obsOnText, obsOnLine, obsOnTimes, obsOnString)
    updateTextAndLine(wiskText, wiskLine, wiskContactTimes)
    updateTextAndLine(rewardText, rewardLine, rewardTimes)
  
    % get frames for trials
    for j = trialInds'

        % get run frame
        frame = uint8(zeros(frameDim));
        topBotFrame = cat(1, rgb2gray(read(vidTop, j)), rgb2gray(read(vidBot, j)));
        topBotFrame = imadjust(topBotFrame, s.contrastLims, [0 1]); % adjust contrast
        frame(:,1:vidBot.Width) = topBotFrame;
        
        % get wisk frame
        if s.includeWiskCam
            wiskFrameInd = find(frameTimeStampsWisk==frameTimeStamps(j), 1, 'first');

            if ~isempty(wiskFrameInd)

                % get wisk frame
                frameWisk = rgb2gray(read(vidWisk, wiskFrameInd));

                % resize, adjust contrast, and draw border
                frameWisk = imresize(frameWisk, wiskScaling);
                frameWisk = imadjust(frameWisk, [.5 1], [0 1]);
                frameWisk = 255 - frameWisk;
                frameWisk([1:s.wiskBorder, end-s.wiskBorder:end], :) = 255;
                frameWisk(:, [1:s.wiskBorder, end-s.wiskBorder:end]) = 255;
                
                % incorporate into frame
                frame(yWiskPos:yWiskPos+size(frameWisk,1)-1, xWiskPos:xWiskPos+size(frameWisk,2)-1, :) = frameWisk;
            end
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
        elseif i==length(trials) && j==trialInds(end)
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
    [~,~] = system([fileName(1) ': & ffmpeg -i ' fileName ' -vcodec mpeg4 -vb 10M -y ' baseDir '\temp.avi']); % run ffmpeg to compress file
    delete(fileName)
    movefile([baseDir '\temp.avi'], fileName)
end

disp(' all done!')





% update position of lines and text marking events of interest
function updateTextAndLine(text, line, eventTimes, textString)
    eventInd = find(eventTimes>=timeEpochs(trials(i),1)-s.voltageWindow & ...
                    eventTimes<=timeEpochs(trials(i),2),1,'first');
    if ~isempty(eventInd)
        eventTime = eventTimes(eventInd);
        set(text, 'Position', [eventTime+s.voltageWindow*.01 s.yLims(2)])
        set(line, 'XData', [eventTime eventTime])
    end
    if exist('textString', 'var'); set(text, 'String', textString); end
end

end


