function makeUnitVid(session, unit_id, timeEpochs)


% notes: make video with mouse on top and simultanneous activity of unit on
% bottom // timeEpochs is nx2 matrix, where each row is begenning and end
% time of a trial // plots evenly spaced trialsToShow of all possible
% trials // if timeEpochs not provided, then assumed to be obsOnTimes and
% obsOffTimes

% TO DO:
% accept inds for specific trials
% accept time epochs to use intead of obs on/off (so i can lot reward times, lick epochs, etc)
% fade in/out sound
% gate sound of spikes
% add paw markers?

% temp
% session = '181019_002';
% unit_id = 67;

% settings
trialsToShow = 5;
timePrePost = [-.2 0]; % time before obsOn and after obOff to show
maxTrialTime = 2.5; % trials exceeding maxTrialTime will be trimmed to this duration (s)
contrastLims = [.1 .9]; % pixels at these proportional values are mapped to 0 and 255
playbackSpeed = .15;
voltageWindow = .8;
audioGain = 15;
yLims = [-800 800]; %[-400 400];
wiskBorder = 2;

includeWiskCam = true;
compressVideo = true;

spkScatterColor = [1 1 0];
lineColors = [.8 .4 1];


% set up video readers / writer
disp('initializing...')
vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
if includeWiskCam; vidWisk = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runWisk.mp4']); end
initialFs = vidTop.FrameRate;
fileName = fullfile(getenv('OBSDATADIR'), 'editedVid', 'vidsWithNeurons', [session 'unit' num2str(unit_id) '.avi']);
vidWriter = vision.VideoFileWriter(fileName, ...
    'AudioInputPort', true, ...
    'FrameRate', round(initialFs*playbackSpeed));
vidWriter.VideoCompressor = 'MJPEG Compressor';

% load spike data
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'obsOnTimes', 'obsOffTimes', 'rewardTimes', 'wiskContactFrames', 'isLightOn', 'obsPixPositions', 'frameTimeStampsWisk')

% get position where wisk frame should overlap with runTop frame
if includeWiskCam    
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
if includeWiskCam
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
% voltage = getVoltage(data, bestChannel, 1:ephysInfo.smps); % !!! this is too slow

% other initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'openEphysToSpikeMapping', 'spkRates', 'timeStamps', 'unit_ids', 'spkTimes', 'rewardTimes')
timeStampsMapped = polyval(openEphysToSpikeMapping, ephysInfo.timeStamps);
audioSmpsPerFrame = round((1/initialFs) * ephysInfo.fs);
wiskContactTimes = frameTimeStampsWisk(wiskContactFrames(wiskContactFrames>0));


% set up figure
fig = figure('color', [0 0 0], 'position', [1925, 50, frameDim(2), frameDim(1)], 'menubar', 'none');
traceLength = voltageWindow*ephysInfo.fs;

axes('position', [0 .2 1 .8], 'CLim', [0 255]); colormap gray
im = image(uint8(zeros(frameDim)), 'CDataMapping', 'scaled'); hold on;
set(gca, 'Visible', 'off')

plotAxis = axes('position', [0 0 1 .2], 'Color', 'black');
tracePlot = plot(plotAxis, 1:traceLength, nan(1,traceLength), 'color', 'white'); set(gca, 'color', 'black'); hold on
set(gca, 'Visible', 'off', 'YLimMode', 'manual', 'YLim', yLims)

unitInd = find(unit_ids==unit_id);
scatter(spkTimes{unitInd}, ...
    repmat(yLims(1),1,length(spkTimes{unitInd})) + range(yLims)*.1, ... % y values
    10, spkScatterColor, 'filled');
obsOnLine = line(plotAxis, [0 0], yLims, 'linewidth', 2, 'color', lineColors);
obsOnText = text(plotAxis, 0, yLims(2), 'obstacle on', 'Color', 'white');
wiskLine = line(plotAxis, [0 0], yLims, 'linewidth', 2, 'color', lineColors);
wiskText = text(plotAxis, 0, yLims(2), 'whisker contact', 'Color', 'white');
rewardLine = line(plotAxis, [0 0], yLims, 'linewidth', 2, 'color', lineColors);
rewardText = text(plotAxis, 0, yLims(2), 'reward', 'Color', 'white');

if ~exist('timeEpochs', 'var'); timeEpochs = cat(2, obsOnTimes, obsOffTimes); end



% get min and max time for unit
minTime = timeStamps(find(~isnan(spkRates(unitInd,:)),1,'first'));
maxTime = timeStamps(find(~isnan(spkRates(unitInd,:)),1,'last'));
validTrials = find(timeEpochs(:,1)>minTime & ...
                   timeEpochs(:,2)<maxTime & ...
                   (diff(timeEpochs,1,2)) < maxTrialTime); % make sure trials aren't too long
trials = validTrials(round(linspace(2,length(validTrials)-1,trialsToShow))); % remove first and last inds to avoid catching beginning and end of usable spikes period


% create video
disp('writing video...')
for i = 1:length(trials)
    
    % get trial frame indices
    trialInds = find(frameTimeStamps>timeEpochs(trials(i),1)+timePrePost(1) & ...
                     frameTimeStamps<timeEpochs(trials(i),2)+timePrePost(2));
    
    % update obstacle and whisker contact lines
    if isLightOn(trials(i)); obsOnString = 'obstacle (light on)'; else; obsOnString = 'obstacle (light off)'; end
	updateTextAndLine(obsOnText, obsOnLine, obsOnTimes, obsOnString)
    updateTextAndLine(wiskText, wiskLine, wiskContactTimes)
    updateTextAndLine(rewardText, rewardLine, rewardTimes)
  
    % get frames for trials
    for j = trialInds'

        % get frame
        frame = uint8(zeros(frameDim));
        topBotFrame = cat(1, rgb2gray(read(vidTop, j)), rgb2gray(read(vidBot, j)));
        topBotFrame = imadjust(topBotFrame, contrastLims, [0 1]); % adjust contrast
        frame(:,1:vidBot.Width) = topBotFrame;
        
        % wisk
        if includeWiskCam
            wiskFrameInd = find(frameTimeStampsWisk==frameTimeStamps(j), 1, 'first');

            if ~isempty(wiskFrameInd)

                % get wisk frame
                frameWisk = rgb2gray(read(vidWisk, wiskFrameInd));

                % resize, adjust contrast, and draw border
                frameWisk = imresize(frameWisk, wiskScaling);
                frameWisk = imadjust(frameWisk, [.5 1], [0 1]);
                frameWisk = 255 - frameWisk;
                frameWisk([1:wiskBorder, end-wiskBorder:end], :) = 255;
                frameWisk(:, [1:wiskBorder, end-wiskBorder:end]) = 255;
                
                % incorporate into frame
                frame(yWiskPos:yWiskPos+size(frameWisk,1)-1, xWiskPos:xWiskPos+size(frameWisk,2)-1, :) = frameWisk;
            end
        end
        
        
        set(im, 'CData', frame);

        % get voltage
        traceStartInd = find(timeStampsMapped>(frameTimeStamps(j)-voltageWindow), 1, 'first');
        traceInds = traceStartInd:traceStartInd+traceLength-1;
    %     trace = voltage(traceStartInd:traceStartInd+traceLength-1);
        trace = double(getVoltage(data, bestChannel, traceInds)) * ephysInfo.bitVolts;
        times = timeStampsMapped(traceInds)';
        set(tracePlot, 'xdata', times, 'ydata', trace);
        set(gca, 'xlim', [times(1) times(end)])

        % get audio
        audioStartInd = find(timeStampsMapped>=frameTimeStamps(j), 1, 'first');
        audio = getVoltage(data, bestChannel, audioStartInd:audioStartInd+audioSmpsPerFrame-1)';
    %     audio = voltage(audioStartInd:audioStartInd+audioSmpsPerFrame-1)';
        audio = audio * audioGain;

        % write to file
        frame = getframe(fig);
        vidWriter(frame.cdata, audio);
    end
end
release(vidWriter);
close(fig)


% compress video
if compressVideo
    baseDir = fileparts(fileName);
    [~,~] = system([fileName(1) ': & ffmpeg -i ' fileName ' -vcodec mpeg4 -vb 10M -y ' baseDir '\temp.avi']); % run ffmpeg to compress file
    delete(fileName)
    movefile([baseDir '\temp.avi'], fileName)
end

disp('all done!')





% update position of lines and text marking events of interest
function updateTextAndLine(text, line, eventTimes, textString)
    eventInd = find(eventTimes>=timeEpochs(trials(i),1)+timePrePost(1)-voltageWindow & ...
                    eventTimes<=timeEpochs(trials(i),2)+timePrePost(2),1,'first');
    if ~isempty(eventInd)
        eventTime = eventTimes(eventInd);
        set(text, 'Position', [eventTime+voltageWindow*.01 yLims(2)])
        set(line, 'XData', [eventTime eventTime])
    end
    if exist('textString', 'var'); set(text, 'String', textString); end
end

end


