function makeUnitVid(session, unit_id, fileName)


% TO DO:
% consider if it's better to incorporate this into makeVidWisk or makeMatrixVid
% compress audio (matlab, or ffmpeg as backup)
% figure out how to select frames
% gate sound of spikes
% trace or otherwise mark spike times

% temp
startFrame = 250*60*15;
totalFrames = 2000;

% settings
initialFs = 250;
playbackSpeed = .25;
voltageWindow = .8;
audioGain = 15;
yLims = [-1000 800];

% set up video readers / writer
disp('initializing...')
vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
if ~exist('fileName', 'var')
    fileName = fullfile(getenv('OBSDATADIR'), 'editedVid', 'vidsWithNeurons', [session 'unit' num2str(unit_id) '.avi']);
end
vidWriter = vision.VideoFileWriter(fileName, ...
    'AudioInputPort', true, ...
    'FrameRate', round(initialFs*playbackSpeed));
vidWriter.VideoCompressor = 'MJPEG Compressor';


% get neural data
ephysInfo = getSessionEphysInfo(session);
[bestChannels, unit_ids] = getBestChannels(session, ephysInfo);
bestChannel = bestChannels(unit_id==unit_ids);
getVoltage = @(data, channel, inds) data.Data.Data(channel,inds);
data = memmapfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysInfo.ephysFolder, [ephysInfo.fileNameBase '_CHs.dat']), ...
    'Format', {'int16', [ephysInfo.channelNum, ephysInfo.smps], 'Data'}, 'Writable', false);
% voltage = getVoltage(data, bestChannel, 1:ephysInfo.smps); % !!! this is too slow

% other initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps')
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'openEphysToSpikeMapping')
timeStampsMapped = polyval(openEphysToSpikeMapping, ephysInfo.timeStamps);
audioSmpsPerFrame = round((1/initialFs) * ephysInfo.fs);


% set up figure
fig = figure('color', [0 0 0], 'position', [1925, 50, vidTop.width, (vidBot.Height + vidTop.Height)], 'menubar', 'none');
traceLength = voltageWindow*ephysInfo.fs;

imAxis = axes('position', [0 .2 1 .8], 'CLim', [0 255]); colormap gray
frame = cat(1, rgb2gray(read(vidTop,1)), rgb2gray(read(vidBot,1)));
im = image(frame, 'CDataMapping', 'scaled'); hold on;
set(gca, 'Visible', 'off')

plotAxis = axes('position', [0 0 1 .2], 'Color', 'black');
tracePlot = plot(plotAxis, 1:traceLength, nan(1,traceLength), 'color', 'white'); set(gca, 'color', 'black')
set(gca, 'Visible', 'off', 'YLimMode', 'manual', 'YLim', yLims / ephysInfo.bitVolts)



% create video
for i = startFrame:startFrame+totalFrames
    
    % get frame
    frame = cat(1, rgb2gray(read(vidTop, i)), rgb2gray(read(vidBot, i)));
    set(im, 'CData', frame);
    
    % get voltage
    traceStartInd = find(timeStampsMapped>(frameTimeStamps(i)-voltageWindow), 1, 'first');
%     trace = voltage(traceStartInd:traceStartInd+traceLength-1);
    trace = getVoltage(data, bestChannel, traceStartInd:traceStartInd+traceLength-1);
    set(tracePlot, 'ydata', trace);
    
    % get audio
    audioStartInd = find(timeStampsMapped>=frameTimeStamps(i), 1, 'first');
    audio = getVoltage(data, bestChannel, audioStartInd:audioStartInd+audioSmpsPerFrame-1)';
%     audio = voltage(audioStartInd:audioStartInd+audioSmpsPerFrame-1)';
    audio = audio * audioGain;
    
    % write to file
    frame = getframe(gcf);
    vidWriter(frame.cdata, audio);
    
end

release(vidWriter);
close(fig)
disp('all done!')



