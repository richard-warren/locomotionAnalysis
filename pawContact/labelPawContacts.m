function labelPawContacts(session)

% settings
vidSizeScaling = 1.25;
classNames = {'fore_ventral', 'fore_dorsal', 'hind_ventral_low', 'hind_ventral_high', 'hind_dorsal', 'no_touch', 'skip_frame'};
classNumbers = [1 4 2 5 8 3 9]; % these are the numbers to press on keypad to assign a frame to a particular class

% initializations
noTouchBin = strcmp(classNames, 'no_touch');
vidDelay = .02;
playing = true;
paused = false;
currentFrameInd = 1;
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
sampleFrame = rgb2gray(read(vid,currentFrameInd));
sampleFrameBot = rgb2gray(read(vidBot,currentFrameInd));
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'frameTimeStamps', 'obsPixPositions', 'obsOnTimes', 'obsOffTimes');
keypadLayout = [7 8 9; 4 5 6; 1 2 3];


% load previous data if it exists
fileName = [getenv('OBSDATADIR') 'tracking\trainingData\pawContact\' session '_contacts.mat'];
if exist(fileName, 'file')
    load(fileName, 'classes', 'classNames')
    fprintf('%s: loaded previous data\n', session)
else
    classes = nan(length(classNames), length(frameTimeStamps));
end

% get trial inds
% frameInds = getFramesToShow(session, true);
frameInds = [];
for i = 1:length(obsOnTimes)
    trialInds= find(frameTimeStamps>obsOnTimes(i) & frameTimeStamps<obsOffTimes(i) & ...
                     obsPixPositions'>0 & obsPixPositions'<=vid.Width);
    if ~isempty(trialInds)
        frameInds = [frameInds; trialInds];
    end
end
trialStartInds = find([0; diff(frameTimeStamps(frameInds))]>1);


% prepare progress bar figure
figProgress = figure('name', 'progres...', 'units', 'pixels', 'position', [601 326 898 40], 'menubar', 'none');
progressIm = imagesc(zeros(1,length(frameInds)));
progressLine = line([1 1], get(gca,'ylim'), 'linewidth', 5);
set(gca, 'Position', [0 0 1 1], 'Visible', 'off');

% prepare category images figures
imgDim = 200;
tileImg = zeros(imgDim*3, imgDim*3, 3);
tileImgMask = ones(length(classNames), imgDim*3, imgDim*3, 3);

for i = 1:length(classNames)
    img = imread([getenv('OBSDATADIR') 'tracking\trainingData\pawContact\exampleImgs\' classNames{i} '.jpeg']);
    img = double(imresize(img, [imgDim imgDim]))/255;
    img = insertText(img, [0,0], classNames{i}, 'boxcolor', [1 1 1], 'boxopacity', 1, 'fontsize', 18);
    
    [row, col] = ind2sub(size(keypadLayout), find(keypadLayout==classNumbers(i)));
    rowInds = (row-1)*imgDim+1:(row)*imgDim;
    colInds = (col-1)*imgDim+1:(col)*imgDim;
    tileImg(rowInds, colInds, :) = img;
    tileImgMask(i, rowInds, colInds, 1) = 0;
end

figSamples = figure('name', 'classes', 'units', 'pixels', 'position', [1100 403 400 400], ...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
samplesIm = imagesc(tileImg); hold on;
set(gca, 'visible', 'off', 'position', [0 0 1 1]);

% prepare main figure
figMain = figure('name', session, 'units', 'pixels', 'position', [600 400 vid.Width*vidSizeScaling (vid.Height+vidBot.Height)*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray
preview = image([sampleFrame; sampleFrameBot], 'CDataMapping', 'scaled'); hold on;
set(gca, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vid.Width*vidSizeScaling (vid.Height+vidBot.Height)*vidSizeScaling]);






% main loop
while playing
    while paused; pause(.001); end
    updateFrame(1);
end
close(figMain, figProgress, figSamples)



% ---------
% FUNCTIONS
% ---------

% keypress controls
function changeFrames(~,~)
    
    key = double(get(figMain, 'currentcharacter'));
    
    if ischar(key) || isscalar(key)
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

            % ESCAPE: close window
            case 27
                playing = false;
                paused = false;

            % 'f': select frame
            case 102
                pause(.001);
                paused = true;
                input = inputdlg('enter frame number');
                currentFrameInd = find(frameInds>=str2double(input{1}),1,'first');
                updateFrame(1);
            
            % 't': skip to next trial
            case 116
                pause(.001);
                paused = true;
                currentFrameInd = trialStartInds(find(trialStartInds>currentFrameInd,1,'first'));
                updateFrame(0);

            % 's': save current progress
            case 115
                save(fileName, 'classes', 'classNames')
                m = msgbox('saving...'); pause(.5); close(m)

            % SPACEBAR: pause or resume playback
            case 32
                paused = ~paused;

            % asign to category
            otherwise
                keypadNum = str2double(char(key));
                if any(classNumbers==keypadNum)
                    classBin = (classNumbers == keypadNum);
                    classes(:,frameInds(currentFrameInd)) = classBin;
                    updateFrame(0);
                end
        end
    end
end



% update frame preview
function updateFrame(frameStep)
    
    % update frame index
    currentFrameInd = currentFrameInd + frameStep;
    if currentFrameInd > length(frameInds); currentFrameInd = length(frameInds);
    elseif currentFrameInd < 1; currentFrameInd = 1; end
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid, frameInds(currentFrameInd)));
    frameBot = rgb2gray(read(vidBot, frameInds(currentFrameInd)));
    frame = [frame; frameBot];
    frame = insertText(frame, [size(frame,2) size(frame,1)], ...
        sprintf('frame %i', frameInds(currentFrameInd)),...
	    'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
    
    % update main figure
    set(preview, 'CData', frame);
    
    % update progress bar
    set(progressIm, 'Cdata', any(isnan(classes(:,frameInds)),1));
    set(progressLine, 'XData', [currentFrameInd currentFrameInd])
    
    % set current frame to 'no touch' if it has not yet been analyzed
    if all(isnan(classes(:, frameInds(currentFrameInd))))
        classes(:, frameInds(currentFrameInd)) = noTouchBin;
    end
    
    % update sample img
    currentClass = find(classes(:, frameInds(currentFrameInd)));
    set(samplesIm, 'Cdata', tileImg .* squeeze(tileImgMask(currentClass,:,:,:)));
    
    % pause to reflect on the little things...
    pause(vidDelay);
end



end