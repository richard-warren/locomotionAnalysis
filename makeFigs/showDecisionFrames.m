function showDecisionFrames(session, varargin)

% for session, shows frames at whisker contact and at the end of the
% modified step for one big and one little step trial // picks trials where
% predictedDistanceToObs is similar (state of body is similar), but choice
% is different


% settings
s.controlColor = [.8 .8 .8];
s.stepColors = lines(2);  % little and big step colors
s.yMax = [];  % (pixels) use to crop out the bottom view if desired (by setting yMax to the height of the top view in pixels)
s.unravel = true;  % whether to add the movement of the wheel to the kinematics
s.contrastLims = [.2 .8]; % pixels at these proportional values are mapped to 0 and 255
s.contactColor = [1 0 0];  % color of scatter point at whisker contact
s.leftPadding = 0;  % add black to the left, which helps all the kinematics fit in for long steps
s.rightPadding = 0;  % add black to the right, which helps all the kinematics fit in for long steps
s.edgeFading = 50; % fading at the edges of images
s.lineWidth = 5;
s.drawKinematics = true;  % whether to draw kinematics on frame


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
imgNames = {'smallContact', 'smallFinal', 'bigContact', 'bigFinal'};
vars = {'isBigStep', 'firstModPaw', 'modPawKin', 'trial', 'isLightOn'};
data = getExperimentData(session, vars);
flat = flattenData(data, vars);
if ~exist(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4')); concatTopBotVids(session); end; pause(1);
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), 'kinData');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'frameTimeStamps', 'frameTimeStampsWisk');
if isempty(s.yMax); hgt = vid.Height; else; hgt = s.yMax; end
[frame, yWiskPos, xWiskPos, wiskScaling] = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, 10000);  % grab arbitrary frame with whisker added
wid = size(frame,2) + s.leftPadding + s.rightPadding;
imgDims = [hgt, wid];



% find distance travelled at moment of contact, then find big and small
% step trials that are as close as possible on this metric
distances = nan(1,length(kinData));
for i = find([kinData.isTrialAnalyzed])
    distances(i) = kinData(i).locationsPix(kinData(i).contactInd,1,3);
end
bigInds = find([flat.isBigStep]==1 & [flat.firstModPaw]==3 & ~[flat.isLightOn]);
smallInds = find([flat.isBigStep]==0 & [flat.firstModPaw]==3 & ~[flat.isLightOn]);
[inds, diffs] = knnsearch(distances(bigInds)', distances(smallInds)');
[~, minInd] = min(diffs);
trialBig = bigInds(inds(minInd));
trialSmall = smallInds(minInd);



trials = [trialSmall, trialBig];
figs = cell(1,4);

for t = 1:2
    
    trial = trials(t);
    
    % get inds for step (wrt trial frames)
    startInd = find(kinData(trial).modifiedStepIdentities(:,3)==1,1,'first');
    contactInd = kinData(trial).contactInd;
    finalInd = find(kinData(trial).modifiedStepIdentities(:,3)==1,1,'last');
    
    if s.unravel
        kinData(trial).locationsPix(:,1,:) = kinData(trial).locationsPix(:,1,:) - kinData(trial).obsPixPos';
        contactOffset = kinData(trial).obsPixPos(contactInd) + s.leftPadding;
        finalOffset = kinData(trial).obsPixPos(finalInd) + s.leftPadding;
    else
        contactOffset = s.leftPadding;
        finalOffset = s.leftPadding;
    end
    
    
    % read frames
    frameNum = kinData(trial).trialInds(finalInd);
    contactFrame = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, kinData(trial).trialInds(contactInd), ...
        'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, 'isPaddingWhite', false, 'edgeFading', s.edgeFading);
    finalFrame = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, kinData(trial).trialInds(finalInd), ...
        'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, 'isPaddingWhite', false, 'edgeFading', s.edgeFading);
    
    if ~isempty(s.yMax)
        contactFrame = contactFrame(1:s.yMax,:,:);
        finalFrame = finalFrame(1:s.yMax,:,:);
    end
    
    if s.leftPadding>0
        contactFrame = cat(2, zeros(imgDims(1), s.leftPadding), contactFrame);
        finalFrame = cat(2, zeros(imgDims(1), s.leftPadding), finalFrame);
    end
    
    if s.rightPadding>0
        contactFrame = cat(2, contactFrame, zeros(imgDims(1), s.rightPadding));
        finalFrame = cat(2, finalFrame, zeros(imgDims(1), s.rightPadding));
    end
    
    
    contactFrame = imadjust(contactFrame, s.contrastLims, [0 1]);
    finalFrame = imadjust(finalFrame, s.contrastLims, [0 1]);
    
    
    % whisker contact frame
    figs{(t-1)*2+1} = figure('name', imgNames{(t-1)*2+1}, 'Color', 'black', 'MenuBar', 'none');
    inds = startInd:contactInd;
    imshow(contactFrame); hold on
    
    if s.drawKinematics
        plot(kinData(trial).locationsPix(inds,1,3)+contactOffset, kinData(trial).locationsPix(inds,3,3), ...
             'lineWidth', s.lineWidth, 'Color', s.controlColor)
        plot(kinData(trial).locationsPix(inds,1,3)+contactOffset, kinData(trial).locationsPix(inds,2,3), ...
             'lineWidth', s.lineWidth, 'Color', s.controlColor)
        scatter(repmat(kinData(trial).locationsPix(contactInd,1,3)+contactOffset, 1, 2), ...
            kinData(trial).locationsPix(contactInd,2:3,3), 100, s.contactColor, 'filled')
    end
    set(gcf, 'Position', [2000 400 imgDims(2) imgDims(1)])
        set(gca, 'Position', [0 0 1 1])
    
    
    % step over frame
    figs{(t-1)*2+2} = figure('name', imgNames{(t-1)*2+2}, 'Color', 'black', 'MenuBar', 'none');
    inds = startInd:finalInd;
    imshow(finalFrame); hold on
    
    if s.drawKinematics
        plot(kinData(trial).locationsPix(inds,1,3)+finalOffset, kinData(trial).locationsPix(inds,3,3), ...
             'lineWidth', s.lineWidth, 'Color', s.stepColors(t,:))
        plot(kinData(trial).locationsPix(inds,1,3)+finalOffset, kinData(trial).locationsPix(inds,2,3), ...
             'lineWidth', s.lineWidth, 'Color', s.stepColors(t,:))

        % draw pre contact kinematics in different color over the initial part of the trace
        inds = startInd:contactInd;
        plot(kinData(trial).locationsPix(inds,1,3)+finalOffset, kinData(trial).locationsPix(inds,3,3), ...
             'lineWidth', s.lineWidth, 'Color', s.controlColor)
        plot(kinData(trial).locationsPix(inds,1,3)+finalOffset, kinData(trial).locationsPix(inds,2,3), ...
             'lineWidth', s.lineWidth, 'Color', s.controlColor)
    
        % add circle at moment of whisker contact
        scatter(repmat(kinData(trial).locationsPix(contactInd,1,3)+finalOffset, 1, 2), ...
            kinData(trial).locationsPix(contactInd,2:3,3), 100, s.contactColor, 'filled')
    end
    set(gcf, 'Position', [2000 400 imgDims(2) imgDims(1)])
    set(gca, 'Position', [0 0 1 1])
end


% save images
for i = 1:4
    file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', ['decision_' imgNames{i} '.png']);
    saveas(figs{i}, file)
end
disp('all done!')






