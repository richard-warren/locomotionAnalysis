%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% global settings



% global initializations
mice = {data.data.mouse};


%% ----------
% PLOT THINGS
%  ----------

%% GET DISTANCE AND TIME TO CONTACT DATA

% settings
trialSmps = 100;

% initializations
[distances, times] = deal(cell(1,length(mice)));

for i = 1:length(mice)
    
    fprintf('%s: collectin data, ya heard...\n', mice{i})
    sessions = {data.data(i).sessions.session}; % sessions for mouse
    
    [distances{i}, times{i}] = deal([]);
    for j = 1:length(sessions)
        
        % load session data
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{j}, 'kinData.mat'), 'kinData')
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{j}, 'runAnalyzed.mat'), 'frameTimeStamps', 'isLightOn');
        secondsPerFrame = nanmedian(diff(frameTimeStamps)); % seconds per frame
        
        % set conditionals here
        modPawOnlySwing = false(size(kinData)); % only trials where mod paw is in swing and non mod paw is in stance
        modPawOnlySwing([kinData.isTrialAnalyzed]) = ([kinData.isLeftSwingAtContact] + [kinData.isRightSwingAtContact]) == 1;
        bins = [kinData.isTrialAnalyzed] & modPawOnlySwing & ~isLightOn(:)';
        
        for k = find(bins)
            
            % get distance of leading paw at contact
            contactInd = find(frameTimeStamps(kinData(k).trialInds) >= kinData(k).wiskContactTimes, 1, 'first'); % ind in trial at which contact occurs
            distances{i}(end+1) = abs(max([kinData(k).locations(contactInd,1,:)]))*1000;

            trialX = max(kinData(k).locations(contactInd-trialSmps+1:contactInd,1,:), [], 3);
            linFit = polyfit(trialX', 1:trialSmps, 1);
            predictedAtObsInd = polyval(linFit, 0);
            times{i}(end+1) = abs((predictedAtObsInd-trialSmps) * secondsPerFrame)*1000; % frame until contact * (seconds/frame)
            
        end
        clear kinData frameTimeStamps isLightOn
    end
end
disp('all done!')

%% PLOT DISTANCE AND TIME TO CONTACT

% settings
xLims = [15 50];
yLims = [0 150];
scatSize = 4;
scatAlpha = .5;
mouseColors = true;
scatPlotSize = .7;
border = .15;


% initializations
d = abs(cat(2, distances{:}));
t = abs(cat(2, times{:}));
mouseIds = repelem(1:length(mice), cellfun(@length, distances));
medDistances = cellfun(@nanmedian, distances);
medTimes = cellfun(@nanmedian, times);

% plot that shit
figure('Color', 'white', 'Position', [2000 400 450 350], 'MenuBar', 'none');
scatterHistoRick(d,t, ...
    {'groupId', mouseIds, 'colors', 'jet', ...
    'xlabel', 'distance to contact (mm)', 'ylabel', 'time to contact (ms)', ...
    'xLims', xLims, 'yLims', yLims, 'showCrossHairs', true, 'scatSize', scatSize, 'scatAlpha', scatAlpha});


% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_distanceTimeToContact');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

%% SCATTER VEL AND BODY ANGLE AT CONTACT

flat = flattenData(data, {'mouse', 'session', 'trial', 'velAtWiskContact', 'angleAtWiskContact', 'isLightOn', 'modPawOnlySwing'});
flat = flat(~[flat.isLightOn]);
flat = flat([flat.modPawOnlySwing]);
[~,~,mouseIds] = unique({flat.mouse});

% settings
xLims = [0 1];
yLims = [-30 30];

% initializations
figure('Color', 'white', 'Position', [2000 400 450 350], 'MenuBar', 'none');
scatterHistoRick([flat.velAtWiskContact], [flat.angleAtWiskContact], ...
    {'groupId', mouseIds, 'colors', 'jet', ...
    'xlabel', 'velocity (m/s)', 'ylabel', ['body angle (' char(176) ')'], ...
    'xLims', xLims, 'yLims', yLims, 'showCrossHairs', true, 'scatSize', scatSize, 'scatAlpha', scatAlpha});

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_phaseVelVariability');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% OVERLAY IMAGES

% settings
session = '180628_004';
trialsToOverlay = 10;

vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'obsOnTimes', 'wiskContactTimes');
wiskContactTimes = wiskContactTimes(~isnan(wiskContactTimes));
frameTimes = sort(wiskContactTimes(randperm(length(wiskContactTimes), trialsToOverlay)));

imgsTop = uint8(zeros(vidTop.Height, vidTop.Width, trialsToOverlay));
imgsBot = uint8(zeros(vidBot.Height, vidBot.Width, trialsToOverlay));

for i = 1:trialsToOverlay
    imgsTop(:,:,i) = rgb2gray(read(vidTop, knnsearch(frameTimeStamps, frameTimes(i))));
    imgsBot(:,:,i) = rgb2gray(read(vidBot, knnsearch(frameTimeStamps, frameTimes(i))));
end

%
overlayTop = overlayImgs(imgsTop, {'colors', 'jet', 'contrastLims', [.3 .75], 'cutoff', 100, 'projection', 'mean'});
overlayBot = overlayImgs(imgsBot, {'colors', 'jet', 'contrastLims', [.3 .5], 'cutoff', 100, 'projection', 'mean'});
overlay = cat(1,overlayTop, overlayBot);
figure; imshow(overlay); pimpFig


% write image to desk
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', ...
        'variability_omg.png');
fprintf('writing %s to disk...\n', file)
imwrite(overlay, file);

%% KINEMATICS WITH LANDING POSITION DISTRIBUTIONS

% settings
trialsToShow = 50;
xLims = [-.1 .07];
stepTypeColors = flipud([0.850 0.325 0.098; 0 0.447 0.741]); % first entry is big step, second is small step
ctlColor = [.2 .2 .2];
histoFillAlpha = .2;



% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', ...
    'modPawKinInterp', 'preModPawKinInterp', 'obsHgt', 'isLightOn', 'modPawOnlySwing', 'isBigStep'});
flat = flat(~[flat.isLightOn]);
flat = flat([flat.modPawOnlySwing]);
[~,~,mouseIds] = unique({flat.mouse});
kinData = permute(cat(3, flat.modPawKinInterp), [3,1,2]);
kinDataCtl = permute(cat(3, flat.preModPawKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1);
condition = ~[flat.isBigStep] + 1;
figure('Color', 'white', 'Position', [2000 400 900 250], 'MenuBar', 'none');




% plot kinematics
subplot(2,1,1)
plotKinematics(kinData(:,[1 3],:), [flat.obsHgt], condition, ...
    {'colors', stepTypeColors, 'trialsToOverlay', trialsToShow, 'trialAlpha', .4, 'lineAlpha', 0, 'yLimZero', false})
plotKinematics(kinDataCtl(:,[1 3],:), [flat.obsHgt], ones(1,length(flat)), ...
    {'colors', [.2 .2 .2], 'lineWidth', 5, 'yLimZero', false})


set(gca, 'XLim', xLims)


% plot pdfs
subplot(2,1,2); hold on;

xGrid = linspace(xLims(1), xLims(2), 500);
longShortRatio = sum(condition==1)/length(condition);

kdCtl = ksdensity(kinDataCtl(:,1,end), xGrid);
kdLong = ksdensity(kinData(condition==1,1,end), xGrid) * longShortRatio;
kdShort = ksdensity(kinData(condition==2,1,end), xGrid) * (1-longShortRatio);


% plot that shit
fill([xGrid xGrid(1)], [kdCtl kdCtl(1)], ctlColor, 'FaceAlpha', histoFillAlpha)
plot(xGrid, kdCtl, 'Color', ctlColor, 'LineWidth', 4)

fill([xGrid xGrid(1)], [kdLong kdLong(1)], stepTypeColors(1,:), 'FaceAlpha', histoFillAlpha)
plot(xGrid, kdLong, 'Color', stepTypeColors(1,:), 'LineWidth', 4)

fill([xGrid xGrid(1)], [kdShort kdShort(1)], stepTypeColors(2,:), 'FaceAlpha', histoFillAlpha)
plot(xGrid, kdShort, 'Color', stepTypeColors(2,:), 'LineWidth', 4)

set(gca, 'XLim', xLims, 'YDir', 'reverse', 'Visible', 'off')



% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_modStepVariability');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');












