%% INITIALIZATIONS

% load global settings
clear all; tcm190910_config;

% settings
referenceModPaw = 2;
velSmps = 10; % how many samples to compute velocity over
frameRate = 250;
useReferencePaw = true; % if true, flip everything with respect to referenceModPaw

% load experiment data
fprintf('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% flatten data and compute predictors
flat = flattenData(data, {'mouse', 'session', 'trial', 'isTrialSuccess', 'modPawOnlySwing', ...
    'firstModPaw', 'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'isBigStep', ...
    'modPawKinInterp', 'preModPawKinInterp', 'isLightOn', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs'});
flat = struct2table(flat);

% flip predictors in flat so everything is relative to firstModPaw
flipBins = [flat.firstModPaw]==referenceModPaw;
flat.angleAtWiskContact(flipBins) = -flat.angleAtWiskContact(flipBins);


% add predictors not in flat already
sessions = unique(flat.session);

for i = 1:length(sessions)
    
    fprintf('%s: computing extra predictor variables\n', sessions{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'kinData.mat'), 'kinData')
    sessionBins = strcmp(flat.session, sessions{i});
    if sum(sessionBins)~=length(kinData); disp('WTF! mismatch in number of trials! fucking shit!!!'); end
    
    % loop over paws
    for j = find([kinData.isTrialAnalyzed])
        
        % flip everything relative to first modified paw
        pawSequence = [1 2 3 4];
        if useReferencePaw && kinData(j).firstModPaw~=referenceModPaw; pawSequence = [4 3 2 1]; end
        flatBin = sessionBins & [flat.trial]==j;
        
        for k = 1:4
            
            % starting paw position
            flat.(['modStepStart_paw' num2str(k)])(flatBin) = kinData(j).modifiedLocations{pawSequence(k)}(1,1,1);
            
            % is in stance
            flat.(['isStance_paw' num2str(k)])(flatBin) = kinData(j).stanceBins(kinData(j).contactInd,pawSequence(k));
            
            % x and z positions
            flat.(['x_paw' num2str(k)])(flatBin) = kinData(j).locations(kinData(j).contactInd, 1, pawSequence(k));
            flat.(['z_paw' num2str(k)])(flatBin) = kinData(j).locations(kinData(j).contactInd, 3, pawSequence(k));
            
            % x and z velocity
            kin = kinData(j).locations(kinData(j).contactInd-velSmps+1:kinData(j).contactInd, :, pawSequence(k));
            flat.(['xVel_paw' num2str(k)])(flatBin) = (kin(end,1)-kin(1,1)) / velSmps / frameRate;
            flat.(['zVel_paw' num2str(k)])(flatBin) = (kin(end,3)-kin(1,3)) / velSmps / frameRate;
        end
        
        % ind in modPawKinInterp at which whiskers contact obs
        flat.contactInd(flatBin) = kinData(j).pawObsPosIndInterp(kinData(j).firstModPaw);
    end
    
    % remove unanalyzed trials
    validBins = ~(sessionBins & ismember(flat.trial, find(~[kinData.isTrialAnalyzed]))); % all trial not in ssession OR in session but analyzed
    flat = flat(validBins,:);
end
disp('all done!')



% BUILD MODELS

% settings
iterations = 20;
normalizeData = true;
barColors = [0 .24 .49; .6 .75 .23; .2 .2 .2];

% prepare training data
[X, y, predictorNames, isCategorical] = ...
    prepareDecisionModelData(flat, 'all', 'isBigStep', true, referenceModPaw, normalizeData, {'balanceClasses', true});
predDistBin = ismember(predictorNames, 'modPawPredictedDistanceToObs'); % this var will ONLY be included in the GLM, not the neural net
    
    
accuracies = nan(3,iterations);
for i = 1:iterations
    
    % NEURAL NET
    net = patternnet(100);
    net.divideParam.trainRatio = .7;
    net.divideParam.valRatio = .15;
    net.divideParam.testRatio = .15;
    [net, tr] = train(net, X(:,~predDistBin)', y');
    outputs = net(X(:,~predDistBin)');
    accuracies(1,i) = mean(round(outputs(tr.testInd))==y(tr.testInd)');
    
    % NN SHUFFLED
    shuffleInds = randperm(length(tr.testInd));
    accuracies(3,i) = mean(round(outputs(tr.testInd))==y(tr.testInd(shuffleInds))'); % NN shuffled
    
    % GLM
    glmPredictorBins = true(1,size(X,2));  % bins of predictors to include in GLM (columns of X)
    glm = fitglm(X([tr.trainInd tr.valInd],glmPredictorBins), y([tr.trainInd tr.valInd]), ...
        'Distribution', 'binomial', 'CategoricalVars', isCategorical(glmPredictorBins));
    accuracies(2,i) = mean(round(predict(glm, X(tr.testInd,glmPredictorBins)))==y(tr.testInd));
    
end


%% MODEL ACCURACY BARS

figure('units', 'inches', 'position', [31.27 2.36 2.8 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
barFancy(accuracies([1 3],:), 'levelNames', {{'ANN', 'shuffled'}}, 'ylabel', 'accuracy', 'colors', modelColors, ...
    'YTick', [0 .5 1], 'YLim', [0 1], barProperties{:})
set(gca, 'position', [.2 .1 .77 .81]);
print -clipboard -dmeta


%% MOD PAW TRIAL OVERLAYS

% settings
trialsToShow = 50;
xLims = [-.1 .07];
histoFillAlpha = .2;
histoHgt = .3;

% initializations
close all
figure('units', 'inches', 'position', [22.76 1.51 7.82 2.31], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
flat_sub = flattenData(data, {'mouse', 'session', 'trial', ...
    'modPawKinInterp', 'preModPawKinInterp', 'obsHgt', 'isLightOn', 'modPawOnlySwing', 'isBigStep'});
flat_sub = flat_sub(~[flat.isLightOn]);
flat_sub = flat_sub([flat_sub.modPawOnlySwing]);
[~,~,mouseIds] = unique({flat_sub.mouse});
kinData = permute(cat(3, flat_sub.modPawKinInterp), [3,1,2]);
kinDataCtl = permute(cat(3, flat_sub.preModPawKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1);
condition = ~[flat_sub.isBigStep] + 1;

% plot kinematics
axes('position', [0 histoHgt 1 1-histoHgt]); hold on;

plotKinematics(kinData(:,[1 3],:), [flat_sub.obsHgt], condition, ...
    {'colors', stepTypeColors, 'trialsToOverlay', trialsToShow, 'trialAlpha', .4, 'lineAlpha', 0, 'yLimZero', false, 'lineColor', axisColor})
plotKinematics(kinDataCtl(:,[1 3],:), [flat_sub.obsHgt], ones(1,length(flat_sub)), ...
    {'colors', axisColor, 'lineWidth', 5, 'yLimZero', false, 'lineColor', 'none'})
set(gca, 'color', 'black', 'xlim', xLims)

% plot pdfs
axes('position', [0 .02 1 histoHgt-.02]); hold on;

xGrid = linspace(xLims(1), xLims(2), 500);
longShortRatio = sum(condition==1)/length(condition);

kdCtl = ksdensity(kinDataCtl(:,1,end), xGrid);
kdLong = ksdensity(kinData(condition==1,1,end), xGrid) * longShortRatio;
kdShort = ksdensity(kinData(condition==2,1,end), xGrid) * (1-longShortRatio);

fill([xGrid xGrid(1)], [kdCtl kdCtl(1)], axisColor, 'FaceAlpha', histoFillAlpha, 'EdgeColor', 'none')
plot(xGrid, kdCtl, 'Color', axisColor, 'LineWidth', 2)
fill([xGrid xGrid(1)], [kdLong kdLong(1)], stepTypeColors(1,:), 'FaceAlpha', histoFillAlpha, 'EdgeColor', 'none')
plot(xGrid, kdLong, 'Color', stepTypeColors(1,:), 'LineWidth', 2)
fill([xGrid xGrid(1)], [kdShort kdShort(1)], stepTypeColors(2,:), 'FaceAlpha', histoFillAlpha, 'EdgeColor', 'none')
plot(xGrid, kdShort, 'Color', stepTypeColors(2,:), 'LineWidth', 2)

set(gca, 'XLim', xLims, 'YDir', 'reverse', 'Visible', 'off', 'color', 'black')
print -clipboard -dmeta

%% BINNED KINEMATICS

% settings
close all
binNum = 5;
pctileBins = true;

% initializations
kin = permute(cat(3,flat.modPawKinInterp{:}), [3,1,2]);
kinCtl = permute(cat(3,flat.preModPawKinInterp{:}), [3,1,2]);
[X, y, predictorNames, isCategorical] = ...
        prepareDecisionModelData(flat, 'all', 'isBigStep', true, referenceModPaw, normalizeData, ...
        {'balanceClasses', false, 'removeNans', false});  % must recom 

    % choose binning variable
binVar = net(X(:,~predDistBin)'); % neural network output

% perctentile bins    
if pctileBins
    binEdges = prctile(binVar, linspace(0, 100, binNum+1));
% evenly spaced bins
else
    binEdges = linspace(min(binVar), max(binVar), binNum+1);
end
condition = discretize(binVar, binEdges);

figure('units', 'inches', 'position', [22.76 1.51 7.82 4.64], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
plotBigStepKin(kin(:,[1,3],:), kinCtl(:,[1,3],:), flat.obsHgt, condition, flat.isBigStep, ...
    {'colors', stepTypeColors, 'xLims', xLims, 'addHistos', false, 'lineWid', 3, 'controlColor', axisColor, ...
    'contactInds', flat.contactInd, 'histoHgt', .5, 'showSmpNum', false, 'obsColor', axisColor})
print -clipboard -dmeta


%% HEATMAP

% settings
xLims = [-.03 .015]*1000;
yLims = [-.03 .03]*1000;

close all
figure('units', 'inches', 'position', [22.76 1.51 7.80 4.64], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
heatmapRick(flat.modPawPredictedDistanceToObs*1000, flat.modPawDistanceToObs*1000, ...
    {'xLims', xLims, 'yLims', yLims, 'xlabel', 'predicted distance (mm)', 'ylabel', 'distance (mm)', 'colormap', 'hot'})
set(gca, 'DataAspectRatio', [1 1 1])
line(xLims, xLims, 'color', axisColor, 'lineWidth', 3)
set(gca, 'color', 'black', 'xcolor', axisColor, 'ycolor', axisColor, 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta

%% NEURAL NET DIAGRAM

% settings
layers = [5 15 1]; % number of units in each layer
dots = {[ceil(layers(1)/2)], [1 0 1]+[ceil(layers(2)/2)], []};
layerHgts = [.3 .8 .5];  % height of layers, expressed as a fraction of figure height
color = modelColors(1,:);
accentColor = axisColor;
nodeSize = 300;  % expressed as fraction of y axis
dotSize = 20;
lineWidth = 2;
lineAlpha = .4;
x = linspace(.1, .9, length(layers));

close all
figure('units', 'inches', 'position', [7.82 3.84 3.22 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
hold on

% get centers of all nodes
centers = cell(1,length(layers));
for i = 1:length(layers)
    if layers(i)>1
        centers{i} = linspace(0, layerHgts(i), layers(i))+.5-layerHgts(i)/2;
    else
        centers{i} = .5;
    end
end


for i = 1:length(layers)
    for j = 1:layers(i)
        
        % connect the dots
        if i<length(layers)
            for k = 1:layers(i+1)
                plot([x(i) x(i+1)], [centers{i}(j) centers{i+1}(k)], ...
                    'lineWidth', lineWidth, 'color', [accentColor lineAlpha]);
            end
        end
        
        % draw node
        if ismember(j, dots{i})
            dy = mode(diff(centers{i}))/3;
            s = scatter(repelem(x(i),3), [centers{i}(j)-dy centers{i}(j) centers{i}(j)+dy], ...
                dotSize, 'markerfacecolor', accentColor, 'MarkerEdgeColor', 'none');
        else
            scatter(x(i), centers{i}(j), nodeSize, 'markerfacecolor', color*.6, 'markeredgecolor', color, 'LineWidth', 1.5)
        end
    end
end

set(gca, 'color', 'black', 'visible', 'off', 'xlim', [0 1], 'ylim', [0 1], 'position', [0 0 1 1])
print -clipboard -dmeta

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

overlayTop = overlayImgs(imgsTop, {'colors', 'jet', 'contrastLims', [.3 .75], 'cutoff', 100, 'projection', 'mean'});
overlayBot = overlayImgs(imgsBot, {'colors', 'jet', 'contrastLims', [.3 .6], 'cutoff', 100, 'projection', 'mean'});
overlay = cat(1,overlayTop, overlayBot);

close all
figure('position', [200 200 size(overlay,2) size(overlay,1)], 'menubar', 'none');
set(gca, 'position', [0 0 1 1])
imshow(overlay);
print -clipboard -dmeta



