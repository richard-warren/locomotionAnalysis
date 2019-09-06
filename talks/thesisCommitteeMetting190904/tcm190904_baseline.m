%% GLOBAL SETTINGS

% first load data // each dataset has a corresponding flattened data set //
% incidual plots may take subsets of the flattened dataset

axisColor = [.95 .95 .95];
pawColors = hsv(4);
mouseColors = 'jet';
ctlStepColor = [1 1 1]*.9;
smallBigStepColors = [0.850 0.325 0.098; 0 0.447 0.741]; % first entry is small step, second is big step;

%% INITIALIZATIONS

% load data
fprintf('loading data... '); 
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data');
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition', 'obsHgt', 'isPawSuccess', ...
    'isLeading', 'isFore', 'stepOverKinInterp', 'paw', 'controlStepKinInterp', 'preObsHgt', 'stepType', ...
    'isBigStep'});
disp('data loaded!')

% conditions and conditionals
conditionNames = {{'fore', 'hind'}, {'leading', 'trailing'}};
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});

conditionals.lightOn = struct('name', 'isLightOn', 'condition', @(x) x==1);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.none = struct('name', '', 'condition', @(x) x); % no conditionals



% build decision making model
iterations = 20;
normalizeData = true;
barColors = [0 .24 .49; .6 .75 .23; .2 .2 .2];

accuracies = nan(3,iterations);

for i = 1:iterations
    
    predDistBin = ismember(predictorNames, 'modPawPredictedDistanceToObs'); % this var will ONLY be included in the GLM, not the neural net
    
    % NEURAL NET
    [X, y, predictorNames, isCategorical] = ...
        prepareDecisionModelData(flat, 'all', 'isBigStep', true, referenceModPaw, normalizeData, {'balanceClasses', true});
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
    glm = fitglm(X([tr.trainInd tr.valInd], predDistBin), y([tr.trainInd tr.valInd]), ...
        'Distribution', 'binomial', 'CategoricalVars', isCategorical(predDistBin));
    accuracies(2,i) = mean(round(predict(glm, X(tr.testInd,colBins)))==y(tr.testInd));
end

disp('all done!')


%% SPEED VS POSITION

% initializations
figure('units', 'inches', 'position', [3.83 3.61 6.67 3.09], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
yLims = [0 .8];

% add obstacle rectangle and lines
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
line([0 0], yLims, 'linewidth', 2, 'color', axisColor)
line([x(1) x(1)], yLims, 'linewidth', 2, 'color', axisColor)
text(x(1), yLims(2), '\itobstacle on', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
text(0, yLims(2), '\itobstacle reached', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

% plot
plotDvPsth(flat, 'velVsPosition', [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [1 1 1], 'mouseColors', mouseColors, ...
    'errorFcn', @(x) nanstd(x), 'xlim', [-.5 .2], 'mouseAlpha', .5, 'errorAlpha', .4})

% pimp fig
set(gca, 'YLim', [0 .8]);
xlabel('distance (m)')
ylabel('velocity (m/s)')
set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black')
print -clipboard -dmeta

%% BASELINE SUCCESS

figure('units', 'inches', 'position', [3.83 3.61 3.92 3.09], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)

logPlotRick([flat.obsHgt]*1000, [flat.isPawSuccess], ...
    {'colors', pawColors, 'conditions', [flat.stepType], 'xlabel', 'obstacle height', 'ylabel', 'success rate', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'computeVariance', false, 'ylim', [0 1], 'plotMouseErrors', false, 'lineWidth', 3})
set(gca, 'xlim', [3.4 10], 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor)

print -clipboard -dmeta

%% BASELINE KINEMATICS

% settings
obsHgtBins = 1; % discretize obstacle heights into this many bins

figure('units', 'inches', 'position', [4.27 6.41 10.71 1.99], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
colNames = {'hind', 'fore'};
rowNames = {'trailing', 'leading'};
conditionSequence = [4 2 3 1]; % mapping between plot index and condition sequence, which is specified above
xLims = [-.05 .05];
yLims = [0 .016];
fading = .5; % within a plot, each paw's color fades from fading*color -> color
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLightOn]);


if obsHgtBins==1; lineWidth = 4; else lineWidth = 3; end

% get kin data
kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);
kinData = cat(1, kinData, permute(cat(3, flat.controlStepKinInterp), [3,1,2])); % append with control steps temporarily

% flip y values s.t. leading is always right and lagging always left
bins = [flat.paw]==1 & [flat.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat.paw]==2 & [flat.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat.paw]==3 & ~[flat.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat.paw]==4 & ~[flat.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);

% split real and control steps
kinDataCtl = kinData(length(flat)+1:end,:,:);
kinData = kinData(1:length(flat),:,:);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1);  % change the x starting x position of ctl steps to match steps over

plotInd = 1;
buf = .03;
wid = (1-3*buf)/2;
hgt = (1-buf)/2;
positions = [buf hgt; wid+2*buf hgt; buf 0; wid+2*buf 0];
for i = conditionSequence
    subplot(2,2,plotInd)
    plotColor = repmat(pawColors(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    bins = [flat.stepType]==i;

    % plot control step
    plotKinematics(kinDataCtl(bins,[1,3],:), [flat(bins).obsHgt], ones(1,sum(bins)), ...
        {'colors', ctlStepColor, 'obsAlpha', 0, 'lineAlpha', .8, 'lineWidth', lineWidth}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], [flat(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat(bins).mouse}, 'lineColor', axisColor, 'lineWidth', lineWidth}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims, 'color', 'black')
    set(gca, 'position', [positions(plotInd,:) wid hgt])
    
    if plotInd==1
        text(0, yLims(2), colNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', axisColor);
        text(xLims(1), mean(yLims), rowNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'color', axisColor);
    elseif plotInd==2
        text(0, yLims(2), colNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', axisColor);
    elseif plotInd==3
        text(xLims(1), mean(yLims), rowNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'color', axisColor);
    end

    plotInd = plotInd+1;
end

print -clipboard -dmeta


%% PRE OBS HEIGHT

figure('units', 'inches', 'position', [4.46 1.36 3.8 2.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
figVars = [vars.isFore; vars.isLeading];
figConditionals = [conditionals.none];

matMod = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlPreObsHgt', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlStepColor,8,1);
colorsWithBl(2:2:8,:) = pawColors;
barFancy(mat, 'axisColor', axisColor, 'levelNames', {figVars.levelNames}, 'ylabel', 'paw height (mm)', ...
    'colors', colorsWithBl, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .1, 'lineThickness', 2)
print -clipboard -dmeta

%% STEP OVER LENGTH

figure('units', 'inches', 'position', [4.46 1.36 3.8 2.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
figVars = [vars.isFore; vars.isLeading];
figConditionals = [conditionals.none];

matMod = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlStepLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlStepColor,8,1);
colorsWithBl(2:2:8,:) = pawColors;
barFancy(mat, 'axisColor', axisColor, 'levelNames', {figVars.levelNames}, 'ylabel', 'step length (mm)', ...
    'colors', colorsWithBl, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .1, 'lineThickness', 2)

print -clipboard -dmeta

%% SHOW LEADING LAGGING SCHEMATIC

showLeadingLaggingImg('190318_000', 43, pawColors)
print -clipboard -dmeta

%% SHAPING CORRELATIONS
figure('units', 'inches', 'position', [13.22 3.65 3.19 3.45], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)

figVars = [vars.isFore; vars.isLeading];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, figVars, {'mouse'}, {'session'}, [conditionals.isLe], 'corr');
barFancy(dv, 'axisColor', axisColor, 'levelNames', {figVars.levelNames}, 'ylabel', 'paw-obstacle correlation', ...
    'colors', pawColors, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .15, 'lineThickness', 2, 'barWidth', .75)
set(gca, 'position', [.15 .11 .77 .81])
print -clipboard -dmeta

%% SHAPING MOVING AVGS

figure('units', 'inches', 'position', [4.46 1.36 3.8 2.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
logPlotRick([flat.obsHgt]*1000, [flat.preObsHgt]*1000, ...
    {'colors', pawColors, 'conditions', [flat.stepType], 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)',  ...
    'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'lineWidth', 2})
set(gca, 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor)
print -clipboard -dmeta

%% DECISION MAKING MODEL BARS

figure('units', 'inches', 'position', [13.22 3.65 3.19 3.45], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
barPlotRick(accuracies, {'conditionNames', {{'NN', 'GLM', 'shuffled'}}, ...
    'lineThickness', 2, 'addBars', true, 'scatColors', 'hsv', 'scatAlpha', .5, 'showStats', false, ...
    'ylim', [0 1], 'ytick', [0 .5 1], 'ylabel', 'accuracy', 'conditionColors', barColors})
set(gca, 'Position', [0.15 0.2 0.75 0.79])
print -clipboard -dmeta

%% DECISION MAKING KINEMATICS

% settings
netNum = 1;
binNum = 5;
pctileBins = true;

% initializations
kin = permute(cat(3,flat.modPawKinInterp{:}), [3,1,2]);
kinCtl = permute(cat(3,flat.preModPawKinInterp{:}), [3,1,2]);
[X, y, predictorNames, isCategorical] = ...
        prepareDecisionModelData(flat, predictors{netNum}, outcome, useAllPaws(netNum), referenceModPaw, normalizeData, ...
        {'balanceClasses', false, 'removeNans', false});

% choose binning variable
binVar = nnets{netNum}(X'); % neural network output
% binVar = -flat.x_paw2; % position of mod paw at moment of contact
% binVar = flat.velAtWiskContact; % vel
% binVar = flat.modStepStart_paw2; % starting position of mod paw
% binVar = flat.modPawPredictedDistanceToObs;

% perctentile bins    
if pctileBins
    binEdges = prctile(binVar, linspace(0, 100, binNum+1));
% evenly spaced bins
else
    binEdges = linspace(min(binVar), max(binVar), binNum+1);
end
condition = discretize(binVar, binEdges);

figure('units', 'inches', 'position', [13.22 3.65 3.19 3.45], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
plotBigStepKin(kin(:,[1,3],:), kinCtl(:,[1,3],:), flat.obsHgt, condition, flat.isBigStep, ...
    {'colors', stepTypeColors, 'xLims', [-.09 .04], 'addHistos', false, 'lineWid', 3, ...
    'contactInds', flat.contactInd, 'histoHgt', .5, 'showSmpNum', false})
print -clipboard -dmeta


%% HEATMAP

% settings
xLims = [-.03 .015];
yLims = [-.03 .03];

figure('units', 'inches', 'position', [13.22 3.65 3.19 3.45], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
heatmapRick(flat.modPawPredictedDistanceToObs, flat.modPawDistanceToObs, ...
    {'xLims', xLims, 'yLims', yLims, ...
    'xlabel', 'predicted distance to obstalce (m)', 'ylabel', 'distance to obstalce (m)'})
set(gca, 'DataAspectRatio', [1 1 1])
line(xLims, xLims, 'color', [.5 .5 1 .8], 'lineWidth', 3)
print -clipboard -dmeta









