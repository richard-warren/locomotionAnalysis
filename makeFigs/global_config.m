% global settings for paper


% global
axisColor = [.15 .15 .15];  % use this for black
% obsColor = [188 125 181] / 255;
obsColor = [1 .7 0];
obsAlpha = .15;
waterColor = [48 135 227]*.75 / 255;
ctlStepColor = [1 1 1] * .4;
barProperties = {'showBars', false, 'showErrorBars', false, 'lineThickness', 4, 'connectDots', true, ...
                 'lineAlpha', .1, 'showScatter', true, 'scatterColors', 'lines', 'scatterSize', 40, ...
                 'scatterAlpha', .8, 'scatterCondColor', true, 'constantEdgeColor', [], 'barWidth', .8};
sessionPlotProperties = {'colors', 'lines', 'alpha', .3, 'scatSize', 0, 'meanColor', axisColor};


% step type colors (leading, lagging, fore, hind)
stepSaturation = .9;
stepColors = hsv2rgb([.65 stepSaturation 1;
                      .55 stepSaturation 1;
                      .02 stepSaturation 1;
                      .12 stepSaturation 1]);  % LF, TF, LH, TH


% sensory dependence colors
colorWisk = [188 125 181] / 255;
colorVision = obsColor;
colorNone = [.2 .2 .2];
colorBoth = hsv2rgb(mean([rgb2hsv(colorWisk); rgb2hsv(colorVision)],1));
sensColors = [colorBoth; colorWisk; colorVision; colorNone];


% decision making
m.deltaMin = .005;  % (m) minimum change in step length for inclusion in model
m.lightOffOnly = false;
m.modPawOnlySwing = true;
m.successOnly = false;  % must set to false for sensoryDependence, bc not enough good trials with no whiskers and no light
m.modSwingContactsMax = 4;  % first swing of first modified paw cannot have more than this many frames of contact with the obstacle // overall success is LESS THAN 5 frames

m.predictorsAll = {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'modPawX', 'modPawXVel', 'modPawZ', 'modPawZVel'};
m.predictorsNamedAll = {{'wheel velocity'}, {'body angle'}, {'obstacle height'}, {'obstacle position';'(horizontal)'}, {'paw position';'(horizontal)'}, {'paw velocity'; '(horizontal)'}, {'paw position';'(vertical)'}, {'paw velocity'; '(vertical)'}};
m.predictors = {'wiskContactPosition', 'modPawX', 'velAtWiskContact', 'modPawXVel', 'modPawZ', 'modPawZVel', 'obsHgt', 'angleAtWiskContact'};  % the values in this array are determined via forward selection in the scipt baselineDecision.m
[~, inds] = ismember(m.predictors, m.predictorsAll);
m.predictorsNamed = m.predictorsNamedAll(inds);

m.heatmapNormalize = 'col';  % normalize row or colum to sum to 1
decisionColors = flipud(colorme(2, 'offset', .2, 'showSamples', false)); % first entry is small step, second is big step
preDecisionColor = hsv2rgb(mean(rgb2hsv(decisionColors),1));
modelColor = [56, 222, 53] / 255;


% manipulations
lesionColor = [177 13 4]/255;
muscimolColor = [13 177 111]/255;


% other
contrast = [0 .75];  % for videos and images
contactColor = [1 .2 .2];  % contact color (paw, whiskers)


