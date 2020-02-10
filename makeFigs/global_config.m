% global settings for paper


% global
axisColor = [.15 .15 .15];  % use this for black
obsColor = [188 125 181] / 255;
obsAlpha = .15;
waterColor = [48 135 227]*.75 / 255;
ctlStepColor = [.5 .5 .5];
barProperties = {'scatterAlpha', .5, 'barAlpha', .4, 'labelSizePerFactor', .1, ...
                 'lineThickness', 2, 'scatterColors', 'lines', 'connectDots', true, ...
                 'lineAlpha', .05, 'showBars', true};


% step type colors (leading, lagging, fore, hind)
stepSaturation = .9;
stepColors = hsv2rgb([.65 stepSaturation 1;
                      .55 stepSaturation 1;
                      .02 stepSaturation 1;
                      .12 stepSaturation 1]);  % LF, TF, LH, TH


% sensory dependence colors
% colorWisk = [51 204 255]/255;
% colorVision = [255 221 21]/255;
colorWisk = hsv2rgb([.05 1 1]);
colorVision = obsColor;
colorNone = [.2 .2 .2];
% both = mean([colorWisk;colorVision],1);
both = hsv2rgb(mean([rgb2hsv(colorWisk); rgb2hsv(colorVision)],1));
sensColors = [both; colorWisk; colorVision; colorNone];


% decision making
m.deltaMin =.5;
m.lightOffOnly = false;
m.modPawOnlySwing = true;
m.successOnly = true;
m.predictorsAll = {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'modPawX', 'modPawXVel', 'modPawZ', 'modPawZVel'};
m.predictors = {'modPawX', 'obsHgt', 'velAtWiskContact', 'wiskContactPosition', 'modPawXVel'};
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


