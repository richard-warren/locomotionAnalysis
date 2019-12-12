% global settings for paper

% sizes
% figHgt = 3.75;  % inches
% font = 'Calibri';
% fontSize = 12;

% global
axisColor = [.15 .15 .15];  % use this for black
obsColor = [188 125 181] / 255;
obsAlpha = .15;
waterColor = [48 135 227]*.75 / 255;
ctlStepColor = [.5 .5 .5];
barProperties = {'scatterAlpha', .5, 'barAlpha', .4, 'labelSizePerFactor', .1, ...
                 'lineThickness', 2, 'scatterColors', 'lines', 'connectDots', true, ...
                 'lineAlpha', .05};

% set step colors
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