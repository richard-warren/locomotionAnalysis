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

% set step colors
stepSaturation = .9;
stepColors = hsv2rgb([.65 stepSaturation 1;
                      .55 stepSaturation 1;
                      .02 stepSaturation 1;
                      .12 stepSaturation 1]);  % LF, TF, LH, TH