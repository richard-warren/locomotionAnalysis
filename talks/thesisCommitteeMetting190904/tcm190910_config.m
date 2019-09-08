% global settings for thesis committee meeting presentation

% sizes
figHgt = 3.75;  % inches
font = 'Calibri';
fontSize = 12;

% global, global
mouseColors = 'jet';
axisColor = [.95 .95 .95];  % this is the 'white' to be used
barProperties = {'axisColor', axisColor, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .1, ...
                 'lineThickness', 2, 'scatterColors', mouseColors, 'FontName', font, 'FontSize', fontSize};

% colors, by experiment
pawColors = colorme(4, 'offset', .65, 'showSamples', false);
sensoryColors = [colorme(3, 'offset', .70, 'showSamples', false); axisColor*.5];
stepTypeColors = colorme(2, 'offset', .9, 'showSamples', false);
modelColors = [.5 1 0; axisColor];
manipColor = [colorme(1, 'offset', .50, 'showSamples', false)];
optoColors = winter(2);