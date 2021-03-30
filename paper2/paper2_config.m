% define global parameters for paper2


cfg.font = 'Arial';
cfg.fontsize = 8;
cfg.axArgs = {'TickDir', 'out', 'box', 'off', ...
    'FontName', cfg.font, 'FontSize', cfg.fontsize};

% colors
cfg.nucleusColors = lines(3);  % dentate, interpositus, fastigial
cfg.heatmapColors = customcolormap([0 .5 1], [1 .2 .2; 1 1 1; .2 .2 1]);  % blue -> white -> red
cfg.lickColor = [48 135 227]*.75 / 255;  % from paper1
cfg.velColor = [1 .3 0];
cfg.obsColor = [1 .7 0];  % from paper1
cfg.wiskColor = cfg.obsColor * .5;
cfg.diiColor = [.99 .33 1];  % dii color for probe traces