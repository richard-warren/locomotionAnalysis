% overlaid tracking

session = '200709_000';
contrast = [.2 .8];
trial = 40;
% colorsTemp = stepColors([4 2 1 3], :);  % this is a hack that makes the colors align with leading, lagging, fore, hind conditions

showSingleFrameTracking(session, trial, ...
    'contrastLims', contrast, 'addWiskCam', true, 'pawColors', lines(4));

% file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'tracking.png');
% fprintf('writing %s to disk...\n', file);
% saveas(gcf, file)

%%