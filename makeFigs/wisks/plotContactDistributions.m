function plotContactDistributions

% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};

nosePix = 197; % these two parameters have to be manually defined... this is a bit of a hack
obsPix = 270;
histoYMax = .4;

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
data(length(sessions)) = struct(); % stores trial data for all sessions



% iterate through sessions

for i = 1:length(sessions)
    
    % load session data
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\wiskContactTimes.mat'], ...
        'contactTimes', 'contactPositions', 'contactFrames')
    
    % store data
    sessionInfoInd = find(strcmp(sessionInfo.session, sessions{i}),1,'first');
    data(i).mouse = sessionInfo.mouse{sessionInfoInd};
    data(i).contactPositions = contactPositions;
    data(i).contactFrames = contactFrames;
    
end



% plot results
close all; figure('color', [1 1 1], 'menubar', 'none');

% get all contact positions
allContactPositions = {data.contactPositions};
allContactPositions = cat(2, allContactPositions{:});

% get average contact frame
allContactFrames = {data.contactFrames};
allContactFrames = cat(3, allContactFrames{:});
meanFrame = uint8(nanmean(allContactFrames,3));

% trim frame
meanFrame = meanFrame(1:220, :);

% set up conversion from image to histogram x scale
pixToMmMapping = polyfit([nosePix obsPix], [0 nanmedian(allContactPositions)], 1);
histoXMin = 1*pixToMmMapping(1) + pixToMmMapping(2);
histoXMax = size(data(1).contactFrames,2)*pixToMmMapping(1) + pixToMmMapping(2);
histoLims = [histoXMin histoXMax] * 1000;

% plot (in mm)
colormap(gray)
imagesc(flipud(meanFrame(:,:)), 'XData', histoLims, 'YData', [0 histoYMax]); hold on
histogram(allContactPositions*1000, 10, 'normalization', 'probability', 'linewidth', 1, ...
    'facealpha', .6, 'facecolor', [244 66 66] / 255); hold on

% pimp fig
set(gca, 'box', 'off', 'xdir', 'reverse', 'ydir', 'normal', 'tickdir', 'out')
set(gca, 'dataaspectratio', [(size(meanFrame,1)/size(meanFrame,2)) * (abs(diff(histoLims))/histoYMax) 1 1])
ax = gca;
ax.YAxis.Visible = 'off';
xlabel('distance from nose (mm)')

saveas(gcf, [getenv('OBSDATADIR') 'figures\wiskContactDistributions.png']);






