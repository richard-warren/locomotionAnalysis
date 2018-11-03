%% plot frame dimensions vs accuracy


% COLLECT DATA

% settings
folder = 'D:\dlcBenchmarking\sizevsspeed\';

% initializations
files = dir(fullfile(folder, '*.csv'));
dims = cell2mat(cellfun(@(x) [str2num(x(5:8)); str2num(x(10:13))], ...
    {files.name}, 'UniformOutput', false)); % parse file name to get video dims
baseColor = [100 100 255] / 255;
darkening = .0;



data = readtable(fullfile(folder, files(1).name));
frameNum = height(data);
features = data(:,2:3:width(data)).Properties.VariableNames;


locations = nan(length(dims), length(features), 2, frameNum); % vid dims X feature number X x/y X frame number
scores = nan(length(dims), length(features), frameNum);

for i = 1:length(dims)
    
    data = readtable(fullfile(folder, files(i).name)); % load table for vid

    for j = 1:length(features)
        xy = table2array(data(:, (j-1)*3+2 : (j-1)*3+3));
        xy = xy .* (dims(1,end)/dims(1,i)); % rescale to original size
        locations(i,j,:,:) = xy';
        scores(i,j,:) = table2array(data(:, (j-1)*3+4));
    end
end

%% PLOT TRACE EXAMPLES

close all;
figure('Color', 'white', 'MenuBar', 'none', 'Position', [2020 400 650 350]);
offset = 100;
colorMap = interp2(1:3, [0;1], cat(1,baseColor,baseColor*darkening), 1:3, linspace(0,1,256)');
colors = interp2(1:3, [0;1], cat(1,baseColor,baseColor*darkening), 1:3, linspace(0,1,length(dims))');

part = 1;

for i = 1:length(dims)
    trace = squeeze(locations(i,part,1,:));
    plot(trace + offset*(length(dims)-i), 'Color', colors(i,:), 'LineWidth', 2); hold on;
end

% pimp fig
set(gca, 'Visible', 'off')
colormap(colorMap)
bar = colorbar('Direction', 'reverse', 'Box', 'off', 'Ticks', [0 1], ...
    'TickLabels', cellfun(@(x) sprintf('%ix%i', x(1), x(2)), num2cell(dims(:,[1 end]),1), 'UniformOutput', false));
saveas(gcf, 'F:\Google Drive\other\dlcSpeedPaper\sizeVsAccuracyTraces.pdf');

%% PLOT ACCURACY VS FRAME DIMS


bl = locations(end,:,:,:); % 'ground truth' is highest resoultion image
distances = sqrt(squeeze(sum((locations-bl).^2,3)));
distancesMean = mean(distances,3);
tranparency = 0.6;

close all;
figure('Color', 'white', 'MenuBar', 'none', 'Position', [2020 400 520 320]);
shapes = ['o','+','x','h'];

for i = 1:size(distancesMean,2)
    if ismember(shapes(i), {'o', 'h'})
        scatter(dims(1,:), distancesMean(:,i), ...
            100, colors, shapes(i), 'filled', 'MarkerFaceAlpha', tranparency); hold on
    else
        scatter(dims(1,:), distancesMean(:,i), ...
            100, colors, shapes(i), 'MarkerFaceAlpha', tranparency); hold on
    end
end

% pimp fig
xlabel('frame size [pixels]')
ylabel('mean euclidian error [pixels]')
set(gca, 'XTick', dims(1,:), 'XTickLabelRotation', 45, ...
    'XTickLabel', cellfun(@(x) sprintf('%ix%i', x(1), x(2)), num2cell(dims,1), 'UniformOutput', false))
saveas(gcf, 'F:\Google Drive\other\dlcSpeedPaper\sizeVsAccuracySummary.pdf');
%%
% OTHER THINGS
%

%% show tracking results
vid = VideoReader('D:\dlcBenchmarking\sizevsspeed\run0400.avi');


parts = 4;
inds = find(trace>250,10,'first')-1;

for frameToShow = inds'
%     frameToShow = 2571;
    

    colors = hsv(length(parts));
    frame = read(vid, frameToShow);
    figure;
    imshow(frame); hold on

    scatter(locations(end,parts,1,frameToShow), locations(end,parts,2,frameToShow), 50, colors, 'filled');
    pause(1)
end



%% test speed as a function of squareness of image


% settings
fileName = 'D:\dlcBenchmarking\results-ricksystem\padtests\data.csv';
varyingDim = 'width';

% initializations
height = unique(data.height);
data = readtable(fileName);
dims = unique(data.(varyingDim));
repetitions = sum(data.(varyingDim)==dims(1)); % number of reps for each condition
xLabels = cellfun(@(x) sprintf('%ix%i', height, x), num2cell(dims), 'UniformOutput', false);
colors = winter(repetitions);

close all; figure('color', 'white', 'menubar', 'none');

speeds = nan(length(dims), repetitions);
for j = 1:length(dims)
    bins = data.(varyingDim)==dims(j);
    speeds(j,:) = data.nframes(bins) ./ data.run_duration(bins);
    
    scatter(repmat(dims(j),1,repetitions), speeds(j,:), 50, colors, 'filled', 'MarkerFaceAlpha', .6); hold on
end

% plot means
line(dims, mean(speeds,2), 'color', [.2 .2 .2], 'linewidth', 2);


% pimp fig
set(gca, 'XTick', dims, 'XTickLabel', xLabels, 'XTickLabelRotation', 45)
xlabel('frame dimensions')
ylabel('fps')


