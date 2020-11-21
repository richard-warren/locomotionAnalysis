close all

% settings
paper1_config;
brainRegion = 'mtc';  % mtc, sen
rostralHippoPos = 77/158;  % vertical position of the rostral edge of the hippocampus, expressed as percentage from top of image (this is used to vertically align lesions)
mmWidth = 11;  % width of cartoon in mm
interpPoints = 200;
smoothSmps = 5;
roiFolder = 'Y:\loco\obstacleData\papers\hurdles_paper1\figures\histology\lesionRois';
ventricalScale = 4.4;  % (mm)
baseScale = 2.9089;  % (mm)
figWidth = 400;
alpha = .6;
lineWidth = .5;


% initializations
cartoon = imread(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'histology', 'brainCartoon', ['cartoon_' brainRegion]), 'png');
cartoon = imgaussfilt(cartoon, 1);
whRatio = size(cartoon,2) / size(cartoon,1);

% find x and y grid in mm
mmPerPix = mmWidth/size(cartoon,2);
x = linspace(-mmWidth/2, mmWidth/2, size(cartoon,2));
y = fliplr(1:size(cartoon,1)) - (size(cartoon,1)*rostralHippoPos);
y = y*mmPerPix;

% show image
figure('menubar', 'none', 'Position', [2023.00 436.00 figWidth figWidth/whRatio]);
colormap gray;
image(x, y, cartoon);
set(gca, 'Position', [0 0 1 1], 'Visible', 'off');
set(gca, 'YDir', 'Normal')
hold on

% % show lines (for troubleshooting)
% line([x(1) x(end)], [0 0])
% line([0 0], [y(1) y(end)])


% get list of mice
if strcmp(brainRegion, 'mtc'); sheet='mtcLesionNotes'; elseif strcmp(brainRegion, 'sen'); sheet='senLesionNotes'; end
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', sheet);
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);
mice = cellfun(@upper, mice, 'UniformOutput', false);
colors = lines(length(mice));

% get names of roi .mat files
files = dir(fullfile(roiFolder, '*.mat'));


edges = cell(1,length(mice));
for i = 1:length(mice)
    
    % load mouse roi
    fileBin = contains({files.name}, mice{i});
    if sum(fileBin)>1
        fileBin = fileBin & contains({files.name}, 'Contra');
    end
    load(fullfile(roiFolder, files(fileBin).name))

    vScale = ventricleDistance;
    bScale = bottomDistance;
    shrinkage = mean([vScale/ventricalScale, bScale/baseScale]);
    
    
    x = table2array(ExportPoints(:,1));
    if mean(x)>0; x=-x; end  % display on lesions on the left side of figure
    y = table2array(ExportPoints(:,2));
    y = y - HippocampusY;
    x = x/shrinkage;
    y = y/shrinkage;
    
    % get segment length
    l = [0 cumsum(sqrt(sum(diff([x' x(1); y' y(1)],1,2).^2,1)))];
    x = smooth(x, smoothSmps);
    y = smooth(y, smoothSmps);
    
    % interpolation parameterized by arc length :)
    method = 'makima';
    x = interp1(l, [x; x(1)], linspace(l(1), l(end), interpPoints), method);
    y = interp1(l, [y; y(1)], linspace(l(1), l(end), interpPoints), method);

    patch(x, y, lesionColor, 'FaceAlpha', alpha, 'LineWidth', lineWidth, 'EdgeColor', axisColor);
%     if showLine; edges{i} = plot(x, y, 'LineWidth', lineWidth, 'Color', colors(i,:)); end
end
for i = 1:length(edges); uistack(edges{i}, 'top'); end

file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'histology', ['lesions_', brainRegion '.svg']);
saveas(gcf, file)

