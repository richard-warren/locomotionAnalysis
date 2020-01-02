close all

% settings
brainRegion = 'sen';  % mtc, sen
rostralHippoPos = .58; %50/101;  % vertical position of the rostral edge of the hippocampus, expressed as percentage from top of image (this is used to vertically align lesions)
mmWidth = 11;  % width of cartoon in mm
interpPoints = 200;
smoothSmps = 5;
roiFolder = 'Y:\obstacleData\papers\hurdles_paper1\figures\histology\lesionRois';
ventricalScale = 4.4;  % (mm)
baseScale = 2.9089;  % (mm)
figWidth = 400;


% initializations
cartoon = imread(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'histology', 'brainCartoon', ['cartoon_' brainRegion]), 'png');
% cartoon = uint8(cartoon(:,:,1)~=0) * 255;
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



for i = 1:length(mice)
    
    % load mouse roi
    fileBin = contains({files.name}, mice{i});
    if sum(fileBin)>1
        fileBin = fileBin & contains({files.name}, 'Contra');
    end
    load(fullfile(roiFolder, files(fileBin).name))

    vScale = table2array(Combine(end,1));
    bScale = table2array(Combine(end-1,1));
    shrinkage = mean([vScale/ventricalScale, bScale/baseScale]);
    
    
    x = table2array(Combine(1:end-2,1));
    disp(mean(x))
    if mean(x)>0; x=-x; end  % display on lesions on the left side of figure
    y = table2array(Combine(1:end-2,2));
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

    % close all; figure('Position', [1980.00 428.00 560.00 420.00])
    fill(x, y, colors(i,:), ...
        'FaceAlpha', .4, 'LineWidth', 1)
end

file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'histology', ['lesions_', brainRegion '.svg']);
saveas(gcf, file)

