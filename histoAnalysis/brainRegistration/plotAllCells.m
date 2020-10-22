function plotAllCells()
% plot all recorded cells on the allen brain common coordinate framework (ccf)


% load allen brain common coordinate framework
ccf = loadCCF();

% load all cells
files = dir(fullfile(getenv('SSD'), 'paper2', 'histo', 'registration', '*_registration.mat'));
cellLocations = cell(1, length(files));
for i = 1:length(files)
    temp = load(fullfile(files(i).folder, files(i).name), 'cellLocationsCcfMm');
    cellLocations{i} = temp.cellLocationsCcfMm;
end

% make color map and concat all cell locations
cellLocations = cat(1, cellLocations{:});
colors = repelem(lines(3),2,1);

% 2D
figure('color', 'white', 'position', [106.00 251.00 996.00 928.00]);
views = {'ap', 'dv', 'ml'};
inds = {[1 3], [1 2], [2 3]};  % cellLocations inds for 'ap' and 'ml' views

for i = 1:3
    % zoomed out
    subplot(3,2,(i-1)*2 + 1); hold on
    plotLabels2D(ccf.labels, 'dim', views{i}, ...
        'colors', colors, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml)
    plotLabels2D(ccf.coarseLabels, 'dim', views{i}, ...
        'patchArgs', {'facecolor', [0 0 0], 'facealpha', .05, 'edgecolor', [0 0 0], 'edgealpha', .4}, ...
        'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml, 'smoothing', 10)
    scatter(cellLocations(:,inds{i}(1)), cellLocations(:,inds{i}(2)), ...
        30, 'black', 'filled', 'MarkerFaceAlpha', .6);

    % registered
    subplot(3,2,(i-1)*2 + 2); hold on
    plotLabels2D(ccf.labels, 'dim', views{i}, ...
        'colors', colors, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml)  
    scatter(cellLocations(:,inds{i}(1)), cellLocations(:,inds{i}(2)), ...
        30, 'black', 'filled', 'MarkerFaceAlpha', .6);
end
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'brainRegistration', 'allCells2D.png'))


% 3D
figure('color', 'white', 'position', [106.00 251.00 996.00 928.00]);
plotLabels3D(ccf.coarseLabels==2, 'method', 'contours', 'colors', zeros(3,3), 'smoothing', 5, 'contourAlpha', .25, ...
    'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml, 'downSampling', 2, 'slices', {'dv', 'ap'});
plotLabels3D(ccf.labels, 'method', 'boundary', 'colors', colors, 'surfArgs', {'edgecolor', 'none', 'facealpha', .4}, ...
    'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml);
scatter3(cellLocations(:,1), cellLocations(:,2), cellLocations(:,3), ...
        30, 'black', 'filled', 'MarkerFaceAlpha', .6);
set(gca, 'visible', 'off')
savefig(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'brainRegistration', 'allCells3D.fig'))
